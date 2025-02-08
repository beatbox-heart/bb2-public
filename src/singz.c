/**
 * Copyright (C) (2010-2025) Vadim Biktashev, Irina Biktasheva et al. 
 * (see ../AUTHORS for the full list of contributors)
 *
 * This file is part of Beatbox.
 *
 * Beatbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beatbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Beatbox.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Finding "tips" - intersections of isolines in a z-section */

#if MPI
#include <mpi.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#ifdef X11
#include "windraw.h"
#endif

#if MPI
extern int idev;
#endif


#if MPI
#define debug(fmt,...) printf (fmt,t,mpi_rank, __VA_ARGS__)
#else
#define debug(fmt,...) printf (fmt,t,0, __VA_ARGS__)
#endif


typedef struct {
  int v0, v1;			/* source layers */
  real c0, c1;			/* isoline constants */
  int K[2];			/* which is which */
  int precise_orientation;	/* if nonzero, phi and psi are done via higher-order grads */
  char filename[MAXPATH];	/* output file name */
  int append;			/* file mode */
  FILE *file;			/* file handle for output; NULL in MPI in non-root procs */
  int fileout;			/* file output flag; may be non-zero even in non-root procs */
  int everypoint;		/* if nonzero, put time stamp at every point */
  int everysection;		/* if nonzero, put time stamp at every z-section */
  int everyrecord;		/* if nonzero, put time stamp only at every t-record */
  char pointsep[80];		/* separator between points */
  char sectionsep[80];		/* separator between z-sections */
  char recordsep[80];		/* separator between t-records */
  #define SINGZBUFMAX (4*1024*1024-1) /* 4M-1, reserve one byte in the end for trailing zero */
  char buf[SINGZBUFMAX+1];	/* output buffer: should be long enough to print whole record incl trailing zero */
  char *p;			/* pointer showing after the end of the record in the buffer so far */
  #if MPI
  MPI_Comm comm;		/* communicator for this device instance */
  int tag;			/* tag for messages of this device instance (don't need if private communicator?) */
  int root;			/* the rank of the instance that does the file output */
  #endif
  #define MAXSTACK 10
  real xtip[MAXSTACK],ytip[MAXSTACK],ztip[MAXSTACK]; /* array to store singular points found at previous/next step */
  int ntip;

  /* Addresses of k-variables to store statistics of found tips */
  KREAL(Npoints); /* better be made integer but might clash with older scripts */
  KREAL(xmean);
  KREAL(ymean);
  KREAL(zmean);
  KREAL(xstddev);
  KREAL(ystddev);
  KREAL(zstddev);
  /* Integer flag showing whether k-variables are used */
  int k_output;
#ifdef X11
  BGIWindow wnd;
#endif
} STR;


/* Add the message to the buffer of this process */
static int overflow_reported=0;
static void output(STR *S,char *fmt, ...) {
  va_list ap;
  int l, n, nmax;
/*   #if MPI */
/*   printf("t=%d rank=%d output fmt=%s\n",t,mpi_rank,fmt); */
/*   #endif */
  if NOT(S->fileout) return;
  
  l = (S->p) - (&(S->buf[0])) ;			/* how many characters already are in the buffer */
  nmax = SINGZBUFMAX - l;			/* how many characters are still available */
  if (nmax>0) {
    va_start(ap, fmt);
    n=vsnprintf(S->p,nmax,fmt,ap);		/* how many characters would be needed */
    va_end(ap);
    if (n<nmax) {				/* enough space, proceed as normal */
      S->p += n;
    } else {
      S->p += nmax;				/* only output what's fit in the buffer */
      if (overflow_reported==0) {
	MESSAGE("singz buffer overflow\n");	/* and repeated the buffer overflow */
	overflow_reported=1;			/* only report it once */
      }
    }
    *(S->p)='\0';				/* don't forget the trailing zero */
  }
}

/* Output buffers of all processes to the file, via the root process */
static void flush_buffer(STR *S) {
  if NOT(S->fileout) return;
#if MPI
  int rank;
  if (mpi_rank == S->root) {			/* If I am the reporting process */
    fputs(S->buf,S->file); 			/* first print out my own buffer */
    for (rank=0; rank<num_active_procs; rank++) { /* then get all other buffers */
      if (rank==S->root) continue;
      MPIDO(MPI_Recv(S->buf,SINGZBUFMAX,MPI_CHAR,rank,S->tag,S->comm,MPI_STATUS_IGNORE),
	    "singz buffer receipt failure");
      fputs(S->buf,S->file); 			/* and output them too. */
    }
    fflush(S->file);				/* always flush the file after writing to it. */
  } else { 					/* Else just send my buffer to the reporting proc */
    MPIDO(MPI_Send(S->buf,(S->p-S->buf+1),MPI_CHAR,S->root,S->tag,S->comm),
	  "singz buffer send failure");
  }
  S->p=&(S->buf[0]);				/* In any case reset my buffer. */
  *(S->p)='\0';
#else
  fputs(S->buf,S->file);			/* Sequential: just output the buffer */
  fflush(S->file);				/* always flush the file after writing to it */
  S->p=&(S->buf[0]);				/* and reset the buffer. */	
  *(S->p)='\0';
#endif
}

#define u(k,x,y) New[ind(x,y,z,S->K[k])]
static void util(int x, int y, real dx, real dy, STR *S,Space *s,int z) {
  if (S->fileout) {
    real ux, uy, vx, vy, phi, psi, cosangle;
    if (S->precise_orientation) {
      real ux00=u(0,x+1,y)-u(0,x-1,y);
      real ux10=u(0,x+2,y)-u(0,x  ,y);
      real ux01=u(0,x+1,y+1)-u(0,x-1,y+1);
      real ux11=u(0,x+2,y+1)-u(0,x  ,y+1);
      real uy00=u(0,x,y+1)  -u(0,x  ,y-1);
      real uy01=u(0,x,y+2)  -u(0,x  ,y  );
      real uy10=u(0,x+1,y+1)-u(0,x+1,y-1);
      real uy11=u(0,x+1,y+2)-u(0,x+1,y  );
      real vx00=u(1,x+1,y)-u(1,x-1,y);
      real vx10=u(1,x+2,y)-u(1,x  ,y);
      real vx01=u(1,x+1,y+1)-u(1,x-1,y+1);
      real vx11=u(1,x+2,y+1)-u(1,x  ,y+1);
      real vy00=u(1,x,y+1)  -u(1,x  ,y-1);
      real vy01=u(1,x,y+2)  -u(1,x  ,y  );
      real vy10=u(1,x+1,y+1)-u(1,x+1,y-1);
      real vy11=u(1,x+1,y+2)-u(1,x+1,y  );
      ux=ux00*(1-dx)*(1-dy)+ux10*dx*(1-dy)+ux01*(1-dx)*dy+ux11*dx*dy;
      uy=uy00*(1-dx)*(1-dy)+uy10*dx*(1-dy)+uy01*(1-dx)*dy+uy11*dx*dy;
      vx=vx00*(1-dx)*(1-dy)+vx10*dx*(1-dy)+vx01*(1-dx)*dy+vx11*dx*dy;
      vy=vy00*(1-dx)*(1-dy)+vy10*dx*(1-dy)+vy01*(1-dx)*dy+vy11*dx*dy;
    } else {
      real u00=u(0,x,y);
      real u10=u(0,x+1,y);
      real u01=u(0,x,y+1);
      real u11=u(0,x+1,y+1);
      real v00=u(1,x,y);
      real v10=u(1,x+1,y);
      real v01=u(1,x,y+1);
      real v11=u(1,x+1,y+1);
      ux=(u10-u00)*(1-dy)+(u11-u01)*dy;
      uy=(u01-u00)*(1-dx)+(u11-u10)*dx;
      vx=(v10-v00)*(1-dy)+(v11-v01)*dy;
      vy=(v01-v00)*(1-dx)+(v11-v10)*dx;
    }
    phi=atan2(uy, ux);
    psi=atan2(vy, vx);
    cosangle=(ux*vx+uy*vy)/sqrt((ux*ux+uy*uy)*(vx*vx+vy*vy));
    if (S->everypoint) output(S,"%10ld",t);
    if (s->z0==s->z1) 
      output(S," %f %f %f %f %f",
	      (float)(x+dx),(float)(y+dy),(float)phi,(float)psi,(float)cosangle);
    else 
      output(S," %f %f %f %f %f %f",
	      (float)(x+dx),(float)(y+dy),(float)(z),(float)phi,(float)psi,(float)cosangle);
    output(S,"%s",S->pointsep);
  }
  if (S->ntip < MAXSTACK) {
    S->xtip[S->ntip]=x+dx;
    S->ytip[S->ntip]=y+dy;
    S->ztip[S->ntip]=z;
    S->ntip++;
  }
}

static int isoline (int k, int x, int y, real fc, real *fx, real *fy, STR *S, int z) {
  static int dx[5]   = {0,  1,  1,  0,  0}; /* x's of cell corners */
  static int dy[5]   = {0,  0,  1,  1,  0}; /* y's of cell corners */
  int i;				    /* num of corner 0..4, 4==0 */
  int j;				    /* num of found crossections 0..2 */
  real fprev, fnext;

  for (i=1,j=0;i<5&&j<2;i++) {
    fprev = u(k,x+dx[i-1],y+dy[i-1])-fc;
    fnext = u(k,x+dx[i],y+dy[i])-fc;
    if (!fnext) continue;     		/* zero corner - next cycle will do */
    if (fnext*fprev > 0) continue;	/* no touches on this rib */
    fx[j] = (fnext*dx[i-1]-fprev*dx[i])/(fnext-fprev);
		if (fx[j]<0 || fx[j]>1)
		  ABORT("ogo");
    fy[j] = (fnext*dy[i-1]-fprev*dy[i])/(fnext-fprev);
		if (fy[j]<0 || fy[j]>1)
		  ABORT("aga");
    j++;
  }
  if (j==2) return 1;
  return 0;		/* j==1 also possible and not interesting */
/*
  if (j==0) return 0;
  ABORT("strange cell: j=%d",j);
*/
}

static int isocross(STR *S,Space *s, int z) {
  int x, y;
  real ux[2], uy[2], vx[2], vy[2];
  real D,tu,tv,dx,dy,a11,a12,a21,a22,d1,d2;
  int uok, vok;

  /* TODO: make this listpointable? */
  /* The space is adjusted by prior choice of x1, y1 */
  for (y=s->y0;y<=s->y1;y++)
  for (x=s->x0;x<=s->x1;x++) {
    if (!isoline(0,x,y,S->c0,ux,uy,S,z)) continue;
    if (!isoline(1,x,y,S->c1,vx,vy,S,z)) continue;
    #define Det(a,b,c,d) ((a)*(d)-(b)*(c))
    a11 = ux[0]-ux[1];
    a12 = vx[1]-vx[0];
    a21 = uy[0]-uy[1];
    a22 = vy[1]-vy[0];
    d1 = vx[1]-ux[1];
    d2 = vy[1]-uy[1];
    D  = Det (a11,a12,a21,a22);
    if (D==0.0) continue;
    tu = Det (d1,a12,d2,a22)/D;
    tv = Det (a11,d1,a21,d2)/D;
    #undef Det
    uok = ( tu>=0 && tu<=1 );
    vok = ( tv>=0 && tv<=1 );
    #if 0
    if (uok && vok) {
      dx = ux[0]*tu + ux[1]*(1-tu);
      dy = uy[0]*tu + uy[1]*(1-tu);
      util (x, y, dx, dy, S, s, z);
    } else if (uok || vok) {
      ABORT("strange cell %d %d",x,y);		/* error in algorithm */
    }
    #else
    if (uok) {
      dx = ux[0]*tu + ux[1]*(1-tu);
      dy = uy[0]*tu + uy[1]*(1-tu);
      util (x, y, dx, dy, S, s, z);
    } else if (vok) {
      dx = ux[0]*tu + ux[1]*(1-tu);
      dy = uy[0]*tu + uy[1]*(1-tu);
      util (x, y, dx, dy, S, s, z);
    }
    #endif
  }
  return 1; 		/* successfull */
}

RUN_HEAD(singz)
{
  int i;
  int z;
#ifdef X11
  DEVICE_CONST(BGIWindow,wnd);
  int graphics=(wnd.row1>wnd.row0) && (wnd.col1>wnd.col0);
#endif
#if MPI
  int imroot=(mpi_rank==S->root);
#else
#define imroot 1
#endif
  
  /* Graphic output: repaint the old tip(s) to the foreground color */
#ifdef X11
  if (graphics) {
    SetWindow(wnd);
    if NOT(SetLimits(s.x0, s.x1, s.y0, s.y1)) return 0;
    for (i=0;i<S->ntip;i++)
      Pixel (S->xtip[i],S->ytip[i],wnd.color%16);
  }
#endif
  S->ntip=0;

  /* File output while finding new tip(s) */
  if (imroot && S->everyrecord) output(S,"%10ld ",t);
  for (z=s.z0;z<=s.z1;z++) {
    if (imroot && S->everysection) output(S,"%10ld ",t);
    isocross(S,&s,z);
    if (imroot) output(S,"%s",S->sectionsep);
  }
  if (imroot) output(S,"%s",S->recordsep);
  flush_buffer(S);

  /* Graphic output: paint the new tip(s) to the border color */
#ifdef X11
  if (graphics) {
    for (i=0;i<S->ntip;i++) 
      Pixel (S->xtip[i],S->ytip[i],wnd.color/16);
  }
#endif

  if (S->k_output) {
    int N=S->ntip;
    if (N>0) {
      double X=0; double Y=0; double Z=0;
      double X2=0; double Y2=0; double Z2=0;
      for (i=0;i<N;i++) {
	X+=S->xtip[i]; X2+=S->xtip[i]*S->xtip[i];
	Y+=S->ytip[i]; Y2+=S->ytip[i]*S->ytip[i];
	Z+=S->ztip[i]; Z2+=S->ztip[i]*S->ztip[i];
      }
      X/=N;Y/=N;Z/=N;
      X2/=N;Y2/=N;Z2/=N;
      if (S->Npoints) *(S->Npoints)=N;
      if (S->xmean)   *(S->xmean)  =X;
      if (S->ymean)   *(S->ymean)  =Y;
      if (S->zmean)   *(S->zmean)  =Z;
      if (S->xstddev) *(S->xstddev)=sqrt(X2-X*X);
      if (S->ystddev) *(S->ystddev)=sqrt(Y2-Y*Y);
      if (S->zstddev) *(S->zstddev)=sqrt(Z2-Z*Z);
    } else {
      if (S->Npoints) *(S->Npoints)=N;
    }
  }
#undef s
}
RUN_TAIL(singz)

DESTROY_HEAD(singz)
  SAFE_CLOSE(S->file);
DESTROY_TAIL(singz)

/********************************************/
CREATE_HEAD(singz)
{
  DEVICE_REQUIRES_SYNC;
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 != v0 );
  ACCEPTR(c0,RNONE,RNONE,RNONE);
  ACCEPTR(c1,RNONE,RNONE,RNONE);
  S->K[0] = v0;  S->K[1] = v1;
  ACCEPTI(precise_orientation,0,0,1);
#if MPI
  if (S->precise_orientation) {
    MESSAGE("\nprecise_orientation=1 is not allowed in MPI. "
	    "Your simulation will continue with precise_orientation=0. */");
    S->precise_orientation=0;
  }
#endif
#if MPI
  /* TODO: make it work for arbitrary user's choice of everys. */
  /* For now, it does not work in any other way.               */
  ACCEPTS(pointsep,"\n");
  ACCEPTS(sectionsep,"");
  ACCEPTS(recordsep,"");
  ACCEPTI(everypoint,1,1,1);
  ACCEPTI(everysection,0,0,0);
  ACCEPTI(everyrecord,0,0,0);
#else
  ACCEPTS(pointsep," ");
  ACCEPTS(sectionsep,"\t");
  ACCEPTS(recordsep,"\n");
  ACCEPTI(everypoint,0,0,1);
  ACCEPTI(everysection,0,0,1);
  ACCEPTI(everyrecord,1,0,1);
#endif
  S->ntip = 0;
  {
    int everypoint=S->everypoint;
    int everysection=S->everysection;
    int everyrecord=S->everyrecord;
    ASSERT(everypoint+everysection+everyrecord==1);
  }
#define ACCEPTIFGIVEN(a) ACCEPTKR(a,1,"")
  ACCEPTIFGIVEN(Npoints);
  ACCEPTIFGIVEN(xmean);
  ACCEPTIFGIVEN(ymean);
  ACCEPTIFGIVEN(zmean);
  ACCEPTIFGIVEN(xstddev);
  ACCEPTIFGIVEN(ystddev);
  ACCEPTIFGIVEN(zstddev);
#undef ACCEPTIFGIVEN
  S->k_output=
    (S->Npoints!=NULL) ||
    (S->xmean!=NULL) ||
    (S->ymean!=NULL) ||
    (S->zmean!=NULL) ||
    (S->xstddev!=NULL) ||
    (S->ystddev!=NULL) ||
    (S->zstddev!=NULL);
  
  
  S->p = &(S->buf[0]);
  
  ACCEPTI(append,1,0,1);
#if MPI
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    EXPECTED_ERROR("Could not create communicator.\n");
  S->tag=idev;
  /* the choice below does not work if global_x0... point falls in the void                */
  /* S->root = getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0); */
  if (mpi_rank==S->root) {
    ACCEPTF(file,S->append?"at":"wt","");
    S->fileout=(S->file!=NULL);
  } else {
    S->fileout=accepts("file=",&(S->filename[0]), "", w) && (S->filename[0]!='\0');
    S->file=NULL;
  }
  /* TODO: make this listpointable? */
  /* Need to cover the cleft space */
  if (mpi_ix == mpi_nx-1) dev->s.x1--;
  if (mpi_iy == mpi_ny-1) dev->s.y1--;
#else
  dev->s.x1--;
  dev->s.y1--;
  ACCEPTF(file,S->append?"at":"wt","");
  S->fileout=(S->file!=NULL);
#endif

#ifdef X11
  ACCEPT_WINDOW(wnd);
#endif
}
CREATE_TAIL(singz,1)
