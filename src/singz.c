/**
 * Copyright (C) (2010-2026) Vadim Biktashev, Irina Biktasheva et al.
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
/*									           */
/* File output: each process prints its tips into its own text buffer, every	   */
/* tip tagged by the global position of its cell. The root collects the buffers,  */
/* sorts the tips and assembles the record, with all time stamps and separators,  */
/* so the file does not depend on the domain decomposition and is the same as in  */
/* the sequential mode (see textbuf.h).					   */
/* Tip statistics k-variables are reduced over all processes.		 	   */

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
#include "textbuf.h"
#ifdef X11
#include "windraw.h"
#endif

/* runHere is only set in MPI, despite the comment in device.h */
#if MPI
#define RUNS_HERE(s) ((s).runHere)
#else
#define RUNS_HERE(s) 1
#endif

/* Tip statistics: number of tips, sums of coordinates and of their squares */
enum {SUM_N, SUM_X, SUM_Y, SUM_Z, SUM_X2, SUM_Y2, SUM_Z2, NSUM};

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
  textpieces mine;		/* texts of this process's tips in the current record */
  #if MPI
  textpieces all;		/* ..., of all processes (root only) */
  #endif
  textbuf out;			/* the record as it goes to the file (root only) */
  #if MPI
  MPI_Comm comm;		/* communicator for this device instance */
  int root;			/* rank of the reporting process in ALL_ACTIVE_PROCS */
  int commroot;			/* ..., and in comm */
  #endif
  #define MAXSTACK 10
  real xtip[MAXSTACK],ytip[MAXSTACK],ztip[MAXSTACK]; /* tips found at previous/next step, for graphics */
  int ntip;
  double sum[NSUM];		/* statistics of all tips found by this process */

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

/* Sort the tips of the record and assemble it in S->out, with time stamps and  */
/* separators. Sequentially, the tips are those of this process; in MPI, of all. */
static void assemble_record(STR *S, textpieces *rec, long z0, long z1) {
  size_t i=0;
  long z;
  tp_sort(rec);
  S->out.n=0;
  if (S->everyrecord) txt_printf(&(S->out),"%10ld ",t);
  for (z=z0;z<=z1;z++) {
    if (S->everysection) txt_printf(&(S->out),"%10ld ",t);
    for (;i<rec->n && tp_key_z(rec->pc[i].key)==z;i++)
      txt_write(&(S->out),rec->txt.s+rec->pc[i].off,(size_t)rec->pc[i].len);
    txt_puts(&(S->out),S->sectionsep);
  }
  txt_puts(&(S->out),S->recordsep);
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
    /* One tip per cell, so the cell defines the order of output */
    tp_begin(&(S->mine),tp_key(x,y,z));
    if (S->everypoint) txt_printf(&(S->mine.txt),"%10ld",t);
    /* 2D or 3D is decided by the whole device space, not this process's part */
    if (s->global_z0==s->global_z1)
      txt_printf(&(S->mine.txt)," %f %f %f %f %f",
		 (float)(x+dx),(float)(y+dy),(float)phi,(float)psi,(float)cosangle);
    else
      txt_printf(&(S->mine.txt)," %f %f %f %f %f %f",
		 (float)(x+dx),(float)(y+dy),(float)(z),(float)phi,(float)psi,(float)cosangle);
    txt_puts(&(S->mine.txt),S->pointsep);
    tp_end(&(S->mine));
  }
  /* Statistics count all tips; only the first MAXSTACK are kept for graphics */
  S->sum[SUM_N]  += 1;
  S->sum[SUM_X]  += x+dx;  S->sum[SUM_X2] += (x+dx)*(x+dx);
  S->sum[SUM_Y]  += y+dy;  S->sum[SUM_Y2] += (y+dy)*(y+dy);
  S->sum[SUM_Z]  += z;     S->sum[SUM_Z2] += (double)z*z;
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
    if (uok) {
      dx = ux[0]*tu + ux[1]*(1-tu);
      dy = uy[0]*tu + uy[1]*(1-tu);
      util (x, y, dx, dy, S, s, z);
    } else if (vok) {
      dx = ux[0]*tu + ux[1]*(1-tu);
      dy = uy[0]*tu + uy[1]*(1-tu);
      util (x, y, dx, dy, S, s, z);
    }
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

  /* In MPI, this may be an instance with nothing to do but the reduction of statistics */
  if (RUNS_HERE(s)) {
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
    for (i=0;i<NSUM;i++) S->sum[i]=0;
    tp_reset(&(S->mine));

    /* Find new tip(s), printing them into the buffer */
    for (z=s.z0;z<=s.z1;z++) isocross(S,&s,z);

    /* File output of the whole record by the root */
    if (S->fileout) {
#if MPI
      if (!tp_gather(&(S->mine),&(S->all),S->commroot,S->comm)) return 0;
      if (imroot) assemble_record(S,&(S->all),s.global_z0,s.global_z1);
#else
      assemble_record(S,&(S->mine),s.global_z0,s.global_z1);
#endif
      if (imroot) {
	fwrite(S->out.s,1,S->out.n,S->file);
	fflush(S->file);
      }
    }

    /* Graphic output: paint the new tip(s) to the border color */
#ifdef X11
    if (graphics) {
      for (i=0;i<S->ntip;i++)
	Pixel (S->xtip[i],S->ytip[i],wnd.color/16);
    }
#endif
  } /* if s.runHere */

  if (S->k_output) {
    double sum[NSUM];
    double N;
    for (i=0;i<NSUM;i++) sum[i]=RUNS_HERE(s)?S->sum[i]:0;
#if MPI
    /* All processes need the k-variables, hence this device always runs */
    MPIDO(MPI_Allreduce(MPI_IN_PLACE,sum,NSUM,MPI_DOUBLE,MPI_SUM,ALL_ACTIVE_PROCS),
	  "singz: could not reduce tip statistics");
#endif
    N=sum[SUM_N];
    if (S->Npoints) *(S->Npoints)=N;
    if (N>0) {
      double X=sum[SUM_X]/N,   Y=sum[SUM_Y]/N,   Z=sum[SUM_Z]/N;
      double X2=sum[SUM_X2]/N, Y2=sum[SUM_Y2]/N, Z2=sum[SUM_Z2]/N;
      if (S->xmean)   *(S->xmean)  =X;
      if (S->ymean)   *(S->ymean)  =Y;
      if (S->zmean)   *(S->zmean)  =Z;
      if (S->xstddev) *(S->xstddev)=sqrt(fmax(X2-X*X,0));
      if (S->ystddev) *(S->ystddev)=sqrt(fmax(Y2-Y*Y,0));
      if (S->zstddev) *(S->zstddev)=sqrt(fmax(Z2-Z*Z,0));
    }
  }
#undef imroot
}
RUN_TAIL(singz)

DESTROY_HEAD(singz)
  tp_free(&(S->mine));
#if MPI
  tp_free(&(S->all));
#endif
  txt_free(&(S->out));
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
  /* It needs two layers of halo, and we have only one */
  if (S->precise_orientation) {
    MESSAGE("\nprecise_orientation=1 is not allowed in MPI. "
	    "Your simulation will continue with precise_orientation=0. */");
    S->precise_orientation=0;
  }
#endif
  /* Since the root assembles the records, these work the same in MPI */
  ACCEPTS(pointsep," ");
  ACCEPTS(sectionsep,"\t");
  ACCEPTS(recordsep,"\n");
  ACCEPTI(everypoint,0,0,1);
  ACCEPTI(everysection,0,0,1);
  ACCEPTI(everyrecord,1,0,1);
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
  /* The statistics are global k-variables, so all processes take part in computing them */
  if (S->k_output) {DEVICE_ALWAYS_RUNS}

  ACCEPTI(append,1,0,1);
#if MPI
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    EXPECTED_ERROR("Could not create communicator.\n");
  if (dev->s.runHere) {
    /* Where is the root in the device's communicator */
    int mine=-1, me;
    MPI_Comm_rank(S->comm,&me);
    if (mpi_rank==S->root) mine=me;
    MPIDO(MPI_Allreduce(&mine,&(S->commroot),1,MPI_INT,MPI_MAX,S->comm),
	  "singz: could not identify the root");
  }
  /* the choice below does not work if global_x0... point falls in the void                */
  /* S->root = getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0); */
  if (mpi_rank==S->root) {
    ACCEPTF(file,S->append?"at":"wt","");
  } else {
    if (!accepts("file=",&(S->filename[0]),"",w)) return 0;
    S->file=NULL;
  }
#else
  ACCEPTF(file,S->append?"at":"wt","");
#endif
  S->fileout=(S->filename[0]!='\0');

  /* TODO: make this listpointable? */
  /* Need to cover the cleft space: cells are [x,x+1]*[y,y+1], so the last column  */
  /* and row of the device space are not cells' origins. In MPI, only the process  */
  /* owning that column/row drops it; others use their halo for x+1, y+1.	   */
  if (dev->s.x1==dev->s.global_x1) dev->s.x1--;
  if (dev->s.y1==dev->s.global_y1) dev->s.y1--;

#ifdef X11
  ACCEPT_WINDOW(wnd);
#endif
}
CREATE_TAIL(singz,1)
