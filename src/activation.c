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

/* This device calculates activation  and recovery times, as defined in */
/* defined as positive and negative crossings of the field in layer v0,	*/
/* through the threshold 'threshold',					*/
/* the time is measured by the value of main loop counter t,		*/
/* the result is put in layer v1 and output to file 'file' if given.	*/
/* An extra layer vd stores the previous values of the field.		*/

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
#include "bikt.h"
#include "qpp.h"

#if MPI
extern int idev;
#endif

typedef struct {
  int begun;			/* already started: previous is defined */
  real threshold;		/* activation/recovery threshold */
  int sign;			/* +1/-1/0 = up/down/both	 */
  int v0;			/* layer for the source field */
  int v1;			/* layer for the answer */
  int vd;			/* layer for previous values 	 */
  /* hereon calc from singz */
  #define SINGZBUFMAX (4*1024*1024-1) /* 4M-1, reserve one byte in the end for trailing zero */
  char buf[SINGZBUFMAX+1];	/* output buffer: should be long enough to print whole record incl trailing zero */
  char *p;			/* pointer showing after the end of the record in the buffer so far */
  int append;			/* 1 if append to the file	 */
  FILE *file;			/* file descriptor 		 */
  char filename[MAXPATH];	/* name of the output file 	 */
  int fileout;			/* file output flag; may be non-zero even in non-root procs */
  #if MPI
  MPI_Comm comm;		/* communicator for this device instance */
  int tag;			/* tag for messages of this device instance (don't need if private communicator?) */
  int root;			/* the rank of the instance that does the file output */
  #endif
} STR;

/*---------------------------------------------------------------*/
/* Output to text file through the root process: calc from singz */
/* Add the message to the buffer of this process */
static int overflow_reported=0;
static void output(STR *S,char *fmt, ...) {
  va_list ap;
  int l, n, nmax;
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
  /* S->p=&(S->buf[0])+snprintf(S->buf,SINGZBUFMAX+1,"\n"); */	/* In any case reset my buffer. DO WE REALLY NEED EXTRA EMPTY LINES THERE??*/
  S->p=&(S->buf[0]);				/* In any case reset my buffer. */
  *(S->p)='\0';
#else
  fputs(S->buf,S->file);			/* Sequential: just output the buffer */
  fflush(S->file);				/* always flush the file after writing to it */
  S->p=&(S->buf[0]);				/* and reset the buffer. */	
  *(S->p)='\0';
#endif
}
/*---------------------------------------------------------------*/


/********************/
RUN_HEAD(activation)
{
  DEVICE_CONST(real,threshold)
  DEVICE_CONST(int,sign)
  DEVICE_CONST(int,v0)
  DEVICE_CONST(int,v1)
  DEVICE_CONST(int,vd)
  DEVICE_CONST(int,fileout)
  int x, y, z;
  real y0, y1, p, q, thetime;
  real *u;
  int prev=DV*(vd-v0);
  int res=DV*(v1-v0);
  for(z=s.z0;z<=s.z1;z++) {
    for(y=s.y0;y<=s.y1;y++) {
      for(x=s.x0;x<=s.x1;x++) {
	u=New+ind(x,y,z,v0);
	y0=u[prev]-threshold;
	y1=u[0]-threshold;
	if( ((sign>=0)&&(y0<=0)&&(y1>0))	/* crossing in positive direction expected and happened */
	    || ((sign<=0)&&(y0>=0)&&(y1<0))) { 	/* crossing in negative direction expected and happened */
	  p=y1/(y1-y0); q=-y0/(y1-y0); 		/* linear interpolation coefficients */
	  assert(p>=0); assert(p<=1);
	  assert(q>=0); assert(q<=1);
	  thetime=p*(t-1)+q*t;  		/* interpolate the t variable */
	  u[res]=thetime;			/* remember it */
	  if (fileout&&S->begun) 		/* and output it */
	    output(S,"%3d %3d %3d %7.3f %2d\n",x,y,z,thetime,((y1>y0)?1:(-1)));
	}
	u[prev]=u[0];
      }
    }
  }
  flush_buffer(S);
  S->begun=1;
}
RUN_TAIL(activation)

/**********************/
DESTROY_HEAD(activation)
DESTROY_TAIL(activation)

/*********************/
CREATE_HEAD(activation)
{
  ACCEPTR(threshold,RNONE,RNONE,RNONE);
  ACCEPTI(sign,INONE,-1,1);
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ACCEPTI(vd,INONE,0,(int)vmax-1);
  ASSERT( v1 != v0 );
  ASSERT( vd != v0 );
  ASSERT( vd != v1 );
  S->begun=0;

  /*-----------------*/
  /* calc from singz */
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
  /* Need to cover the cleft space */
  if (mpi_ix < mpi_nx-1) dev->s.x1++;
  if (mpi_iy < mpi_ny-1) dev->s.y1++;
  if (mpi_iz < mpi_nz-1) dev->s.z1++;
#else
  ACCEPTF(file,S->append?"at":"wt","");
  S->fileout=(S->file!=NULL);
#endif
  /*-----------------*/
}
CREATE_TAIL(activation,1)

