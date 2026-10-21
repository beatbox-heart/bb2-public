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

/* Modify the on-grid diff weights for an Ohmic inner boundary, */
/* with the correction by boundary slope. */
/* Works in sequential and MPI modes; void points are neutral in the */
/* slope estimate. Slope correction is implemented only in 2D. */

#if MPI
#include <mpi.h>
#endif

#include <math.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "qpp.h"

/* Constants in stencils named as in popular formulas */
static const real sixth      = 1./6.;
static const real foursixths = 4./6.;

typedef struct {
  int advance;				/* this instance runs first when created, before time starts */
  int v0;				/* indicator layers */
  real c0, c1;				/* indicator constants */
  PAR(real,Sigma);			/* Conductivity across the boundary (SI: m/s) */
  real hx;				/* Space step (SI: m) */
  int manypoint; 			/* == NINEPOINT, NINETEENPT in ezdiff < ezspiral, ezscroll */ 
  int w0, w1;				/* Layers range to keep the stencil weights */
  int check_positivity;			/* Report points with positive central weight */
  int enforce_positivity;		/* Remove their negative links to make them non-positive */
} STR;

enum i1_nbrs {i10, i1p, i1m, i1nnb};
enum i2_nbrs {i200, i2p0, i2m0, i20p, i20m, i2nnb};
enum m2_nbrs {m200, m2p0, m2m0, m20p, m20m,
	      m2pp, m2pm, m2mp, m2mm, m2nnb};
enum i3_nbrs {i3000, i3p00, i3m00, i30p0, i30m0, i300p, i300m, i3nnb};
enum m3_nbrs {m3000, m3p00, m3m00, m30p0, m30m0, m300p, m300m,
	      m3pp0, m3pm0, m3mp0, m3mm0, m3p0p, m3p0m,
	      m3m0p, m3m0m, m30pp, m30pm, m30mp, m30mm, m3nnb};

/* Neighbours' displacements, in the order of the enums above */
static const int off_i1[i1nnb][3]={{0,0,0},{+1,0,0},{-1,0,0}};
static const int off_i2[i2nnb][3]={{0,0,0},{+1,0,0},{-1,0,0},{0,+1,0},{0,-1,0}};
static const int off_m2[m2nnb][3]={{0,0,0},{+1,0,0},{-1,0,0},{0,+1,0},{0,-1,0},
				   {+1,+1,0},{+1,-1,0},{-1,+1,0},{-1,-1,0}};

/* Report a point with positive central weight, with its stencil laid out as on the grid */
static void report_bad(INT x, INT y, INT z, int w0, int manypoint)
{
#define W(qq) ((double)New[ind(x,y,z,w0+qq)])
  URGENT_MESSAGE("\n/* #%d t=%ld inner_bc: bad point (%ld,%ld,%ld) W=\n",
		 mpi_rank,(long)t,(long)x,(long)y,(long)z);
  switch (dim) {
  case 1:
    URGENT_MESSAGE("%g\t%g\t%g\n",W(i1m),W(i10),W(i1p));
    break;
  case 2:
    if (manypoint)
      URGENT_MESSAGE("%g\t%g\t%g\n%g\t%g\t%g\n%g\t%g\t%g\n",
		     W(m2mp),W(m20p),W(m2pp),
		     W(m2m0),W(m200),W(m2p0),
		     W(m2mm),W(m20m),W(m2pm));
    else
      URGENT_MESSAGE("%g\t%g\t%g\n%g\t%g\t%g\n%g\t%g\t%g\n",
		     0.0,W(i20p),0.0,
		     W(i2m0),W(i200),W(i2p0),
		     0.0,W(i20m),0.0);
    break;
  }
  URGENT_MESSAGE(" The central weight is positive: %g, which will cause numerical instability. */\n",W(0));
#undef W
}

static unsigned char indicator_class(real value, real c0, real c1)
{
  if (value==c0) return 1;
  if (value==c1) return 2;
  return 0;
}

/* Void points are neutral in the slope estimate, whatever their indicator */
static unsigned char point_class(INT x, INT y, INT z, int v0, real c0, real c1)
{
  if (!istissue(x,y,z)) return 0;
  return indicator_class(New[ind(x,y,z,v0)],c0,c1);
}

#if MPI
/* kslope reaches two cells beyond a local subdomain, including diagonals. */
static size_t indicator_index(INT x, INT y, INT xmin, INT ymin, size_t ystride)
{
  return (size_t)(x-(xmin-2))*ystride + (size_t)(y-(ymin-2));
}

/* Absolute boundary points are stored, as halo, by the process owning
   the adjacent interior points, so clamp them into the interior first. */
static INT clamp_to_interior(INT q, INT qmin, INT qmax)
{
  if (q<qmin) return qmin;
  if (q>qmax-1) return qmax-1;
  return q;
}

static int indicator_owner(INT x, INT y, INT z)
{
  long ix, iy, iz;
  ix=iy=iz=-1;
  if (x<0 || y<0 || z<0) return NORANK;
  if ((size_t)x>=xmax || (size_t)y>=ymax || (size_t)z>=zmax) return NORANK;
  x=clamp_to_interior(x,ONE,(INT)xmax-ONE);
  y=clamp_to_interior(y,TWO,(INT)ymax-TWO);
  z=clamp_to_interior(z,TRI,(INT)zmax-TRI);
  for (long i=0;i<mpi_nx;i++) {
    if (all_xmin[i]<=(size_t)x && (size_t)x<all_xmax[i]) {ix=i; break;}
  }
  for (long i=0;i<mpi_ny;i++) {
    if (all_ymin[i]<=(size_t)y && (size_t)y<all_ymax[i]) {iy=i; break;}
  }
  for (long i=0;i<mpi_nz;i++) {
    if (all_zmin[i]<=(size_t)z && (size_t)z<all_zmax[i]) {iz=i; break;}
  }
  if (ix<0 || iy<0 || iz<0) return NORANK;
  return ProcessRank[sind(ix,iy,iz)];
}

static void prepare_indicator_halo(int v0, real c0, real c1,
                                   unsigned char **indicator,
                                   size_t *ystride)
{
  int nranks, total_send=0, total_recv=0;
  int *send_counts, *recv_counts, *send_displs, *recv_displs;
  int *coord_send_counts, *coord_recv_counts;
  int *coord_send_displs, *coord_recv_displs, *cursor;
  long *coords_send, *coords_recv;
  size_t *destinations;
  unsigned char *values_send, *values_recv;
  unsigned char *field;
  size_t nx=(size_t)(local_xmax-local_xmin);
  size_t ny=(size_t)(local_ymax-local_ymin);
  size_t field_x, field_y, field_size;
  INT xmin=(INT)local_xmin, ymin=(INT)local_ymin, z=(INT)local_zmin;

  if (nx>SIZE_MAX-4 || ny>SIZE_MAX-4)
    ABORT("inner_bc indicator halo dimensions overflow\n");
  field_x=nx+4;
  field_y=ny+4;
  if (field_x && field_y>SIZE_MAX/field_x)
    ABORT("inner_bc indicator halo size overflows\n");
  field_size=field_x*field_y;
  field=calloc(field_size,sizeof(*field));
  if (!field) ABORT("Could not allocate inner_bc indicator halo\n");
  *ystride=field_y;
  *indicator=field;

  for (INT x=xmin;x<(INT)local_xmax;x++) {
    for (INT y=ymin;y<(INT)local_ymax;y++) {
      field[indicator_index(x,y,xmin,ymin,field_y)] =
        point_class(x,y,z,v0,c0,c1);
    }
  }

  MPIDO(MPI_Comm_size(ALL_ACTIVE_PROCS,&nranks),
        "Couldn't get active process count for inner_bc.");
  send_counts=calloc(nranks,sizeof(*send_counts));
  recv_counts=calloc(nranks,sizeof(*recv_counts));
  send_displs=calloc(nranks,sizeof(*send_displs));
  recv_displs=calloc(nranks,sizeof(*recv_displs));
  coord_send_counts=calloc(nranks,sizeof(*coord_send_counts));
  coord_recv_counts=calloc(nranks,sizeof(*coord_recv_counts));
  coord_send_displs=calloc(nranks,sizeof(*coord_send_displs));
  coord_recv_displs=calloc(nranks,sizeof(*coord_recv_displs));
  cursor=calloc(nranks,sizeof(*cursor));
  if (!send_counts || !recv_counts || !send_displs || !recv_displs ||
      !coord_send_counts || !coord_recv_counts ||
      !coord_send_displs || !coord_recv_displs || !cursor)
    ABORT("Could not allocate inner_bc halo exchange metadata\n");

  for (INT x=xmin-2;x<(INT)local_xmax+2;x++) {
    for (INT y=ymin-2;y<(INT)local_ymax+2;y++) {
      if (x>=xmin && x<(INT)local_xmax &&
          y>=ymin && y<(INT)local_ymax) continue;
      int owner=indicator_owner(x,y,z);
      if (owner==mpi_rank) { /* absolute boundary point in own halo */
        field[indicator_index(x,y,xmin,ymin,field_y)] =
          point_class(x,y,z,v0,c0,c1);
      } else if (owner>=0) {
        if (send_counts[owner]==INT_MAX)
          ABORT("Too many inner_bc halo requests for one rank\n");
        send_counts[owner]++;
      }
    }
  }

  MPIDO(MPI_Alltoall(send_counts,1,MPI_INT,recv_counts,1,MPI_INT,
                     ALL_ACTIVE_PROCS),
        "Couldn't exchange inner_bc halo request counts.");
  for (int rank=0;rank<nranks;rank++) {
    if (send_counts[rank]>INT_MAX/2 || recv_counts[rank]>INT_MAX/2)
      ABORT("Too many inner_bc halo requests for MPI_Alltoallv\n");
    if (total_send>INT_MAX-send_counts[rank] ||
        total_recv>INT_MAX-recv_counts[rank])
      ABORT("inner_bc halo request count exceeds MPI limits\n");
    send_displs[rank]=total_send;
    recv_displs[rank]=total_recv;
    coord_send_displs[rank]=2*total_send;
    coord_recv_displs[rank]=2*total_recv;
    total_send+=send_counts[rank];
    total_recv+=recv_counts[rank];
    coord_send_counts[rank]=2*send_counts[rank];
    coord_recv_counts[rank]=2*recv_counts[rank];
  }

  coords_send=calloc(total_send?2*(size_t)total_send:1,sizeof(*coords_send));
  coords_recv=calloc(total_recv?2*(size_t)total_recv:1,sizeof(*coords_recv));
  destinations=calloc(total_send?total_send:1,sizeof(*destinations));
  values_send=calloc(total_recv?total_recv:1,sizeof(*values_send));
  values_recv=calloc(total_send?total_send:1,sizeof(*values_recv));
  if (!coords_send || !coords_recv || !destinations ||
      !values_send || !values_recv)
    ABORT("Could not allocate inner_bc halo exchange buffers\n");
  memcpy(cursor,send_displs,(size_t)nranks*sizeof(*cursor));

  for (INT x=xmin-2;x<(INT)local_xmax+2;x++) {
    for (INT y=ymin-2;y<(INT)local_ymax+2;y++) {
      if (x>=xmin && x<(INT)local_xmax &&
          y>=ymin && y<(INT)local_ymax) continue;
      int owner=indicator_owner(x,y,z);
      if (owner>=0 && owner!=mpi_rank) {
        int pos=cursor[owner]++;
        coords_send[2*(size_t)pos]=x;
        coords_send[2*(size_t)pos+1]=y;
        destinations[pos]=indicator_index(x,y,xmin,ymin,field_y);
      }
    }
  }

  MPIDO(MPI_Alltoallv(coords_send,coord_send_counts,coord_send_displs,MPI_LONG,
                      coords_recv,coord_recv_counts,coord_recv_displs,MPI_LONG,
                      ALL_ACTIVE_PROCS),
        "Couldn't exchange inner_bc halo coordinates.");
  for (int i=0;i<total_recv;i++) {
    INT x=coords_recv[2*(size_t)i];
    INT y=coords_recv[2*(size_t)i+1];
    if (indicator_owner(x,y,z)!=mpi_rank)
      ABORT("Invalid remote coordinate in inner_bc halo request: (%ld,%ld)\n",
            (long)x,(long)y);
    values_send[i]=point_class(x,y,z,v0,c0,c1);
  }
  MPIDO(MPI_Alltoallv(values_send,recv_counts,recv_displs,MPI_UNSIGNED_CHAR,
                      values_recv,send_counts,send_displs,MPI_UNSIGNED_CHAR,
                      ALL_ACTIVE_PROCS),
        "Couldn't exchange inner_bc halo values.");
  for (int i=0;i<total_send;i++) field[destinations[i]]=values_recv[i];

  free(values_recv);
  free(values_send);
  free(destinations);
  free(coords_recv);
  free(coords_send);
  free(cursor);
  free(coord_recv_displs);
  free(coord_send_displs);
  free(coord_recv_counts);
  free(coord_send_counts);
  free(recv_displs);
  free(send_displs);
  free(recv_counts);
  free(send_counts);
}
#endif

static unsigned char kslope_sample(INT x, INT y, INT z, INT v0,
                                   real ca, real cb,
                                   const unsigned char *indicator,
                                   size_t ystride, INT xmin, INT ymin)
{
#if MPI
  if (indicator)
    return indicator[indicator_index(x,y,xmin,ymin,ystride)];
#else
  (void)indicator;
  (void)ystride;
  (void)xmin;
  (void)ymin;
#endif
  return point_class(x,y,z,v0,ca,cb);
}

/* Currently only 1D and 2D implemented so z1=z2 always */
static real kslope (INT x1, INT y1, INT z1, INT x2, INT y2, INT z2,
                    INT v0, real ca, real cb, int manypoint,
                    const unsigned char *indicator, size_t ystride,
                    INT xmin, INT ymin)
{
  unsigned char c1, c2, thisv;
  unsigned char v1=kslope_sample(x1,y1,z1,v0,ca,cb,indicator,ystride,xmin,ymin);
  unsigned char v2=kslope_sample(x2,y2,z2,v0,ca,cb,indicator,ystride,xmin,ymin);
  real Q;
  real Sx, Sy, S, dx, dy;
  long x, y, X, Y, N;
  real alpha, p, q, k, cosalpha, sinalpha;
  
  ASSERT(dim==2);	/* do 3D later */
  ASSERT(z1==z2);
  ASSERT((x1!=x2) || (y1!=y2));
  
  if (v1==1 || v1==2) c1=v1;
  else ABORT("New[(x1=%ld,y1=%ld,z1=%ld,v0=%ld)] is neither %g nor %g\n",
	     (long)x1,(long)y1,(long)z1,(long)v0,ca,cb);
  
  if (v2==1 || v2==2) c2=v2;
  else ABORT("New[(x2=%ld,y2=%ld,z2=%ld,v0=%ld)] is neither %g nor %g\n",
	     (long)x2,(long)y2,(long)z2,(long)v0,ca,cb);

  S=Sx=Sy=0.0;
  N=X=Y=0;
  
#define account(dx,dy)				\
  x=x1+dx; y=y1+dy;				\
  thisv=kslope_sample(x,y,z1,v0,ca,cb,indicator,ystride,xmin,ymin); \
  Q=(thisv==c1)?1:(thisv==c2)?(-1):0;		\
  N++; X+=x; Y+=y; S+=Q; Sx+=Q*x; Sy+=Q*y;
  account(+1,0); 
  account(-1,0); 
  account(0,+1); 
  account(0,-1); 
  account(+1,+1); 
  account(+1,-1); 
  account(-1,+1); 
  account(-1,-1); 
#undef account
#define account(dx,dy)				\
  x=x2+dx; y=y2+dy;				\
  thisv=kslope_sample(x,y,z2,v0,ca,cb,indicator,ystride,xmin,ymin); \
  Q=(thisv==c1)?1:(thisv==c2)?(-1):0;		\
  N++; X+=x; Y+=y; S+=Q; Sx+=Q*x; Sy+=Q*y;
  account(+1,0); 
  account(-1,0); 
  account(0,+1); 
  account(0,-1); 
  account(+1,+1); 
  account(+1,-1); 
  account(-1,+1); 
  account(-1,-1); 
#undef account

  ASSERT(N>0);
  dx=fabs(Sx-X*S/N);
  dy=fabs(Sy-Y*S/N);
  /* ASSERT(dx!=0 || dy!=0); */
  if (dx==0 && dy==0) { /* fallback to the normal as the connecting line */
    MESSAGE("/* (x1,y1)=(%ld,%ld), (x2,y2): (%ld,%ld): zero dipole; fallback to straight line */\n",
	    (long)x1,(long)y1,(long)x2,(long)y2);
    dx=1.0*labs(x2-x1);
    dy=1.0*labs(y2-y1);
    alpha=atan2(dy,dx);
  } else {
    alpha=atan2(dy,dx);
  }

  /* alpha=atan2(dy,dx); */
  if (manypoint) {
    p=foursixths;
    q=sixth;
  } else {
    p=1.0;
    q=0.0;
  }
  cosalpha=cos(alpha);
  sinalpha=sin(alpha);
  k=1.0/(p*(cosalpha+sinalpha)+2*q*fmax(cosalpha,sinalpha));
  Debug("kslope:\t%d\t%d\t%d\t%d\t%g\t%g\t%g\t%g\n",x1,y1,x2,y2,dx,dy,alpha,k);
  return k;
}

/****************/
RUN_HEAD(inner_bc)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(real,c0);
  DEVICE_CONST(real,c1);
  DEVICE_CONST(int,w0);
  DEVICE_CONST(int,w1);
  DEVICE_CONST(real,hx);
  DEVICE_PAR(real,Sigma);
  DEVICE_CONST(int,manypoint);
  DEVICE_CONST(int,check_positivity);
  DEVICE_CONST(int,enforce_positivity);

#if MPI
  unsigned char *indicator=NULL;
  size_t indicator_ystride=0;
  prepare_indicator_halo(v0,c0,c1,&indicator,&indicator_ystride);
#else
  unsigned char *indicator=NULL;
  size_t indicator_ystride=0;
#endif

  /* The halo exchange is collective; only ranks with device work update weights. */
#if MPI
  if (s.runHere) {
#else
  { /* s.runHere is not set in sequential mode */
#endif
  int nn=0;
  /* Currently only 1D and 2D implemented so dz=0 always */
#define adjust(dx,dy,dz,qq,wgt)				    \
  if (isTissue(*x+dx,*y+dy,*z+dz)) {                        \
    /* not New: its halo corners may be stale next to empty subdomains */ \
    unsigned char cls0=kslope_sample(*x,*y,*z,v0,c0,c1,     \
      indicator,indicator_ystride,local_xmin,local_ymin);   \
    unsigned char cls1=kslope_sample(*x+dx,*y+dy,*z+dz,v0,c0,c1, \
      indicator,indicator_ystride,local_xmin,local_ymin);   \
    if ((cls0==1 && cls1==2) || (cls0==2 && cls1==1)) {     \
      real *pwgt=New+ind(*x,*y,*z,w0+qq);		    \
      real k=kslope(*x,*y,*z,*x+dx,*y+dy,*z+dz, v0, c0, c1, manypoint, \
                   indicator,indicator_ystride,local_xmin,local_ymin); \
      real awgt=wgt*k;					    \
      New[ind(*x,*y,*z,w0+oo)]+=(*pwgt-awgt);		    \
      *pwgt=awgt;				    	    \
      nn++;						    \
    }							    \
  }

#define w(qq) (real)(New[ind(*x,*y,*z,w0+qq)])

  switch(dim) {
  case 1: {
      real wgt=Sigma/hx;
#define oo i10
#define po i1p /* won't use p to avoid clash with already defined one */
#define mo i1m
#define COMMANDS \
      nn=0; \
      adjust(+1,0,0,po,wgt); \
      adjust(-1,0,0,mo,wgt); \
      if(nn) Debug("%s.%d: (%ld,%ld,%ld): %d\n",__FILE__,__LINE__,(long)(*x),(long)(*y),(long)(*z),nn);
      DO_FOR_ALL_POINTS;
#undef COMMANDS
#undef oo
#undef po
#undef mo
  } break; /* case 1 */
  case 2: {
    if (manypoint==0) {
      real wgt=Sigma/hx;
#define oo i200
#define po i2p0
#define mo i2m0
#define op i20p
#define om i20m
#define COMMANDS \
      nn=0; \
      adjust(+1,0,0,po,wgt); \
      adjust(-1,0,0,mo,wgt); \
      adjust(0,+1,0,op,wgt); \
      adjust(0,-1,0,om,wgt); \
      if (nn) Debug("%s.%d: (%ld,%ld,%ld): %d\n",__FILE__,__LINE__,(long)(*x),(long)(*y),(long)(*z),nn); \
      if (debug && nn>1) fprintf(debug,				\
	      "\t%.6f\t%.6f\t%.6f\n"				\
	      "\t%.6f\t%.6f\t%.6f\n"				\
	      "\t%.6f\t%.6f\t%.6f\n",				\
              0.0,     w(op),0.0,				\
              w(mo),   w(oo), w(po),  				\
              0.0,     w(om), 0.0				\
      );
      DO_FOR_ALL_POINTS;
#undef COMMANDS
#undef oo
#undef po
#undef mo
#undef op
#undef om
    } else {
      real swgt=foursixths*Sigma/hx;
      real dwgt=sixth*Sigma/hx;
#define oo m200
#define po m2p0
#define mo m2m0
#define op m20p
#define om m20m
#define pp m2pp
#define pm m2pm
#define mp m2mp
#define mm m2mm
#define COMMANDS \
      nn=0; \
      adjust(+1,0,0,po,swgt); \
      adjust(-1,0,0,mo,swgt); \
      adjust(0,+1,0,op,swgt); \
      adjust(0,-1,0,om,swgt); \
      adjust(+1,+1,0,pp,dwgt); \
      adjust(+1,-1,0,pm,dwgt); \
      adjust(-1,+1,0,mp,dwgt); \
      adjust(-1,-1,0,mm,dwgt); \
      if(nn) Debug("%s.%d: (%ld,%ld,%ld): %d\n",__FILE__,__LINE__,(long)(*x),(long)(*y),(long)(*z),nn); \
      if (debug && nn>1) fprintf(debug,				\
	      "\t%.6f\t%.6f\t%.6f\n"				\
	      "\t%.6f\t%.6f\t%.6f\n"				\
	      "\t%.6f\t%.6f\t%.6f\n",				\
	      w(mp), w(op), w(pp),				\
	      w(mo), w(oo), w(po),				\
	      w(mm), w(om), w(pm)				\
      );
      DO_FOR_ALL_POINTS;
#undef COMMANDS
#undef oo
#undef po
#undef mo
#undef op
#undef om
#undef pp
#undef pm
#undef mp
#undef mm
#undef adjust
    }
  } break; /* case 2 */
  case 3: {
    /* TBD */
    EXPECTED_ERROR("unexpected dimension %d",dim);
  } break; /* case 3 */
  default:
    EXPECTED_ERROR("unexpected dimension %d",dim);
    break;
  } /* switch dim */
  }

  /* Positivity of the central weights.						*/
  /* An anisotropic stencil has negative (cross-derivative) links, normally	*/
  /* outweighed by the positive links to the orthogonal neighbours. Replacing	*/
  /* the links across the boundary by the weak Ohmic ones may leave a point	*/
  /* with negative links uncompensated, so its central weight becomes positive	*/
  /* and the scheme locally anti-diffusive, i.e. unstable.			*/
  /* Enforcement zeroes all negative links of such points, from both ends:	*/
  /* a link is removed if either of its ends is bad. Each row keeps zero sum,	*/
  /* the bad rows become diagonally dominant, and the operator stays symmetric	*/
  /* if it was. Removing a negative link from a symmetric operator can only	*/
  /* decrease its eigenvalues, so this never makes stability worse.		*/
  if (check_positivity || enforce_positivity) {
    int nnb=0;				/* stencil size */
    const int (*off)[3]=NULL;		/* neighbours' displacements */
    unsigned char *bad=NULL;		/* flags of bad points in local memory, halo included */
    long nbad=0, nlinks=0;		/* counters of bad points and zeroed weights */
    INT xx, yy, zz;
    switch (dim) {
    case 1: nnb=i1nnb; off=off_i1; break;
    case 2: if (manypoint) {nnb=m2nnb; off=off_m2;} else {nnb=i2nnb; off=off_i2;} break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    }
#if MPI
    haloSwap();				/* central weights of neighbours, as adjusted by their ranks */
#endif
#define pind(x,y,z) (ind(x,y,z,0)/vmax)	/* number of the point in local memory */
    /* First pass: flag the bad points before any enforcement, */
    /* so that both ends of a link see the same flags */
    CALLOC(bad,xlen*ylen*zlen,sizeof(unsigned char));
    for (xx=(INT)local_xmin-ONE;xx<(INT)local_xmax+ONE;xx++)
      for (yy=(INT)local_ymin-TWO;yy<(INT)local_ymax+TWO;yy++)
	for (zz=(INT)local_zmin-TRI;zz<(INT)local_zmax+TRI;zz++)
	  bad[pind(xx,yy,zz)]=(isTissue(xx,yy,zz) && New[ind(xx,yy,zz,w0)]>0);
    /* Second pass: report and enforce at the device's points */
#if MPI
    if (s.runHere) {
#else
    {
#endif
#define COMMANDS							\
      if (bad[pind(*x,*y,*z)]) {					\
	nbad++;								\
	if (check_positivity) report_bad(*x,*y,*z,w0,manypoint);	\
      }									\
      if (enforce_positivity) {						\
	int nb;								\
	for (nb=1;nb<nnb;nb++) {					\
	  real *pw=New+ind(*x,*y,*z,w0+nb);				\
	  if (*pw<0 && (bad[pind(*x,*y,*z)] ||				\
			bad[pind(*x+off[nb][0],*y+off[nb][1],*z+off[nb][2])])) { \
	    New[ind(*x,*y,*z,w0)]+=*pw;	/* keep zero row sum */		\
	    *pw=0;							\
	    nlinks++;							\
	  }								\
	}								\
      }
      DO_FOR_ALL_POINTS;
#undef COMMANDS
    }
#undef pind
    FREE(bad);
#if MPI
    {
      long nloc[2]={nbad,nlinks}, nall[2];
      MPIDO(MPI_Allreduce(nloc,nall,2,MPI_LONG,MPI_SUM,ALL_ACTIVE_PROCS),
	    "Couldn't sum up inner_bc positivity counters.");
      nbad=nall[0];
      nlinks=nall[1];
    }
#endif
    if (nbad)
      MESSAGE("\n/* t=%ld inner_bc: %ld points with positive central weight%s */\n",
	      (long)t,nbad,enforce_positivity?" found":"; consider enforce_positivity=1");
    if (nlinks)
      MESSAGE("/* t=%ld inner_bc: %ld negative weights zeroed to enforce positivity */\n",
	      (long)t,nlinks);
  } /* if check_positivity || enforce_positivity */

#if MPI
  free(indicator);
#endif
}
RUN_TAIL(inner_bc)

DESTROY_HEAD(inner_bc)
DESTROY_TAIL(inner_bc)

CREATE_HEAD(inner_bc)
{
  DEVICE_REQUIRES_SYNC;
  DEVICE_ALWAYS_RUNS;
  ACCEPTI(advance,0,0,1);
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTR(c0,RNONE,RNONE,RNONE);
  ACCEPTR(c1,RNONE,RNONE,RNONE);
  ASSERT(c1!=c0);
  ACCEPTR(hx,RNONE,RSUCC(0),RNONE);
  ACCEPTRK(Sigma,RNONE,0,RNONE);  /* Sigma=0 is full disconnect; makes sense, too */
  int nnb;
  ACCEPTI(manypoint,1,0,1); /* this choice is independent on diffusivity of the bulk? */
  switch (dim) {
  case 1:
    nnb=i1nnb;
    break;
  case 2:
    nnb=(manypoint)?m2nnb:i2nnb;
    break;
  case 3:
    nnb=(manypoint)?m3nnb:i3nnb;
    break;
  }
  ACCEPTI(w0,INONE,0,vmax-1);
  ACCEPTI(w1,INONE,w0,vmax-1);  
  if (w1-w0+1!=nnb) {
    MESSAGE("w0=%d w1=%d vmax=%d nnb=%d\n",w0,w1,vmax,nnb);
    ASSERT(w1-w0+1==nnb);
  }
  ACCEPTI(check_positivity,1,0,1);	/* report points with positive central weight */
  ACCEPTI(enforce_positivity,0,0,1);	/* zero their negative links */
  if (dim==3 && (S->check_positivity || S->enforce_positivity))
    EXPECTED_ERROR("check_positivity and enforce_positivity are not implemented in 3D");
  #if MPI
  if (advance) run_inner_bc (dev->s,S,dev->sync,dev->alwaysRun);
  #else
  if (advance) run_inner_bc (dev->s,S);
  #endif
}
CREATE_TAIL(inner_bc,1)
