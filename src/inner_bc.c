/**
 * >Copyright (C) (2010-2026) Vadim Biktashev, Irina Biktasheva et al. 
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
/* This version is not yet for MPI. */

#include <math.h>
#include <assert.h>
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
} STR;

enum i1_nbrs {i10, i1p, i1m, i1nnb};
enum i2_nbrs {i200, i2p0, i2m0, i20p, i20m, i2nnb};
enum m2_nbrs {m200, m2p0, m2m0, m20p, m20m,
	      m2pp, m2pm, m2mp, m2mm, m2nnb};
enum i3_nbrs {i3000, i3p00, i3m00, i30p0, i30m0, i300p, i300m, i3nnb};
enum m3_nbrs {m3000, m3p00, m3m00, m30p0, m30m0, m300p, m300m,
	      m3pp0, m3pm0, m3mp0, m3mm0, m3p0p, m3p0m,
	      m3m0p, m3m0m, m30pp, m30pm, m30mp, m30mm, m3nnb};

static real kslope (INT x1, INT y1, INT z1, INT x2, INT y2, INT z2, INT v0, real ca, real cb,int manypoint) {
  real c1, c2;
  real v1=New[ind(x1,y1,z1,v0)];
  real v2=New[ind(x2,y2,z2,v0)];
  real thisv, Q;
  real Sx, Sy, S, dx, dy;
  long x, y, X, Y, N;
  real alpha, p, q, k, cosalpha, sinalpha;
  
  ASSERT(dim==2);	/* do 3D later */
  ASSERT(z1==z2);
  
  if (v1==ca) c1=ca;
  else if (v1==cb) c1=cb;
  else ABORT("New[(x1=%ld,y1=%ld,z1=%ld,v0=%ld)]=%g neither %g nor %g\n",
	     (long)x1,(long)y1,(long)z1,(long)v0,ca,cb);
  
  if (v2==ca) c2=ca;
  else if (v2==cb) c2=cb;
  else ABORT("New[(x2=%ld,y2=%ld,z2=%ld,v0=%ld)]=%g neither %g nor %g\n",
	     (long)x2,(long)y2,(long)z2,(long)v0,ca,cb);

  S=Sx=Sy=0.0;
  N=X=Y=0;
  
#define account(dx,dy)				\
  x=x1+dx; y=y1+dy;				\
  thisv=New[ind(x,y,z1,v0)];			\
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
  thisv=New[ind(x,y,z2,v0)];			\
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

  ASSERT(N!=0);
  dx=fabs(Sx-X*S/N);
  dy=fabs(Sy-Y*S/N);
  ASSERT(dx!=0 || dy!=0);

  alpha=atan2(dy,dx);
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

  int nn=0;
#define adjust(dx,dy,dz,qq,wgt)				    \
  if (isTissue(*x+dx,*y+dy,*z+dz)) {                        \
    real val0=New[ind(*x,*y,*z,v0)];                        \
    real val1=New[ind(*x+dx,*y+dy,*z+dz,v0)];               \
    if ((val0==c0 && val1==c1) || (val0==c1 && val1==c0)) { \
      real *pwgt=New+ind(*x,*y,*z,w0+qq);		    \
      real k=kslope(*x,*y,*z,*x+dx,*y+dy,*z+dz, v0, c0, c1, manypoint);	\
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
RUN_TAIL(inner_bc)

DESTROY_HEAD(inner_bc)
DESTROY_TAIL(inner_bc)

CREATE_HEAD(inner_bc)
{
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
  #if MPI
  if (advance) run_inner_bc (dev->s,S,dev->sync,dev->alwaysRun);
  #else
  if (advance) run_inner_bc (dev->s,S);
  #endif
}
CREATE_TAIL(inner_bc,1)
