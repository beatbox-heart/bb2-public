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

/* diff.c: 
 * "Universal" diffusion device to supersede many (all) previous devices of the kind. 
 * Key new features: 
 * # The weights are in device memory OR in the dedicated layers of the main grid (# depending on stencil);
 * # The list of points defined by the same mechanism as generic device space definition;
 * # (Re)calculation of the weights by demand when the corresponding flag is raised, not just at the start;
 * # Several different methods of defining the weights, including some old and some new: 
 * # - isotropic without diagonal connections ("manypoint=0") (diff, ANISO=0);
 * # - isotropic with diagonal connections ("manypoint=1") (diff, ANISO=0);
 * # - ..., variable diffusivity that may be space and/or time-dependent;
 * # - anisotropic with diagonal connections based on fibre directions with Dl and Dt fixed pars;
 * # - ..., ..., Dl and Dt may be time and/or space-dependent; 
 * # All of these work in MPI, too. 
 * 
 */

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
#include "geometry.h"

#include "diffweights.h"		/* stencils, shared with other devices */

typedef struct {
#include "k_code.h"			
#include "diffweights_str.h"		/* the fields used by the weights code */
  real ht;				/* Time step or 0 */
  int v0;				/* Layer for the diffusive field */
  int v1;				/* Layer for the Laplacian */
  ssize_t H[a3nnb];			/* Vector of displacements in the stencil; allocate for max */
  Condition reweight_when;		/* When to recompute the weights: variable address */
  char reweight[maxname];		/* When to recompute the weights: variable name */
  int unstable_warned;			/* Whether instability warning has been issued after a reweight */
} STR;

#include "diffweights_code.h"		/* the weights code, shared with other devices */


/* Maximal stable time step for the current weights */
static real diff_htmax (STR *S,Space *s);


/****************/
RUN_HEAD(diff)
{
  DEVICE_CONST(real,ht);		/* time step */
  DEVICE_CONST(int,v0);			/* source layer (diffusive field) */
  DEVICE_CONST(int,v1);			/* target layer (laplacian field) */
  DEVICE_CONST(int,nnb);		/* number of neighbours in the stencil */
  DEVICE_ARRAY(ssize_t,H);		/* vector of neighbours' displacements */
  DEVICE_CONST(Condition,reweight_when);/* flag to recompute weights */
  size_t np=s.np;			/* number of points in the device's space */
  int V0=v0*DV;				/* displacement from anchor to source value */
  int V1=v1*DV;				/* displacement from anchor to target value */
  long ip;				/* points counter */
  devicePoint P;			/* point structure */
  real *X;				/* vector of weights */
  real *u, *un;				/* pointers to source and target */
  ssize_t H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,
    H10,H11,H12,H13,H14,H15,H16,H17,H18;/* neighbours' displacements (max no) */

  /* Recompute weights if requested explicitly or the space has been relisted */
  /* NB the original, not the copy s, to keep the reset 'relisted' and serve the 'where' code */
  if (*reweight_when || s.relisted) {
    if (SUCCESS != diff_connect(S,&(dev[idev].s))) EXPECTED_ERROR("could not (re)make the connection weights");
    /* New weights may violate the stability limit; warn only once per device */
    if (ht>0 && !S->unstable_warned) {
      real htmax=diff_htmax(S,&s);
      if (ht>htmax) {
	MESSAGE("\n/* WARNING: t=%ld device '%s': after reweighting, ht=%g exceeds "
		"the stability limit %g of the forward Euler step */\n",t,dev[idev].n,ht,htmax);
	S->unstable_warned=1;
      }
    }
  }

  /* Simple variables are faster than array elements */
  switch (nnb) {
  case  3: H0=H[0]; H1=H[1]; H2=H[2]; break; /* i1 */
  case  5: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; break; /* i2 */
  case  7: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; H5=H[5]; H6=H[6]; break; /* i3 */
  case  9: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; H5=H[5]; H6=H[6]; H7=H[7]; H8=H[8]; break; /* m2=a2 */
  case 19: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; H5=H[5]; H6=H[6]; H7=H[7]; H8=H[8]; H9=H[9];
    H10=H[10]; H11=H[11]; H12=H[12]; H13=H[13]; H14=H[14]; H15=H[15]; H16=H[16]; H17=H[17]; H18=H[18]; break; /* m3=a3 */
  default: EXPECTED_ERROR("unknown stencil size %d",nnb);
  }

#if 1
  /* TODO: swap nnb switch with ip loop for speed? */
  /* Paradoxically, in debug compile, loop-outside is faster */
  for (ip=0;ip<np;ip++) {
    P=s.p[ip];
    u=P.u+V0;
    un=P.u+V1;
    X=(real *)P.X;
    switch (nnb) {
    case  3: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]; break;
    case  5: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]; break;
    case  7: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]; break;
    case  9: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8]; break;
    case 19: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8]
	+u[H9]*X[9]+u[H10]*X[10]+u[H11]*X[11]+u[H12]*X[12]+u[H13]*X[13]+u[H14]*X[14]+u[H15]*X[15]+u[H16]*X[16]
	+u[H17]*X[17]+u[H18]*X[18]; break;
    default: EXPECTED_ERROR("unknown stencil size %d",nnb);
    }
  }
#else
  switch (nnb) {
  case  3:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2];
    }
    break;
  case  5:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4];
    }
    break;
  case 7:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6];
    }
    break;
  case  9:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8];
    }
    break;
  case 19:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8]
	+u[H9]*X[9]+u[H10]*X[10]+u[H11]*X[11]+u[H12]*X[12]+u[H13]*X[13]+u[H14]*X[14]+u[H15]*X[15]+u[H16]*X[16]
	+u[H17]*X[17]+u[H18]*X[18];
    }
    break;
  default: EXPECTED_ERROR("unknown stencil size %d",nnb);
  } /* switch nnb */
#endif

  /* If nonzero timestep set, do the forward Euler step for diffusion equation */
  if (ht) {
#if MPI
    MPIDO(MPI_Barrier(ALL_ACTIVE_PROCS),"could not sync threads before applying timestep");
#endif
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u;
      u[V0] += ht*u[V1];
    } /* for ip */
  } /* if ht */
}
RUN_TAIL(diff)

/****************/
DESTROY_HEAD(diff)
{
  FREE(S->Xpool);
  FREE(S->Ypool);
  FREE(S->u);
  FREE(S->geom);
  FREE(S->loctb);
  FREE(S->Dsrc);
  FREE(S->Dlsrc);
  FREE(S->Dtsrc);
}
DESTROY_TAIL(diff)

CREATE_HEAD(diff)
{
  Space *s=&(dev->s);
  DEVICE_REQUIRES_SYNC;
  DEVICE_ALWAYS_RUNS;			/* inactive instances must join haloSwap() and MPI_Barrier() */
  if (s->nowhere) EXPECTED_ERROR("'nowhere' is not allowed\n");
  if (find_key("listpoints=",w)) {	/* say something only if the user has a say */
    if (s->listpoints==0)
      MESSAGE("/* listpoints=1 is mandatory in this device and enforced */\n");
    else
      MESSAGE("/* listpoints=1 is mandatory anyway */\n");
  }
  if (!s->listpoints) {			/* otherwise already listed by accept_space() */
    s->listpoints=1;
    s->p=NULL; /* signal that it has not be alloc'd yet */
    if (SUCCESS != list_space(s))
      ABORT("could not make list of points");
  }
    
  size_t np=s->np;			/* number of listed points of this device */
  size_t ip;				/* listed points counter */
  int reweight_entry;			/* deftb entry of 'reweight' k-var */

  ACCEPTI(v0,INONE,0,vmax-1);		/* source layer: diffusive field */
  ACCEPTI(v1,INONE,0,vmax-1);		/* target layer: laplacian field */
  ASSERT( v1 != v0 );			/* source and target must be distinct */

  ACCEPTR(hx,RNONE,RSUCC(0.),RNONE);	/* diffusion makes no sense without space step */
  ACCEPTR(ht,0.,0.,RNONE);		/* make forward Euler step if ht>0 */
  int nb;				/* small range counter */
  int nnb, nD;				/* stencil size, num of tensor components */
  int hasexe;				/* whether any diffusivity is a =~ expression */
  /* Parameters of the weights, common with other devices, see diffweights_code.h */
  if (!diff_weights_accept(S,w,&hasexe)) return 0;
  nnb=S->nnb;
  nD=S->nD;
  int manypoint=S->manypoint;

  ACCEPTS(reweight,"never");		/* global k-var to flag that weights to be recomputed */
  if (0==(reweight_entry=tb_find(deftb,S->reweight)))
    EXPECTED_ERROR("unknown symbol %s for 'when' parameter",S->reweight);
  if ((tb_flag(deftb,reweight_entry) & f_vb)==0)
    EXPECTED_ERROR("symbol %s used for 'when' parameter is not a variable",S->reweight);
  S->reweight_when = (Condition)tb_addr(deftb,reweight_entry); 


  /* Border of the device space: no-flux, or connected to the tissue beyond?	*/
  /* The latter is safe if the diff tensor is constant, so default then.	*/
  /* Any =~ expression, even of global variables only, defaults to no-flux,	*/
  /* for compatibility with the earlier versions.				*/
  ACCEPTI(isolated,(S->loctb || hasexe || ANISOTROPY_ON)?1:0,0,1);

  /* The stability limit of the time step is checked against the actual weights, see diff_htmax() */

  /* Range of layers to keep stencil weights */
  ACCEPTI(w0,-1,-1,vmax-1);
  ACCEPTI(w1,-1,-1,vmax-1);
  if (!((w0==-1 && w1==-1) || (0<=w0 && w0<=w1 && w1<vmax && w1-w0+1==nnb))) {
    MESSAGE("w0=%d w1=%d vmax=%d nnb=%d\n",w0,w1,vmax,nnb);
    ASSERT((w0==-1 && w1==-1) || (0<=w0 && w0<=w1 && w1<vmax && w1-w0+1==nnb));
  }
 

  /* Determine displacements to neighbours == stencil elements */
  ssize_t *H=S->H;			/* filled in directly; S->H is allocated for max size */
  if (ANISOTROPY_ON) {			/* anisotropic */
    switch (dim) {
    case 1:
      EXPECTED_ERROR("don't know how to make 1D anisotropic");
      break;
    case 2:
      H[a200]=0;
      H[a2p0]=+DX; H[a2m0]=-DX;
      H[a20p]=+DY; H[a20m]=-DY;
      H[a2pp]=+DX+DY; H[a2pm]=+DX-DY; H[a2mp]=-DX+DY; H[a2mm]=-DX-DY;
      break;
    case 3:
      H[a3000]=0;
      H[a3p00]=+DX; H[a3m00]=-DX;
      H[a30p0]=+DY; H[a30m0]=-DY;
      H[a300p]=+DZ; H[a300m]=-DZ;
      H[a3pp0]=+DX+DY; H[a3pm0]=+DX-DY; H[a3mp0]=-DX+DY; H[a3mm0]=-DX-DY;
      H[a3p0p]=+DX+DZ; H[a3p0m]=+DX-DZ; H[a3m0p]=-DX+DZ; H[a3m0m]=-DX-DZ;
      H[a30pp]=+DY+DZ; H[a30pm]=+DY-DZ; H[a30mp]=-DY+DZ; H[a30mm]=-DY-DZ;
      break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    } /* switch dim */
  } else if (manypoint==0) {		/* isotropic, small stencil */
    switch (dim) {
    case 1:
      H[i10]=0;
      H[i1p]=+DX; H[i1m]=-DX;
      break;
    case 2:
      H[i200]=0;
      H[i2p0]=+DX; H[i2m0]=-DX;
      H[i20p]=+DY; H[i20m]=-DY;
      break;
    case 3:
      H[i3000]=0;
      H[i3p00]=+DX; H[i3m00]=-DX;
      H[i30p0]=+DY; H[i30m0]=-DY;
      H[i300p]=+DZ; H[i300m]=-DZ;
      break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    } /* switch dim */
  } else { 				/* isotropic, large stencil */
    switch (dim) {
    case 2:
      H[m200]=0;
      H[m2p0]=+DX; H[m2m0]=-DX;
      H[m20p]=+DY; H[m20m]=-DY;
      H[m2pp]=+DX+DY; H[m2pm]=+DX-DY; H[m2mp]=-DX+DY; H[m2mm]=-DX-DY;
      break;
    case 3:
      H[m3000]=0;
      H[m3p00]=+DX; H[m3m00]=-DX;
      H[m30p0]=+DY; H[m30m0]=-DY;
      H[m300p]=+DZ; H[m300m]=-DZ;
      H[m3pp0]=+DX+DY; H[m3pm0]=+DX-DY; H[m3mp0]=-DX+DY; H[m3mm0]=-DX-DY;
      H[m3p0p]=+DX+DZ; H[m3p0m]=+DX-DZ; H[m3m0p]=-DX+DZ; H[m3m0m]=-DX-DZ;
      H[m30pp]=+DY+DZ; H[m30pm]=+DY-DZ; H[m30mp]=-DY+DZ; H[m30mm]=-DY-DZ;
      break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    } /* switch dim */
  } /* if ANISOTROPY_ON else if manypoint==0 else */

  /* Now compute and assign connection weights: this is also potentially doable at each step */
  if (SUCCESS != diff_connect(S,s)) ABORT("could not make the connection weighs");

  /* Connection to adjacent tissue is dangerous unless that is deliberate */
  long nadjacent=S->nadjacent;
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,&nadjacent,1,MPI_LONG,MPI_SUM,ALL_ACTIVE_PROCS),"could not sum up the adjacent tissue points");
#endif
  if (nadjacent)
    MESSAGE("\n/* WARNING: the device space borders on other tissue, and with isolated=0 is connected to it.\n"
	    "   This is correct only if the tissue beyond the border is taken care of, e.g. by another\n"
	    "   diff device or by boundary conditions; otherwise use isolated=1. */\n");

  /* Stability limit of the forward Euler step, from the actual weights */
  if (ht>0) {
    real htmax=diff_htmax(S,s);
    if (ht>htmax) {
      MESSAGE("\n/* WARNING: ht=%g exceeds the stability limit %g of the forward Euler step */\n",ht,htmax);
      S->unstable_warned=1;
    } else {
      MESSAGE("\x01\n/* ht/htmax=%g, htmax=%g */\n",ht/htmax,htmax);
    }
  }
}
CREATE_TAIL(diff,1);

/****************/
/* Maximal stable time step of the forward Euler step with the current weights. 	*/
/* By Gershgorin's theorem, the eigenvalues of the discrete operator lie within	*/
/* max_i (|W_i0| + sum_nb |W_i,nb|) from zero, as W_i0 = -sum_nb W_i,nb; hence 	*/
/* ht <= 2/max_i(...) is sufficient. With nonnegative peripheral weights this is	*/
/* ht*max_i|W_i0| <= 1, which is also the condition of the maximum principle and	*/
/* the exact stability limit for the standard stencils, e.g. D*ht/hx^2 <= 1/(2*dim).	*/
/* Collective in MPI: all instances must call it.					*/
static real diff_htmax (STR *S,Space *s)
{
  int nnb=S->nnb, nb;
  size_t ip;
  double wmax=0.0, w;
  real *X;
  for (ip=0;ip<s->np;ip++) {
    X=(real *)(s->p[ip].X);
    for (w=0.0,nb=0;nb<nnb;nb++) w+=fabs(X[nb]);
    if (w>wmax) wmax=w;
  }
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,&wmax,1,MPI_DOUBLE,MPI_MAX,ALL_ACTIVE_PROCS),"could not find the maximal weight");
#endif
  return (wmax>0.0)?(2.0/wmax):HUGE_VAL;
} /* diff_htmax */

