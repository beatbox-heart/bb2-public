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

/*
 * Modified hybrid of Cushing-Horwood 1994 and Truscott-Brindley 1994
 * models for phyto+zoo and zoo+fish larva
 * TB rescaled from liters to m^3.
 * Modification: regularisation to avoid negative concentrations/
 * weights, provide smoothness etc.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"
#include "mchtb.h"

				/* indices of dynamic variables */
enum {
  #define _(n,v,d,unit) var_##n, 
  #include "mchtbvar.h"
  #undef _ 
  Nvar /* total number of variables */
};

typedef struct {
				/* parameters */
  #define _(n,v,d,unit) real n;
  #include "mchtbpar.h"
  #undef _ 
				/* external currents */
  #define _(n,v,d,unit) real I##n;
  #include "mchtbvar.h"
  #undef _
} STR;

		    /* volume searched by individual larva of weight w per day */
#define vol(W) (k*pow(W,nu))

		    /* Maximal larvae growth rate (foraging capacity) */
#if VARGMAX
  #define Gl(A) (Gmax*pow(1+(A),-Glexp))
#else
  #define Gl(A) (Gmax)
#endif

#define theta(x) (0.5*(1+tanh(x)))

			    /* metamorphosis rate depends on age and weight */
#define metamrate(A,W) ( \
  metagecoeff*theta(((A)-metage)/metagealw) \
 +metwgtcoeff*theta(((W)-metwgt)/metwgtalw) \
)
			    /* death of starvation rate depends on
			     * weight gain deficit agains metabolic demand
			     */
			     /* replace exp -> theta ? */
#define starverate(dW,Spent) (Spent?starvcoeff*exp(-starvexp*(dW)/(Spent)):0)

#define predationrate(A) (Mmax/(1+b*(A)))

real nonnegative(real w);
void combine(real N1, real W1, real N2, real W2, real *Nn, real *Wn);
real Holling(real Rmax,real p,real z,real type);

RHS_HEAD(mchtb,Nvar) {
  #define _(n,v,d,unit) DEVICE_CONST(real,n)
  #include "mchtbpar.h"
  #undef _ 
  #define _(n,v,d,unit) DEVICE_CONST(real,I##n)
  #include "mchtbvar.h"
  #undef _ 
  #define _(n,v,d,unit) real n = u[var_##n];
  #include "mchtbvar.h"
  #undef _ 

  real alp2;	/* TB just alpTB^2  */
  real P2;	/* TB just P^2 */
  real G;	/* TB zoo grazing (fraction of phyto eaten this day) */

  real W;	/* CH larvae weight = B/N */
  #if VARKCH
    real Kch;	/* CH metabolic wgt-dep coeff */
  #endif
  real bet;	/* CH metabolic wgt-dep coeff */
  real pch;	/* CH max fraction of volume searched */
  real R;	/* CH larvae ration per ug larvae per ug zoo per day */
  real Rmax;	/* CH max larvae ration (needs) */
  real GmaxA;	/* CH max larvae growth rate with acct of age */
  real Spent, Obtained; /* CH metabolic balance */

  real dWmax;	/* CH max wgt-dep growth */
  real dW;	/* CH de facto weight gain */
  real dN;	/* CH number gain */

  real eaten;	/* CH larvae predation mortality rate */
  real starved; /* new larvae starvation mortality rate */
  real metamorphosed; /* CH larvae metamorphosis rate */
  
  real incr;	/* TB zoo relative increment */
  real decr;	/* TB+CH zoo relative decrement */


  #if NONNEG_P
    if (P<0) P=0;
  #endif
  #if NONNEG_Z
    if (Z<0) Z=0;
  #endif
  #if NONNEG_B
    if (B<0) B=0;
  #endif
  #if NONNEG_N
    if (N<0) N=0;
  #endif

  if (N) W=B/N; else W=Wi;

  du[var_A]=1;

  /************************/
  /* phyto dynamics by TB */
  /************************/
  P2=P*P;
  alp2=alpTB*alpTB;
  G = Rm*P2/(alp2+P2);
  du[var_P] = TB*(r*P*(1-P/Ktb) - G*Z) + IP;

  /********************************/
  /* larvae weight dynamics by CH */
  /********************************/
  GmaxA = Gl(A);			/* maximal growth rate */
  dWmax = W*GmaxA;			/* max wgt-dep growth */
  bet = betmax-dbet*exp(-j*W);		/* conversion efficiency */
  #if VARKCH
    Kch = (W<1000)?2.57:2.64;		/* wgt-dep metabolism */
  #endif
  Spent = Kch*pow(W,n);		        /* metabolic costs */ 
  Rmax = (dWmax + Spent)/bet;		/* maximal ration */
  pch = k*pow(W,nu);			/* max fraction of volume searched */
  R = Holling(Rmax,pch,Z,Hollingtype);	/* Holling's I/II larvae/zoo "encounter rate" */
  Obtained=bet*R*Z;			/* actual biomass ackquired */ 
  dW=CH*(Obtained-Spent);		/* actual weight gain per day */

  /********************************/
  /* larvae number dynamics by CH */
  /********************************/
  eaten=predationrate(A);
  starved=starverate(dW,Spent);    
  metamorphosed=metamrate(A,W); 
				    if(metamorphosed>1e10) 
				    ABORT("metamorphosed=%lg too big\nA=%lg W=%lg\n",
				      metamorphosed, A, W
				    );
  dN = - CH*N*(eaten+starved+metamorphosed);
  
  du[var_N] = dN + IN;
  du[var_B] = N*dW + W*dN + IB;

  /************************/
  /* metamorphosed larvae */
  /************************/
  du[var_FN] = CH*N*metamorphosed + IFN;
  du[var_FB] = CH*N*W*metamorphosed + IFB;

  /**********************************************/
  /* Zoo: increment by TB and decrement by both */
  /**********************************************/
  /* CH formulation: z' = z*( g - P );              */
  /* TB formulation: Z' = Z*( gam*G - mu );         */
  /* CH+TB formulation: Z' = Z*( gam*G - RN - mu ); */
  incr = TB*gam*G + (1-TB)*g;
  decr = CH*R*N + TB*mu;
  du[var_Z] = Z*(incr - decr) + IZ;

} RHS_TAIL(mchtb)

static real sqr(real a) {return a*a;}

RHS_CREATE_HEAD(mchtb) {
  int i;
  #define _(n,v,d,unit) ACCEPTP(n,v,RNONE,RNONE);
  #include "mchtbpar.h"
  #undef _ 


  #define _(n,v,d,unit) ACCEPTP(I##n,0,RNONE,RNONE);
  #include "mchtbvar.h"
  #undef _ 
  
  MALLOC(*u,Nvar*sizeof(real));

  for(i=0;i<Nvar;i++) {
    #define _(n,v,d,unit) (*u)[var_##n]=v;
    #include "mchtbvar.h"
    #undef _ 
  }

  /* these does not seem to be correct formulae? - 2001/08/20 */
  /* Pure TB model, find exact equilibirum */
/*   if (S->CH==0) { } */
/*   if (S->TB!=0) {  */
/*     real Peq = (S->alpTB)/( (S->gam)*(S->Rm)/(S->mu) - 1. ); */
/*     real Zeq = (1.+sqr(S->alpTB/Peq))*(S->r)*Peq*(1.-Peq/(S->Ktb))/(S->Rm); */
/*     (*u)[var_P]=Peq; */
/*     (*u)[var_Z]=Zeq; */
/*   } */

  /* added 2001/08/20 */
  /* coexistence phyto+zoo w/o fish */
   if (S->TB!=0) {  
     real Peq = (S->alpTB)/sqrt( (S->gam)*(S->Rm) / (S->mu) - 1.0 ); 
     real Zeq = (S->gam)*(S->r)*Peq*((S->Ktb)-Peq)/((S->Ktb)*(S->mu)); 
     (*u)[var_P]=Peq; 
     (*u)[var_Z]=Zeq; 
     (*u)[var_N]=0; 
     (*u)[var_B]=0; 
     (*u)[var_N]=0; 
     (*u)[var_FN]=0; 
     (*u)[var_FB]=0; 
     (*u)[var_A]=0; 
   } 

} RHS_CREATE_TAIL(mchtb,Nvar)

real Holling(real Rmax,real p,real z,real type) {
  if (type>1.5) return (Rmax*p)/(Rmax + p*z);
  else return min(p, Rmax/z);
}

