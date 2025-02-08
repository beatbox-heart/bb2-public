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
 * Hybrid of Cushing-Horwood 1994 and Truscott-Brindley 1994
 * models for phyto+zoo and zoo+fish larva
 * (TB rescaled from liters to m^3)
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

				/* indices of dynamic variables */
enum {
  #define _(n,v,d,unit) var_##n, 
  #include "chtbvar.h"
  #undef _
  
  Nvar /* total number of variables */
};

typedef struct {
				/* parameters */
  #define _(n,v,d,unit) real n;
  #include "chtbpar.h"
  #undef _ 
				/* external currents */
  #define _(n,v,d,unit) real I##n;
  #include "chtbvar.h"
  #undef _
} rhspar;

RHS_HEAD(chtb,Nvar) {
  #define _(n,v,d,unit) DEVICE_CONST(real,n)
  #include "chtbpar.h"
  #undef _ 
  #define _(n,v,d,unit) DEVICE_CONST(real,I##n)
  #include "chtbvar.h"
  #undef _ 
  #define _(n,v,d,unit) real n = u[var_##n];
  #include "chtbvar.h"
  #undef _ 
  #define _(n,v,d,unit) real n##new;
  #include "chtbvar.h"
  #undef _ 

  real P2;	/* TB just P^2 */
  real G;	/* TB zoo grazing (fraction of phyto eaten this day) */
  real alp2;	/* TB just alpTB^2  */
  real Kch;	/* CH metabolic wgt-dep coeff */
  real bet;	/* CH metabolic wgt-dep coeff */
  real pch;	/* CH max fraction of volume searched */
  real R;	/* CH max larvae ration */
  real dW1;	/* CH max growth from mass/energy balance */
  real dW2;	/* CH max wgt-dep growth */
  real dW;	/* CH de facto weight gain */
  real incr;	/* TB zoo relative increment */
  real decr;	/* TB+CH zoo relative decrement */
			/* Relative grazing intensity as function of age */
  #define app(W) (pow((W)/appwgt,-appexp))

  du[var_a]=1;

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
  Kch = (W<1000)?2.57:2.64;	        /* wgt-dep metabolism */
  bet = betmax-(betmax - betmin)*exp(-j*W);
  pch = k*pow(W,nu);			/* max fraction of volume searched */
  R = pch*h*Z;				/* max larvae ration */
  dW1 = (1-alpCH)*bet*R - Kch*pow(W,n); /* max growth from mass/energy balance */
  dW2 = W*Gmax*app(W);			/* max wgt-dep growth */
  dW =  CH*min(dW1,dW2);			/* de facto weight gain */
  du[var_W] = dW + IW;

  /********************************/
  /* larvae number dynamics by CH */
  /********************************/
  du[var_N] = -CH*N*Mmax/(1+b*a) + IN;	/* NB: age-dep instead of time-dep */

  /*******************************/
  /* How much zoo has been eaten */
  /*******************************/
  R = (dW + Kch*pow(W,n))/((1-alpCH)*bet);  /* de facto larvae ration */

  /**********************************************/
  /* Zoo: increment by TB and decrement by both */
  /**********************************************/
      
  /* CH formulation: z' = z + dt*z*( g - P ); */
  /* TB formulation: Z' = Z + dt*Z*( gam*G - mu ); */
  incr = (TB*gam*G + (1-TB)*g)*Z;
  decr = CH*R*N + TB*mu*Z;
				/*eaten+=CH*h*Pchtot;*/
				/*died+=TB*mu;*/
  du[var_Z] = incr - decr + IZ;

} RHS_TAIL(chtb)

RHS_CREATE_HEAD(chtb) {
  int i;
  #define _(n,v,d,unit) ACCEPTP(n,v,RNONE,RNONE);
  #include "chtbpar.h"
  #undef _ 

  #define _(n,v,d,unit) ACCEPTP(I##n,0,RNONE,RNONE);
  #include "chtbvar.h"
  #undef _ 
  
  MALLOC(*u,Nvar*sizeof(real));
  for(i=0;i<Nvar;i++) {
    #define _(n,v,d,unit) (*u)[var_##n]=v;
    #include "chtbvar.h"
    #undef _ 
  }
} RHS_CREATE_TAIL(chtb,Nvar)
