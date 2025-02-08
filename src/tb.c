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
 * Truscott-Brindley 1994 model for phyto+zoo
 * (rescaled from liters to m^3)
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
  #include "tbvar.h"
  #undef _
  Nvar /* total number of variables */
};

typedef struct {
				/* parameters */
  #define _(n,v,d,unit) real n;
  #include "tbpar.h"
  #undef _
				/* external currents */
  #define _(n,v,d,unit) real I##n;
  #include "tbvar.h"
  #undef _
} STR;

RHS_HEAD(tb,Nvar) {
  #define _(n,v,d,unit) DEVICE_CONST(real,n)
  #include "tbpar.h"
  #undef _
  #define _(n,v,d,unit) DEVICE_CONST(real,I##n)
  #include "tbvar.h"
  #undef _
  #define _(n,v,d,unit) real n = u[var_##n];
  #include "tbvar.h"
  #undef _
/*   #define _(n,v,d,unit) real n##new; */
/*   #include "tbvar.h" */
/*   #undef _ */

  real P2;	/* TB just P^2 */
  real G;	/* TB zoo grazing (fraction of phyto eaten this day) */
  real alp2;	/* TB just alpTB^2  */
  real incr;	/* TB zoo relative increment */
  real decr;	/* TB zoo relative decrement */
			/* Relative grazing intensity as function of age */
  /*   #define app(W) (pow((W)/appwgt,-appexp)) */

  /************************/
  /* phyto dynamics by TB */
  /************************/
  P2=P*P;
  alp2=alpTB*alpTB;
  G = Rm*P2/(alp2+P2);
  du[var_P] = (r*P*(1-P/Ktb) - G*Z) + IP;

  /**********************************************/
  /* Zoo: increment by TB and decrement by both */
  /**********************************************/
      
  /* TB formulation: Z' = Z + dt*Z*( gam*G - mu ); */
  incr = gam*G*Z;
  decr = mu*Z;
				/*eaten+=CH*h*Pchtot;*/
				/*died+=TB*mu;*/
  du[var_Z] = incr - decr + IZ;

} RHS_TAIL(tb)

RHS_CREATE_HEAD(tb) {

  #define _(n,v,d,unit) ACCEPTP(n,v,RNONE,RNONE);
  #include "tbpar.h"
  #undef _

  #define _(n,v,d,unit) ACCEPTP(I##n,0,RNONE,RNONE);
  #include "tbvar.h"
  #undef _
  
/*   MALLOC(*u,Nvar*sizeof(real)); */
/*   for(i=0;i<Nvar;i++) { */
/*     #define _(n,v,d,unit) (*u)[var_##n]=v; */
/*     #include "tbvar.h" */
/*     #undef _ */
/*   } */

  MALLOC(*u,Nvar*sizeof(real)); 
  (*u)[var_P]=(S->alpTB)/sqrt(1+(S->gam)*(S->Rm)/(S->mu));
  (*u)[var_Z]=(S->r)/(S->Rm)*(*u)[var_P];
  MESSAGE("/*P0=%g Z0=%g*/",(*u)[var_P],(*u)[var_Z]);

} RHS_CREATE_TAIL(tb,Nvar)
