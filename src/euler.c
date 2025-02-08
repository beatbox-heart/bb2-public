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

/* General pointwise ODE stepper */
/* The specific ODE rhs (e.g. cell model module) is attached as a separate function, or cell model module */

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

typedef struct {
  real ht;		/* time step */
  Name ode;		/* name of the ODE system */
  rhs_str rhs;		/* right-hand side description */
  int rest;		/* how many steps to do to find resting values */
  int v0, v1;		/* layer range */
  real *u;		/* vector of vars for steady state if any */
  real *du;		/* vector of derivatives */
} STR;

/****************/
RUN_HEAD(euler)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_ARRAY(real,du);
  DEVICE_ARRAY(char,ode);
  DEVICE_CONST(real,ht);
  RHS_CONST(RhsProc *,f);
  RHS_CONST(Par,p);
  RHS_CONST(Var,var);
  int nv=v1-v0+1;
  int iv;

#define COMMANDS {	      						\
  if (!f(u+v0*DV,du,p,var,nv)) { /* u start from layer 0, base layer is v0 */ \
    URGENT_MESSAGE("error calculating %s at t=%ld point %ld %ld %ld: u=(",ode,t,xof(u),yof(u),zof(u)); \
    for(iv=0;iv<nv;iv++) URGENT_MESSAGE(" %lg",u[v0+iv]);		\
    ABORT(")\n");							\
  } /*  if !f... */							\
  for (iv=0;iv<nv;iv++) {						\
    u[v0+iv] += ht*du[iv]; /* u start from layer 0, base layer is v0 */	\
    if (!isfinite(u[v0+iv])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) detected at t=%ld x=%d y=%d z=%d v=%d\n",t,xof(u),yof(u),zof(u),v0+iv); \
      URGENT_MESSAGE("This happened after an increment by %lg*",ht);	\
      for(iv=0;iv<nv;iv++) URGENT_MESSAGE("%c%lg",iv?',':'(',du[iv]);	\
      URGENT_MESSAGE(")\n");							\
      ABORT("");							\
    }									\
  }									\
}
  DO_FOR_ALL_POINTS;
#undef COMMANDS
}
RUN_TAIL(euler);

#undef CHECK
#undef F

/****************/
DESTROY_HEAD(euler)
{
  FREE(S->u);
  FREE(S->du);
  if (S->rhs.var.n) {
    FREE(S->rhs.var.src);
    FREE(S->rhs.var.dst);
  }
  FREE(S->rhs.p);
}
DESTROY_TAIL(euler)

/* Declare all available ODE's */
#define D(a) RhsProc a;
#include "rhslist.h"
#undef D
#define D(a) RhsCreate create_##a;
#include "rhslist.h"
#undef D

/****************************************/
CREATE_HEAD(euler)
{
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 >= v0 );
  int nv=v1-v0+1;
  int step, iv, ix, iy, iz, NV;
  real *ufull, *u;

  ACCEPTR(ht,RNONE,0.,RNONE); /* ht=0 can make sense in some cases */
  ACCEPTS(ode,NULL);
  {
    char *pars;
    MALLOC(pars,(long)MAXSTRLEN);
    BEGINBLOCK("par=",pars);
    S->u=NULL;
    #define D(a)							\
    if (0==stricmp(S->ode,#a)) {					\
      NV=create_##a(&(S->rhs.p),&(S->rhs.var),pars,&(S->u),v0);		\
      if (!NV) EXPECTED_ERROR("reading parameters for %s in \"%s\"",S->ode,pars);	\
      if (NV!=nv) EXPECTED_ERROR("%d layers instead of %d for %s",nv,NV,S->ode);	\
      S->rhs.f=a;							\
    } else
    #include "rhslist.h"
    #undef D
    EXPECTED_ERROR("unknown ODE %s",S->ode);
    ENDBLOCK;
    FREE(pars);
  }
  CALLOC(S->du,nv,sizeof(real));

  ACCEPTI(rest,0,0,INONE);
  if (S->rest) {
    /* Need full vector, in case variable parameters use other layers.   */
    /* NB other layers may be above as well as below this device's space.*/
    CALLOC(ufull,vmax,sizeof(real)); 
    u=&(ufull[v0]);
    if (S->u) {
      MESSAGE("\n/* NOTICE: finding resting state while already defined by the rhs */");
      /* in that case, ode-defined state will be initial condition for iterations */
      for (iv=0;iv<nv;iv++) ufull[v0+iv]=S->u[iv];
    } else {
      CALLOC(S->u,nv,sizeof(real));
    }

    if (S->rhs.var.n) MESSAGE("\n/* ATTENTION: you asked to find the resting state while "
			      "some parameters are potentially variable */");
    
    for(step=0;step<S->rest;step++) {
      if( !S->rhs.f(u,S->du,S->rhs.p,S->rhs.var,nv) ) {
	URGENT_MESSAGE("error calculating %s at step %d when initiating device %s: u=",S->ode,step,dev->n);
	for(iv=0;iv<nv;iv++) URGENT_MESSAGE(" %lg",u[iv]);
	ABORT("\n");
	/* FFLUSH(stdout); - never executed! */
      }
      for (iv=0;iv<nv;iv++) {
	u[iv] += S->ht * S->du[iv];
	if (0==isfinite(u[iv])) {
	  URGENT_MESSAGE("\nNAN (not-a-number) detected at step %d when initiating device %s: v=%d",step,iv);
	  URGENT_MESSAGE("This happened after an increment by %lg*",S->ht);
	  for(iv=0;iv<nv;iv++) URGENT_MESSAGE("%c%lg",iv?',':'(',S->du[iv]);
	  URGENT_MESSAGE(")\n");
	  return 0;
	}
      }
    } /* for step */

    /* copy what is needed */
    for (iv=0;iv<nv;iv++) S->u[iv]=u[iv];

    /* don't need this one any more */
    FREE(ufull);
  } /* if S->rest */

  if (S->u) {
    MESSAGE0("\n/* Resting state: "); for(iv=0;iv<nv;iv++) MESSAGE1("%lg ",(S->u)[iv]); MESSAGE0("*/");
    /* Fill up the whole of the grid with the resting state */
    #if MPI
    if (dev->s.runHere) {
    #endif
      for (ix=dev->s.x0;ix<=dev->s.x1;ix++) { 
	for (iy=dev->s.y0;iy<=dev->s.y1;iy++) {
	  for (iz=dev->s.z0;iz<=dev->s.z1;iz++) {
	    for (iv=v0;iv<=v1;iv++) {
	      New[ind(ix,iy,iz,iv)]=(S->u)[iv-v0];
	    } /* for iv */
	  } /* for iz */
	} /* for iy */
      } /* for ix */
    #if MPI
    } /* if runHere */
    #endif
  } else {
    /* No resting state was defined */
    CALLOC(S->u,nv,sizeof(real));
  }
}
CREATE_TAIL(euler,1)
