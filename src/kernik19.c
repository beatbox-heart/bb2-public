/**
 * Copyright (C) (2010-2024) Vadim Biktashev, Irina Biktasheva et al. 
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
  IONIC description of the Kernik et al. 2019 model:

  Kernik, D.C., Morotti, S., Wu, H., Garg, P., Duff, H.J., Kurokawa, J., Jalife, J., Wu, J.C., Grandi, E. and Clancy, C.E., 2019. A computational model of induced pluripotent stem‐cell derived cardiomyocytes incorporating experimental variability from multiple data sources. The Journal of physiology, 597(17), pp.4533-4564.
  This code is based on authors `error-free' Matlab code, reinterpreted to C and converted to the Beatbox 'ionic' format by VNB. 
  NB all essential contents are stoved away into kernik19_*.h files. 
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
#include "ionic.h"
#include "qpp.h"

/* No Markov chains in this model */
#define NMC 0

static inline real sqr(real x) {return x*x;}
static inline real cub(real x) {return x*x*x;}
static inline real eighth(real x) {real x2=x*x; real x4=x2*x2; return x4*x4;}
#define iif(a,b,c) (a)?(b):(c)

/* Transmembrane voltage is the first element of state vector */
static int V_index=0;

/* Constants (unchangeable parameters) */
/* that can be used in tabulated functions. */
#include "kernik19_const.h"

/* Enumerate all dynamic variables */
enum {
  #define _(n,i) var_##n,
  #include "kernik19_other.h"
  #include "kernik19_ngate.h"
  #include "kernik19_tgate.h"
  #undef _
  NV		    		/* total number of variables */
};

/* Enumerate the other (non-gate) variables */
enum {
  #define _(n,i) other_##n,
  #include "kernik19_other.h"
  #undef _
  NO		    		/* total number of other variables */
};

/* Enumerate the non-tabulated gates */
enum {
  #define _(n,i) ngate_##n,
  #include "kernik19_ngate.h"
  #undef _
  NN		    		/* total number of non-tabulated gate variables */
};

/* Enumerate the tabulated gates */
enum {
  #define _(n,i) tgate_##n,
  #include "kernik19_tgate.h"
  #undef _
  NT		    		/* total number of tabulated gate variables */
};

/* Enumerate the tabulated transition rates and other tabulated functions */
enum {
  #define _(n,i) _alp_##n,
  #include "kernik19_tgate.h"
  #undef _
  #define _(n,i) _bet_##n,
  #include "kernik19_tgate.h"
  #undef _
  #define _(n) _##n,
  #include "kernik19_fun.h"
  #undef _
  NTAB				/* total number of tabulated functions */
};

/* The structure containing the parameter values for this instance of the model */
typedef struct {		/* First go the cell parameters as listed in M.C.'s code */
  #define _(name,default) real name;
  #include "kernik19_par.h"
  #undef _
  real IV; 			/* Then the external current (mV/ms) */
} STR;

/* Functions of V to be tabulated, including gates' transition rates */
IONIC_FTAB_HEAD(kernik19) {
  #include "kernik19_ftab.h"
  /* Copy the results into the output array values[]. */
  /* Care is taken that all, and only, tabulated functions are attended here. */
  #define _(n,i) values[_alp_##n]=alp_##n;
  #include "kernik19_tgate.h"
  #undef _
  #define _(n,i) values[_bet_##n]=bet_##n;
  #include "kernik19_tgate.h"
  #undef _
  #define _(n) values[_##n]=n;
  #include "kernik19_fun.h"
  #undef _
  /* if (t>=0) printf("+++++++++++++++++ t=%ld f2_inf=%g, tau_f2=%g\n", t, f2_inf, tau_f2); */
  /* if (t>=0) printf("----------------- alpha[f2]=%g, beta[f2]=%g\n", values[_alp_f2], values[_bet_f2]); */
} IONIC_FTAB_TAIL;

IONIC_FDDT_HEAD(kernik19,NV,NTAB,NO,NN) {
  /* Declare the const pars and take their values from struct S==par (a formal parameter) */
  #define _(name,default) DEVICE_CONST(real,name);
  #include "kernik19_par.h"
  #undef _
  DEVICE_CONST(real,IV);
  /* Declare and assign local variables for dynamic variables from state vector */
  /* ..., first for non-gate variables */
  #define _(name,initial) real name=u[var_##name];
  #include "kernik19_other.h"
  #undef _
  /* ..., for the non-tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "kernik19_ngate.h"
  #undef _
  /* ..., and then for tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "kernik19_tgate.h"
  #undef _
  /* Get the values of the tabulable functions passed through the formal parameter */
  #define _(name) real name=values[_##name];
  #include "kernik19_fun.h"
  #undef _
  /* Calculate the rates of non-gate variables */
  #include "kernik19_fddt.h"
  /* Copy the calculated rates into the output array du[].  */
  /* Care is taken that all, and only, non-gating variables are attended here */
  #define _(name,initial) du[other_##name]=diff_##name;
  #include "kernik19_other.h"
  #undef _
  /* And copy the non-tab transition rates into the output arraysn nalp[] and nbet[] */
  #define _(name,initial) nalp[ngate_##name]=alp_##name; nbet[ngate_##name]=bet_##name; 
  #include "kernik19_ngate.h"
  #undef _
  /* Finally add the "external current" parameter values */
  du[V_index]+=IV;
} IONIC_FDDT_TAIL;

/*
 **********************************************************************************
 * The function initialising an instance of the model.                            *
 * This includes: assigning values for all model parameters,                      *
 * i.e. reading from the script or keeping the default values otherwise;      *
 * and creating the vector of initial (steady-state) values of dynamic variables, *
 * which will be used by euler to initialise the whole medium.                    *
 **********************************************************************************
 */
IONIC_CREATE_HEAD(kernik19) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "kernik19_par.h"
  #undef _

  /* Assign the initial values as given the following *.h files */
  #define _(name,initial) (*u)[var_##name]=initial;
  #include "kernik19_other.h"
  #include "kernik19_ngate.h"
  #include "kernik19_tgate.h"
  #undef _
} IONIC_CREATE_TAIL;



real nalp[NN];
real nbet[NN];
real values[NTAB];
int kernik19rhs (real *u, real *du, Par par, Var var, int ln)
{
  STR *S = (STR *)par;
  int ivar;
  assert(ln==NV);
  if (var.n) for(ivar=0;ivar<var.n;ivar++) *(var.dst[ivar])=u[var.src[ivar]];
  real V;			/* voltage used in tabulated gates definitions */
  int in, it, iv;		/* vector components counters */

  /* printf("#### t=%ld",t); for (iv=0;iv<NV;iv++) printf("\t%g",u[iv]); printf("\n"); */
  
  V=u[V_index];
  /* In: V, ntab */
  /* Out: values */
  if (!ftab_kernik19(V,values,NTAB)) {
    ABORT("\nerror calculating ftab(%s) at t=%ld: V=%g\n","kernik19",t,V);
  } /*  if !ftab... */

  /* In: u, nv, ntab, nn, no, values, par, var */
  /* Out: du;  nalp, nbet - potentially, not in this model */
  if (!fddt_kernik19(u,NV,values,NTAB,par,var,du,NO,nalp,nbet,NN)) {
    URGENT_MESSAGE("\nerror calculating fddt(%s) at t=%ld: u=","kernik19",t);
    for(iv=0;iv<NV;iv++) URGENT_MESSAGE(" %lg",u[iv]);
    ABORT("\n");
  } /*  if !faddy... */

  /* values::  0..NT-1: alp; NT..2*NT-1: bet; 2*NT..NTAB-1: fun */
  /* u, du:: 0..NO-1: other, NO..NO+NG-1: ngate, NO+NG+NT-1=NV-1: tgate */
  /* ngate=u+no; */
  /* tgate=u+no+nn; */
  /* markov=u+no+nn+nt; */
  
  /* nontab gates */
  if (NN>0) for (in=0;in<NN;in++) {
    real a=nalp[in];
    real b=nbet[in];
    du[NO+in]=a-(a+b)*u[NO+in];
  }

  /* tab gates */
  if (NT>0) for (it=0;it<NT;it++) {
    real a=values[it];
    real b=values[NT+it];
    du[NO+NN+it]=a-(a+b)*u[NO+NN+it];
  }
  printf("================= t=%ld a[7]=%g b[7]=%g\n",t,values[7],values[NT+7]);
  return 1;
}


int create_kernik19rhs (Par *par, Var *var, char *w, real **u, int v0)
{
  /* printf("xxxxxxxxxxxxxxxx create_kernik19rhs entered\n"); */
  STR *S = (STR *)Calloc(1,sizeof(STR));
  char *ptr=w;
  int ivar=0;
  if (!S) ABORT("cannot create %s","kernik19rhs");
  for(var->n=0;*ptr;var->n+=(*(ptr++)==AT));
  if(var->n){CALLOC(var->dst,var->n,sizeof(real *));
  CALLOC(var->src,var->n,sizeof(int));}
  else{var->src=NULL;(var->dst)=NULL;}
  
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "kernik19_par.h"
  #undef _
  ACCEPTP(IV,0,RNONE,RNONE);

  CALLOC(*u,NV,sizeof(real));

  /* Assign the initial values as given the following *.h files */
  #define _(name,initial) (*u)[var_##name]=initial;
  #include "kernik19_other.h"
  #include "kernik19_ngate.h"
  #include "kernik19_tgate.h"
  #undef _
  
  {int i; printf("xxxx t=%ld",t); for (i=0;i<NV;i++) printf("\t%g",(*u)[i]); printf("\n");}
  var->n=ivar;
  if(ivar){REALLOC(var->dst,1L*ivar*sizeof(real *));
  REALLOC(var->src,1L*ivar*sizeof(int));}
  else{FREE(var->dst);FREE(var->src);}
  *par = S;
  return NV;
}

