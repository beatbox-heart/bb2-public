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
  IONIC description of the Courtemanche et al. 1998 model:
  Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying
  human atrial action potential properties: insights from a mathematical
  model. Am J Physiol (1998);275:H301-H321. 
  This code is based on the M.Courtemanche's `error-free' electronic 
  version of equations and parameters. 
  Changes made: the forward Euler updates replaced with ODE descriptions, 
  and Rush-Larsen updates replaced with alpha-beta descriptions
  when those depend only on V, and ODE descriptions otherwise. 
  NB all essential contents are stoved away into crn98_*.h files. 
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

/* number of Markov chain models */
#define NMC 0

static inline real sqr(real x) {return x*x;}
static inline real cub(real x) {return x*x*x;}

/* In this model, transmembrane voltage is the first element of state vector */
static int V_index=0;

/* Define constants (unchangeable parameters) */
/* that can be used in tabulated functions. */
#include "crn98_const.h"

/* Enumerate all dynamic variables */
enum {
  #define _(n,i) var_##n,
  #include "crn98_other.h"
  #include "crn98_ngate.h"
  #include "crn98_tgate.h"
  #undef _
  NV		    		/* total number of variables */
};

/* Enumerate the other (non-gate) variables */
enum {
  #define _(n,i) other_##n,
  #include "crn98_other.h"
  #undef _
  NO		    		/* total number of other variables */
};

/* Enumerate the non-tabulated gates */
enum {
  #define _(n,i) ngate_##n,
  #include "crn98_ngate.h"
  #undef _
  NN		    		/* total number of non-tabulated gate variables */
};

/* Enumerate the tabulated gates */
enum {
  #define _(n,i) tgate_##n,
  #include "crn98_tgate.h"
  #undef _
  NT		    		/* total number of tabulated gate variables */
};

/* Enumerate the tabulated transition rates and other tabulated functions */
enum {
  #define _(n,i) _alp_##n,
  #include "crn98_tgate.h"
  #undef _
  #define _(n,i) _bet_##n,
  #include "crn98_tgate.h"
  #undef _
  #define _(n) _##n,
  #include "crn98_fun.h"
  #undef _
  NTAB				/* total number of tabulated functions */
};

/* The structure containing the parameter values for this instance of the model */
typedef struct {		/* First go the cell parameters as listed in M.C.'s code */
  #define _(name,default) real name;
  #include "crn98_par.h"
  #undef _
  real IV; 			/* Then the external current. */
  /* 
  ! NB: IV is in millivolt/millisecond rather than picoampere(/picofarad??),			!
  ! and is positive outword raher than positive inward, as the cellular current densities,	!
  ! so e.g. IV can be = D \nabla^2 V, where D is the voltage diffusion coefficient !!		!
  */
} STR;

/* Define all the functions to be tabulated, including gates' transition rates, */
/* for given function of the transmebrane voltage V. */
IONIC_FTAB_HEAD(crn98) {
  #include "crn98_ftab.h"
  /* Copy the results into the output array values[]. */
  /* Care is taken that all, and only, tabulated functions are attended here. */
  #define _(n,i) values[_alp_##n]=alp_##n;
  #include "crn98_tgate.h"
  #undef _
  #define _(n,i) values[_bet_##n]=bet_##n;
  #include "crn98_tgate.h"
  #undef _
  #define _(n) values[_##n]=n;
  #include "crn98_fun.h"
  #undef _
} IONIC_FTAB_TAIL;	  

IONIC_FDDT_HEAD(crn98,NV,NTAB,NO,NN) {
  /* Declare the const pars and take their values from struct S==par (a formal parameter) */
  #define _(name,default) DEVICE_CONST(real,name);
  #include "crn98_par.h"
  #undef _
  DEVICE_CONST(real,IV);
  /* Declare and assign local variables for dynamic variables from state vector */
  /* ..., first for non-gate variables */
  #define _(name,initial) real name=u[var_##name];
  #include "crn98_other.h"
  #undef _
  /* ..., for the non-tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "crn98_ngate.h"
  #undef _
  /* ..., and then for tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "crn98_tgate.h"
  #undef _
  /* Get the values of the tabulable functions passed through the formal parameter */
  #define _(name) real name=values[_##name];
  #include "crn98_fun.h"
  #undef _
  /* Calculate the rates of non-gate variables */
  #include "crn98_fddt.h"
  /* Copy the calculated rates into the output array du[].  */
  /* Care is taken that all, and only, non-gating variables are attended here */
  #define _(name,initial) du[other_##name]=d_##name;
  #include "crn98_other.h"
  #undef _
  /* And copy the non-tab transition rates into the output arraysn nalp[] and nbet[] */
  #define _(name,initial) nalp[ngate_##name]=alp_##name; nbet[ngate_##name]=bet_##name; 
  #include "crn98_ngate.h"
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
IONIC_CREATE_HEAD(crn98) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "crn98_par.h"
  #undef _

  /* Assign the initial values as given the following *.h files */
  #define _(name,initial) (*u)[var_##name]=initial;
  #include "crn98_other.h"
  #include "crn98_ngate.h"
  #include "crn98_tgate.h"
  #undef _
} IONIC_CREATE_TAIL;


/* /\* */
/*  ******************************************************************************* */
/*  * The function "computing the right-hand sides".                              * */
/*  * Actually, simply recalculates the value of the variables for the next step. * */
/*  * Made as an rhs module to take advantage of parameter substitution mechanism.* */
/*  ******************************************************************************* */
/*  *\/ */
/* RHS_HEAD(crn,N) { */
/*   /\* Declare the constant parameters and take their values from structure S==par (a formal parameter) *\/ */
/*   #define _(name,default,unit) DEVICE_CONST(real,name) */
/*   #include "crn98_par.h" */
/*   #undef _ */
/*   DEVICE_CONST(real,IV) */
/*   DEVICE_CONST(real,ht) */
/*   /\* Declare the dynamic variables and give them values from the vector u[] (a formal parameter) *\/ */
/*   #define _(name,initial) real name=u[crn98_##name]; */
/*   #include "crn98_var.h" */
/*   #undef _ */
/*   /\* The time step is called ht everywhere in Beatbox and dt in M.C.'s code *\/ */
/*   real dt=ht; */

/*   /\* Here go all actual computations of M.C.'s code *\/ */
/*   #include "crn98_step.h" */

/*   /\* Add the external current's effect *\/ */
/*   V += ht*IV; */

/*   /\* Copy the updated values of dynamic variables back to the vector u[].  *\/ */
/*   /\* Assign the derivatives to zero as the values are already fully updated *\/ */
/*   #define _(name,initial) u[crn98_##name]=name; du[crn98_##name]=0; */
/*   #include "crn98_var.h" */
/*   #undef _ */
/* } RHS_TAIL(crn) */
