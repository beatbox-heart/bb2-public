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

/* TenTusscher-Panfilov 2006 human ventricular model.  */
/* Based on TenTusscher's own code from http://www-binf.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/ */
/* retrieved 2013/02/08 */
/* This model was published in "Alternans and spiral breakup in a human ventricular tissue model",
   by K. H. W. J. ten Tusscher and A. V. Panfilov, 
   J Physiol Heart Circ Physiol, 291: pp. H1088 –H1100, 2006 
*/

/* Changes made: the forward Euler updates replaced with ODE descriptions,  */
/* and Rush-Larsen updates replaced with alpha-beta descriptions */
/* when those depend only on V, and ODE descriptions otherwise.  */
/* NB all essential contents are stoved away into ttp2006_*.h files.  */

/* the calcium concentration is now calculated based on
   dynamical formulas for all calcium, rather than free and buffered
   calcium as it was in author's code
 */

/* MAKE THE CODE MRL COMPATIBLE */
/* Turn on Markov chain for mrl -- requires dummy macros  as shown bellow */
/* macro to switch on MRL specific parts (mainly in imported files) */
/* Please, also modify the ioniclist.h accordingly. */

/* To do: make this into a switch controlled from bbs cript */
#define EPI
/* #define ENDO */
/* #define MCELL */


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


/* dummy Macros required by mrl device */
/* number of Markov chain models */
#define NMC 0

/* In this model, transmembrane voltage is the first element of state vector */
static int V_index=0;

/* Define constants (unchangeable parameters) */
/* that can be used in tabulated functions. */
#include "ttp2006_const.h"

/* Enumerate the dynamic ("state") variables within the vector of dynamic variables */
enum {
  #define _(n,i) var_##n,
  #include "ttp2006_other.h"
  #include "ttp2006_ngate.h"
  #include "ttp2006_tgate.h"
  #undef _
  NV		    		/* total number of variables */
};

/* Enumerate the other (non-gate) variables */
enum {
  #define _(n,i) other_##n,
  #include "ttp2006_other.h"
  #undef _
  NO		    		/* total number of other variables */
};

/* Enumerate the non-tabulated gates */
enum {
  #define _(n,i) ngate_##n,
  #include "ttp2006_ngate.h"
  #undef _
  NN		    		/* total number of non-tabulated gate variables */
};

/* Enumerate the tabulated gates */
enum {
  #define _(n,i) tgate_##n,
  #include "ttp2006_tgate.h"
  #undef _
  NT		    		/* total number of tabulated gate variables */
};

/* Enumerate the tabulated transition rates and other tabulated functions */
enum {
  #define _(n,i) _alp_##n,
  #include "ttp2006_tgate.h"
  #undef _
  #define _(n,i) _bet_##n,
  #include "ttp2006_tgate.h"
  #undef _
  #define _(n) _##n,
  #include "ttp2006_fun.h"
  #undef _
  NTAB				/* total number of tabulated functions */
};

/* The structure containing the parameter values for this instance of the model */
typedef struct {                	/* First go the cell parameters as listed in K.tT's code */
  #define _(name,default) real name;
  #include "ttp2006_par.h"
  #undef _
  real IV; 			/* The external current. */
  /* 
  ! NB: IV is in millivolt/millisecond,								!
  ! and is positive outword raher than positive inward, as the cellular current densities,	!
  ! so e.g. IV can be = D \nabla^2 V, where D is the voltage diffusion coefficient !!		!
  */
  real ht;			/* the time step, needed for RL steps */
} STR;

  /* Define all the functions to be tabulated, including gates' transition rates, */
/* for given function of the transmebrane voltage V. */
IONIC_FTAB_HEAD(ttp2006) {
  #include "ttp2006_ftab.h"
  /* Copy the results into the output array values[]. */
  /* Care is taken that all, and only, tabulated functions are attended here. */
  #define _(n,i) values[_alp_##n]=alp_##n;
  #include "ttp2006_tgate.h"
  #undef _
  #define _(n,i) values[_bet_##n]=bet_##n;
  #include "ttp2006_tgate.h"
  #undef _
  #define _(n) values[_##n]=n;
  #include "ttp2006_fun.h"
  #undef _
} IONIC_FTAB_TAIL;

IONIC_FDDT_HEAD(ttp2006,NV,NTAB,NO,NN) {
  /* Declare the const pars and take their values from struct S==par (a formal parameter) */
  #define _(name,default) DEVICE_CONST(real,name);
  #include "ttp2006_par.h"
  #undef _
  DEVICE_CONST(real,IV);
  /* Declare and assign local variables for dynamic variables from state vector */
  /* ..., first for non-gate variables */
  #define _(name,initial) real name=u[var_##name];
  #include "ttp2006_other.h"
  #undef _
  /* ..., for the non-tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "ttp2006_ngate.h"
  #undef _
  /* ..., and then for tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "ttp2006_tgate.h"
  #undef _
  /* Get the values of the tabulable functions passed through the formal parameter */
  #define _(name) real name=values[_##name];
  #include "ttp2006_fun.h"
  #undef _
  /* Calculate the rates of non-gate variables */
  #include "ttp2006_fddt.h"
  /* Copy the calculated rates into the output array du[].  */
  /* Care is taken that all, and only, non-gating variables are attended here */
  #define _(name,initial) du[other_##name]=d##name;
  #include "ttp2006_other.h"
  #undef _
  /* And copy the non-tab transition rates into the output arraysn nalp[] and nbet[] */
  #define _(name,initial) nalp[ngate_##name]=alp_##name; nbet[ngate_##name]=bet_##name; 
  #include "ttp2006_ngate.h"
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
/*
 ********************************************************************************** 
 * NB the initial values of gate variables as defined in ttp2006_{n,t}gate.h    *
 * might be replaced by their steady-state values by the rushlarsen device	  *
 * (rushlarsen.c line 346-377).							  *
 **********************************************************************************
 */

IONIC_CREATE_HEAD(ttp2006) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "ttp2006_par.h"
  #undef _

  /* Assign the initial values as given the following *.h files */
  #define _(name,initial) (*u)[var_##name]=initial;
  #include "ttp2006_other.h"
  #include "ttp2006_ngate.h"
  #include "ttp2006_tgate.h"
  #undef _
} IONIC_CREATE_TAIL;
