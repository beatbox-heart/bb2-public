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

/**
 * IONIC description of the Beeler-Reuter 1977 model.
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

/* In this model, transmembrane voltage is the first element of state vector */
static int V_index=0;

/* Define constants (unchangeable parameters) */
/* that can be used in tabulated functions. */
#include "br77_const.h"

/* Enumerate all dynamic variables */
enum {
  #define _(n,i) var_##n,
  #include "br77_other.h"
  #include "br77_ngate.h"
  #include "br77_tgate.h"
  #undef _
  NV		    		/* total number of variables */
};

/* Enumerate the other (non-gate) variables */
enum {
  #define _(n,i) other_##n,
  #include "br77_other.h"
  #undef _
  NO		    		/* total number of other variables */
};

/* Enumerate the non-tabulated gates */
enum {
  #define _(n,i) ngate_##n,
  #include "br77_ngate.h"
  #undef _
  NN		    		/* total number of non-tabulated gate variables */
};

/* Enumerate the gates */
enum {
  #define _(n,i) gate_##n,
  #include "br77_tgate.h"
  #undef _
  NT		    		/* total number of tabulated gate variables */
};

/* Enumerate the tabulated transition rates and other tabulated functions */
enum {
#define _(n,i) _alp_##n,
  #include "br77_tgate.h"
  #undef _
#define _(n,i) _bet_##n,
  #include "br77_tgate.h"
  #undef _
  #define _(n) _##n,
  #include "br77_fun.h"
  #undef _
  NTAB				/* total number of tabulated functions */
};

/* The structure containing the parameter values for this instance of the model */
typedef struct {                /* First go the canonical cell parameters */
  #define _(name,default) real name;
  #include "br77_par.h"
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
IONIC_FTAB_HEAD(br77) {
  #include "br77_ftab.h"
  /* Copy the results into the output array values[]. */
  /* Care is taken that all, and only, tabulated functions are attended here. */
  #define _(n,i) values[_alp_##n]=alp_##n;
  #include "br77_tgate.h"
  #undef _
  #define _(n,i) values[_bet_##n]=bet_##n;
  #include "br77_tgate.h"
  #undef _
  #define _(n) values[_##n]=n;
  #include "br77_fun.h"
  #undef _
} IONIC_FTAB_TAIL;	  

IONIC_FDDT_HEAD(br77,NV,NTAB,NO,NN) {
  /* Declare the const pars and take their values from struct S==par (a formal parameter) */
  #define _(name,default) DEVICE_CONST(real,name);
  #include "br77_par.h"
  #undef _
  DEVICE_CONST(real,IV);
  /* Declare and assign local variables for dynamic variables from state vector */
  /* ..., first for non-gate variables */
  #define _(name,initial) real name=u[var_##name];
  #include "br77_other.h"
  #undef _
  /* ..., for the non-tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "br77_ngate.h"
  #undef _
  /* ..., and then for tabulated gate variables */
  #define _(name,i) real name=u[var_##name];
  #include "br77_tgate.h"
  #undef _
  /* Get the values of the tabulable functions passed through the formal parameter */
  #define _(name) real name=values[_##name];
  #include "br77_fun.h"
  #undef _
  /* Calculate the rates of non-gate variables */
  #include "br77_fddt.h"
  /* Copy the calculated rates into the output array du[].  */
  /* Care is taken that all, and only, non-gating variables are attended here */
  #define _(name,initial) du[other_##name]=d_##name;
  #include "br77_other.h"
  #undef _
  /* And copy the non-tab transition rates into the output arrays nalp[] and nbet[] */
  #define _(name,initial) nalp[ngate_##name]=alp_##name; nbet[ngate_##name]=bet_##name; 
  #include "br77_ngate.h"
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
IONIC_CREATE_HEAD(br77) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "br77_par.h"
  #undef _
  
  /* Assign the initial values as given the following *.h files */
  #define _(name,initial) (*u)[var_##name]=initial;
  #include "br77_other.h"
  #include "br77_ngate.h"
  #include "br77_tgate.h"
  #undef _
} IONIC_CREATE_TAIL;

