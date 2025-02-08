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
#include "rhs.h"
#include "qpp.h"

/* Enumerate the dynamic ("state") variables within the vector of dynamic variables */
enum {
  #define _(n,i) ttp06_##n,
  #include "ttp06_var.h"
  #undef _
  N		    		/* total number of variables */
};


/* The structure containing the parameter values for this instance of the model */
typedef struct {                	/* First go the cell parameters as listed in K.tT's code */
  #define _(name,default) real name;
  #include "ttp06_par.h"
  #undef _
  real IV; 			/* The external current. */
  /* 
  ! NB: IV is in millivolt/millisecond,								!
  ! and is positive outword raher than positive inward, as the cellular current densities,	!
  ! so e.g. IV can be = D \nabla^2 V, where D is the voltage diffusion coefficient !!		!
  */
  real ht;			/* the time step, needed for RL steps */
} STR;


/*
 **********************************************************************************
 * The function initialising an instance of the model.                            *
 * This includes: assigning values for all model parameters,                      *
 * i.e. reading from the script or keeping the default values otherwise;      *
 * and creating the vector of initial (steady-state) values of dynamic variables, *
 * which will be used by euler to initialise the whole medium.                    *
 **********************************************************************************
 */
RHS_CREATE_HEAD(ttp06) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "ttp06_par.h"
  #undef _

  /* Now the values of the local variables are assigned to the structure S elements */
  #define _(name,default) S->name=name;
  #include "ttp06_par.h"
  #undef _

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */
  ACCEPTP(ht,RNONE,RSUCC(0),RNONE); 	/* CELLML recommended step is 0.001 ms */
  /*  but we need to ensure this value is the same as used in all other modules,  */
  /* so don't provide any default here but require explicit specification in the script. */

  /* Create the vector of initial (steady-state) values of dynamic (state) variables. */
  MALLOC(*u,N*sizeof(real));

  /* Assign the initial values: piece of code from Cellml  */
  #define _(name,initial) (*u)[ttp06_##name]=initial;
  #include "ttp06_var.h"
  #undef _

} RHS_CREATE_TAIL(ttp06,N)

/*
 *******************************************************************************
 * The function "computing the right-hand sides".                              *
 * Actually, simply recalculates the value of the variables for the next step. *
 * Made as an rhs module to take advantage of parameter substitution mechanism.*
 *******************************************************************************
 */
RHS_HEAD(ttp06,N) {
  /* Declare the constant parameters and take their values from structure S==par (a formal parameter) */
  #define _(name,default) DEVICE_CONST(real,name)
  #include "ttp06_par.h"
  #undef _
  DEVICE_CONST(real,IV)
  DEVICE_CONST(real,ht)
    
  #define       HT          ht

  /* Here go all actual computations for 1 time step */
  #include "ttp06_step.h"
  (u[ttp06_Volt]) += ht*IV;

  /* Assign the derivatives to zero as the values are already fully updated */
  #define _(name,initial) du[ttp06_##name]=0;
  #include "ttp06_var.h"
  #undef _

} RHS_TAIL(ttp06)
