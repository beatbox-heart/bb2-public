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
  Human atrial cell model as described in: Courtemanche M, Ramirez RJ,
  Nattel S. Ionic mechanisms underlying human atrial action potential
  properties: insights from a mathematical model. Am J Physiol
  (1998);275:H301-H321.  This code is based on the M.Courtemanche's
  `error-free' electronic version of equations and parameters. See the
  copyright notice in crn_step.h for detail.
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

/* Enumerate the dynamic ("state") variables within the vector of dynamic variables */
enum {
  #define _(n,i) crn_##n,
  #include "crn_var.h"
  #undef _
  N		    		/* total number of variables */
};


/* The structure containing the parameter values for this instance of the model */
typedef struct {                	/* First go the cell parameters as listed in M.C.'s code */
  #define _(name,default,unit) real name;
  #include "crn_par.h"
  #undef _
  real IV; 			/* The external current. */
  /* 
  ! NB: IV is in millivolt/millisecond rather than picoampere(/picofarad??),			!
  ! and is positive outword raher than positive inward, as the cellular current densities,	!
  ! so e.g. IV can be = D \nabla^2 V, where D is the voltage diffusion coefficient !!		!
  */
  real ht;			/* the time step - will be equal to dt in the M.C.'s code */
} rhspar;


/*
 **********************************************************************************
 * The function initialising an instance of the model.                            *
 * This includes: assigning values for all model parameters,                      *
 * i.e. reading from the script or keeping the default values otherwise;      *
 * and creating the vector of initial (steady-state) values of dynamic variables, *
 * which will be used by euler to initialise the whole medium.                    *
 **********************************************************************************
 */
RHS_CREATE_HEAD(crn) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default,unit) ACCEPTP(name,default,0,RNONE);
  #include "crn_par.h"
  #undef _

  /* Now the values of the local variables are assigned to the structure S elements */
  #define _(name,default,unit) S->name=name;
  #include "crn_par.h"
  #undef _

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */
  ACCEPTP(ht,RNONE,RSUCC(0),RNONE);	/* M.C. suggests a 5 microseconds time step, */
  /*  but we need to ensure this value is the same as used in all other modules,  */
  /* so don't provide any default here but require explicit specification in the script. */

  /* Create the vector of initial (steady-state) values of dynamic (state) variables. */
  MALLOC(*u,N*sizeof(real));
  /* Assign the initial values as given in crn_var.h */
  #define _(name,initial) (*u)[crn_##name]=initial;
  #include "crn_var.h"
  #undef _
} RHS_CREATE_TAIL(crn,N);


/*
 *******************************************************************************
 * The function "computing the right-hand sides".                              *
 * Actually, simply recalculates the value of the variables for the next step. *
 * Made as an rhs module to take advantage of parameter substitution mechanism.*
 *******************************************************************************
 */
RHS_HEAD(crn,N) {
  /* Declare the constant parameters and take their values from structure S==par (a formal parameter) */
  #define _(name,default,unit) DEVICE_CONST(real,name)
  #include "crn_par.h"
  #undef _
  DEVICE_CONST(real,IV)
  DEVICE_CONST(real,ht)
  /* Declare the dynamic variables and give them values from the vector u[] (a formal parameter) */
  #define _(name,initial) real name=u[crn_##name];
  #include "crn_var.h"
  #undef _
  /* The time step is called ht everywhere in Beatbox and dt in M.C.'s code */
  real dt=ht;

  /* Here go all actual computations of M.C.'s code */
  #include "crn_step.h"

  /* Add the external current's effect */
  V += ht*IV;

  /* Copy the updated values of dynamic variables back to the vector u[].  */
  /* Assign the derivatives to zero as the values are already fully updated */
  #define _(name,initial) u[crn_##name]=name; du[crn_##name]=0;
  #include "crn_var.h"
  #undef _
} RHS_TAIL(crn)
