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

/* the KYLZ mouse SAN model driver. */

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
  #define _(n,i) mouse_##n,
  #include "mouse_init.h"
  #undef _
  N		    		/* total number of variables */
};

/* The structure containing the parameter values for this instance of the model */
typedef struct {                	/* First go the cell parameters as listed in M.C.'s code */
  #define _(name,default,unit) real name;
  #include "mouse_par.h"
  #undef _
  real ht;			/* the time step as given by user */
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
RHS_CREATE_HEAD(mouse) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default,unit) ACCEPTP(name,default,-1000.0,RNONE); name=S->name;
  #include "mouse_par.h"
  #undef _

  /* Now the values of the local variables are assigned to the structure S elements */
  #define _(name,default,unit) S->name=name;
  #include "mouse_par.h"
  #undef _

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(ht,RNONE,RSUCC(0),RNONE);
  /*  but we need to ensure this value is the same as used in all other modules,  */
  /* so don't provide any default here but require explicit specification in the script. */

  /* Create the vector of initial (steady-state) values of dynamic (state) variables. */
  MALLOC(*u,N*sizeof(real));
  /* Assign the initial values as given in mouse_init.h */
  #define _(name,initial) (*u)[mouse_##name]=initial;
  #include "mouse_init.h"
  #undef _
} RHS_CREATE_TAIL(mouse,N)

/*
 *******************************************************************************
 * The function "computing the right-hand sides".                              *
 * Actually, simply recalculates the value of the variables for the next step. *
 * Made as an rhs module to take advantage of parameter substitution mechanism.*
 *******************************************************************************
 */
RHS_HEAD(mouse,N) {
  /* Declare the constant parameters and take their values from structure S==par (a formal parameter) */
  #define _(name,default,unit) DEVICE_CONST(real,name)
  #include "mouse_par.h"
  #undef _
  DEVICE_CONST(real,ht)
  /* Declare the dynamic variables and give them values from the vector u[] (a formal parameter) */
  #define _(name,initial) real name=u[mouse_##name];
  #include "mouse_init.h"
  #undef _
  /* The time step is called ht everywhere in Beatbox and dt in M.C.'s code */
  real ddt=ht;

  /* Here go all actuall computations of M.C.'s code */
  #include "mouse_step.h"

  /* Copy the updated values of dynamic variables back to the vector u[].  */
  /* Assign the derivatives to zero as the values are already fully updated */
  #define _(name,initial) u[mouse_##name]=name; du[mouse_##name]=0;
  #include "mouse_init.h"
  #undef _
} RHS_TAIL(mouse)
