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
  #define _(n,i) mousevent_##n,
  #include "mousevent_init.h"
  #undef _
  N		    		/* total number of variables */
};

typedef struct {                	/* First go the cell parameters as listed in M.C.'s code */
  #define _(name,default,unit) real name;
  #include "mousevent_par.h"
  #undef _
  real IV; 			/* The external current. */
  real ht;			/* the time step - will be equal to dt in the M.C.'s code */
} STR;

RHS_CREATE_HEAD(mousevent) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default,unit) ACCEPTP(name,default,-1000,RNONE); name=S->name;
  #include "mousevent_par.h"
  #undef _

  /* Now the values of the local variables are assigned to the structure S elements */
  #define _(name,default,unit) S->name=name;
  #include "mousevent_par.h"
  #undef _

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */
  ACCEPTP(ht,RNONE,RSUCC(0),RNONE);	/* M.C. suggests a 5 microseconds time step, */
  MALLOC(*u,N*sizeof(real));
  /* Assign the initial values as given in crn_var.h */
  #define _(name,initial) (*u)[mousevent_##name]=initial;
  #include "mousevent_init.h"
  #undef _
} RHS_CREATE_TAIL(mousevent,N)

/*
 *******************************************************************************
 */
RHS_HEAD(mousevent,N) {
  /* Declare the constant parameters and take their values from structure S==par (a formal parameter) */
  #define _(name,default,unit) DEVICE_CONST(real,name)
  #include "mousevent_par.h"
  #undef _
  DEVICE_CONST(real,IV)
  DEVICE_CONST(real,ht)
  /* Declare the dynamic variables and give them values from the vector u[] (a formal parameter) */
  #define _(name,initial) real name=u[mousevent_##name];
  #include "mousevent_init.h"
  #undef _
  /* The time step is called ht everywhere in Beatbox and dt in M.C.'s code */
  real dt=ht;

  /* Here go all actual computations of Bondarenko paper code */
  #include "mousevent_step.h"

  /* Add the external current's effect */
  vm += ht*IV;
  ki += ht*IV*acap*cm/(vmyo*F);

  /* Copy the updated values of dynamic variables back to the vector u[].  */
  /* Assign the derivatives to zero as the values are already fully updated */
  #define _(name,initial) u[mousevent_##name]=name; du[mousevent_##name]=0;
  #include "mousevent_init.h"
  #undef _
} RHS_TAIL(mousevent)
