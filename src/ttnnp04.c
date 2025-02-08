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

/* TenTusscher-Noble-Noble-Panfilov 2004 human ventricular model.
   The model is published as: "A model for human ventricular tissue", by K. H. W. J. ten Tusscher, D. Noble, P. J. Noble and A. V. Panfilov,
   Am J Physiol Heart Circ Physiol, 286: pp. H1573-H1589, 2004,
   The code was obtained from TenTusscher's own code from http://www-binf.bio.uu.nl/khwjtuss/SourceCodes/HVM/Source/
   as of retrieved 2013/02/02 */

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
  #define _(n,i) ttnnp04_##n,
  #include "ttnnp04_var.h"
  #undef _
  N		    		/* total number of variables */
};


/* The structure containing the parameter values for this instance of the model */
typedef struct {                	/* First go the cell parameters as listed in K.tT's code */
  #define _(name,default) real name;
  #include "ttnnp04_par.h"
  #undef _
  real IV; 			/* The external current. */
  /* 
  ! NB: IV is in millivolt/millisecond,								!
  ! and is positive outword raher than positive inward, as the cellular current densities,	!
  ! so e.g. IV can be = D \nabla^2 V, where D is the voltage diffusion coefficient !!		!
  */
  real ht;			/* the time step, needed for RL steps */
  char variant[16];		/* epi/mcell/endo : string read */
  int model; 			/* the variant's number 0..2 */
} STR;

enum {
  model_EPI,
  model_ENDO,
  model_MCELL
};


/*
 **********************************************************************************
 * The function initialising an instance of the model.                            *
 * This includes: assigning values for all model parameters,                      *
 * i.e. reading from the script or keeping the default values otherwise;      *
 * and creating the vector of initial (steady-state) values of dynamic variables, *
 * which will be used by euler to initialise the whole medium.                    *
 **********************************************************************************
 */
RHS_CREATE_HEAD(ttnnp04)
{
  /* Model should be selected first as the parameter defaults depend on that */
  int model; /* 0..2: the variant of the model */
  ACCEPTS(variant, "EPI");
  STRSWITCH(S->variant);
  STRCASE("EPI")   model=model_EPI;
  STRCASE("epi")   model=model_EPI;
  STRCASE("ENDO")  model=model_ENDO;
  STRCASE("endo")  model=model_ENDO;
  STRCASE("MCELL") model=model_MCELL;
  STRCASE("mcell") model=model_MCELL;
  STRDEFAULT	MESSAGE("/* unknown model; using EPI */"); model=model_EPI;
  STRENDSW;
  S->model=model;
  
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #define _ver(a1,a2,a3) (model==model_EPI?(a1):model==model_ENDO?(a2):model==model_MCELL?(a3):RNONE)
  #include "ttnnp04_par.h"
  #undef _

  /* Now the values of the local variables are assigned to the structure S elements */
  /*   - this is obsolete with the current definition pf ACCEPTP */
  /* #define _(name,default) S->name=name; */
  /* #include "ttnnp04_par.h" */
  /* #undef _ */

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */
  ACCEPTP(ht,RNONE,RSUCC(0),RNONE); 	/* CELLML recommended step is 0.001 ms */
  /*  but we need to ensure this value is the same as used in all other modules,  */
  /* so don't provide any default here but require explicit specification in the script. */

  /* Create the vector of initial (steady-state) values of dynamic (state) variables. */
  MALLOC(*u,N*sizeof(real));

  /* Assign the initial values: piece of code from Cellml  */
  #define _(name,initial) (*u)[ttnnp04_##name]=initial;
  #include "ttnnp04_var.h"
  #undef _

}
RHS_CREATE_TAIL(ttnnp04,N)

/*
 *******************************************************************************
 * The function "computing the right-hand sides".                              *
 * Actually, simply recalculates the value of the variables for the next step. *
 * Made as an rhs module to take advantage of parameter substitution mechanism.*
 *******************************************************************************
 */
RHS_HEAD(ttnnp04,N)
{
  /* Declare the constant parameters and take their values from structure S==par (a formal parameter) */
#define _(name,default) DEVICE_CONST(real,name)
#include "ttnnp04_par.h"
#undef _
  DEVICE_CONST(real,IV);
  DEVICE_CONST(real,ht);
  DEVICE_CONST(int, model);
  
#define       HT          ht

  /* Here go all actual computations for 1 time step */
#include "ttnnp04_step.h"
  (u[ttnnp04_Volt]) += ht*IV;

  /* Assign the derivatives to zero as the values are already fully updated */
  #define _(name,initial) du[ttnnp04_##name]=0;
  #include "ttnnp04_var.h"
  #undef _

} RHS_TAIL(ttnnp04)
