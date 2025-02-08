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
Human ventricle cell model.
This model was published in "Alternans and spiral breakup in a human
ventricular tissue model", by K. H. W. J. ten Tusscher and A. V. Panfilov, 
J Physiol Heart Circ Physiol, 291: pp. H1088–H1100, 2006.

The source code was obtained from CELLML:
http:/www.cellml.org/
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

/*
Do everything here.
*/

/*
this is standard format for Beatbox variables declarations. It is a enum
of all variables. Note that in previously implemented cell models, it is
an explicit list of variables, while in this cell model, it is an array
of size 41. However, this array of Y (from CellML) is going to be assigned
to the Beatbox array u. So no use of keeping it here.
*/

/* Enumerate the dynamic ("state") variables within the vector of dynamic variables */
/*
enum {
#define _(n,i) ord_##n,
#include "ord_var.h"
#undef _
  N	*/	    		/* total number of variables */
// };

#define N 19 // number of variables.

/* The structure containing the parameter values for this instance of the model */
/* The constants in CellML are declared as extern double xxxx = dddd; and so that could potentially be another thing to alter in Beatbox format */
typedef struct {                	/* First go the cell parameters as list in the cellML to COR C/C++ output */
  #define _(name,default,unit) real name;
  #include "ttnp_par.h"
  #undef _
  real IV; 			/* The external current. */
  real ki_IV; /* stimulus also goes into Ki */
  real ht; /* the time step -The basal dt is 0.001 ms according to COR/CellML */
} STR;

/*
 **********************************************************************************
 * The function initialising an instance of the model.                            *
 * This includes: assigning values for all model parameters,                      *
 * i.e. reading from the script or keeping the default values otherwise;      *
 * and creating the vector of initial (steady-state) values of dynamic variables, *
 * which will be used by euler to initialise the whole medium.                    *
This part of the code simply involves replacing ord by whatever new model you are
implementing. It assumes that your model can be solved using the Euler method.
Correction: If the vector u and the derivative du are coming out of the RHS of the
ODEs, then any solver can be used. As of now, we have the Euler device solver,
but new solvers can be implemented.
 **********************************************************************************
 */
RHS_CREATE_HEAD(ttnp) {
  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default,unit) ACCEPTP(name,default,0,RNONE);
  #include "ttnp_par.h"
  #undef _

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */
  ACCEPTP(ht,RNONE,RSUCC(0),RNONE); /* default step is 0.001 ms */
  ACCEPTP(ki_IV,0,RNONE,RNONE);
  /*  but we need to ensure this value is the same as used in all other modules,  */
  /* so don't provide any default here but require explicit specification in the script. */


  /* Create the vector of initial (steady-state) values of dynamic (state) variables. */
  MALLOC(*u,N*sizeof(real));
  /* Cellml code assumes variables are called Y[n], n=0..40 */
  real *Y=*u;
  /* Assign the initial values: piece of code from Cellml  */

  /* Assign the initial values as given in ord_init.h */
  #include "ttnp_init.h"

} RHS_CREATE_TAIL(ttnp,N)

/*
 *******************************************************************************
 * The function "computing the right-hand sides".                              *
 * Actually, simply recalculates the value of the variables for the next step. *
 * Made as an rhs module to take advantage of parameter substitution mechanism.*
 *******************************************************************************
 */
RHS_HEAD(ttnp,N) {
  /* Declare the constant parameters and take their values from structure S==par (a formal parameter) */
  #define _(name,default,unit) DEVICE_CONST(real,name)
  #include "ttnp_par.h"
  #undef _
  DEVICE_CONST(real,IV)
  DEVICE_CONST(real,ki_IV);
  DEVICE_CONST(real,ht)
  /* Declare the dynamic variables and give them values from the vector u[] (a formal parameter) */

real *Y=u; // note how 
real *dY=du; // this is the important step.

/* Here go all actual computations for 1 time step */
#include "ttnp_step.h"

  /* Add the external current's effect. This assumes that the IV is already normalised/multiplied by capacitance. */
  u[11] += ht*IV;
  u[12] += ht*ki_IV;

  /* Copy the updated values of dynamic variables back to the vector u[].  */
  /* Assign the derivatives to zero as the values are already fully updated. Not true any more. */
  /* maybe for analyses of cell model properties and the solver, you want du from the cell model calculations */
/* In the ORD/CellML setting, this part is not required any more? */

} RHS_TAIL(ttnp)
