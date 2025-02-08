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

/* Hodkin & Huxley model (J Physiol 117:500-544,1952) */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

/* flag AUTHENTIC to switch between original Hodgkin-Huxley
   1952 and a new model */
/* 
   the new model has addiabaticaly exluded m-gate from the 
   vector of dynamical variables; and the convention of 
   measuring the voltage is reversed.
*/
/* AUTHENTIC 1 for authors model */
/* AUTHENTIC 0 for new implementation */
#define AUTHENTIC 1

/* number of layers of dynamical variables */
#if AUTHENTIC
#define N 4
#else
#define N 3
#endif /* AUTHENTIC */

#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

static double cub(double x) {return x*x*x;}
static double qrt(double x) {double q=x*x; return q*q;}
#if (AUTHENTIC == 0)
/* Goldman(?) - H&H exponent */
static double hhx(double x){ 
  if(-1.e-4<x && x<1.e-4)
    return 1+x/2.+x*x/12.; 
  else 
    return x/(1-exp(-x));
}
#endif /* AUTHENTIC */

typedef struct {
  #define _(n,v) real n;
  #include "hh.h"
  #undef _
  real I; /* current/area = mV/ms*capacitance-trmbr current */
  real IV;/* mV/ms - intercellular current */
} STR;

/* RHS_HEAD expands to full right hand side function for the model
   (including all the dynamical equations):

   int hh(real *u, real *du, Par par, Var var, int ln)

   INPUT ARGUMENTS
   ---------------
   real *u      | pointer to array of states variables
   real *du     | pointer to array of increments of *u
   Par par	| parameter structure
   Var var	| variable structure 
   int ln	| number of variables
*/
RHS_HEAD(hh,N) {
  /* DEVICE_CONST copies elements of device structure S */
  #define _(n,v) DEVICE_CONST(real,n)
  #include "hh.h"
  #undef _
  DEVICE_CONST(real,I)
  DEVICE_CONST(real,IV)
  real V = u[0];
  real h = u[1];
  real n = u[2];
  real *dV = du+0;
  real *dh = du+1;
  real *dn = du+2;
  real INa, IK, Il;
#if AUTHENTIC
  real m = u[3];
  real *dm = du +3;

  /* define macros for identical names as in ionic model */
  /* dynamical variables */
  # define V_m V
  /* transition rates */
  #define alph alpha_h
  #define beth beta_h
  #define alpn alpha_n
  #define betn beta_n

  real alpha_m = 0.1 * (-V_m + 25.0) /
    (exp ((-V_m + 25.0) / 10.0) - 1.0);
  real beta_m = 4.0 * exp (-V_m / 18.0);
  real alpha_h = 0.07 * exp (-V_m / 20.0);
  real beta_h = 1.0 / (exp ((-V_m + 30.0) / 10.0) + 1.0);
  real alpha_n = 0.01 * (-V_m + 10.0) /
    (exp ((-V_m + 10.0) / 10.0) - 1.0);
  real beta_n = 0.125 * exp (-V_m / 80.0);

  /* this is corresponding dynamical equation for m gate */
  *dm = alpha_m*(1-m)-beta_m*m;
#else
  real m;
  #define alpm (km*hhx(-(V+25.)/10.))
  #define betm (km*4.*exp(V/18.))
  #define alph (kh*0.07*exp(V/20.))
  #define beth (kh*1./(exp((V+30)/10)+1))
  #define alpn (kn*0.1*hhx(-(V+10.)/10.))
  #define betn (kn*0.125*exp(V/80.))
  m=alpm/(alpm+betm);
#endif /* AUTHENTIC */

  IK = gK*qrt(n)*(V-VK);
  INa = gNa*cub(m)*h*(V-VNa);
  Il = gl*(V-Vl);

  *dV = -1./C*( IK + INa + Il + I ) + IV;
  *dh = alph*(1-h)-beth*h;
  *dn = alpn*(1-n)-betn*n;
} RHS_TAIL(hh)

/* RHS_CREATE_HEAD expands to intitialise the model. This 
   includes assigning values for all model parameters, 
   i.e. reading from the script or keeping the default
   values otherwise; and creating the vector of initial
   (steady-state) values of dynamic variables, which will
   be used by time-stepping device to initialise the whole
   medium.

   int 
   create_hh(Par *par, Var *var, char *w, real **u, int v0)

   INPUT ARGUMENTS
   ---------------
   Par par	| parameter structure
   Var var	| variable structure 
   char *w      | parameters to be assigned from script
   real **u     | pointer to array of states variables
   int v0       | number of entries in states array
*/
RHS_CREATE_HEAD(hh) {
  /* ACCEPTP reads value of named parameter from BBS script */
  #define _(n,v) ACCEPTP(n,v,RNONE,RNONE);
  #include "hh.h"
  #undef _
  ACCEPTP(I,0,RNONE,RNONE);
  ACCEPTP(IV,0,RNONE,RNONE);
  MALLOC(*u,N*sizeof(real));
#if AUTHENTIC
/* V_m  - membrane potential (mV) */
  (*u)[0] = 7.0;
/* h - inactivation gate of sodium current I_Na  */
  (*u)[1] =0.5960;
/* n - activation gate of potassium current I_K  */
  (*u)[2] =0.3177;
/* m - activation gate of sodium current I_Na  */
  (*u)[3] =0.0530;
#else
  (*u)[0]=0;
  (*u)[0]=1;
  (*u)[0]=0;
#endif /* AUTHENTIC */
} RHS_CREATE_TAIL(hh,N)
