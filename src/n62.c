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

/* Noble-62 model of Purkinje cell (J Physiol 160:317-352, 1962)
 * m adiabatically excluded
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

/* Temporarily: to mark undefined parameters */
#define UNKNOWN 0

#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

static double cub(double x) {return x*x*x;}
static double qrt(double x) {double q=x*x; return q*q;}
/* Goldman(?) - H&H exponent */
static double hhx(double x){ if(-1.e-4<x && x<1.e-4) return 1+x/2.+x*x/12.; else return x/(1-exp(-x));}

typedef struct {
  #define _(n,v) real n;
  #include "n62.h"
  #undef _
  real I;
} STR;

RHS_HEAD(n62,3) {
  #define _(n,v) DEVICE_CONST(real,n)
  #include "n62.h"
  #undef _
  DEVICE_CONST(real,I)
  real E = u[0];
  real h = u[1];
  real n = u[2];
  real *dE = du+0;
  real *dh = du+1;
  real *dn = du+2;
  real m, INa, IK, IAn;

  #define alpm (1.5*hhx((E+48.)/15.))
  #define betm (0.6*hhx(-(E+8.)/5.))
  #define alph (0.17*exp(-(E+90.)/20.))
  #define beth (1./(exp(-(E+42.)/10.)+1))
  #define alpn (0.001*hhx((E+50.)/10.))
  #define betn (0.002*exp(-(E+90.)/80.))
  #define gK1  (1.2*exp(-(E+90.)/50.)+0.015*exp((E+90.)/60.))

  m=alpm/(alpm+betm);
  IK = (gK1 + gK*qrt(n))*(E-EK);
  INa = (gNa1 + gNa*cub(m)*h)*(E-ENa);
  IAn = gAn*(E-EAn);

  *dE = -1./C*( IK + INa + IAn + I );
  *dh = alph*(1-h)-beth*h;
  *dn = alpn*(1-n)-betn*n;

} RHS_TAIL(n62)

RHS_CREATE_HEAD(n62) {
  #define _(n,v) ACCEPTP(n,v,RNONE,RNONE);
  #include "n62.h"
  #undef _
  ACCEPTP(I,0,RNONE,RNONE);
  MALLOC(*u,3L*sizeof(real));
  (*u)[0]=S->EK;
  (*u)[1]=1;
  (*u)[2]=0;
} RHS_CREATE_TAIL(n62,3)
