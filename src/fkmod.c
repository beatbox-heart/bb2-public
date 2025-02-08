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
  Modified Fenton & Karma model of cardiac AP.
  This code is based on formulas given in the paper:
  F Fenton and A Karma, "Vortex dynamics in three-dimensional 
  continuous myocardium with fiber rotation: Filament instability and 
  fibrillation". Chaos, 8(1): pp. 20-47, 1998.
  The cell model code was developed using Beatbox by modifying the fk cell model module.
 */

/* 
SK January 2013. Modified FK with Ih current.
This is the normalised system where u0 (the membrane potential)
is set up to go from 0 to 1. To transform to biophysical units,
invert this to get V if required:
u = (V - Vo)/(Vfi - Vo)
where Vo is resting
vfi is INa nernest potential, around +40 mV
The currents are obtained by inverting:
J = I/(Cm*(Vfi - Vo))
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

/* Temporally: to mark undefined parameters */
#define UNKNOWN 0

#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

typedef struct {
  char basic[80];
  #define _(n,v1,v2,v3,v4) real n;
  #include "fkmod.h"
  #undef _
  real Iu;
} STR;

RHS_HEAD(fkmod,4) {
  #define _(n,v1,v2,v3,v4) DEVICE_CONST(real,n)
  #include "fkmod.h"
  #undef _
  DEVICE_CONST(real,Iu)
  real U = u[0];
  real V = u[1];
  real W = u[2];
  real Y = u[3];
  real Jfi, Jso, Jsi, tau_vm, tauh;
  real Jh, yinf;

  if (U<u_c) {
    Jfi=0;
    Jso = U/tau_0;
    tau_vm=(U>u_v)?(tau_vm1):(tau_vm2);
    du[1] = (1.-V)/tau_vm;
    du[2] = (1.-W)/tau_wm;
  } else {
    Jfi = - V/tau_d*(1.-U)*(U-u_c);
    Jso = 1/tau_r;
    du[1] = (-V)/tau_vp;
    du[2] = (-W)/tau_wp;
  }
  Jsi = -gsi*W/tau_si*(1+tanh(k*(U-usi_c)));

/*
The hyperpolarisation activated current without kinetics.
u_c is low, so that is the reversal potential of the outward current.
1 or usi_c is reversal potential of the sodium part.
Uh is usually way off.
the time constant for Ih is not small, and cannot be a constant number either.
If I take no time kinetics for Ih, there is no influence of Ih on the AP.
If I take it to be a constant, then it gives a lopsided AP which is not
the desired share.
So since my APD for atrial is now ~ 20 t.u. (same as mouse atrial APD of 20 ms),
I will put the time constant of Ih to be the same as in my mouse SAN model.
*/

yinf = 1.0/(1.0+exp((U+Uh)/kh)); // U-(-Uh) = U+Uh, and then |Uh| >> |Umin|. In mouse, |Uh| = 106 and |Umin| = 63
tauh = p5/(exp(-(U+p1)*p2)+ exp((U-p3)/p4)); // all parameters taken to give lopsided bell shape with paeak of 500 t.u. close to resting.
du[3] = (yinf - Y)/tauh;
Jh = ifswitch*(ghna*Y*(U-1.0) + ghk*Y*(U-u_c));

  du[0] = - Jfi - Jso - Jsi + Iu - Jh;
} RHS_TAIL(fk)

RHS_CREATE_HEAD(fkmod) {
  int model;
  #define Cm 1.0
  #define _(n,v1,v2,v3,v4) real std_##n[4]={v1,v2,v3,v4};
  #include "fkmod.h"
  #undef _
  ACCEPTS(basic, "GP");
  STRSWITCH(S->basic);
    STRCASE("BR")   model=0;
    STRCASE("MBR")  model=1;
    STRCASE("MLR-I")model=2;
    STRCASE("GP")   model=3;
    STRDEFAULT	MESSAGE("/* unknown model; using BR */"); model=0;
  STRENDSW
  #define _(n,v1,v2,v3,v4) ACCEPTP(n,std_##n[model],0,RNONE);
  #include "fkmod.h"
  #undef _
  ACCEPTP(Iu ,   0,RNONE,RNONE);
  MALLOC(*u,4*sizeof(real));
  (*u)[0]=0;
  (*u)[1]=(*u)[2]=(*u)[3]=1;
} RHS_CREATE_TAIL(fkmod,4)
