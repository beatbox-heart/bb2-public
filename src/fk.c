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
  Generic AP cell model as described in:
  F Fenton and A Karma, "Vortex dynamics in three-dimensional 
  continuous myocardium with fiber rotation: Filament instability and 
  fibrillation". Chaos, 8(1): pp. 20-47, 1998.
  This code is based on formulas given in the above paper.
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
  #include "fk.h"
  #undef _
  real Iu;
} STR;

RHS_HEAD(fk,3)
  #define _(n,v1,v2,v3,v4) DEVICE_CONST(real,n)
  #include "fk.h"
  #undef _
  DEVICE_CONST(real,Iu)
  real U = u[0];
  real V = u[1];
  real W = u[2];
  real Jfi, Jso, Jsi, tau_vm;

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
  Jsi = -0.5*W/tau_si*(1+tanh(k*(U-usi_c)));
  du[0] = -Jfi - Jso - Jsi + Iu;
RHS_TAIL(fk)

RHS_CREATE_HEAD(fk) {
  int model;
  #define Cm 1.0
  #define _(n,v1,v2,v3,v4) real std_##n[4]={v1,v2,v3,v4};
  #include "fk.h"
  #undef _
  ACCEPTS(basic, "GP");
  STRSWITCH(S->basic);
    STRCASE("BR")   model=0;
    STRCASE("MBR")  model=1;
    STRCASE("MLR-I")model=2;
    STRCASE("GP")   model=3;
    STRDEFAULT	MESSAGE("/* unknown model; using BR */"); model=0;
  STRENDSW
  #define _(n,v1,v2,v3,v4) ACCEPTP(n,std_##n[model],RNONE,RNONE);
  #include "fk.h"
  #undef _
  ACCEPTP(Iu ,   0,RNONE,RNONE);
  MALLOC(*u,3*sizeof(real));
  (*u)[0]=0;
  (*u)[1]=(*u)[2]=1;
} RHS_CREATE_TAIL(fk,3)
