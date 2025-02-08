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

/* Caricature model of the INa-driven front, described in      */
/* V.N. Biktashev, "Dissipation of the excitation wavefronts", */
/* Phys. Rev. Lett., 89(16): 168102, 2002                      */

#include <assert.h>
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

typedef struct {
  real Em, Eh, tau;
  real IE,  Ih;
} rhspar;

#define theta(x) ((x)>0?1.0:0.0)


RHS_HEAD(b02,2) {
  DEVICE_CONST(real,Em);
  DEVICE_CONST(real,Eh);
  DEVICE_CONST(real,tau);
  DEVICE_CONST(real,IE);
  DEVICE_CONST(real,Ih);

  real E = u[0];
  real h = u[1];

  du[0] = theta(E-Em)*h + IE;
  du[1] = (theta(Eh-E)-h)/tau + Ih;

} RHS_TAIL(b02);

RHS_CREATE_HEAD(b02) {
  ACCEPTP(Em,1.0,RNONE,RNONE);
  ACCEPTP(Eh,0.0,RNONE,RNONE);
  ACCEPTP(tau,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(IE,0.0,RNONE,RNONE);
  ACCEPTP(Ih,0.0,RNONE,RNONE);
} RHS_CREATE_TAIL(b02,2);

