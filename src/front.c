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

/* Caricature of the Na current */

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
  real ENa, EK, tau;
  real IE,  Ih;
} STR;

#define theta(x) ((x)>0?1.0:0.0)


RHS_HEAD(front,2) {

  DEVICE_CONST(real,ENa)
  DEVICE_CONST(real,EK)
  DEVICE_CONST(real,tau)
  DEVICE_CONST(real,IE)
  DEVICE_CONST(real,Ih)

  real E = u[0];
  real h = u[1];

  du[0] = theta(E-ENa)*h + IE;
  du[1] = (theta(EK-E)-h)/tau + Ih;

} RHS_TAIL(front)

RHS_CREATE_HEAD(front)
  ACCEPTP(ENa,RNONE,RNONE,RNONE);
  ACCEPTP(EK,RNONE,RNONE,RNONE);
  ACCEPTP(tau,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(IE,0,RNONE,RNONE);
  ACCEPTP(Ih,0,RNONE,RNONE);
RHS_CREATE_TAIL(front,2)
