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

/* Barkley's [Physica D 49:61-70, 1991] model */

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
  real
    a, b, eps;			/* Barkleys parameters */
  real 
    Iu,  Iv; 			/* External influence */
} STR;

RHS_HEAD(fhnbkl,2) {
  DEVICE_CONST(real,a)
  DEVICE_CONST(real,b)
  DEVICE_CONST(real,eps)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)
  real U, Uth, V;

  U = u[0];
  V = u[1];
  Uth = (V+b)/a;

  du[0] = U*(1.0-U)*(U-Uth)/eps + Iu;
  du[1] = U-V + Iv;

} RHS_TAIL(fhnbkl)

RHS_CREATE_HEAD(fhnbkl)
  ACCEPTP(a,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(b,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(eps,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
  ACCEPTP(Iv,0,RNONE,RNONE);
RHS_CREATE_TAIL(fhnbkl,2)
