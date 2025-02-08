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

/* Hopf normal form - reaction part of the CGLE */
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
    alp,bet,Omg;
  real 
    Iu,  Iv; 			/* Laplacians not just currents! */
} rhspar;

RHS_HEAD(cgle,2) {
  DEVICE_CONST(real,alp)
  DEVICE_CONST(real,bet)
  DEVICE_CONST(real,Omg)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)

  real mod2, U, V;

  U = u[0];
  V = u[1];
  mod2 = U*U+V*V;

  du[0] = U+Omg*V-(U+alp*V)*mod2 + Iu-bet*Iv;
  du[1] = V-Omg*U-(V-alp*U)*mod2 + Iv+bet*Iu;

} RHS_TAIL(cgle)

RHS_CREATE_HEAD(cgle)
  ACCEPTP(alp,RNONE,RNONE,RNONE);
  ACCEPTP(bet,RNONE,RNONE,RNONE);
  ACCEPTP(Omg,RNONE,RNONE,RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
  ACCEPTP(Iv,0,RNONE,RNONE);
RHS_CREATE_TAIL(cgle,2)
