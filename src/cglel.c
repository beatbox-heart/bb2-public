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

/* linearised Hopf normal form - reaction part of the CGLE */
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
    U, V;                       /* u- and v-component of the unperturbed solution */
  real 
    Iu,  Iv; 			/* Laplacians not just currents! */
} rhspar;

RHS_HEAD(cglel,2) {
  DEVICE_CONST(real,alp)
  DEVICE_CONST(real,bet)
  DEVICE_CONST(real,Omg)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)
  DEVICE_CONST(real,U)
  DEVICE_CONST(real,V)

  real u1 = u[0];
  real v1 = u[1];
  real mod2 = U*U+V*V;

  du[0] = (1-mod2-2*U*(U+alp*V))        *u1 + ( Omg-alp*mod2-2*V*(U+alp*V)) *v1 + (   1)*Iu + (-bet)*Iv;
  du[1] = (-Omg+alp*mod2-2*U*(V-alp*U)) *u1 + (1-mod2-2*V*(V-alp*U))        *v1 + ( bet)*Iu + (   1)*Iv;

} RHS_TAIL(cglel)

RHS_CREATE_HEAD(cglel)
  ACCEPTP(alp,RNONE,RNONE,RNONE);
  ACCEPTP(bet,RNONE,RNONE,RNONE);
  ACCEPTP(Omg,RNONE,RNONE,RNONE);
  ACCEPTP(U,RNONE,RNONE,RNONE);
  ACCEPTP(V,RNONE,RNONE,RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
  ACCEPTP(Iv,0,RNONE,RNONE);
RHS_CREATE_TAIL(cglel,2)
