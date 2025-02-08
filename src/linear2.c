/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
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

/* Linear two-component kinetics */
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
  real a11, a12, a21, a22, b1, b2; /* coeffs on linear ODE */
} STR;

RHS_HEAD(linear2,2)
{
  DEVICE_CONST(real,a11);
  DEVICE_CONST(real,a12);
  DEVICE_CONST(real,a21);
  DEVICE_CONST(real,a22);
  DEVICE_CONST(real,b1);
  DEVICE_CONST(real,b2);
  du[0] = a11*u[0]+a12*u[1]+b1;
  du[1] = a21*u[0]+a22*u[1]+b2;
}
RHS_TAIL(linear2)

RHS_CREATE_HEAD(linear2)
{
  ACCEPTP(a11,RNONE,RNONE,RNONE);
  ACCEPTP(a12,RNONE,RNONE,RNONE);
  ACCEPTP(a21,RNONE,RNONE,RNONE);
  ACCEPTP(a22,RNONE,RNONE,RNONE);
  ACCEPTP(b1,RNONE,RNONE,RNONE);
  ACCEPTP(b2,RNONE,RNONE,RNONE);
}
RHS_CREATE_TAIL(linear2,2)
