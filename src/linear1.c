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

/* Linear one-component kinetics */
#include <assert.h>
#include <math.h>
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
  real a, b; /* coeffs on linear ODE */
} STR;

RHS_HEAD(linear1,1)
{
  DEVICE_CONST(real,a);
  DEVICE_CONST(real,b);
  du[0] = a*u[0]+b;
  if (!isfinite(a) || !isfinite(b) || !isfinite(u[0]) || !isfinite(du[0])) 
    printf("%s:%d t=%ld idev=%d a=%g b=%g u=%g du=%g\n",__FILE__,__LINE__,t,idev,a,b,u[0],du[0]);
}
RHS_TAIL(linear1)

RHS_CREATE_HEAD(linear1)
{
  ACCEPTP(a,RNONE,RNONE,RNONE);
  ACCEPTP(b,RNONE,RNONE,RNONE);
}
RHS_CREATE_TAIL(linear1,1)
