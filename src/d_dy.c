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

/* Central Y derivative */

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
#include "qpp.h"

typedef struct {
  int v0, v1;		/* Range of layers */
  real hx;		/* Space step */
} STR;

/****************/
RUN_HEAD(d_dy)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(real,hx);
  real one_o_2h=0.5/hx;
  int V0=v0*DV;
  int V1=v1*DV;
#define COMMANDS u[V1]=(u[V0+DY]-u[V0-DY])*one_o_2h;
  DO_FOR_ALL_POINTS;
}
RUN_TAIL(d_dy)

DESTROY_HEAD(d_dy)
DESTROY_TAIL(d_dy)

CREATE_HEAD(d_dy)
{
  DEVICE_REQUIRES_SYNC;
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT(v1>=v0);
  ACCEPTR(hx,RNONE,RSUCC(0),RNONE);
}
CREATE_TAIL(d_dy,1)
