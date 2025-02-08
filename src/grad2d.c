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

/* Absolute value of 2D gradient, central difference */

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
RUN_HEAD(grad2d)
{
  DEVICE_CONST(real,hx);
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  int x, y, z;
  real *u;
  size_t V1=DV*(v1-v0);
  for(x=s.x0;x<=s.x1;x++) {
    for(y=s.y0;y<=s.y1;y++){
      for(z=s.z0;z<=s.z1;z++){
	u=New+ind(x,y,z,v0);
	u[V1] = hypot(u[DX]-u[-DX],u[DY]-u[-DY])/hx;
      } /* for z */
    } /* for y */
  } /* for x */
}
RUN_TAIL(grad2d)

DESTROY_HEAD(grad2d)
DESTROY_TAIL(grad2d)

CREATE_HEAD(grad2d)
{
  DEVICE_REQUIRES_SYNC;
  DEVICE_IS_RECTANGULAR;
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT(v1>=v0);
  ACCEPTR(hx,RNONE,RSUCC(0),RNONE);
}
CREATE_TAIL(grad2d,1)

