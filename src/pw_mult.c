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

/* Copy from one layer to another */

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
  int v0, v1;	/* argument layers */
  int vres;	/* result layer    */
} STR;

/****************/
RUN_HEAD(pw_mult)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(int,vres);
  int DV0=DV*v0;
  int DV1=DV*v1;
  int DVR=DV*vres;
  int x, y, z;
  real *u;
  for(z=s.z0;z<=s.z1;z++) for(y=s.y0;y<=s.y1;y++) for(x=s.x0;x<=s.x1;x++) {
	u=New+ind(x,y,z,0);
	u[DVR] = u[DV1] * u[DV0];
      }
}
RUN_TAIL(pw_mult)

DESTROY_HEAD(pw_mult)
DESTROY_TAIL(pw_mult)

CREATE_HEAD(pw_mult)
{
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ACCEPTI(vres,INONE,0,(int)vmax-1);
}	
CREATE_TAIL(pw_mult,1)

