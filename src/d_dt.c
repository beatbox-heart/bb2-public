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

/* Differentiation by time device */

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
  real ht;			/* time step for differentiation */
  int v0;			/* layer for the source field */
  int v1;			/* layer for the answer */
  int vd;			/* layer for previous values 	 */
} STR;

/****************/
RUN_HEAD(d_dt)
{
  DEVICE_CONST(real,ht);
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(int,vd);
  int x, y, z;
  real *u;
  for(z=s.z0;z<=s.z1;z++) for(y=s.y0;y<=s.y1;y++) for(x=s.x0;x<=s.x1;x++) {
	u=New+ind(x,y,z,v0);
	u[DV*(vd-v0)] = (u[0]-u[DV*(v1-v0)])/ht;
	u[DV*(v1-v0)] = u[0];
      }
}
RUN_TAIL(d_dt)

DESTROY_HEAD(d_dt)
DESTROY_TAIL(d_dt)

CREATE_HEAD(d_dt)
{
  DEVICE_REQUIRES_SYNC;
  ACCEPTR(ht,RNONE,RSUCC(0),RNONE);
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ACCEPTI(vd,INONE,0,(int)vmax-1);
  ASSERT( S->v1 != S->v0 );
  ASSERT( S->vd != S->v0 );
  ASSERT( S->vd != S->v1 );
}
CREATE_TAIL(d_dt,1)

