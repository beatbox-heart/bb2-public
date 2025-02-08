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

/* 
 * Periodic boundary conditions along x axis.
 * Only for box geometry, sequential mode.
 */
#include <assert.h>
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
  int v0, v1;
} STR;

/****************/
RUN_HEAD(torx)
{
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  int y, z, v;
  real *u;

  for(y=s.y0;y<=s.y1;y++) for(z=s.z0;z<=s.z1;z++) for(v=v0;v<=v1;v++) {
    New[ind(s.x0-1,y,z,v)] = New[ind(s.x1,y,z,v)];
    New[ind(s.x1+1,y,z,v)] = New[ind(s.x0,y,z,v)];
  }
}
RUN_TAIL(torx)

#include "argsused.h"
DESTROY_HEAD(torx)
DESTROY_TAIL(torx)

#include "argsused.h"
CREATE_HEAD(torx)
{
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT(v1 >= v0);
  Space s=dev->s;
  ASSERT(s.x0>1)
  ASSERT(s.x1<xmax-1)
}
CREATE_TAIL(torx,1)

