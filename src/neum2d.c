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
 * Experimental module, not for distribution. 
 * Provide Neumann boundary conditions for the 
 * experimental diffusion devices diff2dv, diff2vh. 
 *
 * Currently only implemented for box geometry
 * and tested in sequential runs. 
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
  int v0;
  int v1;
} STR;

/****************/
RUN_HEAD(neum2d){
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  int x, y, z, v;
  real *u;

  for(v=v0;v<=v1;v++) for(z=s.z0;z<=s.z1;z++) {
    for(x=s.x0;x<=s.x1;x++) {
      if (mpi_iy==0)        {y=s.y0; u=New+ind(x,y,z,v); u[-DY]=u[0];} /* bottom */
      if (mpi_iy==mpi_ny-1) {y=s.y1; u=New+ind(x,y,z,v); u[+DY]=u[0];} /* top */
    }
    /* limits extended by 1 to include corner cells */
    for(y=s.y0-1;y<=s.y1+1;y++) {
      if (mpi_ix==0)        {x=s.x0; u=New+ind(x,y,z,v); u[-DX]=u[0];} /* left */
      if (mpi_ix==mpi_nx-1) {x=s.x1; u=New+ind(x,y,z,v); u[+DX]=u[0];} /* right */
    }
  }

}
RUN_TAIL(neum2d)

/****************/
DESTROY_HEAD(neum2d)
DESTROY_TAIL(neum2d)

/****************/
CREATE_HEAD(neum2d)
{
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT(v1 >= v0);
  if (dim!=2) MESSAGE("/* Warning: 2D device is used in %d dimensional simulation */",dim);
  if (GEOMETRY_ON) MESSAGE("/* Warning: this device is intended for boxes only */");
}
CREATE_TAIL(neum2d,1)

