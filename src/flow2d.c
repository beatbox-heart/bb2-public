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

/* 2D diffusion-advection term */

/* d_x(C_xu+D_{xx}d_xu+D_{xy}d_yu)+d_y(C_yu+D_{yx}d_xu+D_{yy}d_yu) */
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
  real Dxx, Dxy, Dyx, Dyy, Cx, Cy;
  real hx;
  int vi;
} STR;

RHS_HEAD(flow2d,1)
{
  DEVICE_CONST(real,Dxx)
  DEVICE_CONST(real,Dxy)
  DEVICE_CONST(real,Dyx)
  DEVICE_CONST(real,Dyy)
  DEVICE_CONST(real,Cx)
  DEVICE_CONST(real,Cy)
  DEVICE_CONST(real,hx)
  DEVICE_CONST(int,vi)
  real jxp, jxm, jyp, jym;
  real bet=1./hx;

  jxp=Cx*(u[0]+u[+DX]) + bet*(Dxx*(u[+DX]-u[0]) + Dxy*0.25*(u[+DX+DY]+u[DY]-u[+DX-DY]-u[-DY]));
  jxm=Cx*(u[0]+u[-DX]) + bet*(Dxx*(u[0]-u[-DX]) + Dxy*0.25*(u[-DX+DY]+u[DY]-u[-DX-DY]-u[-DY]));
  jyp=Cy*(u[0]+u[+DY]) + bet*(Dyy*(u[+DY]-u[0]) + Dyx*0.25*(u[+DX+DY]+u[DX]-u[-DX+DY]-u[-DX]));
  jym=Cy*(u[0]+u[-DY]) + bet*(Dyy*(u[0]-u[-DY]) + Dyx*0.25*(u[+DX-DY]+u[DX]-u[-DX-DY]-u[-DX]));
  u[vi] = bet*(jxp-jxm+jyp-jym);

}
RHS_TAIL(flow2d)

RHS_CREATE_HEAD(flow2d)
{
  ACCEPTP(Dxx,0.,0.,RNONE);
  ACCEPTP(Dxy,0.,RNONE,RNONE);
  ACCEPTP(Dyx,S->Dxy,RNONE,RNONE);
  ACCEPTP(Dyy,0,0.,RNONE);
  ACCEPTP(Cx,0,0,RNONE);
  ACCEPTP(Cy,0,0,RNONE);
  ACCEPTR(hx,RNONE,0.,RNONE); ASSERT(hx != 0);
  ACCEPTI(vi,INONE,0,INONE);
}
RHS_CREATE_TAIL(flow2d,1)


