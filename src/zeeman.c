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

/* Extended Zeeman's model */

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
    p, q, eps, ooe; 
  real 
    Iu, Iv, Iw; 			/* External influence */
} STR;

RHS_HEAD(zeeman,3) {
  DEVICE_CONST(real,p)
  DEVICE_CONST(real,q)
  DEVICE_CONST(real,ooe)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)
  DEVICE_CONST(real,Iw)
  real U, U3, V, W;
  U = u[0];
  V = u[1];
  W = u[2];
  U3=U*U; U3*=U;
  du[0] = - ooe*(U3+U*V+W) + Iu;
  du[1] = 2*(p+q-U-V) + Iv;
  du[2] = - V + q + Iw;
} RHS_TAIL(zeeman)

RHS_CREATE_HEAD(zeeman)
  ACCEPTP(eps,0.01,RSUCC(0.),RNONE);
  ACCEPTP(p,1,RNONE,RNONE);
  ACCEPTP(q,-1,RNONE,RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
  ACCEPTP(Iv,0,RNONE,RNONE);
  ACCEPTP(Iw,0,RNONE,RNONE);
  S->ooe = 1./(S->eps);

  MALLOC(*u,3*sizeof(real));
  (*u)[0]=S->p; 
  (*u)[1]=S->q;
  (*u)[2]=-(S->p)*((S->p)*(S->p)+(S->q));

RHS_CREATE_TAIL(zeeman,3)
