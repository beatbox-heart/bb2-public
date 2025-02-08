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

/* Zeldovich-Frank-Kamenetsky equation = Nagumo equation = Schloegl model            */
/* Notations as in Keener and Sneyd, Mathematical Physiology, page 231, except a->A  */

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
  real A;			/* magnitude */
  real alpha;			/* threshold */
  real Iu; 			/* connection term */
} STR;


RHS_HEAD(zfk,1) {
  DEVICE_CONST(real,A)
  DEVICE_CONST(real,alpha)
  DEVICE_CONST(real,Iu)
  real U = u[0];

  du[0] = A*U*(U-alpha)*(1.0-U) + Iu;

} RHS_TAIL(zfk)

RHS_CREATE_HEAD(zfk)
  ACCEPTP(A,1.0,0.,RNONE);
  ACCEPTP(alpha,RNONE,RNONE,RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
RHS_CREATE_TAIL(zfk,1)
