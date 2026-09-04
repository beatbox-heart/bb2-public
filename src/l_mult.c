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

/* Layer arithmetics: y=a*x, same as l_scal */

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
  int v0;			/* source layer */
  int v1;			/* target layer */
  PAR(real,multiplier);		/* by what to multiply, mutiplicand, factor, coefficient, scalar */
} STR;

/****************/
RUN_HEAD(l_mult)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_PAR(real,multiplier);
#define COMMANDS u[v1]=u[v0]*multiplier;
  DO_FOR_ALL_POINTS;
#undef COMMANDS
}
RUN_TAIL(l_mult)

DESTROY_HEAD(l_mult)
DESTROY_TAIL(l_mult)

CREATE_HEAD(l_mult)
{
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ACCEPTRK(multiplier,RNONE,RNONE,RNONE);
  /* TODO: synonyms? */
}
CREATE_TAIL(l_mult,1)

