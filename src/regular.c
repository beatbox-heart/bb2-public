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

/* global result=eq(mod(t-t0,dt),0)*ge(t,t0)*le(t,t1) only quicker */

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"

#include "k_.h"

typedef struct {
  INT t0;		/* earliest loop count */
  INT t1;		/* latest loop count */
  INT dt;		/* interval between */
  KREAL(result); 	/* k_var to which result will be assigned. */
  NFILE(debug);
  int debugWriter;
} STR;

RUN_HEAD(regular)
{
  DEVICE_CONST(INT,t0);
  DEVICE_CONST(INT,t1);
  DEVICE_CONST(INT,dt);
  /* DEVICE_CONST(REAL *,result); */
  DEVICE_KREAL (result);
  DEVICE_CONST(int, debugWriter)
  DEVICE_CONST(FILE *,debug)
  *result
    = (t<t0)? 0.0
    : (t>t1)? 0.0
    : (0==(t-t0)%dt) ? 1.0
    : 0.0;
  /* printf("t=%ld [%ld:%ld:%ld] -> %lg\n",(long)t,(long)t0,(long)dt,(long)t1,(double)(*result)); */
  if (debug && debugWriter) {
    fprintf(debug, " t=%ld t0=%ld t1=%ld dt=%ld %s=%lg\n", (long)t, (long)t0, (long)t1, (long)dt, resultname, (double)(*result));
    fflush(debug);
  }
}
RUN_TAIL(regular)

DESTROY_HEAD(regular)
{
} 
DESTROY_TAIL(regular)

CREATE_HEAD(regular)
{
  DEVICE_IS_SPACELESS;
  ACCEPTL(t0,0,0,LNONE);
  ACCEPTL(dt,1,1,LNONE);
  ACCEPTL(t1,LONG_MAX,t0,LNONE);
  ACCEPTKR(result,1,NULL);
  ACCEPTF(debug,"wt","");
#if MPI
  /* Identifying the root process allows only one process to write debug output. */
  if (dev->s.nowhere) {
    S->debugWriter = (mpi_rank == 0);
  } else {
    S->debugWriter = (mpi_rank == getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0));
  }
#else
  S->debugWriter = 1;
#endif
}
CREATE_TAIL(regular,1)

