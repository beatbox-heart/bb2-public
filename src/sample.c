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
 * Assigns the value of a point in the mesh to a given k_variable.
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "system.h"
#include "beatbox.h"
#include "error.h"
#include "state.h"
#include "device.h"
#include "qpp.h"

typedef struct {
  int v0;		/* layer of the source point */
#if MPI
  int the_rank;	/* rank of process holding the point to be sampled */
#endif
  KREAL(result); 	/* k_var to which result will be assigned. */
  NFILE(debug);
} STR;

/******************/
RUN_HEAD(sample)
{
  DEVICE_CONST(int,v0);
#if MPI
  DEVICE_CONST(int,the_rank);
#endif
  DEVICE_CONST(real *,result);
  DEVICE_CONST(char *,resultname);
#if MPI
  if (mpi_rank==the_rank) {
#endif
    DEVICE_CONST(FILE *,debug)
      *result = New[ind(s.x0,s.y0,s.z0,v0)];
    if (debug) fprintf(debug, "t=%ld : %s=%lf\n",t,resultname,*result);
#if MPI
  } /* if the_rank */
  MPI_Bcast(S->result, 1, MPI_DOUBLE, the_rank, ALL_ACTIVE_PROCS);
#endif
}
RUN_TAIL(sample)

/******************/
DESTROY_HEAD(sample)
DESTROY_TAIL(sample)

/******************/
CREATE_HEAD(sample)
{
  DEVICE_ALWAYS_RUNS; /* unless stated otherwise, see below */
  DEVICE_OPERATES_ON_A_SINGLE_POINT;
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTKR(result,1,NULL);
#if MPI
  S->the_rank = getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0);
  if (S->the_rank<0) {
    MESSAGE("/* Warning: point (%d,%d,%d) is void; this sample device will be idle */",
	    dev->s.global_x0,dev->s.global_y0,dev->s.global_z0);
    dev->s.runHere=0; /* NB this is zeroed out in all processes */
    dev->alwaysRun=0; /* here stated otherwise ! */
  } 

  /* Only the_rank process writes to debug */
  if (mpi_rank == S->the_rank) {
#endif
    ACCEPTF(debug,"wt","");
#if MPI
  }
#endif
}
CREATE_TAIL(sample,1)
