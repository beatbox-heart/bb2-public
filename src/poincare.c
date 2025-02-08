/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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
 *	Finds a crossing (can be -ve, +ve or both, in time
 *	for a given variable at a given point.
 *
 *	Assigns the result to a given k_variable.
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "system.h"
#include "beatbox.h"
#include "error.h"
#include "device.h"
#include "state.h"
#include "qpp.h"

#define POSITIVE 1
#define BOTH 0
#define NEGATIVE -1

typedef struct {
#if MPI
  int root;		/* rank of process holding the point to be sampled */
#endif
  int first;		/* flags the first use, so we don't detect a crossing from the default value */
  real then;		/* value of result when last run. */
  real cross;		/* value at the crossing point. */
  REAL *result; 	/* k_var to which result will be assigned. */
  char *resultname;
  REAL *timestep; 	/* k_var to which the timestep at which the crossing occurs will be assigned. */
  char *timestepname; 
  int sign;		/* POSITIVE/NEGATIVE/BOTH: direction in which crossings should be detected */
  int lasttime;	     /* the value of t this device was last called */
} STR;

/******************/
RUN_HEAD(poincare)
#if MPI
  DEVICE_CONST(int,root);
  real buffer[2];
#endif
  DEVICE_CONST(real,then);
  DEVICE_CONST(real,cross);
  DEVICE_CONST(real *,result);
  DEVICE_CONST(real *,timestep);
  DEVICE_CONST(int,sign);
  DEVICE_CONST(int,lasttime);
  int crossingDetected = 0;
  real now = 0;
#if MPI
  if (mpi_rank == root) {
#endif
    now = New[ind(s.x0,s.y0,s.z0,s.v0)];	/* Current value */
#if MPI
  } /* if root */
#endif
  if (S->first) {
    S->first = 0;
  } else {
#if MPI
    if (mpi_rank == root) {
#endif
      /* Detect Crossing */
      if( (sign != POSITIVE && (then > cross && cross >= now)) || /* Negative crossing */
	  (sign != NEGATIVE && (then < cross && cross <= now)) ){ /* Positive crossing */
	crossingDetected = 1;
	/* printf("Crossing detected at t=%ld.\n",t); */
	if(timestep!=NULL) *timestep=(t*then+lasttime*now)/(then+now);
      }
      *result = crossingDetected;
      
#if MPI
    }/* if root */
    /* ------------- DISTRIBUTE RESULT TO ALL PROCESSES ---------------*/
    buffer[0] = *result;
    if(timestep) buffer[1] = *timestep;
    MPI_Bcast(&buffer, timestep?2:1, MPI_DOUBLE, root, ALL_ACTIVE_PROCS);
    *result = buffer[0];
    if(timestep) *timestep = buffer[1];
#endif
  }/* first=0 */
  S->then = now;
  S->lasttime = t;
RUN_TAIL(poincare)

/******************/
DESTROY_HEAD(poincare)
DESTROY_TAIL(poincare)

/******************/
CREATE_HEAD(poincare)
  DEVICE_ALWAYS_RUNS
  DEVICE_OPERATES_ON_A_SINGLE_POINT
  ACCEPTV(result);
  S->timestep = NULL;
  if (find_key("timestep=",w)) {
    ACCEPTV(timestep);
  }
  ACCEPTR(cross,RNONE,RNONE,RNONE);
  ACCEPTI(sign,BOTH,NEGATIVE,POSITIVE);
  S->first = 1;
#if MPI
  S->root = getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0);
#endif
CREATE_TAIL(poincare,1)

#undef POSITIVE
#undef BOTH
#undef NEGATIVE

