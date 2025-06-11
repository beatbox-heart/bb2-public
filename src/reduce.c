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
 *	Carry out reduction operation over the space and assign result to a given k_variable.
 *	Available operations, selected in the user script are sum, product ('prod'), min and max.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


#include "system.h"
#include "beatbox.h"
#include "error.h"
#include "state.h"
#include "device.h"
#include "qpp.h"

typedef enum {r_max,r_min,r_sum,r_prod} reduce_t;

#define REDUCE(name) void name(real new_value, real *running)
typedef REDUCE(Reduce);

void reduce_max(real new_value, real *running){
  if(new_value > *running){
    *running = new_value;
  }
}

void reduce_min(real new_value, real *running){
  if(new_value < *running){
    *running = new_value;
  }
}

void reduce_sum(real new_value, real *running){
  *running += new_value;
}

void reduce_prod(real new_value, real *running){
  *running *= new_value;
}

typedef struct {
  int v0, v1;			/* Range of layers					*/
#if MPI
  MPI_Comm reduce_comm;		/* Communicator for reduction operation.		*/
  int bcast_root;		/* Rank of device root in ALL_ACTIVE_PROCS		*/
  MPI_Op global_reduce;		/* Global reduction operation to use.			*/
#endif
  char operation[10];		/* Name of reduction operation to use.			*/
  reduce_t code;                /* code (0..3) of the operation				*/
  Reduce *local_reduce;		/* Function performing the local reduction		*/
  KREAL(result); 		/* Address of k_var to which result will be assigned	*/
  NFILE(debug);
} STR;

/******************/
RUN_HEAD(reduce)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(reduce_t, code);
  DEVICE_VAR(Reduce *, local_reduce);
  DEVICE_CONST(FILE *, debug);

#if MPI
  DEVICE_CONST(MPI_Comm, reduce_comm);
  DEVICE_CONST(int, bcast_root);
  DEVICE_CONST(MPI_Op, global_reduce);

  /* All active processes need to know the S->result, */
  /* but only those within the space contribute to it. */
  if (s.runHere) {
#endif
    int x, y, z, v;
    int firstPass = 1;
    double this_value;
    double running;
    switch (code) {
    case (r_max): running=-MAXREAL; break;
    case (r_min): running=MAXREAL; break;
    case (r_sum): running=0; break;
    case (r_prod): running=1.0; break;
    }
    for (x=s.x0;x<=s.x1;x++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (z=s.z0;z<=s.z1;z++) {
	  for (v=v0;v<=v1;v++) {
	    this_value = New[ind(x,y,z,v)];
	    (*local_reduce)(this_value, &running);
	  } /* for v */
	} /* for z */
      } /* for y */
    } /* for x */
#if MPI
    MPIDO(MPI_Reduce(&running, S->result, 1, MPI_DOUBLE, global_reduce, 0, reduce_comm),
	  "Couldn't carry out reduction operation.");
    if (debug) fprintf(debug, "t=%ld : %s=%lf\n",t,S->resultname,*(S->result));
  } /* if(s.runHere) */

  /* ------------- DISTRIBUTE RESULT TO ALL PROCESSES ---------------*/
  MPI_Bcast(S->result, 1, MPI_DOUBLE, bcast_root, ALL_ACTIVE_PROCS);
#else
  *(S->result) = running;
  if (debug) fprintf(debug, "t=%ld : %s=%lf\n",t,S->resultname,*(S->result));
#endif
} RUN_TAIL(reduce)

/******************/
DESTROY_HEAD(reduce)
DESTROY_TAIL(reduce)

/******************/
CREATE_HEAD(reduce)
{
  DEVICE_ALWAYS_RUNS;
  DEVICE_IS_RECTANGULAR;
    
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 >= v0 );
  ACCEPTS(operation,"");
  if (strcmp(S->operation, "min")==0) {
    S->code = r_min;
    S->local_reduce = reduce_min;
#if MPI
    S->global_reduce = MPI_MIN;
#endif
  } else if (strcmp(S->operation, "max")==0) {
    S->code = r_max;
    S->local_reduce = reduce_max;
#if MPI
    S->global_reduce = MPI_MAX;
#endif
  } else if (strcmp(S->operation, "sum")==0) {
    S->code = r_sum;
    S->local_reduce = reduce_sum;
#if MPI
    S->global_reduce = MPI_SUM;
#endif
  } else if (strcmp(S->operation, "prod")==0) {
    S->code = r_prod;
    S->local_reduce = reduce_prod;
#if MPI
    S->global_reduce = MPI_PROD;
#endif
  } else {
    EXPECTED_ERROR("The given operation, \"%s\", was not recognised.\n"
		   "Accepted values are \"min\", \"max\", \"sum\" and \"prod\".\n",
		   S->operation);
  }
  ACCEPTKR(result,1,NULL);
  ACCEPTF(debug,"wt","");
#if MPI
  int reduce_rank, i_am_root, bcast_root, i;
  int rootFlags[num_active_procs];
  if (!deviceCommunicator(dev->s.runHere, &(S->reduce_comm))) 
    EXPECTED_ERROR("Could not create communicator.\n");
	
  /*	MPI_Reduce operates on the device communicator, defined where s.runHere is true.
   *	As the MPI_Broadcast that distributes the result must operate on ALL_ACTIVE_PROCS.
   *	The same process shall work as root in reduce communicator, 
   *    to gather the reduction result from the participating processes, 
   *	and in broadcast communicator, to distribute the result among all active processes. 
   *
   *	This process, as MPI_Reduce root, has rank 0 in the device communicator.
   *	Its rank, as of the MPI_Broadcast root, in the ALL_ACTIVE_PROCS communicator
   *    is computed here and stored in S->bcast_root.
   */

  i_am_root = 0;
  if (dev->s.runHere) {
    MPI_Comm_rank(S->reduce_comm,&reduce_rank);
    if (reduce_rank == 0) i_am_root = 1;
  }
  MPIDO(MPI_Allgather(&i_am_root, 1, MPI_INT, &rootFlags, 1, MPI_INT, ALL_ACTIVE_PROCS),
	"Couldn't gather root flags.");
  for (i=0; i<num_active_procs; i++) {
    if(rootFlags[i]) bcast_root = i;
  }
  S->bcast_root = bcast_root;
  
#endif
} CREATE_TAIL(reduce,1)

