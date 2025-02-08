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

/*  Utility functions for devices. */

#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include "dynamic.h"
#include "beatbox.h"
#include "state.h"
#define OWN
#include "device.h"
#undef OWN

#if MPI
/* 	Define new group and communicator for I/O, 	*/
/*	including only instances with runHere = true.	*/
int deviceCommunicatorWithFirstRank (int runHere, MPI_Comm *new_comm, int *first)
{
  MPI_Group old_group;		/* Group of existing communicator (ALL_ACTIVE_PROCS) */
  MPI_Group new_group;		/* Group for new communicator */
  int runHeres[num_active_procs];	/*  runHere flags from all active processes */
  int group_count = 0;		/* Number of processes in the new group */
  unsigned long count;		/* An occupation number */
  int laziest_rank;		/* Rank of least occupied process */
  unsigned long laziest = ULONG_MAX;	/* Its occupation number */
  int i, next, success;		/* Handy integers */

	
  MPIDO(MPI_Allgather(&runHere, 1, MPI_INT, &runHeres, 1, MPI_INT, ALL_ACTIVE_PROCS),"Couldn't gather runHere flags.");
	
  /*  Count the running instances. */
  for (i=0; i<num_active_procs; i++) if(runHeres[i]) group_count++;

  /* Empty group can occur if a device's space is entirely in the void. */
  /* No communicator then. Let the caller decide if this is fatal.      */
  if (group_count==0) {
    MESSAGE("Zero communicator group. The caller device's space may be entirely in the void.\n");
    return 0;
  }
	
  /*  Identify ranks of running instances to include in new group. */
  int group_ranks[group_count];

  next = 0;
  for (i=0; i<num_active_procs; i++) {
    if (runHeres[i]) {
      group_ranks[next] = i;
      count=occupy_rank[i];
      if (count<laziest) {
	laziest_rank=i;
	laziest=count;
      }
      next++;
    }
  }

  /* *first = group_ranks[0]; */
  *first = laziest_rank;
	
  MPIDO(MPI_Comm_group(ALL_ACTIVE_PROCS, &old_group),"Couldn't get old communicator's group.");
  MPIDO(MPI_Group_incl(old_group, group_count, group_ranks, &new_group),"Couldn't make new group.");
  MPIDO(MPI_Comm_create(ALL_ACTIVE_PROCS, new_group, new_comm),"Couldn't create new communicator.")
	
  return 1;
}

/*  Overloaded for the majority of occasions where 'first' isn't needed. */
int deviceCommunicator (int runHere, MPI_Comm *new_comm) {
  int first;
  return deviceCommunicatorWithFirstRank(runHere, new_comm, &first);
}

/* Here it is assumed that NORANK is negative, so neighbour>=0 means real process. */
/* For dim<3, neighbours in non-existing dimensions will have NORANK and automatically ignored, */
/* so no need to check for dim separately. */
int haloSwap (void) {
  int nb; /* rank of a neighbour */

  /* X AXIS */
  if (mpi_ix % 2 == 0) {
    if (0<=(nb=XP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,XP_Type,nb,0,New,1,XP_Halo_Type,XP_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers."); 
    }
  } else { 
    if (0<=(nb=XN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,XN_Type,nb,0,New,1,XN_Halo_Type,XN_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers."); 
    }
  }

  if (mpi_ix % 2 == 1) { 
    if (0<=(nb=XP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,XP_Type,nb,0,New,1,XP_Halo_Type,XP_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers."); 
    } 
  } else { 
    if (0<=(nb=XN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,XN_Type,nb,0,New,1,XN_Halo_Type,XN_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    } 
  }

  /* Y AXIS */ 
  if (mpi_iy % 2 == 0) { 
    if (0<=(nb=YP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,YP_Type,nb,0,New,1,YP_Halo_Type,YP_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    } 
  } else { 
    if (0<=(nb=YN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,YN_Type,nb,0,New,1,YN_Halo_Type,YN_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    }
  } 

  if (mpi_iy % 2 == 1) { 
    if (0<=(nb=YP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,YP_Type,nb,0,New,1,YP_Halo_Type,YP_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    } 
  } else { 
    if (0<=(nb=YN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,YN_Type,nb,0,New,1,YN_Halo_Type,YN_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    } 
  }

  /* Z AXIS */
  if (mpi_iz % 2 == 0) { 
    if (0<=(nb=ZP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,ZP_Type,nb,0,New,1,ZP_Halo_Type,ZP_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    } 
  } else { 
    if (0<=(nb=ZN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,ZN_Type,nb,0,New,1,ZN_Halo_Type,ZN_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    } 
  }

  if (mpi_iz % 2 == 1) { 
    if (0<=(nb=ZP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,ZP_Type,nb,0,New,1,ZP_Halo_Type,ZP_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    }
  } else { 
    if (0<=(nb=ZN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(New,1,ZN_Type,nb,0,New,1,ZN_Halo_Type,ZN_NEIGHBOUR,0,ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),"Failed to swap buffers.");
    }
  }

  return 1;
}
#else
int haloSwap (void) {return 1;}
#endif

