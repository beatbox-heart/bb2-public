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

/* Domain decomposition definitions */

#if MPI
#ifndef _DECOMP_H_
#define _DECOMP_H_
#ifdef _OWN
#undef _OWN
#define OWN
#endif
#include "extern.h"

/*! Computes a preferred domain decomposition for a given number of processes.
 *	\param num_procs The number of available processes.
 *	\param nx Address at which the x dimension of the supergrid should 
 *                be stored.
 *	\param ny Address at which the y dimension of the supergrid should 
 *                be stored.
 *	\param nz Address at which the z dimension of the supergrid should 
 *                be stored.
 *	\return The number of active processes in the preferred decomposition.
 */
int decompose (int num_procs, long *nx, long *ny, long *nz);

/* Allocates subdomains to processes, omitting empty ones if any */
int decomp_allocateSubdomains(void);

/*! Computes the superindices of a subdomain, given its rank. 
 *	\param rank Rank of the process for which superindices are to be 
 *                  calculated.
 *	\param ix Address at which the x supergrid coordinate should be stored.
 *	\param iy Address at which the y supergrid coordinate should be stored.
 *	\param iz Address at which the z supergrid coordinate should be stored.
 */
void decomp_getSuperindices(int rank, long *ix, long *iy, long *iz);

/*! Computes the local subdomain's dimensions. 
 *	Results for all subdomains are stored in subdomain_limits.
 *	Results for local subdomain are stored in local_xmin, local_xmax, etc.
 */
int decomp_getSubdomainDimensions(void);

/*! Computes the local subdomain's size (including boundary/halo points). 
 *	Results are stored in xlen, ylen and zlen.
*/
void decomp_getSubdomainSize(void);

/*! Defines datatypes for halo swapping.
*/
void decomp_defineHaloTypes(void);

/*! Defines a global communicator over all active processes. 
 *	\param old_rank The process's rank in MPI_COMM_WORLD.
 *	\param new_comm Address at which the new communicator will be created.
 *	\param new_rank Address at which the process's rank in new_comm will 
 *             be stored. Will be equal to old_rank if the process is idle.
 *	\return 1 if the process is idle in the decomposition.
*/
int decomp_globalCommunicator(int old_rank, MPI_Comm *new_comm, int *new_rank);

#endif /* end of include guard: _DECOMP_H_ */

#endif /* MPI */
