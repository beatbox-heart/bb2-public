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

/* Domain decomposition */

#include <mpi.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "error.h"
#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#define OWN
#include "decomp.h"
#undef OWN

#define DB(...) DEBUG(__VA_ARGS__)

/* Imported from beatbox.c */
extern int Decomp_version;
extern char *Decomp_string;
extern int Verbose;

/* Imported from geometry.c */
extern unsigned char *gpoints;

/****************************************/
/* Declarations of locally used objects */

/* Tuning parameters for the aim function */
const static int DECOMP_VOLUME_ADJUST = 1;
const static int DECOMP_SURFACE_AREA_ADJUST = 100;


/* This fully defines a decomposition */
typedef struct{
  int nx, ny, nz;
} Decomposition;


/*! Computes the dimensions of the worst (i.e. largest) subdomain of a given 
 *  decomposition.
 *	\param d    Decomposition in which to find the subdomain.
 *	\param xlen Address at which the x dimension of the subdomain 
 *                  should be stored.
 *	\param ylen Address at which the y dimension of the subdomain 
 *                  should be stored.
 *	\param zlen Address at which the z dimension of the subdomain 
 *                  should be stored.
 */
void decomp_getWorstSubdomainDimensions(Decomposition d, size_t *xlen, size_t *ylen, size_t *zlen);

/*! Computes the largest volume of the any subdomain in a given decomposition.
 *	\param d Decomposition in which to find the subdomain.
 *	\return Volume of the subdomain (pt^3).
 */
long decomp_getWorstSubdomainVolume(Decomposition d);

/*! Computes the largest surface area of the any subdomain in a given 
 *  decomposition.
 *	\param d Decomposition in which to find the subdomain.
 *	\return Surface area of the subdomain (pt^2).
 */
long decomp_getWorstSubdomainSurfaceArea(Decomposition d);

/*! Computes metric for comparison of decompositions.
 *	\param d Decomposition to be evaluated.
 *	\return Metric value. Smaller is better.
 */
double decomp_getMetric(Decomposition d);

/*! Compares two decompositions using decomp_getMetric().
 *	\param a First decomposition.
 *	\param b Second decomposition.
 *	\return The better of a and b, according to decomp_getMetric().
 *	\sa decomp_getMetric
 */
Decomposition decomp_getBest(Decomposition a, Decomposition b);


/************************/
/* Function definitions */

void decomp_getWorstSubdomainDimensions (Decomposition d, size_t *xlen, size_t *ylen, size_t *zlen) {
/*
 *xlen = (xmax / d.nx) + (xmax % d.nx == 0)?0:1;
 *ylen = (ymax / d.ny) + (ymax % d.ny == 0)?0:1;
 *zlen = (zmax / d.nz) + (zmax % d.nz == 0)?0:1;
*/

  *xlen = (xmax / d.nx) + (xmax % d.nx);
  *ylen = (ymax / d.ny) + (ymax % d.ny);
  *zlen = (zmax / d.nz) + (zmax % d.nz);
}

long decomp_getWorstSubdomainVolume (Decomposition d) {
  size_t xlen, ylen, zlen;
  decomp_getWorstSubdomainDimensions(d, &xlen, &ylen, &zlen);
  return xlen * ylen * zlen;
}

long decomp_getWorstSubdomainSurfaceArea (Decomposition d) {
  size_t xlen, ylen, zlen;
  decomp_getWorstSubdomainDimensions(d, &xlen, &ylen, &zlen);
  return 2*((xlen*ylen)+(ylen*zlen)+(xlen*zlen));
}

double decomp_getMetric (Decomposition d) {
  return 
    (float) (decomp_getWorstSubdomainVolume(d) * DECOMP_VOLUME_ADJUST) 
    +
    (decomp_getWorstSubdomainSurfaceArea(d) * DECOMP_SURFACE_AREA_ADJUST);
}

Decomposition decomp_getBest(Decomposition a, Decomposition b){	
  return (decomp_getMetric(a) < decomp_getMetric(b)) ? a : b;
}

/********************************************************/ 
/*! Computes how many processes will be allocated along each of the
 *  cartesian directions.
 *
 * \param num_procs Number of processes involved in the calculation.
 * \param nx Number of processes assigned along the x-direction.
 * \param ny Number of processes assigned along the y-direction.
 * \param nz Number of processes assigned along the z-direction.
 * \return Number of processes participating in the domanin decomposition.
 *
 */
int decompose (int num_procs, long *nx, long *ny, long *nz)
{
  int nx_max,ny_max,nz_max;
  int nx_tmp_max,ny_tmp_max,nz_tmp_max;
  Decomposition best, trial;
	
  /* Factorization variables. */
  int no,x;
  int factors[10000],max_factors = 0,factor_counter = 0;

  /*	Set maxima for nx, ny and nz
   *	At most one processor along an inactive axis.
   *	No more than one processor per point.
   */
  nx_max = (ONE) ? min(xmax, num_procs) : 1;
  ny_max = (TWO) ? min(ymax, num_procs) : 1;
  nz_max = (TRI) ? min(zmax, num_procs) : 1;

  /*
   * Different decomposition strategies (default set in beatbox.c):
   *
   * 0 - decomposition stated at the command-line using the -decomp flag.
   * 1 - RMF based decomposition.
   * 2 - SRK based decomposition.
   */
  switch (Decomp_version) {
  /****************************************************/
  case 0: /* explicitly stated in the command-line option using -decomp */
    if (Verbose)
      MESSAGE("/* User-defined domain decomposition */\n");
    if (3!=sscanf(Decomp_string,"%dx%dx%d",&(best.nx),&(best.ny),&(best.nz))) 
      EXPECTED_ERROR("could not interpret %s as NXxNYxNZ",Decomp_string);
    /* This check is only applicable for boxes as in geometries empty subdoms 
     * won't need procs 
     */
    if (GEOMETRY_ON==0 && (best.nx*best.ny*best.nz>num_procs))
      EXPECTED_ERROR("given NX*NY*NZ=%d*%d*%d=%d exceeds requested number of processes %d",
	    best.nx,best.ny,best.nz,best.nx*best.ny*best.nz,num_procs);
    if (best.nx>nx_max)
      EXPECTED_ERROR("given NX=%d is greater than allowed maximum of %d",best.nx,nx_max);
    if (best.ny>ny_max)
      EXPECTED_ERROR("given NY=%d is greater than allowed maximum of %d",best.ny,ny_max);
    if (best.nz>nz_max)
      EXPECTED_ERROR("given NZ=%d is greater than allowed maximum of %d",best.nz,nz_max);

    break;
  /****************************************************/
  case 1: /* RMF version */
    if (Verbose)
      MESSAGE("/* Automatic domain decomposition variant 1 */\n");
    best.nx = 1;
    best.ny = 1;
    best.nz = 1;
    for (trial.nx = 1; trial.nx <= nx_max; trial.nx++) {
      ny_tmp_max = min (ny_max, num_procs/trial.nx);
      for (trial.ny = 1; trial.ny <= ny_tmp_max; trial.ny++) {
	nz_tmp_max = min (nz_max, num_procs/(trial.nx*trial.ny));
	for(trial.nz = 1; trial.nz <= nz_tmp_max; trial.nz++) {
	  best = decomp_getBest(trial, best);
	}
      }
    }
    break;
  /****************************************************/
  case 2: /* SRK version */
    if (Verbose)
      MESSAGE("/* Automatic domain decomposition variant 2 */\n");
    /* TEST: Hard code the 4 threads to see if it helps any. 
     * Answer: it does help! 
     */
    /*
      best.nx = 4;
      best.ny = 1;
      best.nz = 1;
    */
    
    /* this is 1D case specific. */
    /*
      The conditions to check for this are:
      1. get the num_procs: already done.
      2. No geometry should be on.
      3. It should be 1D along the x axis
      4. Always in the MPI #if: this is already true by the time the code 
         gets here.
      
      Some more conditions have to be considered. The way the num_procs are
      being partitioned, it assumes that num_procs is at least size
      dim. This need not always be true.
    */
    
    if (dim==1) { /* this condition alone makes sure that we have satisfied conditions 1-4. */
      best.nx = num_procs;
      best.ny = 1;
      best.nz = 1;
    }  /* end of dim == 1 */
    
    if (dim==2 /* &&!GEOMETRY_ON */) { /* at this stage, there could be a slice. so dont include geometry yet. */
      
      if (num_procs==1||num_procs==2) {
	
	best.nx = num_procs;
	best.ny = 1;
	best.nz = 1;
	
      } else {
	
	no = num_procs;
	
	for(x=1; x<=no; x++) {
	  if (no%x==0) {
	    factors[factor_counter] = x;
	    factor_counter++;
	  }
	}
	
	/* if it is not a perfect square, then the number of factors is even by now. */
	if(factor_counter%2==0){
	  best.nx = factors[(factor_counter-1)/2];
	  best.ny = factors[(factor_counter+1)/2];
	  best.nz = 1;
	}
	
	/* a perfect square. */
	if (factor_counter%2!=0) {
	  best.nx = factors[(factor_counter-1)/2];
	  best.ny = factors[(factor_counter-1)/2];
	  best.nz = 1;
	}
	
      }	 /* end of numproc>2 case */
      
    } /* end of dim==2 */
    
    /* I have not done the 3D partition correctly yet. The 3D geometry also works. */
    /* This code is geared towards 1,2,3, and perfect cube num_procs. */
    /* The code for other values of num_procs is required to cover all bases. */
    if (dim==3 /*&&!GEOMETRY_ON */) { 
      
      if (num_procs==3||num_procs==2||num_procs==1) {
	best.nx = num_procs;
	best.ny = 1;
	best.nz = 1;
      } else if( /* this should be ok for now. */
		num_procs==2*2*2 || num_procs==3*3*3 || num_procs==4*4*4 ||
		num_procs==5*5*5 || num_procs==6*6*6 || num_procs==7*7*7 ||
		num_procs==8*8*8 || num_procs==9*9*9 || num_procs==10*10*10 ||
		num_procs==11*11*11 ) { 
	/* this is only true for full cubes as of now.  */
	/* Improve the cube root finding, or hard code it in for the first  */
        /* several cubes. */
	best.nx = (int)ceil(pow((double)num_procs,1.0/3.0));
	best.ny = (int)ceil(pow((double)num_procs,1.0/3.0));
	best.nz = (int)ceil(pow((double)num_procs,1.0/3.0));
      } else {
	no = num_procs;
	
	for (x=1; x<=no; x++) {
	  if (no%x==0) {
	    factors[factor_counter] = x;
	    factor_counter++;
	  }
	}
	
	/* if it is not a perfect square, then the number of factors is 
         * even by now. 
         */
	if (factor_counter%2==0) {
	  best.nx = factors[(factor_counter-1)/2];
	  best.ny = factors[(factor_counter+1)/2];
	  best.nz = 1;
	}
	
	/* a perfect square. */
	if (factor_counter%2!=0) {
	  best.nx = factors[(factor_counter-1)/2];
	  best.ny = factors[(factor_counter-1)/2];
	  best.nz = 1;
	}
	
      } /* end of numproc>3 case */
      
    } /* end of dim==3 */
    break;
  default:
    EXPECTED_ERROR("Unknown decomposition method %d\n", Decomp_version);
  } /* switch Decomp_version */

  /* Bounds of the subdomains, however determined. */
  size_t ix, iy, iz;
  /* Approximate length of subdomain along each axis.	*/
  size_t x_slice_size = (xmax-(ONE*2)) / best.nx;
  size_t y_slice_size = (ymax-(TWO*2)) / best.ny;
  size_t z_slice_size = (zmax-(TRI*2)) / best.nz;
  
  /* Remainders from i[xyz] * [xyz]_slice_size */
  /* The remainder will be shared amongst the first processes on the axis */
  int x_remainder  = (xmax-(ONE*2)) % best.nx;
  int y_remainder  = (ymax-(TWO*2)) % best.ny;
  int z_remainder  = (zmax-(TRI*2)) % best.nz;
  
  /* Calculate upper and lower bounds along x axis */
  CALLOC(all_xmin,best.nx,sizeof(size_t));
  CALLOC(all_xmax,best.nx,sizeof(size_t));
  for (ix=0;ix<best.nx;ix++) {
    all_xmin[ix] = ONE + ( ix    * x_slice_size) + (ix < x_remainder?  ix    :x_remainder);
    all_xmax[ix] = ONE + ((ix+1) * x_slice_size) + (ix < x_remainder? (ix+1) :x_remainder);
  }
  
  /* Calculate upper and lower bounds along y axis */
  CALLOC(all_ymin,best.ny,sizeof(size_t));
  CALLOC(all_ymax,best.ny,sizeof(size_t));
  for (iy=0;iy<best.ny;iy++) {
    all_ymin[iy] = TWO + ( iy    * y_slice_size) + (iy < y_remainder?  iy    :y_remainder);
    all_ymax[iy] = TWO + ((iy+1) * y_slice_size) + (iy < y_remainder? (iy+1) :y_remainder);
  }
  
  /* Calculate upper and lower bounds along z axis */
  CALLOC(all_zmin,best.nz,sizeof(size_t));
  CALLOC(all_zmax,best.nz,sizeof(size_t));
  for (iz=0;iz<best.nz;iz++) {
    all_zmin[iz] = TRI + ( iz    * z_slice_size) + (iz < z_remainder?  iz    :z_remainder);
    all_zmax[iz] = TRI + ((iz+1) * z_slice_size) + (iz < z_remainder? (iz+1) :z_remainder);
  }

  /* Return the results */
  *nx = best.nx;
  *ny = best.ny;
  *nz = best.nz;
  return best.nx * best.ny * best.nz;
} /* decompose */


/* Allocates subdomains to processes, omitting any empty ones. */
/* Takes num_subdoms and defines num_active_procs.             */
int decomp_allocateSubdomains (void)
{
  int ix,iy,iz;	/* looping superindices */
  size_t x,y,z;	/* indices within subdomain */
  unsigned long count;	/* num of points in the current subdomain */
  unsigned long laziest;	/* .., minimal of */
  int n;	/* number of nonempty subdomains -> new num_active_procs */
  int rank;	/* nonempty subdomain counter */
  int corank;	/* empty subdomain counter */

  ASSERT(num_subdoms == mpi_nx * mpi_ny * mpi_nz);
  CALLOC(occupy_grid,num_subdoms,sizeof(unsigned long));
  CALLOC(occupy_rank,num_subdoms,sizeof(unsigned long));
    
  if (GEOMETRY_ON) {

    /* Geometry version                                                      */
    /* Allocate only non-empty subdomains, for processes with smaller ranks. */
    /* Check that there are enough processes available,                      */
    /* but tolerates if there are too many.                                  */

    /* This is where the geometry skeleton is required */
    ASSERT(gpoints!=NULL);
    
    /* Mark and count non-empty subdomains */
    n=0;
    for (ix=0;ix<mpi_nx;ix++) {
      for (iy=0;iy<mpi_ny;iy++) {
	for (iz=0;iz<mpi_nz;iz++) {
	  count=0;
	  for (x=all_xmin[ix];x<all_xmax[ix];x++) {
	    for (y=all_ymin[iy];y<all_ymax[iy];y++) {
	      for (z=all_zmin[iz];z<all_zmax[iz];z++) {
		long int index=gind(x,y,z);
		long int volume=xmax*ymax*zmax;
		if (0>index || index>=volume)
		  DB("domain(%d,%d,%d)=[%ld:%ld]:[%ld:%ld]:[%ld:%ld] "
		     "xyz=%ld,%ld,%ld -> index=%ld vs volume=%ld=%ld*%ld*%ld\n",
		     (int)ix,(int)iy,(int)iz,
		     (long)all_xmin[ix],(long)all_xmax[ix],
		     (long)all_ymin[iy],(long)all_ymax[iy],
		     (long)all_zmin[iz],(long)all_zmax[iz],
		     (long)x,(long)y,(long)z,(long)index,
		     (long)volume,(long)xmax,(long)ymax,(long)zmax);
		if (gpoints[index]) count++;
		/* if (gpoints[gind(x,y,z)]) count++; */
	      } /* for z */
	    } /* for y */
	  } /* for x */
	  occupy(ix,iy,iz)=count;
	  if (count>0) n++;
	} /* for iz */
      } /* for iy */
    } /* for ix */
    if (!n) EXPECTED_ERROR("No non-empty subdomains found in this geometry.");

    /* Allocate memory for the tables */
    num_active_procs = n;
    num_empty_subdoms = num_subdoms - n;
    num_inactive_procs = mpi_size - n;
    CALLOC(ProcessRank,num_subdoms,sizeof(int));
    CALLOC(rank_ix,num_subdoms,sizeof(int));
    CALLOC(rank_iy,num_subdoms,sizeof(int));
    CALLOC(rank_iz,num_subdoms,sizeof(int));

    /* Allocate non-empty subdomains to processes, */
    /* and identify empty ones just in case */
    rank=0;
    corank=num_active_procs;
    for (ix=0;ix<mpi_nx;ix++) {
      for (iy=0;iy<mpi_ny;iy++) {
	for (iz=0;iz<mpi_nz;iz++) {
	  if ((count=occupy(ix,iy,iz))>0) {
	    ASSERT(rank < num_active_procs);
	    ProcessRank[sind(ix,iy,iz)]=rank;
	    rank_ix[rank]=ix;
	    rank_iy[rank]=iy;
	    rank_iz[rank]=iz;
	    rank++;
	  } else { /* not occupy */
	    ASSERT(corank < num_subdoms);
	    ProcessRank[sind(ix,iy,iz)]=NORANK;
	    rank_ix[corank]=ix;
	    rank_iy[corank]=iy;
	    rank_iz[corank]=iz;
	    corank++;
	  } /* not occupy */
	} /* for iz */
      } /* for iy */
    } /* for ix */

  } else { /* not GEOMETRY_ON */

    /* Box version, no empty subdomains.                                  */
    /* All subdomains are allocated to processes identified by decompose, */
    /* and ordered straightforwardly as defined by sind() macro.          */
    num_active_procs = num_subdoms;
    num_empty_subdoms = 0;
    num_inactive_procs = mpi_size - num_subdoms;

    /* Allocate memory for the tables */
    CALLOC(ProcessRank,num_subdoms,sizeof(int));
    CALLOC(rank_ix,num_subdoms,sizeof(int));
    CALLOC(rank_iy,num_subdoms,sizeof(int));
    CALLOC(rank_iz,num_subdoms,sizeof(int));

    /* Allocate subdomains to the processes */
    for(ix=0;ix<mpi_nx;ix++) {
      for(iy=0;iy<mpi_ny;iy++) {
	for(iz=0;iz<mpi_nz;iz++) {
	  rank=sind(ix,iy,iz);
	  rank_ix[rank]=ix;
	  rank_iy[rank]=iy;
	  rank_iz[rank]=iz;
	  ProcessRank[rank]=rank;
	  occupy(ix,iy,iz)
	    =(all_xmax[ix]-all_xmin[ix]+1)
	    *(all_ymax[iy]-all_ymin[iy]+1)
	    *(all_zmax[iz]-all_zmin[iz]+1);
	} /* for iz */
      } /* for iy */
    } /* for ix */
  } /* not GEOMETRY_ON */

  /* Mark occupation numbers */
  laziest=ULONG_MAX;
  laziest_rank=NORANK;
  for (rank=0;rank<num_active_procs;rank++) {
    count=occupy(rank_ix[rank],rank_iy[rank],rank_iz[rank]);
    occupy_rank[rank]=count;
    if (count<laziest) {
      laziest_rank=rank;
      laziest=count;
    }	
  }

  return num_active_procs;
} /* decomp_allocateSubdomains */

/*  Define communicator over active processes */
int decomp_globalCommunicator (int old_rank, MPI_Comm *new_comm, int *new_rank)
{
  int i,is_idle;
  MPI_Group old_group, new_group;
  int group_ranks[num_active_procs];

  is_idle = (old_rank >= num_active_procs) ? 1 : 0;

  /* Ranks 0..(num_active_procs-1) are included in new communicator. */
  for (i=0; i<num_active_procs; i++)  group_ranks[i] = i;

  MPIDO( MPI_Comm_group(MPI_COMM_WORLD, &old_group), "Couldn't get old communicator's group.");
  MPIDO( MPI_Group_incl(old_group, num_active_procs, group_ranks, &new_group), "Couldn't make new group.");
  MPIDO( MPI_Comm_create(MPI_COMM_WORLD, new_group, new_comm), "Couldn't create new communicator.");

  if (!is_idle) {
    /*	Redefine rank based on new communicator
     *	This shouldn't have changed, but we can't guarantee it.	*/
    MPIDO(MPI_Comm_rank(*new_comm, new_rank),"Couldn't get my new rank.");
  } else {
    /* If the process is idle, we keep the old rank so that it can report. */
    *new_rank = old_rank;
  }
	
  return is_idle;
} /* decomp_globalCommunicator */


/* Map rank -> superindices (report pre-computed values from the tables) */
void decomp_getSuperindices (int rank, long *ix, long *iy, long *iz)
{
  *ix=rank_ix[rank];
  *iy=rank_iy[rank];
  *iz=rank_iz[rank];
} /* decomp_getSuperindices */

/* Compute bounds of all subdomains:                                     */
/* Each process contains active points from local_*min to local_*max -1. */
/* Repackage the bounds of partitions computed by decompose.             */
int decomp_getSubdomainDimensions (void)
{
  int rank;
  long ix, iy, iz;
  CALLOC(subdomain_limits,num_active_procs,sizeof(Subdomain));

  for (rank=0;rank<num_active_procs;rank++) {
    decomp_getSuperindices(rank,&ix,&iy,&iz);
    subdomain_limits[rank].local_xmin = all_xmin[ix];
    subdomain_limits[rank].local_xmax = all_xmax[ix];
    subdomain_limits[rank].local_ymin = all_ymin[iy];
    subdomain_limits[rank].local_ymax = all_ymax[iy];
    subdomain_limits[rank].local_zmin = all_zmin[iz];
    subdomain_limits[rank].local_zmax = all_zmax[iz];
  } /* for(rank) */
	
  /*	Assign own subdomain boundaries */
  local_xmin = subdomain_limits[mpi_rank].local_xmin;
  local_xmax = subdomain_limits[mpi_rank].local_xmax;
  local_ymin = subdomain_limits[mpi_rank].local_ymin;
  local_ymax = subdomain_limits[mpi_rank].local_ymax;
  local_zmin = subdomain_limits[mpi_rank].local_zmin;
  local_zmax = subdomain_limits[mpi_rank].local_zmax;

  return SUCCESS;
} /* decomp_getSubdomainDimensions */

void decomp_getSubdomainSize (void)
{
  xlen = (local_xmax-local_xmin)+(ONE*2);
  ylen = (local_ymax-local_ymin)+(TWO*2);
  zlen = (local_zmax-local_zmin)+(TRI*2);
  DEBUG("\n#%d:"
	 "local_(x,y,z):(max,min)=(%ld:%ld,%ld:%ld,%ld:%ld) (x,y,z)len=(%ld,%ld,%ld)\n",
	 mpi_rank,
	 (long)local_xmin,(long)local_xmax,
	 (long)local_ymin,(long)local_ymax,
	 (long)local_zmin,(long)local_zmax,
	 (long)xlen,(long)ylen,(long)zlen);
} /* decomp_getSubdomainSize */

void decomp_defineHaloTypes (void)
{
  /*	Types for Halo Swapping
   *
   *	Labelled as XN, XP, YN, YP, ZN, ZP
   *	where, e.g. XN is the YZV hyperplane at local_xmin.
   *	ZP is the XYV hyperplane at local_zmax.
   *
   *	ZP_Halo_Type is the halo at local_xmax.
   *	YN_Halo_Type is the halo at local_ymin-1.
   *
   *	Hyperplanes exchanged in the y axis are widened in
   *	the x axis, and those exchanged in the z axis are widened in
   *	in the x and y axes. This includes points from diagonally 
   *	neighbouring subdomains, which prevents the need for exchanging 
   *	corner points explicitly between these processes.
   *	-----------------------------------	*/
  #define NDIMS 4
  int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
  sizes[0] = xlen;
  sizes[1] = ylen;
  sizes[2] = zlen;
  sizes[3] = vmax;

  if (mpi_nx > 1) {
    /*	XN	*/
    subsizes[0] = 1;
    subsizes[1] = local_ymax-local_ymin;
    subsizes[2] = local_zmax-local_zmin;
    subsizes[3] = vmax;
	
    starts[0] = ONE;
    starts[1] = TWO;
    starts[2] = TRI;
    starts[3] = 0;
	
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&XN_Type), "Couldn't define XN_Type.");
    MPIDO(MPI_Type_commit(&XN_Type),"Couldn't commit XN_Type.");

    /*	XN Halo	*/
    starts[0] -= 1;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&XN_Halo_Type), "Couldn't define XN_Halo_Type at mpi_iq=(%ld,%ld,%ld) our of (%ld,%ld,%ld)",mpi_ix,mpi_iy,mpi_iz,mpi_nx,mpi_ny,mpi_nz);
    MPIDO(MPI_Type_commit(&XN_Halo_Type),"Couldn't commit XN_Halo_Type.");

    /*	XP	*/
    starts[0] = (local_xmax-1 + ONE) - local_xmin;
    starts[1] = TWO;
    starts[2] = TRI;
    starts[3] = 0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&XP_Type),"Couldn't define XP_Typea.");
    MPIDO(MPI_Type_commit(&XP_Type),"Couldn't commit XP_Type.");

    /*	XP Halo	*/
    starts[0] += 1;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&XP_Halo_Type),"Couldn't define XP_Halo_Type.");
    MPIDO(MPI_Type_commit(&XP_Halo_Type),"Couldn't commit XP_Halo_Type.");
  } /*	if(mpi_nx > 1)	*/
  
  if (mpi_ny > 1) {
    /*	YN	*/
    subsizes[0] = xlen; /*	Widened to include magic corners.	*/
    subsizes[1] = 1;
    subsizes[2] = local_zmax-local_zmin;
    subsizes[3] = vmax;
    
    starts[0] = 0;		/*	0 to include magic corner.	*/
    starts[1] = TWO;
    starts[2] = TRI;
    starts[3] = 0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&YN_Type),"Couldn't define YN_Type.");
    MPIDO(MPI_Type_commit(&YN_Type),"Couldn't commit YN_Type.");
	
    /*	YN Halo		*/
    starts[1] -= 1;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&YN_Halo_Type),"Couldn't define YN_Halo_Type.");
    MPIDO(MPI_Type_commit(&YN_Halo_Type),"Couldn't commit YN_Halo_Type.");
      
    /*	YP	*/
    starts[0] = 0;		/*	0 to include magic corner.	*/
    starts[1] = (local_ymax-1 + TWO) - local_ymin;
    starts[2] = TRI;
    starts[3] = 0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&YP_Type),"Couldn't define YP_Type.");
    MPIDO(MPI_Type_commit(&YP_Type),"Couldn't commit YP_Type.");

    /*	YP Halo	*/
    starts[1] += 1;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&YP_Halo_Type),"Couldn't define YP_Type.");
    MPIDO(MPI_Type_commit(&YP_Halo_Type),"Couldn't commit YP_Type.");
  } /*	if(mpi_ny > 1)	*/

  if (mpi_nz > 1) {
    /*	ZN	*/
    subsizes[0] = xlen;	/*	Widened to include magic corners	*/
    subsizes[1] = ylen;	/*	Widened to include magic corners	*/
    subsizes[2] = 1;
    subsizes[3] = vmax;

    starts[0] = 0;		/*	0 to include magic corner.	*/
    starts[1] = 0;		/*	0 to include magic corner.	*/
    starts[2] = TRI;
    starts[3] = 0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&ZN_Type),"Couldn't define ZN_Type.");
    MPIDO(MPI_Type_commit(&ZN_Type),"Couldn't commit ZN_Type.");
	
    /*	ZN Halo	*/
    starts[2] -= 1;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&ZN_Halo_Type),"Couldn't define ZN_Halo_Type.");
    MPIDO(MPI_Type_commit(&ZN_Halo_Type),"Couldn't commit ZN_Halo_Type.");

    /*	ZP	*/
    starts[0] = 0;		/*	0 to include magic corner.	*/
    starts[1] = 0;		/*	0 to include magic corner.	*/
    starts[2] = (local_zmax-1 + TRI) - local_zmin;
    starts[3] = 0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&ZP_Type),"Couldn't define ZP_Type.");
    MPIDO(MPI_Type_commit(&ZP_Type),"Couldn't commit ZP_Type.");
      
    /*	ZP Halo	*/
    starts[2] += 1;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&ZP_Halo_Type),"Couldn't define ZP_Halo_Type.");
    MPIDO(MPI_Type_commit(&ZP_Halo_Type),"Couldn't commit ZP_Halo_Type.");
  } /*	if(mpi_nz > 1)	*/

  #undef NDIMS
} /* decomp_defineHaloTypes */
