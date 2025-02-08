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

#ifndef _STATE_H_
#define _STATE_H_ 

#ifdef OWN
#undef OWN
#define _OWN
#endif
#include "extern.h"
#include "k_.h"
#if MPI
#include <mpi.h>
#endif /* MPI */
#ifdef _OWN
#undef _OWN
#define OWN
#endif
#include "extern.h"

EXTERN real *New;			/* the main computational grid */

/*! Description of a Device point */
EXTERN real *Geom;			/* the geometry grid */
typedef struct {
  real *u;	/* anchor point in the grid */
  void *X;	/* other device-specific useful load, e.g. diff template weights */
  void *Y;	/* some other useful load, e.g. diff tensor components */
  INT x,y,z;	/* grid coordinates */
} devicePoint;
EXTERN long numTissuePoints;		/* number of tissue points in subdomain */
EXTERN devicePoint *TissuePoints;	/* list of those points */

EXTERN INT xmax,ymax,zmax,vmax;		/* the sizes of the main grid, global */
EXTERN size_t xlen,ylen,zlen;		/* lengths of axes in memory incl halos, local */
EXTERN INT geom_vmax;			/* number of geom variables */
EXTERN INT t;				/* device loops counter */

#if MPI
EXTERN int *ProcessRank;		/* table subdomain -> process rank */
EXTERN int *rank_ix;			/* table process rank -> subdomain x-superindex */
EXTERN int *rank_iy;			/* table process rank -> subdomain y-superindex */
EXTERN int *rank_iz;			/* table process rank -> subdomain z-superindex */
#define NORANK (-1)			/* empty subdomains are "allocated" to this */
#endif /* MPI */

EXTERN INT GEOMETRY_FILE;		/* geometry is defined by file	*/
EXTERN INT GEOMETRY_PGM;		/* geometry is defined by k-pgm */
EXTERN INT GEOMETRY_ON;			/* geometry is defined by file or by pgm */
EXTERN INT ANISOTROPY_ON;		/* anisotropy is on, by file or by choice */

EXTERN size_t vmax_zmax, vmax_zmax_ymax;		/* for main grid index computation */
EXTERN size_t geom_vmax_zmax, geom_vmax_zmax_ymax;	/* for geometry grid index computation */
EXTERN int gx, gy, gz;

EXTERN size_t DX, DY, DZ, DV;		/* main grid index shifts */
EXTERN int dim,ONE,TWO,TRI;		/* used to control dimension */

/*
 * Defaults for a device's Space box.
 * These depend on dimensionality of the problem:
 * long dimensions require boundary points at both ends. 
 */
#define SPACE_DEFAULT_X0 ONE
#define SPACE_DEFAULT_X1 ((int)xmax-(1+ONE))
#define SPACE_DEFAULT_Y0 TWO
#define SPACE_DEFAULT_Y1 ((int)ymax-(1+TWO))
#define SPACE_DEFAULT_Z0 TRI
#define SPACE_DEFAULT_Z1 ((int)zmax-(1+TRI))


/* SUMMARY OF BOX LIMITS OF VARIOUS KINDS 	*/
/* Global: pertaining to the whole grid.  	*/
/* Local: pertaining to this subdomain.   	*/
/* Brutto: including halo points of either kind.*/
/* Netto: excluding halo points of either kind.	*/
/* For seq mode, local==global. .		*/
/* All intervals are lo-inclusive, hi-exclusive.*/
/* Aim at using only these in new code. 	*/
/* TODO: replace legacy boxes with these.	*/
/* VNB 2024/08/13				*/
/* NB this is distinct from device space coords which are lo-inclusive, hi-inclusive !!! */
EXTERN size_t BX0, BX1, BY0, BY1, BZ0, BZ1; /* brutto global  */
EXTERN size_t NX0, NX1, NY0, NY1, NZ0, NZ1; /* netto global  */
EXTERN size_t bx0, bx1, by0, by1, bz0, bz1; /* brutto local */
EXTERN size_t nx0, nx1, ny0, ny1, nz0, nz1; /* netto local  */

/**************************************************************************
 ************************** GEOMETRY SEMANTICS ****************************
 **************************************************************************/

/*	Layers for geometry data.	*/
#define GEOM_STATUS  0
#define GEOM_FIBRE_1 1
#define GEOM_FIBRE_2 2
#define GEOM_FIBRE_3 3

/*	Point statuses.	*/
#define GEOM_VOID   0
#define GEOM_TISSUE 1
#define GEOM_TISSUE2 2

/**************************************************************************
 ************************** MPI ONLY **************************************
 **************************************************************************/
#ifdef MPI
EXTERN int mpi_size;			/* Total number of processes in MPI_COMM_WORLD */
EXTERN int mpi_rank; 			/* Rank of this process among MPI_COMM_WORLD and later among ALL_ACTIVE_PROCS */
EXTERN int num_subdoms;			/* Number of subdomains, == mpi_nx*mpi_ny*mpi_nz <= mpi_size */
EXTERN int num_active_procs;		/* Number of active (allocated AND nonempty) processes, <= num_subdoms */
EXTERN int num_inactive_procs;		/* Number of idle processes, == mpi_size-num_active_procs */
EXTERN int num_dangling_procs;		/* Number of unallocated processes, == mpi_size-num_subdoms */
EXTERN int num_empty_subdoms;		/* Number of subdomains without points, == num_subdoms-num_active_procs */
EXTERN MPI_Comm ALL_ACTIVE_PROCS;	/* Communicator between active processes */
EXTERN int I_AM_IDLE;			/* This process is idle (not allocated by decomp) or empty (without tissue points) */
EXTERN long mpi_nx, mpi_ny, mpi_nz;	/* Number of partitions for each axis. Make them long for convenience of ctlpoint */
EXTERN long mpi_ix, mpi_iy, mpi_iz;	/* Superindices of this process. Make them long for convenience of ctlpoint */
EXTERN size_t *all_xmin, *all_xmax;	/* Local minima and maxima */
EXTERN size_t *all_ymin, *all_ymax;	/*   of all coordinate     */
EXTERN size_t *all_zmin, *all_zmax;	/*   axes' partitions      */
typedef struct {			/* Limits of a subdomain, "netto"; */
  size_t local_xmin, local_xmax;	/*   allowed upper values are local_?max-1 */
  size_t local_ymin, local_ymax;
  size_t local_zmin, local_zmax;
} Subdomain;
EXTERN Subdomain *subdomain_limits;	/* Local minima and maxima for every subdomain */

EXTERN unsigned long *occupy_grid;	/* Number of tissue points in each subdomain */
#define occupy(ix,iy,iz) (occupy_grid[sind((ix),(iy),(iz))]) /* Access to the same via superindices */
EXTERN unsigned long *occupy_rank;	/* Number of tissue points in each subprocess */
EXTERN int laziest_rank;		/* Rank of the least busy of active processes */

/***********************************/
/* This subdomain's details        */

/* Local bounds for this parallel process */
EXTERN size_t local_xmin, local_xmax, local_ymin, local_ymax, local_zmin, local_zmax;
/* Datatypes for sending hyperplanes */
EXTERN MPI_Datatype XN_Type, XP_Type, YN_Type, YP_Type, ZN_Type, ZP_Type;
/* Datatypes for receiving hyperplanes */
EXTERN MPI_Datatype XN_Halo_Type, XP_Halo_Type, YN_Halo_Type, YP_Halo_Type, ZN_Halo_Type, ZP_Halo_Type;

/* Macros for ranks of neighbouring processes */
#define XN_NEIGHBOUR xn_neighbour /* neighbour to LEFT	*/
#define XP_NEIGHBOUR xp_neighbour /* neighbour to RIGHT	*/
#define YN_NEIGHBOUR yn_neighbour /* neighbour ABOVE	*/
#define YP_NEIGHBOUR yp_neighbour /* neighbour BELOW	*/
#define ZN_NEIGHBOUR zn_neighbour /* neighbour in FRONT	*/
#define ZP_NEIGHBOUR zp_neighbour /* neighbour BEHIND	*/

/* Ranks of neighbours, pre-calculated at start time */
EXTERN int xn_neighbour;
EXTERN int xp_neighbour;
EXTERN int yn_neighbour;
EXTERN int yp_neighbour;
EXTERN int zn_neighbour;
EXTERN int zp_neighbour;

/* Map supercoords -> process rank */
int getProcessRank (int ix,int iy,int iz); 

/* end of this subdomain's details */
/***********************************/

/* Get the rank of the process containing a      */
/*    given internal or absolute boundary point. */
int getRankContainingPoint(int x,int y,int z);

/* Macro to compute the index of a variable in New by its coordinates.
 * The MPI version has to be adjusted to match the continuous coordinate 
 * system to the local subdomain stored in New.
 *
 * Takes account of local bounds and whether or not boundary halos exist 
 * on each axis.
 */
/* To do: optimize this */
#define ind(x,y,z,v) ( \
  ((x + ONE) - local_xmin)*vmax_zmax_ymax + \
  ((y + TWO) - local_ymin)*vmax_zmax + \
  ((z + TRI) - local_zmin)*vmax + \
  (v) \
)

/*! Reverse of ind() */
#define xind(ind) (((ind)/vmax_zmax_ymax)+local_xmin-ONE)
#define yind(ind) (((ind)%(vmax_zmax_ymax))+local_ymin-TWO)
#define zind(ind) (((ind)%(vmax_zmax))+local_zmin-TRI)
#define vind(ind) ((ind)%(vmax))

static size_t indfun (size_t x,size_t y,size_t z,int v)
{
  return ((x + ONE) - local_xmin)*vmax_zmax_ymax +	
  ((y + TWO) - local_ymin)*vmax_zmax + 
  ((z + TRI) - local_zmin)*vmax + 
  (v);
}

/* TODO: check this is equivalent to, and which is faster */
/* #define ind(x,y,z,v) ((v)+vmax*((z)+TRI-local_zmin+zmax*((y)+TWO-local_ymin+ymax*((x)+ONE-local_xmin)))) */

/* Macro to compute the index in the gpoints array */
/* used by decompose (decomp.c) and make_gpoints (geometry.c) */
#define gind(x,y,z) (((x)*ymax+(y))*zmax+(z))

/* Macro to compute the index in the subdomains' supergrid, by supercoordinates */
#define sind(ix,iy,iz) (((ix)*mpi_ny+(iy))*mpi_nz+(iz))

/* Macro to compute the index of a variable in Geom by its coordinates.
 * The MPI version has to be adjusted to match the continuous coordinate system to
 * the local subdomain stored in New.
 * Takes account of local bounds and whether or not boundary halos exist on each axis.
 */
/* To do: optimize this by hand */
#define geom_ind(x,y,z,v) ( \
  ((x + ONE) - local_xmin)*geom_vmax_zmax_ymax + \
  ((y + TWO) - local_ymin)*geom_vmax_zmax + \
  ((z + TRI) - local_zmin)*geom_vmax + \
  (v) \
)

#else
/**************************************************************************
 ************************** SEQUENTIAL ONLY *******************************
 **************************************************************************/
#define ind(x,y,z,v) ((x)*vmax_zmax_ymax+(y)*vmax_zmax+(z)*vmax+(v))
#define xind(ind) ((ind)/vmax_zmax_ymax)
#define yind(ind) (((ind)%(vmax_zmax_ymax))/vmax_zmax)
#define zind(ind) (((ind)%(vmax_zmax))/vmax)
#define vind(ind) ((ind)%(vmax))
static size_t indfun (size_t x,size_t y,size_t z,int v)
{
  return ((x)*vmax_zmax_ymax+(y)*vmax_zmax+(z)*vmax+(v));
}

#define geom_ind(x,y,z,v) ((x)*geom_vmax_zmax_ymax+(y)*geom_vmax_zmax+(z)*geom_vmax+(v))

/* and these definition may help avoid some #if(MPI) clauses: */
#define mpi_size 1
#define mpi_rank 0
EXTERN long mpi_nx,mpi_ny,mpi_nz;	/* ==1 */
#define mpi_ix 0L
#define mpi_iy 0L
#define mpi_iz 0L
#define local_xmin (ONE)
#define local_xmax (xmax-ONE)
#define local_ymin (TWO)
#define local_ymax (ymax-TWO)
#define local_zmin (TRI)
#define local_zmax (zmax-TRI)
#endif

/**************************************************************************
 ************************** COMMON ****************************************
 **************************************************************************/

/*! Tests the geometry status (TISSUE/VOID) of a point.
 * Almost duplicates Geom, BUT: is protected against offsite indices 
 */
#define isTissue(x,y,z)	( \
  ( \
    x>=SPACE_DEFAULT_X0 && x<=SPACE_DEFAULT_X1 && \
    y>=SPACE_DEFAULT_Y0 && y<=SPACE_DEFAULT_Y1 && \
    z>=SPACE_DEFAULT_Z0 && z<=SPACE_DEFAULT_Z1 \
  ) && ( \
    geom_vmax==0 || Geom[geom_ind(x,y,z,GEOM_STATUS)] != (real)GEOM_VOID \
  ) \
)
/*! Tests the geometry status (TISSUE/VOID) of a point.
 * Simplified: relies on Geom, assuming indices are within limits
 */
#define istissue(x,y,z)	( \
  geom_vmax==0 || Geom[geom_ind(x,y,z,GEOM_STATUS)] != (real)GEOM_VOID \
)


/*! Shortcuts to find 4D coords of a pointer to a grid point */
#define xof(u) xind((u)-New)
#define yof(u) yind((u)-New)
#define zof(u) zind((u)-New)
#define vof(u) vind((u)-New)

/* Macro for calculating index of screen arrays by 2D pixel coordinates */
#define indxy(x,y) ((x)+xmax*(y))

/* Function defining the state */
int state(char *w);

/* Function destroying the state */
void state_free(void);

#undef EXTERN

#endif /* _STATE_H_ */
