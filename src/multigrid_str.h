/**
 * Copyright (C) (2010-2026) Vadim Biktashev, Irina Biktasheva et al. 
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

/* multigrid_str.h: the fields of a device's STR used by the multigrid code
 * of multigrid_code.h; include inside the STR typedef, after diffweights_str.h,
 * whose fields (conservative, isolated, w0, w1, ...) the multigrid code also uses. */

  int v0;				/* layer for the diffusive field */
  int v1;				/* auxiliary layer */

  char debugname[MAXPATH];		/* Name of the file to write debug information to */
  FILE *debug;				/*	its file descriptor */
  int debugWriter;			/* Nonzero if this process is to write debug information */

  list_t numPoints;		/* Number of tissue points in subdomain. */
  AnisoPoint *points;			/* List of tissue points. ENUMERATION STARTS WITH 1 !!! */
  /* These points use global coordinates */

#if MPI
  MPI_Comm comm; 			/* Communicator for reduction operation. */
  int *active;				/* runHere flags of all processes: only the active ones */
					/* take part in the multigrid's halo exchanges          */
  int root; 				/* Rank of device root in ALL_ACTIVE_PROCS */
#endif

  /* The operator: L=div(D grad) discretized as in diff if tau==0, otherwise I-tau*L, */
  /* e.g. tau=theta*ht for an implicit step of diffusion. Set by the device before   */
  /* mg_connect(). STR is calloc'd, so tau==0 unless the device sets it.             */
  real tau;

  /*****************************************************************/
  /* The singular problem: no-flux at all borders, no Dirichlet values */
  int singular;				/* nonzero if no point is linked to values beyond the device space */
  real fmean;				/* plain GS: mean of the RHS, subtracted to make it compatible */
  /* Normalization of the solution of the singular problem: zero mean, or if pin=1, */
  /* shifted by a constant so that u(xpin,ypin,zpin)=upin. NB this is applied to   */
  /* the solution, not imposed as an equation, so does not affect convergence.     */
  int pin;				/* 0/1 flag */
  int xpin, ypin, zpin;			/* global coords of the pinning point */
  real upin;				/* pinning value */
  list_t pinpoint;			/* its number in S->points, 0 if not in this subdomain */

  /***********************************/
  /* Base iteration parameters (TBE) */
  /* int vp;				/\* Layer to remember the previous solution - redundant if FMG used *\/ */
  /* int started;				/\* Should not use the previous solution before it has been remembered *\/ */

  real tolerance;			/* C0 norm of change in one iteration that is considered small enough */
  real resnorm;				/* Residual norm achieved so far */
  int maxiter;				/* Max num of iterations after multigrid precondition */
  int preiter;				/* Max num of iter before coarsening */
  int postiter;				/* Max num of iter after coarsening */
  int coarseiter;			/* Max num of iter at the coarsest level (0: preiter+postiter as elsewhere) */
  real coarsereduce;			/* ... iterate there until the residual is reduced by this factor */
 

  /*******************/
  /* Multigrid stuff */

  Space s;				/* Device space, needed by a sub-function */
  Space *sp;				/* The device's own Space, for in_space(): its where-code uses its fields */
  int upper_level;			/* The number of levels is (upper_level+1) */

  /* Iteration parameters */
  real delta;				/* Tolerance factor per multigrid level */
  real damping;				/* Damping for smoothing */
  real *tol;				/* The level-dependent tolerances */
  real *rsn;				/* The level-dependent residual norms */
  int FMGlevel;				/* Depth of the Full Multigrid ("nested iterations") */
  int vcycles;				/* Number of V-cycles */
  int *iter;				/* Continuous iteration counters for each level */
  int *preIter;				/* Max number of iteration for each level */
  int *postIter;			/* Max number of iteration for each level */
  char smoother[32];			/* Name of the smoother */
  smoother_t smootherchoice;		/* Choice of the smoother; currently only Jacobi/GS */
  char transfer[32];			/* Name of the transfer method */
  transfer_t transferchoice;		/* Choice of the transfer method: correct/replace */
  int direction;			/* current GS direction, +1 or -1, when alternate */
  int alternate;			/* If nonzero, GS direction alternates (DOES THIS WORK???) */

  /* Grid structure */
  list_t *nPoints;			/* Number of non-empty points in each grid */
  AnisoPoint **gridPoints;		/* List of non-empty points in each grid. ENUMERATION STARTS WITH 1 !!! */
  /* AnisoPoint gridPoints[numLevels][nPoints+1], these use local netto coordinates */
  coord_t *nx, *ny, *nz;		/* Netto dimensions of the grids (without halos) */
  coord_t *xlen, *ylen, *zlen;		/* Brutto dimensions of the grids (with halos) */
  coord_t *levelsize;			/* ... convenient to have the product of the three */
  coord_t *xmin, *ymin, *zmin;		/* Netto minimals in brutto (local) coordinates. Not used for now */
  coord_t *xmax, *ymax, *zmax;		/* Netto maximals in brutto (local) coordinates (actually accessible) */
  index_t **pointIndex;			/* 1D indexes of points in the lists */
  /* index_t pointIndex[numLevels][nPoints+1] */
  /* The grids are level arrays of "flattened" 1D arrays representing 1D, 2D or 3D meshes */
  real **ugrid;				/* The solution grid. */
  real **rgrid;				/* The residual grid. */
#define ggrid rgrid			/* At start time, same arrays used as "geometry" grid (nonzeros at tissue points) */
  real **fgrid;				/* The RHS (free terms, "force") grid */
  /* TODO: allocate only point lists, to save memory? */

  /* System matrices */
  neighbour_t numNeighbours;		/* max number of neighbours: depends on dimensionality */
  sindex_t **d;				/* 1D index displacements for neighbours */
  /* sindex_t d [numLevels][numNeihbours] */
  list_t numBClinks;			/* Links to values outside this device's space */
  bclink *BClinks;			/*    e.g. for non-homogeneous boundary conditions */

  /* Transition matrices */
  contr_t numContributors;		/* max number of contributors: depends on dimensionality */
  sindex_t **Dc;			/* 1D index displacements for contributors */
  /* sindex_t D[numLevels][numContributors] */
  index_t ***contributors;		/* 1D indices of such contributor nodes for each higher level point */
  /* index_t contributors[numLevels][nPoints+1][numContributorsMax] */
  real ***P;				/* P matrix entries */
  real ***R;				/* R matrix entries */
  /* real P/R[numLevels][nPoints+1][numContributorsMax] */
  /* The same P matrix by rows, for prolongation by the fine points:  */
  /* point i of level l-1 gets parN[l][i] terms parW[l][8*i+n]*u_l[parIdx[l][8*i+n]] */
  int **parN;
  index_t **parIdx;
  real **parW;

  /* Coarsening geometry, for each dimension d (0=x,1=y,2=z) and level l,	*/
  /* in global coordinates measured in the mesh steps of that level.		*/
  int *off[3];				/* l>0: coarse local k <-> fine local 2k+off[d][l] at level l-1 */
  long *g0[3];				/* global coordinate of the local netto point 0 */
  long *glo[3], *ghi[3];		/* global limits of the device space (images of the level-0 ones) */
  long *anchor[3];			/* level l+1 points are at level-l global positions of this parity */
  int dlo[3], dhi[3];			/* nonzero if the low/high face of the device space has Dirichlet links */

#if MPI
  /* Halo exchange MPI datatypes for multigrid levels */
  MPI_Datatype *xn_type, *xp_type;
  MPI_Datatype *yn_type, *yp_type;
  MPI_Datatype *zn_type, *zp_type;
  MPI_Datatype *xn_halo_type, *xp_halo_type;
  MPI_Datatype *yn_halo_type, *yp_halo_type;
  MPI_Datatype *zn_halo_type, *zp_halo_type;
#endif

  /* Profiling stuff */

  char profilename[MAXPATH];		/* Name of the file to write profile information to */
  FILE *profile;				/*	its file descriptor */
  int profileWriter;			/* Nonzero if this process is to write profile information */

  int profile_lastlevel;
  double profile_lasttime;
  double *profile_totalspent;	/* [TOTALSPENTNUM] */
