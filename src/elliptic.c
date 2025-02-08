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
 * This device solves the elliptic equation of the bidomain model 
 * \[ \nabla( \hat{D} \nabla \phi_e ) = f \]
 * Take the source field $f$ prepared in layer v0
 * and calculate the $\phi_e$ field in layer v1.
 * Implemented only for the anisotropic case for now. 
 *
 * Main features: Full Multigrid with vertex-centered
 * restrition/prolongation operators with bi/tri-linear interpolation,
 * and (multicoloured) Gauss-Seidel or Jacobi smoother.
 * Parallel version seems to work in 2D but there are still some issues in 3D.
 * 
 */

#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "qpp.h"

/*==============================================*/
/* Important defined/derived types              */

/* Choice of the smoother */
typedef enum {
  Jacobi,
  GS
} smoother_t;

/* Choice of inter-grid transfer method */
typedef enum {
  correct,
  replace
} transfer_t;

/* Choice of the exact analytical solution.  */
/* Obviously, this is only for debugging !!! */
typedef enum {
  zero,
  sinv,
  ones,
  zfk
} exact_t;

/* Enumerated neighbours.    */
/* 1D precede 2D precede 3D. */
typedef enum {
  ooo,			/* 0D */
  poo, moo,		/* 1D 3-point stencil */
  opo, omo,		/* 2D 5-point stencil */
  ppo, pmo, mpo, mmo,	/* 2D 9-point stencil */
  oop, oom,		/* 3D 7-point stencil */
  pop, pom, mop, mom,	/* 3D 19-point stencil */
  opp, opm, omp, omm,	/* 3D 19-point stencil */
  ppp, ppm, pmp, pmm,	/* 3D 27-point stencil */
  mpp, mpm, mmp, mmm,	/* 3D 27-point stencil */
  maxnb
} neighbour_t;

/* "Systems of coordinates":									*/
typedef size_t coord_t;		/* 3D local (netto) x,y,z coordinates: x=0 <=> local_xmin etc   */
typedef ssize_t scoord_t;	/* Signed version of the same: for displacements		*/
typedef size_t index_t;		/* Flattened 1D array index: 0..xmax_ymax_zmax_vmax             */
typedef ssize_t sindex_t;	/* Signed version of the same: for displacements	        */
typedef size_t list_t;		/* Number of a point in the point list: numPoints         */

/* Conversions between coordinates:								*/
/* gridIndex(level,x,y,z):        coord_t^3 |-> index_t						*/
/*: (z+S->zmin[level]+S->nz[level]*(y+S->ymin[level]+S->ny[level]*(x+S->xmin[level]))) - same only longer  */
/* S->gridPoints[level][i].x etc: list_t    |-> coord_t					        */
/* pointIndex[level][i]:          list_t    |-> index_t						*/
/* gridX(level,ind) etc:          index_t   |-> coord_t				        	*/
/* pointNumber(level,ind):        index_t   |-> list_t						*/

/* Enumerated contributors.  */
/* 1D precede 2D precede 3D. */
typedef enum {
  OOO,			/* 0D */
  POO, MOO,		/* 1D */
  OPO, OMO,		
  PPO, PMO, MPO, MMO,	/* 2D */
  OOP, OOM,		
  POP, POM, MOP, MOM,	
  OPP, OPM, OMP, OMM,	
  PPP, PPM, PMP, PMM,
  MPP, MPM, MMP, MMM,	/* 3D */
  maxcn
} contr_t;

/* Enumerator of colors (as in multicolor GS) */
typedef int color_t; /* 0..7 actually */

/* An entry of the list of points, RMcF style */
typedef struct {
  real *K;				/* Neibours connection strengths (entries of A matrices) */
  coord_t x,y,z;			/* Coordinates */
} AnisoPoint;      

/* A link to value outside space */
typedef struct {
  index_t to;				/* Target point in the grid  */
  real *from;				/* Pointer to outside value  */
  real K;				/* Connection strength       */
} bclink;


/*****************************************************************************/
/* The device parameters structure - in this device a bit longer than usual! */
typedef struct {                

  /***********************/
  /* The usual diff stuff */

  real Dpar;				/* Parallel diffusion coefficient. ANISO only */
  real Dtrans;				/* Transverse diffusion coefficient. ANISO only. */
  real hx;				/* Space step */
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
  int root; 				/* Rank of device root in ALL_ACTIVE_PROCS */
#endif

  /*****************************************************/
  /* The pinning condition, needed for Neumann problem */
  int pin;				/* 0/1 flag */
  int xpin, ypin, zpin;			/* grid coord of pinning point */
  real upin;				/* pinning value */
#define PINCOEFF (10.0)			/* to make pinning equation commensurate with others */

  /***********************************/
  /* Base iteration parameters (TBE) */
  /* int vp;				/\* Layer to remember the previous solution - redundant if FMG used *\/ */
  /* int started;				/\* Should not use the previous solution before it has been remembered *\/ */

  real tolerance;			/* C0 norm of change in one iteration that is considered small enough */
  real resnorm;				/* Residual norm achieved so far */
  int maxiter;				/* Max num of iterations after multigrid precondition */
  int preiter;				/* Max num of iter before coarsening */
  int postiter;				/* Max num of iter after coarsening */
 

  /*******************/
  /* Multigrid stuff */

  Space s;				/* Device space, needed by a sub-function */
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
#define Pgrid ugrid			/* At start time, same arrays used to calc row sums of P matrices */
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
  sindex_t **D;				/* 1D index displacements for contributors */
  /* sindex_t D[numLevels][numContributors] */
  index_t ***contributors;		/* 1D indices of such contributor nodes for each higher level point */
  /* index_t contributors[numLevels][nPoints+1][numContributorsMax] */
  real ***P;				/* P matrix entries */
  real ***R;				/* R matrix entries */
  /* real P/R[numLevels][nPoints+1][numContributorsMax] */

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

} STR;
/*****************************************************************************/

/*==============================================*/
/* Local static objects                         */

/* Non-normalized stencils for prolongation/restriction operators */
static real wgt[maxcn]={
  8,			/* OOO, 		0D */
  4, 4, 		/* POO, MOO,		1D */
  4, 4,			/* OPO, OMO,		   */
  2, 2, 2, 2,		/* PPO, PMO, MPO, MMO,	2D */
  4, 4,			/* OOP, OOM,		   */
  2, 2, 2, 2, 		/* POP, POM, MOP, MOM,	   */
  2, 2, 2, 2,		/* OPP, OPM, OMP, OMM,	   */
  1, 1, 1, 1,		/* PPP, PPM, PMP, PMM,     */
  1, 1, 1, 1,		/* MPP, MPM, MMP, MMM,	3D */
};


/* Coord displacements of neighbours */
static scoord_t ndx[maxnb]={
  0  ,			/* ooo, */
  1  , -1 ,		/* poo, moo, */
  0  , 0  ,		/* opo, omo, */
  1  , 1  , -1 , -1 ,	/* ppo, pmo, mpo, mmo, */
  0  , 0  ,		/* oop, oom, */
  1  , 1  , -1 , -1 ,	/* pop, pom, mop, mom, */
  0  , 0  , 0  , 0  ,	/* opp, opm, omp, omm, */
  1  , 1  , 1  , 1  ,	/* ppp, ppm, pmp, pmm, */
  -1 , -1 , -1 , -1  	/* mpp, mpm, mmp, mmm  */
};
static scoord_t ndy[maxnb]={
  0  ,			/* ooo, */
  0  , 0  ,		/* poo, moo, */
  1  , -1 ,		/* opo, omo, */
  1  , -1 , 1  , -1 ,	/* ppo, pmo, mpo, mmo, */
  0  , 0  ,		/* oop, oom, */
  0  , 0  ,  0 ,  0 ,	/* pop, pom, mop, mom, */
  1  , 1  , -1 , -1 ,	/* opp, opm, omp, omm, */
  1  , 1  , -1 , -1 ,	/* ppp, ppm, pmp, pmm, */
  1  , 1  , -1 , -1  	/* mpp, mpm, mmp, mmm  */
};
static scoord_t ndz[maxnb]={
  0  ,			/* ooo, */
  0  , 0  ,		/* poo, moo, */
  0  , 0  ,		/* opo, omo, */
  0  , 0  , 0  , 0  ,	/* ppo, pmo, mpo, mmo, */
  1  , -1 ,		/* oop, oom, */
  1  , -1 , 1  , -1 ,	/* pop, pom, mop, mom, */
  1  , -1 , 1  , -1 ,	/* opp, opm, omp, omm, */
  1  , -1 , 1  , -1 ,	/* ppp, ppm, pmp, pmm, */
  1  , -1 , 1  , -1  	/* mpp, mpm, mmp, mmm  */
};

/* Coord displacements of contributors */
static scoord_t cdx[maxcn]={
  0  ,			/* OOO, */
  1  , -1 ,		/* POO, MOO, */
  0  , 0  ,		/* OPO, OMO, */
  1  , 1  , -1 , -1 ,	/* PPO, PMO, MPO, MMO, */
  0  , 0  ,		/* OOP, OOM, */
  1  , 1  , -1 , -1 ,	/* POP, POM, MOP, MOM, */
  0  , 0  , 0  , 0  ,	/* OPP, OPM, OMP, OMM, */
  1  , 1  , 1  , 1  ,	/* PPP, PPM, PMP, PMM, */
  -1 , -1 , -1 , -1  	/* MPP, MPM, MMP, MMM  */
};
static scoord_t cdy[maxcn]={
  0  ,			/* OOO, */
  0  , 0  ,		/* POO, MOO, */
  1  , -1 ,		/* OPO, OMO, */
  1  , -1 , 1  , -1 ,	/* PPO, PMO, MPO, MMO, */
  0  , 0  ,		/* OOP, OOM, */
  0  , 0  ,  0 ,  0 ,	/* POP, POM, MOP, MOM, */
  1  , 1  , -1 , -1 ,	/* OPP, OPM, OMP, OMM, */
  1  , 1  , -1 , -1 ,	/* PPP, PPM, PMP, PMM, */
  1  , 1  , -1 , -1  	/* MPP, MPM, MMP, MMM  */
};
static scoord_t cdz[maxcn]={
  0  ,			/* OOO, */
  0  , 0  ,		/* POO, MOO, */
  0  , 0  ,		/* OPO, OMO, */
  0  , 0  , 0  , 0  ,	/* PPO, PMO, MPO, MMO, */
  1  , -1 ,		/* OOP, OOM, */
  1  , -1 , 1  , -1 ,	/* POP, POM, MOP, MOM, */
  1  , -1 , 1  , -1 ,	/* OPP, OPM, OMP, OMM, */
  1  , -1 , 1  , -1 ,	/* PPP, PPM, PMP, PMM, */
  1  , -1 , 1  , -1  	/* MPP, MPM, MMP, MMM  */
};

/*==============================================*/
/* Local functions                              */

/* Make list of points, with their connections */
static int list_points(STR *S,Space s);

/* Prepare multigrid stuff for given level (and above, recursively) */
static int make_level(int level,STR *S);

/* Straight iterations (outside multigrid construct) */
static real iterate_base(int it, STR *S);

/* Transfer data between New and the base grid */
static void to_level0(STR *S);
static void from_level0(STR *S);

/* Do multigrid computations for given level, recursively */
static int do_level(int level,STR *S);

/* Transfer the residual from level to level+1 */
static int up_level(int level,STR *S);

/* Transfer the results from level+1 to level */
static int down_level(int level,STR *S);

/* Do the higher level =  up + there + down */
static int higher_level(int level,STR *S);

/* Do one iteration, or calc residuals if relax=0 */
static real iterate_level(int level,real relax,STR *S);

/* Just calculate the residuals - for debugging only */
static real residual_level(int level,STR *S);

#if MPI
/* Swap the halos at the given level */
static int halo_level(int level,real *grid,STR *S);
/* Define the halo MPI datatypes */
static void define_halo_level(int level,STR *S);
#endif

/*-------------------------------------------------------*/
/* Start profiling episode */
static void profile_start(STR *S);
/* Profiling milestone */
static void profile_level(int level,STR *S);
/* Index of the projection operator totalspent */
#define UPLEVEL (S->upper_level+1)
/* Index of the reduction operator totalspent */
#define DOWNLEVEL (S->upper_level+2)
/* Index of the to-grid conversion totalspent */
#define TOGRID (S->upper_level+3)
/* Index of the from-grid conversion totalspent */
#define FROMGRID (S->upper_level+4)
/* Index for "everything else" totalspent */
#define ALLELSE (S->upper_level+5)
#define TOTALSPENTNUM (S->upper_level+6)
/* End profiling episode and report results */
static void profile_report(STR *S);
/*-------------------------------------------------------*/

/* Coordinates converters */
/* These converters understand (x,y,z) as local netto coordinates */
/* TODO: optimize this ? */
#define gridIndex(level,x,y,z) ((z)+TRI+S->zlen[level]*((y)+TWO+S->ylen[level]*((x)+ONE))) 
#define gridZ(level,index) ((index)                                    %(S->zlen[level])-TRI)
#define gridY(level,index) ((index)/((S->zlen[level])                 )%(S->ylen[level])-TWO)
#define gridX(level,index) ((index)/((S->zlen[level])*(S->ylen[level]))                 -ONE)

/* 0 return means this point is not in the list, ie void */
/* #define pointNumber(level,index) (S->ggrid[level][index]) */
#define pointNumber(level,index) (index<0?0:index>=S->levelsize[level]?0:S->ggrid[level][index]) 

/* Macros to make higher-level nodes coords from lower-level coords */
/* This way it is used only for the dimensions of the coarser grid */
#define coarsex(xl) ((xl+2)/2-1)
#define coarsey(yl) ((yl+2)/2-1)
#define coarsez(zl) ((zl+2)/2-1)
/* And vice-versa */
#define finex(xh) (2*(xh))
#define finey(xh) (2*(yh))
#define finez(xh) (2*(zh))

#undef DEBUG
#define DEBUG(...) if (S->debug&&S->debugWriter) {fprintf(S->debug,__VA_ARGS__);fflush(S->debug);}
#define PROFILE(...) if (S->profile&&S->profileWriter) {profile_level(__VA_ARGS__,S);}

/* END OF DEFINITIONS                                    */
/*-------------------------------------------------------*/

/******************/
RUN_HEAD(elliptic)
{
  DEVICE_CONST(real,hx)
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(real,tolerance)
  DEVICE_CONST(int,FMGlevel)
  DEVICE_CONST(int,vcycles)
  DEVICE_CONST(int,maxiter)
  /* DEVICE_CONST(int,vp) */
  /* DEVICE_VAR(int,started) */
  DEVICE_CONST(int, debugWriter)
  DEVICE_CONST(FILE *,debug)

  DEVICE_ARRAY(AnisoPoint, points)
  DEVICE_CONST(list_t, numPoints)

#if MPI
  DEVICE_CONST(MPI_Comm, comm)
  DEVICE_CONST(int, root)
#endif

  AnisoPoint P;
  coord_t x,y,z;
  list_t i;
  real *u, *ucurr, *uprev;
  real unew, f, sum, absres, rsn0, rsn1, rsn2;
  int it;

  int level=0;
  
  DEBUG("\nelliptic #%d: t=%ld\n", mpi_rank, (long)t);
  profile_start(S);
  
  if (!ANISOTROPY_ON) EXPECTED_ERROR("this device is implemented only for the anisotropic case\n");
  
  /* Apply the linear predictor */
  /* if (vp>=0) { */
  /*   for (i=1;i<=numPoints;i++) { /\* ENUMERATION STARTS WITH 1 !!! *\/ */
  /*     P=points[i]; */
  /*     x = P.x; */
  /*     y = P.y; */
  /*     z = P.z; */
  /*     ucurr=New+ind(x,y,z,v1); */
  /*     uprev=New+ind(x,y,z,vp); */
  /*     if (*started) { (*ucurr) += (*ucurr)-(*uprev);} */
  /*     (*uprev)=(*ucurr); */
  /*   } */
  /*   *started=1; */
  /* } */

  PROFILE(ALLELSE);
  /* Apply the multi-grid step */
  if (S->upper_level>=0) {
    static int cyc;
    static int level=0;

    for (level=0;level<=S->upper_level;level++) 
      S->iter[level]=0;				/* Reset iteration counters */

    PROFILE(TOGRID);
    to_level0(S);				/* New -> grid[0] conversion */
    PROFILE(ALLELSE);

    for (level=0;level<FMGlevel;level++)	/* Setup RHS for the FMG iterations */
      if NOT(up_level(level,S)) return 0;

    /* The Full Multi-Grid (FMG) iterations: the preconditioner */
    for (level=FMGlevel;level>=0;level--) {	/* The nested iterations loop */
      DEBUG("FMG loop from level %d\n",level);
      for (cyc=0;cyc<vcycles;cyc++) {		/* The recurvsive V-cycles */
	DEBUG("V-cycle %d\n",cyc);
	if NOT(do_level(level,S)) return 0;	/* This is a V-cycle from this level    */
	DEBUG("V-cycle %d done, rsn[%d]=%g vs tol=%g\n",
	      cyc,level,S->rsn[level],S->tol[level]);
	if (S->rsn[level]<tolerance) break;	/* Check if this will suffice           */
      }
      if (level>0)				/* Interpolate for next nested iteration */
	down_level(level-1,S);
    }

    /* "Post-conditioning" finalizing iterations */
    for (it=0;it<maxiter && S->rsn[0]>=tolerance;it++)
	iterate_level(0,1,S);
    
    DEBUG("final result: level=0 iter=%ld resnorm=%lg\n", (long)(S->iter[0]), (double)S->rsn[0]);
    if (debug) {
      int l, i;
      real factor=1.0;
      real total=0.0;
      DEBUG("more accurately: iter");
      for (l=0;l<=S->upper_level;l++) {
	i=S->iter[l];
	DEBUG("%c%3d",l?'+':'=',i);
	total+=i*factor;
	factor*=0.25;
      }
      DEBUG("=%.3g resnorm=%lg\n", total, (double)S->rsn[0]);
    }

    PROFILE(FROMGRID);
    from_level0(S);				/* New -> grid[0] conversion */
    PROFILE(ALLELSE);

  } 
  /* Straight GS iterations by the old code. */
  /* This is to be phased out when all done? */
  else {
    for (it=0;it<maxiter;it++)
      if (tolerance>iterate_base(it,S))
	break;
    DEBUG("result: iter=%ld resnorm=%lg\n", (long)it, (double)S->resnorm);
  }

  profile_report(S);
}
RUN_TAIL(elliptic)

  /*******************/
DESTROY_HEAD(elliptic) {
  FREE(S->points);
  SAFE_CLOSE(S->debug);
  SAFE_CLOSE(S->profile);
  /* TODO: free all multigrid objects ... */
} DESTROY_TAIL(elliptic)
  
CREATE_HEAD(elliptic)
{
  DEVICE_REQUIRES_SYNC;
  int level=-1;

  if (! ANISOTROPY_ON) EXPECTED_ERROR("elliptic is only implemented for the anisotropic case\n");
  
  ACCEPTF(debug,"wt","");
  ACCEPTF(profile,"wt","");
#if MPI
  if (S->profile) MESSAGE("\n/* profiling not implemented for MPI; will proceed without profiling */");
#endif	
  ACCEPTR(Dpar,RNONE,0.,RNONE);
  ACCEPTR(Dtrans,RNONE,0.,RNONE);
  if (find_key("D=",w)) 
    MESSAGE("The isotropic diffusion coefficient 'D' is not needed when anisotropy is active."
	    "This parameter will be ignored.\n");
  ACCEPTR(hx,RNONE,0.,RNONE);
  ACCEPTR(tolerance,RNONE,0.,RNONE);
  ACCEPTI(maxiter,INONE,1,INONE);
  ACCEPTI(pin,0,0,1);
  if (S->pin) {
    ACCEPTI(xpin,INONE,dev->s.x0,dev->s.x1);
    ACCEPTI(ypin,INONE,dev->s.y0,dev->s.y1);
    ACCEPTI(zpin,INONE,dev->s.z0,dev->s.z1);
    ACCEPTR(upin,RNONE,RNONE,RNONE);
  }
  /* ACCEPTI(vp,-1,-1,vmax-1); */
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 != v0 );
  /* ASSERT( v0 != vp ); */
  /* ASSERT( v1 != vp ); */
  S->s=dev->s;
  
  /* S->started=0; */

  if (ANISOTROPY_ON) /* always true for now, but retain this if for future development */
    list_points(S,dev->s);
  
#if MPI
  /* We need one process that is guaranteed to be engaged in this device */
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    ABORT("Could not create communicator.\n");
  /* Identifying the root process allows only one process to write debug output. */
  S->debugWriter = (mpi_rank == S->root);
  S->profileWriter = (mpi_rank == S->root);
#else
  S->debugWriter = 1;
  S->profileWriter = 1;
#endif /* MPI */

  ACCEPTI(upper_level,-1,-1,xmax);
  /* upper_level<0 means no multigrid, just straightforward GS iterations by old code. */
  if (S->upper_level>=0) {
    int numLevels=S->upper_level+1;	/* the number of levels: from 0 to upper_level incl */
    int l; 				/* level counter */

    int level=0;

    ACCEPTS(smoother,"GS");
    if (0);
    #define CASE(a) else if (0==stricmp(S->smoother,a)) 
    CASE("GS") {S->smootherchoice=GS; ACCEPTI(alternate,0,0,1);}
    CASE("Jacobi") S->smootherchoice=Jacobi;
    #undef CASE
    else EXPECTED_ERROR("\nUnknown smoother name %s",S->smoother);

    ACCEPTS(transfer,"correct");
    if (0);
    #define CASE(a) else if (0==stricmp(S->transfer,a)) 
    CASE("correct") S->transferchoice=correct;
    CASE("replace") S->transferchoice=replace;
    #undef CASE
    else EXPECTED_ERROR("\nUnknown transfer method name %s",S->transfer);


    ACCEPTR(damping,1.,0.,RNONE);
    ACCEPTR(delta,1.,0.,1.);
    ACCEPTI(FMGlevel,S->upper_level,0,S->upper_level);
    ACCEPTI(vcycles,1,1,INONE);
    ACCEPTI(preiter,S->maxiter,0,INONE);
    ACCEPTI(postiter,S->maxiter,0,INONE);

    /* Allocate level arrays */
    CALLOC(S->tol,numLevels,sizeof(real));
    CALLOC(S->rsn,numLevels,sizeof(real));
    CALLOC(S->iter,numLevels,sizeof(int));
    CALLOC(S->preIter,numLevels,sizeof(int));
    CALLOC(S->postIter,numLevels,sizeof(int));
    CALLOC(S->nPoints,numLevels,sizeof(list_t));
    CALLOC(S->gridPoints,numLevels,sizeof(AnisoPoint *));
    CALLOC(S->contributors,numLevels,sizeof(index_t **));
    CALLOC(S->P,numLevels,sizeof(real **));
    CALLOC(S->R,numLevels,sizeof(real **));

    CALLOC(S->nx,numLevels,sizeof(coord_t));
    CALLOC(S->ny,numLevels,sizeof(coord_t));
    CALLOC(S->nz,numLevels,sizeof(size_t));

    CALLOC(S->xlen,numLevels,sizeof(coord_t));
    CALLOC(S->ylen,numLevels,sizeof(coord_t));
    CALLOC(S->zlen,numLevels,sizeof(coord_t));
    CALLOC(S->levelsize,numLevels,sizeof(coord_t));

    CALLOC(S->xmax,numLevels,sizeof(coord_t));
    CALLOC(S->ymax,numLevels,sizeof(coord_t));
    CALLOC(S->zmax,numLevels,sizeof(coord_t));

    CALLOC(S->xmin,numLevels,sizeof(coord_t)); /* these ones are not actually used for now */
    CALLOC(S->ymin,numLevels,sizeof(coord_t)); /* these ones are not actually used for now */
    CALLOC(S->zmin,numLevels,sizeof(coord_t)); /* these ones are not actually used for now */

    CALLOC(S->pointIndex,numLevels,sizeof(list_t *));

    CALLOC(S->ugrid,numLevels,sizeof(real *));
    CALLOC(S->rgrid,numLevels,sizeof(real *));
    CALLOC(S->fgrid,numLevels,sizeof(real *));
    CALLOC(S->d,numLevels,sizeof(sindex_t *));
    CALLOC(S->D,numLevels,sizeof(sindex_t *));
    
    /* Set tolerances iteration numbers at all levels */
    for (l=0;l<numLevels;l++) {
      ACCEPTRE("tol%d",S->tol,l,l?(S->tol[l-1]*S->delta):(S->tolerance),0.,RNONE);
      ACCEPTIE("preiter%d",S->preIter,l,S->preiter,0,INONE);
      ACCEPTIE("postiter%d",S->postIter,l,S->postiter,0,INONE);
    }

#if MPI
    CALLOC(S->xn_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->xp_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->yn_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->yp_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->zn_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->zp_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->xn_halo_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->xp_halo_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->yn_halo_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->yp_halo_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->zn_halo_type,numLevels,sizeof(MPI_Datatype));
    CALLOC(S->zp_halo_type,numLevels,sizeof(MPI_Datatype));
#endif

    /* Now make level 0 which recursively entails all others */
    make_level(0,S);

    /* Allocate profiling array */
    if (S->profile)
      CALLOC(S->profile_totalspent,TOTALSPENTNUM,sizeof(double));

    /* For debug purposes: check that R, P and A-operators map all-ones to all-1 (R,P) and all-0 (A) */
    if (Verbose) {
      transfer_t transfer_safe=S->transferchoice;
      list_t il, ih;
      index_t indl, indh;
      int level;
      real err,maxerr;
      S->transferchoice=replace;
#define One (1.0)
#define Zero (0.0)
      for (level=0;level<=S->upper_level;level++) {
	/*---------------------------------------------*/
	if (level>0) { /* Checking R: level-1 -> level maps 1->1 */
	  for (il=1;il<=S->nPoints[level-1];il++) 	/* set 1s thru lower level */
	    S->fgrid[level-1][ S->pointIndex[level-1][il] ]=One;
	  up_level(level-1,S);				/* apply R operator */
	  maxerr=0;
	  for (ih=1;ih<=S->nPoints[level];ih++){ 	/* check thru higher level */
	    err=fabs(S->fgrid[level][ S->pointIndex[level][ih] ]-One);
	    if (err>maxerr) maxerr=err;
	  } 
	  MESSAGE("\n/* R matrix (%d->%d) is 1:1 to within %g. */",level-1,level,maxerr);
	} /* if level */
	
	/*---------------------------------------------*/
	if (level>0) { /* Checking P: level -> level-1 maps 1->1 */
	  index_t iwst;
	  for (ih=1;ih<=S->nPoints[level];ih++) 	/* set 1s thru higher level */
	    S->ugrid[level][ S->pointIndex[level][ih] ]=One;
	  down_level(level-1,S);			/* apply P operator */
	  maxerr=0;
	  for (il=1;il<=S->nPoints[level-1];il++){ 	/* check thru lower level */
	    indl=S->pointIndex[level-1][il];
	    err=fabs(S->ugrid[level-1][indl]-One);
	    if (err>maxerr) { 
	      maxerr=err;
	      iwst=indl;
	    }
	  } 
	  MESSAGE("\n/* P matrix (%d->%d) is 1:1 to within %g. */",level,level-1,maxerr);
	  if (maxerr>1.e-13) 
	    MESSAGE("\n/*   worst offender: %d (%d,%d,%d)",
		   (int)iwst,(int)gridX(level-1,iwst),
		   (int)gridY(level-1,iwst),(int)gridZ(level-1,iwst));
	} /* if level */
	
	/*--------------------------------------------*/
	if (level>=0) { /* Checking A: level -> level maps 1->0 */
	  for (ih=1;ih<=S->nPoints[level];ih++) {
	    indh=S->pointIndex[level][ih];
	    S->ugrid[level][indh]=One;			/* set 1 thru u-grid */
	    S->fgrid[level][indh]=Zero;			/* and 0 thru f-grid */
	  }
	  residual_level(level,S);			/* apply A operator */
	  maxerr=0;
	  for (ih=1;ih<=S->nPoints[level];ih++) {
	    indh=S->pointIndex[level][ih];
	    err=fabs(S->rgrid[level][indh]-Zero);		/* check thru r-grid */
	    if (err>maxerr) maxerr=err;
	  }
	} /* if level */
	MESSAGE("\n/* A matrix (%d->%d) is 1:1 to within %g. */",level,level,maxerr);
      } /* for level */
#undef One
#undef Zero
      S->transferchoice=transfer_safe;
    } /* if Verbose debugging of 1-transferrals */
  } /* if upper_level >= 0 */
}
CREATE_TAIL(elliptic,1)
/*=============================================*/

/* Make the list of tissue points with all their connections */
/* as defined by the anisotropic geometry, */
/* by the same code as in diff.c */
static int list_points(STR *S,Space s) {
  DEVICE_CONST(real,Dpar)
  DEVICE_CONST(real,Dtrans)
  DEVICE_CONST(real,hx)
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  coord_t x,y,z,x_shift,y_shift,z_shift;	/* Cartesian grid coords */
  int j,k;					/* dimensions enumerators */
  list_t i;					/* list counter */
  index_t indh;
  color_t color, the_color;
    
  /* Array dimensions are 4 to make indices match paper */
  real f[4];				/* Fibre direction vector */
  real D[4][4];			/* Diffusion tensor */
  real dD[4][4];			/* Derivatives of diffusion tensor */
  real c[4];				/* Derivatives of diffusion coefficients */
  
  /* Neighbours' enumerator */
  neighbour_t nb;
  /* Neighbours' vectors */
  coord_t nx,ny,nz; /* Coords of negative neighbour */
  coord_t px,py,pz; /* Coords of positive neighbour */
  
  /* The current point's structure */
  AnisoPoint P;
  
  /***************************************************/
  
  /* Count tissue points and allocate array of AnisoPoints */
  S->numPoints = 0;
  for(x=s.x0;x<=s.x1;x++){
    for(y=s.y0;y<=s.y1;y++){
      for(z=s.z0;z<=s.z1;z++){
	if(isTissue(x,y,z)) S->numPoints++;
      }
    }
  }
  CALLOC(S->points,S->numPoints+1,sizeof(AnisoPoint)); /* ENUMERATION STARTS WITH 1 !!! */
  
  i = 1; /* ENUMERATION STARTS WITH 1 !!! */
  for (color=0;color<8;color++) {
    for (x=s.x0;x<=s.x1;x++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (z=s.z0;z<=s.z1;z++) {
	  the_color=(x%2)+2*(y%2)+4*(z%2);
	  if (the_color==color) {
	    if (isTissue(x,y,z)) {
	      
	      /*********************************/
	      /* The diff tensor of this point */
	      /* We shell treat dim<3 cases as dim=3 cases where the problem */
	      /*   simply does not depend on z, or z and y, coords.          */
	      
	      /* Fibre direction cosines of this point */
	      f[1] =  Geom[geom_ind(x,y,z,GEOM_FIBRE_1)];
	      f[2] =  Geom[geom_ind(x,y,z,GEOM_FIBRE_2)];
	      f[3] =  Geom[geom_ind(x,y,z,GEOM_FIBRE_3)];
	      /* set up the tensor */
	      for (j=1;j<=3;j++) {
		for (k=1;k<=3;k++) {
		  /* \mathbf{D} = 
		   * D_\perp \mathbf{I} + (D_\parallel - D\perp)
		   * \mathbf{A}\mathbf{A}^T
		   */
		  D[j][k] = ((Dpar-Dtrans)*f[j]*f[k]) + ((j==k)?Dtrans:0);
		} /* for k */
	      } /* for j */
	      
	      for (j=1;j<=3;j++) {
		if(
		   !isTissue( x+(j==1?1:0), y+(j==2?1:0), z+(j==3?1:0) ) ||
		   !isTissue( x-(j==1?1:0), y-(j==2?1:0), z-(j==3?1:0) )
		   ){
		  for (k=1;k<=3;k++) {dD[j][k]=0.0;}
		} else {
		  for (k=1;k<=3;k++) {
		    /* d^{11}_{+..} - d^{11}_{-..}
		     * = (D_\parallel - D_\perp)(a_j a_k)_{+..} -(a_j a_k)_{-..}
		     */
		    nx = x-(j==1?1:0);
		    ny = y-(j==2?1:0);
		    nz = z-(j==3?1:0);
		    px = x+(j==1?1:0);
		    py = y+(j==2?1:0);
		    pz = z+(j==3?1:0); /* spatial derivatives of the tensor */
		    dD[j][k] =
		      (Dpar - Dtrans) * 
		      (
		       (Geom[geom_ind(px,py,pz,(GEOM_FIBRE_1+j-1) )] *
			Geom[geom_ind(px,py,pz,(GEOM_FIBRE_1+k-1) )] )-
		       (Geom[geom_ind(nx,ny,nz,(GEOM_FIBRE_1+j-1) )] *
			Geom[geom_ind(nx,ny,nz,(GEOM_FIBRE_1+k-1) )] )
		       );
		  } /* for k */
		} /* else */
	      } /* for j */
	      
	      for (j=1;j<=3;j++) {
		c[j] = 0;
		for (k=1;k<=3;k++) {
		  c[j] += dD[k][j];
		}
	      }
	      
	      /******************************************/
	      /* Compute the connections to the neighbours. */
	      
	      /* For dim<3 cases, some connections will be calculated but not used. */
	      /* This can be optimized one day, to save memory.                 */
	      
	      CALLOC(P.K,maxnb,sizeof(real));
	      
	      /* Transform tensor coefficients into the sum weights. */
	      /* For dim<3, if (isTissue) will automatically ignore some of these lines. */
	      if (isTissue(x+1,y	 ,z  )) { P.K[poo]+=4*D[1][1]; P.K[ooo]-=4*D[1][1]; }
	      if (isTissue(x-1,y	 ,z  ))	{ P.K[moo]+=4*D[1][1]; P.K[ooo]-=4*D[1][1]; }
	      if (isTissue(x  ,y+1,z  ))	{ P.K[opo]+=4*D[2][2]; P.K[ooo]-=4*D[2][2]; }
	      if (isTissue(x  ,y-1,z  ))	{ P.K[omo]+=4*D[2][2]; P.K[ooo]-=4*D[2][2]; }
	      if (isTissue(x  ,y  ,z+1))	{ P.K[oop]+=4*D[3][3]; P.K[ooo]-=4*D[3][3]; }
	      if (isTissue(x  ,y	 ,z-1))	{ P.K[oom]+=4*D[3][3]; P.K[ooo]-=4*D[3][3]; }
	      if (isTissue(x+1,y  ,z  ) && isTissue(x-1,y  ,z  )) { P.K[poo]+=c[1]; P.K[moo]-=c[1]; }
	      if (isTissue(x  ,y+1,z  ) && isTissue(x  ,y-1,z  )) { P.K[opo]+=c[2]; P.K[omo]-=c[2]; }
	      if (isTissue(x  ,y  ,z+1) && isTissue(x  ,y  ,z-1)) { P.K[oop]+=c[3]; P.K[oom]-=c[3]; }
	      if (isTissue(x+1,y+1,z  ) && isTissue(x+1,y-1,z  )) { P.K[ppo]+=2*D[1][2]; P.K[pmo]-=2*D[1][2]; }
	      if (isTissue(x-1,y-1,z  ) && isTissue(x-1,y+1,z  )) { P.K[mmo]+=2*D[1][2]; P.K[mpo]-=2*D[1][2]; }
	      if (isTissue(x+1,y  ,z+1) && isTissue(x+1,y  ,z-1)) { P.K[pop]+=2*D[1][3]; P.K[pom]-=2*D[1][3]; }
	      if (isTissue(x-1,y  ,z-1) && isTissue(x-1,y  ,z+1)) { P.K[mom]+=2*D[1][3]; P.K[mop]-=2*D[1][3]; }
	      if (isTissue(x  ,y+1,z+1) && isTissue(x  ,y+1,z-1)) { P.K[opp]+=2*D[2][3]; P.K[opm]-=2*D[2][3]; }
	      if (isTissue(x  ,y-1,z-1) && isTissue(x  ,y-1,z+1)) { P.K[omm]+=2*D[2][3]; P.K[omp]-=2*D[2][3]; }
	      
	      /* Far corners maybe involved in coarsegrained connection matrices ? */
	      P.K[ppp]=P.K[ppm]=P.K[pmp]=P.K[pmm]=P.K[mpp]=P.K[mpm]=P.K[mmp]=P.K[mmm]=0.0;

	      /* Otherwise it will be division by zero full stop */
	      if (P.K[ooo]==0) ABORT("P.K[ooo]=0 at x,y,z=%d,%d,%d\n",x,y,z);

	      /* Normalization chosen to make u and f of the same scale */
	      for (nb=0;nb<maxnb;nb++) P.K[nb]/=(4*hx*hx);
	      
	      /* And the coordinates of this point */
	      P.x = x;
	      P.y = y;
	      P.z = z;
	      
	      /******************************/
	      /* Add this point to the list */
	      S->points[i] = P;
	      i++;
	    } /* if isTissue */
	  } /* if the_color */
	} /* for z */
      } /* for y */
    } /* for x */
  } /* for color */

  return 1;
} /* list_points */
/*=========================================*/

/* Prepare the hierarchy of grids, recursively */
static int make_level(int level,STR *S) {
  Space s =S->s;
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  coord_t x0, x1, y0, y1, z0, z1;
  coord_t xlen, ylen, zlen;
  coord_t xl, yl, zl, xh, yh, zh;
  coord_t xng, yng, zng;
  coord_t nx, ny, nz;
  coord_t x, y, z;
  scoord_t xpin, ypin, zpin; /* these may happen to be negative */
  scoord_t xcnl, ycnl, zcnl;

  coord_t levelsize;
  index_t indl, indh, the_point;
  sindex_t the_contributor; /* may happen to be negative */
  index_t the_low_neighbour, the_suspect, the_suspects_contributor;
  sindex_t the_neighbour; /* may happen to be negative */
  list_t npoints, i, j, l, i1;
  contr_t contributorNum, numContributorsMax;
  contr_t ic, jc, ic1, nc; 
  neighbour_t nb, nbl, nbh;
  color_t color, the_color;
  real total_weight, the_weight, the_suspects_weight, the_connection, the_contribution;

  if (level==0) {

    /* The grids cover the (local) device space */
    x0=S->s.x0; x1=S->s.x1; 
    y0=S->s.y0; y1=S->s.y1; 
    z0=S->s.z0; z1=S->s.z1; 
    /* "Netto" dimensions */
    nx=S->nx[0]=x1-x0+1;
    ny=S->ny[0]=y1-y0+1;
    nz=S->nz[0]=z1-z0+1;
    /* "Brutto" dimensions */
    xlen=S->xlen[0]=nx+2*ONE;
    ylen=S->ylen[0]=ny+2*TWO;
    zlen=S->zlen[0]=nz+2*TRI;
    /* "Netto" limits (upper inclusive) */
    S->xmin[level]=ONE;    S->xmax[level]=nx-1+ONE;
    S->ymin[level]=TWO;    S->ymax[level]=ny-1+TWO;
    S->zmin[level]=TRI;    S->zmax[level]=ny-1+TRI;
    S->levelsize[level]=levelsize=xlen*ylen*zlen;
    CALLOC(S->ugrid[0],levelsize,sizeof(real));
    CALLOC(S->rgrid[0],levelsize,sizeof(real));
    CALLOC(S->fgrid[0],levelsize,sizeof(real));

    /* List of points */
    S->nPoints[0]=S->numPoints;
    /* S->gridPoints[0]=S->points; - not good as points are in global coords and gridPoints in local netto coords */
    CALLOC(S->gridPoints[0],S->nPoints[0]+1,sizeof(AnisoPoint)); /* ENUMERATION STARTS WITH 1 !!! */
    CALLOC(S->pointIndex[0],S->nPoints[0]+1,sizeof(index_t)); /* ENUMERATION STARTS WITH 1 !!! */

    /* Level 0 arrays */
    S->contributors[0]=NULL;	/* level 0 nodes are 'atomic'    */
    S->R[0]=NULL;		/* i.e. do not have contributors */
    S->P[0]=NULL;		/* i.e. do not have contributors */
    CALLOC(S->d[0],maxnb,sizeof(sindex_t));
    S->d[0][ooo] = 0;
    S->d[0][poo] = gridIndex(0,+1, 0, 0) - gridIndex(0,0,0,0);
    S->d[0][moo] = gridIndex(0,-1, 0, 0) - gridIndex(0,0,0,0);
    S->d[0][opo] = gridIndex(0, 0,+1, 0) - gridIndex(0,0,0,0);
    S->d[0][omo] = gridIndex(0, 0,-1, 0) - gridIndex(0,0,0,0);
    S->d[0][oop] = gridIndex(0, 0, 0,+1) - gridIndex(0,0,0,0);
    S->d[0][oom] = gridIndex(0, 0, 0,-1) - gridIndex(0,0,0,0);
    S->d[0][ppo] = gridIndex(0,+1,+1, 0) - gridIndex(0,0,0,0);
    S->d[0][pmo] = gridIndex(0,+1,-1, 0) - gridIndex(0,0,0,0);
    S->d[0][mpo] = gridIndex(0,-1,+1, 0) - gridIndex(0,0,0,0);
    S->d[0][mmo] = gridIndex(0,-1,-1, 0) - gridIndex(0,0,0,0);
    S->d[0][pop] = gridIndex(0,+1, 0,+1) - gridIndex(0,0,0,0);
    S->d[0][pom] = gridIndex(0,+1, 0,-1) - gridIndex(0,0,0,0);
    S->d[0][mop] = gridIndex(0,-1, 0,+1) - gridIndex(0,0,0,0);
    S->d[0][mom] = gridIndex(0,-1, 0,-1) - gridIndex(0,0,0,0);
    S->d[0][opp] = gridIndex(0, 0,+1,+1) - gridIndex(0,0,0,0);
    S->d[0][opm] = gridIndex(0, 0,+1,-1) - gridIndex(0,0,0,0);
    S->d[0][omp] = gridIndex(0, 0,-1,+1) - gridIndex(0,0,0,0);
    S->d[0][omm] = gridIndex(0, 0,-1,-1) - gridIndex(0,0,0,0);
    /* Far corners maybe involved in coarsegrained connection matrices ? */
    S->d[0][ppp] = gridIndex(0,+1,+1,+1) - gridIndex(0,0,0,0);
    S->d[0][ppm] = gridIndex(0,+1,+1,-1) - gridIndex(0,0,0,0);
    S->d[0][pmp] = gridIndex(0,+1,-1,+1) - gridIndex(0,0,0,0);
    S->d[0][pmm] = gridIndex(0,+1,-1,-1) - gridIndex(0,0,0,0);
    S->d[0][mpp] = gridIndex(0,-1,+1,+1) - gridIndex(0,0,0,0);
    S->d[0][mpm] = gridIndex(0,-1,+1,-1) - gridIndex(0,0,0,0);
    S->d[0][mmp] = gridIndex(0,-1,-1,+1) - gridIndex(0,0,0,0);
    S->d[0][mmm] = gridIndex(0,-1,-1,-1) - gridIndex(0,0,0,0);
    CALLOC(S->D[0],maxcn,sizeof(sindex_t));
    /* For the chosen vertex-centered grid, */
    /* neighbours' and contributors' vectors coincide */
    for (nb=0;nb<maxnb;nb++) S->D[0][nb]=S->d[0][nb];
    /* S->D[0][PPP] = gridIndex(0,+1,+1,+1) - gridIndex(0,0,0,0); */
    /* S->D[0][PPM] = gridIndex(0,+1,+1,-1) - gridIndex(0,0,0,0); */
    /* S->D[0][PMP] = gridIndex(0,+1,-1,+1) - gridIndex(0,0,0,0); */
    /* S->D[0][PMM] = gridIndex(0,+1,-1,-1) - gridIndex(0,0,0,0); */
    /* S->D[0][MPP] = gridIndex(0,-1,+1,+1) - gridIndex(0,0,0,0); */
    /* S->D[0][MPM] = gridIndex(0,-1,+1,-1) - gridIndex(0,0,0,0); */
    /* S->D[0][MMP] = gridIndex(0,-1,-1,+1) - gridIndex(0,0,0,0); */
    /* S->D[0][MMM] = gridIndex(0,-1,-1,-1) - gridIndex(0,0,0,0); */

    switch (dim) {
    case 1: S->numNeighbours=3;  nc=S->numContributors=3;  break;
    case 2: S->numNeighbours=9;  nc=S->numContributors=9;  break;
    case 3: S->numNeighbours=27/* 19 */; nc=S->numContributors=27; break; 
      /* Far corners maybe involved in coarsegrained connection matrices ? */
    default: ABORT("wrong dimensionality %d\n",dim);
    }

    /* We want ggrid containt list numbers rather than 0/1. */
    /* Also convert global to local netto coordinates. */
    S->numBClinks=0;
    for (i=1;i<=S->nPoints[0];i++) { /* ENUMERATION STARTS WITH 1 !!! */
      AnisoPoint Pg=S->points[i];	/* this has global coordinates */
      AnisoPoint *P=&(S->gridPoints[0][i]); /* this has local netto coordinates */
      /* P->K=Pg.K;  - this copies the address of array K - not good as we need to modify them */
      x=P->x = Pg.x-x0;			/* (x0,y0,z0) are global coordinates */
      y=P->y = Pg.y-y0;			/*        of the local point (0,0,0) */
      z=P->z = Pg.z-z0;
      indh = gridIndex(0,x,y,z);	/* 1D index of the point */
      S->pointIndex[0][i]=indh;		/* remember it in the point list */
      S->ggrid[0][indh]=i;		/* ggrid will contain the list number of this point */

      /* The treatment of K below takes care of non-homogeneous          */
      /*     boundary conditions (links to values beyond device space)   */
      CALLOC(P->K,S->numNeighbours,sizeof(real));	/* Start the list of connections from scratch */
      for (nbh=0;nbh<S->numNeighbours;nbh++) {	/* Look at all neighbours */
	the_neighbour=indh+S->d[0][nbh];/* Neighbour's 1D index,  */
	xng=gridX(0,the_neighbour)+x0;	/*   and its              */
	yng=gridY(0,the_neighbour)+y0;	/*     global             */
	zng=gridZ(0,the_neighbour)+z0;	/*       coordinates      */

	if (s.global_x0<=xng && xng <= s.global_x1 && 
	    s.global_y0<=yng && yng <= s.global_y1 && 
	    s.global_z0<=zng && zng <= s.global_z1) {
	  P->K[nbh]=Pg.K[nbh];		/* the neighbour is within space, valid connection */
	} else if (Pg.K[nbh]) {		/* its outside space, neighbour's value is nonzero */
	  P->K[nbh]=0;			/* no connection within the grid */
	  S->numBClinks++;		/* but count links to outside values */
	} else {			/* outside and no connection outside */
	  P->K[nbh]=0;			/* no connection within the grid either */
	}
      }
    }

    /* Now make the links to the outide values. */
    /* TODO: estimate max numBClinks so could do all work in one go? */
    if (S->numBClinks>0) {
      if (S->pin) {
	MESSAGE("/* Pinning was specified while non-homogenoues boundary conditions available.\n"
		"   Pinning will be honoured but this may lead to trouble. */\n");
      }

      CALLOC(S->BClinks,sizeof(bclink),S->numBClinks);
      l=0;
      for (i=1;i<=S->nPoints[0];i++) { /* ENUMERATION STARTS WITH 1 !!! */
	AnisoPoint Pg=S->points[i];	  /* the global list */
	AnisoPoint P=S->gridPoints[0][i]; /* this device's list */
	x=P.x;
	y=P.y;
	z=P.z;
	indh = gridIndex(0,x,y,z);		/* 1D index of the point */
	for (nbh=0;nbh<S->numNeighbours;nbh++) {		/* Look at all neighbours */
	  the_neighbour=indh+S->d[0][nbh];	/* Neighbour's 1D index,  */
	  xng=gridX(0,the_neighbour)+x0;	/*   and its              */
	  yng=gridY(0,the_neighbour)+y0;	/*     global             */
	  zng=gridZ(0,the_neighbour)+z0;	/*       coordinates      */
	  
	  if (!(
		s.global_x0<=xng && xng<=s.global_x1 && 
		s.global_y0<=yng && yng<=s.global_y1 && 
		s.global_z0<=zng && zng<=s.global_z1)) {
	    if (Pg.K[nbh]) { /* No point remembering the link if the connection is zero */
	      S->BClinks[l].to=indh;
	      S->BClinks[l].from=&(New[ind(xng,yng,zng,v1)]); /* link to SOLUTION outside the space */
	      S->BClinks[l].K=Pg.K[nbh];
	      l++;
	    } /* if Pg.K[nbh] */
	  } /* if in space */
	} /* for nbh */
      } /* for i */
    } /* if numBClinks */ 

    /* Impose the pinning condition by discarding one of the equations             */
    /* (they are linearly dependent anyway, if solvability condition is satisfied) */
    if (S->pin) {
      xpin=S->xpin-x0;			/* (S->xpin,...) are global coordinates,*/
      ypin=S->ypin-y0;			/* xpin, ... are local coordinates      */
      zpin=S->zpin-z0;
      if (0<=xpin && xpin<nx && 0<=ypin && ypin<ny && 0<=zpin && zpin<nz) {
	real *K;  /* the pointer to the array of connections==coefficients of the linear system */
	indh = gridIndex(0,x,y,z);	/* 1D index of the pin point */
	i=pointNumber(0,indh);		/* index of that point in the list */
	if (i) {			/* make sure it is tissue */
	  DEVICE_CONST(real,hx);
	  K=S->gridPoints[0][i].K;	/* the array of coefficients */
	  K[0]=PINCOEFF/(hx*hx);	/* we fix the value at the point; make coeff commensurate with others */
	  for (nbh=1;nbh<S->numNeighbours;nbh++) K[nbh]=0; /* and ignore the neighbours */
	} else {
	  ABORT("pinning point (global %d,%d,%d, local %d,%d,%d) is not tissue\n",
		S->xpin,S->ypin,S->zpin,xpin,ypin,zpin);
	}
	/* And the fixed value upin has to be applied at run time, not build time */
      } /* otherwise the pinning point is in another subdomain so we don't care here */
    }

    if (!S->numBClinks && !S->pin) {
      MESSAGE("/* Warning: Neumann problem without pinning condition detected. Expect troubles!!! */\n");
    }

    /* end of level==0 case */
  } else {
    /* Sizes of the subgrid. Take care of the halos: they are not coarsened! */
    /* "Netto" dimensions */
    nx=(S->nx[level])=coarsex(S->nx[level-1]+1);
    ny=(S->ny[level])=coarsey(S->ny[level-1]+1);
    nz=(S->nz[level])=coarsez(S->nz[level-1]+1);
    /* "Brutto" dimensions */
    S->xlen[level]=nx+2*ONE;
    S->ylen[level]=ny+2*TWO;
    S->zlen[level]=ny+2*TRI;
    /* "Netto" limits (upper inclusive) */
    S->xmin[level]=ONE; S->xmax[level]=nx-1+ONE;
    S->ymin[level]=TWO; S->ymax[level]=ny-1+TWO;
    S->zmin[level]=TRI; S->zmax[level]=ny-1+TRI;
    S->levelsize[level]=levelsize=S->xlen[level]*S->ylen[level]*S->zlen[level];
    CALLOC(S->ugrid[level],levelsize,sizeof(real));
    CALLOC(S->rgrid[level],levelsize,sizeof(real));
    CALLOC(S->fgrid[level],levelsize,sizeof(real));
    /* Index displacements for neighbours */
    CALLOC(S->d[level],maxnb,sizeof(sindex_t));
    S->d[level][ooo] = 0;
    S->d[level][poo] = gridIndex(level,+1, 0, 0) - gridIndex(level,0,0,0);
    S->d[level][moo] = gridIndex(level,-1, 0, 0) - gridIndex(level,0,0,0);
    S->d[level][opo] = gridIndex(level, 0,+1, 0) - gridIndex(level,0,0,0);
    S->d[level][omo] = gridIndex(level, 0,-1, 0) - gridIndex(level,0,0,0);
    S->d[level][oop] = gridIndex(level, 0, 0,+1) - gridIndex(level,0,0,0);
    S->d[level][oom] = gridIndex(level, 0, 0,-1) - gridIndex(level,0,0,0);
    S->d[level][ppo] = gridIndex(level,+1,+1, 0) - gridIndex(level,0,0,0);
    S->d[level][pmo] = gridIndex(level,+1,-1, 0) - gridIndex(level,0,0,0);
    S->d[level][mpo] = gridIndex(level,-1,+1, 0) - gridIndex(level,0,0,0);
    S->d[level][mmo] = gridIndex(level,-1,-1, 0) - gridIndex(level,0,0,0);
    S->d[level][pop] = gridIndex(level,+1, 0,+1) - gridIndex(level,0,0,0);
    S->d[level][pom] = gridIndex(level,+1, 0,-1) - gridIndex(level,0,0,0);
    S->d[level][mop] = gridIndex(level,-1, 0,+1) - gridIndex(level,0,0,0);
    S->d[level][mom] = gridIndex(level,-1, 0,-1) - gridIndex(level,0,0,0);
    S->d[level][opp] = gridIndex(level, 0,+1,+1) - gridIndex(level,0,0,0);
    S->d[level][opm] = gridIndex(level, 0,+1,-1) - gridIndex(level,0,0,0);
    S->d[level][omp] = gridIndex(level, 0,-1,+1) - gridIndex(level,0,0,0);
    S->d[level][omm] = gridIndex(level, 0,-1,-1) - gridIndex(level,0,0,0);
    S->d[level][ppp] = gridIndex(level,+1,+1,+1) - gridIndex(level,0,0,0);
    /* Far corners maybe involved in coarsegrained connection matrices ? */
    S->d[level][ppm] = gridIndex(level,+1,+1,-1) - gridIndex(level,0,0,0);
    S->d[level][pmp] = gridIndex(level,+1,-1,+1) - gridIndex(level,0,0,0);
    S->d[level][pmm] = gridIndex(level,+1,-1,-1) - gridIndex(level,0,0,0);
    S->d[level][mpp] = gridIndex(level,-1,+1,+1) - gridIndex(level,0,0,0);
    S->d[level][mpm] = gridIndex(level,-1,+1,-1) - gridIndex(level,0,0,0);
    S->d[level][mmp] = gridIndex(level,-1,-1,+1) - gridIndex(level,0,0,0);
    S->d[level][mmm] = gridIndex(level,-1,-1,-1) - gridIndex(level,0,0,0);
    CALLOC(S->D[level],maxcn,sizeof(sindex_t));
    /* For the chosen vertex-centered grid, */
    /* neighbours' and contributors' vectors coincide */
    for (nb=0;nb<maxnb;nb++) S->D[level][nb]=S->d[level][nb];

    nc=S->numContributors;

    /* Count contributors to coarse grid points */
    for (xh=0;xh<nx;xh++) {
      xl=finex(xh);
      for (yh=0;yh<ny;yh++) {
	yl=finey(yh);
	for (zh=0;zh<nz;zh++) {
	  zl=finez(zh);
	  indl=gridIndex(level-1,xl,yl,zl);
	  indh=gridIndex(level,xh,yh,zh);
	  for (ic=0;ic<nc;ic++) {			/* loop through all contributors */
	    xcnl=xl+cdx[ic];				/* check that the 3D coords of the contributor */
	    ycnl=yl+cdy[ic];
	    zcnl=zl+cdz[ic];
	    if (
		-ONE<=xcnl && xcnl<=S->nx[level-1]+ONE &&	/* are within the extended box */
		-TWO<=ycnl && ycnl<=S->ny[level-1]+TWO &&
		-TRI<=zcnl && zcnl<=S->nz[level-1]+TRI 
		) {
	      the_contributor=indl+S->D[level-1][ic];	/* contributor's would-be 1D index */
        	      /* debug printout */
        	      if (the_contributor<0 || the_contributor >= S->levelsize[level-1]) {
        		ABORT("\nwrong contributor index %ld ic=%ld of point index %ld of level %ld\n",
        		       (long)the_contributor,(long)ic,(long)indl,(long)level-1);
        	      }
	      /* Make sure higher level point is included */
	      /*  if at least one contributor is included */	
	      S->ggrid[level][indh]+=S->ggrid[level-1][the_contributor];
	    } /* if ic is within extended box */
	  } /* for ncl */
	} /* for zh */
      } /* for yh */
    } /* for xh */

    /* Count the non-empty vertices */
    npoints=0;
    for (xh=0;xh<nx;xh++) {
      for (yh=0;yh<=ny;yh++) {
	for (zh=0;zh<nz;zh++) {
	  indh=gridIndex(level,xh,yh,zh);
	  if (S->ggrid[level][indh]>0) npoints++;
	}
      }
    }

    /* Allocate the list of points */
    CALLOC(S->gridPoints[level],npoints+1,sizeof(AnisoPoint)); /* ENUMERATION STARTS WITH 1 !!! */
    CALLOC(S->pointIndex[level],npoints+1,sizeof(index_t)); /* ENUMERATION STARTS WITH 1 !!! */
    S->nPoints[level]=npoints;

    /* Make the enumerated list of points and */
    CALLOC(S->contributors[level],npoints+1,sizeof(contr_t *));
    CALLOC(S->P[level],npoints+1,sizeof(real *));
    CALLOC(S->R[level],npoints+1,sizeof(real *));
    i=1;						/* ENUMERATION STARTS WITH 1 !!! */

    for (color=0;color<8;color++) {
      for (xh=0;xh<nx;xh++) {
	for (yh=0;yh<=ny;yh++) {
	  for (zh=0;zh<nz;zh++) {
	    the_color=(xh%2)+2*(yh%2)+4*(zh%2);
	    if (the_color==color) {
	      indh=gridIndex(level,xh,yh,zh);
	      if (S->ggrid[level][indh]>0) {
		/* ... mutually link them with the grid */
		S->pointIndex[level][i]=indh;
		S->ggrid[level][indh]=i;
		/* ... their coordinates */
		S->gridPoints[level][i].x=xh;
		S->gridPoints[level][i].y=yh;
		S->gridPoints[level][i].z=zh;
		/* ... initialize their connections */
		CALLOC(S->gridPoints[level][i].K,S->numNeighbours,sizeof(real));
		/* ... and initialize their contributors lists */
		CALLOC(S->contributors[level][i],nc,sizeof(size_t));
		CALLOC(S->P[level][i],nc,sizeof(real));
		CALLOC(S->R[level][i],nc,sizeof(real));
		i++;
	      } /* if ggrid */
	    } /* if the_color */
	  } /* for zh */
	} /* for yh */
      } /* for xh */
    } /* for color */

    /* Assign the weights of contributors */
    for (i=1;i<=S->nPoints[level-1];i++) /* zero out the lower Pgrid of P row totals */
	S->Pgrid[level-1][ S->pointIndex[level-1][i] ]=0;

    /* levelsize=S->nx[level-1]*S->ny[level-1]*S->nz[level-1]; */
    levelsize=S->levelsize[level-1];
    for(indl=0;indl<levelsize;indl++) S->ugrid[level-1][indl]=0.0;

    for (xh=0;xh<nx;xh++) {
      xl=finex(xh);
      for (yh=0;yh<ny;yh++) {
	yl=finey(yh);
	for (zh=0;zh<nz;zh++) {
	  zl=finez(zh);
	  indl=gridIndex(level-1,xl,yl,zl);
	  indh=gridIndex(level,xh,yh,zh);
	  i=pointNumber(level,indh);			/* List number of high-level node. */
	  if (i) {
	    total_weight=0;				/* need sum of weights of existing contributors */
	    for (ic=0;ic<nc;ic++) {			/* loop through all contributors */
	      xcnl=xl+cdx[ic];				/* check that the 3D coords of the contributor */
	      ycnl=yl+cdy[ic];
	      zcnl=zl+cdz[ic];
	      if (
		  -ONE<=xcnl && xcnl<=S->nx[level-1]+ONE &&	/* are within the extended box */
		  -TWO<=ycnl && ycnl<=S->ny[level-1]+TWO &&
		  -TRI<=zcnl && zcnl<=S->nz[level-1]+TRI 
		  ) {
		the_contributor=indl+S->D[level-1][ic];		/* contributor's would-be 1D index */
		S->contributors[level][i][ic]=the_contributor;	/* register it for fast access */
		j=pointNumber(level-1,the_contributor);	/* number of that point in the list of points */
		the_weight=(j)?(wgt[ic]):(0.0);		/* a real point has some weight */
		S->R[level][i][ic]=the_weight;		/* register the P-weight */
		S->P[level][i][ic]=the_weight;		/* register the R-weight */
		total_weight+=the_weight;			/* update "row sum" for R matrix */
		S->Pgrid[level-1][the_contributor]+=the_weight; /* update "row sum" for P matrix */
	      } /* if ic is within extended box */
	    } /* for ic */
	    ASSERT(total_weight>0);
	    for(ic=0;ic<nc;ic++)			/* normalize the R matrix to unit row sums */
	      S->R[level][i][ic]/=total_weight;
	  } /* if i */
	} /* for zh */
      } /* for yh */
    } /* for xh */

    /* Normalize the P matrix to unit row sums */
    for (xh=0;xh<nx;xh++) {
      xl=finex(xh);
      for (yh=0;yh<ny;yh++) {
	yl=finey(yh);
	for (zh=0;zh<nz;zh++) {
	  zl=finez(zh);
	  indl=gridIndex(level-1,xl,yl,zl);
	  indh=gridIndex(level,xh,yh,zh);
	  i=pointNumber(level,indh);			/* List number of high-level node. */
	  if (i) {
	    for (ic=0;ic<nc;ic++) {			/* loop through all possible contributors */
	      the_contributor=indl+S->D[level-1][ic];		/* its would-be 1D index */
	      j=pointNumber(level-1,the_contributor);	/* number of that point in the list of points */
	      if (j) {
		total_weight=S->Pgrid[level-1][the_contributor];
		if (total_weight)			/* some corners are never reached, esp in 3D */
		  S->P[level][i][ic]/=total_weight;
	      } /* if j */
	    } /* for ic */
	  } /* if i */
	} /* for zh */
      } /* for yh */
    } /* for xh */

    /* Make the coarse system Q=RAP */
    for (i=1;i<=S->nPoints[level];i++) { /* ENUMERATION STARTS WITH 1 !!! */
                                                        /* Loop through higher level points */
      real Rvalue, Pvalue, Avalue;
      index_t the_main_contributor, the_side_contributor, the_main_contributors_neighbour;
      the_point=S->pointIndex[level][i];		/* point's 1D index */
      for (nbh=0;nbh<S->numNeighbours;nbh++) {		/* Loop over coarse neighbours */
	the_neighbour=the_point+S->d[level][nbh];
	if (0!=(j=pointNumber(level,the_neighbour))) {	/* if it is a real point */	
	  for (ic=0;ic<nc;ic++) {			/* loop over main point contributors */
	    if (0!=(Rvalue=S->R[level][i][ic])) {	/* main contributor's weight */
	      the_main_contributor=S->contributors[level][i][ic];
	      i1=pointNumber(level-1,the_main_contributor);
	      if (!i1) {
		ABORT("level=%d i=%d ic=%d Rvalue=%g "
		      "level-1=%d the_main_contributor=%d i1=%d\n",
		      (int)level,(int)i,(int)ic,(float)Rvalue,
		      (int)(level-1),(int)the_main_contributor,(int)i1
		      );
	      }
	      for (jc=0;jc<nc;jc++) {			/* loop over side point (=neighbour's) contributors */
		if (0!=(Pvalue=S->P[level][j][jc])) {	/* side contributor's weight */
		  the_side_contributor=S->contributors[level][j][jc];
		  for (nbl=0;nbl<S->numNeighbours;nbl++) {/* Are the two contributors neighbours? */
		    the_main_contributors_neighbour=the_main_contributor+S->d[level-1][nbl];
		    if (the_main_contributors_neighbour==the_side_contributor) {	/* yes they are */
		      Avalue=S->gridPoints[level-1][i1].K[nbl];	/* connection = element of A matrix */
		      S->gridPoints[level][i].K[nbh] +=			/* element of RAP matrix */
			Rvalue*Avalue*Pvalue;
		    } /* if the_mcn==th_sc */
		  } /* for nbl */
		} /* if Pvalue */
	      } /* for ic1 */
	    } /* if Rvalue */
	  } /* for ic */
	} /* if j */
      } /* for nbh */
    } /* for i */
  } /* level>0 case */

#if MPI
  define_halo_level(level,S);
  /* We need to exchange halos for ggrid=rgrid to be aware */
  /* of tissue presence in neighbouring subdomains */
  halo_level(level,S->ggrid[level],S);
#endif

  if (Verbose) MESSAGE("/* level=%d grid=%dx%dx%d */\n",level,nx,ny,nz);

  /* Diagnostics to report to user: how many bad cells     */
  /* NB: bad cells are those with "ratio" smaller than 1/2 */
  if (Verbose) {
    int badcells;
    list_t baddest;
    real *K, Kooo, Knbh;
    real cellsum, maxcellsum, abscellsum, ratio, worstratio, bestratio;
    maxcellsum=0;
    baddest=0;
    worstratio=MAXREAL;
    bestratio=0;
    badcells=0;
    for (i=1;i<=S->nPoints[level];i++) {
      K=S->gridPoints[level][i].K;
      cellsum=K[0];
      abscellsum=Kooo=fabs(K[0]);
      for(nbh=1;nbh<S->numNeighbours;nbh++) {
	Knbh=K[nbh];
	cellsum+=Knbh;
	abscellsum+=fabs(Knbh);
      }
      ratio=Kooo/abscellsum;
      if (ratio<worstratio) {
	worstratio=ratio;
	baddest=i;
      }
      if (ratio>bestratio) {
	bestratio=ratio;
      }
      if (ratio<0.5-2*macheps) 
	badcells+=1;
      cellsum=fabs(cellsum);
      if (cellsum>maxcellsum) maxcellsum=cellsum;
    }
    if (badcells && Verbose) {
      MESSAGE("/* WARNING: at level %d there are %d/%d 'bad cells' that may cause divergence of iterations */\n",
	      level,badcells,S->nPoints[level]);
      K=S->gridPoints[level][baddest].K;
      x=S->gridPoints[level][baddest].x;
      y=S->gridPoints[level][baddest].y;
      z=S->gridPoints[level][baddest].z;
      MESSAGE("/* Level %d ratios are from %.15g to %.15g. The worst cell is #%d (%d,%d,%d) and its 2D stencil is: */\n"
	      "/* %g\t%g\t%g */\n/* %g\t%g\t%g */\n/* %g\t%g\t%g */\n",
	      level,worstratio,bestratio,baddest,x,y,z,
	      K[mpo],K[opo],K[ppo],
	      K[moo],K[ooo],K[poo],
	      K[mmo],K[omo],K[pmo]);
    }
  }

  /* Here goes the recursion */
  if (level<S->upper_level) if (NOT(make_level(level+1,S))) return 0;

  return 1;
}
/*=========================================*/


/* This function does a GS iteration, lexicographic order, */
/*	for the original data in New array. */
/* TODO: consider extending/making version of this for Jacobi iterations */
/*	 and/or for red-black/multicolour order? */
static real iterate_base(int it, STR *S) {      
  DEVICE_CONST(Space,s)
  DEVICE_CONST(int,v0)
  DEVICE_CONST(int,v1)
  DEVICE_CONST(list_t, numPoints)
  DEVICE_ARRAY(AnisoPoint, points)
  AnisoPoint P;
  coord_t x,y,z;
  list_t i;
  real *u;
  real unew,f,sum,absres;
  double resnorm, globalresnorm; /* fixed type for MPI_REDUCE */
#if MPI
  DEVICE_CONST(MPI_Comm, comm)
  DEVICE_CONST(int, root)
#endif

#if MPI
  /* Synchronise the subdomains at this point */
  if (it) haloSwap(); 
#endif	
  
  /* Will measure the progress of the iteration */
  resnorm=0.0;
  
  /* One step, described separately for each dimensions case */
  switch (dim) {
  case 3:
    for(i=1;i<=numPoints;i++){ /* ENUMERATION STARTS WITH 1 !!! */
      P=points[i];
      x = P.x;
      y = P.y;
      z = P.z;
      /* Centre point */
      u=New+ind(x,y,z,v1); 
      /* The righ-hand side */
      f=New[ind(x,y,z,v0)]; 
      /* Sum by proper neighbours */
      sum = /* not including the center point */
	+ u[DX]*P.K[poo] + u[-DX]*P.K[moo] 
	+ u[DY]*P.K[opo] + u[-DY]*P.K[omo]
	+ u[DZ]*P.K[oop] + u[-DZ]*P.K[oom]
	+ u[DX+DY]*P.K[ppo] + u[DX-DY]*P.K[pmo] + u[-DX+DY]*P.K[mpo] + u[-DX-DY]*P.K[mmo]
	+ u[DX+DZ]*P.K[pop] + u[DX-DZ]*P.K[pom] + u[-DX+DZ]*P.K[mop] + u[-DX-DZ]*P.K[mom]
	+ u[DY+DZ]*P.K[opp] + u[DY-DZ]*P.K[opm] + u[-DY+DZ]*P.K[omp] + u[-DY-DZ]*P.K[omm];
      unew=(f-sum)/P.K[ooo]; /* resolve wrt the center point */
      absres=fabs((unew-u[0])*P.K[ooo]);
      if (absres>resnorm) resnorm=absres;
      u[0]=unew;
    } /* for i */
    break;
    /* shortened version for lesser dimensions */
  case 2:
    for(i=1;i<=numPoints;i++){ /* ENUMERATION STARTS WITH 1 !!! */
      P=points[i];
      x = P.x;
      y = P.y;
      z = P.z;
      u=New+ind(x,y,z,v1); 
      f=New[ind(x,y,z,v0)]; 
      sum =
	+ u[DX]*P.K[poo] + u[-DX]*P.K[moo] 
	+ u[DY]*P.K[opo] + u[-DY]*P.K[omo]
	+ u[DX+DY]*P.K[ppo] + u[DX-DY]*P.K[pmo] + u[-DX+DY]*P.K[mpo] + u[-DX-DY]*P.K[mmo];
      unew=(f-sum)/P.K[ooo];
      absres=fabs((unew-u[0])*P.K[ooo]);
      if (absres>resnorm) resnorm=absres;
      u[0]=unew;
    }
    break;
  case 1:
    for(i=1;i<=numPoints;i++){ /* ENUMERATION STARTS WITH 1 !!! */
      P=points[i];
      x = P.x;
      y = P.y;
      z = P.z;
      u=New+ind(x,y,z,v1); 
      f=New[ind(x,y,z,v0)]; 
      sum =
	+ u[DX]*P.K[poo] + u[-DX]*P.K[moo];
      unew=(f-sum)/P.K[ooo];
      absres=fabs((unew-u[0])*P.K[ooo]);
      if (absres>resnorm) resnorm=absres;
      u[0]=unew;
    }
    break;
  default:
    ABORT("impossible dimension %d\n",dim);
  } /* switch dim */
  
  /* Test if iterations have converged */
#if MPI
  MPIDO(MPI_Reduce(&resnorm,&globalresnorm,1,MPI_DOUBLE,MPI_MAX, root, comm),
	"Couldn't carry out reduction operation.");
  MPIDO(MPI_Bcast(&globalresnorm, 1, MPI_DOUBLE, root, comm),
	"Couldn't broadcast reduction result.");
#else
  globalresnorm=resnorm;
#endif
  S->resnorm=globalresnorm;
  DEBUG(" iter=%ld resnorm=%lg\n", (long)it, (double)globalresnorm);
  return(globalresnorm); 
} /* iterate_base */
/*=========================================*/

/* Do a recursive V-cycle of the residual systems */
static int do_level(int level,STR *S) {
  DEVICE_CONST(Space, s)
  real tolerance=S->tol[level];
  int preiter, postiter, it;

  preiter=S->preIter[level];
  postiter=S->postIter[level];

  /* Pre-iterations */
  if (preiter>0) {
    for (it=0;it<preiter;it++)
      if (tolerance>iterate_level(level,S->damping,S)) break;
  } else { /* need residuals even of no presmoothing */
    iterate_level(level,0,S);
  }

  if (S->rsn[level]>=tolerance) {
    /* Solve system at higher level (this is the recursion) */
    if NOT(higher_level(level,S)) return 0;

    /* Post-iterations: */
    /* at least once after the higher-level */
    for (it=0;it<postiter;it++) {
      if (tolerance>iterate_level(level,S->damping,S)) break;
    }
  }
  DEBUG("result: level=%d iter=%ld rsn=%lg\n", level, (long)(S->iter[level]), (double)S->rsn[level]);
  return 1;
}
/*=========================================*/

/* Transfer data from New to the base of the multigrid */
/* TODO: consider, again, identifying ugrid[0] and fgrid[0] with subsets of New */
/* to avoid this unnecessary copying                                            */
static void to_level0(STR *S) {
  DEVICE_CONST(Space, s)
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  coord_t nx=S->nx[0];
  coord_t ny=S->ny[0];
  coord_t nz=S->nz[0];
  coord_t x,y,z;
  scoord_t xpin,ypin,zpin;
  index_t indh;
  list_t l;
  int x0=s.x0; int x1=s.x1;
  int y0=s.y0; int y1=s.y1;
  int z0=s.z0; int z1=s.z1;

  int level=0;

  for (x=0;x<nx;x++) {
    for (y=0;y<ny;y++) {
      for (z=0;z<nz;z++) {
	indh=gridIndex(0,x,y,z);
	S->ugrid[0][indh]=New[ind(x0+x,y0+y,z0+z,v1)];
	S->fgrid[0][indh]=New[ind(x0+x,y0+y,z0+z,v0)];
	S->rgrid[0][indh]=0; /* ??? let it contain stuff from prev use ? */
      }
    }
  }
  /* Now take the outside values into account */
  if (S->numBClinks) {
    for (l=0;l<S->numBClinks;l++) {
      bclink L=S->BClinks[l];
      S->fgrid[0][L.to]-=L.K*(*L.from);
    }
  } 

  /* and/or apply the pinning condition */
  if (S->pin) {
    xpin=S->xpin-x0;			/* (S->xpin,...) are global coordinates,*/
    ypin=S->ypin-y0;			/* xpin, ... are local coordinates      */
    zpin=S->zpin-z0;
    if (0<=xpin && xpin<nx && 0<=ypin && ypin<ny && 0<=zpin && zpin<nz) {
      DEVICE_CONST(real,hx);
      indh = gridIndex(0,x,y,z);	/* 1D index of the pin point */
      S->fgrid[0][indh]=S->upin*PINCOEFF/(hx*hx); /* scaling commensurate with other equations */
    } /* otherwise the pinning point is in another domain so we don't care here */
  } 
#if MPI
  halo_level(0,S->ugrid[0],S); /* "just in case" - do we really need it here? */
  halo_level(0,S->fgrid[0],S); /* "just in case" - do we really need it here? */
  halo_level(0,S->rgrid[0],S); /* "just in case" - do we really need it here? */
#endif
}
/*=========================================*/

/* Transfer data from the base of the multigrid to New */
/* TODO: consider, again, identifying ugrid[0] and fgrid[0] with subsets of New */
/* to avoid this unnecessary copying                                            */
static void from_level0(STR *S) {
  DEVICE_CONST(Space, s)
  DEVICE_CONST(int, v0)
  DEVICE_CONST(int, v1)
  coord_t nx=S->nx[0];
  coord_t ny=S->ny[0];
  coord_t nz=S->nz[0];
  coord_t x,y,z;
  index_t indh;
  int x0=s.x0; int x1=s.x1;
  int y0=s.y0; int y1=s.y1;
  int z0=s.z0; int z1=s.z1;
  for (x=0;x<nx;x++) {
    for (y=0;y<ny;y++) {
      for (z=0;z<nz;z++) {
	indh=gridIndex(0,x,y,z);
	New[ind(x0+x,y0+y,z0+z,v1)]=S->ugrid[0][indh];
	New[ind(x0+x,y0+y,z0+z,v0)]=S->fgrid[0][indh];
      }
    }
  }
}
/*=========================================*/

/* Do the higher level =  up, there, down */
static int higher_level(int level,STR *S) {
  if (level<S->upper_level) {
    if NOT(up_level(level,S)) return 0;		/* Form the coarse system. */
    if NOT(do_level(level+1,S)) return 0;	/* Solve the coarse system (this is recursion!!). */
    if NOT(down_level(level,S)) return 0;	/* Prolong and add the resulting correction */
  }
  return 1;
}
/*=========================================*/

/* Implement the Restriction Operator: pass residual to coarser system */
static int up_level(int level,STR *S) {
  list_t number_of_coarse_points, i;
  contr_t nc, ic;
  index_t *the_contributors, the_contributor, indh;
  real *the_weights;
  real the_weight;
  real usum, fsum, rsum;
  /* By-hand optimization */
  real *ul=S->ugrid[level];
  real *fl=S->fgrid[level];
  real *rl=S->rgrid[level];
  real *uh=S->ugrid[level+1];
  real *fh=S->fgrid[level+1];
  real *rh=S->rgrid[level+1];
  index_t **ch=S->contributors[level+1];
  real **Rh=S->R[level+1];

  if (S->transferchoice==correct)
    iterate_level(level,0,S);		/* Need residuals to move up. */

  PROFILE(UPLEVEL);
  number_of_coarse_points=S->nPoints[level+1];
  for (i=1;i<=number_of_coarse_points;i++) { /* ENUMERATION STARTS WITH 1 !!! */
    nc = S->numContributors;
    if (nc<=0)
      ABORT("level=%d i=%d number_of_contributors=%d\n",level+1,i,nc);
    indh = S->pointIndex[level+1][i];
    the_contributors = ch[i];
    the_weights = Rh[i];	/* entries of the Restriction matrix */
    switch (S->transferchoice) {
    case correct: /* upper level calculates correction to the lower level */
      rsum = 0;
      for (ic=0;ic<nc;ic++) {
	the_contributor = the_contributors[ic];
	the_weight = the_weights[ic];
	rsum += the_weight*rl[the_contributor];
      }
      fh[indh] = rsum;	/* higher-level free term from lower level's residual */
      uh[indh] = 0;	/* initial guess for correction is zero */
      rh[indh] = 0;	/* this should not be needed ? */
      break;
    case replace: /* upper level calculates coarser version of the lower level */
      usum = 0;
      fsum = 0;
      for (ic=0;ic<nc;ic++) {
	the_contributor = the_contributors[ic];
	the_weight = the_weights[ic];
	fsum += the_weight*S->fgrid[level][the_contributor];
	usum += the_weight*S->ugrid[level][the_contributor];
      }
      S->fgrid[level+1][indh] = fsum;	/* higher-level free term from lower level free term */
      S->ugrid[level+1][indh] = usum;	/* higher-level initial guess from lower level solution */
      S->rgrid[level+1][indh] = 0;	/* this should not be needed ? */
      break;
    default:
      ABORT("Unknown transfer method choice %d\n",S->transferchoice);
    } /* switch transfer */
  } /* for i */

  PROFILE(ALLELSE);

  return 1;
}
/*=========================================*/

/* Implement the Prolongation Operator: get the correction from coarse system and add it */
/* P: (level+1) -> (level) */
int down_level(int level,STR *S) {
  list_t number_of_coarse_points, number_of_fine_points, i;
  contr_t nc, ic;
  index_t *the_contributors, the_contributor;
  sindex_t indl, indh; /* may happen to be negative */
  real *the_weights;
  real the_weight;
  real increment;

  PROFILE(DOWNLEVEL);

  /* Prolongation of the upper level will overwrite the lower level */
  if (S->transferchoice==replace) {
    number_of_fine_points=S->nPoints[level];
    for (i=1;i<=number_of_fine_points;i++) { /* ENUMERATION STARTS WITH 1 !!! */
      indl = S->pointIndex[level][i];
              if (indl<0 || indl>=S->levelsize[level]) {
        	ABORT("\nwrong index %ld of point %ld at level %ld\n",
        	       (long)indl, (long)i, (long)level);
              }
      S->ugrid[level][indl]=0;
    }
  }

  number_of_coarse_points=S->nPoints[level+1];
  for (i=1;i<=number_of_coarse_points;i++) { /* ENUMERATION STARTS WITH 1 !!! */
    indh = S->pointIndex[level+1][i];
              if (indh<0 || indh>=S->levelsize[level+1]) {
        	ABORT("\nwrong index %ld of point %ld at level %ld\n",
        	       (long)indh, (long)i, (long)level+1);
              }
    increment = S->ugrid[level+1][indh];
    nc = S->numContributors;
    ASSERT(nc>0);
    the_contributors=S->contributors[level+1][i];
    the_weights = S->P[level+1][i];	/* entries of the Prolongation matrix */
    for(ic=0;ic<nc;ic++) {
      the_contributor = the_contributors[ic];
              if (/* the_contributor<0 && */ the_contributor>=S->levelsize[level]) {
		ABORT("\nwrong contributor index %ld of contributor %ld "
		       "of point %ld of level %ld\n",
        	       (long)the_contributor,(long)ic,(long)i,(long)level+1);
              }
      the_weight = the_weights[ic];
      /* This adds to 0 (replace) or to previous contents (correct) */
      S->ugrid[level][the_contributor] += increment*the_weight;
    }
  }

  PROFILE(ALLELSE);
  return 1;
}
/*=========================================*/


/* Perform one iteration loop and return the max norm of the residual.         */
/* Normally relax=1; if relax=0, just measure the residuals;                   */
/* Take relax>1 for SOR - but NR do not recommend it for smoothing operators.  */
static int direction=1;
static real iterate_level(int level,real relax,STR *S) {
  /* Rules of thumb for coarsegraining	:
     solution (u) is intensive variable;
     forcing term (f) is extensive variable;
     matrix elements==capacitance and connections (K) are extensive variables.
     Interpolation/restriction operators: 2/3-linear.
  */
  index_t indh;
  list_t i;
  neighbour_t numNeighbours=S->numNeighbours;
  sindex_t *d=S->d[level];
  real *u, *r, *K;
  real f, sum, unew, change, usum, avg;
  neighbour_t nb;
  AnisoPoint P;
  double absres, resnorm, globalresnorm; /* fixed type for MPI_REDUCE */
  int istart, inotreach, istep;
  /* By-hand optimization */
  real *ul=S->ugrid[level];
  real *fl=S->fgrid[level];
  real *rl=S->rgrid[level];
  real rvalue;
  index_t *pil=S->pointIndex[level];
  AnisoPoint *gpl=S->gridPoints[level];

  PROFILE(level);
#if MPI
  /* Synchronise the subdomains at this point */
  halo_level(level,S->ugrid[level],S);
#endif	

  /* Will measure the progress of the iteration */
  resnorm=0.0; 

  /* Will count the average in case we need to zero it */
  usum=0.0;

  /* TODO: special case if only residuals are required */

  switch (S->smootherchoice) {
  /*=======================*/
  case GS: /* Gauss-Seidel */
    if (S->alternate) direction=-direction;

    if (direction>0) {	/* ENUMERATION STARTS WITH 1 !!! */
      istart=1;
      inotreach=S->nPoints[level]+1;
      istep=1;
    } else {		/* ENUMERATION STARTS WITH 1 !!! */
      istart=S->nPoints[level];
      inotreach=0;
      istep=-1;
    }

    /* Loop through all points */
    for (i=istart;i!=inotreach;i+=istep) {
      P = gpl[i];				/* This point */
      K = &(P.K[0]);				/* This point's connections */
      indh = pil[i];				/* This point's index in the grid */
      u=ul+indh;				/* The pointer to solution */
      r=rl+indh;				/* The pointer to residual */
      f=fl[indh];				/* The value of the force */
      sum=0;
      for (nb=0;nb<numNeighbours;nb++)
	sum += u[d[nb]]*K[nb]; 			/* the LHS of the system */
      rvalue=f-sum;				/* the residual */ 
      if (rvalue!=rvalue) {
	URGENT_MESSAGE("\nNAN (not-a-number) detected in #%d at t=%ld level=%ld indh=%ld",mpi_rank,t,(long)level,(long)indh);
	URGENT_MESSAGE("\nThis happened to rvalue for f=%g sum=%g\n",f,sum);
	ABORT("");
      }
      unew=(*u)+relax*rvalue/K[0];		/* updated value */
      if (unew!=unew) {
	URGENT_MESSAGE("\nNAN (not-a-number) detected in #%d at t=%ld level=%ld indh=%ld",mpi_rank,t,(long)level,(long)indh);
	URGENT_MESSAGE("\nThis happened to unew for uold=%g relax=%g rvalue=%g K[0]=%g\n",*u,relax,rvalue,K[0]);
	ABORT("");
      }
      absres=fabs(rvalue);			/* c_0 norm of */
      if (absres>resnorm) resnorm=absres;	/*     the residual */
      (*u)=unew;
      (*r)=rvalue;
      usum+=unew;
    } /* for i */

    break; /* end of case GS */

  /*=====================*/
  case Jacobi: /* Jacobi */
    /* First calc residuals */
    for (i=1;i<=S->nPoints[level];i++) { /* ENUMERATION STARTS WITH 1 !!! */
      P = S->gridPoints[level][i];		/* This point */
      K = &(P.K[0]);				/* This point's connections */
      indh = S->pointIndex[level][i];		/* This point's index in the grid */
      u=&(S->ugrid[level][indh]);		/* The pointer to solution */
      r=&(S->rgrid[level][indh]);		/* The pointer to residual */
      f=S->fgrid[level][indh];			/* The value of the force */
      for (sum=nb=0;nb<numNeighbours;nb++) sum += u[d[nb]]*K[nb]; /* the LHS of the system */
      (*r)=f-sum;				/* the residual */ 
      absres=fabs(*r);				/* c_0 norm of */
      if (absres>resnorm) resnorm=absres;	/*     the residual */
    } /* for i */
    /* Then update the solution */
    if (relax) {
      for (i=1;i<=S->nPoints[level];i++) { /* ENUMERATION STARTS WITH 1 !!! */
	P = S->gridPoints[level][i];		/* This point */
	K = &(P.K[0]);				/* This point's connections */
	indh = S->pointIndex[level][i];		/* This point's index in the grid */
	u=&(S->ugrid[level][indh]);		/* The pointer to solution */
	r=&(S->rgrid[level][indh]);		/* The pointer to residual */
	(*u)=(*u)+relax*(*r)/K[0];		/* updated value */
	usum+=*u;
      } /* for i */
    } /* if relax */
    break; /* end of case Jacobi */

  /*=====*/
  default:
    EXPECTED_ERROR("Unknown smoother choice %d\n",S->smootherchoice);
  }

  /* Go on calculating max norm of residual across subdomains */
#if MPI
  MPIDO(MPI_Reduce(&resnorm,&globalresnorm,1,MPI_DOUBLE,MPI_MAX, S->root, S->comm),
	"Couldn't carry out reduction operation.");
  MPIDO(MPI_Bcast(&globalresnorm, 1, MPI_DOUBLE, S->root, S->comm),
	"Couldn't broadcast reduction result.");
#else
  globalresnorm=resnorm;
#endif
  S->rsn[level]=globalresnorm;
  if (relax!=0) S->iter[level]++;

  PROFILE(ALLELSE);
  if (S->profile) {
    DEBUG(" level=%d iter=%d%s [%gms] resnorm=%lg\n",
	  level, S->iter[level], relax?"":"*",
	  S->profile_totalspent[level],
	  (double)globalresnorm);
  } else {
    DEBUG(" level=%d iter=%d%s resnorm=%lg\n", level, S->iter[level], relax?"":"*",(double)globalresnorm);
  }

  return globalresnorm;
}
/*=========================================*/

/* Just calculate the residuals and save and return their max norm. */
/* Made a separate function, to be used for debugging the main code. */
static real residual_level(int level,STR *S) {
  index_t indh;
  list_t i;
  neighbour_t numNeighbours=S->numNeighbours;
  sindex_t *d=S->d[level];
  real *u, *r, *K;
  real f, sum;
  neighbour_t nb;
  AnisoPoint P;
  double absres, resnorm, globalresnorm; /* fixed type for MPI_REDUCE */

  PROFILE(level);
#if MPI
  /* Synchronise the subdomains at this point */
  halo_level(level,S->ugrid[level],S);
#endif	

  /* Will measure the progress of the iteration */
  resnorm=0.0; 

  for (i=1;i<=S->nPoints[level];i++) { /* ENUMERATION STARTS WITH 1 !!! */
    P = S->gridPoints[level][i];		/* This point */
    K = &(P.K[0]);				/* This point's connections */
    indh = S->pointIndex[level][i];		/* This point's index in the grid */
    u=&(S->ugrid[level][indh]);			/* The pointer to solution */
    r=&(S->rgrid[level][indh]);			/* The pointer to residual */
    f=S->fgrid[level][indh];			/* The value of the force */
    for (sum=nb=0;nb<numNeighbours;nb++) 
      sum += u[d[nb]]*K[nb]; 			/* the LHS of the system */
    (*r)=f-sum;					/* the residual */ 
    absres=fabs(*r);				/* c_0 norm of */
    if (absres>resnorm) resnorm=absres;		/*     the residual */
  } /* for i */

  /* Go on calculating max norm of residual across subdomains */
#if MPI
  MPIDO(MPI_Reduce(&resnorm,&globalresnorm,1,MPI_DOUBLE,MPI_MAX, S->root, S->comm),
	"Couldn't carry out reduction operation.");
  MPIDO(MPI_Bcast(&globalresnorm, 1, MPI_DOUBLE, S->root, S->comm),
	"Couldn't broadcast reduction result.");
#else
  globalresnorm=resnorm;
#endif
  S->rsn[level]=globalresnorm;

  return globalresnorm;
}
/*=========================================*/

#if MPI
static int halo_level(int level,real *grid,STR *S) {
  int nb; /* rank of a neighbour */
  /* X AXIS */
  if (mpi_ix % 2 == 0) {
    if (0<=(nb=XP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->xp_type[level],nb,0,
			 grid,1,S->xp_halo_type[level],XP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XP neighbour."); 
    }
  } else { 
    if (0<=(nb=XN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->xn_type[level],nb,0,
			 grid,1,S->xn_halo_type[level],XN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XN neighbour."); 
    }
  }

  if (mpi_ix % 2 == 1) { 
    if (0<=(nb=XP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->xp_type[level],nb,0,
			 grid,1,S->xp_halo_type[level],XP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XP neighbour."); 
    } 
  } else { 
    if (0<=(nb=XN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->xn_type[level],nb,0,
			 grid,1,S->xn_halo_type[level],XN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XN neighbour.");
    } 
  }

  /* Y AXIS */ 
  if (mpi_iy % 2 == 0) { 
    if (0<=(nb=YP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->yp_type[level],nb,0,
			 grid,1,S->yp_halo_type[level],YP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YP neighbour.");
    } 
  } else { 
    if (0<=(nb=YN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->yn_type[level],nb,0,
			 grid,1,S->yn_halo_type[level],YN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YN neighbour.");
    }
  } 

  if (mpi_iy % 2 == 1) { 
    if (0<=(nb=YP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->yp_type[level],nb,0,
			 grid,1,S->yp_halo_type[level],YP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YP neighbour.");
    } 
  } else { 
    if (0<=(nb=YN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->yn_type[level],nb,0,
			 grid,1,S->yn_halo_type[level],YN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YN neighbour.");
    } 
  }

  /* Z AXIS */
  if (mpi_iz % 2 == 0) { 
    if (0<=(nb=ZP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->zp_type[level],nb,0,
			 grid,1,S->zp_halo_type[level],ZP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZP neighbour.");
    } 
  } else { 
    if (0<=(nb=ZN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->zn_type[level],nb,0,
			 grid,1,S->zn_halo_type[level],ZN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZN neighbour.");
    } 
  }

  if (mpi_iz % 2 == 1) { 
    if (0<=(nb=ZP_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->zp_type[level],nb,0,
			 grid,1,S->zp_halo_type[level],ZP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZP neighbour.");
    }
  } else { 
    if (0<=(nb=ZN_NEIGHBOUR)) {
      MPIDO(MPI_Sendrecv(grid,1,S->zn_type[level],nb,0,
			 grid,1,S->zn_halo_type[level],ZN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZN neighbour.");
    }
  }

  return 1;

}
/*=========================================*/
#endif

#if MPI
static void define_halo_level(int level,STR *S) {
  /*	Datatypes for halo swapping for the given level.
   *    As decomp_defineHaloTypes() in decomp.c, 
   *	adjusted for the multigrid layers, 
   *    except it is just for u-grid, i.e. vmax->1.
   *
   *	Correspondence: basic (local coord) vs grid
   *	xlen			vs	S->xlen[level]
   *	ylen			vs	S->ylen[level]
   *	zlen			vs	S->zlen[level]
   *	vmax			vs	1
   *
   *	local_xmax-local_xmin	vs	S->nx[level]
   *	local_ymax-local_ymin	vs	S->ny[level]
   *	local_zmax-local_zmin	vs	S->nz[level]
   *
   *	ONE			vs	ONE
   *	TWO			vs	TWO
   *	TRI			vs 	TRI
   *
   *	1			vs	1
   *
   *	(local_xmax-1+ONE)-local_xmin	vs S->nx[level]-1+ONE=S->xmax[level]
   *	(local_ymax-1+TWO)-local_ymin	vs S->ny[level]-1+TWO=S->ymax[level]
   *	(local_zmax-1+TRI)-local_zmin	vs S->nz[level]-1+TRI=S->zmax[level]
   *
   *	-----------------------------------	*/
  #define NDIMS 3
  int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
  sizes[0] = S->xlen[level];
  sizes[1] = S->ylen[level];
  sizes[2] = S->zlen[level];
  /* sizes[3] = 1; */

  /*=============*/
  if (mpi_nx > 1) {

    /*	XN	*/
    starts[0] = 1; 	subsizes[0] = 1;		/* We know ONE=1 here */
    starts[1] = TWO; 	subsizes[1] = S->ny[level];
    starts[2] = TRI;	subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
    				   MPI_DOUBLE, &(S->xn_type[level])),
    	  "Couldn't define xn_type.");
    MPIDO(MPI_Type_commit(&(S->xn_type[level])),
	  "Couldn't commit xn_type.");

    /*	XN Halo	*/
    starts[0] = 0; 	subsizes[0] = 1;
    starts[1] = TWO; 	subsizes[1] = S->ny[level];
    starts[2] = TRI;	subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->xn_halo_type[level])),
	  "Couldn't define xn_halo_type.");
    MPIDO(MPI_Type_commit(&(S->xn_halo_type[level])),
	  "Couldn't commit xn_halo_type.");

    /*	XP	*/
    starts[0] = S->xmax[level];		subsizes[0] = 1;
    starts[1] = TWO; 			subsizes[1] = S->ny[level];
    starts[2] = TRI;			subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->xp_type[level])),
	  "Couldn't define xp_type.");
    MPIDO(MPI_Type_commit(&(S->xp_type[level])),
	  "Couldn't commit xp_type.");

    /*	XP Halo	*/
    starts[0] = S->xmax[level]+1;	subsizes[0] = 1;
    starts[1] = TWO; 			subsizes[1] = S->ny[level];
    starts[2] = TRI;			subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->xp_halo_type[level])),
	  "Couldn't define xp_halo_type.");
    MPIDO(MPI_Type_commit(&(S->xp_halo_type[level])),
	  "Couldn't commit xp_halo_type.");
  } /*	if(mpi_nx > 1)	*/
  
  /*=============*/
  if (mpi_ny > 1) {
    /*	YN	*/
    starts[0] = 0; 	subsizes[0] = S->xlen[level];	/* Widened in x to include magic corners */
    starts[1] = 1; 	subsizes[1] = 1;		/* We know TWO=1 here */
    starts[2] = TRI;	subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
    				   MPI_DOUBLE, &(S->yn_type[level])),
    	  "Couldn't define yn_type.");
    MPIDO(MPI_Type_commit(&(S->yn_type[level])),
	  "Couldn't commit yn_type.");
	
    /*	YN Halo		*/
    starts[0] = 0; 	subsizes[0] = S->xlen[level];	/* Widened in x to include magic corners */
    starts[1] = 0; 	subsizes[1] = 1;		/* We know TWO=1 here */
    starts[2] = TRI;	subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->yn_halo_type[level])),
	  "Couldn't define yn_halo_type.");
    MPIDO(MPI_Type_commit(&(S->yn_halo_type[level])),
	  "Couldn't commit yn_halo_type.");
      
    /*	YP	*/
    starts[0] = 0; 			subsizes[0] = S->xlen[level];	/* Widened in x to include magic corners */
    starts[1] = S->ymax[level];		subsizes[1] = 1;
    starts[2] = TRI;			subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->yp_type[level])),
	  "Couldn't define yp_type.");
    MPIDO(MPI_Type_commit(&(S->yp_type[level])),
	  "Couldn't commit yp_type.")

    /*	YP Halo	*/
    starts[0] = 0; 			subsizes[0] = S->xlen[level];	/* Widened in x to include magic corners */
    starts[1] = S->ymax[level]+1;	subsizes[1] = 1;		
    starts[2] = TRI;			subsizes[2] = S->nz[level];
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->yp_halo_type[level])),
	  "Couldn't define yp_type.");
    MPIDO(MPI_Type_commit(&(S->yp_halo_type[level])),
	  "Couldn't commit yp_type.");
  } /*	if(mpi_ny > 1)	*/

  /*=============*/
  if (mpi_nz > 1) {
    /*	ZN	*/
    starts[0] = 0;	subsizes[0] = S->xlen[level];		/* Widened in x to include magic corners */
    starts[1] = 0;	subsizes[1] = S->ylen[level];		/* Widened in y to include magic corners */
    starts[2] = 1;	subsizes[2] = 1;			/* We know TRI=1 here */
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C, 
				   MPI_DOUBLE, &(S->zn_type[level])),
	  "Couldn't define zn_type.");
    MPIDO(MPI_Type_commit(&(S->zn_type[level])),
	  "Couldn't commit zn_type.");
	
    /*	ZN Halo	*/
    starts[0] = 0;	subsizes[0] = S->xlen[level];		/* Widened in x to include magic corners */
    starts[1] = 0;	subsizes[1] = S->ylen[level];		/* Widened in y to include magic corners */
    starts[2] = 0;	subsizes[2] = 1;			/* We know TRI=1 here */
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->zn_halo_type[level])),
	  "Couldn't define zn_halo_type.");
    MPIDO(MPI_Type_commit(&(S->zn_halo_type[level])),
	  "Couldn't commit zn_halo_type.");

    /*	ZP	*/
    starts[0] = 0;			subsizes[0] = S->xlen[level];	/* Widened in x to include magic corners */
    starts[1] = 0;			subsizes[1] = S->ylen[level];	/* Widened in y to include magic corners */
    starts[2] = S->zmax[level];		subsizes[2] = 1;		/* We know TRI=1 here */
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->zp_type[level])),
	  "Couldn't define zp_type.");
    MPIDO(MPI_Type_commit(&(S->zp_type[level])),
	  "Couldn't commit zp_type.");
      
    /*	ZP Halo	*/
    starts[0] = 0;			subsizes[0] = S->xlen[level];	/* Widened in x to include magic corners */
    starts[1] = 0;			subsizes[1] = S->ylen[level];	/* Widened in y to include magic corners */
    starts[2] = S->zmax[level]+1;	subsizes[2] = 1;		/* We know TRI=1 here */
    MPIDO(MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C,
				   MPI_DOUBLE, &(S->zp_halo_type[level])),
	  "Couldn't define zp_halo_type.");
    MPIDO(MPI_Type_commit(&(S->zp_halo_type[level])),
	  "Couldn't commit zp_halo_type.");
  } /*	if(mpi_nz > 1)	*/

  #undef NDIMS
}
/*=========================================*/
#endif

static inline double Wtime(void) {
  #if MPI
    return MPI_Wtime();
  #else
    return (double)clock()/(double)CLOCKS_PER_SEC;
  #endif
}
/*=========================================*/

#if MPI
static void profile_start(STR *S) {};
static void profile_level(int level,STR *S) {};
static void profile_report(STR *S) {};
#else

static void profile_start(STR *S) {
  int level;
  if (!(S->profile)) return;
  for (level=0;level<TOTALSPENTNUM;level++)
    S->profile_totalspent[level]=0.0;
  S->profile_lastlevel=-1;
}
/*=========================================*/

static void profile_level(int level,STR *S) {
  double now, spent;
  int last=S->profile_lastlevel;
  if (level==last) return;
  now=Wtime();
  spent=now-S->profile_lasttime;
  /* level<0 used to signal begin and end of profiling session */
  if (last>=0) S->profile_totalspent[last]+=spent;
  S->profile_lastlevel=level;
  S->profile_lasttime=now;
}
/*=========================================*/

static void profile_report(STR *S) {
  DEVICE_CONST(FILE *,profile);
  DEVICE_CONST(int,profileWriter);
  DEVICE_CONST(int,upper_level);
  DEVICE_ARRAY(double,profile_totalspent);
  int level;
  double k=1.e3, spent, total;
  #define FMT "\t%5.3f"
  if (!profile) return;
  PROFILE(-1);
  if (profileWriter) {
    total=0.0;
    fprintf(profile,"%-5ld: elliptic profile:",(long)t);
    for (level=upper_level+1;level<TOTALSPENTNUM;level++) {
      spent=profile_totalspent[level];
      fprintf(profile,FMT,k*spent);
      total+=spent;
    }
    fprintf(profile," | ");
    for (level=0;level<=upper_level;level++) {
      spent=profile_totalspent[level];
      fprintf(profile,FMT,k*spent);
      total+=spent;
    }
    fprintf(profile," | ");
    fprintf(profile,FMT,k*total);
    fprintf(profile,"\n");
  }

  for (level=0;level<TOTALSPENTNUM;level++)
    profile_totalspent[level]=0.0;

  S->profile_lastlevel=-1;
}
#endif
/*=========================================*/
