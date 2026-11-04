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

/* multigrid_code.h: the multigrid solver of the elliptic device, shared with
 * other devices that use it. Like diffweights_code.h, it is textually
 * included in the device's source, after its STR typedef (which includes
 * diffweights_str.h and multigrid_str.h) and after diffweights_code.h.
 *
 * The device's create function calls, in this order:
 *   mg_accept(dev,w,S,isodef) - the space checks and the solver's parameters,
 *                          including the connection weights' ones (diff_weights_accept);
 *                          isodef is the default of isolated: 0 or 1, or -1 for the
 *                          default of diff (1 if anisotropic or variable diffusivity);
 *   [accept v0 (layer of the RHS) and v1 (layer of the solution and initial guess),
 *    and set S->tau if the operator is I-tau*L rather than L, see multigrid_str.h]
 *   mg_connect(dev,w,S)  - the connection weights, the list of points, the communicator,
 *                          the singular problem detection;
 *   mg_build(dev,w,S)    - the multigrid parameters and the hierarchy of grids;
 * its run function calls mg_solve(S), and its destroy function calls mg_destroy(S)
 * and frees the fields of diffweights_str.h.
 *
 * The functions return 0 on failure, so that the create function can return 0;
 * NB this file redefines DEBUG.
 */

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
/* Whether a point (global coords) is in the device space: the box, and also restricted */
/* and where=, so that the tissue points excluded by where= give Dirichlet values.       */
static int mg_in_space(STR *S,long x,long y,long z) {
  Space *s=S->sp;
  if (x<s->global_x0 || x>s->global_x1 || y<s->global_y0 || y>s->global_y1 ||
      z<s->global_z0 || z>s->global_z1) return 0;
  return in_space(s,(INT)x,(INT)y,(INT)z);
}
/* The singular problem */
static int is_singular(STR *S);
static real points_mean(real *grid,STR *S);
static real layer_mean(int v,STR *S);
static void apply_pin(STR *S);

/* Prepare multigrid stuff for given level (and above, recursively) */
static int make_level(int level,STR *S);

/* Find the parity of the coarsening of the given level */
static void set_anchor(int level,STR *S);
/* Coarse parents of a fine point, and their kinds */
static int fine_parents(int level,const long f[3],long k[8][3],real w[8],STR *S);
static int parent_kind(int level,const long k[3],STR *S);
static real fine_norm(int level,const long f[3],STR *S);

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
/* And vice-versa */
/* Fine (level-1) local coordinate of a coarse (level) point; uses 'level' in scope */
#define finex(xh) (2*(xh)+S->off[0][level])
#define finey(yh) (2*(yh)+S->off[1][level])
#define finez(zh) (2*(zh)+S->off[2][level])
/* Integer floor division by 2 and remainder mod 2, correct for negative numbers */
#define mod2(a) ((((a)%2)+2)%2)
#define fdiv2(a) (((a)-mod2(a))/2)

#undef DEBUG
#define DEBUG(...) if (S->debug&&S->debugWriter) {fprintf(S->debug,__VA_ARGS__);fflush(S->debug);}
#define PROFILE(...) if (S->profile&&S->profileWriter) {profile_level(__VA_ARGS__,S);}

/* END OF DEFINITIONS                                    */

/*==============================================*/
/* The interface functions                      */

/* Solve the equation with the RHS in layer v0, the solution and its initial guess in layer v1 */
static int mg_solve(STR *S)
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
    S->fmean=(S->singular)?layer_mean(v0,S):0.0; /* the singular problem: compatible RHS */
    for (it=0;it<maxiter;it++)
      if (tolerance>iterate_base(it,S))
	break;
    DEBUG("result: iter=%ld resnorm=%lg\n", (long)it, (double)S->resnorm);
    if (S->singular) {			/* ... and zero-mean solution */
      real mean=layer_mean(v1,S);
      for (i=1;i<=numPoints;i++) New[ind(points[i].x,points[i].y,points[i].z,v1)]-=mean;
    }
  }

  /* The singular problem with pin=1: shift the solution so that u(pin)=upin */
  if (S->singular && S->pin) apply_pin(S);

  profile_report(S);
  return 1;
}
/*=========================================*/

/* Free what was allocated by mg_connect and mg_build */
static void mg_destroy(STR *S)
{
  if (S->points) {
    list_t i;
    for (i=1;i<=S->numPoints;i++) FREE(S->points[i].K);
  }
  FREE(S->points);
  SAFE_CLOSE(S->debug);
  SAFE_CLOSE(S->profile);
#if MPI
  FREE(S->active);
  /* Halo types of the levels that have been made */
  {
    MPI_Datatype **types[] = {
      &S->xn_type, &S->xp_type, &S->yn_type, &S->yp_type, &S->zn_type, &S->zp_type,
      &S->xn_halo_type, &S->xp_halo_type, &S->yn_halo_type, &S->yp_halo_type,
      &S->zn_halo_type, &S->zp_halo_type};
    int i, l;
    for (i=0; i<sizeof(types)/sizeof(types[0]); i++) if (*types[i]) {
      for (l=0; l<=S->upper_level; l++)
	if ((*types[i])[l]!=MPI_DATATYPE_NULL) MPI_Type_free(&((*types[i])[l]));
      FREE(*types[i]);
    }
  }
#endif
  {
    int d;
    for (d=0;d<3;d++) {
      FREE(S->off[d]); FREE(S->g0[d]); FREE(S->glo[d]); FREE(S->ghi[d]); FREE(S->anchor[d]);
    }
    if (S->parN) {
      int l;
      for (l=0;l<=S->upper_level;l++) {FREE(S->parN[l]); FREE(S->parIdx[l]); FREE(S->parW[l]);}
      FREE(S->parN); FREE(S->parIdx); FREE(S->parW);
    }
  }
  /* TODO: free all other multigrid objects ... */
}
/*=========================================*/

/* The space checks and the parameters accepted before the layers */
static int mg_accept(Device *dev,char *w,STR *S,int isodef)
{
  Space *s=&(dev->s);
  int hasexe;				/* whether any diffusivity is a =~ expression */

  if (s->nowhere) EXPECTED_ERROR("'nowhere' is not allowed\n");
  if (find_key("listpoints=",w)) {	/* say something only if the user has a say */
    if (s->listpoints==0)
      MESSAGE("/* listpoints=1 is mandatory in this device and enforced */\n");
    else
      MESSAGE("/* listpoints=1 is mandatory anyway */\n");
  }
  if (!s->listpoints) {			/* otherwise already listed by accept_space() */
    s->listpoints=1;
    s->p=NULL;				/* signal that it has not be alloc'd yet */
    if (SUCCESS != list_space(s))
      ABORT("could not make list of points");
  }
  
  ACCEPTF(debug,"wt","");
  ACCEPTF(profile,"wt","");
#if MPI
  if (S->profile) MESSAGE("\n/* profiling not implemented for MPI; will proceed without profiling */");
#endif	
  /* The connection weights: the same parameters and code as in diff, see diffweights_code.h */
  ACCEPTR(hx,RNONE,RSUCC(0.),RNONE);
  if (!diff_weights_accept(S,w,&hasexe)) return 0;
  /* The border of the device space: by default, connected to the tissue beyond */
  /* it, whose values are then the Dirichlet boundary conditions; with          */
  /* isolated=1, it is a no-flux boundary, as is the border of the tissue.      */
  /* The default is the device's choice; isodef<0 means the default of diff.    */
  if (isodef<0) isodef=(S->loctb || hasexe || ANISOTROPY_ON)?1:0;
  ACCEPTI(isolated,isodef,0,1);
  S->w0=S->w1=-1;			/* the weights are kept in the device's pool */
  ACCEPTR(tolerance,RNONE,0.,RNONE);
  ACCEPTI(maxiter,INONE,1,INONE);
  /* Normalization of the solution of the singular problem (no-flux at all borders): */
  /* zero mean, or, if pin=1, u(xpin,ypin,zpin)=upin. Applied after the solution.    */
  ACCEPTI(pin,0,0,1);
  if (S->pin) {
    /* Global limits: in MPI, the pin point is in one subdomain only */
    ACCEPTI(xpin,INONE,dev->s.global_x0,dev->s.global_x1);
    ACCEPTI(ypin,INONE,dev->s.global_y0,dev->s.global_y1);
    ACCEPTI(zpin,INONE,dev->s.global_z0,dev->s.global_z1);
    ACCEPTR(upin,RNONE,RNONE,RNONE);
  }
  return 1;
}
/*=========================================*/

/* The connections, the points and the communicator; NB v0, v1 accepted by now */
static int mg_connect(Device *dev,char *w,STR *S)
{
  Space *s=&(dev->s);
  S->s=dev->s;
  S->sp=&(dev->s);
  
  /* S->started=0; */

  /* The connection weights at the listed points, and the device's list of points */
  if (SUCCESS != diff_connect(S,s)) ABORT("could not make the connection weights");
  if (!list_points(S,*s)) return 0;
  /* The operator I-tau*L instead of L. NB after list_points(), which leaves  */
  /* out the points without connections: there the solution stays unchanged. */
  if (S->tau) {
    list_t i;
    neighbour_t nb;
    for (i=1;i<=S->numPoints;i++) {
      for (nb=0;nb<maxnb;nb++) S->points[i].K[nb]*=-S->tau; /* NB numNeighbours not yet known */
      S->points[i].K[ooo]+=1.0;
    }
  }

#if MPI
  /* We need one process that is guaranteed to be engaged in this device */
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    ABORT("Could not create communicator.\n");
  /* Identifying the root process allows only one process to write debug output. */
  S->debugWriter = (mpi_rank == S->root);
  S->profileWriter = (mpi_rank == S->root);
  /* Which processes are active: the multigrid's halos are exchanged only between */
  /* them, as an inactive one has no points of the device space, and does not run. */
  CALLOC(S->active,num_active_procs,sizeof(int));
  MPIDO(MPI_Allgather(&(dev->s.runHere),1,MPI_INT,S->active,1,MPI_INT,ALL_ACTIVE_PROCS),
	"Couldn't gather runHere flags.");
#else
  S->debugWriter = 1;
  S->profileWriter = 1;
#endif /* MPI */

  /* No-flux conditions at all borders? Then the solution is up to a constant. */
  /* NB in MPI this is collective over the active instances only; runHere is   */
  /* not set in sequential mode, where the instance is always active.          */
#if MPI
  int active=dev->s.runHere;
#else
  int active=1;
#endif
  if (active) S->singular=S->tau?0:is_singular(S); /* I-tau*L is never singular */
  if (S->pin && active) {
    long found=0;
    list_t i;
    if (!S->singular)
      EXPECTED_ERROR("pin=1 is only for the singular problem, with no-flux conditions at all borders;\n"
		     "here the solution is unique, and is fixed by the Dirichlet boundary values.\n");
    for (i=1;i<=S->numPoints;i++)
      if ((long)S->points[i].x==S->xpin && (long)S->points[i].y==S->ypin && (long)S->points[i].z==S->zpin)
	S->pinpoint=i;
    found=(S->pinpoint!=0);
#if MPI
    MPIDO(MPI_Allreduce(MPI_IN_PLACE,&found,1,MPI_LONG,MPI_SUM,S->comm),"Couldn't find the pinning point.");
#endif
    if (!found)
      EXPECTED_ERROR("the pinning point (%d,%d,%d) is not a point of the device space\n",
		     S->xpin,S->ypin,S->zpin);
  }
  if (S->singular) {
    if (S->pin)
      MESSAGE("\n/* No Dirichlet boundary values: the problem is singular. The RHS will be made compatible\n"
	      "   by subtracting its mean, and the solution shifted so that u(%d,%d,%d)=%g. */\n",
	      S->xpin,S->ypin,S->zpin,S->upin);
    else
      MESSAGE("\n/* No Dirichlet boundary values: the problem is singular. The RHS will be made compatible\n"
	      "   by subtracting its mean, and the solution normalized to zero mean. */\n");
    if (!S->conservative)
      MESSAGE("/* WARNING: with conservative=0 the operator is not symmetric, so subtracting the mean\n"
	      "   of the RHS may not make it compatible, and the iterations may stall. */\n");
  }
  return 1;
}
/*=========================================*/

/* The multigrid parameters and the hierarchy of grids */
static int mg_build(Device *dev,char *w,STR *S)
{
  S->direction=1;
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
    ACCEPTI(coarseiter,100,0,INONE);
    ACCEPTR(coarsereduce,1e-2,0.,1.);

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
    CALLOC(S->parN,numLevels,sizeof(int *));
    CALLOC(S->parIdx,numLevels,sizeof(index_t *));
    CALLOC(S->parW,numLevels,sizeof(real *));

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
    CALLOC(S->Dc,numLevels,sizeof(sindex_t *));
    for (l=0;l<3;l++) {
      CALLOC(S->off[l],numLevels,sizeof(int));
      CALLOC(S->g0[l],numLevels,sizeof(long));
      CALLOC(S->glo[l],numLevels,sizeof(long));
      CALLOC(S->ghi[l],numLevels,sizeof(long));
      CALLOC(S->anchor[l],numLevels,sizeof(long));
    }
    
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
    /* so that destroy frees only the types that have been created */
    for (l=0;l<numLevels;l++) {
      S->xn_type[l]=S->xp_type[l]=S->yn_type[l]=S->yp_type[l]=S->zn_type[l]=S->zp_type[l]=MPI_DATATYPE_NULL;
      S->xn_halo_type[l]=S->xp_halo_type[l]=S->yn_halo_type[l]=MPI_DATATYPE_NULL;
      S->yp_halo_type[l]=S->zn_halo_type[l]=S->zp_halo_type[l]=MPI_DATATYPE_NULL;
    }
#endif

#if MPI
    /* An inactive instance has no multigrid to make: it is not collective over */
    /* all processes, only over the active ones (S->comm and halo exchanges).   */
    if (!dev->s.runHere) goto inactive;
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
	if (level>0) { /* Checking R: level-1 -> level maps 1->1 in the interior (R=P^T/2^dim) */
	  for (il=1;il<=S->nPoints[level-1];il++) 	/* set 1s thru lower level */
	    S->fgrid[level-1][ S->pointIndex[level-1][il] ]=One;
	  up_level(level-1,S);				/* apply R operator */
	  maxerr=0;
	  for (ih=1;ih<=S->nPoints[level];ih++){ 	/* check thru higher level */
	    err=fabs(S->fgrid[level][ S->pointIndex[level][ih] ]-One);
	    if (err>maxerr) maxerr=err;
	  } 
	  MESSAGE("\n/* R matrix (%d->%d) is 1:1 to within %g (deviations are expected at borders). */",level-1,level,maxerr);
	} /* if level */
	
	/*---------------------------------------------*/
	if (level>0) { /* Checking P: level -> level-1 maps 1->1, except next to Dirichlet faces */
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
	  MESSAGE("\n/* P matrix (%d->%d) is 1:1 to within %g (deviations are expected next to Dirichlet faces). */",level,level-1,maxerr);
	  if (maxerr>1.e-13) 
	    MESSAGE("\n/*   worst offender: %d (%d,%d,%d)",
		   (int)iwst,(int)gridX(level-1,iwst),
		   (int)gridY(level-1,iwst),(int)gridZ(level-1,iwst));
	} /* if level */
	
	/*--------------------------------------------*/
	if (level>=0) { /* Checking A: level -> level maps 1->0 (1->1 for I-tau*L) */
	  for (ih=1;ih<=S->nPoints[level];ih++) {
	    indh=S->pointIndex[level][ih];
	    S->ugrid[level][indh]=One;			/* set 1 thru u-grid */
	    S->fgrid[level][indh]=S->tau?One:Zero;	/* and A1 thru f-grid */
	  }
	  residual_level(level,S);			/* apply A operator */
	  maxerr=0;
	  for (ih=1;ih<=S->nPoints[level];ih++) {
	    indh=S->pointIndex[level][ih];
	    err=fabs(S->rgrid[level][indh]-Zero);		/* check thru r-grid */
	    if (err>maxerr) maxerr=err;
	  }
	} /* if level */
	MESSAGE("\n/* A matrix (%d->%d) maps 1 to %d to within %g. */",level,level,S->tau?1:0,maxerr);
      } /* for level */
#undef One
#undef Zero
      S->transferchoice=transfer_safe;
    } /* if Verbose debugging of 1-transferrals */
  } /* if upper_level >= 0 */
#if MPI
 inactive: ;
#endif
  return 1;
}

/*=============================================*/

/* Make the list of points of the device space with all their connections:  */
/* the weights computed by diff_connect() at the listed points s.p, as in    */
/* the diff device, rearranged into the 27-neighbour arrays K of this device, */
/* in the multicolour order. Points without connections are left out, as     */
/* their equations would be 0=f.                                             */
static int list_points(STR *S,Space s) {
  int off[a3nnb][3];			/* displacements of the diff stencil */
  int nbtab[3][3][3];			/* this device's neighbour of a displacement */
  int nnb, nb, a, b, c;
  size_t ip;
  list_t i, ndropped=0;
  color_t color;

  nnb=diff_stencil_offsets(ANISOTROPY_ON,S->manypoint,dim,off);
  if (nnb!=S->nnb) ABORT("stencil size mismatch: %d vs %d\n",nnb,S->nnb);
  for (a=0;a<3;a++) for (b=0;b<3;b++) for (c=0;c<3;c++) nbtab[a][b][c]=-1;
  for (nb=0;nb<maxnb;nb++) nbtab[ndx[nb]+1][ndy[nb]+1][ndz[nb]+1]=nb;

  CALLOC(S->points,s.np+1,sizeof(AnisoPoint)); /* ENUMERATION STARTS WITH 1 !!! */
  i=1;
  for (color=0;color<8;color++) {
    for (ip=0;ip<s.np;ip++) {
      devicePoint *Q=s.p+ip;
      real *X=(real *)(Q->X);
      AnisoPoint P;
      if ((Q->x%2)+2*(Q->y%2)+4*(Q->z%2)!=color) continue;
      if (X[0]==0) {			/* no connections */
	ndropped++;
	continue;
      }
      CALLOC(P.K,maxnb,sizeof(real));
      for (nb=0;nb<nnb;nb++)
	P.K[nbtab[off[nb][0]+1][off[nb][1]+1][off[nb][2]+1]]=X[nb];
      P.x=Q->x;
      P.y=Q->y;
      P.z=Q->z;
      S->points[i++]=P;
    } /* for ip */
  } /* for color */
  S->numPoints=i-1;
  if (ndropped)
    MESSAGE("/* %ld point(s) of the device space have no connections and are left out */\n",(long)ndropped);

  /* The weights are now in the list; the pools of diff_connect() are not needed */
  FREE(S->Xpool);
  FREE(S->Ypool);
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
    S->zmin[level]=TRI;    S->zmax[level]=nz-1+TRI;
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
    CALLOC(S->Dc[0],maxcn,sizeof(sindex_t));
    /* For the chosen vertex-centered grid, */
    /* neighbours' and contributors' vectors coincide */
    for (nb=0;nb<maxnb;nb++) S->Dc[0][nb]=S->d[0][nb];
    /* S->Dc[0][PPP] = gridIndex(0,+1,+1,+1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][PPM] = gridIndex(0,+1,+1,-1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][PMP] = gridIndex(0,+1,-1,+1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][PMM] = gridIndex(0,+1,-1,-1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][MPP] = gridIndex(0,-1,+1,+1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][MPM] = gridIndex(0,-1,+1,-1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][MMP] = gridIndex(0,-1,-1,+1) - gridIndex(0,0,0,0); */
    /* S->Dc[0][MMM] = gridIndex(0,-1,-1,-1) - gridIndex(0,0,0,0); */

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

	if (mg_in_space(S,(long)xng,(long)yng,(long)zng)) {
	  P->K[nbh]=Pg.K[nbh];		/* the neighbour is within space, valid connection */
	} else if (Pg.K[nbh]) {		/* its outside space, neighbour's value is nonzero */
	  P->K[nbh]=0;			/* no connection within the grid */
	  S->numBClinks++;		/* but count links to outside values */
	  /* and remember which faces of the space have Dirichlet links */
	  if ((long)xng<(long)s.global_x0) S->dlo[0]=1;
	  if ((long)xng>(long)s.global_x1) S->dhi[0]=1;
	  if ((long)yng<(long)s.global_y0) S->dlo[1]=1;
	  if ((long)yng>(long)s.global_y1) S->dhi[1]=1;
	  if ((long)zng<(long)s.global_z0) S->dlo[2]=1;
	  if ((long)zng>(long)s.global_z1) S->dhi[2]=1;
	} else {			/* outside and no connection outside */
	  P->K[nbh]=0;			/* no connection within the grid either */
	}
      }
    }

    /* Now make the links to the outide values. */
    /* TODO: estimate max numBClinks so could do all work in one go? */
    if (S->numBClinks>0) {
      CALLOC(S->BClinks,S->numBClinks,sizeof(bclink));
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
	  
	  if (!mg_in_space(S,(long)xng,(long)yng,(long)zng)) {
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


    /* Coarsening geometry of level 0 */
    S->g0[0][0]=x0; S->glo[0][0]=s.global_x0; S->ghi[0][0]=s.global_x1;
    S->g0[1][0]=y0; S->glo[1][0]=s.global_y0; S->ghi[1][0]=s.global_y1;
    S->g0[2][0]=z0; S->glo[2][0]=s.global_z0; S->ghi[2][0]=s.global_z1;
#if MPI
    MPIDO(MPI_Allreduce(MPI_IN_PLACE,S->dlo,3,MPI_INT,MPI_MAX,S->comm),"Couldn't reduce dlo.");
    MPIDO(MPI_Allreduce(MPI_IN_PLACE,S->dhi,3,MPI_INT,MPI_MAX,S->comm),"Couldn't reduce dhi.");
#endif

    /* end of level==0 case */
  } else {
    /* Sizes of the subgrid. Take care of the halos: they are not coarsened! */
    /* The coarse points are at the fine global positions with the parity of  */
    /* the anchor, the boundary node of the domain (see set_anchor): so a     */
    /* Dirichlet wall, or the outermost points at a no-flux border, are coarse */
    /* grid positions on every level. The parity is global, which keeps the   */
    /* coarse grids of neighbouring MPI subdomains consistent.                 */
    {
      int dd, nocoarse=0;
      coord_t nn[3];
      /* Can this level be made? Every (active) process must have coarse points, */
      /* otherwise the nearest-neighbour halo exchanges would not work. If not,  */
      /* the previous level is the coarsest, in all processes.                   */
      for (dd=0;dd<dim;dd++) {
	coord_t nf=(dd==0)?S->nx[level-1]:(dd==1)?S->ny[level-1]:S->nz[level-1];
	int o=mod2(S->anchor[dd][level-1]-S->g0[dd][level-1]);
	if ((long)nf-o+1<2) nocoarse=1;
      }
#if MPI
      MPIDO(MPI_Allreduce(MPI_IN_PLACE,&nocoarse,1,MPI_INT,MPI_MAX,S->comm),"Couldn't reduce nocoarse.");
#endif
      if (nocoarse) {
	MESSAGE("\n/* upper_level=%d is too large for this grid%s: some subdomain would have no points at level %d.\n"
		"   The multigrid will use levels 0..%d. */\n",S->upper_level,mpi_size>1?" and decomposition":"",level,level-1);
	S->upper_level=level-1;
	if (S->FMGlevel>S->upper_level) S->FMGlevel=S->upper_level;
	return 1;
      }
      for (dd=0;dd<3;dd++) {
	coord_t nf=(dd==0)?S->nx[level-1]:(dd==1)?S->ny[level-1]:S->nz[level-1];
	if (dd<dim) {
	  long a=S->anchor[dd][level-1];
	  int o=S->off[dd][level]=mod2(a-S->g0[dd][level-1]);
	  nn[dd]=(nf-o+1)/2;
	  S->g0[dd][level]=(S->g0[dd][level-1]+o-a)/2;	/* exact: numerator is even */
	  S->glo[dd][level]=-fdiv2(a-S->glo[dd][level-1]); /* ceil((glo-a)/2) */
	  S->ghi[dd][level]=fdiv2(S->ghi[dd][level-1]-a);	/* floor((ghi-a)/2) */
	} else {
	  S->off[dd][level]=0;
	  nn[dd]=1;
	  S->g0[dd][level]=S->glo[dd][level]=S->ghi[dd][level]=0;
	}
      }
      /* "Netto" dimensions */
      nx=(S->nx[level])=nn[0];
      ny=(S->ny[level])=nn[1];
      nz=(S->nz[level])=nn[2];
    }
    /* "Brutto" dimensions */
    S->xlen[level]=nx+2*ONE;
    S->ylen[level]=ny+2*TWO;
    S->zlen[level]=nz+2*TRI;
    /* "Netto" limits (upper inclusive) */
    S->xmin[level]=ONE; S->xmax[level]=nx-1+ONE;
    S->ymin[level]=TWO; S->ymax[level]=ny-1+TWO;
    S->zmin[level]=TRI; S->zmax[level]=nz-1+TRI;
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
    CALLOC(S->Dc[level],maxcn,sizeof(sindex_t));
    /* For the chosen vertex-centered grid, */
    /* neighbours' and contributors' vectors coincide */
    for (nb=0;nb<maxnb;nb++) S->Dc[level][nb]=S->d[level][nb];

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
	      the_contributor=indl+S->Dc[level-1][ic];	/* contributor's would-be 1D index */
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
      for (yh=0;yh<ny;yh++) {
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
	for (yh=0;yh<ny;yh++) {
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

    /* Prolongation weights. A fine point gets the bi/tri-linear weights wgt[] */
    /* of its coarse parents (1D weight 2 from a coincident parent, 1 from each  */
    /* of two neighbouring ones; 8 in total), normalized by the total weight of  */
    /* the parents that count (see parent_kind): the coarse points in the device */
    /* space, and the phantom ones beyond a Dirichlet face of it, which carry    */
    /* zero correction. So P is the plain linear interpolation towards a         */
    /* Dirichlet wall, and is renormalized next to a no-flux border, where it    */
    /* still reproduces constants. The normalization is geometric rather than    */
    /* accumulated over the local coarse points, so that it is the same on both  */
    /* sides of an MPI subdomain interface, including at the fine halo points.   */
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
	    total_weight=0;
	    for (ic=0;ic<nc;ic++) {			/* loop through all contributors */
	      xcnl=xl+cdx[ic];				/* check that the 3D coords of the contributor */
	      ycnl=yl+cdy[ic];
	      zcnl=zl+cdz[ic];
	      if (
		  -ONE<=xcnl && xcnl<=S->nx[level-1]+ONE &&	/* are within the extended box */
		  -TWO<=ycnl && ycnl<=S->ny[level-1]+TWO &&
		  -TRI<=zcnl && zcnl<=S->nz[level-1]+TRI 
		  ) {
		the_contributor=indl+S->Dc[level-1][ic];		/* contributor's would-be 1D index */
		S->contributors[level][i][ic]=the_contributor;	/* register it for fast access */
		j=pointNumber(level-1,the_contributor);	/* nonzero if tissue (in the halo: the owner's number) */
		if (j) {
		  long f[3];
		  f[0]=(long)xcnl; f[1]=(long)ycnl; f[2]=(long)zcnl;
		  the_weight=wgt[ic]/fine_norm(level,f,S);
		} else {
		  the_weight=0.0;
		}
		S->P[level][i][ic]=the_weight;		/* R is made from P below */
		total_weight+=the_weight;
	      } /* if ic is within extended box */
	    } /* for ic */
	    ASSERT(total_weight>0);
	  } /* if i */
	} /* for zh */
      } /* for yh */
    } /* for xh */

    /* The same P matrix by rows, for prolongation (see down_level): the parents */
    /* of each fine point of this subdomain, which may be in the coarse halo.     */
    CALLOC(S->parN[level],S->nPoints[level-1]+1,sizeof(int));
    CALLOC(S->parIdx[level],8*(S->nPoints[level-1]+1),sizeof(index_t));
    CALLOC(S->parW[level],8*(S->nPoints[level-1]+1),sizeof(real));
    for (i=1;i<=S->nPoints[level-1];i++) {
      AnisoPoint *Q=&(S->gridPoints[level-1][i]);
      long f[3], k[8][3];
      real w[8], norm;
      int n, npar, m=0;
      f[0]=(long)Q->x; f[1]=(long)Q->y; f[2]=(long)Q->z;
      npar=fine_parents(level,f,k,w,S);
      norm=fine_norm(level,f,S);
      for (n=0;n<npar;n++) {
	if (parent_kind(level,k[n],S)>0) {
	  S->parIdx[level][8*i+m]=gridIndex(level,k[n][0],k[n][1],k[n][2]);
	  S->parW[level][8*i+m]=w[n]/norm;
	  m++;
	}
      }
      ASSERT(m>0);
      S->parN[level][i]=m;
    }


    /* Restriction is the transpose of prolongation, R = P^T/2^dim (full       */
    /* weighting in the interior). Then the Galerkin coarse operators RAP      */
    /* inherit symmetry and definiteness, otherwise they can degrade at        */
    /* irregular borders. The scaling must be uniform: a point-dependent       */
    /* scaling of the rows of R, although harmless for this level's coarse    */
    /* correction, acts as a column scaling in the next level's restriction.  */
    {
      real scale=1.0/(1<<dim);
      for (i=1;i<=S->nPoints[level];i++)
	for (ic=0;ic<nc;ic++)
	  S->R[level][i][ic]=S->P[level][i][ic]*scale;
    }

    /* Make the coarse system Q=RAP:                                            */
    /*   Q[I][J] = sum_p R[I][p] sum_q A[p][q] P[q][J],                          */
    /* over the contributors p of the coarse point I and their neighbours q.     */
    /* With MPI, some p are fine halo points: their rows of A are taken from     */
    /* their owners by a halo exchange. Some q may lie beyond the fine halo, but */
    /* P[q][J] is geometric (see fine_norm) and needs no data at q, since A[p][q]*/
    /* is nonzero only if q is tissue. Q[I][J] is stored as K[nbh] with          */
    /* J=I+(ndx,ndy,ndz)[nbh]; |J-I|<=1 in each direction.                        */
    {
      neighbour_t nnb=S->numNeighbours;
      int nbtab[3][3][3], a, b, c;
      real **Ag;
      CALLOC(Ag,nnb,sizeof(real *));
      for (nbl=0;nbl<nnb;nbl++) CALLOC(Ag[nbl],S->levelsize[level-1],sizeof(real));
      for (i=1;i<=S->nPoints[level-1];i++)
	for (nbl=0;nbl<nnb;nbl++)
	  Ag[nbl][S->pointIndex[level-1][i]]=S->gridPoints[level-1][i].K[nbl];
#if MPI
      for (nbl=0;nbl<nnb;nbl++) halo_level(level-1,Ag[nbl],S);
#endif
      for (a=0;a<3;a++) for (b=0;b<3;b++) for (c=0;c<3;c++) nbtab[a][b][c]=-1;
      for (nbh=0;nbh<nnb;nbh++) nbtab[ndx[nbh]+1][ndy[nbh]+1][ndz[nbh]+1]=nbh;

      for (i=1;i<=S->nPoints[level];i++) { /* ENUMERATION STARTS WITH 1 !!! */
	AnisoPoint *Q=&(S->gridPoints[level][i]);
	long I[3];
	I[0]=(long)Q->x; I[1]=(long)Q->y; I[2]=(long)Q->z;
	for (ic=0;ic<nc;ic++) {
	  real Rvalue=S->R[level][i][ic];
	  index_t p;
	  long f[3];
	  if (Rvalue==0) continue;
	  p=S->contributors[level][i][ic];
	  f[0]=(long)finex(Q->x)+cdx[ic];
	  f[1]=(long)finey(Q->y)+cdy[ic];
	  f[2]=(long)finez(Q->z)+cdz[ic];
	  for (nbl=0;nbl<nnb;nbl++) {
	    real Avalue=Ag[nbl][p];
	    long q[3], k[8][3];
	    real w[8], norm;
	    int n, npar;
	    if (Avalue==0) continue;
	    q[0]=f[0]+ndx[nbl]; q[1]=f[1]+ndy[nbl]; q[2]=f[2]+ndz[nbl];
	    npar=fine_parents(level,q,k,w,S);
	    norm=fine_norm(level,q,S);
	    for (n=0;n<npar;n++) {
	      long dx=k[n][0]-I[0], dy=k[n][1]-I[1], dz=k[n][2]-I[2];
	      if (parent_kind(level,k[n],S)<=0) continue;
	      if (dx<-1 || dx>1 || dy<-1 || dy>1 || dz<-1 || dz>1 || nbtab[dx+1][dy+1][dz+1]<0)
		ABORT("level %d point %ld: coarse connection (%ld,%ld,%ld) out of stencil\n",
		      level,(long)i,dx,dy,dz);
	      Q->K[nbtab[dx+1][dy+1][dz+1]]+=Rvalue*Avalue*w[n]/norm;
	    }
	  } /* for nbl */
	} /* for ic */
      } /* for i */
      for (nbl=0;nbl<nnb;nbl++) FREE(Ag[nbl]);
      FREE(Ag);
    }
  } /* level>0 case */

  /* Parity of the coarsening of this level into the next one */
  set_anchor(level,S);

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


/* The anchor of a level: next level's points are at positions of its parity. */
/* At level 0, in each direction, it is the boundary node of the domain, as  */
/* in the classic vertex-centred coarsening: the Dirichlet wall (the position */
/* next to the first or last point), or else the outermost point itself at a */
/* no-flux border. Coarser levels keep the same alignment: the image of the  */
/* fine anchor, which is position 0 in coarse global coordinates.            */
static void set_anchor(int level,STR *S) {
  long amin[3], amax[3];
  list_t i;
  int d;
  if (level>0) {
    for (d=0;d<3;d++) S->anchor[d][level]=(d<dim)?0:S->g0[d][level];
    return;
  }
  for (d=0;d<3;d++) {amin[d]=LONG_MAX; amax[d]=-LONG_MAX;}
  for (i=1;i<=S->nPoints[level];i++) {
    AnisoPoint *Q=&(S->gridPoints[level][i]);
    long g[3];
    g[0]=S->g0[0][level]+(long)Q->x;
    g[1]=S->g0[1][level]+(long)Q->y;
    g[2]=S->g0[2][level]+(long)Q->z;
    for (d=0;d<3;d++) {
      if (g[d]<amin[d]) amin[d]=g[d];
      if (g[d]>amax[d]) amax[d]=g[d];
    }
  }
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,amin,3,MPI_LONG,MPI_MIN,S->comm),"Couldn't reduce anchors.");
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,amax,3,MPI_LONG,MPI_MAX,S->comm),"Couldn't reduce anchors.");
#endif
  for (d=0;d<3;d++) {
    if (d>=dim || amin[d]==LONG_MAX) S->anchor[d][level]=S->g0[d][level];
    else if (S->dlo[d]) S->anchor[d][level]=amin[d]-1;	/* Dirichlet wall at the low end */
    else if (S->dhi[d]) S->anchor[d][level]=amax[d]+1;	/* Dirichlet wall at the high end */
    else S->anchor[d][level]=amin[d];			/* no-flux: the outermost point */
  }
}
/*=========================================*/

/* The coarse parents (at the given level) of the point of level-1 with local */
/* coordinates f: their local coordinates k, possibly in the halo or beyond,  */
/* and their non-normalized weights w, the products of the 1D weights (2 for  */
/* a coincident parent, 1 for each of two neighbouring ones), as in wgt[];    */
/* the weights of all parents sum to 8. Returns the number of parents.        */
static int fine_parents(int level,const long f[3],long k[8][3],real w[8],STR *S) {
  long pk[3][2];
  real pw[3][2];
  int np[3], dd, a, b, c, n=0;
  for (dd=0;dd<3;dd++) {
    if (dd<dim) {
      long r=f[dd]-S->off[dd][level];
      if (mod2(r)==0) {np[dd]=1; pk[dd][0]=r/2; pw[dd][0]=2;}
      else {np[dd]=2; pk[dd][0]=fdiv2(r); pk[dd][1]=fdiv2(r)+1; pw[dd][0]=pw[dd][1]=1;}
    } else {np[dd]=1; pk[dd][0]=0; pw[dd][0]=2;}
  }
  for (a=0;a<np[0];a++) for (b=0;b<np[1];b++) for (c=0;c<np[2];c++) {
    k[n][0]=pk[0][a]; k[n][1]=pk[1][b]; k[n][2]=pk[2][c];
    w[n]=pw[0][a]*pw[1][b]*pw[2][c];
    n++;
  }
  return n;
}
/*=========================================*/

/* The kind of a potential coarse point, at local coordinates k of the level:  */
/*  1: within the device space, so a real coarse point if it is a parent of a  */
/*     tissue point (it belongs to the subdomain that owns its fine position); */
/* -1: beyond a Dirichlet face of the device space: a phantom parent with zero */
/*     correction, which counts in the normalization of P;                     */
/*  0: beyond a no-flux border: does not count.                                */
static int parent_kind(int level,const long k[3],STR *S) {
  int dd, out=0, dir=0;
  for (dd=0;dd<dim;dd++) {
    long G=S->g0[dd][level]+k[dd];	/* global coarse coordinate */
    if (G<S->glo[dd][level]) {out=1; if (S->dlo[dd]) dir=1;}
    else if (G>S->ghi[dd][level]) {out=1; if (S->dhi[dd]) dir=1;}
  }
  return out?(dir?-1:0):1;
}
/*=========================================*/

/* The normalization of the row of P of the point of level-1 at local        */
/* coordinates f: the total weight of its parents that count (kind!=0).      */
static real fine_norm(int level,const long f[3],STR *S) {
  long k[8][3];
  real w[8], norm=0;
  int n, npar=fine_parents(level,f,k,w,S);
  for (n=0;n<npar;n++) if (parent_kind(level,k[n],S)) norm+=w[n];
  return norm;
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
      f=New[ind(x,y,z,v0)]-S->fmean; 
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
      f=New[ind(x,y,z,v0)]-S->fmean; 
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
      f=New[ind(x,y,z,v0)]-S->fmean; 
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

  /* The coarsest level: solve, i.e. iterate until the residual is reduced    */
  /* by the factor coarsereduce (or is below the tolerance), but at most       */
  /* coarseiter times. A fixed small number of sweeps here is not enough when  */
  /* the coarsest grid is not very coarse, or is split between MPI subdomains, */
  /* where the smoother is Jacobi-like across their interfaces.                */
  if (level>0 && level==S->upper_level && S->coarseiter>0) {
    real target=S->coarsereduce*iterate_level(level,0,S);
    if (target<tolerance) target=tolerance;
    for (it=0;it<S->coarseiter;it++)
      if (target>iterate_level(level,S->damping,S)) break;
    DEBUG("result: level=%d (coarsest) sweeps=%d rsn=%lg\n", level, it, (double)S->rsn[level]);
    return 1;
  }

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

  /* The singular problem: make the RHS compatible, i.e. orthogonal to the     */
  /* null space of the (symmetric) operator, the constants: subtract its mean.  */
  if (S->singular) {
    list_t i;
    real mean=points_mean(S->fgrid[0],S);
    for (i=1;i<=S->nPoints[0];i++) S->fgrid[0][S->pointIndex[0][i]]-=mean;
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
  if (S->singular) {			/* the singular problem: normalize to zero mean */
    list_t i;
    real mean=points_mean(S->ugrid[0],S);
    for (i=1;i<=S->nPoints[0];i++) S->ugrid[0][S->pointIndex[0][i]]-=mean;
  }
  for (x=0;x<nx;x++) {
    for (y=0;y<ny;y++) {
      for (z=0;z<nz;z++) {
	indh=gridIndex(0,x,y,z);
	New[ind(x0+x,y0+y,z0+z,v1)]=S->ugrid[0][indh];
	/* NB: do not copy fgrid back to v0: it contains the boundary-link and singular-problem
	   modifications, which would corrupt the RHS layer for the next call or other devices */
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
#if MPI
  /* The contributors of a coarse point include fine halo points */
  if (S->transferchoice==correct) {
    halo_level(level,rl,S);
  } else {
    halo_level(level,ul,S);
    halo_level(level,fl,S);
  }
#endif

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
/* P: (level+1) -> (level). Each fine point of this subdomain takes the contributions   */
/* of its coarse parents (see make_level), which may be in the coarse halo.             */
int down_level(int level,STR *S) {
  list_t number_of_fine_points, i;
  int n, npar;
  index_t indl;
  real sum;
  real *ul=S->ugrid[level];
  real *uh=S->ugrid[level+1];
  int *parN=S->parN[level+1];
  index_t *parIdx=S->parIdx[level+1];
  real *parW=S->parW[level+1];

  PROFILE(DOWNLEVEL);

#if MPI
  halo_level(level+1,uh,S);
#endif

  number_of_fine_points=S->nPoints[level];
  for (i=1;i<=number_of_fine_points;i++) { /* ENUMERATION STARTS WITH 1 !!! */
    indl=S->pointIndex[level][i];
    npar=parN[i];
    for (sum=0,n=0;n<npar;n++)
      sum+=parW[8*i+n]*uh[parIdx[8*i+n]];
    if (S->transferchoice==replace)
      ul[indl]=sum;			/* the upper level overwrites the lower level */
    else
      ul[indl]+=sum;			/* the upper level corrects the lower level */
  }

  PROFILE(ALLELSE);
  return 1;
}
/*=========================================*/


/* Perform one iteration loop and return the max norm of the residual.         */
/* Normally relax=1; if relax=0, just measure the residuals;                   */
/* Take relax>1 for SOR - but NR do not recommend it for smoothing operators.  */
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
    if (S->alternate) S->direction=-S->direction;

    if (S->direction>0) {	/* ENUMERATION STARTS WITH 1 !!! */
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
    if (0<=(nb=XP_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->xp_type[level],nb,0,
			 grid,1,S->xp_halo_type[level],XP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XP neighbour."); 
    }
  } else { 
    if (0<=(nb=XN_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->xn_type[level],nb,0,
			 grid,1,S->xn_halo_type[level],XN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XN neighbour."); 
    }
  }

  if (mpi_ix % 2 == 1) { 
    if (0<=(nb=XP_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->xp_type[level],nb,0,
			 grid,1,S->xp_halo_type[level],XP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XP neighbour."); 
    } 
  } else { 
    if (0<=(nb=XN_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->xn_type[level],nb,0,
			 grid,1,S->xn_halo_type[level],XN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with XN neighbour.");
    } 
  }

  /* Y AXIS */ 
  if (mpi_iy % 2 == 0) { 
    if (0<=(nb=YP_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->yp_type[level],nb,0,
			 grid,1,S->yp_halo_type[level],YP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YP neighbour.");
    } 
  } else { 
    if (0<=(nb=YN_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->yn_type[level],nb,0,
			 grid,1,S->yn_halo_type[level],YN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YN neighbour.");
    }
  } 

  if (mpi_iy % 2 == 1) { 
    if (0<=(nb=YP_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->yp_type[level],nb,0,
			 grid,1,S->yp_halo_type[level],YP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YP neighbour.");
    } 
  } else { 
    if (0<=(nb=YN_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->yn_type[level],nb,0,
			 grid,1,S->yn_halo_type[level],YN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with YN neighbour.");
    } 
  }

  /* Z AXIS */
  if (mpi_iz % 2 == 0) { 
    if (0<=(nb=ZP_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->zp_type[level],nb,0,
			 grid,1,S->zp_halo_type[level],ZP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZP neighbour.");
    } 
  } else { 
    if (0<=(nb=ZN_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->zn_type[level],nb,0,
			 grid,1,S->zn_halo_type[level],ZN_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZN neighbour.");
    } 
  }

  if (mpi_iz % 2 == 1) { 
    if (0<=(nb=ZP_NEIGHBOUR) && S->active[nb]) {
      MPIDO(MPI_Sendrecv(grid,1,S->zp_type[level],nb,0,
			 grid,1,S->zp_halo_type[level],ZP_NEIGHBOUR,0,
			 ALL_ACTIVE_PROCS,MPI_STATUS_IGNORE),
	    "Failed to swap buffers with ZP neighbour.");
    }
  } else { 
    if (0<=(nb=ZN_NEIGHBOUR) && S->active[nb]) {
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

/*=========================================*/
/* The singular problem: is no point of the device space linked to a value  */
/* beyond it, in any subdomain? The links are the nonzero weights of the     */
/* neighbours outside the device space (isolated=0, tissue beyond).          */
static int is_singular(STR *S) {
  Space s=S->s;
  list_t i;
  int nb;
  long nlinks=0;
  for (i=1;i<=S->numPoints;i++) {
    AnisoPoint *P=&(S->points[i]);
    for (nb=1;nb<maxnb;nb++) {
      long x=(long)P->x+ndx[nb], y=(long)P->y+ndy[nb], z=(long)P->z+ndz[nb];
      if (P->K[nb]!=0 && !mg_in_space(S,x,y,z)) nlinks++;
    }
  }
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,&nlinks,1,MPI_LONG,MPI_SUM,S->comm),"Couldn't sum up the boundary links.");
#endif
  return nlinks==0;
}
/*=========================================*/

/* Mean of a level-0 grid over the points of the device space, in all subdomains */
static real points_mean(real *grid,STR *S) {
  double sum[2]={0.0,0.0};
  list_t i;
  for (i=1;i<=S->nPoints[0];i++) sum[0]+=grid[S->pointIndex[0][i]];
  sum[1]=S->nPoints[0];
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,sum,2,MPI_DOUBLE,MPI_SUM,S->comm),"Couldn't sum up a mean.");
#endif
  return (sum[1]>0)?(sum[0]/sum[1]):0.0;
}
/*=========================================*/

/* Mean of layer v over the points of the device space, in all subdomains */
static real layer_mean(int v,STR *S) {
  double sum[2]={0.0,0.0};
  list_t i;
  for (i=1;i<=S->numPoints;i++) sum[0]+=New[ind(S->points[i].x,S->points[i].y,S->points[i].z,v)];
  sum[1]=S->numPoints;
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,sum,2,MPI_DOUBLE,MPI_SUM,S->comm),"Couldn't sum up a mean.");
#endif
  return (sum[1]>0)?(sum[0]/sum[1]):0.0;
}
/*=========================================*/

/* The singular problem with pin=1: add a constant to the solution in layer v1 */
/* so that its value at the pinning point is upin.                              */
static void apply_pin(STR *S) {
  int v1=S->v1;
  double value=0.0;
  list_t i;
  AnisoPoint *P;
  if (S->pinpoint) {
    P=&(S->points[S->pinpoint]);
    value=New[ind(P->x,P->y,P->z,v1)];
  }
#if MPI
  MPIDO(MPI_Allreduce(MPI_IN_PLACE,&value,1,MPI_DOUBLE,MPI_SUM,S->comm),"Couldn't get the value at the pinning point.");
#endif
  value=S->upin-value;
  for (i=1;i<=S->numPoints;i++) {
    P=&(S->points[i]);
    New[ind(P->x,P->y,P->z,v1)]+=value;
  }
}
/*=========================================*/
