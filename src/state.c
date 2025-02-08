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

#if MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "system.h"
#include "beatbox.h"
#include "k_.h"
#include "bikt.h"
#include "decomp.h"
#include "error.h"
#include "geometry.h"
#define OWN
#include "state.h"
#undef OWN
#include "device.h"
#include "qpp.h"
#include "mpi_io_choice.h"

#define DB DEBUG

/* macro for accepting STR pars eponymous to global vars, unlike local ones in qpp.h */
#define ACCEPTLG(b,c,d,e) if (!acceptl(#b"=",&(S->b),c,d,e,w)) return(0); b=S->b

extern int Verbose;
static int already_created=0;

/* This is for ACCEPT macros and compatibility with k_comp.h contents. */
typedef struct {
  /* Geometry file */
  NFILE(geometry);
  /* /\* Geometry options *\/ */
  int normaliseVectors;
  int checkVectors;
  int anisotropy;
  int padding;
  /* Output normalized geometry */
  NFILE(out);
  /* Debugging */
  NFILE(debug);		/* if and where to print debug messages */
  int debugWriter;	/* this thread is entrusted with common debug messages */
  /* Global bounds: "allow sign" to suppress compiler warnings */
  ssize_t xmax, ymax, zmax, vmax;
  /* Stuff for the k-pgm */
  #include "k_code.h"			/* compiled k-pgm */
} STR;

/* This function is performed once just when reading from parameter file */
int state (char *w)
{
  INT x,y,z;		/* local k-vars accessible to the k-pgm */
  real *u;		/*                           [vmax], .. */
  real *geom;		/*                      [geom_vmax], .. */
  int base;		/* index from which shifts are measured. */
  int iv;		/* loop iterator */
  char name[80];	/* buffer for names of local k-vars */
  int hasu;		/* flag that k-pgm assigns u[] */
  p_tb loctb = tb_new(); /* local symbol table */
  STR *S= (STR *) Calloc(1,sizeof(STR));
  
  if (!S) ABORT("not enough memory");
	
  if (already_created) EXPECTED_ERROR("state must be defined only once");

  ACCEPTF(debug,"wt","");
  DEVICE_CONST(int, debugWriter);
  /* DEVICE_CONST(FILE *, debug); */
  
  ACCEPTF(geometry,"r","");
  GEOMETRY_FILE = (S->geometry!=NULL);
  GEOMETRY_PGM = (find_key("pgm=",w)!=NULL);
  ACCEPTI(anisotropy,0,0,1);
  ANISOTROPY_ON = S->anisotropy;
  GEOMETRY_ON = GEOMETRY_FILE || GEOMETRY_PGM || ANISOTROPY_ON;
  geom_vmax = ANISOTROPY_ON? 4: GEOMETRY_ON? 1: 0; /* STATUS + 3 VECTOR COMPONENTS */
	
  if (GEOMETRY_FILE && ANISOTROPY_ON) {
    ACCEPTI(normaliseVectors,0,0,1);
    ACCEPTI(checkVectors,1,0,1);
  }

  if (stateDimensionsExist(w)) {
    ACCEPTLG(xmax,LNONE,1L,LNONE);
    ACCEPTLG(ymax,LNONE,1L,LNONE);
    ACCEPTLG(zmax,   1L,1L,LNONE);
  } else if (GEOMETRY_FILE) {
    ACCEPTI(padding,1,0,INONE);

    if (Verbose) MESSAGE("\n/* Finding minimal enclosing box for geometry from %s ... */\n",
			 S->geometryname);
    getMeshDimensions(S->geometry,&xmax,&ymax,&zmax,S->padding);
    S->xmax = xmax;
    S->ymax = ymax;
    S->zmax = zmax;
    if (Verbose) MESSAGE("\n/* ... %ld x %ld x %ld */\n",xmax,ymax,zmax);
  } else {
    EXPECTED_ERROR("state dimensions are not specified, neither explicitly nor via geometry");
  }
  ACCEPTLG(vmax,2L,2L,LNONE);
  if (Verbose) MESSAGE("\n/* Grid size (%d x %d x %d x %d) */\n",xmax,ymax,zmax,vmax);
  
  /* The geometry may be defined or restricted by a k-code */
  if (GEOMETRY_PGM) {
    int iv;			/* loop iterator */
    Space fake_space;		/* for           */
    Device fake_dev;		/*   the         */
    Device *dev=&fake_dev;	/*   benefit     */
    fake_space.nowhere=0;	/*   of          */
    fake_dev.s=fake_space;	/*   k_comp.h    */
    k_on();						CHK(NULL);
    memcpy(loctb,deftb,sizeof(*deftb));
    tb_insert_int (loctb,"x",&x);			CHK("x");
    tb_insert_int (loctb,"y",&y);			CHK("y");
    tb_insert_int (loctb,"z",&z);			CHK("z");

    CALLOC(u,vmax,sizeof(real));
    for (iv=0;iv<vmax;iv++) {
      snprintf(name,80,"u%d",iv);
      tb_insert_real(loctb,name,&(u[iv]));		CHK(name);
    } /* for iv<vmax */
    
    CALLOC(geom,geom_vmax,sizeof(real));
    for (iv=0;iv<geom_vmax;iv++) {
      snprintf(name,maxname,"geom%d",iv);
      tb_insert_real(loctb,name,&(geom[iv]));  		CHK(name);
    } /* for iv<geom_vmax */
    
    #include "k_comp.h"

    /* geom0 must be assigned */
    if (!used(S->data,S->ncode,&(geom[0])))
      EXPECTED_ERROR("Block pgm={} specified but variable 'geom0' never assigned in it!");

    /* u0..u[nv] may be assigned, need to know if they were */
    hasu=0;
    for (iv=0;iv<vmax;iv++) {
      if (used(S->data,S->ncode,&(u[iv]))) {
	hasu=1;
	break;
      }
    } /* for iv<vmax */
  } /* if GEOM_PGM */

  if (GEOMETRY_FILE) {
    if (mpi_rank==0) {
      ACCEPTF(out,"w","");
    } else {
      S->out=NULL;
      S->outname[0]='\0';
    }
  }
 
  /* If any dimension is exactly 2, complain.
   * Since it's not a single point, we'll need diffusion,
   * but we don't have enough points for boundaries.
   */
  if (xmax==2||ymax==2||zmax==2) 
    EXPECTED_ERROR("One or more dimensions (xmax,ymax,zmax) is equal to 2. Dimensions must be either 1 or >=3.");
  
  /* If the dimension is in use, flag it so we can 
     treat the first and last points as boundaries.*/	
  TRI = (zmax >=3) ? 1 : 0;
  TWO = (ymax >=3) ? 1 : 0;
  ONE = (xmax >=3) ? 1 : 0;
  dim = ONE + TWO + TRI;
  /* safeone=&ONE; */
  /* DEBUG("ONE=%d=%d TWO=%d TRI=%d dim=%d\n",ONE,*safeone,TWO,TRI,dim); */

  if(dim==1 && !ONE)		EXPECTED_ERROR("1-Dimensional simulations must be defined on the x axis.");
  if(dim==1 && GEOMETRY_ON)	EXPECTED_ERROR("Geometry requires a 2- or 3-dimensional simulation medium.");
  if(dim==2 && !(ONE && TWO))	EXPECTED_ERROR("2-Dimensional simulations must be defined on the x and y axes.");

  /* Global boxes */
  BX0=0; BX1=xmax; BY0=0; BY1=ymax; BZ0=0; BZ1=zmax;
  NX0=BX0+ONE; NX1=BX1-ONE; NY0=BY0+TWO; NY1=BY1-TWO; NZ0=BZ0+TRI; NZ1=BZ1-TRI;
  

/**************************************************************************
 ************************** MPI ONLY **************************************
 **************************************************************************/
#if MPI
  /*  Get total number of processes. */
  MPIDO(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size),"Could not get number of MPI processes.");

  /* Geometry skeleton is required for decomp_allocaeSubdomains */
  /* and may be useful for decompose in the future */
  if (GEOMETRY_ON) make_gpoints(S->geometry);

  /* Compute domain decomposition.                 */
  /* Also, define bounds of subdomains along axes. */
  num_subdoms = decompose(mpi_size, &mpi_nx, &mpi_ny, &mpi_nz);
  if (num_subdoms<=0) EXPECTED_ERROR("Could not decompose the domain\n");
  if (Verbose) MESSAGE("\n/* Domain decomposed as (%03d,%03d,%03d). */",mpi_nx,mpi_ny,mpi_nz);
  num_dangling_procs = mpi_size - num_subdoms;

  /* Allocate subdomains to processes */
  num_active_procs = decomp_allocateSubdomains();
  if (num_active_procs<=0) 
    EXPECTED_ERROR("could not allocate subdomains per processes: num_active_procs=%d\n",num_active_procs);
  if (Verbose) MESSAGE("\n/* %d processes will be active, %d will be idle. */",num_active_procs,mpi_size-num_active_procs);
  num_empty_subdoms = num_subdoms-num_active_procs;
  if (Verbose) MESSAGE("\n/* %d subdomains are empty. */",num_empty_subdoms);

  /*  Define communicator over active processes. */
  I_AM_IDLE = decomp_globalCommunicator (mpi_rank, &ALL_ACTIVE_PROCS, &mpi_rank);

  /*  If we're idle, we're done. */
  if (I_AM_IDLE) {
    if (Verbose) URGENT_MESSAGE("\n/* Process %d will be idle */",mpi_rank);
    FREE(S);
    return (1);
  }
	
  decomp_getSuperindices(mpi_rank, &mpi_ix, &mpi_iy, &mpi_iz);
  decomp_getSubdomainDimensions();
  decomp_getSubdomainSize();
  decomp_defineHaloTypes();

  /* Local boxes as determined by decomp_getSubdomainDimensions */
  nx0=local_xmin; nx1=local_xmax;
  ny0=local_ymin; ny1=local_ymax;
  nz0=local_zmin; nz1=local_zmax;
  bx0=nx0-ONE; bx1=nx1+ONE;
  by0=ny0-TWO; by1=ny1+TWO;
  bz0=nz0-TRI; bz1=nz1+TRI;

  /* Locate processes of neighbouring subdomains. */
  /* These will return NORANK for empty neighbours. */
  xn_neighbour=getProcessRank(mpi_ix-1,mpi_iy  ,mpi_iz  ); /* neighbour to LEFT	*/
  xp_neighbour=getProcessRank(mpi_ix+1,mpi_iy  ,mpi_iz  ); /* neighbour to RIGHT	*/
  yn_neighbour=getProcessRank(mpi_ix  ,mpi_iy-1,mpi_iz  ); /* neighbour ABOVE	*/
  yp_neighbour=getProcessRank(mpi_ix  ,mpi_iy+1,mpi_iz  ); /* neighbour BELOW	*/
  zn_neighbour=getProcessRank(mpi_ix  ,mpi_iy  ,mpi_iz-1); /* neighbour in FRONT	*/
  zp_neighbour=getProcessRank(mpi_ix  ,mpi_iy  ,mpi_iz+1); /* neighbour BEHIND	*/

  /* Report subdomain details of this processes. */
  if (Verbose) {
    FFLUSH(stdout);
    MPI_Barrier(ALL_ACTIVE_PROCS);
    URGENT_MESSAGE(
      "\n/* Process %d=(%d,%d,%d) domain size (%dx%dx%d) limits [%d:%d]x[%d:%d]x[%d:%d] neighbours %2d:%2d, %2d:%2d, %2d:%2d */",
      mpi_rank,mpi_ix,mpi_iy,mpi_iz,xlen,ylen,zlen,
      local_xmin, local_xmax-1, local_ymin, local_ymax-1, local_zmin, local_zmax-1,
      xn_neighbour,xp_neighbour,yn_neighbour,yp_neighbour,zn_neighbour,zp_neighbour
    );
    MPI_Barrier(ALL_ACTIVE_PROCS);
  }

#else /* not MPI */
/**************************************************************************
 ************************** SEQUENTIAL ONLY *******************************
 **************************************************************************/
  /*  No partitioning. */
  mpi_nx=1;
  mpi_ny=1;
  mpi_nz=1;
  /*  Size of the medium, including boundaries. */
  xlen = xmax;
  ylen = ymax;
  zlen = zmax;
  /* Local boxes same as global boxes */
  nx0=NX0; nx1=NX1;
  ny0=NY0; ny1=NY1;
  nz0=NZ0; nz1=NZ1;
  bx0=BX0; bx1=BX1;
  by0=BY0; by1=BY1;
  bz0=BZ0; bz1=BZ1;  
#endif /* not MPI */

/**************************************************************************
 ************************** COMMON ****************************************
 **************************************************************************/

  /* This is the actual declaration of memory (used by this subdomain if MPI) */
  CALLOC(New,xlen*ylen*zlen*vmax,sizeof(real)); DEBUG("qlen=%ld*%ld*%ld vmax=%ld New=%p\n",xlen,ylen,zlen,vmax,New);
	
  /*  Get constants for index computation. */
  vmax_zmax = vmax * zlen;
  vmax_zmax_ymax = vmax * zlen * ylen;
	
  /* If non-box geometry (for either reason) */
  if (geom_vmax>0) {
    /* Create the geometry grid */
    CALLOC(Geom,xlen*ylen*zlen*geom_vmax,sizeof(real));
    /*  Get constants for geometry index computation. */
    geom_vmax_zmax= geom_vmax * zlen;
    geom_vmax_zmax_ymax= geom_vmax * zlen * ylen;
    gx=geom_ind(1,0,0,0);
    gy=geom_ind(0,1,0,0);
    gz=geom_ind(0,0,1,0);
      
    /* Fill geom grid, with fake data if no geom file provided */
    if (!populateMesh(S->geometry,S->geometryname,S->normaliseVectors,S->checkVectors,S->out,S->outname)) return 0;
    if (geometry) {fclose(geometry); geometry=S->geometry=NULL;}
    if (S->out) {fclose(S->out); S->out=NULL;}
  } /* if geom_vmax */

  /* If k-code amendments are requested */
  if (GEOMETRY_PGM) {
    if (Verbose) MESSAGE("/* Geometry amended by k-pgm */\n");
    DEVICE_CONST(int,ncode);
    DEVICE_ARRAY(pp_fn,code);
    DEVICE_ARRAY(p_real,data);
    assert(geom_vmax>0);
    /* */				DEBUG("program geometry in [%ld:%ld]x[%ld:%ld]x[%ld:%ld]\n",
					      (long)local_xmin,local_xmax-1,
					      (long)local_ymin,local_ymax-1,
					      (long)local_zmin,local_zmax-1);
    /* brutto local */
#if MPI
    size_t x0 = local_xmin - ONE;
    size_t x1 = local_xmax + ONE;
    size_t y0 = local_ymin - TWO;
    size_t y1 = local_ymax + TWO;
    size_t z0 = local_zmin - TRI;
    size_t z1 = local_zmax + TRI;
#else
    size_t x0 = 0;
    size_t x1 = xmax;
    size_t y0 = 0;
    size_t y1 = ymax;
    size_t z0 = 0;
    size_t z1 = zmax;
#endif
    /* We need to designate tissue points in the box which is	*/
    /* brutto box minus outer halos,	   			*/
    /* or equivalently, netto box plus interface halos,		*/
    /* or by set algebra, NETTO*brutto.				*/
    size_t xlo=max(SPACE_DEFAULT_X0,x0), xhi=min(SPACE_DEFAULT_X1+1,x1);
    size_t ylo=max(SPACE_DEFAULT_Y0,y0), yhi=min(SPACE_DEFAULT_Y1+1,y1);
    size_t zlo=max(SPACE_DEFAULT_Z0,z0), zhi=min(SPACE_DEFAULT_Z1+1,z1);

    for (x=xlo; x<xhi; x++) {
      for (y=ylo; y<yhi; y++) {
	for (z=zlo; z<zhi; z++) {
	  if (GEOM_VOID != Geom[geom_ind(x,y,z,GEOM_STATUS)]) {
            #include "k_exec.h"
	    memcpy(Geom+geom_ind(x,y,z,0),geom,geom_vmax*sizeof(real));
	    if (hasu) memcpy(New+ind(x,y,z,0),u,vmax*sizeof(real));
	    /* if (x==9 && y==6) */
	      /* DEBUG("xyz=(%ld,%ld,%ld) geom=%g &geom=%p\n", */
	      /* 	    x,y,z,Geom[geom_ind(x,y,z,0)],Geom+geom_ind(x,y,z,0)); */
	  } /* not GEOM_VOID */
	} /* for z */
      } /* for y */
    } /* for x */
  } /* GEOMETRY_PGM */

  /* Make the list of 'tissue' points on which we can operate.	*/
  /* Here we are restricted to the netto box:			*/
  /* interface halo may be inquired, but not operated on	*/
  if (GEOMETRY_ON) {
    assert(geom_vmax>0);
    numTissuePoints=0;
    CALLOC(TissuePoints,
	   (local_xmax-local_xmin)*(local_ymax-local_ymin)*(local_zmax-local_zmin),
	   sizeof(devicePoint));
    /* DEBUG("TissuePoints=%p\n",TissuePoints); */
    /* DEBUG("(xyz):[%ld:%ld[*[%ld:%ld[*[%ld:%ld[\n", */
    /* 	  (long)local_xmin, (long)local_xmax, */
    /* 	  (long)local_ymin, (long)local_ymax, */
    /* 	  (long)local_zmin, (long)local_zmax); */
    for (z=local_zmin; z<local_zmax; z++) {
      for (y=local_ymin; y<local_ymax; y++) {
	for (x=local_xmin; x<local_xmax; x++) {
	  if (GEOM_VOID != Geom[geom_ind(x,y,z,GEOM_STATUS)]) {
	    TissuePoints[numTissuePoints].x=x;
	    TissuePoints[numTissuePoints].y=y;
	    TissuePoints[numTissuePoints].z=z;
	    TissuePoints[numTissuePoints].u=New+ind(x,y,z,0);
	    TissuePoints[numTissuePoints].X=NULL;
	    TissuePoints[numTissuePoints].Y=NULL;
	    /* DEBUG("%ld:(%ld,%ld,%ld)\n",numTissuePoints,x,y,z); */
	    numTissuePoints++;
	  } /* not GEOM_VOID */
	} /* for z */
      } /* for y */
    } /* for x */
    REALLOC(TissuePoints,numTissuePoints*sizeof(devicePoint));
    /* DEBUG("TissuePoints=%p\n",TissuePoints); */
    MESSAGE("/* Geometry is ON, %ld points listed */\n", numTissuePoints);
  } /* GEOMETRY_ON */
	
 /* Shift increments: should be computed after ind is fully defined */
#if MPI
  /*  Computed from local points to avoid upsetting ind() */
  base = ind(local_xmin,   local_ymin,   local_zmin,   0);
  DX   = ind(local_xmin+1, local_ymin,   local_zmin,   0) - base;
  DY   = ind(local_xmin,   local_ymin+1, local_zmin,   0) - base;
  DZ   = ind(local_xmin,   local_ymin,   local_zmin+1, 0) - base;
  DV   = ind(local_xmin,   local_ymin,   local_zmin,   1) - base;
#else /*  Sequential */
  base = ind(0,0,0,0);
  DX   = ind(1,0,0,0) - base;
  DY   = ind(0,1,0,0) - base;
  DZ   = ind(0,0,1,0) - base;
  DV   = ind(0,0,0,1) - base;
#endif

  DEBUG("\n"
	"BRUTTO:[%ld:%ld)x[%ld:%ld)x[%ld:%ld)\n"
	"NETTO:[%ld:%ld)x[%ld:%ld)x[%ld:%ld)\n"
	"brutto:[%ld:%ld)x[%ld:%ld)x[%ld:%ld)\n"
	"netto:[%ld:%ld)x[%ld:%ld)x[%ld:%ld)\n",
	BX0,BX1,BY0,BY1,BZ0,BZ1,
	NX0,NX1,NY0,NY1,NZ0,NZ1,
	bx0,bx1,by0,by1,bz0,bz1,
	nx0,nx1,ny0,ny1,nz0,nz1);
  
  FREE(S);
  return SUCCESS;
} /*  state() */

void state_free(void)
{
  FREE(Geom);
  FREE(New);
#if MPI
  FREE(occupy_grid);
  FREE(occupy_rank);
#endif
} /* state_free () */


/**************************************************************************
 ************************** MPI FUNCTIONS *********************************
 **************************************************************************/
#if MPI

/* This is a slow and careful version intended for use */
/* at start time rather than run time                  */
int getProcessRank (int ix,int iy,int iz) {
  /* These checks here are unnecessary if done in the calling function */
  if (ix<0) return NORANK; if (ix>=mpi_nx) return NORANK;
  if (iy<0) return NORANK; if (iy>=mpi_ny) return NORANK;
  if (iz<0) return NORANK; if (iz>=mpi_nz) return NORANK;
  /* The line below is enough if checks done in the calling function. */
  /* If the subdomain is empty it still will return NORANK. */
  return ProcessRank[sind(ix,iy,iz)];
}
/*------------*/

/* This version exploits axes partitions rather than going through all subdomains. */
/* To do: optimize this using details of the partitioning algorithm */
int getRankContainingPoint(int x,int y,int z) {
  int i, ix, iy, iz;
  ix=iy=iz=-1;
  for (i=0;i<mpi_nx;i++) {if (all_xmin[i]<=x && x<all_xmax[i]) {ix=i; break;};} if (ix<0) return NORANK;
  for (i=0;i<mpi_ny;i++) {if (all_ymin[i]<=y && y<all_ymax[i]) {iy=i; break;};} if (iy<0) return NORANK;
  for (i=0;i<mpi_nz;i++) {if (all_zmin[i]<=z && z<all_zmax[i]) {iz=i; break;};} if (iz<0) return NORANK;
  /* No need for extra checks now. */
  /* If the subdomain is empty it still will return NORANK. */
  return ProcessRank[sind(ix,iy,iz)];
}
/*------------*/

#endif /* MPI */
