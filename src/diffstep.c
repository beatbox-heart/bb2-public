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
 * Diffusion timestep device for operator-splitting schemes.
 * Can supersede diff (just set ht=0). 
 */

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "qpp.h"

/* Point descriptions for different cases */

/* Anisotropic case, 3D */
typedef struct {
  real K000,				/* Neibours weights */
    Kp00, Km00,
    K0p0, K0m0,
    K00p, K00m,
    Kpp0, Kpm0, Kmp0, Kmm0,
    Kp0p, Kp0m, Km0p, Km0m,
    K0pp, K0pm, K0mp, K0mm;
  int x,y,z;				/* Coordinates */
} AnisoPoint;      

/* Isotropic 3D */
typedef struct {
  real K000,				/* Neibours weights */
    Kp00, Km00,
    K0p0, K0m0,
    K00p, K00m;
  int x,y,z;				/* Coordinates */
} Iso3Point;      

/* Isotropic 2D */
typedef struct {
  real K000,				/* Neibours weights */
    Kp00, Km00,
    K0p0, K0m0;
  int x,y,z;				/* Coordinates */
} Iso2Point;      

/* Isotropic 1D */
typedef struct {
  real K000,				/* Neibours weights */
    Kp00, Km00;
  int x,y,z;				/* Coordinates */
} Iso1Point;      

typedef struct {                
  real D;				/* Isotropic diffusion coefficient. */
  real Dpar;				/* Parallel diffusion coefficient. ANISO only */
  real Dtrans;				/* Transverse diffusion coefficient. ANISO only. */
  real hx;				/* Space step */
  real ht;				/* Time step */
  int v0;				/* layer for the diffusive field */
  int v1;				/* auxiliary layer */
  
  /* ANISOTROPY ONLY */
  int numPoints;			/* Number of tissue apoints in subdomain. */
  AnisoPoint *apoints;			/* List of tissue apoints. */
  Iso3Point *i3points;			/* List of tissue apoints. */
  Iso2Point *i2points;			/* List of tissue apoints. */
  Iso1Point *i1points;			/* List of tissue apoints. */
} STR;

  /* If nonzero timestep set, do the Euler step for diffusion equation */
#define IF_HT_DO_STEP(pointschoice) \
  if (ht) { \
    for(ind=0;ind<numPoints;ind++){ \
      P=S->pointschoice[ind];\
      x = P.x; \
      y = P.y; \
      z = P.z; \
      u=New+ind(x,y,z,v0); \
      u[0] += ht*u[DV*(v1-v0)]; \
    } \
  }

/****************/
RUN_HEAD(diffstep)
{
  int x, y, z, ind;
  real sum, value;
  real *u;
  DEVICE_CONST(real,hx)
  DEVICE_CONST(real,ht)
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(int,numPoints)

  if(!ANISOTROPY_ON){

    /*********************/
    /* ISOTROPIC VERSION */
    DEVICE_CONST(real,D)
    real gam=D/(hx*hx);
    if (dim==3) {
      Iso3Point P;
      for(ind=0;ind<numPoints;ind++){
	P=S->i3points[ind];
	x = P.x;
	y = P.y;
	z = P.z;
	/* Centre point */
	u=New+ind(x,y,z,v0);
	/* Sum by neighbours */
	sum = u[0]*P.K000
	  + u[DX]*P.Kp00 + u[-DX]*P.Km00 
	  + u[DY]*P.K0p0 + u[-DY]*P.K0m0
	  + u[DZ]*P.K00p + u[-DZ]*P.K00m;
	/* Laplacian value */
	u[DV*(v1-v0)] = gam*sum;
      } /* for ind */
      IF_HT_DO_STEP(i3points);
    } /* shortened version for lesser dimensions */
    else if (dim==2) {
      Iso2Point P;
      for(ind=0;ind<numPoints;ind++){
	P=S->i2points[ind];
	x = P.x;
	y = P.y;
	z = P.z;
	/* Centre point */
	u=New+ind(x,y,z,v0);
	/* Sum by neighbours */
	sum = u[0]*P.K000
	  + u[DX]*P.Kp00 + u[-DX]*P.Km00 
	  + u[DY]*P.K0p0 + u[-DY]*P.K0m0;
	/* Laplacian value */
	u[DV*(v1-v0)] = gam*sum;
      } /* for ind */
      IF_HT_DO_STEP(i2points);
    } else if (dim==1) {
      Iso1Point P;
      for(ind=0;ind<numPoints;ind++){
	P=S->i1points[ind];
	x = P.x;
	y = P.y;
	z = P.z;
	/* Centre point */
	u=New+ind(x,y,z,v0);
	/* Sum by neighbours */
	sum = u[0]*P.K000
	  + u[DX]*P.Kp00 + u[-DX]*P.Km00;
	/* Laplacian value */
	u[DV*(v1-v0)] = gam*sum;
      } /* for ind */
      IF_HT_DO_STEP(i1points);
    } /* dimension cases */

  } else {
    
    /***********************/
    /* ANISOTROPIC VERSION */
    DEVICE_ARRAY(AnisoPoint, apoints)
    DEVICE_CONST(int, numPoints)
    int ind;
    real sum;
    /* The current point's structure */
    AnisoPoint P;
    
    if (dim==3) {
      for(ind=0;ind<numPoints;ind++){
	P=S->apoints[ind];
	x = P.x;
	y = P.y;
	z = P.z;
	/* Centre point */
	u=New+ind(x,y,z,v0);
	/* Sum by neighbours */
	sum = u[0]*P.K000
	  + u[DX]*P.Kp00 + u[-DX]*P.Km00 
	  + u[DY]*P.K0p0 + u[-DY]*P.K0m0
	  + u[DZ]*P.K00p + u[-DZ]*P.K00m
	  + u[DX+DY]*P.Kpp0 + u[DX-DY]*P.Kpm0 + u[-DX+DY]*P.Kmp0 + u[-DX-DY]*P.Kmm0
	  + u[DX+DZ]*P.Kp0p + u[DX-DZ]*P.Kp0m + u[-DX+DZ]*P.Km0p + u[-DX-DZ]*P.Km0m
	  + u[DY+DZ]*P.K0pp + u[DY-DZ]*P.K0pm + u[-DY+DZ]*P.K0mp + u[-DY-DZ]*P.K0mm;
	/* Laplacian value */
	u[DV*(v1-v0)] = sum / (4*hx*hx);
      } /* for ind */
    } /* shortened version for lesser dimensions */
    else if (dim==2) {
      for(ind=0;ind<numPoints;ind++){
	P=S->apoints[ind];
	x = P.x;
	y = P.y;
	z = P.z;
	u=New+ind(x,y,z,v0);
	sum = u[0]*P.K000
	  + u[DX]*P.Kp00 + u[-DX]*P.Km00 
	  + u[DY]*P.K0p0 + u[-DY]*P.K0m0
	  + u[DX+DY]*P.Kpp0 + u[DX-DY]*P.Kpm0 + u[-DX+DY]*P.Kmp0 + u[-DX-DY]*P.Kmm0;
	u[DV*(v1-v0)] = sum / (4*hx*hx);
      } 
    } else if (dim==1) {
      for(ind=0;ind<numPoints;ind++){
	P=S->apoints[ind];
	x = P.x;
	y = P.y;
	z = P.z;
	u=New+ind(x,y,z,v0);
	sum = u[0]*P.K000
	  + u[DX]*P.Kp00 + u[-DX]*P.Km00;
	u[DV*(v1-v0)] = sum / (4*hx*hx);
      }
    } /* dimensions cases */
    IF_HT_DO_STEP(apoints);
  } /* if !ANISOTROPY */
}
RUN_TAIL(diffstep)

/****************/
DESTROY_HEAD(diffstep)
  if (S->apoints) FREE(S->apoints);
DESTROY_TAIL(diffstep)

/****************/
CREATE_HEAD(diffstep)
{
  DEVICE_REQUIRES_SYNC;
  /* DEVICE_HAS_DEFAULT_SPACE */

  /* counters */
  int x,y,z,i,j,ind; 

  ACCEPTR(hx,RNONE,0.,RNONE);
  ACCEPTR(ht,0.,0.,RNONE);
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 != v0 );

  /* Count tissue points */
  S->numPoints = 0;
  for (x=dev->s.x0;x<=dev->s.x1;x++) {
    for (y=dev->s.y0;y<=dev->s.y1;y++) {
      for (z=dev->s.z0;z<=dev->s.z1;z++) {
	if (isTissue(x,y,z)) S->numPoints++;
      }
    }
  }

  /***************************************************************/
  /* ANISOTROPIC CASE                                           */
  /* We shell treat dim<3 cases as dim=3 cases where the problem */
  /*   simply does not depend on z, or z and y, coords.          */
  /***************************************************************/
  if (ANISOTROPY_ON) {
    /* Array dimensions are 4 to make indices match paper */
    real f[4];				/* Fibre direction vector */
    real D[4][4];			/* Diffusion tensor */
    real dD[4][4];			/* Derivatives of diffusion tensor */
    real c[4];				/* Derivatives of diffusion coefficients */
    real Dscalar;			/* Back-up scalar diffusivity for singular points */
    /* Neighbours' vectors */
    int nx,ny,nz; /* Coords of negative neighbour */
    int px,py,pz; /* Coords of positive neighbour */

    /* Anisotropic diffusion coefficients */
    ACCEPTR(Dpar,RNONE,0.,RNONE);
    ACCEPTR(Dtrans,RNONE,0.,RNONE);

    /* There should be no single diffusion coefficient normally */
    Dscalar=pow(Dpar*Dpar*Dtrans,1.0/3.0);
    if (find_key("D=",w)) { /* Insulate this scalar D from the tensor D defined in the outside block */
      ACCEPTR(D,Dscalar,0.,RNONE);
      Dscalar=D;
      MESSAGE("/* The isotropic diffusion coefficient 'D' "
	      "is mostly unused when anisotropy is active. "
	      "It will apply only at exceptional 'isotropic' points, if any. */\n");
    }	

    /* The current point's structure */
    AnisoPoint P;
    
    /* Allocate array of AnisoPoints */
    CALLOC(S->apoints,S->numPoints,sizeof(AnisoPoint));

    ind = 0;
    for (x=dev->s.x0;x<=dev->s.x1;x++) {
      for (y=dev->s.y0;y<=dev->s.y1;y++) {
	for (z=dev->s.z0;z<=dev->s.z1;z++) {
	  if (isTissue(x,y,z)) {

	    /*-------------------------------*/
	    /* The diff tensor of this point */

	    /* Fibre direction cosines of this point */
	    f[1] =  Geom[geom_ind(x,y,z,GEOM_FIBRE_1)];
	    f[2] =  Geom[geom_ind(x,y,z,GEOM_FIBRE_2)];
	    f[3] =  Geom[geom_ind(x,y,z,GEOM_FIBRE_3)];
#define UNIT_VECTOR_TOLERANCE 0.01 /* defined in geometry.h; should we rather #include that? */
	    if (f[1]*f[1]+f[2]*f[2]+f[3]*f[3]<UNIT_VECTOR_TOLERANCE) { /* this is an 'isotropic point' */
#undef UNIT_VECTOR_TOLERANCE
	      /* set up the scalar */
	      for (i=1;i<=3;i++) {
		for (j=1;j<3;j++) {
		  D[i][j] = (i==j)?Dscalar:0;
		}
	      }
	      for (i=1;i<=3;i++) c[i] = 0; /* simple and ad hoc; rethink it later */
	    } else {
	      /* set up the tensor */
	      for (i=1;i<=3;i++) {
		for (j=1;j<=3;j++) {
		  /* \mathbf{D} = 
		   * D_\perp \mathbf{I} + (D_\parallel - D\perp)
		   * \mathbf{A}\mathbf{A}^T
		   */
		  D[i][j] = ((Dpar-Dtrans)*f[i]*f[j]) + ((i==j)?Dtrans:0);
		}
	      }
	      /* and the first derivatives' vector */
	      for (i=1;i<=3;i++) {
		if(
		   !isTissue( x+(i==1?1:0), y+(i==2?1:0), z+(i==3?1:0) ) ||
		   !isTissue( x-(i==1?1:0), y-(i==2?1:0), z-(i==3?1:0) )
		   ){
		  for (j=1;j<=3;j++) {dD[i][j]=0.0;}
		} else {
		  for (j=1;j<=3;j++) {
		    /* d^{11}_{+..} - d^{11}_{-..}
		     * = (D_\parallel - D_\perp)(a_i a_j)_{+..} -(a_i a_j)_{-..}
		     */
		    nx = x-(i==1?1:0);
		    ny = y-(i==2?1:0);
		    nz = z-(i==3?1:0);
		    px = x+(i==1?1:0);
		    py = y+(i==2?1:0);
		    pz = z+(i==3?1:0); /* spatial derivatives of the tensor */
		    dD[i][j] =
		      (Dpar - Dtrans) * 
		      (
		       (Geom[geom_ind(px,py,pz,(GEOM_FIBRE_1+i-1) )] *
			Geom[geom_ind(px,py,pz,(GEOM_FIBRE_1+j-1) )] )-
		       (Geom[geom_ind(nx,ny,nz,(GEOM_FIBRE_1+i-1) )] *
			Geom[geom_ind(nx,ny,nz,(GEOM_FIBRE_1+j-1) )] )
		       );
		  } /* for j */
		} /* else */
	      } /* for i */
	      
	      for (i=1;i<=3;i++) {
		c[i] = 0;
		for (j=1;j<=3;j++) {
		  c[i] += dD[j][i];
		}
	      }
	    } /* if (f[1]*...<TOLERANCE) else */

	    /*----------------------------------------*/
	    /* Compute the weights of the neighbours. */
	    
	    /* For dim<3 cases, some weights will be calculated but not used. */
	    /* This can be optimized one day, to save memory.                 */

	    P.K000=
	      P.Kp00=P.Km00=
	      P.K0p0=P.K0m0=
	      P.K00p=P.K00m=
	      P.Kpp0=P.Kpm0=P.Kmp0=P.Kmm0=
	      P.Kp0p=P.Kp0m=P.Km0p=P.Km0m=
	      P.K0pp=P.K0pm=P.K0mp=P.K0mm=0.0;
	    
	    /* Transform tensor coefficients into the sum weights. */
	    /* For dim<3, if (isTissue) will automatically ignore some of these lines. */
	    /* NB this is potentially unsafe in case the "non-existent" neighbours contain NaN !! */
	    
	    if (isTissue(x+1,y	 ,z  )) { P.Kp00+=4*D[1][1]; P.K000-=4*D[1][1]; }
	    if (isTissue(x-1,y	 ,z  ))	{ P.Km00+=4*D[1][1]; P.K000-=4*D[1][1]; }
	    if (isTissue(x  ,y+1,z  ))	{ P.K0p0+=4*D[2][2]; P.K000-=4*D[2][2]; }
	    if (isTissue(x  ,y-1,z  ))	{ P.K0m0+=4*D[2][2]; P.K000-=4*D[2][2]; }
	    if (isTissue(x  ,y  ,z+1))	{ P.K00p+=4*D[3][3]; P.K000-=4*D[3][3]; }
	    if (isTissue(x  ,y	 ,z-1))	{ P.K00m+=4*D[3][3]; P.K000-=4*D[3][3]; }
	    if (isTissue(x+1,y  ,z  ) && isTissue(x-1,y  ,z  )) { P.Kp00+=c[1]; P.Km00-=c[1]; }
	    if (isTissue(x  ,y+1,z  ) && isTissue(x  ,y-1,z  )) { P.K0p0+=c[2]; P.K0m0-=c[2]; }
	    if (isTissue(x  ,y  ,z+1) && isTissue(x  ,y  ,z-1)) { P.K00p+=c[3]; P.K00m-=c[3]; }
	    if (isTissue(x+1,y+1,z  ) && isTissue(x+1,y-1,z  )) { P.Kpp0+=2*D[1][2]; P.Kpm0-=2*D[1][2]; }
	    if (isTissue(x-1,y-1,z  ) && isTissue(x-1,y+1,z  )) { P.Kmm0+=2*D[1][2]; P.Kmp0-=2*D[1][2]; }
	    if (isTissue(x+1,y  ,z+1) && isTissue(x+1,y  ,z-1)) { P.Kp0p+=2*D[1][3]; P.Kp0m-=2*D[1][3]; }
	    if (isTissue(x-1,y  ,z-1) && isTissue(x-1,y  ,z+1)) { P.Km0m+=2*D[1][3]; P.Km0p-=2*D[1][3]; }
	    if (isTissue(x  ,y+1,z+1) && isTissue(x  ,y+1,z-1)) { P.K0pp+=2*D[2][3]; P.K0pm-=2*D[2][3]; }
	    if (isTissue(x  ,y-1,z-1) && isTissue(x  ,y-1,z+1)) { P.K0mm+=2*D[2][3]; P.K0mp-=2*D[2][3]; }

	    /* And the coordinates of this point */
	    P.x = x;
	    P.y = y;
	    P.z = z;
	    
	    /******************************/
	    /* Add this point to the list */
	    S->apoints[ind] = P;
	    ind++;
	  } /* isTissue */
	} /* z */
      } /* y */
    } /* x */
    S->i3points=NULL;
    S->i2points=NULL;
    S->i1points=NULL;
  } else {
  /***************************************************************/
  /* ISOTROPIC CASE                                            */
  /* Here different dimensionalities are treated separately      */
  /***************************************************************/

    /* Single diffusion coefficient */
    ACCEPTR(D,RNONE,0.,RNONE);
    
    /* Should not have anisotropic diffusion coefficients */
    if (find_key("Dpar=",w) || find_key("Dtrans=",w)) {
      MESSAGE("The anisotropic diffusion coefficients 'Dpar' and 'Dtrans' "
	      "are unused when anisotropy is inactive. "
	      "The parameter(s) will be ignored.\n");
    }

    /* Make shorthand copies of important parameters */
    DEVICE_CONST(real,Dpar)
    DEVICE_CONST(real,Dtrans)
    DEVICE_CONST(real,hx)

    /***********/
    /* 3D case */
    if (dim==3) {
      /* The current point's structure */
      Iso3Point P; 

      /* Allocate array of points */
      CALLOC(S->i3points,S->numPoints,sizeof(Iso3Point));

      /* Fill in the array */
      ind = 0;
      for (x=dev->s.x0;x<=dev->s.x1;x++) {
	for (y=dev->s.y0;y<=dev->s.y1;y++) {
	  for (z=dev->s.z0;z<=dev->s.z1;z++) {
	    if (isTissue(x,y,z)) {

	    /* Connection weights */
	      P.K000=
		P.Kp00=P.Km00=
		P.K0p0=P.K0m0=
		P.K00p=P.K00m=0.0;
	      
	      if (isTissue(x+1,y	 ,z  )) { P.Kp00++; P.K000--; }
	      if (isTissue(x-1,y	 ,z  )) { P.Km00++; P.K000--; }
	      if (isTissue(x  ,y+1 ,z  )) { P.K0p0++; P.K000--; }
	      if (isTissue(x  ,y-1 ,z  )) { P.K0m0++; P.K000--; }
	      if (isTissue(x  ,y   ,z+1)) { P.K00p++; P.K000--; }
	      if (isTissue(x  ,y   ,z-1)) { P.K00m++; P.K000--; }
	      
	      /* And the coordinates of this point */
	      P.x = x;
	      P.y = y;
	      P.z = z;
	      
	      /******************************/
	      /* Add this point to the list */
	      S->i3points[ind] = P;
	      ind++;
	    } /* isTissue */
	  } /* z */
	} /* y */
      } /* x */
      S->apoints=NULL;
      S->i2points=NULL;
      S->i1points=NULL;
    } else 
      
    /***********/
    /* 2D case */
    if (dim==2) {
      /* The current point's structure */
      Iso2Point P; 
      
      /* Allocate array of points */
      CALLOC(S->i2points,S->numPoints,sizeof(Iso2Point));
      
      /* Fill in the array */
      ind = 0;
      for (x=dev->s.x0;x<=dev->s.x1;x++) {
	for (y=dev->s.y0;y<=dev->s.y1;y++) {
	  for (z=dev->s.z0;z<=dev->s.z1;z++) {
	    if (isTissue(x,y,z)) {
	      
	      /* Connection weights */
	      P.K000=
		P.Kp00=P.Km00=
		P.K0p0=P.K0m0=0.0;
	      
	      if (isTissue(x+1,y	 ,z  )) { P.Kp00++; P.K000--; }
	      if (isTissue(x-1,y	 ,z  )) { P.Km00++; P.K000--; }
	      if (isTissue(x  ,y+1 ,z  )) { P.K0p0++; P.K000--; }
	      if (isTissue(x  ,y-1 ,z  )) { P.K0m0++; P.K000--; }
	      
	      /* And the coordinates of this point */
	      P.x = x;
	      P.y = y;
	      P.z = z;
	      
	      /******************************/
	      /* Add this point to the list */
	      S->i2points[ind] = P;
	      ind++;
	    } /* isTissue */
	  } /* z */
	} /* y */
      } /* x */
      S->apoints=NULL;
      S->i3points=NULL;
      S->i1points=NULL;
    } else
      
    /***********/
    /* 1D case */
    if (dim==1) {
      /* The current point's structure */
      Iso1Point P; 
      
      /* Allocate array of points */
      CALLOC(S->i1points,S->numPoints,sizeof(Iso1Point));
      
      /* Fill in the array */
      ind = 0;
      for (x=dev->s.x0;x<=dev->s.x1;x++) {
	for (y=dev->s.y0;y<=dev->s.y1;y++) {
	  for (z=dev->s.z0;z<=dev->s.z1;z++) {
	    if (isTissue(x,y,z)) {
	      
	      /* Connection weights */
	      P.K000=
		P.Kp00=P.Km00=0.0;
	      
	      if (isTissue(x+1,y	 ,z  )) { P.Kp00++; P.K000--; }
	      if (isTissue(x-1,y	 ,z  )) { P.Km00++; P.K000--; }
	      
	      /* And the coordinates of this point */
	      P.x = x;
	      P.y = y;
	      P.z = z;
	      
	      /******************************/
	      /* Add this point to the list */
	      S->i1points[ind] = P;
	      ind++;
	    } /* isTissue */
	  } /* z */
	} /* y */
      } /* x */
      S->apoints=NULL;
      S->i3points=NULL;
      S->i2points=NULL;
    } else {
      EXPECTED_ERROR("diffstep device requested for 0D case, which does not make sense\n");
    } /* dimension cases */

  } /* anisotropy cases */
}
CREATE_TAIL(diffstep,1)
