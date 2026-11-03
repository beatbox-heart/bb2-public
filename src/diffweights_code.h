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

/* diffweights_code.h: computation of the connection weights of the "universal"
 * diffusion device diff, shared with other devices (e.g. elliptic). Include
 * after the STR typedef, which must include diffweights_str.h. Provides:
 *   diff_weights_accept(): accept the parameters of the weights,
 *   diff_connect(): compute the weights at the listed points of a space. */

#ifndef DIFFWEIGHTS_CODE_H
#define DIFFWEIGHTS_CODE_H

/* (Re-)Compute the weights */
static int diff_connect (STR *S,Space *s);

/* Averaging of diffusivities of neighbouring nodes for conservative schemes */
static real arithmean(real a, real b)
{
  return 0.5*(a+b);
}
static real harmomean(real a, real b)
{
  real ab=a*b;
  return (ab==0.0)?0.0:(2.0*ab/(a+b));
}
#define MEAN(a,b) (harmonic?(harmomean((a),(b))):(arithmean((a),(b))))

#define DB(...) DEBUG(__VA_ARGS__)
/* #define DB(...)  */

/* Accepting real variable parameter by new name 'b' or its old alias 'a' */
/* TODO: make the concept of aliases to parameters universal for all devices that might benefit */
#define ACCEPTRR(a,b,c,d,e) { \
  int founda=find_key(#a"=",w)?1:0; \
  int foundb=find_key(#b"=",w)?1:0; \
  switch (2*founda+foundb) { \
  case 2*0+0: \
  case 2*0+1: \
    ACCEPTRL(b,c,d,e);	\
    break; \
  case 2*1+0:  \
    ACCEPTRA(a,b,c,d,e); \
    break; \
  case 2*1+1:  \
    EXPECTED_ERROR("\nBoth %s and %s specified, mean the same thing, make up your mind!\n",#a,#b); \
    break; \
  } /* switch */ \
} real b=S->b;

/* Check if given expression depends on local vars, i.e. is space-dependent;	*/
/* hask: yes; hasu: in particular u*; hasg: in particular geom*.		*/
/* An expression of global variables only, e.g. D=~2*Dglob, is not		*/
/* space-dependent, and is evaluated once per (re)weighting.			*/
#define PAR_USE_CHECK(name) { \
  if (S->name##exe) { \
    char varname[maxname]; \
    int iv; \
    hasexe=1; \
    if (k_expr_depends (S->name##src,"x") || \
        k_expr_depends (S->name##src,"y") || \
        k_expr_depends (S->name##src,"z")) hask=1; \
    for (iv=0;iv<vmax;iv++) { \
      snprintf(varname,maxname,"u%d",iv); \
      if (k_expr_depends (S->name##src,varname)) hasu=hask=1; \
    } /* for iv */ \
    for (iv=0;iv<geom_vmax;iv++) { \
      snprintf(varname,maxname,"geom%d",iv); \
      if (k_expr_depends (S->name##src,varname)) hasg=hask=1; \
    } /* for iv */ \
  } else if (S->name##ptr) { \
    for (iv=0;iv<vmax;iv++)\
      if (S->name##ptr == &(S->u[iv]))\
	hasu=hask=1;\
    for (iv=0;iv<geom_vmax;iv++) \
      if (S->name##ptr == &(S->geom[iv]))\
	hasg=hask=1;\
  }\
} /* PAR_USE_CHECK */

/****************/
/* Accept the parameters of the weights: manypoint, conservative, harmonic,	*/
/* check_positivity, the diffusivities D(=Ds), Dl(=Dpar), Dt(=Dtrans), which may	*/
/* be =~ expressions, and d0..d1; and determine the stencil size nnb and nD.	*/
/* Not accepted here: hx (accepted before), isolated (default is up to the	*/
/* device; *hasexe tells if any diffusivity is a =~ expression), w0..w1.	*/
/* Returns 1 on success, 0 on failure of an ACCEPT, as a create function.	*/
static int diff_weights_accept(STR *S, char *w, int *hasexe_out)
{
  int nnb;				/* num of neighbours in the stencil */
  int nD;				/* num of essential components in diff tensor */
  int iv;				/* small range counter */
  char name[maxname];			/* name of local k-var */

  ACCEPTI(manypoint,0,0,1);
  ACCEPTI(conservative,ONE,0,1);	/* default 1: ONE is 1 unless xmax<3 */
  ACCEPTI(harmonic,0,0,1);		/* default is arithmetic, for historical reason  */

  ACCEPTI(check_positivity,0,0,1);
  
  /*-----------------------------*/
  /* Variable parameters section */
  /* Prepare local symtable */
  int hask=0, hasu=0, hasg=0;		/* flags of what k-expressions used */
  int hasexe=0;				/* flag that any coefficient is a =~ expression */
  k_on();					CHK(NULL);
  S->loctb = tb_new();
  memcpy(S->loctb,deftb,sizeof(*deftb));
  if (ONE) {tb_insert_int_ro(S->loctb,"x",&(S->x));	CHK("x");}
  if (TWO) {tb_insert_int_ro(S->loctb,"y",&(S->y));	CHK("y");}
  if (TRI) {tb_insert_int_ro(S->loctb,"z",&(S->z));	CHK("z");}
  CALLOC(S->u,vmax,sizeof(real));
  for (iv=0;iv<vmax;iv++) {
    snprintf (name,maxname,"u%d",iv);
    tb_insert_real (S->loctb,name,&(S->u[iv]));  CHK(name);
  } /* for iv */
  if (geom_vmax) {
    CALLOC(S->geom,geom_vmax,sizeof(real));
    for (iv=0;iv<geom_vmax;iv++) {
      snprintf (name,maxname,"geom%d",iv);
      tb_insert_real (S->loctb,name,&(S->geom[iv]));  CHK(name);
    } /* for iv */
  } else { /* not geom_vmax */
    CALLOC(S->geom,1,sizeof(real));	/* w/o Geom, isTissue is everywhere within NETTO box */
    S->geom[0]=1.0;    
  } /* not geom_vmax */

  /* Accept the recalculatable parameters, with aliases */
  ACCEPTRR(Ds,D,0.,0.,RNONE); 		PAR_USE_CHECK(D);
  ACCEPTRR(Dpar,Dl,0.,0.,RNONE);	PAR_USE_CHECK(Dl);
  ACCEPTRR(Dtrans,Dt,0.,0.,RNONE);	PAR_USE_CHECK(Dt);

  if NOT(hasu) FREE(S->u);		/* need copy of u only if a parameter depends on it */
  if NOT(hasg) FREE(S->geom);		/* need copy of geom only if a parameter depends on it */
  if NOT(hask) FREE(S->loctb);		/* need local symtable only if a parameter depends on a local var */
  k_off();
  /* End of variable parameters section */
  /*------------------------------------*/

  /* Check parameters for relevance and consistency */
  if (ANISOTROPY_ON) {
    /* At least one anisotropic diffusion coefficients is mandatory */
    ASSERT((S->Dlptr!=NULL)||(S->Dlexe!=NULL) || (S->Dtptr!=NULL)||(S->Dtexe!=NULL));
    /* Scalar diffusion is only for singular points */
    if (find_key("D=",w)) {
      MESSAGE("/* Anisotropy is on, so the scalar diffusivity 'D' "
	      "will apply only at exceptional 'isotropic' points, if any. */\n");
    } else {
      S->Dptr=NULL;
      S->Dexe=NULL;
      S->Dsrc=NULL;
    }
    if (manypoint)
      EXPECTED_ERROR("/* Option 'manypoint' in anisotropy is not allowed. */\n");
  } else { /* not ANISOTROPY_on */
    /*------------------------------*/
    /* Scalar diffusion coefficient is mandatory */
    ASSERT((S->Dptr!=NULL)||(S->Dexe!=NULL));
    /* Anisotropic diffusion coefficients don't have role to play */
    if (find_key("Dl=",w) || find_key("Dpar=",w) || 
	find_key("Dt=",w) || find_key("Dtrans=",w)) {
      MESSAGE("The anisotropic diffusion coefficients 'Dl'='Dpar' and 'Dt'='Dtrans' "
	      "are unused when anisotropy is inactive. "
	      "The parameter(s) will be ignored.\n");
    } /* not ANISOTROPY_on */
    if (manypoint && dim==1) {
      MESSAGE("/* Parameter 'manypoint' is not applicable in 1D, and will be ignored.*/\n");
      manypoint=S->manypoint=0;
    }
    if (manypoint && S->loctb)
      EXPECTED_ERROR("Option 'manypoint' only works with constant isotropic diffusivity, "
		     "so is not allowed with space-dependent D.\n");
  } /* if ANISOTROPY_ON else */

  /* Determine number of neighbours and diff tensor components */
  if (ANISOTROPY_ON) {
    switch (dim) {
    case 1: EXPECTED_ERROR("don't know how to make 1D anisotropic"); break;
    case 2: nnb=a2nnb; nD=a2nD; break;
    case 3: nnb=a3nnb; nD=a3nD; break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    }
  } else if (manypoint==0) {
    switch (dim) {
    case 1: nnb=i1nnb; nD=i1nD; break;
    case 2: nnb=i2nnb; nD=i2nD; break;
    case 3: nnb=i3nnb; nD=i3nD; break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    }
  } else {
    switch (dim) {
    case 2: nnb=m2nnb; nD=m2nD; break;
    case 3: nnb=m3nnb; nD=m3nD; break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    }
  }
  S->nnb=nnb;
  S->nD=nD;
  MESSAGE("\x01/*nnb=%d, nD=%d*/ ",nnb,nD);

  /* Range of layers to keep diff tensor components, mandatory if var diff & MPI.	*/
  /* With constant Dl and Dt, the neighbours' tensors are computed from their	*/
  /* fibre directions in Geom, which is available in the halo, so no layers needed.	*/
  ACCEPTI(d0,-1,-1,vmax-1);
  ACCEPTI(d1,-1,-1,vmax-1);
  #if MPI
  if (S->loctb && d0<0)
    EXPECTED_ERROR("space-dependent diffusivity in MPI requires layers d0..d1 to keep the %d diff tensor component(s)",nD);
  #endif
  ASSERT((d0==-1 && d1==-1) || (0<=d0 && d0<=d1 && d1<vmax && d1-d0+1==nD));

  *hasexe_out=hasexe;
  return 1;
} /* diff_weights_accept */

/* Macros to print out the weights of the current point, */
/* for debug purposes only. Remove when not needed.      */
#include "diff_printweights.h"

/****************/
/* Adjacent tissue: tissue points beyond the border of the device space,	*/
/* neighbouring its points. With isolated=0, they are connected to	*/
/* the device space, so their diff tensor is needed, too.		*/
/* List them, with the tensor kept in layers d0..d1 if any, 		*/
/* otherwise in arrays of nD as for the listed points.			*/
/* Duplicates are possible, and are harmless. 				*/
static devicePoint *diff_adjacent (Space *s, int d0, int nD, size_t *nadjacent)
{
  devicePoint *adjacent=NULL;
  size_t ip, n=0;
  INT x, y, z, dx, dy, dz;
  int pass;
  for (pass=0;pass<2;pass++) {		/* count, then fill in */
    n=0;
    for (ip=0;ip<s->np;ip++) {
      for (dx=-ONE;dx<=ONE;dx++) {
	for (dy=-TWO;dy<=TWO;dy++) {
	  for (dz=-TRI;dz<=TRI;dz++) {
	    x=s->p[ip].x+dx;
	    y=s->p[ip].y+dy;
	    z=s->p[ip].z+dz;
	    if (!isTissue(x,y,z) || in_space(s,x,y,z)) continue;
	    if (pass) {
	      adjacent[n].x=x;
	      adjacent[n].y=y;
	      adjacent[n].z=z;
	      adjacent[n].u=New+ind(x,y,z,0);
	      adjacent[n].X=NULL;
	      if (d0>=0) {
		adjacent[n].Y=(void *)(New+ind(x,y,z,d0));
	      } else {
		CALLOC((adjacent[n].Y),nD*sizeof(real),1); /* NB Y is (void *) so pointing to size 1 */
	      }
	    }
	    n++;
	  } /* for dz */
	} /* for dy */
      } /* for dx */
    } /* for ip */
    if (pass==0) {
      if (n==0) break;
      CALLOC(adjacent,n,sizeof(devicePoint));
    }
  } /* for pass */
  *nadjacent=n;
  return adjacent;
} /* diff_adjacent */

/****************/
/* Component (i,j), 0-based, of the scaled diff tensor at point (x,y,z), computed	*/
/* from its fibre direction in Geom, for constant Dl and Dt; the same formulas as	*/
/* for the listed points in diff2d.h, diff3d.h. Geom is available in the halo,	*/
/* so this works for neighbours in other MPI subdomains without any haloSwap.	*/
static real diff_geom_tensor (INT x, INT y, INT z, int i, int j, real Ddh, real Dth, real Dsh)
{
  real f[3];
  f[0]=Geom[geom_ind(x,y,z,GEOM_FIBRE_1)];
  f[1]=Geom[geom_ind(x,y,z,GEOM_FIBRE_2)];
  f[2]=(dim==3)?Geom[geom_ind(x,y,z,GEOM_FIBRE_3)]:0.0; /* as in diff2d.h */
  if (isUnitVector(f[0],f[1],f[2],UNIT_VECTOR_TOLERANCE))
    return Ddh*f[i]*f[j] + ((i==j)?Dth:0.0);
  else /* exceptional 'isotropic' point */
    return (i==j)?Dsh:0.0;
} /* diff_geom_tensor */

/****************/
static int diff_connect (STR *S,Space *s)
{
  DEVICE_CONST(real,hx);		/* space step */
  DEVICE_CONST(int,manypoint);		/* 0/1 flag for small/large stencil */
  DEVICE_CONST(int,conservative);	/* 0/1 flag for conservative scheme */
  DEVICE_CONST(int,harmonic);		/* 0/1 flag for harmonic averaging */
  DEVICE_CONST(int,nnb);		/* number of neighbours in the stencil */
  DEVICE_CONST(int,w0);			/* layers range to keep */
  DEVICE_CONST(int,w1);			/*          the weights */
  DEVICE_CONST(int,nD);			/* number of diff tensor components */
  DEVICE_CONST(int,d0);			/* layers range to keep */
  DEVICE_CONST(int,d1);			/*    tensor components */
  DEVICE_CONST(int,isolated);		/* 0/1 flag for no-flux border of device space */

  size_t x,y,z,ip;			/* grid range counters */
  int i,j,nb;				/* short range counters */
  size_t np=s->np;			/* number of listed points of this device */
  devicePoint P;			/* the current point's structure */
  devicePoint *adjacent=NULL;		/* adjacent tissue, connected to the device space */
  size_t nadjacent=0;			/* .., their number */
  int ipass;				/* 0: listed points, 1: adjacent tissue */

  if (!isolated) adjacent=diff_adjacent(s,d0,nD,&nadjacent);
  S->nadjacent=nadjacent;
  
  assert(GEOM_FIBRE_1==1);		/* shorter */	
  assert(GEOM_FIBRE_2==2);		/* if rely */
  assert(GEOM_FIBRE_3==3);		/* on this */

 /* Create and initialize extra load arrays at all points.			*/
 /* Without the layers, they are slices of the device's own pools, which are	*/
 /* reallocated here, as the list of points may have been renewed by relisting	*/
 /* which does not free the arrays attached to the old points.			*/
  if (w0<0) {
    FREE(S->Xpool);
    if (np) CALLOC(S->Xpool,np*nnb,sizeof(real));
  }
  if (d0<0) {
    FREE(S->Ypool);
    if (np) CALLOC(S->Ypool,np*nD,sizeof(real));
  }
  for (ip=0;ip<np;ip++) {
    devicePoint *P=s->p+ip;

    if (w0<0) { /* no weights layers */
      P->X=(void *)(S->Xpool+ip*nnb);	/* zeroed by CALLOC */
    } else { /* yes weights layers */
      P->X=(void *)(P->u+w0*DV);
      memset(P->X, 0, nnb*sizeof(real));
    }

    if (d0<0) { /* no tensor layers */
      P->Y=(void *)(S->Ypool+ip*nD);	/* zeroed by CALLOC */
    } else { /* yes tensor layers */
      P->Y=(void *)(P->u+d0*DV);
      memset(P->Y, 0, nD*sizeof(real));
    }
    
  } /* for ip */

  /* Macros to access option-dependent and stencil-dependent arrays */
#define XP ((real *)(P.X))		/* array of weights */
#define YP ((real *)(P.Y))		/* array of diff tensor components */
#define WENUM(stc,nb) JOIN2(stc,nb)	/* weights enumerator, according to stencil stc */
#define DENUM(stc,ij) JOIN3(stc,D,ij)	/* tensor component enumerator, according to stencil stc */
#define W(nb) (XP[WENUM(STC,nb)])	/* stencil weight for current stencil STC */
#define D(ij) (YP[DENUM(STC,ij)])	/* diff tensor component for current stencil STC */ 
  /* Tensor component of neighbour computed from its fibre direction: for constant	*/
  /* Dl, Dt without the layers. ij is a two-digit literal, e.g. 12 -> i=0, j=1.	*/
  /* Ddh, Dth, Dsh are the local variables of diff?d*.h; declared dummies in 1D.	*/
#define DNB_GEOM(ij,...) (diff_geom_tensor(_(__VA_ARGS__),(ij)/10-1,(ij)%10-1,Ddh,Dth,Dsh))
  /* Whether the neighbour's tensor is to be computed from Geom */
#define TENSOR_FROM_GEOM (ANISO && !VARIA && d0<0)
#if MPI
  /* tensor component of neighbour: from the layers or from Geom */
#define DNB(ij,...) (TENSOR_FROM_GEOM?						\
  (DNB_GEOM(ij,__VA_ARGS__)):							\
  (New[APPLY(ind,(_(__VA_ARGS__),d0+DENUM(STC,ij)))]))
#else /* not MPI */
  /* size_t ipgrid[xlen*ylen*zlen];	/\* bookkeeping the points *\/ */
  size_t *ipgrid;
  CALLOC(ipgrid,xlen*ylen*zlen,sizeof(size_t));	/* bookkeeping the points */
#define IPNBR(x,y,z) (ipgrid[(z)+zlen*((y)+ylen*(x))])	/* number in the list of point with these coords */
/* tensor component of neighbour: from the layers if any, or from Geom if possible, */
/* otherwise from the listed point or adjacent tissue */
#define DNB(ij,...) ((d0>=0)?						\
  (New[APPLY(ind,(_(__VA_ARGS__),d0+DENUM(STC,ij)))]):			\
  TENSOR_FROM_GEOM?(DNB_GEOM(ij,__VA_ARGS__)):					\
  (((real *)((APPLY(IPNBR,(_(__VA_ARGS__)))<np)?				\
	      (s->p[APPLY(IPNBR,(_(__VA_ARGS__)))].Y):			\
	      (adjacent[APPLY(IPNBR,(_(__VA_ARGS__)))-np].Y)))[DENUM(STC,ij)]))
#endif

#define CLEAR_WEIGHTS memset(P.X, 0, nnb*sizeof(real))	/* zero out before filling them in */
  
  /* Beginning of a 4-fold #include chain that converts variables to CPP macros 		*/
  /* Semantics of these macro flags:  								*/
  /* GEOM: there is Geom array, i.e. geometry file, geometry pgm in 'state', or just 'aniso=1'.	*/
  /* ANISO: anisotropy is on. Must have been defined in 'state'. Implies GEOM. 			*/
  /* MANY: many-point stencil. Incompatible with anisotropy nor with variable diffusivity.	*/
  /* VARIA: space-dependent diffusivities (principal values of the tensor).			*/
  /* Nota Bene: tensor may be variable due to principal values OR fibre directions!		*/
  /* G(0,...) tells if a neighbour is connected: it must be tissue, and if isolated,		*/
  /* also in the device space, so the border of that is a no-flux boundary.			*/
  /* G(v,...), v>0 are the fibre direction components.						*/
#define CONNECTED(x,y,z) (isTissue(x,y,z) && (!isolated || in_space(s,x,y,z)))

  if (Geom==NULL) {			/* Not Geom, so not Aniso */
#define GEOM 0
#define ANISO 0
#define G(v,...) ((v)?(0):(APPLY(CONNECTED,(_(__VA_ARGS__)))))
    /* */				DB("D=%g=%g=*(%p)=(double *)exec(%p)\n",S->D,*(S->Dptr),S->Dptr,S->Dexe);
#include "diff1.h"
#undef G
#undef ANISO
#undef GEOM
  } else if (ANISOTROPY_ON==0) { 	/* Geom but not Aniso */
#define GEOM 1
#define ANISO 0
#define G(v,...) ((v)?(0):(APPLY(CONNECTED,(_(__VA_ARGS__)))))
#include "diff1.h"
#undef G
#undef ANISO
#undef GEOM
  } else {				/* both Geom and Aniso */
#define GEOM 1
#define ANISO 1
#define G(v,...) ((v)?(Geom[APPLY(geom_ind,(_(__VA_ARGS__),v))]):(APPLY(CONNECTED,(_(__VA_ARGS__)))))
#include "diff1.h"
#undef G
#undef ANISO
#undef GEOM
    if (S->check_positivity) {
      for (ip=0;ip<np;ip++) {
	P=s->p[ip];
	/* if (XP[0]>0) {} */
	if (XP[0]>0) {
	  MESSAGE("\n/* #%d t=%ld idev=%d='%s' bad point %ld:(%ld,%ld,%ld) W=\n" ,
		  mpi_rank,t,idev,dev[idev].n, ip,P.x,P.y,P.z);
	  switch (dim) {
	  case 2:
	    MESSAGE(
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n", 
		    XP[m2mp],XP[m20p],XP[m2pp], 
		    XP[m2m0],XP[m200],XP[m2p0], 
		    XP[m2mm],XP[m20m],XP[m2pm]); 
	    break;
	  case 3:
	    MESSAGE(
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "-----------\n" 
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "-----------\n" 
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "%g\t%g\t%g\n" 
		    "===========\n", 
		    0.0 ,XP[m30pp],0.0, 
		    XP[m3m0p],XP[m300p],XP[m3p0p], 
		    0.0 ,XP[m30mp],0.0, 
		    XP[m3mp0],XP[m30p0],XP[m3pp0], 
		    XP[m3m00],XP[m3000],XP[m3p00], 
		    XP[m3mm0],XP[m30m0],XP[m3pm0], 
		    0.0 ,XP[m30pm],0.0, 
		    XP[m3m0m],XP[m300m],XP[m3p0m], 
		    0.0 ,XP[m30mm],0.0 
		    ); 
	    break;
	  default: EXPECTED_ERROR("dimensionality %d is cannot be anisotropic so unexpected here",dim);
	  } /* switch dim */
	  MESSAGE(" The central weight is positivie: %lg, which will cause numerical instability. */\n",XP[0]);
	} /* if XP[0]<0 */
      } /* for ip */
    } /* if check_positivity */
  } /* Geom or Aniso? */
  /* End of #include chain */
  if (d0<0) for (ip=0;ip<nadjacent;ip++) FREE(adjacent[ip].Y);
  FREE(adjacent);
#if !MPI
  FREE(ipgrid);
#endif
  
  s->relisted=0;							DB("diff_connect done\n");
  return SUCCESS;
} /* diff_connect */

#endif /* DIFFWEIGHTS_CODE_H */
