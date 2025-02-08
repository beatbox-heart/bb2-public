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

/* diff.c: 
 * "Universal" diffusion device to supersede many (all) previous devices of the kind. 
 * Key new features: 
 * # The weights are in device memory OR in the dedicated layers of the main grid (# depending on stencil);
 * # The list of points defined by the same mechanism as generic device space definition;
 * # (Re)calculation of the weights by demand when the corresponding flag is raised, not just at the start;
 * # Several different methods of defining the weights, including some old and some new: 
 * # - isotropic without diagonal connections ("manypoint=0") (diff, ANISO=0);
 * # - isotropic with diagonal connections ("manypoint=1") (diff, ANISO=0);
 * # - ..., variable diffusivity that may be space and/or time-dependent;
 * # - anisotropic with diagonal connections based on fibre directions with Dl and Dt fixed pars;
 * # - ..., ..., Dl and Dt may be time and/or space-dependent; 
 * # All of these work in MPI, too. 
 * 
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
#include "geometry.h"

/* Handle tokes for CPP prescan */
#define JOIN2(a,b) a##b
#define JOIN3(a,b,c) a##b##c
#define JOIN4(a,b,c,d) a##b##c##d
#define APPLY(a,b) a b

/* Constants in stencils named as in popular formulas */
static const real sixth      = 1./6.;
static const real quarter    = 1./4.;
static const real twosixths  = 2./6.;
static const real half       = 1./2.;
static const real foursixths = 4./6.;

#ifndef UNIT_VECTOR_TOLERANCE
#define UNIT_VECTOR_TOLERANCE 0.01
#endif

/* Neighbourhoods and diff arrays descriptions 		*/
/* for different cases.					*/
/* ..nnb = # of neighbours,				*/
/* ..nD  = # of diff tensor components.			*/
/* NB: unique components listed lexically,		*/
/* then others added as synonyms, 			*/
/* for code convenience and clarity.			*/
/* STD?D??? macros will be used in diff3.h.		*/
/* 3-bit list of stencils is excessive in some cases.	*/

/* 1+3; 1 : 1D iso (can't be otherwise in 1D!) */
enum i1_nbrs {i10, i1p, i1m, i1nnb};
enum i1_D    {i1D11, i1nD};
#define STC1D000 i1
#define STC1D001 i1  
#define STC1D010 i1  
#define STC1D011 i1
#define STC1D100 i1  
#define STC1D101 i1  
#define STC1D110 i1  
#define STC1D111 i1
#define STC1D(a,m,v) JOIN4(STC1D,a,m,v)

/* 1+5; 1 : 2D iso */
enum i2_nbrs {i200, i2p0, i2m0, i20p, i20m, i2nnb};
enum i2_D    {i2D11, i2nD, i2D22=i2D11};
/* 1+9; 1 : 2D manypoint */
enum m2_nbrs {m200, m2p0, m2m0, m20p, m20m,
	      m2pp, m2pm, m2mp, m2mm, m2nnb};
enum m2_D    {m2D11, m2nD, m2D22=m2D11};
/* 1+9; 3 : 2D aniso */
enum a2_nbrs {a200, 
              a2p0, a2m0, a20p, a20m,
	      a2pp, a2pm, a2mp, a2mm, a2nnb};
enum a2_D    {a2D11, a2D12, a2D22, a2nD, a2D21=a2D12};
#define STC2D000 i2
#define STC2D001 i2  
#define STC2D010 m2  
#define STC2D011 m2
#define STC2D100 a2  
#define STC2D101 a2  
#define STC2D110 a2  
#define STC2D111 a2
#define STC2D(a,m,v) JOIN4(STC2D,a,m,v)

/* 1+7; 1 : 3D iso */
enum i3_nbrs {i3000, i3p00, i3m00, i30p0, i30m0, i300p, i300m, i3nnb};
enum i3_D    {i3D11, i3nD, i3D22=i3D11, i3D33=i3D11};
/* 1+19, 1 : 3D manypoint */
enum m3_nbrs {m3000, m3p00, m3m00, m30p0, m30m0, m300p, m300m,
	      m3pp0, m3pm0, m3mp0, m3mm0, m3p0p, m3p0m,
	      m3m0p, m3m0m, m30pp, m30pm, m30mp, m30mm, m3nnb};
enum m3_D    {m3D11, m3nD, m3D22=m3D11, m3D33=m3D11};
/* 1+19, 6 : 3D aniso */
enum a3_nbrs {a3000, a3p00, a3m00, a30p0, a30m0, a300p, a300m,
              a3pp0, a3pm0, a3mp0, a3mm0, a3p0p, a3p0m,
	      a3m0p, a3m0m, a30pp, a30pm, a30mp, a30mm, 
	      a3nnb};
enum a3_D    {a3D11, a3D12, a3D13, a3D22, a3D23, a3D33, a3nD, a3D21=a3D12, a3D31=a3D13, a3D32=a3D23};
#define STC3D000 i3
#define STC3D001 i3  
#define STC3D010 m3  
#define STC3D011 m3
#define STC3D100 a3  
#define STC3D101 a3  
#define STC3D110 a3  
#define STC3D111 a3
#define STC3D(a,m,v) JOIN4(STC3D,a,m,v)
  
typedef struct {
#include "k_code.h"			
  p_tb loctb;				/* Local symtable; NULL if none of PARs are k-expressions */
  real *u;				/* Copy of local u vector; NULL if none of PARs depends on it */
  real *geom;				/* ..., geom vector, ... */
  INT x, y, z;				/* Copy of current coords for (non-NULL only) loctb */
  PAR(real,Dl);				/* Parallel diffusion coefficient (ANISO only) */
  PAR(real,Dt);				/* Transverse diffusion coefficient (ANISO only) */
  PAR(real,D);				/* Scalar diffusivity for isotropic or singular anisotropic */
  real hx;				/* Space step */
  real ht;				/* Time step or 0 */
  int manypoint; 			/* == NINEPOINT, NINETEENPT in ezdiff < ezspiral, ezscroll */
  int v0;				/* Layer for the diffusive field */
  int v1;				/* Layer for the Laplacian */
  int nnb;				/* Num of neighbours in the stencil */
  ssize_t H[a3nnb];			/* Vector of displacements in the stencil; allocate for max */
  int w0, w1;				/* Layers range to keep the stencil weights */
  int nD;				/* Num of diff tensor components */
  int d0, d1;				/* Layers range to keep the diff tensor; must for var diff in MPI */
  Condition reweight_when;		/* When to recompute the weights: variable address */
  char reweight[maxname];		/* When to recompute the weights: variable name */
} STR;

/* (Re-)Compute the weights */
static int diff_connect (STR *S,Space *s); 

/****************/
RUN_HEAD(diff)
{
  DEVICE_CONST(real,ht);		/* time step */
  DEVICE_CONST(int,v0);			/* source layer (diffusive field) */
  DEVICE_CONST(int,v1);			/* target layer (laplacian field) */
  DEVICE_CONST(int,nnb);		/* number of neighbours in the stencil */
  DEVICE_ARRAY(ssize_t,H);		/* vector of neighbours' displacements */
  DEVICE_CONST(Condition,reweight_when);/* flag to recompute weights */
  size_t np=s.np;			/* number of points in the device's space */
  int V0=v0*DV;				/* displacement from anchor to source value */
  int V1=v1*DV;				/* displacement from anchor to target value */
  int ip;				/* points counter */
  devicePoint P;			/* point structure */
  real *X;				/* vector of weights */
  real *u, *un;				/* pointers to source and target */
  ssize_t H0,H1,H2,H3,H4,H5,H6,H7,H8,H9,
    H10,H11,H12,H13,H14,H15,H16,H17,H18;/* neighbours' displacements (max no) */

  /* Recompute weights if requested explicitly or the space has been relisted */
  if (*reweight_when || s.relisted) {
    if (SUCCESS != diff_connect(S,&s)) EXPECTED_ERROR("could not (re)make the connection weights");
  }

  /* Simple variables are faster than array elements */
  switch (nnb) {
  case  3: H0=H[0]; H1=H[1]; H2=H[2]; break; /* i1 */
  case  5: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; break; /* i2 */
  case  7: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; H5=H[5]; H6=H[6]; break; /* i3 */
  case  9: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; H5=H[5]; H6=H[6]; H7=H[7]; H8=H[8]; break; /* m2=a2 */
  case 19: H0=H[0]; H1=H[1]; H2=H[2]; H3=H[3]; H4=H[4]; H5=H[5]; H6=H[6]; H7=H[7]; H8=H[8]; H9=H[9];
    H10=H[10]; H11=H[11]; H12=H[12]; H13=H[13]; H14=H[14]; H15=H[15]; H16=H[16]; H17=H[17]; H18=H[18]; break; /* m3=a3 */
  default: EXPECTED_ERROR("unknown stencil size %d",nnb);
  }

#if 1
  /* TODO: swap nnb switch with ip loop for speed? */
  /* Paradoxically, in debug compile, loop-outside is faster */
  for (ip=0;ip<np;ip++) {
    P=s.p[ip];
    u=P.u+V0;
    un=P.u+V1;
    X=(real *)P.X;
    switch (nnb) {
    case  3: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]; break;
    case  5: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]; break;
    case  7: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]; break;
    case  9: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8]; break;
    case 19: *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8]
	+u[H9]*X[9]+u[H10]*X[10]+u[H11]*X[11]+u[H12]*X[12]+u[H13]*X[13]+u[H14]*X[14]+u[H15]*X[15]+u[H16]*X[16]
	+u[H17]*X[17]+u[H18]*X[18]; break;
    default: EXPECTED_ERROR("unknown stencil size %d",nnb);
    }
  }
#else
  switch (nnb) {
  case  3:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2];
    }
    break;
  case  5:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4];
    }
    break;
  case 7:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6];
    }
    break;
  case  9:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8];
    }
    break;
  case 19:
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u+V0;
      un=P.u+V1;
      X=(real *)P.X;
      *un=u[H0]*X[0]+u[H1]*X[1]+u[H2]*X[2]+u[H3]*X[3]+u[H4]*X[4]+u[H5]*X[5]+u[H6]*X[6]+u[H7]*X[7]+u[H8]*X[8]
	+u[H9]*X[9]+u[H10]*X[10]+u[H11]*X[11]+u[H12]*X[12]+u[H13]*X[13]+u[H14]*X[14]+u[H15]*X[15]+u[H16]*X[16]
	+u[H17]*X[17]+u[H18]*X[18];
    }
    break;
  default: EXPECTED_ERROR("unknown stencil size %d",nnb);
  } /* switch nnb */
#endif

  /* If nonzero timestep set, do the forward Euler step for diffusion equation */
  if (ht) {
#if MPI
    MPIDO(MPI_Barrier(ALL_ACTIVE_PROCS),"could not sync threads before applying timestep");
#endif
    for (ip=0;ip<np;ip++) {
      P=s.p[ip];
      u=P.u;
      u[V0] += ht*u[V1];
    } /* for ip */
  } /* if ht */
}
RUN_TAIL(diff)

/****************/
DESTROY_HEAD(diff)
/* TODO: there have been dynamically allocated objects; free them! */
DESTROY_TAIL(diff)

#define DB(...) DEBUG(__VA_ARGS__)

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

/* Check if given expression depends on local vars; */
/* hask: yes; hasu: in particular u*; hasg: in particular geom* */
#define PAR_USE_CHECK(name) { \
  if (S->name##exe) { \
    char varname[maxname]; \
    int iv; \
    hask=1; \
    for (iv=0;iv<vmax;iv++) { \
      snprintf(varname,maxname,"u%d",iv); \
      if (k_expr_depends (S->name##src,varname)) hasu=1; \
    } /* for iv */ \
    for (iv=0;iv<geom_vmax;iv++) { \
      snprintf(varname,maxname,"geom%d",iv); \
      if (k_expr_depends (S->name##src,varname)) hasg=1; \
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
CREATE_HEAD(diff)
{
  Space *s=&(dev->s);
  DEVICE_REQUIRES_SYNC;
  if (s->nowhere) EXPECTED_ERROR("'nowhere' is not allowed\n");
  if (0==s->listpoints==0) 
    MESSAGE("/* listpoints=1 is mandatory in this device and enforced */\n");
  else if (find_key("listpoints=",w))
    MESSAGE("/* listpoints=1 is mandatory anyway*/\n");
  s->listpoints=1;
  s->p=NULL; /* signal that it has not be alloc'd yet */
  if (SUCCESS != list_space(s))
    ABORT("could not make list of points");
    
  size_t np=s->np;			/* number of listed points of this device */
  int nnb;				/* num of neighbours in the stencil */
  int nD;				/* num of essential components in diff tensor */
  size_t ip;				/* listed points counter */
  int nb,iv;				/* small range counters */
  int reweight_entry;			/* deftb entry of 'reweight' k-var */
  char name[maxname];			/* name of local k-var */

  ACCEPTI(v0,INONE,0,vmax-1);		/* source layer: diffusive field */
  ACCEPTI(v1,INONE,0,vmax-1);		/* target layer: laplacian field */
  ASSERT( v1 != v0 );			/* source and target must be distinct */

  ACCEPTR(hx,RNONE,RSUCC(0.),RNONE);	/* diffusion makes no sense without space step */
  ACCEPTR(ht,0.,0.,RNONE);		/* make forward Euler step if ht>0 */
  ACCEPTI(manypoint,0,0,1);

  ACCEPTS(reweight,"never");		/* global k-var to flag that weights to be recomputed */
  if (0==(reweight_entry=tb_find(deftb,S->reweight)))
    EXPECTED_ERROR("unknown symbol %s for 'when' parameter",S->reweight);
  if ((tb_flag(deftb,reweight_entry) & f_vb)==0)
    EXPECTED_ERROR("symbol %s used for 'when' parameter is not a variable",S->reweight);
  S->reweight_when = (Condition)tb_addr(deftb,reweight_entry); 
  
  /*-----------------------------*/
  /* Variable parameters section */
  /* Prepare local symtable */
  int hask=0, hasu=0, hasg=0;		/* flags of what k-expressions used */
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

  /* Range of layers to keep stencil weights */
  ACCEPTI(w0,-1,-1,vmax-1);
  ACCEPTI(w1,-1,-1,vmax-1);
  ASSERT((w0==-1 && w1==-1) || (0<=w0 && w0<=w1 && w1<vmax && w1-w0+1==nnb));

  /* Range of layers to keep diff tensor components, mandatory if var diff & MPI */
  ACCEPTI(d0,-1,-1,vmax-1);
  ACCEPTI(d1,-1,-1,vmax-1);
  #if MPI
  if (S->loctb) ASSERT(0<=d0 && d0<=d1 && d1<vmax && d1-d0+1==nD);
  #else
  ASSERT((d0==-1 && d1==-1) || (0<=d0 && d0<=d1 && d1<vmax && d1-d0+1==nD));
  #endif
  
  /* Determine displacements to neighbours == stencil elements */
  ssize_t *H;
  CALLOC(H,nnb,sizeof(ssize_t));
  if (ANISOTROPY_ON) {			/* anisotropic */
    switch (dim) {
      EXPECTED_ERROR("don't know how to make 1D anisotropic");
      break;
    case 2:
      H[a200]=0;
      H[a2p0]=+DX; H[a2m0]=-DX;
      H[a20p]=+DY; H[a20m]=-DY;
      H[a2pp]=+DX+DY; H[a2pm]=+DX-DY; H[a2mp]=-DX+DY; H[a2mm]=-DX-DY;
      break;
    case 3:
      H[a3000]=0;
      H[a3p00]=+DX; H[a3m00]=-DX;
      H[a30p0]=+DY; H[a30m0]=-DY;
      H[a300p]=+DZ; H[a300m]=-DZ;
      H[a3pp0]=+DX+DY; H[a3pm0]=+DX-DY; H[a3mp0]=-DX+DY; H[a3mm0]=-DX-DY;
      H[a3p0p]=+DX+DZ; H[a3p0m]=+DX-DZ; H[a3m0p]=-DX+DZ; H[a3m0m]=-DX-DZ;
      H[a30pp]=+DY+DZ; H[a30pm]=+DY-DZ; H[a30mp]=-DY+DZ; H[a30mm]=-DY-DZ;
      break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    } /* switch dim */
  } else if (manypoint==0) {		/* isotropic, small stencil */
    switch (dim) {
    case 1:
      H[i10]=0;
      H[i1p]=+DX; H[i1m]=-DX;
      break;
    case 2:
      H[i200]=0;
      H[i2p0]=+DX; H[i2m0]=-DX;
      H[i20p]=+DY; H[i20m]=-DY;
      break;
    case 3:
      H[i3000]=0;
      H[i3p00]=+DX; H[i3m00]=-DX;
      H[i30p0]=+DY; H[i30m0]=-DY;
      H[i300p]=+DZ; H[i300m]=-DZ;
      break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    } /* switch dim */
  } else { 				/* isotropic, large stencil */
    switch (dim) {
    case 2:
      H[m200]=0;
      H[m2p0]=+DX; H[m2m0]=-DX;
      H[m20p]=+DY; H[m20m]=-DY;
      H[m2pp]=+DX+DY; H[m2pm]=+DX-DY; H[m2mp]=-DX+DY; H[m2mm]=-DX-DY;
      break;
    case 3:
      H[m3000]=0;
      H[m3p00]=+DX; H[m3m00]=-DX;
      H[m30p0]=+DY; H[m30m0]=-DY;
      H[m300p]=+DZ; H[m300m]=-DZ;
      H[m3pp0]=+DX+DY; H[m3pm0]=+DX-DY; H[m3mp0]=-DX+DY; H[m3mm0]=-DX-DY;
      H[m3p0p]=+DX+DZ; H[m3p0m]=+DX-DZ; H[m3m0p]=-DX+DZ; H[m3m0m]=-DX-DZ;
      H[m30pp]=+DY+DZ; H[m30pm]=+DY-DZ; H[m30mp]=-DY+DZ; H[m30mm]=-DY-DZ;
      break;
    default: EXPECTED_ERROR("unexpected dimension %d",dim);
    } /* switch dim */
  } /* if ANISOTROPY_ON else if manypoint==0 else */
  memcpy(S->H,H,nnb*sizeof(ssize_t *));

  /* Create and initialize extra load arrays at all points */
  for (ip=0;ip<np;ip++) {
    devicePoint *P=s->p+ip;
    FREE(P->X);
    if (w0<0) { /* no weights layers */
      CALLOC((P->X),nnb*sizeof(real),1); /* NB P->X is (void *) so pointing to size 1 */
    } else { /* yes weights layers */
      P->X=(void *)(P->u+w0*DV);
      memset(P->X, 0, nnb*sizeof(real));
    }
    FREE(P->Y);
    if (d0<0) { /* no tensor layers */
      CALLOC((P->Y),nD*sizeof(real),1); /* NB P->Y is (void *) so pointing to size 1 */
    } else { /* yes tensor layers */
      P->Y=(void *)(P->u+d0*DV);
      memset(P->Y, 0, nD*sizeof(real));
    }
  }
  
  /* Now compute and assign connection weights: this is potentially doable at each step */
  if (SUCCESS != diff_connect(S,s)) ABORT("could not make the connection weighs");
}
CREATE_TAIL(diff,1);

/* Macros to print out the weights of the current point, */
/* for debug purposes only. Remove when not needed.      */
#include "diff_printweights.h"

/****************/
static int diff_connect (STR *S,Space *s)
{
  DEVICE_CONST(real,hx);		/* space step */
  DEVICE_CONST(int,manypoint);		/* 0/1 flag for small/large stencil */
  DEVICE_CONST(int,nnb);		/* number of neighbours in the stencil */
  DEVICE_CONST(int,w0);			/* layers range to keep */
  DEVICE_CONST(int,w1);			/*          the weights */
  DEVICE_CONST(int,nD);			/* number of diff tensor components */
  DEVICE_CONST(int,d0);			/* layers range to keep */
  DEVICE_CONST(int,d1);			/*    tensor components */
  
  size_t x,y,z,ip;			/* grid range counters */
  int i,j,nb;				/* short range counters */
  size_t np=s->np;			/* number of listed points of this device */
  devicePoint P;			/* the current point's structure */
  
  assert(GEOM_FIBRE_1==1);		/* shorter */	
  assert(GEOM_FIBRE_2==2);		/* if rely */
  assert(GEOM_FIBRE_3==3);		/* on this */

  /* Macros to access option-dependent and stencil-dependent arrays */
#define XP ((real *)(P.X))		/* array of weights */
#define YP ((real *)(P.Y))		/* array of diff tensor components */
#define WENUM(stc,nb) JOIN2(stc,nb)	/* weights enumerator, according to stencil stc */
#define DENUM(stc,ij) JOIN3(stc,D,ij)	/* tensor component enumerator, according to stencil stc */
#define W(nb) (XP[WENUM(STC,nb)])	/* stencil weight for current stencil STC */
#define D(ij) (YP[DENUM(STC,ij)])	/* diff tensor component for current stencil STC */ 
#if MPI
#define DNB(ij,...) (New[APPLY(ind,(_(__VA_ARGS__),d0+DENUM(STC,ij)))]) /* tensor component of neighbour */
#else /* not MPI */
  size_t ipgrid[xlen*ylen*zlen];	/* bookkeeping the points */
#define IPNBR(x,y,z) (ipgrid[(z)+zlen*((y)+ylen*(x))])	/* number in the list of point with these coords */
#define DNB(ij,...) (((real *)(s->p[APPLY(IPNBR,(_(__VA_ARGS__)))].Y))[DENUM(STC,ij)]) /* tensor component of neighbour */
#endif

#define CLEAR_WEIGHTS memset(P.X, 0, nnb*sizeof(real))	/* zero out before filling them in */
  
  /* Beginning of a 4-fold #include chain that converts variables to CPP macros 		*/
  /* Semantics of these macro flags:  								*/
  /* GEOM: there is Geom array, i.e. geometry file, geometry pgm in 'state', or just 'aniso=1'.	*/
  /* ANISO: anisotropy is on. Must have been defined in 'state'. Implies GEOM. 			*/
  /* MANY: many-point stencil. Incompatible with anisotropy nor with variable diffusivity.	*/
  /* VARIA: space-dependent diffusivities (principal values of the tensor).			*/
  /* Nota Bene: tensor may be vary due to principal values OR fibre directions!			*/
  
  if (Geom==NULL) {			/* Not Geom, so not Aniso */
#define GEOM 0
#define ANISO 0
    /* */				/* isTissue are all points in the device box  */
#define IN(x,y,z) (s->global_x0<=x && x<=s->global_x1 && \
		   s->global_y0<=y && y<=s->global_y1 && \
		   s->global_z0<=z && z<=s->global_z1)
#define G(v,...) ((v)?(0):(APPLY(IN,(_(__VA_ARGS__)))))
    /* */				DB("D=%g=%g=*(%p)=(double *)exec(%p)\n",S->D,*(S->Dptr),S->Dptr,S->Dexe);    
#include "diff1.h"
#undef G
#undef ANISO
#undef IN
#undef GEOM 
  } else if (ANISOTROPY_ON==0) { 	/* Geom but not Aniso */
#define GEOM 1
#define ANISO 0
#define G(v,...) Geom[APPLY(geom_ind,(_(__VA_ARGS__),v))] /* isTissue are points in Geom array */
#include "diff1.h"
#undef G
#undef ANISO
#undef GEOM
  } else {				/* both Geom and Aniso */
#define GEOM 1
#define ANISO 1
#define G(v,...) Geom[APPLY(geom_ind,(_(__VA_ARGS__),v))] /* isTissue are points in Geom array */
#include "diff1.h"
#undef G
#undef ANISO
#undef GEOM
  } /* Geom or Aniso? */
  /* End of #include chain */
  
  s->relisted=0;							DB("diff_connect done\n");
  return SUCCESS;
} /* diff_connect */
