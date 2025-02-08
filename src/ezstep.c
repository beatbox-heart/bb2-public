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

/* Wrapper around ezstep & ezstep3d from D. Barkley's ezspiral and ezscroll.  */
/* For now, periodic directions are not allowed.                              */
/* To do: MPI                                                                 */

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

typedef struct {
  /* Front-end, specified in input file */
  #define _(Type,Name,Dflt,miN,maX) Type Name;
  #include "ezstep.h"
  #undef _
  /* Back-end, ezspiral/ezscroll computed constants */
  real dt;
  real one_o_a;
  real b_o_a;
  real dt_o_eps;
  real dt_o_2eps;
  real dt_o_wh2;
  real dtDv_o_wh2;
  int u_diff_on;		/* if 1 then u-field diffuses */
  int v_diff_on;		/* if 1 then v-field diffuses */
  /* Auxiliary variables */
  int s1;
  int s2;
  int v0, v1;
  int first_time;
  int ezdim;
} STR;

/* True constants */
static const real zero       = 0.;
static const real one        = 1.;
static const real two        = 2.;
static const real four       = 4.;
static const real six        = 6.;
static const real twenty     = 20.;
static const real twentyfour = 24.;

/* EZ/Beatbox matching macros */
#define NX (xmax-2)
#define NY (ymax-2)
#define NZ (zmax-2)
#define I_INC (DX)
#define J_INC (DY)
#define K_INC (DZ)

static void Step_ini_2D (Space s, STR *S);
static void Step_ini_3D (Space s, STR *S);
static void Impose_boundary_conditions_2D (Space s, STR *S, real *w, real *sigma_w);
static void Impose_boundary_conditions_3D (Space s, STR *S, real *w, real *sigma_w);

/****************/
RUN_HEAD(ezstep)
{
  DEVICE_CONST(int, v0);
  DEVICE_CONST(real,delta);
  DEVICE_CONST(int,v_diff_on);
  DEVICE_CONST(int,explicit);
  DEVICE_CONST(int,split);
  DEVICE_CONST(int,manypoint);
  DEVICE_CONST(real,dt);
  DEVICE_CONST(real,one_o_a);
  DEVICE_CONST(real,b_o_a);
  DEVICE_CONST(real,dt_o_eps);
  DEVICE_CONST(real,dt_o_2eps);
  DEVICE_CONST(real,dt_o_wh2);
  DEVICE_CONST(real,dtDv_o_wh2);
  DEVICE_VAR(int,s1);
  DEVICE_VAR(int,s2);
  DEVICE_VAR(int,first_time);
  DEVICE_CONST(int,ezdim);
  real u_thresh;
  register int x, y, z;
  register int index, index_1, index_2;
  real *u       = &(New[ind(0,0,0,v0+0)]);
  real *v       = &(New[ind(0,0,0,v0+1)]);
  real *sigma_u = &(New[ind(0,0,0,v0+2)]);
  real *sigma_v = &(New[ind(0,0,0,v0+4)]);

  if (*first_time) {
    /* Initialize arrays for timestepping */
    switch (ezdim) {
    case 0: break;
    case 2: Step_ini_2D(s,S); break;
    case 3: Step_ini_3D(s,S); break;
    default: EXPECTED_ERROR("ezstep is not implemented for %d dimensions\n",dim);
    }
    *first_time=0;
  }

  if (ezdim) {
    /* Interchange s1 and s2 */
    (*s1) = 1-(*s1);
    (*s2) = 1-(*s2);
  }

  /* A 5-fold #include chain that converts integer flags to CPP macros */
  if (v_diff_on) {
    #define V_DIFF_ON 1
    #include "ezstep1.h"
    #undef V_DIFF_ON
  } else {
    #define V_DIFF_ON 0
    #include "ezstep1.h"
    #undef V_DIFF_ON
  } /* if diff_on else */
  
}
RUN_TAIL(ezstep);

DESTROY_HEAD(ezstep)
DESTROY_TAIL(ezstep)

CREATE_HEAD(ezstep)
{
  DEVICE_REQUIRES_SYNC;
  DEVICE_IS_RECTANGULAR;
  
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  real WEIGHT, STABILITY_LIMIT;

  /* Front-end parameters are read from the script */
  #define acceptreal(b,c,d,e) if (!acceptr(#b"=",&(S->b),c,d,e,w)) return(0); real b=S->b
  #define acceptint(b,c,d,e) if (!accepti(#b"=",&(S->b),c,d,e,w)) return(0); int b=S->b 
  #define _(Type,Name,Dflt,miN,maX) accept##Type(Name,Dflt,miN,maX);
  #include "ezstep.h"
  #undef _
  #undef acceptint
  #undef acceptreal

  /* Some checks on the front-end parameters */
  if (D!=0.0) {
    ASSERT(hx>0);
    if (manypoint) {
      WEIGHT=6.;
      STABILITY_LIMIT=3.*hx*hx/8.;
    } else {
      WEIGHT=1.;
      switch (dim) {
      case 2: STABILITY_LIMIT=hx*hx/4.; break;
      case 3: STABILITY_LIMIT=hx*hx/6.; break;
      default: EXPECTED_ERROR("%s not implemented for %d dimensions\n",dev->n,dim);
      }
    }
    ASSERT(ht < STABILITY_LIMIT);
  } else {
    if (find_key("hx=",w)) MESSAGE("/* Warning: D=0 so the given hx value will be ignored */\n");
    if (find_key("split=",w)) MESSAGE("/* Warning: D=0 so the given split value will be ignored */\n");
    if (find_key("manypoint=",w)) MESSAGE("/* Warning: D=0 so the given manypoint value will be ignored */\n");
    if (find_key("pbc_x=",w)) MESSAGE("/* Warning: D=0 so the given pbc_x value will be ignored */\n");
    if (find_key("pbc_y=",w)) MESSAGE("/* Warning: D=0 so the given pbc_y value will be ignored */\n");
    if (find_key("pbc_z=",w)) MESSAGE("/* Warning: D=0 so the given pbc_z value will be ignored */\n");
  }
  ASSERT(a != 0);
  ASSERT(eps != 0);
  
  /* Back-end parameters calculated from the front-end ones */
  #define _(var,expr) int var=expr; S->var=var; printf("%s=%d\n",#var,var)
  _(u_diff_on, (D!=0));
  _(v_diff_on, (Dv!=0));
  _(ezdim , (u_diff_on || v_diff_on)?dim:0);
  #undef _  
  switch (u_diff_on*1 + v_diff_on*2) {
  case 0*1+0*2:
    if (v1 - v0 < 1)
      EXPECTED_ERROR("%s device without either diffusion requires two contiguous layers; "
		     "given %d..%d do not satisfy that\n",dev->n,v0,v1);
    break;
  case 0*1+1*2:
    EXPECTED_ERROR("only second diffusion requested, this is not implemented\n");
    break;
  case 1*1+0*2:
    if (v1 - v0 < 3)
      EXPECTED_ERROR("%s device without second diffusion requires four contiguous layers; "
		     "given %d..%d do not satisfy that\n",dev->n,v0,v1);
    break;
  case 1*1+1*2:
    if (v1 - v0 < 5)
      EXPECTED_ERROR("%s device with second diffusion requires six contiguous layers; "
		     "given %d..%d do not satisfy that\n",dev->n,v0,v1);
    break;
  default: EXPECTED_ERROR("impossible case\n");
  }
  #define _(var,expr) real var=expr; S->var=var; printf("%s=%g\n",#var,var)
  _(dt , ht);
  _(one_o_a , 1.0/a);
  _(b_o_a , b/a);
  _(dt_o_eps , dt/eps);
  _(dt_o_2eps  , dt_o_eps/2.);
  _(dt_o_wh2 , D?(D*dt/(WEIGHT*hx*hx)):0);
  _(dtDv_o_wh2 , D?(Dv/D*dt_o_wh2):0);
  #undef _

  /* Initialize arrays for timestepping at the first run, not now */
  S->first_time=1;
}
CREATE_TAIL(ezstep,1);
/* ========================================================================= */;

static void Step_ini_2D (Space s, STR *S)
{
  DEVICE_CONST(int,manypoint);
  DEVICE_CONST(int,v_diff_on);
  DEVICE_VAR(int,s1);
  DEVICE_VAR(int,s2);
  DEVICE_CONST(int,v0);
  int x, y;
  int index, index_1, index_2;
  real *u=&(New[ind(0,0,0,v0)]);
  real *v=&(New[ind(0,0,0,v0+1)]);
  real *sigma_u=&(New[ind(0,0,0,v0+2)]);
  real *sigma_v=v_diff_on?(&(New[ind(0,0,0,v0+4)])):NULL;

  /* Set initial s1 and s2 */
  (*s1) = 0;
  (*s2) = 1;

  /* Initialize spatial sums. */
  switch (2*v_diff_on+1*manypoint) {
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = zero;
  case 2*0+1*0:
    for (y=s.y0;y<=s.y1;y++) {
      for (x=s.x0;x<=s.x1;x++) {
	index = ind(x,y,0,0);
	index_1 = index + (*s1)*DV;
	index_2 = index + (*s2)*DV;
	ADD_TO_U_SUM5(index,index_2);
	ZERO_USED_SUMS(index_1);
      }
    }
    break;
  case 2*0+1*1:
    for (y=s.y0;y<=s.y1;y++) {
      for (x=s.x0;x<=s.x1;x++) {
	index = ind(x,y,0,0);
	index_1 = index + (*s1)*DV;
	index_2 = index + (*s2)*DV;
	ADD_TO_U_SUM9(index,index_2);
	ZERO_USED_SUMS(index_1);
      }
    }
    break;
  #undef ZERO_USED_SUMS
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = sigma_v[index_1] = zero;
  case 2*1+1*0:
    for (y=s.y0;y<=s.y1;y++) {
      for (x=s.x0;x<=s.x1;x++) {
	index = ind(x,y,0,0);
	index_1 = index + (*s1)*DV;
	index_2 = index + (*s2)*DV;
	ADD_TO_U_SUM5(index,index_2);
	ADD_TO_V_SUM5(index,index_2);
	ZERO_USED_SUMS(index_1);
      }
    }
    break;
  case 2*1+1*1:
    for (y=s.y0;y<=s.y1;y++) {
      for (x=s.x0;x<=s.x1;x++) {
	index = ind(x,y,0,0);
	index_1 = index + (*s1)*DV;
	index_2 = index + (*s2)*DV;
	ADD_TO_U_SUM9(index,index_2);
	ADD_TO_V_SUM9(index,index_2);
	ZERO_USED_SUMS(index_1);
      }
    }
    break;
  #undef ZERO_USED_SUMS
  }/* switch v_diff_on and manypoint */
  Impose_boundary_conditions_2D(s, S, u, sigma_u);
  if (v_diff_on) 
    Impose_boundary_conditions_2D(s, S, v, sigma_v);
} /* Step_ini_2D */
/* ========================================================================= */

static void Step_ini_3D (Space s, STR *S)
{
  DEVICE_CONST(int,manypoint);
  DEVICE_CONST(int,v_diff_on);
  DEVICE_VAR(int,s1);
  DEVICE_VAR(int,s2);
  DEVICE_CONST(int,v0);
  int x, y, z;
  int index, index_1, index_2;
  real *u=&(New[ind(0,0,0,v0)]);
  real *v=&(New[ind(0,0,0,v0+1)]);
  real *sigma_u=&(New[ind(0,0,0,v0+2)]);
  real *sigma_v=v_diff_on?(&(New[ind(0,0,0,v0+4)])):NULL;

  /* Set initial (*s1) and (*s2) */
  (*s1) = 0;
  (*s2) = 1;

  /* Initialize spatial sums. */
  switch (2*v_diff_on+1*manypoint) {
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = zero;
  case 2*0+1*0:
    for (z=s.z0;z<=s.z1;z++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (x=s.x0;x<=s.x1;x++) {
	  index = ind(x,y,z,0);
	  index_1 = index + (*s1)*DV;
	  index_2 = index + (*s2)*DV;
	  ADD_TO_U_SUM7(index,index_2);
	  ZERO_USED_SUMS(index_1);
	}
      }
    }
    break;
  case 2*0+1*1:
    for (z=s.z0;z<=s.z1;z++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (x=s.x0;x<=s.x1;x++) {
	  index = ind(x,y,z,0);
	  index_1 = index + (*s1)*DV;
	  index_2 = index + (*s2)*DV;
	  ADD_TO_U_SUM19(index,index_2);
	  ZERO_USED_SUMS(index_1);
	}
      }
    }
    break;
  #undef ZERO_USED_SUMS
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = sigma_v[index_1] = zero;
  case 2*1+1*0:
    for (z=s.z0;z<=s.z1;z++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (x=s.x0;x<=s.x1;x++) {
	  index = ind(x,y,z,0);
	  index_1 = index + (*s1)*DV;
	  index_2 = index + (*s2)*DV;
	  ADD_TO_U_SUM7(index,index_2);
	  ADD_TO_V_SUM7(index,index_2);
	  ZERO_USED_SUMS(index_1);
	}
      }
    }
    break;
  case 2*1+1*1:
    for (z=s.z0;z<=s.z1;z++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (x=s.x0;x<=s.x1;x++) {
	  index = ind(x,y,z,0);
	  index_1 = index + (*s1)*DV;
	  index_2 = index + (*s2)*DV;
	  ADD_TO_U_SUM19(index,index_2);
	  ADD_TO_V_SUM19(index,index_2);
	  ZERO_USED_SUMS(index_1);
	}
      }
    }
    break;
  #undef ZERO_USED_SUMS
  }/* switch v_diff_on and manypoint */
  Impose_boundary_conditions_3D(s, S, u, sigma_u);
  if (v_diff_on) 
    Impose_boundary_conditions_3D(s, S, v, sigma_v);
} /* Step_ini_3D */
/* ========================================================================= */


static void  Impose_boundary_conditions_2D (Space s, STR *S, real *w, real *sigma_w)
{
  DEVICE_CONST(int,pbc_x);
  DEVICE_CONST(int,pbc_y);
  DEVICE_CONST(int,manypoint);
  DEVICE_CONST(int,s2);
  int x,y;
#define W(x,y)          w[ind(x,y,0,0)]
#define Sigma_w(s,x,y)  sigma_w[ind(x,y,0,s)]

  
  /* First the values of u on the "real" part of the mesh are copied to
   * "fictitious" boundary points.  The way points are copied depend on
   * whether Neumann or periodic boundary conditions are being imposed.
   *
   * Note that the indices i, j do not range over the same values in each
   * case (x, y).  This is correct though a bit subtle. */

  /* Set fictitious points in x-direction */
  if (pbc_x) {
    for (y=s.y0;y<=s.y1;y++) { 
      W(s.x0-1,y) = W(s.x1-1,y);  
      W(s.x1+1,y) = W(s.x0+1,y);  
    }
  } else {
    for (y=s.y0;y<=s.y1;y++) { 
      W(s.x0-1,y) = W(s.x0+1,y);  
      W(s.x1+1,y) = W(s.x1-1,y);  
    }
  }

  /* Set fictitious points in y-direction */
  if (pbc_y) {
    for (x=s.x0-1;x<=s.x1+1;x++) { 
      W(x,s.y0-1) = W(x,s.y1-1);  
      W(x,s.y1+1) = W(x,s.y0+1);  
    }
  } else {
    for (x=s.x0-1;x<=s.x1+1;x++) { 
      W(x,s.y0-1) = W(x,s.y0+1);  
      W(x,s.y1+1) = W(x,s.y1-1);  
    }
  }

  if (manypoint) {
    
    /* -----------------------------
     * Nine-point Laplacian formulas
     * ----------------------------- */
    
    for (x=s.x0;x<=s.x1;x++) {
      Sigma_w(s2,x  ,s.y0) += four*W(x,s.y0-1);
      Sigma_w(s2,x+1,s.y0) +=      W(x,s.y0-1);  
      Sigma_w(s2,x-1,s.y0) +=      W(x,s.y0-1); 
      
      Sigma_w(s2,x  ,s.y1) += four*W(x,s.y1+1); 
      Sigma_w(s2,x+1,s.y1) +=      W(x,s.y1+1);  
      Sigma_w(s2,x-1,s.y1) +=      W(x,s.y1+1); 
    } 
    
    for (y=s.y0;y<=s.y1;y++) {
      Sigma_w(s2,s.x0,y  ) += four*W(s.x0-1,y); 
      Sigma_w(s2,s.x0,y+1) +=      W(s.x0-1,y);  
      Sigma_w(s2,s.x0,y-1) +=      W(s.x0-1,y); 
      
      Sigma_w(s2,s.x1,y  ) += four*W(s.x1+1,y); 
      Sigma_w(s2,s.x1,y+1) +=      W(s.x1+1,y);  
      Sigma_w(s2,s.x1,y-1) +=      W(s.x1+1,y); 
    }
    
    Sigma_w(s2,s.x0,s.y0) += W(s.x0-1,s.y0-1);  
    Sigma_w(s2,s.x1,s.y0) += W(s.x1+1,s.y0-1);  
    Sigma_w(s2,s.x0,s.y1) += W(s.x0-1,s.y1+1);  
    Sigma_w(s2,s.x1,s.y1) += W(s.x1+1,s.y1+1);  

  } else { /* manypoint */

    /* -----------------------------
     * Five-point Laplacian formulas 
     * ----------------------------- */
    
    for (x=s.x0;x<=s.x1;x++) {
      Sigma_w(s2,x,s.y0) += W(x,s.y0-1);  
      Sigma_w(s2,x,s.y1) += W(x,s.y1+1); 
    }
    
    for (y=s.y0;y<=s.y1;y++) {
      Sigma_w(s2,s.x0,y) += W(s.x0-1,y); 
      Sigma_w(s2,s.x1,y) += W(s.x1+1,y); 
    }
  } /* if manypoint else */
#undef W
#undef Sigma_w
} /* Impose_boundary_conditions_2D */
/* ========================================================================= */

static void  Impose_boundary_conditions_3D (Space s, STR *S, real *w, real *sigma_w)
{
  DEVICE_CONST(int,pbc_x);
  DEVICE_CONST(int,pbc_y);
  DEVICE_CONST(int,pbc_z);
  DEVICE_CONST(int,manypoint);
  DEVICE_CONST(int,s2);
  int x,y,z;
#define W(i,j,k)          w[ind(i,j,k,0)]
#define Sigma_w(s,i,j,k)  sigma_w[ind(i,j,k,s)]


  /* First the values of w on the "real" part of the mesh are copied to
   * "fictitious" boundary points.  The way points are copied depend on
   * whether Neumann or periodic boundary conditions are being imposed.
   *
   * Note that the indices i, j, k do not range over the same values in each
   * case (x, y, z).  This is correct though a bit subtle. */

  /* Set fictitious points in x-direction */
  if (pbc_x) {
    for (z=s.z0;z<=s.z1;z++) { 
      for (y=s.y0;y<=s.y1;y++) { 
	W(s.x0-1,y,z) = W(s.x1-1,y,z);  
	W(s.x1+1,y,z) = W(s.x0+1,y,z);  
      }
    }
  } else {
    for (z=s.z0;z<=s.z1;z++) { 
      for (y=s.y0;y<=s.y1;y++) { 
	W(s.x0-1,y,z) = W(s.x0+1,y,z);  
	W(s.x1+1,y,z) = W(s.x1-1,y,z);  
      }
    }
  }

  /* Set fictitious points in y-direction */
  if (pbc_y) {
    for (z=s.z0;z<=s.z1;z++) { 
      for (x=s.x0-1;x<=s.x1+1;x++) { 
	W(x,s.y0-1,z) = W(x,s.y1-1,z);  
	W(x,s.y1+1,z) = W(x,s.y0+1,z);  
      }
    }
  } else {
    for (z=s.z0;z<=s.z1;z++) { 
      for (x=s.x0-1;x<=s.x1+1;x++) { 
	W(x,s.y0-1,z) = W(x,s.y0+1,z);  
	W(x,s.y1+1,z) = W(x,s.y1-1,z);  
      }
    }
  }

  /* Set fictitious points in z-direction */
  if (pbc_z) {
    for (y=s.y0-1;y<=s.y1+1;y++) { 
      for (x=s.x0-1;x<=s.x1+1;x++) { 
	W(x,y,s.z0-1) = W(x,y,s.z1-1);  
	W(x,y,s.z1+1) = W(x,y,s.z0+1);  
      }
    }
  } else {
    for (y=s.y0-1;y<=s.y1+1;y++) { 
      for (x=s.x0-1;x<=s.x1+1;x++) { 
	W(x,y,s.z0-1) = W(x,y,s.z0+1);  
	W(x,y,s.z1+1) = W(x,y,s.z1-1);  
      }
    }
  }

  if (manypoint) {

/* ---------------------------------
 * Nineteen-point Laplacian formulas
 * --------------------------------- */

    for (y=s.y0;y<=s.y1;y++) {
      for (x=s.x0;x<=s.x1;x++) {
	Sigma_w(s2,x  ,y  ,s.z0) += two*W(x,y,s.z0-1);  
	Sigma_w(s2,x+1,y  ,s.z0) +=     W(x,y,s.z0-1);  
	Sigma_w(s2,x-1,y  ,s.z0) +=     W(x,y,s.z0-1); 
	Sigma_w(s2,x  ,y+1,s.z0) +=     W(x,y,s.z0-1);  
	Sigma_w(s2,x  ,y-1,s.z0) +=     W(x,y,s.z0-1); 
	
	Sigma_w(s2,x  ,y  ,s.z1) += two*W(x,y,s.z1+1); 
	Sigma_w(s2,x+1,y  ,s.z1) +=     W(x,y,s.z1+1);  
	Sigma_w(s2,x-1,y  ,s.z1) +=     W(x,y,s.z1+1); 
	Sigma_w(s2,x  ,y+1,s.z1) +=     W(x,y,s.z1+1);  
	Sigma_w(s2,x  ,y-1,s.z1) +=     W(x,y,s.z1+1); 
      }
    }
    
    for (z=s.z0;z<=s.z1;z++) {
      for (x=s.x0;x<=s.x1;x++) {
	Sigma_w(s2,x  ,s.y0,z  ) += two*W(x,s.y0-1,z);  
	Sigma_w(s2,x+1,s.y0,z  ) +=     W(x,s.y0-1,z);  
	Sigma_w(s2,x-1,s.y0,z  ) +=     W(x,s.y0-1,z); 
	Sigma_w(s2,x  ,s.y0,z+1) +=     W(x,s.y0-1,z);  
	Sigma_w(s2,x  ,s.y0,z-1) +=     W(x,s.y0-1,z); 
	
	Sigma_w(s2,x  ,s.y1,z  ) += two*W(x,s.y1+1,z); 
	Sigma_w(s2,x+1,s.y1,z  ) +=     W(x,s.y1+1,z);  
	Sigma_w(s2,x-1,s.y1,z  ) +=     W(x,s.y1+1,z); 
	Sigma_w(s2,x  ,s.y1,z+1) +=     W(x,s.y1+1,z);  
	Sigma_w(s2,x  ,s.y1,z-1) +=     W(x,s.y1+1,z); 
      }
    } 
    
    for (z=s.z0;z<=s.z1;z++) {
      for (y=s.y0;y<=s.y1;y++) {
	Sigma_w(s2,s.x0,y  ,z  ) += two*W(s.x0-1,y,z); 
	Sigma_w(s2,s.x0,y+1,z  ) +=     W(s.x0-1,y,z);  
	Sigma_w(s2,s.x0,y-1,z  ) +=     W(s.x0-1,y,z); 
	Sigma_w(s2,s.x0,y  ,z+1) +=     W(s.x0-1,y,z);  
	Sigma_w(s2,s.x0,y  ,z-1) +=     W(s.x0-1,y,z); 
	
	Sigma_w(s2,s.x1,y  ,z  ) += two*W(s.x1+1,y,z); 
	Sigma_w(s2,s.x1,y+1,z  ) +=     W(s.x1+1,y,z);  
	Sigma_w(s2,s.x1,y-1,z  ) +=     W(s.x1+1,y,z); 
	Sigma_w(s2,s.x1,y  ,z+1) +=     W(s.x1+1,y,z);  
	Sigma_w(s2,s.x1,y  ,z-1) +=     W(s.x1+1,y,z); 
      }
    }
    
    for (x=s.x0;x<=s.x1;x++) {
      Sigma_w(s2,x,s.y0,s.z0) += W(x,s.y0-1,s.z0-1);  
      Sigma_w(s2,x,s.y1,s.z0) += W(x,s.y1+1,s.z0-1);  
      Sigma_w(s2,x,s.y0,s.z1) += W(x,s.y0-1,s.z1+1);  
      Sigma_w(s2,x,s.y1,s.z1) += W(x,s.y1+1,s.z1+1);  
    }
    
    for (y=s.y0;y<=s.y1;y++) {
      Sigma_w(s2,s.x0,y,s.z0) += W(s.x0-1,y,s.z0-1);  
      Sigma_w(s2,s.x1,y,s.z0) += W(s.x1+1,y,s.z0-1);  
      Sigma_w(s2,s.x0,y,s.z1) += W(s.x0-1,y,s.z1+1);  
      Sigma_w(s2,s.x1,y,s.z1) += W(s.x1+1,y,s.z1+1);  
    }
    
    for (z=s.z0;z<=s.z1;z++) {
      Sigma_w(s2,s.x0,s.y0,z) += W(s.x0-1,s.y0-1,z);  
      Sigma_w(s2,s.x1,s.y0,z) += W(s.x1+1,s.y0-1,z);  
      Sigma_w(s2,s.x0,s.y1,z) += W(s.x0-1,s.y1+1,z);  
      Sigma_w(s2,s.x1,s.y1,z) += W(s.x1+1,s.y1+1,z);  
    }
    
  } else { /* if manypoint */
    
    /* ------------------------------
     * Seven-point Laplacian formulas 
     * ------------------------------ */
    
    for (y=s.y0;y<=s.y1;y++) {
      for (x=s.x0;x<=s.x1;x++) {
	Sigma_w(s2,x,y,s.z0) += W(x,y,s.z0-1);  
	Sigma_w(s2,x,y,s.z1) += W(x,y,s.z1+1); 
      }
    }
    
    for (z=s.z0;z<=s.z1;z++) {
      for (x=s.x0;x<=s.x1;x++) {
	Sigma_w(s2,x,s.y0,z) += W(x,s.y0-1,z);  
	Sigma_w(s2,x,s.y1,z) += W(x,s.y1+1,z); 
      }
    }
    
    for (z=s.z0;z<=s.z1;z++) {
      for (y=s.y0;y<=s.y1;y++) {
	Sigma_w(s2,s.x0,y,z) += W(s.x0-1,y,z); 
	Sigma_w(s2,s.x1,y,z) += W(s.x1+1,y,z); 
      }
    }
    
  } /* if manypoint else */
#undef W
#undef Sigma_w
} /* Impose_boundary_conditions_3D */
/* ========================================================================= */
