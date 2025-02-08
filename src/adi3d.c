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
A serial Beatbox device for the 3D implicit diffusion step
using alternating directions implicit method.
Neumann boundary condition. 
Brian algorithm from Carnahan et al. book for implicit
solution of heat equation (eqn. 7.56 on page 453). 
*/

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "qpp.h"
#include "rhs.h"
#include "state.h" // for the DX, DY, DZ


static void solveMatrix3d (
  real *x, int n0, int inc0,
  real *dd1, int inc1, 
  real *dd2, int inc2,
  real gam
) {
  /*
   * x - output, the implicitly resolved array
   * n0 - number of points in the sweep direction
   * inc0 - increment in sweep direction

   * dd1 - input across in direction 1, the array which is also present in the LHS
   * inc1 - increment in direction 1

   * dd2 - input across in direction 2, the array which is not present in the LHS
   * inc2 - increment in direction 2

   * gam - diffusion*ht/(2*hx*hx)
   */

  int i, imin, imax;
  real m;

  /* Progonka/double sweep/Thomas algorithm. */
  /* Notations as in http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm */
  /* (as retrieved 2012/03/14) */

  /* cp == $c'$, dp == $d'$ */
  real *cp=(real *)calloc(n0,sizeof(real));
  real *dp=(real *)calloc(n0,sizeof(real));

  /* The tridiagonal matrix coefficients */
  /*  - the formulas for the first, intermediate and last rows are different */
  /* The subdiagonal */
  /* #define amin(i) (0.0) - non-existent */
  #define amid(i) (-gam)
  #define amax(i) (-2.0*gam)
  /* The main diagonal */
  #define bmin(i) (1+2*gam)
  #define bmid(i) (1+2*gam)
  #define bmax(i) (1+2*gam)
  /* The superdiagonal */
  #define cmin(i) (-2*gam)
  #define cmid(i) (-gam)
  /* #define cmax(i) (0.0) - non-existent */

  /* The vector of right-hand sides */
  /* - similar structure in all three cases of the Brian method */
  #define DD(i) ( \
    (dd1[(i)*inc0]) \
    +gam*((dd1[(i)*inc0+inc1])+(dd1[(i)*inc0-inc1])-2.0*(dd1[(i)*inc0])) \
    +gam*((dd2[(i)*inc0+inc2])+(dd2[(i)*inc0-inc2])-2.0*(dd2[(i)*inc0])) \
  )

  /* The solution vector */
  #define X(i) (x[(i)*inc0])

  /* C and maths notations differ; we try and avoid any confusion */
  imin=0; imax=n0-1;

  /* The forward sweep */
  cp[imin]=cmin(imin)/bmin(imin);
  dp[imin]=DD(imin)/bmin(imin);
  /* the formulas at i=imax are different so it is excluded */
  for (i=imin+1; i<imax; i++) {
    m = bmid(i) - cp[i-1]*amid(i);
    cp[i] = cmid(i)/m;
    dp[i] = (DD(i) - dp[i-1]*amid(i))/m;
  }

  /* The back sweep */
  dp[imax] = (DD(imax) - dp[imax-1]*amax(imax))/(bmax(imax) - cp[imax-1]*amax(imax));
  X(imax) = dp[imax];
  for (i=imax-1; i>=imin; i--) {
    X(i)=dp[i]-cp[i]*X(i+1);
  }

  free(dp);
  free(cp);
}

typedef struct {
  real D;
  real ht;
  real hx;
  int v0;
  int v1;
} STR;

RUN_HEAD(adi3d)
{
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  DEVICE_CONST(real,D)
  DEVICE_CONST(real,ht)
  DEVICE_CONST(real,hx)
  real gam=D*ht/(2.0*hx*hx); /* This method used two alternative HALF-steps */
  int x, y, z, n0;

  /* Allocate layers for the starting, intermediate and final values of the method. */
  /* Note that voltage layer is likely to be within a contiguous state vector,      */
  /* so the auxiliary layers for the intermediate steps have to be separate         */

  /* First $v_n$ of (7.56); later $v_{n+1}$ */
  #define V0 (New+ind(x,y,z,v0))
  /* $v^*$ of (7.56) */
  #define V1 (New+ind(x,y,z,v1))
  /* $v^{**}$ of (7.56) */
  #define V2 (New+ind(x,y,z,v1+1))

  /* Eqn 7.56(i): implicit along the x axis, build on $v_n$ = V0 */
  /* First set the halo points, to impose the lateral boundary conditions */
  for (x=s.x0;x<=s.x1;x++) {
    for (z=s.z0;z<=s.z1;z++) {y=s.y0;V0[-DY]=V0[DY]; y=s.y1; V0[DY]=V0[-DY];}
    for (y=s.y0;y<=s.y1;y++) {z=s.z0;V0[-DZ]=V0[DZ]; z=s.z1; V0[DZ]=V0[-DZ];}
  }
  /* Now do the double-sweeps on all lines */
  for(z=s.z0;z<=s.z1;z++) {
    for(y=s.y0;y<=s.y1;y++) {
      x=s.x0; n0=s.x1-s.x0+1;
      solveMatrix3d(V1,n0,DX, V0,DY, V0,DZ, gam);
    } /* for y */
  } /* for z */

  /* Eqn 7.56(ii): implicit along the y axis, build on $v_n$ = V0 */
  for (y=s.y0;y<=s.y1;y++) {
    for (x=s.x0;x<=s.x1;x++) {z=s.z0;V0[-DZ]=V0[DZ]; z=s.z1;V0[DZ]=V0[-DZ];}
    for (z=s.z0;z<=s.z1;z++) {x=s.x0;V1[-DX]=V1[DX]; x=s.x1;V1[DX]=V1[-DX];}
  } /* for y */
  for(z=s.z0;z<=s.z1;z++) {
    for(x=s.x0;x<=s.x1;x++) {
      y=s.y0; n0=s.y1-s.y0+1;
      solveMatrix3d(V2,n0,DY, V0,DZ, V1,DX, gam);
    } /* for x */
  } /* for z */

  /* Eqn 7.56(iii): implicit along the z axis, build on $v^{**}$ = V2 */
  for (z=s.z0;z<=s.z1;z++) {
    for (x=s.x0;x<=s.x1;x++) {y=s.y0;V2[-DY]=V2[DY]; y=s.y1; V2[DY]=V2[-DY];}
    for (y=s.y0;y<=s.y1;y++) {x=s.x0;V1[-DX]=V1[DX]; x=s.x1; V1[DX]=V1[-DX];}
  } /* for z */
  for(x=s.x0;x<=s.x1;x++) {
    for(y=s.y0;y<=s.y1;y++) {
      z=s.z0; n0=s.z1-s.z0+1;
      solveMatrix3d(V0,n0,DZ, V2,DY, V1,DX, gam);
    } /* for y */
  } /* for x */
}
RUN_TAIL(adi3d)

/* At some stage, in might be an idea to allocate cp, dp once; */
/* currently no need to destroy anything in this device */
DESTROY_HEAD(adi3d) 
DESTROY_TAIL(adi3d)

CREATE_HEAD(adi3d)
{
  ACCEPTR(D,RNONE,0.,RNONE);
  ACCEPTR(hx,RNONE,RSUCC(0.),RNONE);
  ACCEPTR(ht,RNONE,0.,RNONE);
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT (v1 != v0);
}
CREATE_TAIL(adi3d,1)
