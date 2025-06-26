/**
 * Copyright (C) (2010-2020) Vadim Biktashev, Irina Biktasheva et al. 
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
A serial Beatbox device for the 2D implicit diffusion step
using alternating directions implicit method.
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

#define a(x) (((x)==0)?0.0:((x)==(n-1))?(-2.0*gam):(-gam))
#define b(x) (1+2*gam)
#define c(x) (((x)==0)?(-2*gam):((x)==(n-1))?(0.0):(-gam))


void solveMatrix2d(int n, int ymin, int ymax, int yloc, real *dd, real *x, int incx, int incy, real gam) {
  /*
   * n - number of points in the sweep direction
   * yloc - location for rhs
   * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
   * b - the main diagonal
   * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
   * d - right hand sides
   * x - the answer
   * incx - increment of x and d arrays in the big memory
   */

  int i;
  real m;

  real *cp, *dp;

  cp=(real *)calloc(n,sizeof(real));
  dp=(real *)calloc(n,sizeof(real));

#define top(i)     ((dd[(i)*incx])+gam*(2.0*(dd[(i)*incx-incy])-2.0*(dd[(i)*incx])))
#define middle(i)  ((dd[(i)*incx])+gam*((dd[(i)*incx+incy])+(dd[(i)*incx-incy])-2.0*(dd[(i)*incx])))
#define bottom(i)  ((dd[(i)*incx])+gam*(2.0*(dd[(i)*incx+incy])-2.0*(dd[(i)*incx])))

/* DD is the RHS and computed only from the input, not the output. This will work only for a perfect square. The n here is ny */
#define DD(i,j) (((j)==ymin)?(bottom((i))):((j)==(ymax))?(top((i))):(middle((i))))

// This is the solution.
#define X(i) (x[(i)*incx])


  cp[0]=c(0)/b(0);
  dp[0]=DD(0,yloc)/b(0); 

  for (i = 1; i < n ; i++) {
    m = b(i) - cp[i-1]*a(i);
	if(i<n-1)
    cp[i] = c(i)/m;
    dp[i] = (DD(i,yloc) - dp[i-1]*a(i))/m;
  }

dp[n-1] = (DD(n-1,yloc) - a(n-1)*dp[n-2])/(b(n-1) - a(n-1)*cp[n-2]);

X(n-1) = dp[n-1];
  
  for (i=n-2; i>=0; i--)
    X(i)=dp[i]-cp[i]*X(i+1);


  free(dp);
  free(cp);
}

typedef struct {
  int v0, v1;
  real D;
  real ht;
  real hx;
} STR;

RUN_HEAD(adi2d)
{
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  DEVICE_CONST(real,D)
  DEVICE_CONST(real,ht)
  DEVICE_CONST(real,hx)
  int x, y, z;

  real gam=D*ht/(2.0*hx*hx);
  real *input, *output;

// set up the a, b, c, v, x_sol, and RHS.
  for(z=s.z0;z<=s.z1;z++) {
    for(y=s.y0;y<=s.y1;y++) {
      input  = &(New[ind(s.x0,y,z,v0)]);
      output = &(New[ind(s.x0,y,z,v1)]); // it is not supposed to be populated: this is what comes back from solvematrix2d/1d
      solveMatrix2d((s.x1-s.x0+1), s.y0, s.y1, y, input, output, DX, DY, gam); /* incx=DY for diffusion in y direction etc */
    } /* for y */
  } /* for z */

  for(z=s.z0;z<=s.z1;z++) {
    for(x=s.x0;x<=s.x1;x++) {
      input =  &(New[ind(x,s.y0,z,v1)]);
      output = &(New[ind(x,s.y0,z,v0)]); // it is not supposed to be populated: this is what comes back from solvematrix2d/1d
      solveMatrix2d((s.y1-s.y0+1), s.x0, s.x1, x, input, output, DY, DX, gam); /* incx=DY for diffusion in y direction etc */
    } /* for y */
  } /* for z */
}
RUN_TAIL(adi2d)

DESTROY_HEAD(adi2d) /* no need to destroy anything in adi1d */
DESTROY_TAIL(adi2d)

CREATE_HEAD(adi2d)
{
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT(v1 != v0);
  ACCEPTR(D,RNONE,RNONE,RNONE); /* could be negative if it is cross-diffusion */
  ACCEPTR(hx,RNONE,0.,RNONE);
  ACCEPTR(ht,RNONE,0.,RNONE);
}
CREATE_TAIL(adi2d,1)

