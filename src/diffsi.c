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

/*
 * This device makes a (semi-)implicit step of diffusion, the theta-scheme
 * \[ u^{n+1} - \theta h_t L u^{n+1} = u^n + (1-\theta) h_t L u^n, \quad L u = \nabla(\hat{D}\nabla u), \]
 * with theta=1 (backward Euler, the default) or theta=1/2 (Crank-Nicolson).
 * The operator L is discretized exactly as in the diff device (the same parameters
 * and code, see diffweights_code.h), and the linear system is solved by the multigrid
 * of the elliptic device (the same parameters and code, see multigrid_code.h).
 *
 * As in diff, v0 is the diffusive field, here updated in place, and v1 an auxiliary
 * layer, which receives the right-hand side of the system. NB in the terms of the
 * multigrid code, S->v1 is the solution and S->v0 the RHS, hence the swap below.
 */

#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "qpp.h"
#include "geometry.h"

/* The connection weights are computed by the code of the diff device */
#include "diffweights.h"
/* The multigrid solver of the elliptic device */
#include "multigrid.h"

/*****************************************************************************/
/* The device parameters structure */
typedef struct {                
  /* The parameters of the connection weights, as in diff */
#include "diffweights_str.h"
  /* The multigrid solver; its v1 is the diffusive field and v0 the RHS */
#include "multigrid_str.h"
  real ht;				/* time step */
  real theta;				/* implicitness: 1 backward Euler, 1/2 Crank-Nicolson */
} STR;
/*****************************************************************************/

/* The weights code of the diff device: diff_weights_accept(), diff_connect() */
#include "diffweights_code.h"
/* The multigrid code: mg_accept(), mg_connect(), mg_build(), mg_solve(), mg_destroy() */
#include "multigrid_code.h"

/******************/
RUN_HEAD(diffsi)
{
  DEVICE_CONST(int,v0)			/* the RHS layer (user's v1) */
  DEVICE_CONST(int,v1)			/* the diffusive field (user's v0) */
  DEVICE_CONST(real,theta)
  DEVICE_CONST(list_t,numPoints)
  DEVICE_ARRAY(AnisoPoint,points)
  real a=(1.0-theta)/theta;		/* explicit to implicit ratio */
  sindex_t off[maxnb];			/* neighbours' displacements in New */
  neighbour_t nb;
  list_t i;
  AnisoPoint *P;
  real *u, sum;

  for (nb=0;nb<maxnb;nb++) off[nb]=ndx[nb]*(sindex_t)DX+ndy[nb]*(sindex_t)DY+ndz[nb]*(sindex_t)DZ;

  /* The RHS u + (1-theta)*ht*L u = u + a*(u - (I-tau*L) u), tau=theta*ht, with the  */
  /* weights K of I-tau*L, incl. the links to the tissue beyond the device space.   */
  /* The points without connections are not in the list, and are left unchanged.   */
  for (i=1;i<=numPoints;i++) {
    P=&(points[i]);
    u=New+ind(P->x,P->y,P->z,v1);
    if (a) {
      for (sum=0,nb=0;nb<maxnb;nb++)	/* K has maxnb entries, unused ones 0 */
	if (P->K[nb]) sum+=u[off[nb]]*P->K[nb];
      New[ind(P->x,P->y,P->z,v0)]=u[0]+a*(u[0]-sum);
    } else {
      New[ind(P->x,P->y,P->z,v0)]=u[0];
    }
  }

  /* Solve (I-tau*L) u = RHS, with the current u as the initial guess */
  if NOT(mg_solve(S)) return 0;
}
RUN_TAIL(diffsi)

/*******************/
DESTROY_HEAD(diffsi) {
  mg_destroy(S);
  /* the parameters of the weights, see diffweights_code.h */
  FREE(S->Xpool);
  FREE(S->Ypool);
  FREE(S->u);
  FREE(S->geom);
  FREE(S->loctb);
  FREE(S->Dsrc);
  FREE(S->Dlsrc);
  FREE(S->Dtsrc);
} DESTROY_TAIL(diffsi)

/*******************/
CREATE_HEAD(diffsi)
{
  DEVICE_REQUIRES_SYNC;
  ACCEPTR(ht,RNONE,RSUCC(0.),RNONE);
  ACCEPTR(theta,1.0,RSUCC(0.),1.0);
  if NOT(mg_accept(dev,w,S,-1)) return 0;	/* isolated by default as in diff */
  /* The layers as in diff, swapped into the multigrid's terms */
  if (!accepti("v0=",&(S->v1),INONE,0,(int)vmax-1,w)) return 0; /* the diffusive field */
  if (!accepti("v1=",&(S->v0),INONE,0,(int)vmax-1,w)) return 0; /* auxiliary: the RHS */
  if (S->v0==S->v1) EXPECTED_ERROR("v0 and v1 must be different layers\n");
  /* The operator of the system is I-tau*L */
  S->tau=theta*ht;
  if NOT(mg_connect(dev,w,S)) return 0;
  if NOT(mg_build(dev,w,S)) return 0;
}
CREATE_TAIL(diffsi,1)
