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
/* The multigrid solver, shared with other devices */
#include "multigrid.h"

/*****************************************************************************/
/* The device parameters structure */
typedef struct {                
  /* The parameters of the connection weights, as in diff */
#include "diffweights_str.h"
  /* The multigrid solver, incl. layers v0 (RHS) and v1 (solution) */
#include "multigrid_str.h"
} STR;
/*****************************************************************************/

/* The weights code of the diff device: diff_weights_accept(), diff_connect() */
#include "diffweights_code.h"
/* The multigrid code: mg_accept(), mg_connect(), mg_build(), mg_solve(), mg_destroy() */
#include "multigrid_code.h"

/******************/
RUN_HEAD(elliptic)
{
  if NOT(mg_solve(S)) return 0;
}
RUN_TAIL(elliptic)

/*******************/
DESTROY_HEAD(elliptic) {
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
} DESTROY_TAIL(elliptic)

/*******************/
CREATE_HEAD(elliptic)
{
  DEVICE_REQUIRES_SYNC;
  if NOT(mg_accept(dev,w,S,0)) return 0;	/* isolated=0 by default */
  /* ACCEPTI(vp,-1,-1,vmax-1); */
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 != v0 );
  /* ASSERT( v0 != vp ); */
  /* ASSERT( v1 != vp ); */
  if NOT(mg_connect(dev,w,S)) return 0;
  if NOT(mg_build(dev,w,S)) return 0;
}
CREATE_TAIL(elliptic,1)
