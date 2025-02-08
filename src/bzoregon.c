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
 * `Oregonator' model of BZ reaction, as in
 * J.P.~Keener \& J.J.~Tyson `Spiral waves in the Belousov-Zhabotinskii reaction
 * Physica D 21:307-324, 1986 
 * and
 * W.~Jahnke \& A.~Winfree `A survey of spiral-wave behaviors in the 
 * Oregonator model' IJBC 1(2):445-466, 1991 (p.446)
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

typedef struct {
  real eps, f, q;		/* reaction constants */
  real Iu, Iv; 			/* external currents */
} rhspar;

RHS_HEAD(bzoregon,2) {
  DEVICE_CONST(real,eps)
  DEVICE_CONST(real,f  )
  DEVICE_CONST(real,q  )
  DEVICE_CONST(real,Iu )
  DEVICE_CONST(real,Iv )
  real U = u[0];
  real V = u[1];

  //  if (U<q) U=q;			/* gimmick blessed by J&W p. 447 */
  du[0] = ( U - U*U - f*V*(U-q)/(U+q) ) / eps + Iu; /* K&T p.309, J&W p.446 (1) */
  du[1] = U - V + Iv;				    
} RHS_TAIL(bzoregon)

RHS_CREATE_HEAD(bzoregon) {
  real eps0=1.e-2, f0=3, q0=2.e-4;			/* Keener \& Tyson, p. 309, last eqns */
  ACCEPTP(eps,eps0,RSUCC(0.),RNONE);
  ACCEPTP(f  ,f0  ,RSUCC(0.),RNONE);
  ACCEPTP(q  ,q0  ,RSUCC(0.),RNONE);

  ACCEPTP(Iu ,   0,RNONE,RNONE);
  ACCEPTP(Iv ,   0,RNONE,RNONE);
  
  MALLOC(*u,2*sizeof(real));
  {
    real f=S->f; 
    real q=S->q;
    real b=f+q-1;
    (*u)[0]=(*u)[1]=0.5*(sqrt(b*b+4*q*(1+f))-b);
  }
  
} RHS_CREATE_TAIL(bzoregon,2)
