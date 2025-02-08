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

/* Aliev & Panfilov 2-var model of cardiac AP
 * CSF 7(3):293-301, 1996
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
  real a, eps0, k, mu1, mu2;	/* reaction constants */
  real Iu, Iv; 			/* external currents */
} rhspar;

RHS_HEAD(ap,2) {
  DEVICE_CONST(real,a)	
  DEVICE_CONST(real,eps0)
  DEVICE_CONST(real,k)
  DEVICE_CONST(real,mu1)
  DEVICE_CONST(real,mu2)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)
  real U = u[0];
  real V = u[1];
  if(U<-mu2) ABORT("U=" REALF "<mu2=" REALF, U, mu2);
  du[0] = -U*(V+k*(U-a)*(U-1)) + Iu;
  du[1] = -(eps0+mu1*V/(U+mu2))*(V+k*U*(U-a-1)) + Iv;
} RHS_TAIL(ap)

RHS_CREATE_HEAD(ap) {
  ACCEPTP(a,	0.15,	RSUCC(0.),RNONE);
  ACCEPTP(eps0,	0.002,	RSUCC(0.),RNONE);
  ACCEPTP(k,	8.0,	RSUCC(0.),RNONE);
  ACCEPTP(mu1,  0.2,	RSUCC(0.),RNONE);
  ACCEPTP(mu2,  0.3,	RSUCC(0.),RNONE);
  ACCEPTP(Iu ,   0,RNONE,RNONE);
  ACCEPTP(Iv ,   0,RNONE,RNONE);
} RHS_CREATE_TAIL(ap,2)
