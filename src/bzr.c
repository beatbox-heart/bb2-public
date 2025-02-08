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

/* Rovinsky's model of BZ reaction */
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
  real alp, bet, eps, mu, q;	/* reaction constants */
  real Ix, Iz; 			/* external currents */
} rhspar;

RHS_HEAD(bzr,2) {
  DEVICE_CONST(real,alp)
  DEVICE_CONST(real,bet)
  DEVICE_CONST(real,eps)
  DEVICE_CONST(real,mu )
  DEVICE_CONST(real,q  )
  DEVICE_CONST(real,Ix )
  DEVICE_CONST(real,Iz )
  real x = u[0];
  real z = u[1];
  real v;
  if (x<=mu) return 0;
  if (z>=1) return 0;
  v = alp*z/(1-z);
  du[0] = ( x*(1-x) - (2*q*v + bet)*(x-mu)/(x+mu) ) / eps + Ix;
  du[1] = x - v + Iz;
} RHS_TAIL(bzr)

static real sqrbuf;
#define sqr(a) (sqrbuf=a*a,sqrbuf)

RHS_CREATE_HEAD(bzr) {
  real A=0.33; real B=0.167; real C=0.00167; real h0=0.35;
  real k1=17.8; real k4=3023.; real k5=1.78e6;
  real k7=2.67; real k8=3.56e-6; real k13=1.78e-7;
  real alp0=k4*k8*B/sqr(k1*A*h0);
  real bet0=2*k4*k13*B/(sqr(k1*A)*h0);
  real eps0=k1*A/(k4*C);
  real mu0 =k4*k7/(k1*k5);
  real q0  =0.6;
  ACCEPTP(alp,alp0,RSUCC(0.),RNONE);
  ACCEPTP(bet,bet0,RSUCC(0.),RNONE);
  ACCEPTP(eps,eps0,RSUCC(0.),RNONE);
  ACCEPTP(mu ,mu0 ,RSUCC(0.),RNONE);
  ACCEPTP(q  ,q0  ,RSUCC(0.),RNONE);

  ACCEPTP(Ix ,   0,RNONE,RNONE);
  ACCEPTP(Iz ,   0,RNONE,RNONE);
} RHS_CREATE_TAIL(bzr,2)
