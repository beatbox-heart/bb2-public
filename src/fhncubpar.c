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

/* Cubic FitzHugh-Nagumo model in Winfree's notation [Chaos1(3):305,1991],
 * with some more parameters added.

Generic AP using cubic FitzHugh-Nagumo model as descibed in:
Winfree AT. "Varieties of spiral wave behavior: An experimentalist's approach to the theory of excitable media.", Chaos. 1991 Oct;1(3):303-334.
This code is an implementation of formulas given in the above paper with minimal modifications.

 */
#include <assert.h>
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
  real
    bet, gam,			/* Winfree's Kinetics A parameters */
    ku, epsu, epsv,		/* may be required different */
    gNa, gK;			/* may be required different */
  real 
    Iu,  Iv; 			/* External influence */
} STR;

RHS_HEAD(fhncubpar,2)
{
  DEVICE_CONST(real,bet)
  DEVICE_CONST(real,gam)
  DEVICE_CONST(real,ku)
  DEVICE_CONST(real,gNa)
  DEVICE_CONST(real,gK)
  DEVICE_CONST(real,epsu)
  DEVICE_CONST(real,epsv)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)
  real U, U3, V;
  U = u[0];
  V = u[1];
  U3=U*U; U3*=U;
  du[0] = ku*(gNa*(U-U3/3.0)-gK*V)/epsu + Iu;
  du[1] = (U+bet-gam*V)*epsv + Iv;
} RHS_TAIL(fhncubpar)

RHS_CREATE_HEAD(fhncubpar)
  ACCEPTP(epsu,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(epsv,RNONE,0.,RNONE);
  ACCEPTP(ku,1.0,0.,RNONE);
  ACCEPTP(gNa,1.0,0.,RNONE);
  ACCEPTP(gK,1.0,0.,RNONE);
  ACCEPTP(bet,RNONE,0.,RNONE);
  ACCEPTP(gam,RNONE,0.,RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
  ACCEPTP(Iv,0,RNONE,RNONE);
RHS_CREATE_TAIL(fhncubpar,2)
