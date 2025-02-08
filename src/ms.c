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
  Generic AP cell model as described in:
  C C Mitchell and D G Schaeffer, 
  "A Two-Current Model for the Dynamics of Cardiac Membrane"
  Bull Math Biol 65: 767-793 (2003). 
  This code is based on formulas given in the above paper.
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
  real tau_in;
  real tau_out;
  real J_stim;
  real v_gate;
  real tau_open;
  real tau_close;
} STR;

RHS_HEAD(ms,2)
{
  DEVICE_CONST(real,tau_in);
  DEVICE_CONST(real,tau_out);
  DEVICE_CONST(real,J_stim);
  DEVICE_CONST(real,v_gate);
  DEVICE_CONST(real,tau_open);
  DEVICE_CONST(real,tau_close);
  real v = u[0];
  real h = u[1];
  #define C(v) ((v)*(v)*(1-(v)))
  real J_in = h*C(v)/tau_in;
  real J_out = -v/tau_out;
  real dvdt = J_in + J_out + J_stim;
  real dhdt = (v<v_gate)?((1.0-h)/tau_open):(-h/tau_close);
  du[0] = dvdt;
  du[1] = dhdt;
}
RHS_TAIL(ms)

RHS_CREATE_HEAD(ms)
{
  /* Default parameters as in fig 1 of the paper */
  ACCEPTP(tau_in,   0.3,RSUCC(0),RNONE);
  ACCEPTP(tau_out,    6,RSUCC(0),RNONE);
  ACCEPTP(J_stim, RNONE,RNONE,   RNONE);
  ACCEPTP(v_gate,  0.13,RSUCC(0),RPRED(1));
  ACCEPTP(tau_open, 120,RSUCC(0),RNONE);
  ACCEPTP(tau_close,150,RSUCC(0),RNONE);
  MALLOC(*u,2*sizeof(real));
  (*u)[0]=0.0;
  (*u)[1]=1.0;
}
RHS_CREATE_TAIL(ms,2)
