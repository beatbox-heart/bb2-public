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

/* Show graph of specified variables */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "k_.h"

#if defined(NOX11)
  #include "nograph.h"
  NOGRAPH_DUMMY(k_draw)
#else

#include "windraw.h"

extern int Verbose;            /* defined in main */

typedef struct {
  BGIWindow wnd;
  #include "k_code.h"
  real absmin, absmax, ordmin, ordmax;
  real absold, ordold;
  double abs, ord, col;
  real lines;
} STR;

RUN_HEAD(k_draw)
{
  DEVICE_CONST(BGIWindow,wnd);
  #include "k_def.h"
  DEVICE_CONST(real,lines)
  DEVICE_VAR(real,abs) DEVICE_VAR(real,absold) DEVICE_CONST(real,absmin) DEVICE_CONST(real,absmax)
  DEVICE_VAR(real,ord) DEVICE_VAR(real,ordold) DEVICE_CONST(real,ordmin) DEVICE_CONST(real,ordmax)
  DEVICE_VAR(real,col)
  int c;
  #include "k_exec.h"
  SetWindow(wnd);
  if NOT(SetLimits(absmin,absmax,ordmin,ordmax)) return 0;
  c=((int)*col)%16;
  if( (*absold!=RNONE) && (*ordold!=RNONE)
  &&  fabsl((*abs-*absold)/(absmax-absmin))<=lines
  &&  fabsl((*ord-*ordold)/(ordmax-ordmin))<=lines )
    Line(*absold,*ordold,*abs,*ord,c);
  else
    Pixel(*abs,*ord,c);
  *absold=*abs; *ordold=*ord;
}
RUN_TAIL(k_draw)

DESTROY_HEAD(k_draw)
{
  #include "k_free.h"
}
DESTROY_TAIL(k_draw)

CREATE_HEAD(k_draw)
{
  ACCEPT_WINDOW(wnd);
  k_on();                                       CHK(NULL);
  memcpy(loctb,deftb,sizeof(*deftb));
  tb_insert_real(loctb,"abs",&(S->abs));	CHK("abs");
  tb_insert_real(loctb,"ord",&(S->ord));	CHK("ord");
  tb_insert_real(loctb,"col",&(S->col));	CHK("col");
  #include "k_comp.h"
  if(!used(S->data,S->ncode,&(S->abs)))EXPECTED_ERROR("/*WARNING: variable \'abs\' never assigned!!*/");
  if(!used(S->data,S->ncode,&(S->ord)))EXPECTED_ERROR("/*WARNING: variable \'ord\' never assigned!!*/");
  if(!used(S->data,S->ncode,&(S->col)))EXPECTED_ERROR("/*WARNING: variable \'col\' never assigned!!*/");
  k_off();					CHK(NULL);
  ACCEPTR(lines,1,0,RNONE);
  ACCEPTR(absmin,0.0,RNONE,RNONE);
  ACCEPTR(absmax,1.0,RNONE,RNONE); ASSERT(S->absmin!=S->absmax);
  ACCEPTR(ordmin,0.0,RNONE,RNONE);
  ACCEPTR(ordmax,1.0,RNONE,RNONE); ASSERT(S->ordmin!=S->ordmax);
  S->absold=S->ordold=RNONE;
}
CREATE_TAIL(k_draw,0)

#endif
