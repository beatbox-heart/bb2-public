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

/* Show graph of specified variables (synchronous) */

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
#include "windraw.h"
#include "k_.h"


#if defined(NOX11)
  #include "nograph.h"
  NOGRAPH_DUMMY(k_plot)
#else

extern int Verbose;            /* defined in main */

typedef struct {
  BGIWindow wnd;
  #include "k_code.h"
  INT i;
  int N;
  real absmin, absmax;
  real ordmin, ordmax;
  REAL abs, ord, col;
  real lines;
  int linewidth;
  int clean;
} STR;

RUN_HEAD(k_plot)
{
  DEVICE_CONST(BGIWindow,wnd);
  #include "k_def.h"
  DEVICE_CONST(int,N) 
  DEVICE_VAR(INT,i)
  DEVICE_VAR(double,abs) DEVICE_CONST(double,absmin) DEVICE_CONST(double,absmax)
  DEVICE_VAR(double,ord) DEVICE_CONST(double,ordmin) DEVICE_CONST(double,ordmax)
  DEVICE_VAR(double,col)
  DEVICE_CONST(real,lines) DEVICE_CONST(int,linewidth) DEVICE_CONST(int,clean)
  real absold, ordold;
  int c;

  if (clean) {SetWindow(wnd); Clean();}
  absold=ordold=RNONE;
  SetWindow(wnd);
  if NOT(SetLimits(absmin,absmax,ordmin,ordmax)) return 0;
  setlinewidth(linewidth);
  k_on();
			    /*printf("%d:\n", N);*/
  for((*i)=0;(*i)<N;(*i)++) {
    *abs=(N>1)?absmin+(REAL)(*i)*(absmax-absmin)/(N-1):0.5*(absmin+absmax); /* default abscissa */
    #include "k_exec.h"
    c=((int)*col)%16;
    if( (absold!=RNONE) && (ordold!=RNONE)
    &&  fabs((*abs-absold)/(absmax-absmin))<=lines
    &&  fabs((*ord-ordold)/(ordmax-ordmin))<=lines )
      Line(absold,ordold,*abs,*ord,c);
    /*else
      Pixel(*abs,*ord,c);*/
    absold=*abs; ordold=*ord;
				/*printf("%d %lf %lf\n", *i, *abs, *ord);*/
  }
				/*printf("\n\n", *abs, *ord);*/
  setlinewidth(1);
}
RUN_TAIL(k_plot)

DESTROY_HEAD(k_plot)
{
  #include "k_free.h"
}
DESTROY_TAIL(k_plot)

CREATE_HEAD(k_plot)
{
  DEVICE_IS_SPACELESS;
  ACCEPT_WINDOW(wnd);
  k_on();				    CHK(NULL);
  memcpy(loctb,deftb,sizeof(*deftb));
  tb_insert_int (loctb,  "i",&(S->i));	    CHK("i");
  tb_insert_real(loctb,"abs",&(S->abs));    CHK("abs");
  tb_insert_real(loctb,"ord",&(S->ord));    CHK("ord");
  tb_insert_real(loctb,"col",&(S->col));    CHK("col");
  #include "k_comp.h"
  if(!used(S->data,S->ncode,&(S->col))) EXPECTED_ERROR("variable \'col\' never assigned");
  if(!used(S->data,S->ncode,&(S->ord))) EXPECTED_ERROR("variable \'ord\' never assigned");
  ACCEPTR(lines,1,0,1);
  ACCEPTI(linewidth,1,1,INONE);
  ACCEPTI(clean,1,0,1);
  ACCEPTI(N,INONE,1,INONE);
  ACCEPTR(absmin,0,RNONE,RNONE);
  ACCEPTR(absmax,S->N-1,RNONE,RNONE); ASSERT(S->absmin!=S->absmax);
  ACCEPTR(ordmin,0,RNONE,RNONE);
  ACCEPTR(ordmax,1,RNONE,RNONE); ASSERT(S->ordmin!=S->ordmax);
  /* FREE(loctb); */
  k_off();
}
CREATE_TAIL(k_plot,0)

#endif
