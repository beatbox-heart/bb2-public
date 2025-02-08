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

/* Poincare cross-section of an input data stream.			*/
/* Can be used e.g. to measure APD or CV.         			*/
/* 									*/
/* An N-dimensional "trajectory" is defined by a			*/
/* list={var1=expr1; var2=expr2; ... varN=exprN} block 			*/
/* One of the vars is used to define the crossing: 			*/
/* it happens when this var crosses zero (up, down or either).		*/
/* This selected variable is "magical": after this device,		*/
/* its value is zero unless a crossing has happened.			*/
/* All other variables will contain the values of their expressions	*/
/* linearly interpolated to the moment of the crossing,			*/
/* or remain unchanged if none has happened 				*/

#include <assert.h>
#include <stdarg.h>
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

typedef struct {
  /* The code in "pgm={}" block includes ncode semicolon-separated statements "var=expr" */
  #include "k_code.h"
  /* User-defined parameters: */
  int which;			/* which variable used to define the cross-section (0..ncode-1) */
  int sign;			/* crossing is when which is increasing(+) decreasing(-) or either(0) */
  char debugname[MAXPATH];	/* filename for debug output, if any */
  /* Internal state variables */
  int begun;			/* already started: previous is defined */
  real *prev, *next;		/* previous/next values of the which variable */
  FILE *debug;			/* debug file descriptor */
} STR;

RUN_HEAD(k_poincare)
  DEVICE_CONST(int,ncode)
  DEVICE_ARRAY(pp_fn,code)
  DEVICE_ARRAY(p_real,data)
  DEVICE_ARRAY(real,prev)
  DEVICE_ARRAY(real,next)
  DEVICE_CONST(int,which)
  DEVICE_CONST(int,sign)
  DEVICE_VAR(int,begun)
  DEVICE_CONST(FILE *,debug)
  int icode, n;
  void *the_code;
  REAL the_result;
  real p, q;

  if (mpi_rank==0 && debug) fprintf(debug, "%ld ", t);
  k_on();

  /* for every assignment in the pgm block */
  for(icode=0;icode<ncode;icode++) {
    prev[icode]=next[icode];			/* remember the previous value of the variable */
    if NOT(the_code=code[icode]) continue;	/* skip if the assignment is empty */
    the_result = *(p_real)execute(the_code); CHK(NULL); /* calc the rvalue */
    next[icode]=the_result;			/* assign the next value of the variable */
    if (mpi_rank==0 && debug) fprintf(debug," %s=" REALF,var_name(deftb,data[icode],n), the_result);
  }
  if (mpi_rank==0 && debug) {fprintf(debug,"\n"); fflush(debug);}

  if(*begun) { /* "previous" is already defined */
    real y0=prev[which], y1=next[which];	/* prev/next values of which var */
    if( ((sign>=0)&&(y0<=0)&&(y1>0)) 		/* crossing in positive direction expected and happened */
	|| ((sign<=0)&&(y0>=0)&&(y1<0)) ) { 	/* crossing in negative direction expected and happened */
      if (mpi_rank==0 && debug) {fprintf(debug,"!!"); fflush(debug);}

      p=y1/(y1-y0); q=-y0/(y1-y0); 		/* linear interpolation coefficients */
		    assert(p>=0); assert(p<=1);
		    assert(q>=0); assert(q<=1);

      for(icode=0;icode<ncode;icode++)  {
        *(data[icode])=
	  (icode!=which)?p*prev[icode]+q*next[icode]:  /* interpolate all vars except */
	  (y1>y0)?1:-1;					/* the which var which is given the sign of the crossing */
        if (mpi_rank==0 && debug) fprintf(debug," %s=" REALF,var_name(deftb,data[icode],n), *(data[icode]));
      }

      if (mpi_rank==0 && debug) {fprintf(debug,"\n"); fflush(debug);}
    } else {
      *(data[which])=0;				/* there were no crossing */
      if (mpi_rank==0 && debug) {
	fprintf(debug," no crossing:"); fflush(debug);
	for(icode=0;icode<ncode;icode++) {
	  fprintf(debug," %s=" REALF,var_name(deftb,data[icode],n), *(data[icode]));
	  fflush(debug);
	}
	fprintf(debug,"\n"); fflush(debug);
      }
    }
  }

  k_off();
  *begun=1;
RUN_TAIL(k_poincare)

DESTROY_HEAD(k_poincare)
  #include "k_free.h"
  FREE(S->prev);
  FREE(S->next);
  SAFE_CLOSE(S->debug);
DESTROY_TAIL(k_poincare)

CREATE_HEAD(k_poincare)
{
  DEVICE_IS_SPACELESS;
  /* DEVICE_MUST_BE_NOWHERE */
  int idata;
  ACCEPTF(debug,"wt","");

  k_on();                                          		CHK(NULL);
  memcpy(loctb,deftb,sizeof(*deftb));
  #include "k_comp.h"
  if (!S->ncode) EXPECTED_ERROR("no statements in \"%s\"",buf);
  /*if (S->ncode==1) EXPECTED_ERROR("must be at least two statements");*/
  k_off();

  ACCEPTI(which,0,0,S->ncode-1);
  ACCEPTI(sign,1,-1,1);
  S->begun=0;

  CALLOC(S->prev,S->ncode,sizeof(real));
  CALLOC(S->next,S->ncode,sizeof(real));

}
CREATE_TAIL(k_poincare,0)

