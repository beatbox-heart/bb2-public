/**
 * Copyright (C) (2010-2023) Vadim Biktashev, Irina Biktasheva et al. 
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

/* Generic device executing k-codes for global var data */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"

#include "k_.h"

/*++++++++++++++++++++++++++++++++++++++++*/ extern int idev;

#undef SEPARATORS
#define SEPARATORS ";"
#define BLANK " \n\t\r"

typedef struct {
  int advance;		/* this instance runs first when created, before time starts */
#include "k_code.h"	/* typical k-program components */
  p_tb loctb;		/* local k_table */
  NFILE(debug);		/* if and where to print debug messages */
  int debugWriter;	/* this thread is entrusted with common debug messages */
} STR;

RUN_HEAD(globalfunc)
{
  DEVICE_CONST(int, debugWriter)
  DEVICE_CONST(FILE *,debug)
  DEVICE_CONST(p_tb, loctb)
  int icode, n;
  /* DEBUG("\n"); */

  if (debug && debugWriter) {
    fprintf(debug, " t=%ld", (long)t);
    fflush(debug);
  }
  k_on();
  for(icode=0;icode<S->ncode;icode++) {
    void *code, *result;
    /* DEBUG("icode=%d\n",icode); */
    if NOT(code=(S->code)[icode]) continue;
    result = execute(code); CHK(NULL);
    if (debug && debugWriter) {
      fprintf(debug," %s=%s",var_name(loctb,(S->data)[icode],n), prt(result, res_type(code)));
      fflush(debug);
    }
    memcpy((S->data)[icode],result,sizetable[res_type(code)]);
  } /* for icode */
  k_off();
  if (debug && debugWriter) {
    fprintf(debug,"\n");
    fflush(debug);
  }
  /* DEBUG("\n"); */
}
RUN_TAIL(globalfunc)

DESTROY_HEAD(globalfunc)
{
  #include "k_free.h"
  FREE(S->loctb);
  SAFE_CLOSE(S->debug);
} 
DESTROY_TAIL(globalfunc)

CREATE_HEAD(globalfunc)
{
  /* DEVICE_IS_SPACELESS; */
  DEVICE_MUST_BE_NOWHERE;
  #define DELIM "\\/,; \r\t\n"
  char *p, name[80];
  int iv, ip, nv, np;
  real *array, *pv;

  ACCEPTI(advance,0,0,1);

#if MPI
  /* Identifying the root process allows only one process to write debug output. */
  if (dev->s.nowhere) {
    S->debugWriter = (mpi_rank == 0);
  } else {
    S->debugWriter = (mpi_rank == getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0));
  }
#else
  S->debugWriter = 1;
#endif
  ACCEPTF(debug,"wt","");
  
  k_on();				CHK(NULL);
  S->loctb = tb_new();
  memcpy(S->loctb,deftb,sizeof(*deftb));
  #define loctb (S->loctb)
  #include "k_comp.h"
  #undef loctb
  if (!S->ncode)
    EXPECTED_ERROR("no statements in \"%s\"",buf);
  k_off();

  #if MPI
  if (advance) run_globalfunc (dev->s,S,dev->sync,dev->alwaysRun);
  #else
  if (advance) run_globalfunc (dev->s,S);
  #endif
}
CREATE_TAIL(globalfunc,1)

