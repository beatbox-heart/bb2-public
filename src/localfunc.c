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

/* Generic device executing k-codes for grid data, 
  including "phase" inital conditions */

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
  real *u;		/* [vmax] k-vars */
  real *geom;		/* [geom_vmax] k-vars */
  real x, y, z;		/* k-vars promised to be real in the manual */
  INT nr, nc;		/* size of the table in file */
  real *a;		/* [nr*nc] table from file */
  double *pv;		/* [nc] k-vars assigned by phasep magic */
  double phaseu, phasep; /* magic k-vars */
  NFILE(file);		/* file with the table for magic k-vars */
  NFILE(debug);		/* if and where to print debug messages */
  int debugWriter;	/* this thread is entrusted with common debug messages */
} STR;

RUN_HEAD(localfunc)
{
  DEVICE_CONST(INT, nr);
  DEVICE_CONST(INT, nc);
  DEVICE_ARRAY(real, a);
  DEVICE_ARRAY(double, pv);
  DEVICE_VAR(double,phasep);
  DEVICE_VAR(double,phaseu);
  DEVICE_CONST(int, debugWriter);
  DEVICE_CONST(FILE *, debug);
  DEVICE_CONST(p_tb, loctb);
  SPACE_CONST(int, restricted);
  int ic, ir, ir1;	/* file table counters: column, row, next row */
  double pp, qq;	/* file table linear approximation coeffs */
  int icode;		/* k-code counter */
  int n;		/* aux var in DO_FOR_ALL_POINTS macro */

  assert(DV==1); /* need all layers of the same point be contiguous */

  if (debug && debugWriter) { /* only one thread to report this */
    fprintf(debug, " t=%ld", (long)t);
    fflush(debug);
  }
  k_on();

  if (a) { /* there is a file table and magic k-vars */
  /*********************************************************************/
#define COMMANDS 							\
    if (debug) { /* reported for every point */				\
      fprintf(debug," (%ld %ld %ld)",*x,*y,*z);				\
      fflush(debug); 							\
    } 									\
    									\
    S->x=*x; S->y=*y; S->z=*z; /* R/O k-vars to match current coords */	\
    memcpy(S->u,New+ind((*x),(*y),(*z),0),vmax*sizeof(real)); /* R/W k-vars */ \
    if (geom_vmax)							\
      memcpy(S->geom,Geom+geom_ind((*x),(*y),(*z),0),geom_vmax*sizeof(real)); /* R/W k-vars */ \
    for (icode=0;icode<S->ncode;icode++) { 				\
      void *code, *result; 						\
      if NOT(code=(S->code)[icode]) continue; 				\
      result = execute(code); CHK(NULL);				\
      if (debug) { /* reported for every point */ 			\
	fprintf(debug," %s=%s",var_name(loctb,(S->data)[icode],n), prt(result, res_type(code))); \
	fflush(debug); 							\
      }									\
      memcpy((S->data)[icode],result,sizetable[res_type(code)]);	\
      									\
      if (((S->data)[icode])==phasep) { /* the magic of 'phasep' k-var */ \
	if (nr>1) { 							\
	  (*phasep) *= nr*M_1_PI*0.5;					\
	  while ((*phasep)<  0) (*phasep) += nr;			\
	  while ((*phasep)>=nr) (*phasep) -= nr; 			\
	  ir=floor(*phasep); ir1=(ir+1)%nr; pp=(*phasep-ir); qq=1-pp; 	\
	  assert(ir<nr);						\
	  for (ic=0;ic<nc;ic++) pv[ic]=a[ic+nc*ir]*qq+a[ic+nc*ir1]*pp;	\
	} else {							\
	  for (ic=0;ic<nc;ic++) pv[ic]=a[ic+nc*ir];			\
	} /* if nr else */ 						\
      } /* if phasep */							\
      									\
      if (((S->data)[icode])==phaseu) { /* the magic of 'phaseu' k-var */ \
	if (nr>1) { 							\
	  (*phaseu) *= nr*M_1_PI*0.5;					\
	  while ((*phaseu)<  0) (*phaseu) += nr; 			\
	  while ((*phaseu)>=nr) (*phaseu) -= nr; 			\
	  ir=floor(*phaseu); ir1=(ir+1)%nr; pp=(*phaseu)-ir; qq=1-pp; 	\
	  assert(ir<nr);						\
	  for (ic=0;ic<nc;ic++) S->u[ic]=a[ic+nc*ir]*qq+a[ic+nc*ir1]*pp; \
	} else {							\
	  for (ic=0;ic<nc;ic++) S->u[ic]=a[ic+nc*ir];			\
	} /* if nr else */ 						\
      } /* if phaseu */							\
    } /* for icode */ 							\
    									\
    memcpy(New+ind((*x),(*y),(*z),0),S->u,vmax*sizeof(real)); /* R/W k-vars */ \
    if (geom_vmax)							\
      memcpy(Geom+geom_ind((*x),(*y),(*z),0),S->geom,geom_vmax*sizeof(real)); /* R/W k-vars */   
    DO_FOR_ALL_POINTS;
#undef COMMANDS
  } else { /* not a, no "magic" */
  /*********************************************************************/
#define COMMANDS 							\
    if (debug) { /* reported for every point */				\
      fprintf(debug," (%ld %ld %ld)",*x,*y,*z);				\
      fflush(debug); 							\
    } 									\
    S->x=*x; S->y=*y; S->z=*z; /* R/O k-vars to match current coords */	\
    memcpy(S->u,New+ind((*x),(*y),(*z),0),vmax*sizeof(real)); /* R/W k-vars */ \
    if (geom_vmax)							\
      memcpy(S->geom,Geom+geom_ind((*x),(*y),(*z),0),geom_vmax*sizeof(real)); /* R/W k-vars */ \
    for (icode=0;icode<S->ncode;icode++) { 				\
      void *code, *result; 						\
      if NOT(code=(S->code)[icode]) continue; 				\
      result = execute(code); CHK(NULL);				\
      if (debug) { /* reported for every point */ 			\
	fprintf(debug," %s=%s",var_name(loctb,(S->data)[icode],n), prt(result, res_type(code))); \
	fflush(debug); 							\
      }									\
      memcpy((S->data)[icode],result,sizetable[res_type(code)]);	\
    } /* for icode */ 							\
    									\
    memcpy(New+ind((*x),(*y),(*z),0),S->u,vmax*sizeof(real)); /* R/W k-vars */ \
    if (geom_vmax)							\
      memcpy(Geom+geom_ind((*x),(*y),(*z),0),S->geom,geom_vmax*sizeof(real)); /* R/W k-vars */   
    DO_FOR_ALL_POINTS;
#undef COMMANDS
  } /* not a */
  /*********************************************************************/

  k_off();
  if (debug && debugWriter) { /* only one thread to report this */
    fprintf(debug,"\n");
    fflush(debug);
  }
}
RUN_TAIL(localfunc)

DESTROY_HEAD(localfunc)
{
  #include "k_free.h"
  FREE(S->u);
  FREE(S->geom);
  FREE(S->loctb);
  FREE(S->a);
  FREE(S->pv);
  SAFE_CLOSE(S->file); /* actually it must have been closed after read from */
  SAFE_CLOSE(S->debug);
} 
DESTROY_TAIL(localfunc)

CREATE_HEAD(localfunc)
{
  #define DELIM "\\/,; \r\t\n"
  char *p;		/* char pointer when parsing strings from table file */
  char name[80];	/* k-var name */
  int ic, ir;		/* counters of file table */
  int nc, nr;		/* size of file table */
  real *array;		/* file table */
  real *pv;		/* 'p*' k-vars */
  int iv;		/* layer counter */

  ACCEPTI(advance,0,0,1);
  CALLOC(S->u,vmax,sizeof(real));
  if (geom_vmax) CALLOC(S->geom,geom_vmax,sizeof(real));
	
#if MPI
  /* Identifying the one process to write the common debug output. */
  S->debugWriter = (mpi_rank == getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0));
#else
  S->debugWriter = 1;
#endif
	
  ACCEPTF(file,"rt","");
  if (S->file) {
    fgets(buf,MAXSTRLEN,S->file);		/* Count words in the first line */
    if NOT(p=strtok(buf,DELIM)) EXPECTED_ERROR("first line of %s seems empty, cannot recognize",S->filename);
    for(nc=1;NULL!=(p=strtok(NULL,DELIM));nc++);
    rewind(S->file);				/* Count lines in the file */
    for(nr=0;fgets(buf,MAXSTRLEN,S->file);nr++);
    switch(nr) {
    case 0: EXPECTED_ERROR("no lines in %s",S->filename);
    case 1: MESSAGE("/* WARNING: one line found in %s, phase values will be ingored",S->filename);
    }
    rewind(S->file);				/* read the array from file */
    MESSAGE("\n/* Table of nc=%d x nr=%d is assumed in %s */\n",nc,nr,S->filename);
    CALLOC(array,nc*nr,sizeof(real));
    CALLOC(pv,nc,sizeof(real));
    for (ir=0;ir<nr;ir++) {
      if (ir>0) {				/* prefill the row */
	for(ic=0;ic<nc;ic++) {			/* with previous values */
	  array[ic+ir*nc]=array[ic+(ir-1)*nc];	/* in case this one is */
	}					/* incomplete in the file */
      }
      fgets(buf,MAXSTRLEN,S->file);
      p = strtok(buf,DELIM);
      for (ic=0;ic<nc;ic++) {
	if NOT(p) break;
	sscanf(p,REALF,array+ic+ir*nc);
	p = strtok(NULL,DELIM);
      } /*  for ic */
    } /*  for ir */
    fclose(S->file);
    S->a = array;
    S->nr = nr;
    S->nc = nc;
    S->pv = pv;
  } else {
    S->a=S->pv=NULL;
    S->nr=S->nc=0;
  }
  k_on();				CHK(NULL);
  S->loctb = tb_new();
  memcpy(S->loctb,deftb,sizeof(*deftb));
  tb_insert_real_ro(S->loctb,"x",&(S->x));	CHK("x");
  tb_insert_real_ro(S->loctb,"y",&(S->y));	CHK("y");
  tb_insert_real_ro(S->loctb,"z",&(S->z));	CHK("z");

  /* Read-only k-variables nr and nc count dimensions of the array */
  /* contained in the file; if there is no file, there is no nr nor nc */
  if (S->file) {
    tb_insert_int_ro(S->loctb,"np",&(S->nr)); CHK("np"); /* old name for backward compatibility */
    tb_insert_int_ro(S->loctb,"nr",&(S->nr)); CHK("nr");
    tb_insert_int_ro(S->loctb,"nv",&(S->nc)); CHK("nc"); /* old name for backward compatibility */
    tb_insert_int_ro(S->loctb,"nc",&(S->nc)); CHK("nc");
  }
	
  for(iv=0;iv<vmax;iv++) {
    snprintf(name,80,"u%d",iv);
    tb_insert_real(S->loctb,name,&(S->u[iv]));  CHK(name);
  } /* for iv */
  if (geom_vmax) {
    for(iv=0;iv<geom_vmax;iv++) {
      snprintf(name,80,"geom%d",iv);
      tb_insert_real(S->loctb,name,&(S->geom[iv]));  CHK(name);
    } /* for iv */
  } /* if geom_vmax */
	
  if (S->file) {
    tb_insert_real(S->loctb,"phasep",&(S->phasep));	CHK("phasep");
    tb_insert_real(S->loctb,"phaseu",&(S->phaseu));	CHK("phaseu");
    for(ic=0;ic<nc;ic++) {
      snprintf(name,80,"p%d",ic);
      tb_insert_real(S->loctb,name, &(S->pv[ic]));	CHK(name);
    } /* for ic */
  } /* if file */
  
  #define loctb (S->loctb)
  #include "k_comp.h"
  #undef loctb
  
  if (!S->ncode)
    EXPECTED_ERROR("no statements in \"%s\"",buf);
  if (S->file) {
    if (!used(S->data,S->ncode,&(S->phaseu)) && !used(S->data,S->ncode,&(S->phasep))) {
      MESSAGE("/*WARNING: file was specified but phaseu, phasep not assigned in the expression!!*/");
    }
  }
  k_off();

  #if MPI
  if (advance) run_localfunc (dev->s,S,dev->sync,dev->alwaysRun);
  #else
  if (advance) run_localfunc (dev->s,S);
  #endif
 
  ACCEPTF(debug,"wt","");
}
CREATE_TAIL(localfunc,1)
