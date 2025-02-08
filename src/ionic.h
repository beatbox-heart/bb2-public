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

/**
 * Interface with cardiac cell model description,
 * describing HH-type gates separately. 
 * TS 2016: modified to include Markov Chain models. 
 */

#ifndef _ionic
#define _ionic

typedef struct {		/* description of dependent parameters */
  int n;
  int *src;
  real **dst;
} Var;

#define IONICFTAB(name) int name(real V, real *values, int ntab)
typedef IONICFTAB(IonicFtab);

#define IONIC_FTAB_HEAD(name)	\
IONICFTAB(ftab_##name) {	\
  ASSERT(ntab==NTAB);
#define IONIC_FTAB_TAIL		\
  return 1;             	\
}

/* Type of function definining the right-hand side for non-gate variables */
/* u: vector of dynamic variables */
/* nv: number of elements in v */
/* values: table of tabulated functions */
/* ntab: number of rows in the table (number of voltage values) */
/* Par: the array of parameters of this ionic model */
/* Var: the array of descriptors of variable parameters of this ionic model */
/* du: vector for the derivatives of dynamic variables (output) */
/* no: the number of "other" variables */
/* nalp: the array of non-tabulated alpha-rates */
/* nbet: the array of non-tabulated beta-rates */
/* nn: the number of non-tabulated gates */

#include "channel.h"

#define IONICFDDT(name) int name(real *u,int nv,real *values,int ntab,Par par,Var var,real *du,int no,real *nalp,real *nbet,int nn)
typedef IONICFDDT(IonicFddt);

/* Header of a standard ionic rhs calculator:	*/
/* - nostrify the parameters list,		*/
/* - check dimensions of subvectors,		*/
/* - implement parameter substitution.		*/
#define IONIC_FDDT_HEAD(name,NV,NTAB,NO,NN)	\
IONICFDDT(fddt_##name) {	\
  STR *S = (STR *)par;		\
  int ivar;			\
  ASSERT(nv==NV);		\
  ASSERT(ntab==NTAB);		\
  ASSERT(no==NO);		\
  ASSERT(nn==NN);		\
  if (var.n) for(ivar=0;ivar<var.n;ivar++) *(var.dst[ivar])=u[var.src[ivar]];

#define IONIC_FDDT_TAIL		\
  return 1;             	\
}

/* Solver-independent entities exported by an ionic model description */
typedef struct {
  IonicFtab *ftab;		/* voltage dependent functions that can be tabulated */
  int nmc;			/* number of Markov chain models */
  int nmv;			/* total number of Markov chain variables */
  IonicFddt *fddt;		/* right-hand sides of non-gate equations */
  int no;			/* number of non-gate variables */
  int nn;			/* number of nontab gate variables */
  int nt;			/* number of tab gate variables */
  int ntab;			/* number of tabulated functions (number of columns in the table) */
  int V_index;			/* index of voltage in the state vector */
  Par p;			/* vector of model parameters */
  Var var;			/* description of dependent parameters */
  channel_str * channel;	/* definitions of ion channel */
} ionic_str;

#define IONICCREATE(name) int name(ionic_str *I,char *w,real **u,int v0)
typedef IONICCREATE(IonicCreate);

#define IONIC_CREATE_HEAD(name)				\
IONICCREATE(create_##name) {				\
  STR *S = (STR *)Calloc(1,sizeof(STR));        	\
  char *p=w;                                    	\
  Var *var=&(I->var);					\
  int ivar=0;						\
  int ig; 	/* gates counter */			\
  real V0;	/* initial voltage */			\
  if (!S) ABORT("cannot create %s",#name);		\
  for(var->n=0;*p;var->n+=(*(p++)==AT));		\
  if(var->n){CALLOC(var->dst,var->n,sizeof(real *));	\
  CALLOC(var->src,var->n,sizeof(int));}			\
  else {var->src=NULL;(var->dst)=NULL;}			\
  I->V_index=V_index;					\
  /* Accept the non-cell parameter values by the standard macro */	\
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */\
  /* Create the vector of initial values of dynamic variables. */	\
  MALLOC(*u,(long int)NV*sizeof(real));			\
  int ii;	 /* TODO: rewrite */			\
  /* for (ii=0; ii < NMC; ii++) {nmv+=nm[ii]; nme+=nm[ii]*nm[ii];} */	\
  /* allocate the channel structure */ 			\
  channel_str * ch;					\
  subchain_str * sbch;					\
  CALLOC(I->channel,NMC,sizeof(channel_str));		\
  ch = &(I->channel[0]);				\
  for (ii=0; ii <NMC; ii++)				\
    CALLOC(I->channel[ii].subchain, MAX_SUBCHAINS, sizeof(subchain_str));

#define IONIC_CREATE_TAIL		\
  var->n=ivar;				\
  if(ivar){REALLOC(var->dst,1L*ivar*sizeof(real *));	\
  REALLOC(var->src,1L*ivar*sizeof(int));}		\
  else{FREE(var->dst);FREE(var->src);}	\
  I->p = S;				\
  I->no = NO;				\
  I->nn = NN;				\
  I->nt = NT;				\
  I->ntab = NTAB;			\
  I->nmc = NMC;				\
  int nmv=0;				\
  for (ii=0; ii < NMC; ii++){		\
    nmv += I->channel[ii].dimension;	\
  }					\
  I->nmv=nmv;				\
  ASSERT(NV==NO+NN+NT+nmv);		\
  return NV;				\
}

#define IONIC_CONST(type,name) type name=S->I.name;
#define IONIC_ARRAY(type,name) type *name=&(S->I.name[0]);

#endif
