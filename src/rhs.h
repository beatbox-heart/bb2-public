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

/* Definitions related to rhs-type cell model definitions */

#ifndef _rhs
#define _rhs

typedef struct {		/* description of dependent parameters */
  int n;			/* # of those */
  int *src;			/* array of layer numbers, relative to the rhs base layer */
  real **dst;			/* array of pointers to dependent parameters */
} Var;

/* rhs run definition wrapper */
#define RHSPROC(name) int name(real *u, real *du, Par par, Var var, int ln)
typedef RHSPROC(RhsProc);

/* Header of a generic rhs calculator:	*/
/* - nostrify the parameters list, 	*/
/* - check number of dynamic variables,	*/
/* - implement parameter substitution.	*/
#define RHS_HEAD(name,LN) \
RHSPROC(name) {           \
  rhspar *S = (rhspar *)par;	  \
  int ivar;               \
  ASSERT(ln==LN);         \
  if (var.n) for(ivar=0;ivar<var.n;ivar++) *(var.dst[ivar])=u[var.src[ivar]]; /* u starts with rhs base layer */

#define RHS_TAIL(name)	\
  return 1;             \
}

/* Solver-independent entities exported by an ionic model description */
typedef struct {
  RhsProc  *f;		/* function of right-hand sides */
  Par p;		/* vector of idiomatic ODE parameters */
  Var var;		/* description of variable parameters */
} rhs_str;


/* rhs create definition wrapper */
#define RHSCREATE(name) int name(Par *par, Var *var, char *w, real **u, int v0)
typedef RHSCREATE(RhsCreate);

#define RHS_CREATE_HEAD(name)				\
RHSCREATE(create_##name) {				\
  rhspar *S = (rhspar *)Calloc(1,sizeof(rhspar));        	\
  char *ptr=w;                                    	\
  int ivar=0;						\
  if (!S) ABORT("cannot create %s",#name);		\
  for(var->n=0;*ptr;var->n+=(*(ptr++)==AT));		\
  if(var->n){CALLOC(var->dst,var->n,sizeof(real *));	\
  CALLOC(var->src,var->n,sizeof(int));}			\
  else{var->src=NULL;(var->dst)=NULL;}

#define RHS_CREATE_TAIL(name,rc)    	\
  var->n=ivar;				\
  if(ivar){REALLOC(var->dst,1L*ivar*sizeof(real *));\
  REALLOC(var->src,1L*ivar*sizeof(int));}	\
  else{FREE(var->dst);FREE(var->src);}	\
  *par = S;				\
  return rc;				\
}

/* Make rhs_str components available for the calling device e.g. euler */
#define RHS_CONST(type,name) type name=S->rhs.name;
#define RHS_ARRAY(type,name) type *name=&(S->rhs.name[0]);

/* I used to be indecisive but now I am not so sure ... */
#define rhspar STR
#endif

