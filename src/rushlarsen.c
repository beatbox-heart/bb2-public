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
 * Rush-Larsen solver for cardiac excitability kinetic models. 
 * Utilizes a special kinetics format, describing HH-type gates separately. 
 * TS 2016: add Matrix Rush Larsen solver for Markov chain (MC).
 */


#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_permutation.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "ionic.h"
#include "qpp.h"

/*++++++++++++++++++++++++++++++++++++++++*/ extern int idev;

typedef struct {
  real ht;			/* time step */
  Name ionic;			/* name of the ionic cell model */
  ionic_str I;			/* ionic cell description */
  int rest;			/* how many steps to do to find resting values */
  int v0, v1;			/* layer range */
  real *u;			/* vector of vars for steady state if any */
  real *du;			/* vector of derivatives */
  Name order;			/* text code of the order of execution */
  int whichorder;		/* numeric code of the order of execution */
  int exp_ngate;		/* if nonzero, non-tabulated gates are steppedn by RL, if zero FE */
  Name exp_mc;			/* text code of the MC method */
  int which_exp_mc;		/* numeric code of the MC method */
  real Vmin, Vmax;		/* limits of voltage in the table (Vmax is approximate) */
  real dV;			/* voltage increment in the table */
  int equilibrate_gates;	/* whether to equilibrate gates in the initial state */
  real *nalp;			/* vector of nontab alphas */
  real *nbet;			/* vector of nontab betas */
  int nV;			/* number of rows in the table */
  real one_o_dV;		/* inverse of voltage increment in the table */
  real *tab;			/* the table of values of functions */
  real *adhoc;			/* array of ad hoc values of tabulable functions */
  real *chains;			/* matrix of MC transition rates */
} STR;

/************************************************/
/* Structure of the state vector:		*/
/* 0 .. no-1: no "other variables"		*/
/* no .. no+nn-1: nn non-tab gates      	*/
/* no+nn .. no+nn+nt-1=nv-1":  nt tab gates	*/
/* So nv=no+nn+nt				*/

/****************************************************************/
/* Structure of the vector of tabulated functions: 		*/
/* 0 .. nt-1 		alphas/a-coefficients	   		*/
/* nt .. 2*nt-1 	betas/b-coefficients	   		*/
/* 2*nt .. ntab-1	everything else		   		*/
/* NB: the first 2*nt values are dual purpose (alp/bet->a/b)	*/

/********************************************************/
/* Conversion of alpha anb beta into the		*/
/* Rush-Larsen linear combination coefficients: 	*/
/* a=alpha/(alpha+beta)*(1-exp(-(alpha+beta)*ht)); 	*/
/* b=exp(-(alpha+beta)*ht). 			        */
/* This is a bit tricky; for debugging might be an idea */
/* to do this in a more straightforward way...          */
static inline void calcab (real *a,real *b,real ht) { 
  /* first a=alp, b=bet */
  (*b)+=(*a);	   /* b=alp+bet         */
  (*a)/=(*b);	   /* a=alp/(alp+bet) */
  (*b)*=ht;	   /* b=(alp+bet)*ht) */
  (*b)=exp(-(*b)); /* b=exp(-(alp+bet)*ht) */
  (*a)*=(1-(*b));  /* a=alp/(alp+bet)*(1-exp(-(alp+bet)*ht)) */
}

int
memcpy_gsl_matrix_complex_to_real (real *dest, gsl_matrix_complex *src, int n)
{
  /* 
     Convert the gsl format of eigenvector matrices to the format used
     in Beatbox.

     As the eigenvectors come from MC, the complex part should be zero
     (asserted in the code) and can be ignored.

     This function copies the entries of array "src" into array
     "dest", both arrays should be of size n*n as we operate on square
     matrices of dimension n.
  */
  /* assert the gsl_matrix "src" is square of dimension n */
  ASSERT(src->size1 == n);
  ASSERT(src->size2 == n);
  /* loop to copy real entries of "src" into "dest" */
  int unsigned i, j;
  for (i=0; i < n; i++)
    {
      for (j=0; j < n; j++)
	{
	  /* assert imaginary part is zero */
	  /* the imaginary part has some tolerance */
	  if ( fabs(src->data[2 * (src->tda * i + j)+1]) != 0)
	    ABORT("Error: Imaginary part of eigenvector matrix is not zero (it is %g).\n", fabs(src->data[2 * (src->tda * i + j)+1]));
	  /* copy real part of "src" to "dest"  matrix */
	  dest[n * i + j] = (real) src->data[2 * (src->tda * i + j)];
	}
    }
  return 1;
}

int
memcpy_gsl_vector_complex_to_real (real *dest, gsl_vector_complex *src, int n)
{
  /* 
     Convert the gsl format of eigenvalue vector to the format used
     in Beatbox.

     As the eigenvalues come from MC, the complex part should be zero
     (asserted in the code) and can be ignored.

     This function copies the entries of array "src" into array
     "dest", both arrays should be of size n (number of eigenvectors).
  */
  /* assert the gsl vector "src" is square of dimension n */
  ASSERT(src->size == n);
  /* loop to copy real entries of "src" into "dest" */
  int unsigned i;
  for (i=0; i < n; i++)
    {
	  /* assert imaginary part is zero */
	  /* the imaginary part has some tolerance */
      if ( fabs(src->data[2 * src->stride * i + 1]) != 0)
	ABORT("Error: Imaginary part of eigenvalue vector is not zero (it is %g).\n", src->data[2 * src->stride * i + 1]);
	  /* copy real part of "src" to "dest"  matrix */
	  dest[i] = (real) src->data[2 * src->stride * i];
    }
  return 1;
}

/* function to obtain matrix rush larsen by eigenvalue decomposition and exponentiation */
/* it returns the mrl in Beatbox format -- real */
int
get_matrix_rush_larsen (real *markov_rates, real ht, real *trans_rates, int n)
{
  unsigned int i, j, k;		/* loop counters */
  double factor;		/* factor for eigenvector normalization */
  /* ******************************************************* */
  /* create and allocate memory in gsl format structure */
  /* right eigenvalues complex  */
  gsl_vector_complex *eval_right_gsl = gsl_vector_complex_alloc ((size_t) n);
  /* left eigenvalues complex */
  gsl_vector_complex *eval_left_gsl = gsl_vector_complex_alloc ((size_t) n);
  /* right eigenvector matrix */
  gsl_matrix_complex *evec_right_gsl = 
    gsl_matrix_complex_alloc ((size_t) n, (size_t) n);
  /* left eigenvector matrix */
  gsl_matrix_complex *evec_left_gsl = 
    gsl_matrix_complex_alloc ((size_t) n, (size_t) n);
  /* workspace for nonsymetric eigenvalue problem */
  gsl_eigen_nonsymmv_workspace *workspace = 
    gsl_eigen_nonsymmv_alloc ((size_t) n); 

  real *tr, *tr_transp;                              /* transposed transition rates matrix of INa */
  real *evec_right, *evec_left, *evals;		/* to eigenvalues and eigenvectors in simple format */
  real *exp_evals_ht, *evec_right_eval;	        /* auxilary arrays for MRL computations */

  /* allocate space to eigenvalues and eigenvectors in Beatbox format */
  CALLOC(tr_transp, n*n, sizeof(real));
  CALLOC(tr, n*n, sizeof(real));
  CALLOC(evec_right, n*n, sizeof(real));
  CALLOC(evec_left, n*n, sizeof(real));
  CALLOC(evals, n, sizeof(real));
  CALLOC(exp_evals_ht, n, sizeof(real));
  CALLOC(evec_right_eval, n*n, sizeof(real));


  /* view to gsl structures */
  gsl_matrix_view rates_matrix = 
    gsl_matrix_view_array (tr, (size_t) n, (size_t) n);
  gsl_matrix_view rates_matrix_transp = 
    gsl_matrix_view_array (tr_transp, (size_t) n, (size_t) n);

  /* copy tr and transposed tr matrix (this is seen by gsl_view rates_matrix_transp )*/
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	{
	  tr[n * i + j] = trans_rates[n * i + j];
	  tr_transp[n * i + j] = trans_rates[n * j + i];
	}
    }

  /* compute eigenvalues and eigenvectors */
  /* get the RIGHT evals and evecs */
  gsl_eigen_nonsymmv (&rates_matrix.matrix, eval_right_gsl, evec_right_gsl,
		      workspace);
  /* get the LEFT evals and evecs */
  gsl_eigen_nonsymmv (&rates_matrix_transp.matrix, eval_left_gsl, evec_left_gsl,
		      workspace);	
  /* ****************************** */
  /* process eval and evec */
  /* sort the eigenvalues and eigenvectors in descending order */
  gsl_eigen_nonsymmv_sort (eval_right_gsl, evec_right_gsl, GSL_EIGEN_SORT_ABS_DESC);
  gsl_eigen_nonsymmv_sort (eval_left_gsl, evec_left_gsl, GSL_EIGEN_SORT_ABS_DESC);
  /* convert gsl arrays to real */
  memcpy_gsl_matrix_complex_to_real (evec_left, evec_left_gsl, n);
  memcpy_gsl_matrix_complex_to_real (evec_right, evec_right_gsl, n);
  memcpy_gsl_vector_complex_to_real (evals, eval_left_gsl, n);

  /* reorder eigenvectors in case of multiple (so far only double) eigenvalues */
  #define M 2
  real  eta[M*M];
  real tmp;
  for (i = 0; i < (n - 1); i++)
    {
      /* the "identical" eigenvalues are computed with some tolerance if (evals[i] = evals[i + 1]) */
      if (fabs(evals[i] - evals[i + 1]) < 1e-13 )
	{
	  /* check for triple eigenvalues -- this has not been resolved */
	  if ((i < n - 2) ? (evals[i] == evals[i + 2]) : (0))
	    {
	      URGENT_MESSAGE ("There are three identical eigenvalues.\n");
	      ABORT("");
	    }
	  /* reset the eta's */
	  for (j = 0; j < n; j++)
	    eta[j] = 0.0;
	  for (j = 0; j < n; j++)
	    {
	      eta[0] += evec_right[n * j + i] * evec_left[n * j + i];
	      eta[1] += evec_right[n * j + i + 1] * evec_left[n * j + i];
	      eta[2] += evec_right[n * j + i] * evec_left[n * j + i + 1];
	      eta[3] += evec_right[n * j + i + 1] * evec_left[n * j + i + 1];
	    }
	  /* if the eigenvector are in reversed order swap right eigenvectors */
	  if (eta[0] == 0 && eta[3] == 0)
	    {
	      for (j = 0; j < n; j++)
		{
		  tmp = evec_right[n * j + i];
		  evec_right[n * j + i] = evec_right[n * j + i + 1];
		  evec_right[n * j + i + 1] = tmp;
		}
	    }
	  for (j=0;j<M*M; j++) eta[j] = 0.0;
	  for (j = 0; j < n; j++)
	    {
	      eta[0] += evec_right[n * j + i] * evec_left[n * j + i];
	      eta[1] += evec_right[n * j + i + 1] * evec_left[n * j + i];
	      eta[2] += evec_right[n * j + i] * evec_left[n * j + i + 1];
	      eta[3] += evec_right[n * j + i + 1] * evec_left[n * j + i + 1];
	    }
	  if (eta[0] == 0 || eta[3] == 0 || eta[1] != 0 || eta[2] != 0 )
	    {
	      URGENT_MESSAGE ("The eigenvectors are not biorthogonal.\n");
	      ABORT ("");
	    }
	}
    }

  /* scale left eigenvectors such that the multiplication
     with right eigenvectors gives identity */
  for (i = 0; i < n; i++)
    {
      factor = 0.0;
      for (j = 0; j < n; j++)
	{
	  factor += evec_right[n * j + i] * evec_left[n * j + i];
	}
      factor = 1/factor;
      for (j = 0; j < n; j++)
	{
	  evec_right[n * j + i] *= factor;
	}
    }
  /* exponentiate eigenvalue and time step */
  for (i = 0; i < n; i++)
    exp_evals_ht[i] = exp (evals[i] * ht);
  /* multiply left eigenvalue and eigenvectors and place to auxilary array */
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)		  
	{
	  evec_right_eval[n * i + j] = evec_right[n * i + j] * exp_evals_ht[j];
	}
    }
  /* multiply auxilary array with right eigenvectors so we get the MRL */
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)		  
	{
	  markov_rates[n * i + j] = 0.0;
	  for (k = 0; k < n; k++)
	    {
	      markov_rates[n * i + j] += 
		evec_right_eval[n * i + k] * evec_left[n * j + k];
	    }
	}
    }

  /* free gsl memory */
  FREE(tr_transp);
  FREE(tr);
  FREE(evec_right);
  FREE(evec_left);
  FREE(evals);
  FREE(exp_evals_ht);
  FREE(evec_right_eval);

  gsl_vector_complex_free (eval_right_gsl);
  gsl_vector_complex_free (eval_left_gsl);
  gsl_matrix_complex_free (evec_right_gsl);
  gsl_matrix_complex_free (evec_left_gsl);
  gsl_eigen_nonsymmv_free (workspace);

  return 1;
}

/* print matrix -- for debugging purposes*/
void
print_matrix (real *m, int n)
{
  /* algorithm to print square matrix m of dimension n */
  int i, j;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	{
	  printf ("%10.2g", m[n * i + j]);
	}
      printf ("\n");
    }
  printf ("\n");
}

/* print matrix -- for debugging purposes*/
void
print_complex_matrix (real *m, int n)
{
  /* algorithm to print square complex matrix m of dimension n */
  int i, j;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < 2*n; j=j+2)
	{
	    printf ("(%0.2g", m[2*n * i + j]);
	    printf ("%+0.2gi)\t", m[2*n * i + j + 1 ]);
	}
      printf ("\n");
    }
  printf ("\n");
}

/* Substeps of the algorithm can be done in different order */

#define XYZ xof(u),yof(u),zof(u)

/***********************/
/* TABULATED FUNCTIONS */
/* check voltage */							
#define DOTABLES \
  V=u[V_index];								\
  iV=floor((V-Vmin)*one_o_dV);						\
  /* if outside limits or tabulation step dV is zero calculate */	\
  /* tabulated values by formulas */					\
  if (iV<0 || iV>=nV || S->dV == 0.0 ) {				\
    if (S->dV != 0.0) printf("V=%g outside [%g,%g] at t=%ld at point %ld,%ld,%ld\n",V,Vmin,Vmax,t,XYZ); \
    values=adhoc;							\
    if (!ftab(V,values,ntab))						\
      ABORT("\nerror calculating ftab(%s) at t=%ld point %ld,%ld,%ld: V=%g\n",ionic,t,XYZ,V); \
    for (it=0;it<nt;it++)						\
      calcab(values+it,values+nt+it,ht);				\
    /* TODO: report the value of V, to extend the table in the future */\
  } else {								\
    /* otherwise get them from the table */				\
    values=tab+iV*ntab;							\
  } /* if iV .. else */

/************************/
/* tabulated GATES by Rush-Larsen */
#define DOTGATES \
  a=values;								\
  b=values+nt;								\
  for (it=0;it<nt;it++) {						\
    gateold=tgate[it];							\
    tgate[it]=a[it]+b[it]*gateold; /* this is the RL step */		\
    if (!isfinite(tgate[it])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) at t=%ld x=%ld,y=%ld,z=%ld tab gate %d\n",t,XYZ,it); \
      URGENT_MESSAGE("This happened while multiplying %g by %g and adding %g\n",gateold,b[it],a[it]); \
      ABORT("");							\
    } /* if !isfinite */						\
  } /* for it */

/************************/
/* equilibrating tabulated GATES */
#define DOTGATESEQ \
  a=values;								\
  b=values+nt;								\
  for (it=0;it<nt;it++) {						\
    gateold=tgate[it];							\
    tgate[it]=a[it]/(1-b[it]); /* this is the equilibrium value */	\
    if (!isfinite(tgate[it])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) equilibrating tab gate %d: a=%g, b=%g\n",it,a[it],b[it]); \
      ABORT("");							\
    } /* if !isfinite */						\
  } /* for it */


/*************************************************************/
/* OTHER VARIABLES by forward Euler */
#define DOOTHER \
  if (!fddt(u,nv,values,ntab,p,var,du,no,nalp,nbet,nn)) {		\
    URGENT_MESSAGE("\nerror calculating fddt(%s) at t=%ld point %ld,%ld,%ld: u=",ionic,t,XYZ); \
    for(iv=0;iv<nv;iv++) URGENT_MESSAGE(" %lg",u[iv]);			\
    ABORT("\n");							\
  } /*  if !fddt... */							\
  for (io=0;io<no;io++) {						\
    if (!isfinite(du[io])) {						\
      URGENT_MESSAGE("\nerror calculating fddt(%s) at t=%ld point %ld,%ld,%ld: u=",ionic,t,XYZ); \
      for(iv=0;iv<nv;iv++) URGENT_MESSAGE("%c%lg",iv?',':'(',u[iv]);	\
      URGENT_MESSAGE(") ->");						\
      for(jo=0;jo<no;jo++) URGENT_MESSAGE("%c%lg",jo?',':'(',du[jo]);	\
      URGENT_MESSAGE(")\n");						\
      ABORT("");							\
    } /* if du=NaN */							\
    u[io] += ht*du[io]; /* this is the FE step */			\
    if (!isfinite(u[io])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) detected at t=%ld point %ld,%ld,%ld v=%d\n",t,XYZ,io); \
      URGENT_MESSAGE("This happened while incrementing ");		\
      for(iv=0;iv<nv;iv++) URGENT_MESSAGE("%c%lg",iv?',':'(',u[iv]);	\
      URGENT_MESSAGE(") by %lg*",ht);					\
      for(jo=0;jo<no;jo++) URGENT_MESSAGE("%c%lg",jo?',':'(',du[jo]);	\
      URGENT_MESSAGE(")\n");						\
      ABORT("");							\
    } /* if u=NaN */							\
  } /* for io */

/*************************************************************/
/* NONTAB GATES by Rush-Larsen                               */
#define DONGATESRL \
  for (in=0;in<nn;in++) {						\
    a=nalp+in;								\
    b=nbet+in;								\
    calcab(a,b,ht);							\
    gateold=ngate[in];							\
    ngate[in]=(*a)+(*b)*gateold; /* this is the RL step */		\
    if (!isfinite(ngate[in])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) at t=%ld point %ld,%ld,%ld nontab gate %d\n",t,XYZ,in); \
      URGENT_MESSAGE("This happened while multiplying %g by %g and adding %g\n",gateold,*b,*a); \
      ABORT("");							\
    } /* if !isfinite */						\
  } /* for in */


/*************************************************************/
/* equilibrating NONTAB GATES                                */
#define DONGATESEQ \
  for (in=0;in<nn;in++) {						\
    a=nalp+in;								\
    b=nbet+in;								\
    ngate[in]=(*a)/((*a)+(*b)); /* this is the equilibrium value */	\
    if (!isfinite(ngate[in])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) equilibrating nontab gate %d: alp=%g, bet=%g\n",in,*a,*b); \
      ABORT("");							\
    } /* if !isfinite */						\
  } /* for in */
  
/*************************************************************/
/* NONTAB GATES by forward Euler                             */
#define DONGATESFE \
  for (in=0;in<nn;in++) {						\
    a=nalp+in;								\
    b=nbet+in;								\
    gateold=ngate[in];							\
    ngate[in] = gateold+ht*((*a)-((*a)+(*b))*gateold); /* this is the FE step */ \
    if (!isfinite(ngate[in])) {						\
      URGENT_MESSAGE("\nNAN (not-a-number) at t=%ld point %ld,%ld,%ld nontab gate %d\n",t,XYZ,in); \
      URGENT_MESSAGE("This happened at gateold=%g ht=%g alpha=%g beta=%g\n",gateold,ht,*b,*a); \
      ABORT("");							\
    } /* if !isfinite */						\
  } /* for in */


#define DOMARKOV \
  for (im=0;im<nmc;im++) {						\
    subchain = &(channel[im].subchain[0]);				\
    dimension = channel[im].dimension;					\
    									\
    CALLOC(trm, dimension*dimension, sizeof(real));			\
    CALLOC(mrl, dimension*dimension, sizeof(real));			\
    CALLOC(hchain, dimension, sizeof(real));				\
    									\
    for (jm=0;jm<channel[im].num_sub;jm++) {				\
	/* limits of the tabulation */					\
	dvar = subchain[jm].tincr;					\
	tmax = subchain[jm].tmax;					\
	tmin = subchain[jm].tmin;					\
									\
	if (subchain[jm].i_control >= 0) {				\
	  /* find out the scale of the trans_rates_mat function */	\
	  if (subchain[jm].scale == 0) {				\
	    /* linear scale */						\
	    valscale[0] = u[subchain[jm].i_control];			\
	  } else if (subchain[jm].scale == 1) {				\
	    /* logarithmic scale */					\
	    valscale[0] = log(u[subchain[jm].i_control]);		\
	  } else {							\
	    EXPECTED_ERROR("unknown tabulation scale %d", subchain[jm].scale); \
	  } /* if .. scale .. else */					\
	  ch_exp_mc = which_exp_mc;					\
									\
	  val=&(valscale[0]);						\
	  it=floor((val[0]-tmin)/dvar);					\
	  if (dvar > 0) {						\
	      nT=ceil( (subchain[jm].tmax - subchain[jm].tmin) / dvar );\
	  } else {							\
	      nT = 0;							\
	    } /* if dvar>0 else */					\
	} else {							\
	  /* for channels without unique controling variable */		\
	  /* the values will be computed on fly */			\
	  val = u;							\
	  nT=0;								\
	  it = -1;	/* places the table index out of range */	\
	  ch_exp_mc = 0;	/* computes using forward Euler */	\
	} /* if .. control .. else */					\
									\
	/* if the index is outside of precomputed range  */		\
	/* find the solution on fly */					\
	if ( (it<0 || it>=nT) || ch_exp_mc == ntabmrl) {		\
	  on_fly = 1;							\
	  /* compute transition rates matrix */				\
	  if (!(channel[im].subchain[jm].trans_rates_mat)(val, trm))	\
	    ABORT("\nerror calculating mtab(%s) table for V=%g\n",S->ionic,V); \
	  /* asign matrix for forward Euler calculations */		\
	  markov_adhoc = trm;						\
	  /* exponential integrator -- matrix_rush_larsen */		\
	  if (ch_exp_mc) {						\
	    if (!get_matrix_rush_larsen (mrl, ht, trm, dimension))	\
	      URGENT_MESSAGE("\nerror calculating mrl(%s) table for V=%g\n",S->ionic,V); /* mention the index of mc */ \
	    markov_adhoc = mrl;						\
	  } /* if ch_exp_mc */						\
	} else {							\
	  on_fly = 0;							\
	  /* pointer to precomputed and tabulated matrix */		\
	  markov_adhoc = &(curr_chain[it*dimension*dimension]);		\
	} /* if it<0 .. else */						\
									\
	/* integration method */					\
	mat_vec_mult(hchain, markov_adhoc, markov, dimension);		\
	if (ch_exp_mc) {						\
	  /* MRL step */						\
	  for (ii=0; ii<dimension; ii++)				\
	    markov[ii] = hchain[ii];					\
	} else if (ch_exp_mc == mcfe) {					\
	  /* FE step */							\
	  for (ii=0; ii<dimension; ii++)				\
	    markov[ii] += ht*hchain[ii];				\
	} else {							\
	  EXPECTED_ERROR("The which_exp_mc == %d is not supported",which_exp_mc); \
	} /* if ch_exp_mc .. else */					\
									\
	/* check if the results are finite and within range */		\
	sum = 0.0;							\
	for (ii=0; ii<dimension; ii++) {				\
	  sum += markov[ii];						\
	  if (!isfinite(markov[ii])) {					\
	    URGENT_MESSAGE("\nNAN (not-a-number) detected at t=%ld point %ld,%ld,%ld v=%d\n",t,XYZ,io); \
	    URGENT_MESSAGE("This happened while incrementing ");	\
	    for(iv=0;iv<nv;iv++) URGENT_MESSAGE("%c%lg",iv?',':'(',u[iv]); \
	    URGENT_MESSAGE(") by %lg*",ht);				\
	    for(jo=0;jo<no;jo++) URGENT_MESSAGE("%c%lg",jo?',':'(',du[jo]); \
	    URGENT_MESSAGE(")\n");					\
	    ABORT("");							\
	  } /* if NaN */						\
	} /* for ii */							\
									\
	/* increment channel pointers */				\
	curr_chain += dimension*dimension*nT;				\
	/* reset auxilary arrays */					\
	if (on_fly) {							\
	  for (ii=0; ii<dimension*dimension; ii++) {			\
	    trm[ii]=0.0;						\
	    mrl[ii]=0.0;						\
	  } /* for ii */						\
	} /* on_fly */							\
	/* reset vector with solution increment */			\
	for (ii=0;ii< dimension;ii++) hchain[ii] = 0.0;			\
    } /* for jm */							\
    /* increment markov variable pointer */				\
    markov += dimension;						\
    /* free auxilary arrays */						\
    FREE(trm);								\
    FREE(mrl);								\
    FREE(hchain);							\
  } /* for im */


/* TODO: check if the new u is finite at the end of the time step,
   rather than in each subunit calculation */

enum {
  tgo,
  tog,
  totg,
  numorders
} ordertype;

enum {
  mcfe,
  tabmrl,
  ntabmrl
} mrltype;

/* do matrix vector multiplication of dimension N as dest = mat * vect  */
static inline void mat_vec_mult (real *dest,real *mat,real *vec,int N)
{
  int im, jm;  
  for (im=0; im<N; im++) {
    dest[im] = 0.0;		/* reset destination vector */
    for (jm=0; jm<N; jm++) {
      dest[im] += mat[im*N+jm] * vec[jm];
    }
  }
 } /* mat_vec_mult */

/* compute the rushlarsen step for one cell of the mesh */
/* nv is not used */
static inline int rushlarsen_step (real *u,int nv,STR *S)
{
  IONIC_CONST(Par,p);			/* set of ionic cell parameters */
  IONIC_CONST(Var,var);			/* description of variable parameters */
  IONIC_CONST(IonicFddt *,fddt);	/* the rhs functions except gates */
  IONIC_CONST(IonicFtab *,ftab);	/* the tabulated functions calculator */
  IONIC_CONST(int,no);			/* number of non-gate variables in the state vector */
  IONIC_CONST(int,nn);			/* num of nontabulated gate variables in the state vector */
  IONIC_CONST(int,nt);			/* num of tabulated gate variables in the state vector */
  IONIC_CONST(int,ntab);		/* number of tabulated functions */
  IONIC_CONST(int,nmc);			/* number of Markov chain models */
  IONIC_CONST(int,nmv);			/* number of Markov chain models variables */
  IONIC_ARRAY(channel_str,channel);	/* markov chain structures */
  IONIC_CONST(int,V_index);		/* index of V in state vector */
  DEVICE_ARRAY(char,ionic);		/* name of the ionic cell model */
  DEVICE_CONST(real,ht);		/* the time step */
  DEVICE_CONST(int,whichorder);		/* numeric code of the order of execution */
  DEVICE_CONST(int,which_exp_mc);	/* if nonzero, markov chains are stepped by MRL, otherwise FE */
  DEVICE_CONST(int,nV);			/* number of 'rows' in the table */
  DEVICE_CONST(real,Vmin);		/* minimal value of V in the table */
  DEVICE_CONST(real,Vmax);		/* maximal value of V in the table */
  DEVICE_CONST(real,one_o_dV);		/* the inverse of increment of V in the table */
  DEVICE_ARRAY(real,du/*[no]*/);	/* array of right-hand sides for the FE part */
  DEVICE_ARRAY(real,tab/*[ntab*nV]*/);	/* the table of function values */
  DEVICE_ARRAY(real,nalp/*[nn]*/);	/* array of alphas for the nontab RL part */
  DEVICE_ARRAY(real,nbet/*[nn]*/);	/* array of betas for the nontab RL part */
  DEVICE_ARRAY(real,adhoc/*[ntab]*/);	/* freshly calculated values of functions and matrix rush larsen */
  int io, jo, in, it, it1, iv;		/* vector components counters */
  int iV;				/* table row counter */
  real V;				/* transmembrane voltage value */
  real *ngate;				/* subvector of nontab gate values */
  real *tgate;				/* subvector of tab gate values */
  real *values;				/* vector of values of tabulated functions */
  real *a, *b;				/* pointers to subvectors of Rush-Larsen step coefficients */
  real gateold;				/* aux variable */
  real *markov,*markov_cur;		/* subvector of markov chain values */
  real *dmarkov,*dmarkov_cur;		/* subvector of markov chain values increments */
  real markov_entry;			/* auxilary markov state */
  real *mrl;				/* matrix for rush-larsen computations */
  int im,jm,km;				/* vector components counters */
  int ii;				/* MC counter */
  real *curr_chain;			/* pointer to current chain matrix */
  subchain_str *subchain;		/* pointer to current subchain matrix */
  int dimension;			/* subchain matrix size */
  int nT;				/* num elements in the tabulated array */
  int on_fly;				/* whether the markov_adhoc was done on fly or obtained from tabulated data */
  real *dchain, *hchain, *trm, *markov_adhoc, *val ; 
  real valscale[1];
  real dvar, tmax, tmin;
  real sum;				/* variable to check the sum of states */
  int ch_exp_mc;    


  CALLOC(dmarkov,nmv,sizeof(real));

  ngate=u+no;
  tgate=u+no+nn;
  markov=u+no+nn+nt;
  curr_chain=&(S->chains[0]);
  switch (whichorder) {
  case tgo:
    DOTABLES;
    DOTGATES;
    DOOTHER;
    DOMARKOV;
    DONGATESFE;
    break;
  case tog:
    DOTABLES;
    DOOTHER;
    DOMARKOV;
    DONGATESFE;
    DOTGATES;
    break;
  case totg:
    DOTABLES;
    DOOTHER;
    DOMARKOV;
    DONGATESFE;
    DOTABLES;
    DOTGATES;
    break;
  case tgo+numorders:
    DOTABLES;
    DOTGATES;
    DOOTHER;
    DOMARKOV;
    DONGATESRL;
    break;
  case tog+numorders:
    DOTABLES;
    DOOTHER;
    DOMARKOV;
    DONGATESRL;
    DOTGATES;
    break;
  case totg+numorders:
    DOTABLES;
    DOOTHER;
    DOTABLES;
    DOMARKOV;
    DONGATESRL;
    DOTABLES;
    DOTGATES;
    break;
  default:
    EXPECTED_ERROR("unknown order of execution code %d\n",whichorder);
  } /* swicth whichorder */
  FREE(dmarkov);
  return 1;
}

static inline int equilibration_step(real *u,int nv,STR *S)
{
  IONIC_CONST(Par,p);			/* set of ionic cell parameters */
  IONIC_CONST(Var,var);			/* description of variable parameters */
  IONIC_CONST(IonicFddt *,fddt);	/* the rhs functions except gates */
  IONIC_CONST(IonicFtab *,ftab);	/* the tabulated functions calculator */
  IONIC_CONST(int,no);			/* number of non-gate variables in the state vector */
  IONIC_CONST(int,nn);			/* num of nontabulated gate variables in the state vector */
  IONIC_CONST(int,nt);			/* num of tabulated gate variables in the state vector */
  IONIC_CONST(int,ntab);		/* number of tabulated functions */
  IONIC_CONST(int,V_index);		/* index of V in state vector */
  DEVICE_ARRAY(char,ionic);		/* name of the ionic cell model */
  DEVICE_CONST(real,ht);		/* the time step */
  DEVICE_CONST(int,whichorder);		/* numeric code of the order of execution */
  DEVICE_CONST(int,nV);			/* number of 'rows' in the table */
  DEVICE_CONST(real,Vmin);		/* minimal value of V in the table */
  DEVICE_CONST(real,Vmax);		/* maximal value of V in the table */
  DEVICE_CONST(real,one_o_dV);		/* the inverse of increment of V in the table */
  DEVICE_ARRAY(real,du/*[no]*/);	/* array of right-hand sides for the FE part */
  DEVICE_ARRAY(real,tab/*[ntab*nV]*/);	/* the table of function values */
  DEVICE_ARRAY(real,nalp/*[nn]*/);	/* array of alphas for the nontab RL part */
  DEVICE_ARRAY(real,nbet/*[nn]*/);	/* array of betas for the nontab RL part */
  DEVICE_ARRAY(real,adhoc/*[ntab]*/);	/* freshly calculated values of functions */
  int io, jo, in, it, it1, iv;		/* vector components counters */
  int iV;				/* table row counter */
  real V;				/* transmembrane voltage value */
  real *ngate;				/* subvector of nontab gate values */
  real *tgate;				/* subvector of tab gate values */
  real *values;				/* vector of values of tabulated functions */
  real *a, *b;				/* pointers to subvectors of Rush-Larsen step coefficients */
  real gateold;				/* aux variable */
  int x=-1;				/* these are */
  int y=-1;				/*   required */
  int z=-1;				/*   by DOTABLES */

  ngate=u+no;
  tgate=u+no+nn;

  DOTABLES;
  DOOTHER;
  DONGATESEQ;
  DOTGATESEQ;

  return 1;
}

/******************/
RUN_HEAD(rushlarsen)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  int nv=v1-v0+1;			/* size of the state vector */
  int x, y, z;				/* space grid counters */
  real V;				/* transmembrane voltage value */
/* #define COMMANDS { if NOT(rushlarsen_step(u,nv,S,(*x),(*y),(*z))) return 0; } */
#define COMMANDS { if NOT(rushlarsen_step(u,nv,S)) return 0; }
  DO_FOR_ALL_POINTS;
#undef COMMANDS
}
RUN_TAIL(rushlarsen)

/**********************/
DESTROY_HEAD(rushlarsen)
{
  FREE(S->du);
  FREE(S->nalp);
  FREE(S->nbet);
  FREE(S->adhoc);
  FREE(S->tab);
  FREE(S->u);
  if (S->I.var.n) {
    FREE(S->I.var.src);
    FREE(S->I.var.dst);
  }
  FREE(S->I.p);
}
DESTROY_TAIL(rushlarsen)

/* Declare all available ionic models */
#define D(a) IonicFtab ftab_##a;
#include "ioniclist.h"
#undef D
#define D(a) IonicFddt fddt_##a;
#include "ioniclist.h"
#undef D
#define D(a) IonicCreate create_##a;
#include "ioniclist.h"
#undef D

/*********************/
CREATE_HEAD(rushlarsen)
{
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 >= v0 );
  int nv=v1-v0+1;
  int no, nn, nt, ntab, nV, nmc, nmv;
  channel_str * channel;
  int size_tr=0;
  /* int *nm; */
  int step, iv, ix, iy, iz;
  int in, it, io, jo, iV;
  real *ufull, *u, *values;
  real *tr, *mrl, *adhoc;
  real V, alp, bet;
  int im,jm,ii;

  /* Create tables for Markov chain models */
  /* goes to the top */
  real dvar, tmax, tmin, one_o_dvar;
  real * trm, * markov_adhoc;			/* pointer to the tr_tab */
  subchain_str  * subchain;
  int dimension, nT;
  real var[1];

  /* Accept the time step */
  ACCEPTR(ht,RNONE,0.,RNONE);
  if (ht==0) MESSAGE("/* WARNING: ht=0 is formally allowed but hardly makes sense */");

  /* Accept the execution order */
  ACCEPTS(order,"totg");
  STRSWITCH(order);
  STRCASE("tgo")  S->whichorder=tgo;
  STRCASE("tog")  S->whichorder=tog;
  STRCASE("totg") S->whichorder=totg;
  STRDEFAULT EXPECTED_ERROR("\nrushlarsen: unknown execution order '%s'\n",order); 
  STRENDSW

  ACCEPTI(exp_ngate,0,0,1);
  if (exp_ngate) S->whichorder+=numorders;

  /* Accept the MC integration method */
  /* TODO: since now it will contains only mcfe and tabmrl, perahps should be numeric */
  ACCEPTS(exp_mc,"tabmrl");
  STRSWITCH(exp_mc);
  STRCASE("mcfe")  S->which_exp_mc=mcfe;
  STRCASE("tabmrl")  S->which_exp_mc=tabmrl;
  STRCASE("ntabmrl") S->which_exp_mc=ntabmrl;
  STRDEFAULT EXPECTED_ERROR("\nrushlarsen: unknown MC integration method exp_mc '%s'\n",exp_mc); 
  STRENDSW

  /* Accept the ionic cell model */
  ACCEPTS(ionic,NULL);
  {
    char *pars;
    MALLOC(pars,(long)MAXSTRLEN);
    BEGINBLOCK("par=",pars);
    S->u=NULL;
    #define D(a)							\
    if (0==stricmp(S->ionic,#a)) {					\
      if NOT(create_##a(&(S->I),pars,&(S->u),v0))		\
    	EXPECTED_ERROR("reading parameters for %s in \"%s\"",S->ionic,pars); \
      no=S->I.no; nn=S->I.nn; nt=S->I.nt;				\
      nmc=S->I.nmc;							\
      channel=&(S->I.channel[0]);					\
      nmv=S->I.nmv;							\
      ntab=S->I.ntab;							\
      if (no+nn+nt+nmv!=nv)						\
        EXPECTED_ERROR("no+nn+nt+nmv=%d+%d+%d+%d != nv=%d for %s\n"	\
		       "Please correct the number of layer in the input file (BBScript).\n",no,nn,nt,nmv,nv,S->ionic); \
      S->I.ftab=ftab_##a;						\
      S->I.fddt=fddt_##a;						\
    } else
    #include "ioniclist.h"
    #undef D
    EXPECTED_ERROR("unknown ionic model %s",S->ionic);
    ENDBLOCK;
    FREE(pars);
  }

  /* compute required sizes for MC */
  for (im=0; im<nmc; im++) {
    /* the table is only created for subchains which are tabulated
       i.e. subchains with specified i_control variable (only single
       variable dependent markov chains which are suitable for
       tabulation). */
    /* the corresponding size for each subchain is the
       ntab_entries*dimension^2 */
    
    subchain = &(S->I.channel[im].subchain[0]);
    
    for (jm=0; jm<S->I.channel[im].num_sub; jm++) {
      dvar = subchain[jm].tincr;
      dimension = channel[im].dimension;
      tmax = subchain[jm].tmax;
      tmin = subchain[jm].tmin;
      if (subchain[jm].i_control >= 0 && dvar>0) {
	nT = ceil ((tmax - tmin) / dvar);
	size_tr += dimension * dimension * nT;
      } /* if .. control */
    } /* for jm */
  } /* for im */

  /* Allocated working arrays */
  CALLOC(S->du, no, sizeof (real));
  CALLOC(S->nalp, nn, sizeof (real));
  CALLOC(S->nbet, nn, sizeof (real));
  CALLOC(S->adhoc, ntab, sizeof (real));	/* creates space for mrl although might not be needed */
  CALLOC(S->chains, size_tr, sizeof (real));

  /* tabulate Markov chain transition rates matrices */
  markov_adhoc = &(S->chains[0]);
  for (im=0; im<nmc; im++) {
    subchain = &(S->I.channel[im].subchain[0]);
    dimension = S->I.channel[im].dimension;
    
    CALLOC (trm, dimension * dimension, sizeof (real));

    for (jm=0; jm<channel[im].num_sub; jm++) {
      dvar = subchain[jm].tincr;
      tmax = subchain[jm].tmax;
      tmin = subchain[jm].tmin;
      if (subchain[jm].i_control >= 0 && dvar>0) {
	nT = ceil ((subchain[jm].tmax - subchain[jm].tmin) / dvar);
	one_o_dvar = 1.0 / dvar;
	
	for (it=0; it<nT; it++) {
	  /* tabulate transition rates of gate and Markov chain and other variables */
	  /* TODO: if the scale is logarithmic, should it still be (it + 0.5)? */
	  var[0] = tmin + (it + 0.5) * dvar;
	  if (!(S->I.channel[im].subchain[jm].trans_rates_mat) (var, trm))
	    ABORT ("\nerror calculating mtab(%s) table for V=%g\n", S->ionic, V);
	  /* get exponential integrator -- matrix_rush_larsen */
	  if (S->which_exp_mc == tabmrl) {
	    if (!get_matrix_rush_larsen (&(markov_adhoc[it * dimension * dimension]), ht, trm, dimension))
			URGENT_MESSAGE ("\nerror calculating mrl(%s) table for V=%g\n", S->ionic, V);
	  } else if (S->which_exp_mc == mcfe) {
	    memcpy (&(markov_adhoc[it * dimension * dimension]), trm, dimension * dimension * sizeof (real));
	  }
	  
	  for (ii=0; ii<dimension*dimension; ii++) {
	    trm[ii] = 0.0;
	    if (markov_adhoc[it*dimension*dimension + ii] != markov_adhoc[it*dimension*dimension + ii]) {
	      URGENT_MESSAGE ("\nNaN in calculation of markov_adhoc(%s) table for var=%.15g\n", S->ionic, var[0]);
	      ABORT ("");
	    } /* if markov_adhoc .. */
	  } /* for ii */
	} /* for it */
	markov_adhoc += dimension * dimension * nT;
      } /* if .. control .. */
    } /* for jm */
    FREE (trm);
  } /* for im */


  /* If the ionic model allocated initial state, possibly */
  /* finalize it by calculating stationary values of the gates */
#if 0
  if (S->u) {
    ACCEPTI(equilibrate_gates,0,0,1);
    V=S->u[S->I.V_index];
    values=S->adhoc;
    
    if (!(S->I.fddt(S->u,nv,values,ntab,S->I.p,S->I.var,S->du,no,S->nalp,S->nbet,nn))) {
      URGENT_MESSAGE("\nerror calculating fddt(%s) at parse time: u=",ionic);
      for(iv=0;iv<nv;iv++) URGENT_MESSAGE(" %lg",u[iv]);
      ABORT("\n");
    } /*  if !fddt... */
    for (in=0;in<nn;in++) {
      alp=S->nalp[in];
      bet=S->nbet[in];
      S->u[no+in]=alp/(alp+bet);
    } /* for in */
    if (!(S->I.ftab)(V,values,ntab))
      EXPECTED_ERROR("\nerror calculating ftab(%s) for V=%g\n",S->ionic,V);
    for (it=0;it<nt;it++) {
      calcab(values+it,values+nt+it,ht);
      alp=values[it];
      bet=values[nt+it];
      S->u[no+nn+it]=alp/(alp+bet);
    } /* for it */
  } else {
    if (find_key("equilibrate_gates=",w)) {
      MESSAGE("\n/* Warning: equilibrate_gate parameter specified for a model which "
  	      " does not define a standard state. "
  	      "The parameter will be ignored. */\n");
    } /* if findkey */
  } /* if S->u else */
#endif	/* 0 -- comment */
  
  /* In any case, do the tabulation */
  ACCEPTR(Vmin,-200,RNONE,RNONE);
  ACCEPTR(Vmax,+200,RNONE,RNONE);
  ACCEPTR(dV,0.01,0.,RNONE);
  if (dV == 0.0) {
    MESSAGE("/* NOTE: dV=0 means the transition rates are calculated on fly and not tabulated. */");
    CALLOC(S->tab,1,sizeof(real)); /* for compatibility reasons,
				      otherwise there are problems
				      with destruction of the
				      device */
  } else {
    S->nV=nV=ceil((S->Vmax-S->Vmin)/S->dV);
    S->one_o_dV=1.0/S->dV;
    CALLOC(S->tab,ntab*nV,sizeof(real));
    
    for (iV=0;iV<nV;iV++) {
      /* tabulate transition rates of gate and Markov chain and other variables */
      V=Vmin+(iV+0.5)*dV;
      var[0]=V;
      values=&(S->tab[ntab * iV]);
      if (!(S->I.ftab)(V,values,ntab))
	ABORT("\nerror calculating ftab(%s) table for V=%g\n",S->ionic,V);
      for (it=0;it<nt;it++)
	calcab(values+it,values+nt+it,ht);
    } /* for iV */
  } /* if dV==0 else */

  ACCEPTI(rest,0,0,INONE);
  if (S->rest) {
    /* Need full vector, in case variable parameters use other layers.   */
    /* NB other layers may be above as well as below this device's space */
    CALLOC(ufull,vmax,sizeof(real)); 
    u=&(ufull[v0]);

    if (S->u) {
      MESSAGE("\n/* NOTICE: finding resting state while already defined by the model */");
      /* use that as the initial condition */
      for (iv=0;iv<nv;iv++) u[iv]=S->u[iv];
    } else {
      CALLOC(S->u,nv,sizeof(real));
    }

    if (S->I.var.n) MESSAGE("\n/* NOTICE: finding resting state while parameters are variable */");

    for (step=0;step<S->rest;step++) 
      if NOT(equilibration_step(u,nv,S)) return 0;

    /* copy what is needed */
    for (iv=0;iv<nv;iv++) S->u[iv]=u[iv];

    /* don't need this one any more */
    FREE(ufull);
  } /* if S->rest */

  if (S->u) {
    MESSAGE0("\n/* Initial state: "); for(iv=0;iv<nv;iv++) MESSAGE1("%lg ",(S->u)[iv]); MESSAGE0("*/");
    /* Fill up the whole of the grid with the resting state */
#if MPI
    if (dev->s.runHere) {
#endif
      Space s=dev->s;
#define COMMANDS for (iv=0;iv<=nv;iv++) {u[v0+iv]=(S->u)[iv];}
      DO_FOR_ALL_POINTS;
#undef COMMANDS
      
    #if MPI
    } /* if runHere */
    #endif
  } else {
    /* No resting state was defined */
    MESSAGE("\n/* Notice: no standard state defined for ionic model '%s' */\n",ionic);
    CALLOC(S->u,nv,sizeof(real));
  } /* if S0-> u else */
}
CREATE_TAIL(rushlarsen,1)
