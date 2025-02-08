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

/* EDITING INFORMATION */
/* Please edit only hh52.c, the hh52m.c is created by the command: */
/* sed 's/(hh52m/(hh52mm/g;s/#define GATE 2/#define GATE 2/g;' hh52.c > hh52m.c */

/**
 * TS 2016: IONIC description of the Hodgkin-Huxley 1952 model.
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
#include "ionic.h"
#include "qpp.h"

/* GATE flag to switch between definition of gate variables 
   as other, as tabulated gates or as Markov chains */

/* 0 for definition as other variables */
/* 1 for definition as tabulated gating variables to take 
   advantage of Rush-Larsen method */
/* 2 for definition as Markov chains to use matrix
   Rush-Larsen method */
#define GATE 2

#if (GATE == 2)
/* Enumerate the INa Markov chains variables */
enum {
  #define _RATE(n,i,a,b)
  #define _(n,i) ina_##n,
  #include "hh52m_ina.h"
  #undef _RATE
  #undef _
  NM_ina /* total number of INa Markov variables */
};

/* Enumerate the IK Markov chains variables */
enum {
  #define _RATE(n,i,a,b)
  #define _(n,i) ik_##n,
  #include "hh52m_ik.h"
  #undef _RATE
  #undef _
  NM_ik	/* total number of IK Markov variables */
};

/* Enumerate Markov chains */
enum {
  MC_ina,
  MC_ik,
  NMC
};
#else  /* GATE */
/* number of Markov chain models */
#define NMC 0
#endif /* GATE */

/* possition of membrane voltage in state vector */
static int V_index = 0;

/* Enumerate all dynamic variables */
enum
{
#define _(n,i) var_##n,
#include "hh52_other.h"
#if (GATE == 2)
#define _RATE(a,b,c,d)
#include "hh52m_ina.h"
#include "hh52m_ik.h"
#undef _RATE
#else  /* GATE */
#include "hh52_tgate.h"
#endif /* GATE */
#undef _
  NV	/* total number of variables */
};

/* Enumerate the other (non-gate) variables */
enum
{
#define _(n,i) other_##n,
#include "hh52_other.h"
#if (GATE == 0)
#include "hh52_tgate.h"
#endif /* GATE */
#undef _
  NO	/* total number of other variables */
};

/* there are none of non-tabulated gate variables */
#define NN 0

#if (GATE == 1)
/* Enumerate the gates */
enum
{
#define _(n,i) gate_##n,
#include "hh52_tgate.h"
#undef _
  NT	/* total number of tabulated gate variables */
};

/* Enumerate the tabulated transition rates */
enum
{
#define _(n,i) _alp_##n,
#include "hh52_tgate.h"
#undef _
#define _(n,i) _bet_##n,
#include "hh52_tgate.h"
#undef _
/* there are no other tabulated functions in HH52 */
  NTAB	/* total number of tabulated functions */
};
#else /* GATE */
#define NT 0
#define NTAB 0
#endif /* GATE */

/* The structure containing the parameter values
   for this instance of the model */
typedef struct
{
  /* First go the canonical cell parameters */
#define _(name,default) real name;
#include "hh52_par.h"
#undef _
  /* Then the external current. */
  real IV;
} STR;

/* IONIC_FTAB_HEAD expands to a function defining voltage dependent
   transition rates for tabulation:
   int ftab_hh52(real V, real *values, int ntab)

   INPUT
   -----
   real V       | membrane voltage
   real *values | pointer to array to be filled with Transition
                | rates
   int ntab     | number of tabulated variables

   Returns 1 if succeeds.
*/
IONIC_FTAB_HEAD (hh52m)
{
#if (GATE == 1)
#include "hh52_ftab.h"
  /* Copy the results into the output array values[]. */
#define _(n,i) values[_alp_##n]=alpha_##n;
#include "hh52_tgate.h"
#undef _
#define _(n,i) values[_bet_##n]=beta_##n;
#include "hh52_tgate.h"
#undef _
#endif /* GATE */
} IONIC_FTAB_TAIL;

#if (GATE == 2)
/* CHANNEL_TR_MATRIX expands to a function which fills the transition
   rates matrix:
   int hh52m_{ina,ik}(real *u, real *tr_mat)

   INPUT ARGUMENTS
   ---------------
   real *u      | states vector
   real *tr_mat | matrix to be filled with transtion rates
 */
CHANNEL_TR_MATRIX(hh52m_ina){
  #define V (u[0])
  /* transition rates */
  real alpha_m = 0.1 * (-V + 25.0) /
    (exp ((-V + 25.0) / 10.0) - 1.0);
  real beta_m = 4.0 * exp (-V / 18.0);
  real alpha_h = 0.07 * exp (-V / 20.0);
  real beta_h = 1.0 / (exp ((-V + 30.0) / 10.0) + 1.0);
  #undef V
 
  #define _(n,i)
  #define _RATE(from,to,direct,reverse)		\
    TR_MAT(ina,from,to,direct,reverse)
  #include "hh52m_ina.h"
  #undef _RATE
  #undef _
  return 1;
}

CHANNEL_TR_MATRIX(hh52m_ik){	/* described above */
  #define V (u[0])
  /* transition rates */
  real alpha_n = 0.01 * (-V + 10.0) /
    (exp ((-V + 10.0) / 10.0) - 1.0);
  real beta_n = 0.125 * exp (-V / 80.0);
  #undef V

  #define _(n,i)
  #define _RATE(from,to,direct,reverse)		\
    TR_MAT(ik,from,to,direct,reverse)
  #include "hh52m_ik.h"
  #undef _RATE
  #undef _
  return 1;
}
#endif /* GATE */

/* IONIC_FDDT_HEAD expands to a function of right hand sides for the
   computation of increment of other and non-tabulated gates.

   int fddt_hh52(real *u,int nv,real *values,int ntab,Par par,\
         Var var,real *du,int no,real *nalp,real *nbet,int nn)

   INPUT ARGUMENTS
   ---------------
   real *u      | pointer to array of states variables
   int nv	| total number of variables
   real *values | pointer to array of tab. transition rates
   int ntab     | number of tab. transition rates
   Par par	| parameter structure
   Var var	| variable structure 
   real *du     | pointer to array of increments of *u
   int no	| number of other variables
   real *nalp   | pointer to array of non-tabulated alphas
   real *nbet   | pointer to array of non-tabulated betas
   int nn       | number of non-tab. gates
*/
IONIC_FDDT_HEAD (hh52m, NV, NTAB, NO, NN)
{
  /* Declare the const pars and take their values from struct 
     S==par (a formal parameter) */
#define _(name,default) DEVICE_CONST(real,name);
#include "hh52_par.h"
#undef _
  DEVICE_CONST (real, IV);
  /* Declare and assign local variables for dynamic variables
     from state vector */
  /* ..., first for non-gate variables */
#define _(name,initial) real name=u[var_##name];
#include "hh52_other.h"
#undef _
#if (GATE == 2)
  /* ..., then the Markov chains. */
#define _RATE(a,b,c,d)
#define _(name,i) real name=u[var_##name];
#include "hh52m_ina.h"
#undef _
#define _(name,i) real name=u[var_##name];
#include "hh52m_ik.h"
#undef _
#undef _RATE
#else  /* GATE */
  /* ..., and then for tabulated gate variables */
#define _(name,i) real name=u[var_##name];
#include "hh52_tgate.h"
#undef _
#endif /* GATE */
  /* Calculate the rates of non-gate variables */
#if (GATE == 0)
#include "hh52_ftab.h"
real dot_n = alpha_n * (1 - n) - beta_n * n;
real dot_m = alpha_m * (1 - m) - beta_m * m;
real dot_h = alpha_h * (1 - h) - beta_h * h;
#endif /* GATE */

#if (GATE != 2)
  real O_K = n * n * n * n;
  real O_Na = m * m * m * h;
#endif /* GATE */
/* currents */
real I_K = G_K * O_K * (V - E_K);
real I_Na = G_Na * O_Na * (V - E_Na);
real I_l = G_l * (V - E_l);

/* model equations */
real dot_V = -1.0 / C_m * (I_Na + I_K + I_l);

  /* Copy the calculated rates into the output array du[].  */
  /* Care is taken that all, and only, non-gating variables
     are attended here */
#define _(name,initial) du[other_##name]=dot_##name;
#include "hh52_other.h"
#if (GATE == 0)
#include "hh52_tgate.h"
#endif /* GATE */
#undef _
  /* Finally add the "external current" parameter values */
  du[V_index] += IV;
} IONIC_FDDT_TAIL;

/* IONIC_CREATE_HEAD expands to a function which initialises\
   an instance of the model.

   int create_hh52(ionic_str *I,char *w,real **u,int v0)

   INPUT ARGUMENTS
   ---------------
   ionic_str *I | pointer to ionic structure to be initialised
   char *w      | parameters to be assigned from script
   real **u     | pointer to array of states variables
   int v0       | number of entries in states array
*/
IONIC_CREATE_HEAD (hh52m)
{
  /* Here we assign the parameter values to the structure 
     AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  #include "hh52_par.h"
  #undef _
  
  /* Assign the initial values as given in the *.h files */
  #if (GATE == 2)
  /* 
     the macro SUBCHAIN(fun_tr, index, min, max, incr, sc) intialises
     the parameters of Markov subchain.

     variable | meaning
     --------------------------------------------------
     fun_tr   | function of transition rates as defined 
              | by CHANNEL_TR_MATRIX
     index    | index of control variable for tabulation
              | (negative to avoid tabulation)
     min      | minimal value for tabulation
     max      | maximal value for tabulation
     incr     | increment in the tabulation
     sc       | scale for tabulation
     --------------------------------------------------
  */
  /* intialize INa Markov chain */
  sbch = &(ch->subchain[0]);
  ch->dimension = NM_ina;
  SUBCHAIN(hh52m_ina, -1, -200, 200, 0.01, 0);
  /* intialize IK Markov chains */
  ch += 1;
  sbch = &(ch->subchain[0]);
  ch->dimension = NM_ik;
  SUBCHAIN(hh52m_ik, -1, -200, 200, 0.01, 0);  

  /* asign gates for computation of MCs initial conditions */
  #define _(name, initial) real name = initial;
  #include "hh52_tgate.h"
  #undef _
  #endif /* GATE */

  #define _(name,initial) (*u)[var_##name]=initial;
  #include "hh52_other.h"
  #if (GATE == 2)
  #define _RATE(a,b,c,d)
  #include "hh52m_ina.h"
  #include "hh52m_ik.h"
  #undef _RATE
  #else	 /* GATE */
  #include "hh52_tgate.h"
  #endif /* GATE */
  #undef _
} IONIC_CREATE_TAIL;
