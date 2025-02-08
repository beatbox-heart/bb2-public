/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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

#include <math.h>
#include "beatbox.h"
#include "p2c.h"
#define OWN
#include "noble.h"
#include "nobleini.h"
#include "nobleftn.h"
					#undef _L


/*  Na/Ca exchange model from Hilgemann et al., 1991
    To use this procedure, call it with the appropriate 4 ion concentrations,
    apparent in the PROCEDURE line, followed by a scalor value for output, 'KNACA.'
    A value of '1' for Knaca will give a maximum inward exchange current of
    about 100 units. A value of 0.1 is realistic for a guinea pig
    myocyte of 120 pF, if simulations are done in nA.  This gives 12 æA/æF
    outward current at 0 mV with saturating cytoplasmic Ca and extracellular Na
*/

void INACADH(double ni, double no, double ci, double co, double knaca,
		    double *inaca, /*double *t,*/ double *y) {
  double xci, xco;
  double kn, k1, k2, k3, k4, cosurf, nosurf, kmc, knc, fasymi, fasymo,
	 rnet, femo, x1, x2, x3, x4, d, kmn, km1, km2, ko, kem2, km3;
  double kem, x3ni, x3no, k5, k6, k7, k8, kci, kco /*,e3, e4*/;
  femo = 0.0322;    /*fractional drop of Em to extracellular binding sites*/
  kco = 415.0;	/*transition rates between Eoc & E2c*/
  kci = 517.0;	/*transition rates between Eoc & E1c*/
  kn = 100.0;	/*Em-dependent transition rates between Eonnn & E2nnn*/
  ko = 184.4;	/*transition rates between Eonnn & E1nnn*/
  kmc = 0.031;	/*dissociation constant for cytoplasmic Ca in mM*/
  km1 = 32.23;	/*dissociation constants for 1st, 2nd and 3rd cytoplasmic Na binding in mM*/
  km2 = 11.11;
  km3 = 11.89;  
  fasymo = 328.6;   /*fractional availability of empty binding sites on
		    the extracellular side, allowing low affinity binding of BOTH
		    the first sodium ion and calcium on extracellular side*/
  fasymi = 0.32;    /*fractional availability of empty binding sites on cytoplasmic side*/

  /*=====================end of parameter descriptions=======================*/

  knc = km3;    /*one sodium can bind noncompetetively on cytoplasmic
		side reasons for this assumption are described in NY
		Academy of Sci article, in press*/
  kem = exp(y[0] * 0.5 * (1 - femo) / RTonF);
  kem2 = exp(-(y[0] * femo / RTonF));
  kmn = km1 * km2 * km3;
  cosurf = co * kem2 * kem2;
  nosurf = no * kem2;
  d = knc * (ci * kmn + kmc * (kmn * (1 + fasymi) + km2 * km3 * ni +
			       km3 * ni * ni + ni * ni * ni)) + ci * ni * kmn;
  x3ni = knc * kmc * ni * ni * ni / d;
  xci = ci * knc * km1 * km2 * km3 / d;
  d=cosurf*kmn+kmc*(kmn*(1+fasymo)+km2*km3*nosurf+km3*nosurf*nosurf+nosurf*nosurf*nosurf);
  x3no = kmc * nosurf * nosurf * nosurf / d;
  xco = cosurf * kmn / d;
  k6 = xci * kci;
  k5 = kci;
  k3 = xco * kco;
  k4 = kco;
  k1 = kn * kem;
  k2 = x3no * kn / kem;
  k7 = ko * x3ni;
  k8 = ko;
  x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
  x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
  x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
  x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
  d = x1 + x2 + x3 + x4;
  rnet = (k1 * k3 * k5 * k7 - k2 * k4 * k6 * k8) / d;
  *inaca = rnet * knaca;
}


static void Yset(short n, /*double *T, */double *y, double *f, double *alpha, double *beta) {
  y[n - 1] = alpha[n - 1] / (alpha[n - 1] + beta[n - 1]);
  f[n - 1] = 0.0;
}  /* of procedure Yset */

static void RateSet(short n, /*double *T,*/ double *y, double *f, double *alpha, double *beta) {
  if (Speed[n - 1] >= Infinite || Slow)
    Yset(n,/* T, */y, f, alpha, beta);
  else
    f[n - 1] = Speed[n - 1] *
	       (alpha[n - 1] - y[n - 1] * (alpha[n - 1] + beta[n - 1]));

}  /* of procedure RateSet */


/* controls use of currents, rates and concentrations */
/*========================================================================*/
// #include "used.h"
// #include "warnpar.h"
void FTN(/*double *T, */double *y, double *f, double *alpha, double *beta) {

  /*       this procedure is the driver procedure of a set of procedures
	   (Rates, Currents, Concentrations) that compute alphas and betas,
	   ionic currents, rates of change of ion concentrations and gating
	   variables. The variables are Y[1] - Y[neqn], their derivatives,
	   F[1] - F[neqn]. Voltage shifts of gating reactions are Shift[1 - neqn]
	   and speeds of reactions are scaled by Speed[1 - neqn]
	   This procedure and its set of procedures are called with the
	   variable T explicitly available; currently none of the procedures
	   actually uses T. They would only need to do so if in the future,
	   functions are introduced with explict dependence on time.

calls: Rates, Currents, Concentrations

=========================================================================*/
  #define __AA(type, name) type name;
  #include "nobleloc.h"
  #undef __AA

  short l;

			Fast=0; /* Originally was set automatically by DESOLVE/ADAMS.. */

  for (l = 0; l < neqn; l++)
    f[l] = 0.0;


  /* Call Rates except when instantaneous IV plots are being computed */

  #include "nobleft3.h"

  /* In all cases, compute ionic currents and set voltage clamp if required */

  #include "nobleft1.h"

  /* compute changes in ion concentrations except when when ivplots are being
     computed. Also compute when T=TSTART to obtain values for imna and imca
     for procedures naiss and caiss*/

  #include "nobleft2.h"

}  
