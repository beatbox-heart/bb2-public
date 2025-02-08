/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
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


/* void RATES(double *T, double *y, double *f) */ {

  /*=========================================================================

	  This procedure computes the rate coefficients for the gating
	  variables for HEART.

	  This version  includes Hilgemann formulation of Ca current
	  inactivation (series V and Ca instead of parallel formulation of
	  the DiFrancesco-Noble model).

	  Update 1.4 incorporated separate T & L channels for the "fast"
	  calcium current.

	  Note: the constants delta3 (0.0001), delta4 (0.00001) and
	  delta6 (0.0000001) are defined as such in initial.pas

calls: no other procedures

===========================================================================*/
  double Z1, /*Z2,*/ Z3, /*Z4,*/ E0, FSS;
  short L;





  /* No need to re-compute V-dependent rates if V constant, unless
     IFLAG is set to 0 to force a computation */

  /*if (y[0] != Esav || Iflag == 0 ) */{
    switch (Ymode) {   /* ymode determines rates for i(f) equations */

    case 0:
      /* y rate equations as in MNT model with alpha/beta interchanged */

      E0 = y[0] + 52 - Shift[4];
      alpha[4] = 0.05 * exp(-0.067 * E0);
      if (fabs(E0) < delta4)
	beta[4] = 2.5;
      else
	beta[4] = E0 / (1 - exp(-0.2 * E0));

      break;


    /* ymode 1 is DiFrancesco's theory curve for y */

    case 1:
    case 10:
      E0 = y[0] - Shift[4];
      alpha[4] = 0.085 * exp(E0 / -15);
      beta[4] = 6.6 * exp(E0 / 25);

      break;


    /*ymode 2 is DiFrancesco's experimental result for y*/

    case 2:
    case 20:
      E0 = y[0] - Shift[4];
      alpha[4] = 0.014 * exp(E0 / -16);
      beta[4] = 9.75 * exp(E0 / 19);

      break;


    }/*of ymode cases*/




    switch (Kmode) {   /* Kmode determines iK equations */

    /* x rate equations as in MNT model when kmode = 0 */

    case 0:
      E0 = y[0] + 50 - Shift[5];
      alpha[5] = 0.5 * exp(0.0826 * E0) / (1 + exp(0.057 * E0));
      E0 = y[0] + 20 - Shift[5];
      beta[5] = 1.3 * exp(-0.06 * E0) / (1 + exp(-0.04 * E0));
      break;


    /* x rates  based  on DiFrancesco,  Noma  & Trautwein  and Oxford
             sinus  results  when  kmode = 1 */

    case 1:
      E0 = y[0] + 22 - Shift[5];
      if (fabs(E0) < delta3) {
	alpha[5] = 2.5;
	beta[5] = 2.5;
      } else {
	alpha[5] = 0.5 * E0 / (1 - exp(E0 / -5));
	beta[5] = 0.357 * E0 / (exp(E0 / 7) - 1);
      }
      break;


    /* x rates slower in diastolic range if kmode = 2 */

    case 2:
      E0 = y[0] + 22 - Shift[5];
      if (fabs(E0) < delta3) {
	alpha[5] = 2.5;
	beta[5] = 2.5;
      } else {
	alpha[5] = 0.5 * E0 / (1 - exp(E0 / -5));
	beta[5] = 0.178 * E0 / (exp(E0 / 15) - 1);
      }
      break;


    /* Galveston rate functions used when kmode = 3 */

    case 3:
      E0 = y[0] + 26.488 - Shift[5];
      if (fabs(E0) < delta3)
	alpha[5] = 0.1072;
      else
	alpha[5] = 0.01344 * E0 / (1 - exp(-0.124 * E0));
      beta[5] = 0.2063 * exp(-0.039268 * E0);
      break;


    case 4:
      E0 = y[0] - Shift[5];
      alpha[5] = 2.1 * exp(E0 / 28);
      beta[5] = 0.96 * exp(E0 / -24);

      break;


    case 5:
      E0 = y[0] + 45 - Shift[5];
      alpha[5] = exp(0.0826 * E0) / (1 + exp(0.057 * E0));
      alpha[35] = 0.2 * alpha[5];
      alpha[5] = 0.7 * alpha[5];
      E0 = y[0] - 10 - Shift[5];
      beta[5] = exp(-0.06 * E0) / (1 + exp(-0.04 * E0));
      beta[35] = 1.3 * beta[5];
      beta[5] = 2.8 * beta[5];

      break;


    case 6:
      E0 = y[0] - 20 - Shift[5];
      if (E0 == 0)
	E0 = delta4;
      alpha[5] = 0.4 * E0 / (1 - exp(E0 / -13));
      E0 = y[0] + 80 - Shift[5];
      if (E0 == 0)
	E0 = delta4;
      beta[5] = 1.7 * E0 / (exp(E0 / 15) - 1);

      break;


    }/* of KMode cases */




    switch (Camode) {

    /* CaMode determines iCa equations :

       CaMode  activation (d)  inactivation (f)  Ca inact   slow inact
         1      Oxford           Oxford
         2      Oxford           Oxford                       y[15]
         3      Galveston        Oxford
         4      Galveston        Oxford                       y[15]
         5      Oxford           Oxford            y[15]
         6      Galveston        Oxford            y[15]
         7      Calgary          Rasmussen
         8      Oxford           -----  Hilgemann ------
         9      Oxford           -----  Hilgemann ------      y[15]
         10     Oxford           Oxford            y[15]      y[53..56]

     */

    /* if camode = 1, 2, 5, 8, 9 or 10 use 'Oxford' sinus d & f */

    case 1:
    case 2:
    case 5:
    case 8:
    case 9:
    case 10:
      E0 = y[0] + 24 - Shift[6];
      if (fabs(E0) < delta3)
	alpha[6] = 120.0;
      else
	alpha[6] = 30 * E0 / (1 - exp(E0 / -4));
      if (fabs(E0) < delta3)
	beta[6] = 120.0;
      else
	beta[6] = 12 * E0 / (exp(E0 / 10) - 1);
      break;

    /* if camode = 3, 4 or 6 use Galveston d equations */

    case 3:
    case 4:
    case 6:
      E0 = y[0] + 12.07 - Shift[6];
      if (fabs(E0) < delta3)
	alpha[6] = 10.7;
      else
	alpha[6] = 18.1 * E0 / (1 - exp(-0.169 * E0));
      beta[6] = 40.18 * exp(-0.0459 * E0);
      break;

    /* if CaMode=7 use Calgary d equations */

    case 7:
      E0 = y[0] + 9.1 - Shift[6];
      alpha[6] = 10 * E0 / (1 - exp(E0 / -8.23));
      beta[6] = alpha[6] * exp(E0 / -8.23);
      break;

    }/* of CaMode cases for activation */


    switch (Camode) {

    /* in all camodes except 7, 8 and 9 use 'oxford' f equations */


    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 10:
      E0 = y[0] + 34 - Shift[7];
      if (fabs(E0) < delta3)
	alpha[7] = 25.0;
      else
	alpha[7] = 6.25 * E0 / (exp(E0 / 4) - 1);
      beta[7] = 50 / (1 + exp(E0 / -4));
      break;

    case 7:
      /* in CaMode=7 use Randy's f formulation */

      E0 = 0.03368 * (y[0] + 10);
      E0 *= E0;
      alpha[7] = 19.6978 * exp(-E0) + 19.78848;
      FSS = (1 - 0.01) / (1 + exp((y[0] - -28.05) / 8.58));
      FSS += 0.01 + 1 / (1 + exp((50 - y[0]) / 20));
      beta[7] = 0.0;
      break;


    /* in CaModes 8 or 9 use Hilgemann formulation */

    case 8:
    case 9:
      E0 = y[0] + 34 - Shift[7];
      if (fabs(E0) < delta3)
	alpha[7] = 25.0;
      else
	alpha[7] = 6.25 * E0 / (exp(E0 / 4) - 1);
      beta[7] = 12 / (1 + exp(E0 / -4));

      /* Note: rest of Hilgemann formulation left until
         computation of F[8]*/

      break;


    }/* of CaMode cases for inactivation */



    /* in camode 2, 4, 9  compute slow ca inactivation */

    if (Camode == 2 || Camode == 4 || Camode == 9) {
      E0 = y[0] + 34 - Shift[14];
      if (fabs(E0) < delta3)
	alpha[14] = 2.5;
      else
	alpha[14] = 0.625 * E0 / (exp(E0 / 4) - 1);
      beta[14] = 5 / (1 + exp(E0 / -4));
    }


    /* CaMode 10 must use y[53] since y[15] used for Ca-induced inact */

    if (Camode == 10) {
      E0 = y[0] - Shift[52];
      alpha[52] = 0.0164 * exp(E0 / -18.4) + 0.0746;
      beta[52] = 0.113 / (1 + exp((E0 + 17.9) / -26.1));
      alpha[53] = alpha[52];
      alpha[54] = alpha[53];
      alpha[55] = alpha[54];
      beta[53] = beta[52];
      beta[54] = beta[53];
      beta[55] = beta[54];
    }





    /* ca store repriming equations based on those for f for voltage dep.*/

    E0 = y[0] + 34 - Shift[13];
    if (fabs(E0) < delta3)
      alpha[13] = 2.5;
    else
      alpha[13] = 0.625 * E0 / (exp(E0 / 4) - 1);
    beta[13] = 5 / (1 + exp(E0 / -4));





    /* TO mode determines equations for iTO */

    /* Ito inactivation based on Dario's modification of Fozzard &
              Hiraoka (J.Physiol, 1973) */

    if (TOmode != 0) {
      E0 = y[0] - Shift[17];
      alpha[17] = 0.033 * exp(E0 / -17);
      beta[17] = 33 / (1 + exp((E0 + 10) / -8));

      /* Ito activation based on MNT model equations */

      if (TOmode == 2) {
	E0 = y[0] - Shift[18];
	if (fabs(E0) < delta3)
	  alpha[18] = 80.0;
	else
	  alpha[18] = 8 * E0 / (1 - exp(-0.1 * E0));
	beta[18] = 80 / (0.1 + exp(0.0888 * E0));
      }
    }




    /* slow ca current: KS Lee current */

    if (PCa3 != 0) {
      E0 = y[0] + 55 - Shift[15];
      if (fabs(E0) < delta3)
	alpha[15] = 5.0;
      else
	alpha[15] = 1.25 * E0 / (1 - exp(E0 / -4));
      if (fabs(E0) < delta3)
	beta[15] = 5.0;
      else
	beta[15] = 0.5 * E0 / (exp(E0 / 10) - 1);
      E0 = y[0] + 45 - Shift[16];
      if (fabs(E0) < delta3)
	alpha[16] = 0.84;
      else
	alpha[16] = 0.21 * E0 / (exp(E0 / 4) - 1);
      beta[16] = 1.7 / (1 + exp(E0 / -4));
    }




    /* Low threshold T channel */

    /* provisionally, this channel is given the same equations as for
       the L channel (ICA). You should use shift32, shift33, speed32 and
       speed33 to shift the activation and inactivation curves and to
       increase the rates. When more information is available, these
       equations will be updated */

    if (PCa2 != 0) {
      E0 = y[0] + 24 - Shift[31];
      if (fabs(E0) < delta3)
	alpha[31] = 120.0;
      else
	alpha[31] = 30 * E0 / (1 - exp(E0 / -4));
      if (fabs(E0) < delta3)
	beta[31] = 120.0;
      else
	beta[31] = 12 * E0 / (exp(E0 / 10) - 1);
      E0 = y[0] + 34 - Shift[32];
      if (fabs(E0) < delta3)
	alpha[32] = 25.0;
      else
	alpha[32] = 6.25 * E0 / (exp(E0 / 4) - 1);
      beta[32] = 50 / (1 + exp(E0 / -4));
    }





    /* For NaMode=1  m & h equations based on Colatsky & Brown, Lee & Powell */

    if (Namode == 1) {  /* Namode determines i Na equations */
      E0 = y[0] + 41 - Shift[9];
      if (fabs(E0) < delta4)
	alpha[9] = 2000.0;
      else
	alpha[9] = 200 * E0 / (1 - exp(-0.1 * E0));
      beta[9] = 8000 * exp(-0.056 * (y[0] + 66 - Shift[9]));
      alpha[8] = 20 * exp(-0.125 * (y[0] + 75 - Shift[8]));
      beta[8] = 2000 / (320 * exp(-0.1 * (y[0] + 75 - Shift[8])) + 1);
      /* Correction due to gB fraction of persistently open h subunits
       * in LQT3 syndrome - vnb 95/12/16 */
      alpha[8] += beta[8]*gB;
      beta[8]  -= beta[8]*gB;
    }


    /* For NaMode=2 use Calgary atrial rate equations */

    if (Namode == 2) {
      E0 = y[0] - 31.23 - Shift[9];
      alpha[9] = 757.22 * E0 / (1 - exp(-0.0617 * E0));
      E0 = y[0] + 52.853 - Shift[9];
      beta[9] = 2256.5 * exp(-0.0418 * E0);
      E0 = y[0] + 85.45 - Shift[8];
      alpha[8] = 51.1 * exp(-0.1235 * E0);
      E0 = y[0] + 2.91 - Shift[8];
      beta[8] = 1270.7 / (2.5 + exp(-0.0764 * E0));
    }



    /* Mark Boyett's rate equations for relaxation of iach */

    if (AChgKmax != 0) {
      beta[50] = 5.82 / (1 + exp((y[0] + Shift[50]) / -15));
      beta[51] = 120 / (1 + exp((y[0] + Shift[51]) / -15));
    }

  }  /* of computation of gate rate coefficients */





  /* ------------------------------------------------------------------
      in camodes 5, 6 and 10 inactivation is ca induced using equation:

          df2/dt = alpha[15](1 - f2)   -   beta[15]*cai*f2

      The equations use the relations:

          alpha[15] = 1/TAUF2

      Where TAUF2 is the limiting recovery time constant when [Ca]i->0

          beta[15] = alpha[15]/KMINACT

      where KMINACT is [Ca]i for half inactivation in steady state.

      For convenience in programming the actual beta[15] in the program
      is [Ca]i*beta[15]. This enables the same differential equation to
      be used later in the program regardless of the value of CaMode

  ---------------------------------------------------------------------*/

  if (Camode == 5 || Camode == 6 || Camode == 10) {
    beta[14] = y[3] * alpha[14] / Kminact;
    if (SPvol != 0)
      beta[42] = y[36] * alpha[14] / Kminact;

    /*Use same alpha[15] for Space ICa*/

  }


  /* slow code now omitted: see slow.pas */


  /* compute dy/dt,  dx/dt */

  RateSet(5,/* T,*/ y, f, alpha, beta);
  RateSet(6,/* T,*/ y, f, alpha, beta);
  if (Kmode == 5) {
    Speed[35] = Speed[5];
    RateSet(36,/* T,*/ y, f, alpha, beta);
  }

  /* compute dr/dt and dq/dt */

  if (TOmode != 0) {
    RateSet(18, /*T,*/ y, f, alpha, beta);
    if (TOmode != 5) {
      if (TOmode == 2)
	RateSet(19,/* T,*/ y, f, alpha, beta);
      else {
	f[18] = 0.0;
	y[18] = 1.0;
      }
    }
  } else {
    f[17] = 0.0;
    f[18] = 0.0;
    y[17] = 1.0;
    y[18] = 1.0;
  }

  /* compute dd/dt, df/dt, & dh/dt */

  RateSet(7, /*T,*/ y, f, alpha, beta);

  switch (Camode) {

  /* When CaMode<7 use DiFrancesco-Noble formulation */

  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 10:
    RateSet(8,/* T, */y, f, alpha, beta);
    break;

  /* When CaMode=7 use Randy's formulation */

  case 7:
    f[7] = alpha[7] * (FSS - y[7]);
    break;

  /*When CaMode=8 or 9 use Hilgemann formulation */

  case 8:
  case 9:
	/* !!! This line inserted to eliminated dependence on previous calls !!! VNB*/
			CaCHoff = y[3] / (KCaCHoff + y[3]);

    Z1 = (1 - y[7]) * (1 - CaCHoff);   /* -Ca/not V-inactivated */
    /*Z2 = y[7] * (1 - CaCHoff);   / * -Ca/    V-inactivated */
    Z3 = (1 - y[7]) * CaCHoff;   /* +Ca/not V-inactivated */
    /*Z4 = y[7] * CaCHoff;   / * +Ca/    V-inactivated */
    f[7] = Speed[7] * (120 * Z3 + Z1) * beta[7] - y[7] * alpha[7];

    break;

  }/* of CaMode cases for F[8] */

  if (Camode == 2 || Camode == 4 || Camode == 5 || Camode == 6 ||
      Camode == 9 || Camode == 10)
    RateSet(15,/* T, */y, f, alpha, beta);
  else
    f[14] = 0.0;

  for (L = 53; L <= 56; L++) {
    if (Camode == 10)
      RateSet(L,/* T, */y, f, alpha, beta);
    else
      f[L - 1] = 0.0;
  }

  /* compute dd3/dt and df3/dt */

  if (PCa3 != 0) {
    RateSet(16,/* T, */y, f, alpha, beta);
    RateSet(17,/* T, */y, f, alpha, beta);
  } else {
    f[15] = 0.0;
    f[16] = 0.0;
  }

  /* compute d and f changes for T channel */

  if (PCa2 != 0) {
    RateSet(32,/* T, */y, f, alpha, beta);
    RateSet(33,/* T, */y, f, alpha, beta);
  } else {
    f[31] = 0.0;
    f[32] = 0.0;
  }

  /* test whether ina at steady state */

  if (GNaSS)
    Yset(9,/* T, */y, f, alpha, beta);
  else
    RateSet(9,/* T, */y, f, alpha, beta);

  /* compute dm/dt when fast is true, else use steady state value */

  #if 0
  if (Fast)
    RateSet(10,/* T, */y, f, alpha, beta);
  else
  #endif
    Yset(10,/* T, */y, f, alpha, beta);


  /* Mark Boyett's iAch rates */

  if (AChgKmax != 0) {
    RateSet(51,/* T, */y, f, alpha, beta);
    RateSet(52,/* T, */y, f, alpha, beta);
  } else {
    f[50] = 0.0;
    f[51] = 0.0;
  }

}  /* of procedure RATES */

/* End. */
/* #endif QHEART */
