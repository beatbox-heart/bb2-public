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


/* void CONCENTRATIONS(double *T, double *y, double *f) */ {

  /* This procedure computes ion concentrations

It fully incorporates the developments described in
Hilgemann & Noble, Proceedings of the Royal Society, 1987, B230, 163-205.
It also incorporates updates 1.21 (mitochondrial calcium sequestration)
1.3 (the A-Ch activated K current) and 1.4 (T & L Ca channels).

The following is the list of conversions of variables from their names
in the Hilgemann-Noble paper:


               H-N paper                       HEART

Contraction:

       R1 ('Light chain conformation')         Y[29]

       R2 ('cross bridge reactions')           Y[30]

Ca release:

       F1 (precursor fraction)                 Y[26]

       F2 (activator fraction)                 Y[27]

       F3 (product fraction)                   Y[28]

Ca stores:

              Uptake                           Y[12]

              Release                          Y[13]

Ca buffers:

              Ca,cytot                        Y[31]

	      Calmodulin                      Y[20]

              Troponin & SR buff              Y[21]


Mitochondria:

              Mitochondrial Ca                Y[25]

----------------------------------------------------------------------

Updated on 14th April 1991 to include subsarcolemmal space in
preparation for single ventricular cell model


New variables:


              SPVOL              Fraction of VI occupied by space

              VSpace             Space volume using same units as VI

              VSpaceF            VSpace*F

              SRfract            Fraction of SR sensing Ca in space

              ICAfract           Fraction of ICa flowing into space

              INCfract           Fraction of INaCa using space

              TauSPvol           Time constant for exchange with [Ca]i

              y[37]              [Ca]i in space (equivalent of [[4])

              y[38]              SR Ca (equivalent of y[12])

              y[39]              SR release Ca (equivalent of y[13])

	      y[40]              Buffer 1 in space (equivalent to y[20])

              y[41]              Buffer 2 in space (equivalent to y[21])

              y[42]              [Na] in space (equivalent to y[2])

              y[43]              f2 for ICa sensing space Ca

              y[44]              Ca probe in space (equivalent to y[34])

	      y[45]              Total Ca (equivalent of y[31])

              y[46]              Precursor fraction (equivalent of y[26])

              y[47]              Activator fraction (equivalent of y[27])

              y[48]              Product fraction (equivalent of y[28])



calls: no other procedures

===========================================================================*/
  /* local again */
  double _Z1, _Z2, _Z3, _Z4, _Z5;
					  #undef _L
  short _L;

  short FORLIM;


  /* d[Na]i/dt = net sodium current/(F x internal volume)
     assume fraction of ibna carried by Na proportional to [Na]o */

  imNa = iNa + ibNa * Nao / 140 + ifNa + iCaNa + iP * nNaK / (nNaK - 1);
  if (NNaCa > 2)
    imNa += iNaCa * NNaCa / (NNaCa - 2);
  f[1] = -(imNa / ViF);

  /* NOTE: still to decide what to do about NAIS */


  /* add Na change due to perfusion */

  if (Space == 5)
    f[1] += (NaPE - y[1]) / TNaP;


  /* for 3 compartment case set d[K]c/dt = iK/VF - K diffusion */

  if (Space == 3)
    f[2] = imK / VF - PF * (y[2] - Kb);
  else
    f[2] = 0.0;


  /* when ecs is not restricted, set [K]c = K[b] */

  if (Space == 4 || Space == 5)
    y[2] = Kb;


  /* Compute mitochondrial calcium */

  if (Mmode == 1) {
    _Z1 = PoN(KcMup / y[3], NcMup);
    _Z1 = alpha[24] / (1 + _Z1);
    _Z2 = PoN(KnMRel / y[1], NNMRel);
    _Z3 = 1 / (1 + KcMRel / y[24]);
    _Z2 = beta[24] * _Z3 / (1 + _Z2);
    _Z4 = 0.0639 * _Z3 * beta[24];

    _Z1 += -_Z2 - _Z4;
    f[24] = _Z1 / Vmit;
    CaFmit = _Z1 / Vi;

  } else
    f[24] = 0.0;


  /* compute calcium stores */

  if (!CaSS) {
    /* When BMODE=1 use Calgary buffer equations*/

    switch (Bmode) {

    case 0:   /* D-N formulation */
      for (_L = 20; _L <= 23; _L++)
	f[_L - 1] = 0.0;
      break;

    case 1:  /* Calgary formulation */
      f[19] = alpha[19] * y[3] * (1 - y[19]) - beta[19] * y[19];
      f[20] = alpha[20] * y[3] * (1 - y[20]) - beta[20] * y[20];
      f[21] = alpha[21] * Mg * (1 - y[21] - y[22]) - beta[21] * y[21];
      f[22] = alpha[22] * y[3] * (1 - y[21] - y[22]) - beta[22] * y[22];
      break;

    case 2:  /* Hilgemann formulation */
      if (BufFast) {
	f[30] = 0.0;
	y[30] = y[3] + y[19] + y[20];
	f[20] = (CTrop - y[20]) * alpha[20] * y[3] - y[20] * beta[20];
	f[19] = (MTrop - y[19]) * alpha[19] * y[3] - y[19] * beta[19];
	if (Fura != 0) {  /*2.2*/
	  f[33] = (Fura - y[33]) * alpha[33] * y[3] - y[33] * beta[33];
	  y[30] += y[33];
	} else
	  f[33] = 0.0;
	if (SPvol != 0) {
	  /*Assume buffer characteristics same as main space*/
	  f[44] = 0.0;
	  y[44] = y[36] + y[39] + y[40];
	  f[40] = (CTrop - y[40]) * alpha[20] * y[36] - y[40] * beta[20];
	  f[39] = (MTrop - y[39]) * alpha[19] * y[36] - y[39] * beta[19];
	  if (FuraS != 0) {
	    f[43] = (FuraS - y[43]) * alpha[33] * y[36] - y[43] * beta[33];
	    y[44] += y[43];
	  } else
	    f[43] = 0.0;
	}

      } else {
	/*BUFFSS(T, y, f);*/
	/* static void BUFFSS(double *T, double *y, double *f)  */
	/*This procedure was added to version 2.0. It calculates steady state
	  buffers using the Hilgemann formulation.*/
	do {

	  _Z2 = y[20];
	  _Z3 = y[19];
	  _Z4 = y[33];
	  #define OLDVERSION
	  #ifdef OLDVERSION
	  y[20] = CTrop * (y[19] + y[33] - y[30]) /
		  (y[20] - y[30] - CTrop + y[19] + y[33] - beta[20] / alpha[20]);
	  y[19] = MTrop * (y[20] + y[33] - y[30]) /
		  (y[19] - y[30] - MTrop + y[20] + y[33] - beta[19] / alpha[19]);
	  if (Fura != 0) {
	    y[33] = Fura * (y[20] + y[19] - y[30]) /
		    (y[33] - y[30] - Fura + y[20] + y[19] - beta[33] / alpha[33]);
	    _Z5 = fabs((_Z2 - y[20]) / y[20]) + fabs((_Z3 - y[19]) / y[19]) +
		  fabs((_Z4 - y[33]) / y[33]);
	  } else
	    _Z5 = fabs((_Z2 - y[20]) / y[20]) + fabs((_Z3 - y[19]) / y[19]);
	  #else
	  y[20] = CTrop * (_Z3 + _Z4 - y[30]) /
		  (_Z2 - y[30] - CTrop + _Z3 + _Z4 - beta[20] / alpha[20]);
	  y[19] = MTrop * (_Z2 + _Z4 - y[30]) /
		  (_Z3 - y[30] - MTrop + _Z2 + _Z4 - beta[19] / alpha[19]);
	  if (Fura != 0) {
	    y[33] = Fura * (_Z2 + _Z3 - y[30]) /
		    (_Z4 - y[30] - Fura + _Z2 + _Z3 - beta[33] / alpha[33]);
	    _Z5 = fabs((_Z2 - y[20]) / y[20]) + fabs((_Z3 - y[19]) / y[19]) +
		  fabs((_Z4 - y[33]) / y[33]);
	  } else
	    _Z5 = (fabs(_Z2 - y[20]) + fabs(_Z3 - y[19])) /
		  (fabs(y[19]) + fabs(y[20]) + fabs(y[30]) + fabs(y[33]));
	  #endif
	} while (_Z5 >= 0.001);

	y[3] = y[30] - y[19] - y[20] - y[33];

	f[3] = 0.0;
	f[20] = 0.0;
	f[19] = 0.0;
	f[33] = 0.0;

	if (SPvol != 0) {
	  fprintf(Filout3, " Cannot use BUFFSS with sspace \n");
	  Finish = true;
	}

	/*Temporarily assume that alphas and betas and concentrations of buffers
	  are same as in main space*/

	/*NOTE: This does not appear to work!! so commented out*/

	/*      REPEAT

		_Z2 := Y[41];    _Z3 := Y[40];    _Z4 := Y[44];
		Y[41] := (CTROP*(Y[40]+Y[44]-Y[45]))/
		       (Y[41]-Y[45]-CTROP+Y[40]+Y[44]-BETA[21]/ALPHA[21]);
		Y[40] := (MTROP*(Y[41]+Y[44]-Y[45]))/
		       (Y[40]-Y[45]-MTROP+Y[41]+Y[44]-BETA[20]/ALPHA[20]);
		IF FURAS<>0 THEN begin
		       Y[44] := (FURAS*(Y[41]+Y[40]-Y[45]))/
		       (Y[44]-Y[45]-FURAS+Y[41]+Y[40]-BETA[34]/ALPHA[34]);
		       _Z5 := ABS((_Z2-Y[41])/Y[41])+ABS((_Z3-Y[40])/Y[40])+ABS((_Z4-Y[44])/Y[44]);
		end
		ELSE _Z5 := ABS((_Z2-Y[41])/Y[41])+ABS((_Z3-Y[40])/Y[40]);

	      UNTIL _Z5<0.001;

	      Y[37] := Y[45]-Y[40]-Y[41]-Y[44];


	      F[37] := 0;
	      F[41] := 0;
	      F[40] := 0;
	      F[44] := 0;  */

      }  /* of Procedure BuffSS */


      /* Contraction parameters depend on main space */

      CaCHoff = y[3] / (KCaCHoff + y[3]);
      FCtrop = y[20] / CTrop;
      FMtrop = y[19] / MTrop;

      f[21] = 0.0;
      f[22] = 0.0;
      break;


    }/* of BMODE cases */



    switch (SRmode) {

    case 0:  /* SR calcium release and transfer zero */
      for (_L = 12; _L <= 14; _L++) {
	f[_L - 1] = 0.0;
	y[_L - 1] = 0.0;
      }
      for (_L = 26; _L <= 28; _L++) {
	f[_L - 1] = 0.0;
	y[_L - 1] = 0.0;
      }
      iUP = 0.0;
      iRel = 0.0;
      iTran = 0.0;
      break;


    case 1:  /* DiFrancesco-Noble SR model */
      iUP = alpha[11] * y[3] * (Ca12m - y[11]);
      iTran = alpha[12] * y[13] * (y[11] - y[12]);
      /* when rmode is zero repriming is voltage dependent */
      if (Rmode == 0)
	f[13] = Speed[13] * (alpha[13] * (1 - y[13]) - beta[13] * y[13]);
      else {
	y[13] = 1.0;
	f[13] = 0.0;
      }
      _Z4 = 1.0;
      _Z5 = 1.0;
      FORLIM = Nrel;
      for (_L = 1; _L <= FORLIM; _L++) {
	_Z5 *= y[3];
	_Z4 *= KmCa;
      }
      iRel = beta[12] * y[12] * _Z5 / (_Z5 + _Z4);
      f[11] = (iUP - iTran) / (2 * V12F);
      f[12] = (iTran - iRel) / (2 * V13F);
      for (_L = 26; _L <= 30; _L++) {
	f[_L - 1] = 0.0;
	y[_L - 1] = 0.0;
      }
      break;

    case 2:  /* Hilgemann formulation */
      _Z3 = exp(0.08 * (y[0] - 40));   /* Voltage dependence of release */


      /* mimic subsarcolemmal space by artificially decreasing
         KmCa when ica exceeds 0.5 nA */

      if (spmimic && -iCaCa > 0.5 && SPvol == 0)
	_Z1 = y[3] / (y[3] + 0.1 * KmCa);
      else
	_Z1 = y[3] / (y[3] + KmCa);

      _Z1 *= _Z1;   /* Regulatory Ca binding site */
      _Z4 = beta[25] * _Z1 + alpha[25] * _Z3;   /* Activation rate */
      _Z2 = beta[26] * _Z1 + alpha[26];   /* Inactivation rate */

      y[25] = 1 - y[27] - y[26];
      f[25] = y[27] * beta[27] - y[25] * _Z4;   /* Precursor fraction */
      f[26] = y[25] * _Z4 - y[26] * _Z2;   /* Activator fraction */
      f[27] = y[26] * _Z2 - y[27] * beta[27];   /* Product   fraction */
      beta[12] = y[26] / (y[26] + 0.25);
      beta[12] *= beta[12];   /* Open rel channel fraction*/

      iTran = alpha[12] * (y[11] - y[12]);   /* Change in uptake store */


      /* NOTE that all concentration changes to and from the cytosol
         are calculated as concentration changes with respect to the
         cytosol, and related by fractional volumes to the changes of
         other states */


      /* SR PUMP EQUATIONS */

      _Z1 = KCyCa * Kxcs / KSRCa;   /* K1 in paper */
      _Z2 = y[3] + y[11] * _Z1 + KCyCa * Kxcs + KCyCa;   /* K2 in paper */
      Fcyt = y[3] / _Z2;   /* Fr uptake Ca sites */
      FSR = y[11] * _Z1 / _Z2;   /* Fr back SR sites   */
      /* FcytFree = KCyCa / _Z2;   / * Used for display   */
      /* FSRFree = KCyCa * Kxcs / _Z2; */

      iUP = Fcyt * alpha[11] - FSR * beta[11];   /* SR Ca uptake  */

      /*------------------------------------------------------------*/

      f[11] = iUP * Vi / ((1 - SRfract) * SRvol) - iTran;
	  /* Uptake store  */
      iRel = (beta[12] * KmCa2 + SRleak) * y[12];
	  /* Conc change of rel st */

      f[12] = iTran * V12 / V13 - iRel;   /* Release store */

      /* Note: W[8]=10 = 0.2/0.02 = fsrrel/fsrup - See Hilgemann
         PCIN1.PAS */

      if (SPvol != 0) {
	_Z3 = exp(0.08 * (y[0] - 40));   /* Voltage dependence of release */

	_Z1 = y[36] / (y[36] + KmCa);
	_Z1 *= _Z1;   /* Regulatory Ca binding site */
	_Z4 = beta[25] * _Z1 + alpha[25] * _Z3;   /* Activation rate */
	_Z2 = beta[26] * _Z1 + alpha[26];   /* Inactivation rate */

	y[45] = 1 - y[47] - y[46];
	f[45] = y[47] * beta[47] - y[45] * _Z4;   /* Precursor fraction */
	f[46] = y[45] * _Z4 - y[46] * _Z2;   /* Activator fraction */
	f[47] = y[46] * _Z2 - y[47] * beta[27];   /* Product   fraction */
	beta[12] = y[46] / (y[46] + 0.25);
	beta[12] *= beta[12];   /* Open rel channel fraction*/

	ITranS = alpha[12] * (y[37] - y[38]);   /* Change in uptake store */


	/* NOTE that all concentration changes to and from the cytosol
	   are calculated as concentration changes with respect to the
	   cytosol, and related by fractional volumes to the changes of
	   other states */


	/* SR PUMP EQUATIONS */

	_Z1 = KCyCa * Kxcs / KSRCa;   /* K1 in paper */
	_Z2 = y[36] + y[37] * _Z1 + KCyCa * Kxcs + KCyCa;   /* K2 in paper */
	Fcyt = y[36] / _Z2;   /* Fr uptake Ca sites */
	FSR = y[37] * _Z1 / _Z2;   /* Fr back SR sites   */
	/* FcytFree = KCyCa / _Z2;   / * Used for display   */
 /*	FSRFree = KCyCa * Kxcs / _Z2; */

	IupS = Fcyt * alpha[11] - FSR * beta[11];   /* SR Ca uptake  */

	/*------------------------------------------------------------*/

	f[37] = IupS * Vspace / (SRfract * SRvol) - ITranS;
	    /* Uptake store  */
	IrelS = (beta[12] * KmCa2 + SRleak) * y[38];
	    /* Conc change of rel st */

	f[38] = ITranS * V12 / V13 - IrelS;   /* Release store */


      }


      break;



    }/* of SRMODE cases */



    /* Needs improving: this assumes distribution of ICA and INACA that
     merely divides currents according to Fract ratio. They should
     be allowed to be separate currents with their own Ca and Na dependence*/

    ImCa = iCaCa + ibCa + iCa2 + iCa3 - 2 * iNaCa / (NNaCa - 2);

    /* Assume only ICAS and INACAS flow into space */
    ImCaS = ICaCaS - 2 * INaCaS / (NNaCa - 2);


    if (SRmode < 2) {  /* Pre-Hilgemann formulations */
      if (CaBuff)
	f[3] = 0.0;
      else {
	f[3] = (iRel - ImCa - iUP) / (2 * ViF);
	f[30] = 0.0;
	if (Bmode > 0) {
	  /* when BMODE = 1 use Calgary buffer equations to subtract
	     calcium taken up by buffers and to adjust calcium space */

	  _Z1 = Ncalmod * Calmod * f[19] + Nctrop * CTrop * f[20];
	  _Z1 += Nmtrop * MTrop * f[22];
	  f[3] /= CaVOL * BufVol;
	  f[3] -= _Z1;
	}

      }
    } else {  /* Hilgemann formulation */
      if (SLPumpMode) {
	/*DCNP = ImCa;*/
	ImCa += y[3] / (y[3] + KMSLpump) * KSLpump;
	/*DCP = ImCa;*/

	/*NOTE: still to deal with Space SLpump*/

      }


      if (BufFast) {
	f[3] = iRel * (1 - SRfract) * SRvol * V13 / (Vi * V12) - iUP -
	       ImCa / (2 * ViF) - f[20] - f[19] - f[33];
	if (SPvol != 0)
	  f[36] = IrelS * SRfract * SRvol * V13 / (Vspace * V12) - IupS -
		  ImCaS / (2 * VSpaceF) - f[40] - f[39] - f[43];


      } else
	f[30] = iRel * (1 - SRfract) * SRvol * V13 / (Vi * V12) - iUP -
		ImCa / (2 * ViF);

      if (SPvol != 0)
	f[44] = IrelS * SRfract * SRvol * V13 / (Vspace * V12) - IupS -
		ImCaS / (2 * VSpaceF);


      if (ContractMode) {
	f[28] = (1 - y[28]) * FMtrop * FCtrop * FCtrop * alpha[28] -
		y[28] * beta[28];
	f[29] = (1 - y[29]) * y[28] * alpha[29] - y[29] * beta[29];

	/* when sarcomere length is set compute rising limb of
	   length-tension curve as _Z1, falling limb as _Z2, and
	   resting tension as _Z3. Compute isometric tension
	   and store temporarily in y[60] to allow grafy60
	   to set graphics for tension */

	if (SarcoLength != 0) {
	  if (SarcoLength < 1)
	    _Z1 = 0.0;
	  else
	    _Z1 = 1 - exp(-3 * (SarcoLength - 1));
	  if (SarcoLength < 2)
	    _Z2 = 1.0;
	  else
	    _Z2 = 1 - 0.625 * (SarcoLength - 2);
	  _Z3 = 0.0002 * exp(2 * SarcoLength);
	  CBdenA = _Z1 * _Z2 * CBden;
	  /* IsoTen = CBdenA * y[29] + _Z3; */
	  /* y[59] = IsoTen; */
	  /* f[59] = 0.00000000000001; */
	      /* forces y[60] to be included in neqnn */
	}


      } else {
	f[28] = 0.0;
	f[29] = 0.0;
      }

    }


    if (Mmode == 1)
      f[3] -= CaFmit;


  }



  /* compute extracellular calcium change */

  if (CaOmode == 1) {
    f[23] = 0.5 * ImCa / VF - DiffCa * (y[23] - CaB);
    Cao = y[23];
  } else
    f[23] = 0.0;

  /* compute change in intracellular potassium */

  f[10] = -(imK / ViF);


  /* compute Ca & K changes due to perfusion */

  if (Space == 5) {
    if (SRmode < 2 || BufFast)
      f[3] += (CaPE - y[3]) / TCaP;
    else
      f[30] += (CaPE - y[3]) / TCaP;
    f[10] += (KPE - y[10]) / TKP;
  }

  /* compute change in ATP concentration */


  if (ComputeATP) {
    f[34] = -((iP + 0.25 * iUP) / ViF);
    if (ContractMode) {
      if (SarcoLength != 0)
	_Z1 = CBdenA;
      else
	_Z1 = CBden;
      f[34] -= y[29] * _Z1 / CBturn;
    }
  }


  /* Compute diffusion between cytosol spaces */

  if (SPvol != 0) {
    _Z1 = (y[36] - y[3]) / TauSPvol;
    f[36] -= _Z1;
    f[3] += _Z1 / SPvol;
  }

  /* Compute integrated calcium transfer by iCa (y[49) and INaCa (y[50]
     in terms of the change that would occur in cytosol calcium
     if distributed uniformly and not buffered */

  if (ComputeCaIntegral) {
    f[48] = -0.5 * iCa / ViF;
    f[49] = iNaCa / ViF;
  }

}  /* of procedure Concentrations */

/* End. */

/* #endif QHEART */
