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


#include "extern.h"
/* void CURRENTS(double *T, double *y, double *f) */ {

  /* this procedure computes ionic currents for HEART

note: the constants delta3, delta4 and delta6 are as described in
rates.pas

calls: no other procedures*/
  #define EXTERN
  /* !!!  came from initial.h */
  EXTERN Ardsub EKC;   /* cleft values of EK */
  /*EXTERN Ardsub IMKSAV;   / * saved values of imK */
  EXTERN Ardsub IMKC;   /* cleft values of imK */
  EXTERN Ardsub IKC;   /* cleft values of iK */
  EXTERN Ardsub IK1C;   /* cleft values of iK1 */
  EXTERN Ardsub IPC;   /* cleft values of ip */
/*  EXTERN Ardsub KCSAV;   / * saved values of [K] cleft */
  EXTERN Ardsub ITOC;   /* cleft values of iTO */
  EXTERN Ardsub IFKC;   /* cleft values of K component of i(f) */
  EXTERN Ardsub ICAKC;   /* cleft values of K component of iCa */
  /*  EXTERN Ardsub YSAV; */
  EXTERN Ardsub IFNAC;   /* cleft values of Na component of i(f) */
  EXTERN Ardsub IBKC;   /* cleft values of background K current */
  /*EXTERN Ardsub IACHC;  / * cleft values of iACh */
  #undef EXTERN
  /* previously local variables of parent procedure, now local again */
  short _L;
  double _E0, _Z1, _Z2, _Z3, _Z4, _Z5;

  short FORLIM;
  /* used in integrating cleft currents dependent on [K] */
  #define IntRad(r) (((r)[_L - 1] + (r)[_L]) * _Z1)

  /* CALCIUM CURRENT INACTIVATION */

  /* Hilgemann formulation for i Ca inactivation requires this calculation
     for subsequent use */

  if (Camode == 8 || Camode == 9) {
    CaCHoff = y[3] / (KCaCHoff + y[3]);
    CaCHon = (1 - y[7]) * (1 - CaCHoff);
  } else
    CaCHon = y[7];

  /* when [ACh] not zero incorporate ACh inhibition of calcium channel */

  if (ACh != 0)
    CaCHon *= 1 - ACh / (AChRecCKm + ACh);



  /* SODIUM CURRENT */

  /* assume 12% K leak in sodium (m3h) channel */

  ENa = -RTonF * log(y[1] / Nao);
  /*        				- K was a full substitute for y[2]
  if (Space > 2)                          though K used at Space<=2
    _Z1 = y[2];                           and y[2] used at Space>2
  else                                    What for? - VNB
    _Z1 = K;
  */
  _Z1 = K = y[2];

  Emh = RTonF * log((Nao + PNaK * _Z1) / (y[1] + PNaK * y[10]));
  iNa = gNa * y[8] * y[9] * y[9] * y[9] * (y[0] - Emh);



  /* BACKGROUND CURRENT */

  /* background inward current ---- note that the background channel
  may either have separate values for GBNA & GBK or one may assume equal Na
  and K conductance (set GB only) as in Colqhoun,Neher,Reuter & Stevens */

  ibNa = gbNa * (y[0] - ENa);

  /* SODIUM-CALCIUM EXCHANGE CURRENT */

  /*
           for details of formulation of equations for inaca see Mullins,L.
           1977  J.gen.physiol. 70, 681-695
	   1981  Ion Transport in Heart. Raven press
           and
           Noble, D. (1986) In Cardiac Muscle: Regulation of Exitation
           and Contraction ( Ed R.D. Nathan), Academic Press.

           Note that the procedure will
           hold inaca constant if ivplot is not instantaneous to eliminate
           large transient inaca at very negative potentials */


  if (y[3] < 0.00000001)
    y[3] = 0.00000001;
		/* This line added to remove dependence on previous step VNB */
					Cao=y[23];
  ECa = -0.5 * RTonF * log(y[3] / Cao);
  if (y[1] < delta4)
    y[1] = delta4;
  if (NNaCa == 2)
    iNaCa = 0.0;
  else {   /*if nnaca..else*/
    ENaCa = (NNaCa * ENa - 2 * ECa) / (NNaCa - 2);
    _Z1 = RTonF * 2 / (NNaCa - 2);

    /* when Nmode is 0 use simple sinh formulation for inaca */
    /* Strongly recommended not to use this mode */

    if (Nmode == 0)
      iNaCa = kNaCa * (exp((y[0] - ENaCa) / _Z1) - exp((ENaCa - y[0]) / _Z1)) / 2;

    /* when Nmode is 1 use Ca and Na activated inaca */

    if (Nmode == 1) {   /*nmode=1*/
      _Z4 = 1.0;
      _Z5 = 1.0;
      FORLIM = NNaCa;
      for (_L = 1; _L <= FORLIM; _L++) {   /*for _L*/
	_Z4 *= Nao;
	_Z5 *= y[1];
      }
      _Z2 = 1 + DNaCa * (y[3] * _Z4 + Cao * _Z5);
      _Z3 = _Z5 * Cao * exp(2 * yNaCa * y[0] / _Z1);
      _Z3 -= _Z4 * y[3] * exp(-2 * (1 - yNaCa) * y[0] / _Z1);
      iNaCa = (1 - iNCfract) * kNaCa * _Z3 / _Z2;
    }
    if (fabs(iNaCa) > iNaCam) {   /*abs(inaca)>inacam*/
      if (iNaCa > 0)
	iNaCa = iNaCam;
      else
	iNaCa = -iNaCam;
    }

    if (Nmode == 2) {   /*nmode=2*/
      INACADH(y[1], Nao, y[3], Cao, kNaCa, &iNaCa, /*T, */y);
      iNaCa = (1 - iNCfract) * iNaCa;
    }

    if (kNaCa < 0)
      iNaCa = 0.0;
    if (KMNaCa != 0) {   /*kmnaca>*/
      _Z1 = y[3] / (y[3] + KMNaCa);
      FORLIM = KPNaCa;
      for (_L = 1; _L <= FORLIM; _L++)   /*for _L*/
	iNaCa *= _Z1;
    }

    if (SPvol != 0) {
      /*ECaS = -0.5 * RTonF * log(y[3] / Cao);*/
      /*ENaCas = (NNaCa * ENa - 2 * ECaS) / (NNaCa - 2);   / *Assume same ENA*/
      _Z1 = RTonF * 2 / (NNaCa - 2);
      if (Nmode == 1) {   /*if nmode=1*/
	_Z4 = 1.0;
	_Z5 = 1.0;
	FORLIM = NNaCa;
	for (_L = 1; _L <= FORLIM; _L++) {   /*for _L*/
	  _Z4 *= Nao;
	  _Z5 *= y[1];
	}
	_Z2 = 1 + DNaCa * (y[36] * _Z4 + Cao * _Z5);
	_Z3 = _Z5 * Cao * exp(2 * yNaCa * y[0] / _Z1);
	_Z3 -= _Z4 * y[36] * exp(-2 * (1 - yNaCa) * y[0] / _Z1);
	INaCaS = iNCfract * kNaCa * _Z3 / _Z2;
      }
      if (Nmode == 2) {   /*if nmode=1*/
	INACADH(y[1], Nao, y[3], Cao, kNaCa, &INaCaS, /*T,*/ y);
	INaCaS = iNCfract * INaCaS;
      }
      if (KMNaCa != 0) {   /*if kmnaca>*/
	_Z1 = y[36] / (y[36] + KMNaCa);
	FORLIM = KPNaCa;
	for (_L = 1; _L <= FORLIM; _L++)   /*for _L*/
	  INaCaS *= _Z1;
      }

    }  /*svpol<>0*/
    else
      INaCaS = 0.0;
  }

  /* BACKGROUND CALCIUM CURRENT */

  ibCa = gbCa * (y[0] - ECa);

  /* K AND K-DEPENDENT CURRENTS */

  /* compute K dependent currents for 3 compartment model and for non-restricted ecs */

  if (Space > 2) {
    if (y[2] < delta4)
      y[2] = delta4;
    EK = RTonF * log(y[2] / y[10]);

    /* IK1 rectifier described by blocking particle model with
    K activation.  K activation uses Michaelis Menten - cf Sakmann & Trube */

    _E0 = y[0] - EK + 10 - ShiftK1;
    iK1 = gK1 * (y[2] / (KmK1 + y[2])) * (y[0] - EK) /
	  (1 + exp(steepK1 * _E0 / RTonF));

    switch (TOmode) {

    case 0:
    case 1:
      /* When TOmode=0,1 ito described by outward rectifier with instantaneous
	    voltage activation */

      _E0 = y[0] + 10 - ShiftTO;
      if (fabs(_E0) < delta3)
	iTO = 1.4 * (0.2 + y[2] / (KmTO + y[2]));
      else
	iTO = 0.28 * (0.2 + y[2] / (KmTO + y[2])) * _E0 / (1 - exp(-0.2 * _E0));
      iTO = gTO * iTO * (y[10] * exp(0.02 * y[0]) - y[2] * exp(-0.02 * y[0]));
      if (TOmode == 1)
	iTO = y[17] * iTO;

      break;


    case 2:
      /* When TOmode = 2 use full equations */

      iTO = y[17] * y[18] * 0.28 * (0.2 + y[2] / (KmTO + y[2]));
      iTO = gTO * iTO * (y[10] * exp(0.02 * y[0]) - y[2] * exp(-0.02 * y[0]));

      break;


    case 3:
    case 4:
      /* When TOmode = 3,4 use Hilgemann formulation */

      _E0 = y[0] - ShiftTO;
      if (TOmode == 4)
	f[16] = (1 / (1 + exp((y[0] + 35) / -5)) - y[16]) * alpha[16];
      else {
	f[16] = 0.0;
	y[16] = 1 / (1 + exp((y[0] + 35) / -5));
      }
      iTO = gTO * y[16] * (y[0] - EK) / (1 + exp(_E0 / 25));

      break;

    case 5:   /*Hilgemann-Noble 1987 formulation */
      _E0 = y[0] + 4 - ShiftTO;
      f[18] = 333 * (1 / (1 + exp(_E0 / -5)) - y[18]);
      iTO = gTO * y[17] * y[18] * (y[0] - EK);

      break;
      /*2.2*/

    }/*of TOmode cases*/


    /* IK uses rate theory inward rectifier, no -ve slope, no xover */
    /* except when Kmode = 5, when linear relation assumed */

    if (Kmode == 5)
      iK = gK * (6.5 * y[5] * y[5] * (y[0] - EK) + 1.8 * y[35] * y[35] * (y[0] - EK));
    else {
      iK = y[5] * iKm * (y[10] - y[2] * exp(-(y[0] / RTonF))) / 140;
      if (Kmode == 3)
	iK = y[5] * iK;
    }

    /* Version 1.3 uses Acetyl-choline activated channel with i(V)
            relation based on IK fully-activated relation */

    /* Version 3.7 also uses Mark Boyett's formulation when required*/

    if (AChgKmax != 0) {
      if (ACh != 0) {
	AChgK = AChgKmax * PoN(ACh, 1.4969) /
		(PoN(AChRecKm, 1.4969) + PoN(ACh, 1.4969));
	iACh = AChgK * (y[2] / (10 + y[2])) * (y[0] - EK) /
	       (1 + exp((y[0] - EK - 140) / (2.5 * RTonF)));
	iACh *= y[50] * y[51];
      } else
	iACh = 0.0;
    } else
      iACh = iAChm * (y[10] - y[2] * exp(-(y[0] / RTonF))) / 140;

    /* i f based on DiFrancesco (1981,1982) - K activation with
           linear K and Na components */

    ifK = y[4] * gfK * (y[2] / (Kmf + y[2])) * (y[0] - EK);
    ifNa = y[4] * gfNa * (y[2] / (Kmf + y[2])) * (y[0] - ENa);

    /* Ymode 10 and 20 USE Y SQUARED KINETICS */

    if (Ymode > 9) {
      ifK = y[4] * ifK;
      ifNa = y[4] * ifNa;
    }

    /* K component of ICA assumes 1% K permeability */

    _Z1 = PCaK * PCa * (y[0] - VSurfCa) /
	  (RTonF * (1 - exp((VSurfCa - y[0]) / RTonF)));
    iCaK = y[6] * CaCHon * _Z1 *
	(y[10] * exp(VSurfCa / RTonF) - y[2] * exp((VSurfCa - y[0]) / RTonF));
    if (Camode == 2 || Camode == 4 || Camode == 5 || Camode == 6 ||
	Camode == 9 || Camode == 10)
      iCaK = y[14] * iCaK;
    if (Camode == 10)
      iCaK = y[52] * y[53] * y[54] * y[55] * iCaK;
    if (SPvol != 0) {
      ICaKS = iCafract * iCaK;
      iCaK = (1 - iCafract) * iCaK;
    } else
      ICaKS = 0.0;

    /* NaK pump current activated by Nai and Kc */

    iP = Pump * (y[1] / (KmNa + y[1])) * y[2] / (Km + y[2]);

    /* K contribution to background channel */

    ibK = gbK * (y[0] - EK);

    /* ATP dependent K channel --- Nichols & Lederer, 1990 */

    ikATP = (y[0] + 80) * gkATPm / (1 + y[34] / kATP * (y[34] / kATP));

    /* compute total K flux as a current */

    imK = iK + iK1 + ifK + iCaK + ICaKS + ibK + iTO + ikATP + iP / (1 - nNaK);

  }


  /* all K dependent currents computed as function of [K]c and of
             radial distance for cylinder or sphere. current functions
             otherwise same as in previous section */

  if (Space < 3) {
    for (_L = 0; _L <= depth; _L++) {
      if (KCE[_L] < delta3)
	KCE[_L] = delta3;
      EKC[_L] = RTonF * log(KCE[_L] / y[10]);
      _E0 = y[0] - EKC[_L] + 10 - ShiftK1;
      IK1C[_L] = gK1 * (KCE[_L] / (KmK1 + KCE[_L])) * (y[0] - EKC[_L]);
      IK1C[_L] /= 1 + exp(_E0 / (0.5 * RTonF));
      _E0 = y[0] + 10 - ShiftTO;
      if (fabs(_E0) < 0.001)
	ITOC[_L] = 5 * (0.2 + KCE[_L] / (KmTO + KCE[_L]));
      else
	ITOC[_L] = (0.2 + KCE[_L] / (KmTO + KCE[_L])) * _E0 / (1 - exp(-0.2 * _E0));
      ITOC[_L] = 0.28 * gTO * ITOC[_L] *
		 (y[10] * exp(0.02 * y[0]) - KCE[_L] * exp(-0.02 * y[0]));
      if (TOmode != 0)
	ITOC[_L] = y[17] * ITOC[_L];
      IKC[_L] = y[5] * iKm * (y[10] - KCE[_L] * exp(-(y[0] / RTonF))) / 140;
      /*IACHC[_L] = iAChm * (y[10] - KCE[_L] * exp(-(y[0] / RTonF))) / 140;*/
	  /*Version 1.3*/
      if (Kmode == 3)
	IKC[_L] = y[5] * IKC[_L];
      IFKC[_L] = y[4] * gfK * (KCE[_L] / (Kmf + KCE[_L])) * (y[0] - EKC[_L]);
      IFNAC[_L] = y[4] * gfNa * (KCE[_L] / (Kmf + KCE[_L])) * (y[0] - ENa);
      IPC[_L] = Pump * (y[1] / (KmNa + y[1])) * KCE[_L] / (Km + KCE[_L]);
      _Z1 = PCaK * PCa * (y[0] - VSurfCa) /
	    (RTonF * (1 - exp((VSurfCa - y[0]) / RTonF)));
      ICAKC[_L] = _Z1 * (y[10] * exp(VSurfCa / RTonF) -
			 KCE[_L] * exp((VSurfCa - y[0]) / RTonF));
      ICAKC[_L] = y[6] * CaCHon * ICAKC[_L];
      if (Camode == 2 || Camode == 4 || Camode == 5 || Camode == 6 ||
	  Camode == 9)
	ICAKC[_L] = y[14] * ICAKC[_L];
      /* need to add CaMode 10 */
      IBKC[_L] = gbK * (y[0] - EKC[_L]);
      IMKC[_L] = IKC[_L] + IK1C[_L] + IFKC[_L] + ICAKC[_L] + IBKC[_L] +
		 IPC[_L] / (1 - nNaK);

    }

    K = 0.0;
    iK1 = 0.0;
    iK = 0.0;
    ifK = 0.0;
    ibK = 0.0;
    ifNa = 0.0;
    iP = 0.0;
    imK = 0.0;
    iCaK = 0.0;
    iTO = 0.0;
    iACh = 0.0;
  }

  /* total K dependent currents obtained by integrating radially */

  if (Space == 1) {  /* cylindrical case */
    for (_L = 1; _L <= depth; _L++) {
      _Z1 = (_L - 0.5) / (depth * depth);
      iK1 += IntRad(IK1C);
      iTO += IntRad(ITOC);
      iK += IntRad(IKC);
      /*Version 1.3 --- Set IACH if needed*/
      ifK += IntRad(IFKC);
      ifNa += IntRad(IFNAC);
      iP += IntRad(IPC);
      iCaK += IntRad(ICAKC);
      ibK += IntRad(IBKC);
      imK += IntRad(IMKC);
      K += IntRad(KC);
    }
    y[2]=K;
  }
  if (Space == 2) {  /* spherical case */
    for (_L = 1; _L <= depth; _L++) {
      _Z1 = (_L * _L * _L - (_L - 1.0) * (_L - 1) * (_L - 1)) /
	    (depth * depth * depth);
      iK1 += 0.5 * IntRad(IK1C);
      iTO += 0.5 * IntRad(ITOC);
      iK += 0.5 * IntRad(IKC);
      /*Version 1.3  ------- Set IACH if needed*/
      ifK += 0.5 * IntRad(IFKC);
      ifNa += 0.5 * IntRad(IFNAC);
      iP += 0.5 * IntRad(IPC);
      iCaK += 0.5 * IntRad(ICAKC);
      ibK += 0.5 * IntRad(IBKC);
      imK += 0.5 * IntRad(IMKC);
      K += 0.5 * IntRad(KC);
    }
    y[2]=K;
  }



  /* CALCIUM CURRENT */

  /* for details of formulation of components of ICA see
             Reuter & Scholz, J.Physiol,1977,264,17-47
             The formulation here uses the constant field
             current voltage relations from that paper but does not
             use their (incorrect) equation for the channel
             reversal potential. See Campbell, Giles, Hume & Noble,
             J.Physiol (1987). The reversal potential is computed in
             procedure CAREV.

             Up to three types of calcium channel are available.

             The correspondence with patch clamp work is that ICA
             corresponds best to the l current, ICA2 corresponds to
             the T current and ICA3 corresponds to the ICas of
             Lee et al, Proc Roy Soc 1984 */

  _Z1 = (y[0] - VSurfCa) / (RTonF * (1 - exp((VSurfCa - y[0]) / (0.5 * RTonF))));
  _Z2 = _Z1 * (y[3] * exp(VSurfCa / (0.5 * RTonF)) -
	       Cao * exp((VSurfCa - y[0]) / (0.5 * RTonF)));
  iCaCa = 4 * PCa * y[6] * CaCHon * _Z2;

  if (SPvol != 0) {
    _Z3 = _Z1 * (y[36] * exp(VSurfCa / (0.5 * RTonF)) -
		 Cao * exp((VSurfCa - y[0]) / (0.5 * RTonF)));
    ICaCaS = 4 * PCa * y[6] * CaCHon * _Z3;
  }

  if (PCa2 != 0)
    iCa2 = 4 * PCa2 * y[31] * y[32] * _Z2;
  else
    iCa2 = 0.0;
  if (PCa3 != 0)
    iCa3 = 4 * PCa3 * y[15] * y[16] * _Z2;
  else
    iCa3 = 0.0;
  _Z1 = PCaK * PCa * (y[0] - VSurfCa) /
	(RTonF * (1 - exp((VSurfCa - y[0]) / RTonF)));
  iCaNa = y[6] * CaCHon * _Z1 *
	  (y[1] * exp(VSurfCa / RTonF) - Nao * exp((VSurfCa - y[0]) / RTonF));
  if (Camode == 2 || Camode == 4 || Camode == 5 || Camode == 6 ||
      Camode == 9 || Camode == 10) {
    if (SPvol != 0) {
      ICaNaS = y[42] * iCaNa;   /* !!!!! Need to use y[53] in Camode 10 */
      ICaCaS = y[42] * iCaCa;
    }
    if (Camode == 10) {
      iCaNa = y[52] * y[53] * y[54] * y[55] * iCaNa;
      iCaCa = y[52] * y[53] * y[54] * y[55] * iCaCa;
    }
    iCaNa = y[14] * iCaNa;
    iCaCa = y[14] * iCaCa;
  } else {
    ICaNaS = iCaNa;
    ICaCaS = iCaCa;
  }
  if (SPvol != 0) {
    ICaNaS = iCafract * iCaNa;
    ICaCaS = iCafract * iCaCa;
    iCaNa = (1 - iCafract) * iCaNa;
    iCaCa = (1 - iCafract) * iCaCa;
  }

  if (CaNmode == 0)
    iCaNa = 0.0;
  iCa = iCaCa + iCaNa + iCaK;
  if (SPvol != 0)
    iCas = ICaCaS + ICaNaS + ICaKS;
  else
    iCas = 0.0;



  /* CALCIUM-ACTIVATION OF IONIC CURRENTS */

  /* compute activation of Ca sensitive currents (other than inaca)
             each current is computed according to equation

                i  =  i ( mi  +  (cai/(cai + cact)))

             where mi is minimum value of current, expressed as a fraction
             of the value in amode 0 and cact is the calcium concentration
             for half-activation */

  if (Amode == 1 || Amode == 9)
    iK1 *= MiK1 + y[3] / (y[3] + CACT[0]);
  if (Amode == 2 || Amode == 9)
    iK *= MiK + y[3] / (y[3] + CACT[1]);
  if (Amode == 3 || Amode == 9) {
    ibK *= MiB + y[3] / (y[3] + CACT[2]);
    ibNa *= MiB + y[3] / (y[3] + CACT[2]);
  }
  if (Amode == 4 || Amode == 9)
    iTO *= MiTO + y[3] / (y[3] + CACT[3]);
  if (Amode == 5 || Amode == 9)
    iNaK = 0.5 * (y[3] / (y[3] + CACT[2])) *
	   (gNaK * (y[0] - ENa) + gNaK * (y[0] - EK));
  else
    iNaK = 0.0;




  /* TOTAL CURRENT */

  itot = iNa + iCa + iCa2 + iCa3 + iK + ifNa + ifK + iP + ibNa + ibK + ibCa +
	 iK1 + iTO + iPulse + Stim + iNaCa + iNaK + ikATP;
  itot += iACh;
  if (SPvol != 0)
    itot += iCas + INaCaS;
  iSI = iCa + iCa2 + iCa3 + iNaCa;
  if (SPvol != 0)
    iSI += iCas + INaCaS;

  /* set dV/dt = 0 to give voltage clamp (modes 2 - 8).  in  mode  2
     compute effect of series resistance and of external clamp circuit.
             in action potential mode set
             dV/dt = - total current/capacitance */

  f[0] = -CAP * itot;


  #undef IntRad

}  /* of procedure CURRENTS */
/* End. */

/* #endif  QHEART */





