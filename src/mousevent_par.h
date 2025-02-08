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

/*
Cell morphology: Table 1
*/

_(     acap , 1.534e-4  , 1 ) // cm2 capacitive mambrane area
_(     vmyo  , 25.84e-6  , 1 ) // microL
_(     vjsr  , 0.12e-6   , 1 ) //
_(     vnsr  , 2.098e-6  , 1 ) //
_(     vss  , 1.485e-9  , 1 ) //

/*
SR Parameters table 4
*/

_(       kplusa  , 0.006075 	 , 1 ) // microM^-4/ms 
_(       kminusa  , 0.07125 	 , 1 ) // /ms 
_(       kplusb  , 0.00405  	 , 1 ) // microM^-3/ms 
_(       kminusb  , 0.965 	 , 1 ) // /ms 
_(       kplusc  , 0.009 		 , 1 ) // /ms 
_(       kminusc  , 0.0008 	 , 1 ) // /ms 
_(       n  , 4 			 , 1 ) // cooperativity parameters
_(       m  , 3, 1 )
_(     nu1  , 4.5        , 1 ) //  /ms
_(     nu2  , 1.74e-5    , 1 ) //  /ms
_(     nu3  , 0.45       , 1 ) //  microM/ms
_(     kmup  , 0.5       , 1 ) // microM
_(     tautr  , 20.0     , 1 ) // ms 
_(     tauxfer  , 8.0    , 1 ) // ms 

/*
Buffering parameters Table 6
*/

_(      ltrpntot  , 70.0          , 1 ) // microM *
_(      htrpntot  , 140.0         , 1 ) // microM *
_(      kplushtrpn  , 0.00237     , 1 ) // microM^-1/ms
_(      kminushtrpn  , 3.2e-5     , 1 ) // /ms
_(      kplusltrpn  , 0.0327      , 1 ) // microM^-1/ms
_(      kminusltrpn  , 0.0196     , 1 ) // /ms
_(      cmdntot  , 50.0           , 1 ) // microM *
_(      csqntot  , 15000.0        , 1 ) // microM *
_(      kmcmdn  , 0.238           , 1 ) // microM *
_(      kmcsqn  , 800.0           , 1 ) // microM *
 
/*
Parameters: table 7.
*/

_(     cm  , 1.0 		 , 1 ) // microF/cm2
_(     F  , 96.5 	 , 1 ) // C/mmol
_(     R  , 8.314	  , 1 ) // J/(mol-K)
_(     temp  , 298.0 , 1 )
_(     knaca  , 292.8 	 , 1 ) // pA/pF: scaling factor of NaCa exchanger
_(     kmna  , 87500.0 	 , 1 ) // micro-M: sodium half saturation for NaCa exchanger
_(     kmca  , 1380.0 	 , 1 ) // micro-M: calcium half saturation for NaCa exchanger
_(     ksat  , 0.1 	 , 1 ) // NaCa saturation at very negative potentials
_(     eta  , 0.35 	 , 1 ) // controls voltage dependence of NaCa exchanger
_(     imaxnak  , 0.88 	 , 1 ) // pA/pF  , 1 ) // max. INaK exchange current
_(     kmnai  , 21000.0 	 , 1 ) // micro-M  , 1 ) // K+ half saturation constant for Na/K exchange current
_(     kmko  , 1500.0 	 , 1 ) // micro-M  , 1 ) // K+ half saturation constant for Na/K exchange current
_(     ipcamax  , 1.0 	 , 1 ) // pA/pF Maximum Ca2+ pump current
_(     kmpca  , 0.5 	 , 1 ) // micro-M Ca2+ half-saturation constant for Ca2+ pump current
_(     gcab  , 0.000367 	 , 1 ) // milli-S/micro-F: max. background calcium conductance
_(     gna  , 13.0 	 , 1 ) // milli-S/micro-F: fast Na current conductance.
_(     gnab  , 0.0026 	 , 1 ) // milli-S/micro-F: background Na current conductance
_(     gktof  , 0.4067 	 , 1 ) // milli-S/micro-F: maximum transient outward current conductance (apex) FAST
// _(     gktof  , 0.0798  , 1 ) // milli-S/micro-F: maximum transient outward current conductance (septum) FAST
_(     gks  , 0.00575 	 , 1 ) // milli-S/micro-F: slow delated rectifier current conductance
_(     gktos  , 0.0 	 , 1 ) // milli-S/micro-F: maximum transient outward current conductance (apex) SLOW
// _(     gktos  , 0.0629  , 1 ) // milli-S/micro-F: maximum transient outward current conductance (septum) SLOW
_(     gkur  , 0.16 	 , 1 ) // milli-S/micro-F: ultrarapidly delated rectifier current conductance (apex)
// _(     gkur  , 0.0975 	 , 1 ) // milli-S/micro-F: ultrarapidly delated rectifier current conductance (septum)
_(     gkss  , 0.05 	 , 1 ) // milli-S/micro-F: non-inactivating potassium current conductance (apex)
// _(     gkss  , 0.0324 	 , 1 ) // milli-S/micro-F: non-inactivating potassium current conductance (septum)
_(     gkr  , 0.078 	 , 1 ) // milli-S/micro-F: rapid delayed-rectifier potassium current conductance
_(     kf  , 0.023761 	 , 1 ) // ms^-1: rate constant for rapid delayed-rectifier potassium current Ikr 
_(     kb  , 0.036778 	 , 1 ) // ms^-1: rate constant for rapid delayed-rectifier potassium current Ikr
_(     gclca  , 10.0 	 , 1 ) // milli-S/micro-F: max. ICl(Ca) current conductance
_(     kmcl  , 10.0 	 , 1 ) // micro-M Half-saturation constant for Ca activated Cl- current
_(     ecl  , -40.0 	 , 1 ) // mV: reversal potential  for Ca2+ activated Cl- current


/* L-type calcium channel parameters table 5 */

_(     gcal  , 0.1729 	 , 1 ) // milliS/microF: maximum conductance of the L-type calcium channel.
_(     kpcmax  , 0.23324 	 , 1 ) // ms^-1 Maximum time constant for ca2+ induced inactivation
_(     kpchalf  , 20.0 	 , 1 ) // microM: half saturation constant for ca2+induced inactivation
_(     ecal  , 63.0 	 , 1 ) // mV
_(     kpcb  , 0.0005 	 , 1 ) // ms^-1 voltage-insensitive rate constant for ICaL inactivation
_(     icalmax  , 7.0  , 1 ) // pA/pF

/*
Extracellular concentrations table 3.
*/

_(     ko  , 5400.0   , 1 ) // microM
_(     nao  , 140000.0  , 1 ) // microM
_(     cao  , 1800.0   , 1 ) // microM

// soicr
_(	soicr_switch , 0.0 , 1 ) // soicr usually is switched off.
_(	soicr_max , 0.45 , 1 ) // 10% of nu1
_(	soicr_thresh , 8400.0 , 1 ) // 8.4 mM, i.e. 8400 microM. thresh is comparable to carel. In basal model, carel is ~ 0.1 mM
_(	soicr_slope , 0.0001 ,  1 ) // 0.001 mM
_(	soicr_tau , 125.0 , 1 )

// for ISO
// most of them are to do with ss of ical and ikr - in the ventricle model, that is markovian
