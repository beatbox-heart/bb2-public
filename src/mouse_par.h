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

//---------------------------------------------------------------------------
// Constants. If really required, you can type in the units of the constants
// in place of the , 1 ) as indicated in the comments. formulas are treated otherwise.
//---------------------------------------------------------------------------

_(         R , 8.314472 , 1 ) // mol^-1
_(         T  , 310.5 , 1 )    
_(         F  , 96.4846 , 1 )   // units: C/mol
_(         capacitance  , 0.025 , 1 )  // nF
_(         vcell  , 3.0 , 1 )  // pL
_(         l_cell  , 66.3767257 , 1 )  // microM
_(         r_cell  , 3.792956 , 1 )  // microM
_(         vrel  , 0.0036 , 1 ) 
_(         vsub  , 0.03328117 , 1 ) 
_(         vup   , 0.0348 , 1 ) 
_(         vi  , 1.34671883 , 1 ) 
_(         Mgi  , 2.5 , 1 ) 	 	// Internal Mg2+ concentration
_(         nao  , 140.0 , 1 )        // mM
_(         cao  , 1.8 , 1 )          // mM
_(         ko  , 5.4 , 1 )           // mM
_(         gst  , 0.00006 , 1 )  // it can be lower. Free parameter. Maltsev has 0.000075 microS
_(         eist  , 17.0 , 1 )  // mV reversal of Ist
_(         gbna  , 0.0001215 , 1 )          //uS 
_(         gbca  , 0.000015 , 1 )          //uS 
_(         gbk   , 0.0000025 , 1 )          //uS
_(         gk1  , 0.000808489 , 1 )              // microS: mangoni is 0.0009 microS
_(         gkr  , 0.002364 , 1 )       // micro-S   increased it temporarily to improve AP 0.8*0.002955
_(         gks  , 0.000299 , 1 )     // micro-S ki=140 mM
_(         gcal12  , 0.006 , 1 ) 
_(         gcal13  , 0.018 , 1 ) 
_(         ecal  , 47.0 , 1 ) 
_(         kmfca  , 0.00035 , 1 ) 
_(         alpha_fca  , 0.021 , 1 ) 
_(         ecat  , 45.0 , 1 ) 
_(         enattxr  , 41.5761 , 1 ) 
_(         gsus  , 0.00039060 , 1 )         // micro-S
_(         inakmax  , 0.14245 , 1 )  // 2.88 pA/pF as in Kurata and Maltsev x conversion to give nA
_(         kmnap  , 14.0 , 1 )  // mM AFFECTED BY ISO
_(         kmkp  , 1.4 , 1 )  // mM
_(         kNaCa  , 5.5 , 1 ) 
_(         K1ni  , 395.3 , 1 ) 
_(         K1no  , 1628.0 , 1 ) 
_(         K2ni  , 2.289 , 1 ) 
_(         K2no  , 561.4 , 1 ) 
_(         K3ni  , 26.44 , 1 ) 
_(         K3no  , 4.663 , 1 ) 
_(         Kci  , 0.0207 , 1 ) 
_(         Kco  , 3.663 , 1 ) 
_(         Kcni  , 26.44 , 1 ) 
_(         Qci  , 0.1369 , 1 ) 
_(         Qco  , 0.0 , 1 ) 
_(         Qn  , 0.4315 , 1 ) 
_(         tdifca  , 0.04 , 1 )  // as in Maltsev. This is diffusion from Cai to Casub, this is not to do with the transients.
_(         Krel  , 0.0015 , 1 ) 
_(         nrel  , 2.0 , 1 ) 	//% SR Ca2+ release Km (mM) and Hill coefficient
_(         Kup  , 0.0006 , 1 )  	// 0.0006
_(         nup  , 1.0 , 1 ) 		//% SR Ca2+ uptake Km (mM) and Hill coefficient
_(         Ttr  , 40.0 , 1 )  	// Time constant for Ca2+ transfer from NSR to JSR (ms) 60
_(         ConcTC  , 0.031 , 1 )             //              % Concentration of Troponin-Ca complex (mM)
_(         ConcTMC  , 0.062 , 1 )            //      % Concentration of Troponin-Mg complex (mM)
_(         kfTC  , 88.8 , 1 ) 
_(         kfTMC  , 237.7 , 1 )              //% Rate constant for Ca2+ binding to Troponin(mM/ms)
_(         kbTC  , 0.446 , 1 ) 
_(         kbTMC  , 0.00751 , 1 )    //      % Rate constant for Ca2+ unbinding from Troponin (ms-1)
_(         kfTMM  , 2.277 , 1 )              //              % Rate constant for Mg2+ binding to Troponin(mM/ms)
_(         kbTMM  , 0.751 , 1 )              //              % Rate constant for Mg2+ unbinding from Troponin (ms-1)
_(         ConcCM  , 0.045 , 1 )     //                      % Concentration of Calmodulin (mM)
_(         kfCM  , 237.7 , 1 )               // manip.               % Rate constant for Ca2+ binding to Calmodulin (mM/ms)
_(         kbCM  , 0.542 , 1 )               //              % Rate constant for Ca2+ unbinding from Calmodulin (mM/ms)
_(         ConcCQ  , 10.0 , 1 )              // % Concentration of Calsequestrin (mM)
_(         kfCQ  , 0.534 , 1 )               //              % Rate constant for Ca2+ binding to Calsequestrin (mM/ms)
_(         kbCQ  , 0.445 , 1 )               //              % Rate constant for Ca2+ unbinding from Calsequestrin (mM/ms)
_(         koca  , 10.0 , 1 )   // mM-2 ms-1
_(         kom  , 0.06 , 1 )  // ms-1
_(         kica  , 0.5 , 1 )  // mM-1ms-1
_(         kim  , 0.005 , 1 )  // ms-1
_(         eca50sr  , 0.45 , 1 )  // mM
_(         ks  , 1300000.0 , 1 )  // ms-1 for a frequency of 5 hertz or more
_(         maxsr  , 15.0 , 1 ) 
_(         minsr  , 1.0 , 1 ) 
_(         hsrr  , 2.5 , 1 ) 
_(         pumpkmf  , 0.000246 , 1 )  // mM
_(         pumpkmr  , 3.29 , 1 )  // mM
_(         pumphill  , 2.0 , 1 )  // 2 seems to be alright
_(         gna_ttxs  ,  5.925e-06 , 1 ) 
_(         gna_ttxr  ,  5.925e-06 , 1 ) 
_(         gcat  , 0.013965 , 1 ) 
_(         gh  , 0.0057 , 1 ) 
_(         gto  , 0.00492 , 1 ) 
_(         Pup  , 0.04 , 1 )
_(	   soicr_switch , 0.0 , 1 ) // soicr usually is switched off.
_(	   soicr_max , 25.0 , 1 )
_(	   soicr_thresh , 8.4 , 1 ) // 8.4 mM, i.e. 8400 microM. thresh is comparable to carel. In basal model, carel is ~ 0.1 mM
_(	   soicr_slope , 0.001 ,  1 ) // 0.001 mM
_(	   soicr_tau , 125.0 , 1 )
_(	   if_vhalf , 106.8 , 1 ) // ISO
_(	   dl13vhalf , 13.5 , 1 )
_(	   fl13vhalf , 35.0 , 1 )
_(	   dl12vhalf , 3.0 , 1 )
_(	   fl12vhalf , 36.0 , 1 )
_(	   ikractvhalf , 21.173694 , 1 )
