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

/* Describes the parameters for the LRd model.
 * Each declaration takes the form:
 * _(variable name, default value, unit)
 * If unit is not known,  1 is used
 *
 * It is permissible to use simple expressions for the initial value,
 * including previously declared variables.   */

/*>>>>>>>>>>>>>>>>>> Constants from LRd1994.c <<<<<<<<<<<<<<<<<<*/


_(R, 8314.0, J kmol-1 K-1) /*  The universal gas constant */
_(F, 96500.0, C mol-1) /*  Faraday's constant */
_(T, 310.0, K) /*  Temperature */

/*------------------------------------------------------------------------------
Variables for cell geometry calculations (I)
------------------------------------------------------------------------------*/

_(L, 0.01, cm) /*  The length of a cell in the model */
_(r, 0.0011, cm) /*  The radius of a cell in the model */

/*------------------------------------------------------------------------------
Variables for the fast sodium current calculations (IIIa)
------------------------------------------------------------------------------*/

_(gna, 16.0, mS/uf) /*  Maximum conductance of the Na channel */

/*------------------------------------------------------------------------------
Variables for calculating the current through the L-type calcium channel (IIIb)
------------------------------------------------------------------------------*/

_(kmca, 6.0e-4, mmol/L) /*  Half saturation conc. of L-type Ca channel */
_(pca, 5.4e-4, cm/s) /*  Permeability of the membrane to calcium */
_(pna, 6.75e-7, cm/s) /*  Permeability of the membrane to sodium */
_(pk, 1.93e-7, cm/s) /*  Permeability of the membrane to K */
_(gamma_cai, 1.0, 1) /*  Activity coefficient of calcium */
_(gamma_cao, 0.341, 1) /*  Activity coefficient of calcium */
_(gamma_nai, 0.75, 1) /*  Activity coefficient of sodium */
_(gamma_nao, 0.75, 1) /*  Activity coefficient of sodium */
_(gamma_ki, 0.75, 1) /*  Activity coefficient of potassium */
_(gamma_ko, 0.75, 1) /*  Activity coefficient of potassium */

/*------------------------------------------------------------------------------
Variables for the time dependent potassium current calculations (IIIc)
------------------------------------------------------------------------------*/

_(pnak, 0.01833, 1) /*  Sodium/Potassium permeability ratio */

/*------------------------------------------------------------------------------
Variables for the plateau potassium current caclulations (IIIe)
------------------------------------------------------------------------------*/

_(gkp, 0.0183, ms/uF) /*  Max. channel conductance of ikp */

/*------------------------------------------------------------------------------
Variables for the sodium-calcium exchanger current calculations (IIIf)
------------------------------------------------------------------------------*/

_(knaca, 2000, 1) /*  Scaling factor of inaca */
_(kmna, 87.5, mmol/L) /*  Half saturation channel conc. for Na */
_(kmca2, 1.38, mmol/L) /*  Half saturation channel conc. for Ca  */
_(ksat, 0.1, 1) /*  Saturation factor of inaca at very -ve potentials */
_(eta, 0.35, 1) /*  Position of the energy barrier controlling voltage dependence of inaca */

/*------------------------------------------------------------------------------
Variables for the sodium-potassium pump current calculations (IIIg)
------------------------------------------------------------------------------*/

_(inak_max, 1.5, uA/uF) /*  Maximum current through the Na-K pump */
_(kmnai, 10.0, mmol/L) /*  Half saturation concentration for Na */
_(kmko, 1.5, mmol/L)/*  Half saturation concentration for K */

/*------------------------------------------------------------------------------
Variables for the nonspecific calcium activated current calculations (IIIh)
------------------------------------------------------------------------------*/

_(pnsca, 1.75e-7, cm/s) /*  Permiability of the nonspecfic Ca channel to Na and K */
_(kmnsca, 0.0012, mmol/L) /*  Half saturation conc. of the nonspecfic Ca channel */

/*------------------------------------------------------------------------------
Variables for the sarcolemmal calcium pump calculations (IIIi)
------------------------------------------------------------------------------*/

_(ipca_max, 1.15, uA/uF) /*  Max. Ca current through the sarcolemmal Ca pump */
_(kmpca, 0.0005, mmol/L) /*  Half saturation conc. of the sarcolemmal Ca pump */

/*------------------------------------------------------------------------------
Variables for the calcium background current calculations (IIIj)
------------------------------------------------------------------------------*/

_(gcab, 0.003016, mS/uF) /*  Maximum conductance of Ca background */

/*------------------------------------------------------------------------------
Variables for the sodium background current calculations (IIIk)
------------------------------------------------------------------------------*/

_(gnab, 0.00141, mS/uF) /*  Maximum conductance of Na background */

/*------------------------------------------------------------------------------
Variables for the calcium buffers in the myoplasm calculations (IV)
------------------------------------------------------------------------------*/

_(trpn_max, 0.07, mmol/L) /*  Maximum Ca buffered in troponin */
_(cmdn_max, 0.05, mmol/L) /*  Maximum Ca buffered in calmodulin */
_(kmtrpn, 0.0005, mmol/L) /*  Equilibrium constant of buffering for TRPN */
_(kmcmdn, 0.00238, mmol/L) /*  Equilibrium constant of buffering for CMDN */

/*------------------------------------------------------------------------------
Variables for the calcium induced calcium release of JSR caclulations (Va)
------------------------------------------------------------------------------*/

_(kmrel, 0.0008, mmol/L) /*  Half saturation conc. of Ca induced Ca release of JSR */
_(tau_on, 2.0, ms) /*  Time constant of activating Ca release from JSR */
_(tau_off, 2.0, ms) /*  Time constant of deactivating Ca release from JSR */
_(delta_caith, 0.00018, mmol/L) /*  Threshold for external triggering of Ca release from JSR */

/*------------------------------------------------------------------------------
Variables for the Ca release of JSR under Ca overload conditions calculations (Vb)
------------------------------------------------------------------------------*/

_(csqnth, 7.0, 1) /*  Threshold of Ca bound CSQN for internal triggering of Ca release from JSR under Ca overload conditions */

/*------------------------------------------------------------------------------
Variables for the Ca buffer in JSR and CSQN calcualtion (Vc)
------------------------------------------------------------------------------*/

_(csqn_max, 10.0, mmol/L) /*  Max. CSQN buffered Ca concentration */
_(kmcsqn, 0.8, mmol/L) /*  Equilibrium constant of buffering for CSQN */


/*------------------------------------------------------------------------------
Variables for the calcium uptake and leakage of NSR calculations (Vd)
------------------------------------------------------------------------------*/

_(kmup, 0.00092, mmol/L) /*  Half saturation conc. of Ca uptake of NSR */
_(iup_max, 0.005, mmol L-1 ms-1) /*  Max. current through Ca uptake of NSR */
_(cansr_max, 15.0, mmol/L) /*  Maximum calcium in NSR */
_(caiinit, 1.2e-4, mmol/L) /*  Initial value for cai. */

/*------------------------------------------------------------------------------
Variables for the translocation of calcium ions from NSR to JSR (Ve)
------------------------------------------------------------------------------*/

_(tau_tr, 180.0, ms) /*  Calcium transfer time constant */


/*------------------------------------------------------------------------------
Ionic valences
------------------------------------------------------------------------------*/

_(zna, 1.0, 1) /*  Sodium ion valence */
_(zk, 1.0, 1) /*  Potassium ion valence */
_(zca, 2.0, 1) /*  Calcium ion valence */


/*>>>>>>>>>>>>>>>>>>>>>>> Variables from LRd1994.c <<<<<<<<<<<<<<<<<<<<<<<<<<<*/
/*>>>>>>>>>>>>>>>>>>>>>>> LHS once only - in main() <<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*------------------------------------------------------------------------------
Variables for cell geometry calculations
------------------------------------------------------------------------------*/
_(vcell, M_PI * r * r * L * 1000.0, uL)
_(ageo, 2.0 * M_PI * r * r + 2.0 * M_PI * r * L, cm2)
_(acap, 2.0 * ageo, cm2)
_(vmyo, vcell * 0.68, vmyo)
_(vmito, vcell * 0.26, uL)
_(vsr, vcell * 0.06, uL)
_(vnsr, vcell * 0.0552, uL)
_(vjsr, vcell * 0.0048, uL)
_(vcleft, (vcell / 0.88) * 0.12, uL) /*  Unused,  but included for completeness. */

/*------------------------------------------------------------------------------
Ion concentrations that remain unchanged.
------------------------------------------------------------------------------*/
_(nao, 140.0, mmol/L)
_(ko, 5.4, mmol/L)
_(cao, 1.8, mmol/L)
