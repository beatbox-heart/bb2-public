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

#include "error.h"

/*------------------------------------------------------------------------------
Variables relating to the membrane potential
------------------------------------------------------------------------------*/
double V_old = 0.0; /* The old value of the membrane potential (mV) UTILITY */


/*------------------------------------------------------------------------------
Variables relating to ionic concentrations
------------------------------------------------------------------------------*/
double dcajsr = 0.0; /* The change in JSR concentration due to influx / efflux (mmol/L) */

/*------------------------------------------------------------------------------
Variables for the fast sodium current calculations (IIIa)
------------------------------------------------------------------------------*/ 
#ifndef CELLTEST
double ina = 0.0; /* The fast sodium current (uA/uF)  */
#endif
double ena = 0.0; /* The reversal potential of sodium (mV) */


/*------------------------------------------------------------------------------
Variables for calculating the current through the L-type calcium channel (IIIb)
------------------------------------------------------------------------------*/

#ifndef CELLTEST
	double fca = 0.0; /* Ca dependent inactivation gate of the L-type Ca channel */
	double ica = 0.0; /* Calcium current through the L-type calcium channel (uA/uF) */
	double icak = 0.0; /* Potassium current through the L-type calcium channel (uA/uF) */
	double icana = 0.0; /* Sodium current through the L-type calcium channel (uA/uF) */
#endif

double ica_max = 0.0; /* Max. Ca current through the L-type Ca channel (uA/uF) */
double icak_max = 0.0; /* Max. K current through the L-type Ca channel (uA/uF) */
double icana_max = 0.0; /* Max. Na current through the L-type Ca channel (uA/uF) */

#ifndef CELLTEST
	double icat = 0.0; /* Total current through the L-type calcium channel (uA/uF) */
#endif


/*------------------------------------------------------------------------------
Variables for the time dependent potassium current calculations (IIIc)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double ik = 0.0; /* The time dependent potassium current (uA/uF) */
#endif

double gk = 0.0; /* Maximum channel conductance of ik (mS/uF) */
//double xi = 0.0; /* Inactivation gate of the time dependent potassium current */ REMOVE
double ek = 0.0; /* Reversal potential of the time dependent K current (mV) */
//double ax = 0.0; /* Opening rate constant of gate x (ms^-1) UNUSED */
//double bx = 0.0; /* Closing rate constant of gate x (ms^-1) UNUSED */


/*------------------------------------------------------------------------------
Variables for the time independent potassium current calculations (IIId)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double ik1 = 0.0; /* The time independent potassium curent (uA/uF) */
#endif

double gk1 = 0.0; /* Maximum channel conductance of ik1 (mS/uF) */
double ek1 = 0.0; /* Reversal potential of the time independent K current (mV) */


/*------------------------------------------------------------------------------
Variables for the plateau potassium current caclulations (IIIe)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double ikp = 0.0; /* The plateau potassium current (uA/uF) */
#endif

double kp = 0.0; /* Potassium plateau factor */


/*------------------------------------------------------------------------------
Variables for the sodium-calcium exchanger current calculations (IIIf)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double inaca = 0.0; /* The sodium-calcium exchanger current (uA/mF) */
#endif


/*------------------------------------------------------------------------------
Variables for the sodium-potassium pump current calculations (IIIg)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double inak = 0.0; /* The sodium-potassium pump current (uA/uF) */
#endif

double fnak = 0.0; /* Voltage dependence parameter of the sodium-potassium pump */
double sigma = 0.0; /* [Na]o dependence factor of fnak */


/*------------------------------------------------------------------------------
Variables for the nonspecific calcium activated current calculations (IIIh)
------------------------------------------------------------------------------*/
double insk = 0.0; /* Non-specific potassium current (uA/uF) */
double insk_max = 0.0; /* Max. K current through the nonspecific Ca channel (uA/uF) */
double insna = 0.0; /* Non-specific sodium current (uA/uF) */
double insna_max = 0.0; /* Max. Na current through the nonspecific Ca channel (uA/uF) */
double insca = 0.0; /* Nonspecific calcium activated current (uA/uF) */


/*------------------------------------------------------------------------------
Variables for the sarcolemmal calcium pump calculations (IIIi)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double ipca = 0.0; /* Sarcolemmal calcium pump current (uA/uF) */
#endif


/*------------------------------------------------------------------------------
Variables for the calcium background current calculations (IIIj)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double icab = 0.0; /* Calcium background current (uA/uF) */
#endif

double ecan = 0.0; /* Nernst potential for calcium (mV) */


/*------------------------------------------------------------------------------
Variables for the sodium background current calculations (IIIk)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double inab = 0.0; /* Sodium background current (uA/uF) */
#endif


/*------------------------------------------------------------------------------
Variables for the total time independent current calculations (IIIl)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double iv = 0.0; /* The total time indepedent current (uA/uF) */
#endif


/*------------------------------------------------------------------------------
Variables for the calcium buffers in the myoplasm calculations (IV)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double trpn = 0.0; /* Troponin buffered calcium */
	double cmdn = 0.0; /* Calmodulin buffered calcium */
#endif


/*------------------------------------------------------------------------------
Variables for the calcium induced calcium release of JSR caclulations (Va)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double irelcicr = 0.0; /* Ca release from JSR to myoplasm due to CICR (mmol/L per ms) */
#endif

double grelcicr = 0.0; /* Rate constant of Ca release from JSR due to CICR (ms^-1) */

/*------------------------------------------------------------------------------
Variables for the Ca release of JSR under Ca overload conditions calculations (Vb)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double ireljsrol = 0; /* Ca release from JSR to myoplasm due to JSR overload (mmol/L per ms) */
#endif

double greljsrol = 0; /* Rate constant of Ca release from JSR due to JSR overload (ms^-1) */
double greljsrol_max = 0; /* Rate constant of Ca release from JSR due to JSR overload (ms^-1) */


/*------------------------------------------------------------------------------
Variables for the calcium uptake and leakage of NSR calculations (Vd)
------------------------------------------------------------------------------*/
double kleak; /* Rate constant of calcium leakage from NSR (ms^-1) */

#ifndef CELLTEST
	double iup = 0; /* Calcium uptake myoplasm to NSR (mmol/L per ms) */
	double ileak = 0; /* Ca leakage from NSR to myoplasm (mmol/L per ms) */
#endif

/*------------------------------------------------------------------------------
Variables for the translocation of calcium ions from NSR to JSR (Ve)
------------------------------------------------------------------------------*/
#ifndef CELLTEST
	double itr = 0; /* Transfer current of Ca from the NSR to the JSR (mmol/L per ms) */
#endif


/*------------------------------------------------------------------------------
Gate utility variables
------------------------------------------------------------------------------*/
double SteadyState = 0; /* Steady-state value of activation gate */
double tau = 0; /* Time constant of gate (ms) */
double a = 0; /* The opening rate constant for the gate (ms^-1) */
double b = 0; /* The closing rate constant for the gate (ms^-1) */


/*------------------------------------------------------------------------------
Variables for determining the maximum value of dvdt
------------------------------------------------------------------------------*/
double caiontot = 0; /* The total Ca current (uA/uF) */
double dcaiontotnew = 0; /* The new change in the total Ca current (uA/uF) */

/*------------------------------------------------------------------------------
Variables for the Steffensen algorithm.
------------------------------------------------------------------------------*/
int k, Cond = 0, Max = 20; /* Max is the maximum number of iterations */
double P0, P1, P2, D1, D2, DP, P3, RelErr, value;
double Small = 1.0e-6, Delta = 1.0e-6; /* Set tolerances */
double CaTot_cai; // Temporary variable representing CaTot(cai).
double CaJSRTot_cajsr; // Temporary variable representing CaJSRTot(cajsr).
double gP0 = 0, gP1 = 0; // Temporary variables representing g(P0) and g(P1).

/******************************************************************************
*******************************************************************************
*******************************************************************************
************************		Calculations		***************
*******************************************************************************
*******************************************************************************
*******************************************************************************
******************************************************************************/

/*------------------------------------------------------------------------------
The fast sodium current calculations (IIIa)
------------------------------------------------------------------------------*/
    
	ena = ((R * T) / F) * log(nao / nai);

    /* Update ina m gate */
    if (fabs(V + 47.13) < 1e-10) /* Avoid denominator = 0 */
        a = 3.2;
    else
        a = 0.32 * (V + 47.13) / (1.0 - exp(-0.1 * (V + 47.13)));
    b = 0.08 * exp(-V / 11.0);
	
    if (t*dt == 0)
        m = a / (a + b); /* Set the initial (steady state) value of the gate */
    else
        m = (a / (a + b)) - ((a / (a + b)) - m) * exp(-dt / (1.0 / (a + b)));

    /* Update ina h gate */
    if (V < -40.0) {
        a = 0.135 * exp((80.0 + V) / -6.8);
        b = 3.56 * exp(0.079 * V) + 3.1e5 * exp(0.35 * V);    
    }
    else { 
	a = 0.0;
        b = 1.0 / (0.13 * (1.0 + exp((V + 10.66) / -11.1)));
    } 

    if (t*dt == 0)
        h = a / (a + b); /* Set the initial (steady state) value of the gate */
    else
        h = (a / (a + b)) - ((a / (a + b)) - h) * exp(-dt / (1.0 / (a + b)));

    /* Update ina j gate */
    if (V < -40.0) {
        a = (-1.2714e5 * exp(0.2444 * V) - 3.474e-5 * exp(-0.04391 * V)) *
          (V + 37.78) / (1.0 + exp(0.311 * (V + 79.23)));
        b = 0.1212 * exp(-0.01052 * V) / (1.0 + exp(-0.1378 * (V + 40.14)));
    }
    else {
        a = 0.0;
        b = 0.3 * exp(-2.535e-7 * V) / (1.0 + exp(-0.1 * (V + 32.0)));
    }

    if (t*dt == 0)
        j = a / (a + b); /* Set the initial (steady state) value of the gate */
    else
        j = (a / (a + b)) - ((a / (a + b)) - j) * exp(-dt / (1.0 / (a + b)));

    ina = gna * pow(m, 3) * h * j * (V - ena);


/*------------------------------------------------------------------------------
The calculations for the L-type calcium channel (IIIb)
------------------------------------------------------------------------------*/
    
	/* Update ical d gate */
    SteadyState = 1.0 / (1.0 + exp(-(V + 10.0) / 6.24));
        
        if (fabs(V + 10.0) < 1e-10) /* Avoid denominator = 0 */
        tau = SteadyState * (1.0 / 6.24) / 0.035;
    else
        tau = SteadyState * (1.0 - exp(-(V + 10.0) / 6.24)) /
          (0.035 * (V + 10.0));
		
    if (t*dt == 0)
        d = SteadyState; /* Set the initial (steady state) value of the gate */
    else
        d = SteadyState - (SteadyState - d) * exp(-dt / tau);
	
    /* Update ical f gate */
    SteadyState = (1.0 / (1.0 + exp((V + 32.0) / 8.0))) + (0.6 / (1.0 +
      exp((50.0 - V) / 20.0)));
          
    tau = 1.0 / (0.0197 * exp(-pow(0.0337 * (V + 10.0), 2)) + 0.02);

    if (t*dt == 0)
        f = SteadyState; /* Set the initial (steady state) value of the gate */
    else
        f = SteadyState - (SteadyState - f) * exp(-dt / tau);
	
    if (fabs(V) < 1e-10)  /* Avoid denominator = 0 */
        ica_max = pca * 2 * F * (gamma_cai * cai - gamma_cao * cao);
    else
        ica_max = pca * zca * zca * ((V * F * F) / (R * T)) * ((gamma_cai * cai
          * exp((zca * V * F) / (R * T)) - gamma_cao * cao) / (exp((zca * V * F)
          / (R * T)) - 1.0));

    if (fabs(V) < 1e-10)  /* Avoid denominator = 0 */
        icana_max = pna * F * (gamma_nai * nai - gamma_nao * nao);
    else
        icana_max = pna * zna * zna * ((V * F * F) / (R * T)) * ((gamma_nai
          * nai * exp((zna * V * F) / (R * T)) - gamma_nao * nao) / (exp((zna
          * V * F) / (R * T)) - 1.0));

    if (fabs(V) < 1e-10)  /* Avoid denominator = 0 */
        icak_max = pk * F * (gamma_ki * ki - gamma_ko * ko);
    else
        icak_max = pk * zk * zk * ((V * F * F) / (R * T)) * ((gamma_ki * ki
          * exp((zk * V * F) / (R * T)) - gamma_ko * ko) / (exp((zk * V * F)
          / (R * T)) - 1.0));

    fca = 1.0 / (1.0 + pow(cai / kmca, 2));

    ica = d * f * fca * ica_max;
    
    icana = d * f * fca * icana_max;
    
    icak = d * f * fca * icak_max;

    icat = ica + icana + icak;


/*------------------------------------------------------------------------------
The calculations for the time dependent K current (IIIc)
------------------------------------------------------------------------------*/

    ek = ((R * T) / F) * log((ko + pnak * nao) / (ki + pnak * nai));

    gk = 0.282 * sqrt(ko / 5.4);

    xi = 1.0 / (1.0 + exp((V - 56.26) / 32.1));

    a = 7.19e-5 * (V + 30.0) / (1.0 - exp(-0.148 * (V + 30.0)));
    b = 1.31e-4 * (V + 30.0) / (-1.0 + exp(0.0687 * (V + 30.0)));

    if (t*dt == 0)
        x = a / (a + b); /* Set the initial (steady state) value of the gate */
    else
        x = (a / (a + b)) - ((a / (a + b)) - x) * exp(-dt / (1.0 / (a + b)));

    ik = gk * xi * x * x * (V - ek);


/*------------------------------------------------------------------------------
The calculations for the time independent K current (IIId)
------------------------------------------------------------------------------*/

    ek1 = ((R * T) / F) * log(ko / ki);

    gk1 = 0.75 * sqrt(ko / 5.4);

    a = 1.02 / (1.0 + exp(0.2385 * (V - ek1 - 59.215)));

    b = (0.49124 * exp(0.08032 * (V - ek1 + 5.476)) + exp(0.06175 *
      (V - ek1 - 594.31))) / (1.0 + exp(-0.5143 * (V - ek1 + 4.753)));

    SteadyState = a / (a + b);

    ik1 = gk1 * SteadyState  * (V - ek1);


/*------------------------------------------------------------------------------
Calculate the plateau potassium current (IIIe)
------------------------------------------------------------------------------*/

    kp = 1.0 / (1.0 + exp((7.488 - V) / 5.98));	

    ikp = gkp * kp * (V - ek1);


/*------------------------------------------------------------------------------
Calculate the sodium-calcium exchanger current (IIIf)
------------------------------------------------------------------------------*/
    
    //printf("inaca before = %lf\n",inaca); // DEBUG
    inaca = knaca * (1.0 / (pow(kmna, 3) + pow(nao, 3))) * (1.0 /
      (kmca2 + cao)) * (1.0 / (1.0 + ksat * exp((eta - 1.0) * V * (F /
      (R * T))))) * (exp(eta * V * (F / (R * T))) * pow(nai, 3) * cao
      - exp((eta - 1.0) * V * (F / (R * T))) * pow(nao, 3) * cai);
	//printf("inaca after = %lf\n",inaca); // DEBUG


/*------------------------------------------------------------------------------
Calculate the sodium-potassium pump current (IIIg)
------------------------------------------------------------------------------*/

    sigma = (1.0 / 7.0) * (exp(nao / 67.3) - 1.0);

    fnak = 1.0 / (1.0 + 0.1245 * exp((-0.1 * V * F) / (R * T)) + 0.0365 * sigma
      * exp((-V * F) / (R * T)));

      inak = inak_max * fnak *(1.0 / (1.0 + pow(kmnai / nai, 1.5))) * (ko /
      (ko + kmko));


/*------------------------------------------------------------------------------
Calculate the nonspecific calcium activated current (IIIh)
------------------------------------------------------------------------------*/

    insk_max = pnsca * zk * zk * ((V * F * F) / (R * T)) * ((gamma_ki * ki
      * exp((zk * V * F) / (R * T)) - gamma_ko * ko) / (exp((zk * V * F)
      / (R * T)) - 1.0));

    insna_max = pnsca * zna * zna * ((V * F * F) / (R * T)) * ((gamma_nai * nai
      * exp((zna * V * F) / (R * T)) - gamma_nao * nao) / (exp((zna * V * F)
      / (R * T)) - 1.0));

    insk = insk_max * (1.0 / (1.0 + pow(kmnsca / cai, 3)));

    insna = insna_max * (1.0 / (1.0 + pow(kmnsca / cai, 3)));

    insca = insk + insna;

    /* The following equation is given in the manuscript but unused */
    /* ensca = ((R * T) / F) * log((ko + nao) / (ki + nai)); */


/*------------------------------------------------------------------------------
Calculate the sarcolemmal calcium pump (IIIi)
------------------------------------------------------------------------------*/

    ipca = ipca_max * (cai / (kmpca + cai));

/*------------------------------------------------------------------------------
Calculate the calcium background current (IIIj)
------------------------------------------------------------------------------*/
    ecan = ((R * T) / (2.0 * F)) * log(cao / cai);

    icab = gcab * (V - ecan);

/*------------------------------------------------------------------------------
Calculate the sodium background current (IIIk)
------------------------------------------------------------------------------*/

    inab = gnab * (V - ena);


/*------------------------------------------------------------------------------
Calculate the total time independent current (IIIl)
------------------------------------------------------------------------------*/

    iv = ik1 + ikp + ipca + inab + icab + inak;


/*------------------------------------------------------------------------------
Calculate the calcium induced calcium relase of JSR (Va)
------------------------------------------------------------------------------*/
    caiontot = ica + icab + ipca - 2 * inaca;
    
    dcaiontotnew = (caiontot - caiontotold) / dt;

    /* Approximate the time of the maximum value of dvdt by using Ca entry */
    if (V > -35 && dcaiontotnew > dcaiontot && test == 0) {
        test = 1;
	tdvdtmax = t*dt;
    }

    if (dvdt > 10.0 && t*dt > tdvdtmax + 10.0 && test == 1)
	test = 0;
	
    /* Calculate the change in intracellular Ca concetration for 2 ms after the time of dvdt max. */
    if ((tdvdtmax > 0.0) && (t*dt >= tdvdtmax) && (t*dt < (tdvdtmax + 2.00999999))){
        delta_cai2 += dcai;
    }

    /* Determine if Ca 2 ms after the time of dvdt max. is greater then the threshold */
    if (delta_cai2 > delta_caith && t*dt > tdvdtmax + 1.99999999 && t*dt < tdvdtmax + 2.00999999) {
        grelcicr_max = 60.0;
		tcicr = 0.0;
    }
    else if (delta_cai2 <= delta_caith && t*dt > tdvdtmax + 1.99999999 && t*dt < tdvdtmax + 2.00999999)
        grelcicr_max = 0.0;
	
    grelcicr = grelcicr_max * (delta_cai2 - delta_caith) / (kmrel + delta_cai2
      - delta_caith) * (1.0 - exp(-tcicr / tau_on)) * exp(-tcicr / tau_off);

    irelcicr = grelcicr * (cajsr - cai);
	
    tcicr += dt;


/*------------------------------------------------------------------------------
Calculate the ca release of JSR under ca overload conditions (Vb)
------------------------------------------------------------------------------*/

    if (csqn >= csqnth && tjsrol > 50.0) {
        greljsrol_max = 4.0;
        tjsrol = 0.0;
	//printf("Spontaneous release occured at %g milliseconds!\n", t*dt);
    }
    
    greljsrol = greljsrol_max * (1.0 - exp(-tjsrol / tau_on)) * exp(-tjsrol / tau_off);

    ireljsrol = greljsrol * (cajsr - cai);

    tjsrol += dt;


/*------------------------------------------------------------------------------
Calculate the calcium uptake and leakage of the NSR (Vd)
------------------------------------------------------------------------------*/

    iup = iup_max * cai / (cai + kmup);
	
    kleak = iup_max / cansr_max;
	
    ileak = kleak * cansr;


/*------------------------------------------------------------------------------
Calculate the translocation of ca ions from the NSR to the JSR (Ve)
------------------------------------------------------------------------------*/

    itr = (cansr - cajsr) / tau_tr;


/*------------------------------------------------------------------------------
Update the ion concentrations
------------------------------------------------------------------------------*/

    /* Update the myoplasmic Ca concentration */
    nai += dt * -((ina + icana + inab + insna + inaca * 3.0 + inak * 3.0)
      * acap) / (vmyo * zna * F);

    /* Update the myoplasmic K concentration */
    ki += dt * -((icak + ik + ik1 + ikp + insk - inak * 2.0) * acap)
      / (vmyo * zk * F);

    /* Update the NSR Ca concentration */
    cansr += dt * (iup - ileak - itr * vjsr / vnsr);

/* The calculation of Ca in the myoplasm and the JSR is complicated by the buffers
in these cell compartments. Therefore 2 equations are required in calculating these
concentrations. The first computes the change in Ca due to the influx / efflux of
Ca from the compartment. Following this the buffering of Ca is calcuated using
Steffensen's iterative method asssuming steady-state. */

    /* Change in myoplasmic Ca concentration due to influx / efflux of Ca */
    dcai = dt * -(((caiontot * acap) / (vmyo * zca * F)) + ((iup - ileak)
      * vnsr / vmyo) - (irelcicr * vjsr / vmyo) - (ireljsrol * vjsr / vmyo));
	
    /* Calculate myoplasmic Ca concentration using the Steffensen iterative method */
    // Calculate CaTot(cai)
    trpn = trpn_max * (cai / (cai + kmtrpn));
    cmdn = cmdn_max * (cai / (cai + kmcmdn));
    CaTot_cai  = cai + trpn + cmdn; /* Total myoplasmic Ca concentration = dissolved Ca in myoplasm + buffered Ca in myoplasm */
    
    Cond=0;
    value = CaTot_cai + dcai;
    P0 = cai + dcai * dt;// * 5.0;
	
    for (k = 1; k <= Max; k++) {
    	
    	// Calculate CaTot(P0)
    	trpn = trpn_max * (P0 / (P0 + kmtrpn));
    	cmdn = cmdn_max * (P0 / (P0 + kmcmdn));
    	gP0  = P0 + trpn + cmdn;
    	
		P1 = P0 + value - gP0; /* First new iterate */
		
		// Calculate CaTot(P1)
    	trpn = trpn_max * (P1 / (P1 + kmtrpn));
    	cmdn = cmdn_max * (P1 / (P1 + kmcmdn));
    	gP1  = P1 + trpn + cmdn;
		
		P2 = P1 + value - gP1; /* Second new iterate */

		D1 = (P1 - P0) * (P1 - P0); /* Form the differences */
		D2 = P2 - 2 * P1 + P0;

		if (fabs(D2) <= 1.e-14) { /* Check division by zero */
	    	Cond = 1;
	    	DP = P2 - P1;
	    	P3 = P2;
		} else {
	    	DP = D1 / D2;
	    	P3 = P0 - DP; /* Aitken's improvement */
		}
        RelErr = 2 * fabs(DP) / (fabs(P3) + Small);

		if (RelErr < Delta) /* Check for convergence */
	    	if (Cond != 1)
				Cond = 2;
		P0 = P3; /* Update iterate */
		

        if(Cond != 0) /* If condition is met terminate loop */
            break;
    }

    if (Cond == 0) {
		EXPECTED_ERROR("ERROR: cai Steffensen. The iterative procedure did not converge at time %lf.\n", t*dt);
	}

    cai=P0;

    /* Change in JSR Ca concentration due to influx / efflux of Ca */
    dcajsr = dt * (itr - irelcicr - ireljsrol);

    /* Calculate JSR Ca concentration using the Steffensen iterative method */
    // Calculate CaJSRTot(cajsr)
    csqn = csqn_max * (cajsr / (cajsr + kmcsqn));
	CaJSRTot_cajsr = cajsr + csqn;
    
    Cond = 0;
    value = CaJSRTot_cajsr + dcajsr;
    P0 = cajsr + dcajsr * dt;// * 5.0;


    for (k = 1; k <= Max; k++) {
    
    	// Calculate CaJSRTot(P0)
    	csqn = csqn_max * (P0 / (P0 + kmcsqn));
		gP0 = P0 + csqn;
    
		P1 = P0 + value - gP0; /* First new iterate */
		
		// Calculate CaJSRTot(P1)
    	csqn = csqn_max * (P1 / (P1 + kmcsqn));
		gP1 = P1 + csqn;
		
		P2 = P1 + value - gP1; /* Second new iterate */

		D1 = (P1 - P0) * (P1 - P0); /* Form the differences */
		D2 = P2 - 2 * P1 + P0;

		if (fabs(D2) <= 1.e-14) { /* Check division by zero */
	    	Cond = 1;
	    	DP = P2 - P1;
	    	P3 = P2;
		} else {
	    	DP = D1 / D2;
	    	P3 = P0 - DP; /* Aitken's improvement */
		}
        RelErr = 2 * fabs(DP) / (fabs(P3) + Small);

		if (RelErr < Delta) /* Check for convergence */
	    	if (Cond != 1)
				Cond = 2;
		P0 = P3; /* Update iterate */


        if(Cond != 0) /* If condition is met terminate loop */
            break;
    }

    if (Cond == 0) {
		EXPECTED_ERROR("ERROR: cajsr Steffensen. The iterative procedure did not converge at time %lf.\n", t*dt);
	}

    cajsr = P0;

/*******************************************************************************
Calculate the membrane potential
*******************************************************************************/
	/* Store the value of the previous membrane potential before calculating the new one */ 
    V_old = V;
	
	V += dt * -(ina + icat + ik + ik1 + ikp + inaca + inak + insca + ipca + icab + inab);

    /* Calculate the value of dvdt */
    dvdt = (V - V_old)/dt;

    caiontotold = caiontot;
    dcaiontot = dcaiontotnew;

    
/*------------------------------------------------------------------------------
Output values.
------------------------------------------------------------------------------*/    
//printf("%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",t*dt,V,ina,ica,ik,iv,inaca,cajsr,csqn,cai,trpn,cmdn,cansr,irelcicr,ireljsrol,iup,ileak,itr,d,f,fca,icat,icak,icana,ik1,ikp,inak,ipca,inab,icab);
