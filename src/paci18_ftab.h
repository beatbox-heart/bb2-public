/**
 * Copyright (C) (2010-2024) Vadim Biktashev, Irina Biktasheva et al. 
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

/* PPCPASH2018 V-dependent functions including gate rates */
#define small (1e-10)
#define TI2AB(g) real alp_##g = g##_inf/tau_##g; real bet_##g = (1.0-g##_inf)/tau_##g;

real m_inf	= 1.0/pow((1.0+exp((-V-34.1)/5.9)),(1.0/3.0));	/* (1) */
real al_m	= 1.0/(1.0+exp((-V-60.0)/5.0));			/* (1) */
real be_m	= 0.1/(1.0+exp((V+35.0)/5.0))+0.1/(1.0+exp((V-50.0)/200.0));	/* (1) */
real tau_m	= al_m*be_m;					/* (ms) */
TI2AB(m);

real h_inf	= 1.0/sqrt(1.0+exp((V+72.1)/5.7));		/* (1) */
real al_h	= 0.057*exp(-(V+80.0)/6.8);			/* (1) */
real be_h	= 2.7*exp(0.079*V)+3.1e5*exp(0.3485*V);		/* (1) */
real tau_h	= iif( (V < -38.5) , ( 1.5/((al_h+be_h)) ) , ( 1.5*1.6947 ) );	/* (ms) */
TI2AB(h);

real j_inf	= 1.0/sqrt(1.0+exp((V+72.1)/5.7));		/* (1) */
real al_j	= iif ( (V < -40.0) , ( 
                (-25428.0*exp(0.2444*V)-6.948e-6*exp(-0.04391*V))*(V+37.78)/(1.0+exp(0.311*(V+79.23))) 
              ) , ( 
                0.0 
		    ) );                                    	/* (1) */
real be_j	= iif ( (V < -40.0) , ( 
                 ((0.02424*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))))) 
               ) , ( 
                 0.6*exp((0.057)*V)/(1+exp(-0.1*(V+32))) 
		     ) );					/* (1) */
real tau_j	= 7.0/((al_j+be_j));				/* (ms) */
TI2AB(j);

real mL_inf	= 1/(1+exp(-(V+42.85)/(5.264)));		/* (1) */
real al_mL	= 1/(1+exp((-60-V)/5));				/* (1) */
real be_mL	= 0.1/(1+exp((V+35)/5))+0.1/(1+exp((V-50)/200));/* (1) */
real tau_mL	= al_mL*be_mL;					/* (ms) */
TI2AB(mL);

real hL_inf	= 1/(1+exp((V+Vh_hLate)/(7.488)));		/* (1) */
real tau_hL	= 200;						/* (ms) */
TI2AB(hL);

real Xf_inf	= 1.0/(1.0+exp((V+77.85)/5.0));			/* (1) */
real tau_Xf	= 1900.0/(1.0+exp((V+15.0)/10.0));		/* (ms) */
TI2AB(Xf);

real d_inf	= 1.0/(1.0+exp(-(V+9.1)/7.0));			/* (1) */
real tau_d	= ((0.25+1.4/(1.0+exp((-V-35.0)/13.0)))*1.4/(1.0+exp((V+5.0)/5.0))+1.0/(1.0+exp((-V+50.0)/20.0)))*1.0;	/* (ms) */
TI2AB(d);

real f2_inf	= 0.33+0.67/(1.0+exp((V+32.0)/4.0));		/* (1) */
real constf2	= 1.0;						/* (1) */
real tau_f2	= (600.0*exp(-sqr(V+25.0)/170.0)+31.0/(1.0+exp((25.0-V)/10.0))+16.0/(1.0+exp((30.0+V)/10.0)))*constf2;	/* (ms) */
TI2AB(f2);

real q_inf	= 1.0/(1.0+exp((V+53.0)/13.0));			/* (1) */
real tau_q	= (6.06+39.102/(0.57*exp(-0.08*(V+44.0))+0.065*exp(0.1*(V+45.93))));	/* (ms) */
TI2AB(q);

real r_inf	= 1.0/(1.0+exp(-(V-22.3)/18.75));		/* (1) */
real tau_r	= (2.75352+14.40516/(1.037*exp(0.09*(V+30.61))+0.369*exp(-0.12*(V+23.84))));	/* (ms) */
TI2AB(r);

real Xs_inf	= 1.0/(1.0+exp((-V-20.0)/16.0));		/* (1) */
real tau_Xs	= 1.0*(1100.0/sqrt(1.0+exp((-10.0-V)/6.0)))*(1.0/(1.0+exp((-60.0+V)/20.0)));	/* (ms) */
TI2AB(Xs);

real L0		= 0.025;					/* (1) */
real logL0	= log(L0);					/* (1) */
real Q		= 2.3;						/* (1) */
real V_half	= (RTF/Q*(logL0+4*log((1.0+Cao/0.58)/(1.0+Cao/2.6)))-19.0); /* (mV) */ /* NB this demands that Cao is constant rather than parameter */
real Xr1_inf	= 1.0/(1.0+exp((V_half-V)/4.9));		/* (1) */
real tau_Xr1	= 1.0*(450.0/(1.0+exp((-45.0-V)/10.0)))*(6.0/(1.0+exp((30.0+V)/11.5)));	/* (ms) */
TI2AB(Xr1);

real Xr2_inf	= 1.0/(1.0+exp((V+88.0)/50.0));		/* (1) */
real tau_Xr2	= 1.0*(3.0/(1.0+exp((-60.0-V)/20.0)))*(1.12/(1.0+exp((-60.0+V)/20.0)));	/* (ms) */
TI2AB(Xr2);
