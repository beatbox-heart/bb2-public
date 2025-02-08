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

/* CRN98 V-dependent functions including gate rates */
#define small (1e-10)
#define TI2AB(g) real alp_##g = inf_##g/tau_##g; real bet_##g = (1.0-inf_##g)/tau_##g;

real gkur_rel = 0.005+0.05/(1+exp((V-15)/-13)); 
real gkr_rel  = 1.0/(1+exp((V+15)/22.4)); 
real gk1_rel  = 1.0/(1+exp(0.07*(V+80)));
real fnak1 = 0.1245*exp(-0.1*V*F/(R*T));
real fnak2 = 0.0365*exp(-V*F/(R*T));
real expVgammalrF_RT=exp(V*gammalr*F/(R*T));
real expVgammalr1F_RT=exp(V*(gammalr-1)*F/(R*T));
/* ina m gate */
real alp_m = (fabs(V+47.13) < small)?(3.2):(0.32*(V+47.13)/(1-exp(-0.1*(V+47.13))));
real bet_m = 0.08*exp(-V/11);
/* ina h gate */
real alp_h = (V >= -40.0)?0.0:(0.135*exp((V+80)/-6.8));
real bet_h = (V >= -40.0)?(1.0/(0.13*(1+exp((V+10.66)/-11.1)))):(3.56*exp(0.079*V)+3.1e5*exp(0.35*V));
/* ina j gate */
real alp_j = (V >= -40.0)?0.0:((-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))*(V+37.78)/
			       (1+exp(0.311*(V+79.23))));
real bet_j = (V >= -40.0)?(0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)))):(0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14))));
/* oa ito gate, corrected for 37 deg */
/* voltage shift due to junctional potential and effect of Cd++ */
#define vshift (-10)
real a_oa = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));
real b_oa = 0.65/(2.5+exp((V-vshift+72.0)/17.0)); 
real tau_oa = 1/(a_oa+b_oa)/Tfac; 
real inf_oa = 1/(1+exp((V-vshift+10.47)/-17.54));
TI2AB(oa);
#undef vshift
/* oi and ois ito gate, corrected for 37 deg */
/* voltage shift due to junctional potential and effect of Cd++ */
#define vshift (-10)
real a_oi = 1/(18.53+exp((V-vshift+103.7)/10.95));
real b_oi = 1/(35.56+exp((V-vshift-8.74)/-7.44));
real tau_oi = 1/(a_oi+b_oi)/Tfac; 
real inf_oi = 1/(1+exp((V-vshift+33.1)/5.3));
TI2AB(oi);
#undef vshift
/* ua ikur gate, corrected for 37 deg */
/* voltage shift due to junctional potential and effect of Cd++ */
#define vshift (-10)
real a_ua = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));
real b_ua = 0.65/(2.5+exp((V-vshift+72.0)/17.0)); 
real tau_ua = 1/(a_ua+b_ua)/Tfac; 
real inf_ua = 1/(1+exp((V-vshift+20.3)/-9.6));
TI2AB(ua);
#undef vshift
/* ui ikur gate, corrected for 37 deg */
/* voltage shift due to junctional potential and effect of Cd++ */
#define vshift (-10)
real a_ui = 1/(exp((V-vshift-195)/-28)+21);
real b_ui = 1/(exp((V-vshift-168)/-16));
real tau_ui = 1/(a_ui+b_ui)/Tfac; 
real inf_ui = 1/(1+exp((V-vshift-109.45)/27.48));
TI2AB(ui);
#undef vshift
/* xr ikr gate */
real a_xr = (fabs(V+14.1) < small)?0.0015:(0.0003*(V+14.1)/(1-exp((V+14.1)/-5)));
real b_xr = (fabs(V-3.3328) < small)?3.7836118e-4:(0.000073898*(V-3.3328)/(exp((V-3.3328)/5.1237)-1));
real tau_xr = 1/(a_xr+b_xr); 
real inf_xr = 1/(1+exp((V+14.1)/-6.5));
TI2AB(xr);
/* xs ikr gate */
real a_xs = (fabs(V-19.9) < small)?0.00068:(0.00004*(V-19.9)/(1-exp((V-19.9)/-17)));
real b_xs = (fabs(V-19.9) < small)?0.000315:(0.000035*(V-19.9)/(exp((V-19.9)/9)-1));
real tau_xs = 0.5/(a_xs+b_xs);   /* reduced by 50% as described in manuscript */
real inf_xs = sqrt(1/(1+exp((V-19.9)/-12.7)));
TI2AB(xs);
/* icaL d gate */
real a_d = 1/(1+exp((V+10)/-6.24));
real tau_d = a_d*((fabs(V+10) < small)?4.579:((1-exp((V+10)/-6.24))/(0.035*(V+10))));
real inf_d = 1/(1+exp((V+10)/-8));
TI2AB(d);
/* icaL f gate */
real inf_f = exp(-(V+28)/6.9)/(1+exp(-(V+28)/6.9));
real tau_f = 1.5*2*3/(0.0197*exp(-0.0337*0.0337*(V+10)*(V+10))+0.02);
TI2AB(f);
/* one of SR gates */
real inf_ww = 1-1/(1+exp(-(V-40)/17.0));
real tau_ww = (fabs(V-7.9) < small)?(6*0.2/1.3):(6.0*(1-exp(-(V-7.9)/5.0))/(1+0.3*exp(-(V-7.9)/5.0))/(V-7.9));
TI2AB(ww);

#undef TI2AB
#undef small
