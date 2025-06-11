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

/* kernik2019 V-dependent functions including gate rates */
#define small (1e-10)
#define TI2AB(g) real alp_##g = g##_inf/tau_##g; real bet_##g = (1.0-g##_inf)/tau_##g;

/* INa m gate */
real m1 = 1.08045846384818006e+02; /* && (1/ms) */
real m2 = 1.31070157339409992e+01; /* && (mV) */
real m5 = 2.3269e-03; /* && */
real m6 = -7.91772628951300028; ; /* && (mV) */
real m3 = m5*m1; /* && */
real m4 = 1.0/((1.0/m2)+(1.0/m6));
real tau_m_const = 3.19775803839999970e-02; /* && (ms) */

real al_m = m1*exp(V/m2); /* && (1/ms) */
real be_m = m3*exp(V/m4); /* && (1/ms) */
real tau_m = (1.0/(al_m+be_m))+tau_m_const; /* && (ms) */
real m_inf = al_m/(al_m+be_m); /* && */
TI2AB(m);

/* INa h gate */
real h1 = 3.62659886399999999e-03; /* && (1/ms) */
real h2 = -1.98393588600259996e1; /* && (mV) */
real h5 = 9.66329497711473959e+03; /* && */
real h6 = 7.39550356461299963; /* && (mV) */
real h3 = (h5*h1); /* && (1/ms) */
real h4 = (1.0/((1.0/h2)+(1.0/h6))); /* &&mV */
real tau_h_const = 1.67331502516000014e-01; /* && (ms) */

real al_h = (h1*exp((V/h2))); /* && (1/ms) */
real be_h = (h3*exp((V/h4))); /* && (1/ms) */
real h_inf = (al_h/(al_h+be_h));/* && */
real tau_h = ((1.0/(al_h+be_h))+tau_h_const); /* && (ms) */ 
TI2AB(h);

/* INa j gate */
real j1 = 5.12257182000000044e-04; /* && (1/ms) */
real j2 = -6.65837555026519965e01; /* && (mV) */
real j5 = h5; /* && */
real j6 = h6; /* && (mV) */
real j3 = (j5*j1); /* && (1/ms) */
real j4 = (1.0/((1.0/j2)+(1.0/j6))); /* && (mV) */ 
real tau_j_const = 9.51088724962000032e-01; /* && (ms) */
real ina_j_a = (j1*(exp((V/j2))));

real al_j = (j1*exp((V/j2))); /* && (1/ms) */
real be_j = (j3*exp((V/j4))); /* && (1/ms) */
real j_inf = (al_j/(ina_j_a+be_j)); /* && */
real tau_j = ((1.0/(al_j+be_j))+tau_j_const); /* && (ms) */
TI2AB(j);

/* ICaL d gate */
real d1 = 1.29662941897219994e+01; /* && (1/ms) */
real d2 = 7.07914596471100044; /* && (mV) */
real d5 = 4.49094155069999987e-02; /* && */
real d6 = -6.90988036924199989; /* && (mV) */
real d3 = (d5*d1); /* && (1/ms) */
real d4 = (1.0/((1.0/d2)+(1.0/d6))); /* && (mV) */
real taud_const = 1.65824694683000007; /* && (ms) */

real al_d = (d1*exp((V/d2))); /* && (1/ms) */
real be_d = (d3*exp((V/d4))); /* && (1/ms) */
real d_inf = (al_d/(al_d+be_d)); /* && */
real tau_d = (((1.0/(al_d+be_d))+taud_const)); /* && (ms) */
TI2AB(d);

/* ICaL f gate */
real f1 = 5.1258925999999987e-04; /* && (1/ms) */
real f2 = -4.95057120338699974e1; /* && (mV) */
real f5 = 1.93121122351431995e+03; /* && */
real f6 = 5.73002749969900016; /* && (mV) */
real f3 = (f5*f1); /* && (1/ms) */
real f4 = (1.0/((1.0/f2)+(1.0/f6))); /* && (mV) */
real tauf_const = 1.00462559171102995e+02; /* && (ms) */

real al_f = (f1*exp((V/f2))); /* && (1/ms) */
real be_f = (f3*exp((V/f4))); /* && (1/ms) */
real f_inf = (al_f/(al_f+be_f)); /* && */
real tau_f = (((1.0/(al_f+be_f))+tauf_const)); /* && (ms) */
TI2AB(f); 

/* Ito */
real r1 = 5.53614181712999975e-02; /* && (1/ms) */
real r2 = 1.16842023429669002e+01; /* && (mV) */
real r5 = 3.98918108037750008; /* && */
real r6 = -1.10471393012032006e1; /* && (mV) */
real r3 = (r5*r1); /* && (1/ms) */
real r4 = (1.0/((1.0/r2)+(1.0/r6))); /* && (mV) */
real tau_r_const = 6.96758421171499998e-01; /* && (ms) */

real al_r = (r1*exp((V/r2))); /* && (1/ms) */
real be_r = (r3*exp((V/r4))); /* && (1/ms) */
real r_inf = (al_r/(al_r+be_r)); /* && */
real tau_r = ((1.0/(al_r+be_r))+tau_r_const); /* && (ms) */
TI2AB(r)

real s1 = 3.44230944300000013e-04; /* && (1/ms) */
real s2 = -1.76344722898096009e1; /* && (mV) */
real s5 = 1.86760536909694991e+02;/* && */
real s6 = 8.18093387332270083; /* && (mV) */
real s3 = (s5*s1); /* && (1/ms) */
real s4 = (1.0/((1.0/s2)+(1.0/s6))); /* (mV) */
real tau_s_const = 1.12244577239468999e+01; /* && (ms) */ 

real al_s = (s1*exp((V/s2))); /* && (1/ms) */
real be_s = (s3*exp((V/s4))); /* && (1/ms) */
real s_inf = (al_s/(al_s+be_s)); /* && */
real tau_s = ((1.0/(al_s+be_s))+tau_s_const); /* && (ms) */
TI2AB(s)

/* ICaT*/
real d_icat_inf = (1.0/(1.0+exp(((V+26.3)/-6.0)))); /* && */
real tau_d_icat = (1.0/((1.068*exp(((V+26.3)/30.0)))+(1.068*exp(((V+26.3)/-30.0))))); /* && (ms) */
TI2AB(d_icat)

real f_icat_inf = (1.0/(1.0+exp(((V+61.7)/5.6)))); /* && */
real tau_f_icat = (1.0/((0.0153*exp((-(V+61.7)/83.3)))+(0.015*exp(((V+61.7)/15.38))))); /* && (ms) */
TI2AB(f_icat)

/* real diff_d_icat = ((d_icat_inf - d_icat)/d_icat_tau); /\* && (1/ms) *\/ */
/* real diff_f_icat = ((f_icat_inf - f_icat)/f_icat_tau); /\* && (1/ms) *\/ */

/* Ifunny */
real xF1 = 5.78970000000000002e-07; /* && (1/ms) */
real xF2 = -1.45897121702000003e1; /* && (mV) */
real xF5 = 2.00866502378844016e+04; /* && */
real xF6 = 1.02023528452800001e+01; /* && (mV) */
real xF3 = (xF5*xF1); /* && (1/ms) */
real xF4 = (1.0/((1.0/xF2)+(1.0/xF6))); /* && (mV) */
real xF_const = 2.39452913465299986e+01; /* && (ms) */

/* real Na_frac = (NatoK_ratio/(NatoK_ratio+1.0)); /\* && *\/ */

real al_Xf = (xF1*exp((V/xF2))); /* && (1/ms) */
real be_Xf = (xF3*exp((V/xF4))); /* && (1/ms) */
real Xf_inf = (al_Xf/(al_Xf+be_Xf)); /* && */
real tau_Xf = ((1.0/(al_Xf+be_Xf))+xF_const); /* && (ms) */
TI2AB(Xf)

/* IK1 */
real xK11 = 4.77994972217041014e-01; /* && (1/ms) */
real xK12 = 2.72427558793486995e+01; /* && (mV) */
real xK13 = 4.92502331781412028; /* && (mV) */
real xK14 = 8.72223760006881932; /* && (mV) */
real xK15 = 5.66361974998243980e+01; /* && (mV) */

real al_XK1 = (xK11*exp(((V+xK13)/xK12))); /* && (1/ms) */
real be_XK1 = ((1.0*exp(((V+xK15)/xK14)))); /* && (1/ms) */
real XK1_inf = (al_XK1/(al_XK1+be_XK1)); /* && */

/* IKr */
real Xr1_1 = 5.74885237435000026e-03; /* && (1/ms) */
real Xr1_2 = 1.36234926362576001e+01; /* && (mV) */
real Xr1_5 = 4.76305711818360011e-02; /* && */
real Xr1_6 = -7.06808742965549008; /* && (mV) */
real Xr1_3 = (Xr1_5*Xr1_1); /* && (1/ms) */
real Xr1_4 = (1.0/((1.0/Xr1_2)+(1.0/Xr1_6))); /* && (mV) */
real tau_1_offset = 50.0; /* && (ms) */

real al_Xr1 = (Xr1_1*exp((V/Xr1_2))); /* && (1/ms) */
real be_Xr1 = (Xr1_3*exp((V/Xr1_4))); /* && (1/ms) */
real Xr1_inf = (al_Xr1/(al_Xr1+be_Xr1)); /* && */
real tau_Xr1 = ((1.0/(al_Xr1+be_Xr1))+tau_1_offset); /* && (ms) */
TI2AB(Xr1)
/* real diff_Xr1 = ((Xr1_inf - Xr1)/tau_Xr1); */

real Xr2_1 = 1.24566405268270002e-02; /* && (1/ms) */
real Xr2_2 = -2.59944581644376989e1; /* && (mV) */ 
real Xr2_5 = 3.73426331501040991e+01; /* && */
real Xr2_6 = 2.20919642353902006e+01; /* && (mV) */
real Xr2_3 = (Xr2_5*Xr2_1); /* && (1/ms) */
real Xr2_4 = (1.0/((1.0/Xr2_2)+(1.0/Xr2_6))); /* && (mV) */
real tau_2_offset = 0.0; /* && (ms) */
 
real al_Xr2 = (Xr2_1*exp((V/Xr2_2))); /* && (1/ms) */ 
real be_Xr2 = (Xr2_3*exp((V/Xr2_4))); /* && (1/ms) */
real Xr2_inf = (al_Xr2/(al_Xr2+be_Xr2)); /* && */ 
real tau_Xr2 = ((1.0/(al_Xr2+be_Xr2))+tau_2_offset); /* && (ms) */
TI2AB(Xr2)
/* real diff_Xr2 = ((Xr2_inf - Xr2)/tau_Xr2); */ 

/* IKs */
real ks1 = 1.16558447999999992e-03; /* && (1/ms) */ 
real ks2 = 6.67268386758935958e+04; /* && (mV) */ 
real ks5 = 2.80458908250000027e-01; /* && */ 
real ks6 = -1.88669715729099998e1; /* && (mV) */ 
real ks3 = (ks5*ks1); /* && (1/ms) */ 
real ks4 = (1.0/((1.0/ks2)+(1.0/ks6))); /* && (mV) */ 
real tauks_const = 4.74115000000000034e-06; /* && (ms) */  

real al_Xs = (ks1*exp((V/ks2))); /* && (1/ms) */
real be_Xs = (ks3*exp((V/ks4))); /* && (1/ms) */ 
real Xs_inf = (al_Xs/(al_Xs+be_Xs)); /* && */ 
real tau_Xs = ((1.0/(al_Xs+be_Xs))+tauks_const); /* && (ms) */ 
TI2AB(Xs)
/* diff_Xs = ((Xs_inf - Xs)/tau_Xs); */
