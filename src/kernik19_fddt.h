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

/* ODE right-hand sides for Kernik et al 2019, except V-controlled gates */


#define small (1e-10)
#define TI2AB(g) real alp_##g = g##_inf/tau_##g; real bet_##g = (1.0-g##_inf)/tau_##g;

real E_Ca = ((0.5*RTF)*log((Cao/Cai))); /* && (mV) */
real E_K = (RTF*log((Ko/Ki))); /* && (mV) */
real E_Na = (RTF*log((Nao/Nai))); /* && (mV) */ 

/* cell geometry */
real V_tot_tenT = (Vc_tenT+VSR_tenT); /* && (L) */
real V_SR = (V_tot*(VSR_tenT/V_tot_tenT)); /* && (L) */
real Vc = (V_tot*(Vc_tenT/V_tot_tenT)); /* && (L) */

/* INa */
real i_Na = (g_Na*cub(m)*h*j*(V-E_Na)); /* && (mV/ms) */

/* ICaL*/
real al_fCa = (1.0/(1.0+eighth(((scale*Cai)/0.000325)))); /* && */
real be_fCa = (0.1/(1.0+exp((((scale*Cai) - 0.0005)/0.0001)))); /* && */
real gamma_fCa = (0.2/(1.0+exp((((scale*Cai) - 0.00075)/0.0008))));  /* && */
real fCa_inf = ((((al_fCa+be_fCa)+gamma_fCa)+0.23)/1.46); /* && */
real k_fca = iif(((fCa_inf>fCa) && (V>-60.0)), 0.0, 1.0); /* && */
/* NB can't rewrite this in standard form for TI2AB macro */
/* because k_fca=0 in parts of the phase plane */
/* real diff_fCa = ((k_fca*(fCa_inf-fCa))/tau_fCa); */ /* && (1/ms) */
real alp_fCa = ((k_fca*(fCa_inf))/tau_fCa); /* && (1/ms) */
real bet_fCa = ((k_fca*(1.0-fCa_inf))/tau_fCa); /* && (1/ms) */

real p_CaL_shannonTot = ((p_CaL_shannonCa+p_CaL_shannonNa)+p_CaL_shannonK); /* && */ 
real p_CaL_shannonNap = (p_CaL_shannonNa/p_CaL_shannonTot); /* && */
real p_CaL_Na = (p_CaL_shannonNap*p_CaL); /* && (L/F ms) */
real p_CaL_shannonKp = (p_CaL_shannonK/p_CaL_shannonTot); /* && */
real p_CaL_K = (p_CaL_shannonKp*p_CaL); /* && (L/F ms) */
real p_CaL_shannonCap = (p_CaL_shannonCa/p_CaL_shannonTot); /* && */
real p_CaL_Ca = (p_CaL_shannonCap*p_CaL); /* && (L/F ms) */

real ibarca = (((((p_CaL_Ca*4.0)*V)*FFRT)*(((0.341*Cai)*exp(((2.0*V)*FRT))) - (0.341*Cao)))/(exp(((2.0*V)*FRT)) - 1.0)); /* && (mV/ms) */
real ibark = ((((p_CaL_K*V)*FFRT)*(((0.75*Ki)*exp((V*FRT))) - (0.75*Ko)))/(exp((V*FRT)) - 1.0)); /* && (mV/ms) */
real ibarna = ((((p_CaL_Na*V)*FFRT)*(((0.75*Nai)*exp((V*FRT))) - (0.75*Nao)))/(exp((V*FRT)) - 1.0)); /* && (mV/ms) */
real i_CaL_Ca = (((ibarca*d)*f)*fCa); /* && (mV/ms) */
real i_CaL_K = (((ibark*d)*f)*fCa); /* && (mV/ms) */
real i_CaL_Na = (((ibarna*d)*f)*fCa); /* && (mV/ms) */
real i_CaL = ((i_CaL_Ca+i_CaL_Na)+i_CaL_K); /* && */

/* IbCa */
real i_b_Ca = g_b_Ca*(V - E_Ca); /* && (mV/ms) */

/* IbNa */
real i_b_Na = g_b_Na*(V - E_Na); /* && (mV/ms) */

/* Jleak */ 
real j_leak = ((Ca_SR - Cai)*V_leak); /* && (mM/ms) */

/* INaCa */
real i_NaCa = (kNaCa*(((exp(((gamma*V)*FRT))*cub(Nai)*Cao) - (((exp((((gamma - 1.0)*V)*FRT))*cub(Nao)*Cai)*alpha)))/(((cub(KmNai)+cub(Nao))*(KmCa+Cao))*(1.0+(Ksat*exp((((gamma - 1.0)*V)*FRT))))))); /* && (mV/ms) */ 

/* INak */
real i_NaK = (((PNaK*Ko)*Nai)/(((Ko+Km_K)*(Nai+Km_Na))*((1.0+(0.1245*exp(((-0.1*V)*FRT))))+(0.0353*exp((-V*FRT)))))); /* && (mV/ms) */

/* IpCa */
real i_PCa = ((g_PCa*Cai)/(Cai+KPCa)); 

/* Ito */
real i_to = ((g_to*r*s)*(V - E_K)); /* && (mV/ms) */

/* ICaT */ 
real i_CaT = (((g_CaT*d_icat)*f_icat)*(V - E_Ca)); /* && (mV/ms) */

/* ifunny */
real i_fNa = (((Na_frac*g_f)*Xf)*(V - E_Na)); /* && (mV/ms) */
real i_fK = ((((1.0 - Na_frac)*g_f)*Xf)*(V - E_K)); /* && (mV/ms) */
real i_f = (i_fNa+i_fK); /* && (mV/ms) */

/* IK1 */
real i_K1 = (((g_K1*sqrt((Ko/5.4)))*XK1_inf)*(V - E_K)); /* && (mV/ms) */

/* IKr */
real i_Kr = ((((g_Kr*sqrt((Ko/5.4)))*Xr1)*Xr2)*(V - E_K)); /* && (mV/ms) */

/* IKs */
real i_Ks = ((g_Ks*sqr(Xs))*(V - E_K)); /* && (mV/ms) */

/* Jrel */
real kCaSR = (MaxSR - ((MaxSR - MinSR)/(1.0+pow((ec50SR/Ca_SR),2.5)))); /* && */
real kiSRCa = (kiCa*kCaSR); /* && (1/(mM ms)) */
real koSRCa = (koCa/kCaSR); /* && 1/(mM^2 ms) */

real RI = (((1.0 - R1) - O) - I); /* && */

real diff_I = (kiSRCa*Cai)*O - (kim)*I - (kom)*I + (koSRCa*sqr(Cai))*RI; /* && (1/ms) */
real diff_O = (koSRCa*sqr(Cai))*R1 - (kom)*O - (kiSRCa*Cai)*O + (kim)*I; /* && (1/ms) */
real diff_R1 = (kim)*RI - (kiSRCa*Cai)*R1 - (koSRCa*sqr(Cai))*R1 + (kom)*O; /* && (1/ms) */

real j_rel = ((((ks*O)*(Ca_SR - Cai))*(V_SR))/Vc); /* && (mM/ms) */


/* Jup */
real j_up = (VmaxUp/(1.0+(sqr(Kup)/sqr(Cai)))); /* && (mM/ms) */

/* Nai */
real diff_Nai = ((-Cm/(F*Vc))*((((((i_Na)+i_b_Na)+i_fNa)+(3.0*i_NaK))+(3.0*i_NaCa))+i_CaL_Na)); /* && (mM/ms) */

/* Cai */
real Cai_bufc = (1.0/(1.0+((Buf_C*Kbuf_C)/((Cai+Kbuf_C)*(Cai+Kbuf_C))))); /* && */
real diff_Cai = (Cai_bufc*(((-j_up+j_leak)+j_rel) - ((Cm/((2.0*Vc)*F))*(((((i_CaL_Ca)+i_CaT)+i_b_Ca)+i_PCa) - (2.0*i_NaCa))))); /* && (mM/ms) */

/* Ki */
real diff_Ki = ((-Cm/(F*Vc))*(((((((i_K1)+i_to)+i_Kr)+i_Ks)+i_fK) - (2.0*i_NaK))+i_CaL_K)); /* && (mM/ms) */

/* CaSR */
real Ca_SR_bufSR = (1.0/(1.0+((Buf_SR*Kbuf_SR)/sqr(Ca_SR+Kbuf_SR)))); /* && */
real diff_Ca_SR = (((Ca_SR_bufSR*Vc)/V_SR)*((j_up - j_rel) - j_leak)); /* && (mM/ms) */


/* V */
real diff_V = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_CaT+i_NaK+i_Na+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca);	/* (mV/ms) */

