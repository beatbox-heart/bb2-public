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

/* ODE right-hand sides for Paci et al 2018, except V-controlled gates */

#define small (1e-10)
#define TI2AB(g) real alp_##g = g##_inf/tau_##g; real bet_##g = (1.0-g##_inf)/tau_##g;

real E_Na = RTF*log(Nao/Nai);			/* (mV) */
real E_Ca = 0.5*RTF*log(Cao/Cai);		/* (mV) */
real E_K  = RTF*log(Ko/Ki);			/* (mV) */
real PkNa = 0.03;				/* (1) */
real E_Ks = RTF*log((Ko+PkNa*Nao)/(Ki+PkNa*Nai));	/* (mV) */

real i_Na	= g_Na*cub(m)*h*j*(V-E_Na);	/* (mV/ms) */
real i_NaL	= g_NaL*cub(mL)*hL*(V-E_Na);	/* (mV/ms) */
real i_f	= g_f*Xf*(V-E_f);		/* (mV/ms) */
real i_fNa	= 0.42*g_f*Xf*(V-E_Na);		/* (mV/ms) */
real i_CaL	= g_CaL*4.0*F*V/RTF*(Cai*exp(2.0*V/RTF)-0.341*Cao)/(exp(2.0*V/RTF)-1.0)*d*f1*f2*fCa;	/* (mV/ms) */
real i_to	= g_to*(V-E_K)*q*r;		/* (mV/ms) */
real i_Ks	= g_Ks*(V-E_Ks)*sqr(Xs)*(1.0+0.6/(1.0+pow((3.8e-5/Cai),1.4)));/* (mV/ms)  */
real i_Kr	= (g_Kr*(V-E_K)*Xr1*Xr2*sqrt(Ko/5.4));/* (mV/ms)  */
real al_K1	= 3.91/(1.0+exp(0.5942*(V-E_K-200.0)));		/* (1) */
real be_K1	= (-1.509*exp(0.0002*(V-E_K+100.0))+exp(0.5886*(V-E_K-10.0)))/(1.0+exp(0.4547*(V-E_K)));	/* (1) */
real XK1_inf	= al_K1/(al_K1+be_K1);		/* (1) */
real i_K1	= g_K1*XK1_inf*(V-E_K)*sqrt(Ko/5.4);	/* (mV/ms)  */
real i_NaCa	= kNaCa*
  (exp(gamma*V/RTF)*cub(Nai)*Cao-exp((gamma-1.0)*V/RTF)*cub(Nao)*Cai*alpha)/
  ((cub(KmNai)+cub(Nao))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*(V/RTF))));		/* (mV/ms) */
real i_NaK = PNaK*Ko/(Ko+Km_K)*Nai/(Nai+Km_Na)/(1.0+0.1245*exp(-0.1*V/RTF)+0.0353*exp(-V/RTF));	/* (mV/ms) */
real i_PCa = g_PCa*Cai/(Cai+KPCa);		/* (mV/ms) */
real i_b_Na = g_b_Na*(V-E_Na);			/* (mV/ms) */
real i_b_Ca = g_b_Ca*(V-E_Ca);			/* (mV/ms) */
real j_up   = j_up_max/(1.0+sqr(Kup/Cai));	/* (mM/ms) */
real j_leak = j_leak_max*(Ca_SR-Cai);		/* (mM/ms) */
real RyRSRCass = (1 - 1/(1 +  exp((Ca_SR-0.3)/0.1)));	/* (1)  */
real j_rel = j_rel_max*RyRSRCass*RyRo*RyRc*(Ca_SR-Cai);	/* (mM/ms) */


real f1_inf	= 1.0/(1.0+exp((V+26.0)/3.0));			/* (1) */
real constf1	= iif( (f1_inf-f1 > 0.0) , (1.0+1433.0*(Cai-50.0e-6)), (1.0));	/* (1) */
real tau_f1	= (20.0+1102.5*exp(-sqr(sqr(V+27.0)/15.0))+200.0/(1.0+exp((13.0-V)/10.0))+180.0/(1.0+exp((30.0+V)/10.0)))*constf1;	/* (ms) */
TI2AB(f1);

real fCa_inf	= (1.0/(1.0+eighth(Cai/0.0006)) + 0.1/(1.0+exp((Cai-0.0009)/0.0001)) + 0.3/(1.0+exp((Cai-0.00075)/0.0008)))/1.3156;	/* (1) */
real constfCa	= iif( ((V > -60.0) && (fCa_inf > fCa)) , (0.0) , (1.0) );	/* (1) lit mV */
real tau_fCa	= 2;						/* (ms) */
TI2AB(fCa);

real RyRa_inf = RyRa1-RyRa2/(1 + exp((Cai/(1e-3)-(RyRahalf))/0.0082));	/* (1)  */
real tau_RyRa = 1.0e3;				/* (ms) */
TI2AB(RyRa)

real RyRo_inf = (1 - 1/(1 +  exp((Cai/(1e-3)-(RyRa+ RyRohalf))/0.003)));/* (1)  */
real tau_RyRo = iif( (RyRo_inf>= RyRo), 18.75, 0.1*18.75);		/* (ms)  */
TI2AB(RyRo)

real RyRc_inf = (1/(1 + exp((Cai/(1e-3)-(RyRa+RyRchalf))/0.001)));	/* (1)  */
real tau_RyRc = iif( (RyRc_inf>= RyRc), 2*87.5, 87.5 );			/* (ms)  */
TI2AB(RyRc)

real Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/sqr(Cai+Kbuf_C));		/* (1) */
real Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/sqr(Ca_SR+Kbuf_SR));		/* (1) */

real diff_Nai    = -CFV*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa);	/* (mM/ms) */
real diff_Cai    = Cai_bufc*(j_leak-j_up+j_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*CFV*0.5);	/* (mM/ms) */
real diff_Ca_SR  = Ca_SR_bufSR*(Vc/Vsr)*(j_up-(j_rel+j_leak));		/* (mM/ms) */
real diff_V = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca);	/* (mV/ms) */

