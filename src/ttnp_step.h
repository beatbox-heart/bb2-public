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


//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

double alpha_d;   // dimensionless (in L_type_Ca_current_d_gate)
double beta_d;   // dimensionless (in L_type_Ca_current_d_gate)
double d_inf;   // dimensionless (in L_type_Ca_current_d_gate)
double gamma_d;   // millisecond (in L_type_Ca_current_d_gate)
double tau_d;   // millisecond (in L_type_Ca_current_d_gate)
double f2_inf;   // dimensionless (in L_type_Ca_current_f2_gate)
double tau_f2;   // millisecond (in L_type_Ca_current_f2_gate)
double fCass_inf;   // dimensionless (in L_type_Ca_current_fCass_gate)
double tau_fCass;   // millisecond (in L_type_Ca_current_fCass_gate)
double f_inf;   // dimensionless (in L_type_Ca_current_f_gate)
double tau_f;   // millisecond (in L_type_Ca_current_f_gate)
double i_CaL;   // picoA_per_picoF (in L_type_Ca_current)
double i_b_Ca;   // picoA_per_picoF (in calcium_background_current)
double Ca_i_bufc;   // dimensionless (in calcium_dynamics)
double Ca_sr_bufsr;   // dimensionless (in calcium_dynamics)
double Ca_ss_bufss;   // dimensionless (in calcium_dynamics)
double O;   // dimensionless (in calcium_dynamics)
double i_leak;   // millimolar_per_millisecond (in calcium_dynamics)
double i_rel;   // millimolar_per_millisecond (in calcium_dynamics)
double i_up;   // millimolar_per_millisecond (in calcium_dynamics)
double i_xfer;   // millimolar_per_millisecond (in calcium_dynamics)
double k1;   // per_millimolar2_per_millisecond (in calcium_dynamics)
double k2;   // per_millimolar_per_millisecond (in calcium_dynamics)
double kcasr;   // dimensionless (in calcium_dynamics)
double i_p_Ca;   // picoA_per_picoF (in calcium_pump_current)
double alpha_h;   // per_millisecond (in fast_sodium_current_h_gate)
double beta_h;   // per_millisecond (in fast_sodium_current_h_gate)
double h_inf;   // dimensionless (in fast_sodium_current_h_gate)
double tau_h;   // millisecond (in fast_sodium_current_h_gate)
double alpha_j;   // per_millisecond (in fast_sodium_current_j_gate)
double beta_j;   // per_millisecond (in fast_sodium_current_j_gate)
double j_inf;   // dimensionless (in fast_sodium_current_j_gate)
double tau_j;   // millisecond (in fast_sodium_current_j_gate)
double alpha_m;   // dimensionless (in fast_sodium_current_m_gate)
double beta_m;   // dimensionless (in fast_sodium_current_m_gate)
double m_inf;   // dimensionless (in fast_sodium_current_m_gate)
double tau_m;   // millisecond (in fast_sodium_current_m_gate)
double i_Na;   // picoA_per_picoF (in fast_sodium_current)
double alpha_K1;   // dimensionless (in inward_rectifier_potassium_current)
double beta_K1;   // dimensionless (in inward_rectifier_potassium_current)
double i_K1;   // picoA_per_picoF (in inward_rectifier_potassium_current)
double xK1_inf;   // dimensionless (in inward_rectifier_potassium_current)
// double i_Stim;   // picoA_per_picoF (in membrane)
double i_p_K;   // picoA_per_picoF (in potassium_pump_current)
double alpha_xr1;   // dimensionless (in rapid_time_dependent_potassium_current_Xr1_gate)
double beta_xr1;   // dimensionless (in rapid_time_dependent_potassium_current_Xr1_gate)
double tau_xr1;   // millisecond (in rapid_time_dependent_potassium_current_Xr1_gate)
double xr1_inf;   // dimensionless (in rapid_time_dependent_potassium_current_Xr1_gate)
double alpha_xr2;   // dimensionless (in rapid_time_dependent_potassium_current_Xr2_gate)
double beta_xr2;   // dimensionless (in rapid_time_dependent_potassium_current_Xr2_gate)
double tau_xr2;   // millisecond (in rapid_time_dependent_potassium_current_Xr2_gate)
double xr2_inf;   // dimensionless (in rapid_time_dependent_potassium_current_Xr2_gate)
double i_Kr;   // picoA_per_picoF (in rapid_time_dependent_potassium_current)
double E_Ca;   // millivolt (in reversal_potentials)
double E_K;   // millivolt (in reversal_potentials)
double E_Ks;   // millivolt (in reversal_potentials)
double E_Na;   // millivolt (in reversal_potentials)
double alpha_xs;   // dimensionless (in slow_time_dependent_potassium_current_Xs_gate)
double beta_xs;   // dimensionless (in slow_time_dependent_potassium_current_Xs_gate)
double tau_xs;   // millisecond (in slow_time_dependent_potassium_current_Xs_gate)
double xs_inf;   // dimensionless (in slow_time_dependent_potassium_current_Xs_gate)
double i_Ks;   // picoA_per_picoF (in slow_time_dependent_potassium_current)
double i_b_Na;   // picoA_per_picoF (in sodium_background_current)
double i_NaCa;   // picoA_per_picoF (in sodium_calcium_exchanger_current)
double i_NaK;   // picoA_per_picoF (in sodium_potassium_pump_current)
double r_inf;   // dimensionless (in transient_outward_current_r_gate)
double tau_r;   // millisecond (in transient_outward_current_r_gate)
double s_inf;   // dimensionless (in transient_outward_current_s_gate)
double tau_s;   // millisecond (in transient_outward_current_s_gate)
double i_to;   // picoA_per_picoF (in transient_outward_current)

   // time: time (millisecond)

   i_CaL = g_CaL*Y[0]*Y[3]*Y[1]*Y[2]*4.0*(Y[11]-15.0)*pow(F, 2.0)/(R*T)*(0.25*Y[6]*exp(2.0*(Y[11]-15.0)*F/(R*T))-Ca_o)/(exp(2.0*(Y[11]-15.0)*F/(R*T))-1.0);
   d_inf = 1.0/(1.0+exp((-8.0-Y[11])/7.5));
   alpha_d = 1.4/(1.0+exp((-35.0-Y[11])/13.0))+0.25;
   beta_d = 1.4/(1.0+exp((Y[11]+5.0)/5.0));
   gamma_d = 1.0/(1.0+exp((50.0-Y[11])/20.0));
   tau_d = 1.0*alpha_d*beta_d+gamma_d;
   dY[0] = (d_inf-Y[0])/tau_d;
   f2_inf = 0.67/(1.0+exp((Y[11]+35.0)/7.0))+0.33;
   tau_f2 = 562.0*exp(-pow(Y[11]+27.0, 2.0)/240.0)+31.0/(1.0+exp((25.0-Y[11])/10.0))+80.0/(1.0+exp((Y[11]+30.0)/10.0));
   dY[1] = (f2_inf-Y[1])/tau_f2;
   fCass_inf = 0.6/(1.0+pow(Y[6]/0.05, 2.0))+0.4;
   tau_fCass = 80.0/(1.0+pow(Y[6]/0.05, 2.0))+2.0;
   dY[2] = (fCass_inf-Y[2])/tau_fCass;
   f_inf = 1.0/(1.0+exp((Y[11]+20.0)/7.0));
   tau_f = 1102.5*exp(-pow(Y[11]+27.0, 2.0)/225.0)+200.0/(1.0+exp((13.0-Y[11])/10.0))+180.0/(1.0+exp((Y[11]+30.0)/10.0))+20.0;
   dY[3] = (f_inf-Y[3])/tau_f;
   E_Ca = 0.5*R*T/F*log(Ca_o/Y[5]);
   i_b_Ca = g_bca*(Y[11]-E_Ca);
   kcasr = max_sr-(max_sr-min_sr)/(1.0+pow(EC/Y[4], 2.0));
   k1 = k1_prime/kcasr;
   O = k1*pow(Y[6], 2.0)*Y[7]/(k3+k1*pow(Y[6], 2.0));
   i_rel = V_rel*O*(Y[4]-Y[6]);
   i_up = Vmax_up/(1.0+pow(K_up, 2.0)/pow(Y[5], 2.0));
   i_leak = V_leak*(Y[4]-Y[5]);
   i_xfer = V_xfer*(Y[6]-Y[5]);
   k2 = k2_prime*kcasr;
   dY[7] = -k2*Y[6]*Y[7]+k4*(1.0-Y[7]);
   Ca_i_bufc = 1.0/(1.0+Buf_c*K_buf_c/pow(Y[5]+K_buf_c, 2.0));
   Ca_sr_bufsr = 1.0/(1.0+Buf_sr*K_buf_sr/pow(Y[4]+K_buf_sr, 2.0));
   Ca_ss_bufss = 1.0/(1.0+Buf_ss*K_buf_ss/pow(Y[6]+K_buf_ss, 2.0));
   i_p_Ca = g_pCa*Y[5]/(Y[5]+K_pCa);
   i_NaCa = K_NaCa*(exp(gamma*Y[11]*F/(R*T))*pow(Y[16], 3.0)*Ca_o-exp((gamma-1.0)*Y[11]*F/(R*T))*pow(Na_o, 3.0)*Y[5]*alpha)/((pow(Km_Nai, 3.0)+pow(Na_o, 3.0))*(Km_Ca+Ca_o)*(1.0+K_sat*exp((gamma-1.0)*Y[11]*F/(R*T))));
   dY[5] = Ca_i_bufc*((i_leak-i_up)*V_sr/V_c+i_xfer-1.0*(i_b_Ca+i_p_Ca-2.0*i_NaCa)*Cm/(2.0*1.0*V_c*F));
   dY[4] = Ca_sr_bufsr*(i_up-(i_rel+i_leak));
   dY[6] = Ca_ss_bufss*(-1.0*i_CaL*Cm/(2.0*1.0*V_ss*F)+i_rel*V_sr/V_ss-i_xfer*V_c/V_ss);
   E_Na = R*T/F*log(Na_o/Y[16]);
   i_Na = g_Na*pow(Y[10], 3.0)*Y[8]*Y[9]*(Y[11]-E_Na);
   h_inf = 1.0/pow(1.0+exp((Y[11]+71.55)/7.43), 2.0);

   if (Y[11] < -40.0)
      alpha_h = 0.057*exp(-(Y[11]+80.0)/6.8);
   else
      alpha_h = 0.0;

   if (Y[11] < -40.0)
      beta_h = 2.7*exp(0.079*Y[11])+310000.0*exp(0.3485*Y[11]);
   else
      beta_h = 0.77/(0.13*(1.0+exp((Y[11]+10.66)/-11.1)));

   tau_h = 1.0/(alpha_h+beta_h);
   dY[8] = (h_inf-Y[8])/tau_h;
   j_inf = 1.0/pow(1.0+exp((Y[11]+71.55)/7.43), 2.0);

   if (Y[11] < -40.0)
      alpha_j = (-25428.0*exp(0.2444*Y[11])-6.948e-6*exp(-0.04391*Y[11]))*(Y[11]+37.78)/1.0/(1.0+exp(0.311*(Y[11]+79.23)));
   else
      alpha_j = 0.0;

   if (Y[11] < -40.0)
      beta_j = 0.02424*exp(-0.01052*Y[11])/(1.0+exp(-0.1378*(Y[11]+40.14)));
   else
      beta_j = 0.6*exp(0.057*Y[11])/(1.0+exp(-0.1*(Y[11]+32.0)));

   tau_j = 1.0/(alpha_j+beta_j);
   dY[9] = (j_inf-Y[9])/tau_j;
   m_inf = 1.0/pow(1.0+exp((-56.86-Y[11])/9.03), 2.0);
   alpha_m = 1.0/(1.0+exp((-60.0-Y[11])/5.0));
   beta_m = 0.1/(1.0+exp((Y[11]+35.0)/5.0))+0.1/(1.0+exp((Y[11]-50.0)/200.0));
   tau_m = 1.0*alpha_m*beta_m;
   dY[10] = (m_inf-Y[10])/tau_m;
   E_K = R*T/F*log(K_o/Y[12]);
   alpha_K1 = 0.1/(1.0+exp(0.06*(Y[11]-E_K-200.0)));
   beta_K1 = (3.0*exp(0.0002*(Y[11]-E_K+100.0))+exp(0.1*(Y[11]-E_K-10.0)))/(1.0+exp(-0.5*(Y[11]-E_K)));
   xK1_inf = alpha_K1/(alpha_K1+beta_K1);
   i_K1 = g_K1*xK1_inf*sqrt(K_o/5.4)*(Y[11]-E_K);

/*
   if ((time-floor(time/stim_period)*stim_period >= stim_start) && (time-floor(time/stim_period)*stim_period <= stim_start+stim_duration))
      i_Stim = -stim_amplitude;
   else
      i_Stim = 0.0;
*/


   i_to = g_to*Y[17]*Y[18]*(Y[11]-E_K);
   i_Kr = g_Kr*sqrt(K_o/5.4)*Y[13]*Y[14]*(Y[11]-E_K);
   E_Ks = R*T/F*log((K_o+P_kna*Na_o)/(Y[12]+P_kna*Y[16]));
   i_Ks = g_Ks*pow(Y[15], 2.0)*(Y[11]-E_Ks);
   i_NaK = P_NaK*K_o/(K_o+K_mk)*Y[16]/(Y[16]+K_mNa)/(1.0+0.1245*exp(-0.1*Y[11]*F/(R*T))+0.0353*exp(-Y[11]*F/(R*T)));
   i_b_Na = g_bna*(Y[11]-E_Na);
   i_p_K = g_pK*(Y[11]-E_K)/(1.0+exp((25.0-Y[11])/5.98));
   dY[11] = -1.0/1.0*(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_b_Na+i_NaCa+i_b_Ca+i_p_K+i_p_Ca/*+i_Stim*/);
   dY[12] = -1.0*(i_K1+i_to+i_Kr+i_Ks+i_p_K/*+i_Stim*/-2.0*i_NaK)/(1.0*V_c*F)*Cm;
   xr1_inf = 1.0/(1.0+exp((-26.0-Y[11])/7.0));
   alpha_xr1 = 450.0/(1.0+exp((-45.0-Y[11])/10.0));
   beta_xr1 = 6.0/(1.0+exp((Y[11]+30.0)/11.5));
   tau_xr1 = 1.0*alpha_xr1*beta_xr1;
   dY[13] = (xr1_inf-Y[13])/tau_xr1;
   xr2_inf = 1.0/(1.0+exp((Y[11]+88.0)/24.0));
   alpha_xr2 = 3.0/(1.0+exp((-60.0-Y[11])/20.0));
   beta_xr2 = 1.12/(1.0+exp((Y[11]-60.0)/20.0));
   tau_xr2 = 1.0*alpha_xr2*beta_xr2;
   dY[14] = (xr2_inf-Y[14])/tau_xr2;
   xs_inf = 1.0/(1.0+exp((-5.0-Y[11])/14.0));
   alpha_xs = 1400.0/sqrt(1.0+exp((5.0-Y[11])/6.0));
   beta_xs = 1.0/(1.0+exp((Y[11]-35.0)/15.0));
   tau_xs = 1.0*alpha_xs*beta_xs+80.0;
   dY[15] = (xs_inf-Y[15])/tau_xs;
   dY[16] = -1.0*(i_Na+i_b_Na+3.0*i_NaK+3.0*i_NaCa)/(1.0*V_c*F)*Cm;
   r_inf = 1.0/(1.0+exp((20.0-Y[11])/6.0));
   tau_r = 9.5*exp(-pow(Y[11]+40.0, 2.0)/1800.0)+0.8;
   dY[17] = (r_inf-Y[17])/tau_r;
   s_inf = 1.0/(1.0+exp((Y[11]+20.0)/5.0));
   tau_s = 85.0*exp(-pow(Y[11]+45.0, 2.0)/320.0)+5.0/(1.0+exp((Y[11]-20.0)/5.0))+3.0;
   dY[18] = (s_inf-Y[18])/tau_s;

