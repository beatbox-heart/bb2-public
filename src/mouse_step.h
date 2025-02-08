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
the formulas and the RHS for KLYZMOUSE
*/

/*
Model variables
*/

// double v;
double sstime;
int  current_trace_counter = 0, output_counter = 0;

double vclamp;
double total_current;

/*
Ion concentrations and reversal potentials
*/

//double cai;             // mM
//double casub;           // mM 

double ena,eca,ek,eks;

/*
Ist variables.
*/

double ist;
double qa, qi, tauqa,tauqi,alphaqa,betaqa,alphaqi,betaqi;
//double dst; // ist gating variables
//double fst;

/*
Ib variables
*/

double ibca,ibna,ibk,ib;

/*
IK1 variables
*/

double ik1 ;
double xk1inf;

/*
ICaT variables
*/

double icat;
//double dt;
//double ft;
double dt_inf, tau_dt, ft_inf, tau_ft;

/*
IKr variables
*/

double ikr_act_inf,tau_ikr_act;
double ikr;
double ikr_inact_inf,tau_ikr_inact,tau_ikr_inact2;
// ,ikr_inact,ikr_inact2;

/*
IKs variables
*/

// double iks_act
double iks_act_inf,tau_iks_act;
double iks;

/*
ICaL 1.2 and 1.3 parameters.
*/

double alpha_dl, beta_dl, tau_dl, alpha_fl, beta_fl, tau_fl;
// double fl12, dl12, 
double dl12_inf, fl12_inf,ical12;
// double fl13, dl13, 
double dl13_inf, fl13_inf,ical13;

// double fca; 
double fca_inf ;
double taufca ;

/*
INa Nav1.1 (TTXS) and Nav1.5 (TTXR) variables,
*/

double ina_ttxr, ina_ttxs; // ina: Nav 1.1 and Nav 1.5
double m3_inf_ttxr, h_inf_ttxr;
double m3_inf_ttxs, h_inf_ttxs;
double m_inf_ttxr,m_inf_ttxs, hs,hsr;
//double m_ttxr,h_ttxr,j_ttxr;
//double m_ttxs,h_ttxs,j_ttxs;
double tau_m,tau_h,tau_j,tau_mr,tau_hr,tau_jr;
double delta_ttxr_m,delta_ttxr_h,delta_ttxr_j;
double fna;

/*
If variables
*/

double ih, ihk, ihna, ih_1, ih_2, ih_4, /* y_1_2, */ y_4, tau_y_1_2, tau_y_4, y_inf, y_1_2_inf, y_4_inf,ih_1_k, ih_2_k, ih_4_k, ih_1_na, ih_2_na, ih_4_na;

/*
INaK variables
*/

double inak;

/*
Inaca variables.
*/
double inaca;

double di,doo,k43,k12,k14,k41,k34,k21,k23,k32,x1,x2,x3,x4;

/*
Isus variables
*/

double r_inf, /* r,*/ tau_r, isus;

/*
Ito variables
*/

// double q,
double q_inf,tau_q;
double ito;

/*
Ionic homeostasis
*/

double ca_flux;

/*
Calcium diffusion
*/

double Jcadif;

double Jrel, Jup, Jtr;
//,carel,caup;

double dFtc; // ,Ftc;
double dFtmc; // ,Ftmc,Ftmm;
double dFtmm; // ,Ftmm;
double dFcms; // ,Fcms;
double dFcmi; // ,Fcmi;
double dFcq; // ,Fcq;

//double nai=8.0;         // mM 
//double ki=140;
double FRT, RTF;

/*
RyR markov chain
*/

// double resting, open, inactivated, resting_inactivated;
double kcasr,kosrca,kisrca;

double nai_tot, ki_tot;
double dvdt;

double soicr_ss;

FRT = F/(R*T);

ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);

/*Ist********************************************************************/
/* Shinagawa  */

qa = 1.0/(1.0 + exp(-(v+67.0)/5.0));

alphaqa = 1.0/(0.15*exp(-(v)/11.0)+0.2*exp(-(v)/700.0));
betaqa  =  1.0/(16.0*exp((v)/8.0)+15.0*exp((v)/50.0));
tauqa = 1.0/(alphaqa + betaqa);

alphaqi = 0.15*1.0/(3100.0*exp((v+10.0)/13.0)+700.3*exp((v+10.0)/70.0));
betaqi =  0.15*1.0/(95.7*exp(-(v+10.0)/10.0) + 50.0*exp(-(v+10.0)/700.0)) + 0.000229/(1+exp(-(v+10.0)/5.0));
qi = alphaqi/(alphaqi + betaqi);
tauqi = 1.0/(alphaqi + betaqi);

dst = dst + ddt*((qa-dst)/tauqa);
fst = fst + ddt*((qi-fst)/tauqi);

ist = gst*dst*fst*(v - eist);

/* Ib ************************************************************************/
  ibna = gbna*(v - ena);
  ibca = gbca*(v - eca);
  ibk  =  gbk*(v - ek);
  ib = (ibna + ibca + ibk);
  
/*IK1**********************************************************************/
xk1inf = 1.0/(1.0 + exp(0.070727*(v - ek)));
ik1 = gk1*xk1inf*(ko/(ko + 0.228880))*(v - ek);

/**ICaT Cav3.1**************************************************************/
/* ICaT */

tau_dt = 1.0/(1.068*exp((v + 26.3)/30.0) + 1.068*exp(-(v + 26.3)/30.0)); // Kurata
dt_inf = 1.0/(1.0+exp(-(v + 26.0)/6.0));
dt = dt + ddt*((dt_inf - dt)/tau_dt);

tau_ft = 1.0/(0.0153*exp(-(v+61.7)/83.3)+0.015*exp((v+61.7)/15.38)); // Kurata/Maltsev.
ft_inf = 1.0/(1.0+exp((v + 61.7)/5.6)); // rabbit SS, since the Mangoni data is 10 mV low, and very steep.
ft = ft + ddt*((ft_inf - ft)/tau_ft);

icat = gcat*ft*dt*(v - ecat);           // nA

/*Ikr*************************************************************************/

ikr_act_inf = 1.0/(1.0 + exp(-(v+21.173694)/9.757086)); // sacred // ISO shifts it to the left
tau_ikr_act = 0.699821/(0.003596*exp((v)/15.339290) + 0.000177*exp(-(v)/25.868423)); // this fits the Q10 variation.
ikr_act = ikr_act + ddt*(ikr_act_inf-ikr_act)/tau_ikr_act;
     
/* SS inactivation */

ikr_inact_inf = 1.0/(1.0 + exp((v+20.758474-4.0)/(19.0)));
tau_ikr_inact = 0.2+0.9*1.0/(0.1*exp(v/54.645)+0.656*exp(v/106.157));
ikr_inact = ikr_inact + ddt*(ikr_inact_inf - ikr_inact)/tau_ikr_inact;
ikr = gkr*ikr_act*ikr_inact*(v - ek);

/**IKs********************************************************************/

iks_act_inf = 1.0/(1.0 + exp(-(v-20.876040)/11.852723));
tau_iks_act =  1000.0/(13.097938/(1.0 + exp(-(v-48.910584)/10.630272)) + exp(-(v)/35.316539)); // Zhang model. The other models are similar.
iks_act = iks_act + ddt*(iks_act_inf - iks_act)/tau_iks_act; // Ding/Maatsura
iks = gks*iks_act*iks_act*(v - eks);

/*ICaL*******************************************************************/


if(fabs(v)<=0.001){
 alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)+408.173;
}
else
if(fabs(v+35.0)<=0.001){
 alpha_dl  = 70.975-84.9*v/(exp(-0.208*v)-1.0);
}
else
if(fabs(v)>0.001&&fabs(v+35.0)>0.001){
 alpha_dl  = -28.39*(v+35.0)/(exp(-(v+35.0)/2.5)-1.0)-84.9*v/(exp(-0.208*v)-1.0);
}

if(fabs(v-5.0)<=0.001)
 beta_dl   = 28.575;
else
if(fabs(v-5.0)>0.001)
 beta_dl   = 11.43*(v-5.0)/(exp(0.4*(v-5.0))-1.0);


tau_dl  = 2000.0/(alpha_dl +beta_dl);

/* Cav 1.3 */
/* Cav 1.3 SS values */

dl13_inf = 1.0/(1+exp(-(v+dl13vhalf)/6.0));
fl13_inf = 1.0/(1+exp((v+fl13vhalf)/7.3));
tau_fl = (7.4 + 45.77*exp(-0.5*(v+28.1)*(v+28.1)/(11*11)));

dl13 = dl13 + ddt*(dl13_inf - dl13)/tau_dl;
fl13 = fl13 + ddt*(fl13_inf - fl13)/tau_fl;

/* Cav 1.2 */

dl12_inf = 1.0/(1+exp(-(v+dl12vhalf)/5.0)); // Mangoni 2006: according to his email
fl12_inf = 1.0/(1+exp((v+fl12vhalf)/4.6)); // Mangoni 2006 35.9

dl12 = dl12 + ddt*(dl12_inf - dl12)/tau_dl;
fl12 = fl12 + ddt*(fl12_inf - fl12)/tau_fl;

/* fca */

fca_inf = kmfca/(kmfca+casub);
taufca = fca_inf/alpha_fca;
fca = fca + ddt*(fca_inf - fca)/taufca;

 ical12 = gcal12*fl12*dl12*fca*(v-ecal);
 ical13 = gcal13*fl13*dl13*fca*(v-ecal);
    
/**INa**********************************************************************/

fna = (9.52e-02*exp(-6.3e-2*(v+34.4))/(1+1.66*exp(-0.225*(v+63.7))))+8.69e-2; // 1-fna = proportion of fast inactivation, fna = contribution of slow

m3_inf_ttxr = 1.0/(1.0 + exp(-(v+45.213705)/7.219547));
h_inf_ttxr = 1.0/(1.0 + exp(-(v+62.578120 )/(-6.084036)));
m3_inf_ttxs = 1.0/(1.0 + exp(-(v+36.097331-5.0)/5.0)); // Mangoni takes -29.
h_inf_ttxs = 1.0/(1.0 + exp((v+56.0)/3.0));

m_inf_ttxr = pow(m3_inf_ttxr,0.333);
m_inf_ttxs = pow(m3_inf_ttxs,0.333);

tau_m = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_h = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_j = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);

m_ttxs = m_ttxs + ddt*(m_inf_ttxs - m_ttxs)/tau_m;
h_ttxs = h_ttxs + ddt*(h_inf_ttxs - h_ttxs)/tau_h;
j_ttxs = j_ttxs + ddt*(h_inf_ttxs - j_ttxs)/tau_j;
 
hs = (1.0-fna)*h_ttxs+fna*j_ttxs;

tau_mr = 1000.0*((0.6247e-03/(0.832*exp(-0.335*(v+56.7))+0.627*exp(0.082*(v+65.01))))+0.0000492);
tau_hr = 1000.0*(((3.717e-06*exp(-0.2815*(v+17.11)))/(1+0.003732*exp(-0.3426*(v + 37.76))))+0.0005977);
tau_jr = 1000.0*(((0.00000003186*exp(-0.6219*(v+18.8)))/(1+0.00007189*exp(-0.6683*(v+34.07))))+0.003556);

m_ttxr = m_ttxr + ddt*(m_inf_ttxr - m_ttxr)/tau_mr;
h_ttxr = h_ttxr + ddt*(h_inf_ttxr - h_ttxr)/tau_hr;
j_ttxr = j_ttxr + ddt*(h_inf_ttxr - j_ttxr)/tau_jr;

hsr = (1.0-fna)*h_ttxr+fna*j_ttxr;

if(fabs(v)>0.005)
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*(F*F/(R*T))*((exp((v-ena)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxs= gna_ttxs*m_ttxs*m_ttxs*m_ttxs*hs*nao*F*((exp((v-ena)*F/(R*T))-1.0));


if(fabs(v)>0.005)
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*(F*F/(R*T))*((exp((v-enattxr)*F/(R*T))-1.0)/(exp(v*F/(R*T))-1.0))*v;
else
ina_ttxr = gna_ttxr*m_ttxr*m_ttxr*m_ttxr*hsr*nao*F*((exp((v-enattxr)*F/(R*T))-1.0));


/**If**************************************************************************/

y_inf = 1.0/(1.0 + exp((v+if_vhalf)/16.3));
tau_y_1_2 = 1.5049/(exp(-(v+590.3)*0.01094)+ exp((v-85.1)/17.2));
y_1_2 = y_1_2 + ddt*(y_inf - y_1_2)/tau_y_1_2;

ihk  = 0.6167*gh*y_1_2*(v - ek);
ihna = 0.3833*gh*y_1_2*(v - ena);

ih = (ihk + ihna);

/*Ito*************************************************************************/

q_inf = 1.0/(1.0+exp((v+49.0)/13.0));
tau_q = (6.06 + 39.102/(0.57*exp(-0.08*(v+44.0))+0.065*exp(0.1*(v+45.93))))/0.67; 
q = q + ddt*((q_inf-q)/tau_q);
r_inf = 1.0/(1.0+exp(-(v-19.3)/15.0));
tau_r = (2.75+14.40516/(1.037*exp(0.09*(v+30.61))+0.369*exp(-0.12*(v+23.84))))/0.303;
r = r + ddt*((r_inf-r)/tau_r);
ito = gto*q*r*(v-ek);
 
/*Isus***********************************************************************/

isus = gsus*r*(v-ek);

/*Inak***********************************************************************/

inak = inakmax*(pow(ko,1.2)/(pow(kmkp,1.2)+pow(ko,1.2)))*(pow(nai,1.3)/(pow(kmnap,1.3)+pow(nai,1.3)))/(1.0+exp(-(v-ena+120.0)/30.0));

/****iNaCa*******************************************************************/

// INaCa
       di=1+(casub/Kci)*(1+exp(-Qci*v*FRT)+nai/Kcni)+(nai/K1ni)*(1+(nai/K2ni)*(1+nai/K3ni));
       doo=1+(cao/Kco)*(1+exp(Qco*v*FRT))+(nao/K1no)*(1+(nao/K2no)*(1+nao/K3no));
       k43=nai/(K3ni+nai);
       k12=(casub/Kci)*exp(-Qci*v*FRT)/di;
       k14=(nai/K1ni)*(nai/K2ni)*(1+nai/K3ni)*exp(Qn*v*FRT/2.0)/di;
       k41=exp(-Qn*v*FRT/2.0);
       k34=nao/(K3no+nao);
       k21=(cao/Kco)*exp(Qco*v*FRT)/doo;
       k23=(nao/K1no)*(nao/K2no)*(1+nao/K3no)*exp(-Qn*v*FRT/2.0)/doo;
       k32=exp(Qn*v*FRT/2);
       x1=k34*k41*(k23+k21)+k21*k32*(k43+k41);
       x2=k43*k32*(k14+k12)+k41*k12*(k34+k32);
       x3=k43*k14*(k23+k21)+k12*k23*(k43+k41);
       x4=k34*k23*(k14+k12)+k21*k14*(k34+k32);
   
 inaca = kNaCa*(k21*x2-k12*x1)/(x1+x2+x3+x4);   

/*****************************************************************************/

/* Ion fluxes */
/* complete these equations */
ca_flux = (ical12+ical13+icat-2.0*inaca+ibca)/(2.0*F);

/*****************************************************************************/

 /* Ca2+ difusion */
 Jcadif = (casub - cai)/tdifca;

 /* Ca handling in SR (mM/ms) */
 kcasr = maxsr - (maxsr - minsr)/(1.0 + pow(eca50sr/carel,hsrr));

 kosrca = koca/kcasr;

 kisrca = kica*kcasr;

resting = resting + ddt*(kim*resting_inactivated - kisrca*casub*resting - kosrca*casub*casub*resting + kom*open);
open    = open    + ddt*(kosrca*casub*casub*resting - kom*open - kisrca*casub*open + kim*inactivated);
inactivated = inactivated + ddt*(kisrca*casub*open - kim*inactivated - kom*inactivated + kosrca*casub*casub*resting_inactivated);
resting_inactivated = resting_inactivated + ddt*(kom*inactivated - kosrca*casub*casub*resting_inactivated - kim*resting_inactivated + kisrca*casub*resting);

// Sung et al. 2011
soicr_ss = soicr_max/(1.0 + exp((soicr_thresh - carel)/soicr_slope));
soicr    = soicr_ss - (soicr_ss - soicr)*exp(-ddt/soicr_tau);
Jrel     = (ks+soicr_switch*soicr)*open*(carel - casub);

Jup = Pup*(pow(cai/pumpkmf,pumphill) - pow(caup/pumpkmr,pumphill))/(1.0 + pow(cai/pumpkmf,pumphill) + pow(caup/pumpkmr,pumphill));

Jtr  = (caup - carel)/Ttr;

// Ca buffering flux // these are the derivatives of the F's
dFtc  = kfTC*cai*(1.0-Ftc)-kbTC*Ftc;
dFtmc = kfTMC*cai*(1.0-Ftmc-Ftmm)-kbTMC*Ftmc;
dFtmm = kfTMM*Mgi*(1.0-Ftmc-Ftmm)-kbTMM*Ftmm;
dFcms = kfCM*casub*(1.0-Fcms)-kbCM*Fcms;
dFcmi = kfCM*cai*(1.0-Fcmi)-kbCM*Fcmi;
dFcq  = kfCQ*carel*(1.0-Fcq)-kbCQ*Fcq;

Ftc = Ftc + ddt*dFtc;                   //	    % dfTnCa/dt
Ftmc = Ftmc + ddt*dFtmc;                  //	% dfTnMgCa/dt
Ftmm = Ftmm + ddt*dFtmm;                  //	% dfTnMgMg/dt
Fcms = Fcms + ddt*dFcms;                 //	    % dfCMs/dt
Fcmi = Fcmi + ddt*dFcmi;                 //	    % dfCMi/dt
Fcq = Fcq + ddt*dFcq;                    //	    % dfCQ/dt

 casub = casub + ddt*((-ca_flux+Jrel*vrel)/vsub-Jcadif-ConcCM*dFcms);
 cai = cai + ddt*((Jcadif*vsub-Jup*vup)/vi - (ConcCM*dFcmi + ConcTC*dFtc + ConcTMC*dFtmc)); 



carel = carel + ddt*(Jtr - Jrel - ConcCQ*dFcq);
caup = caup + ddt*(Jup-Jtr*vrel/vup);

/***************************************************************************/

total_current = ih+ina_ttxr+ina_ttxs+ical12+ical13+iks+ikr+ik1+ist+ib+icat+inak+isus+inaca+ito;
dvdt = - total_current/capacitance;

v = v  + ddt*dvdt;
 
ena = (R*T/F)*log(nao/nai);
ek  = (R*T/F)*log(ko/ki);
eks = ((R*T)/F)*log((ko+0.12*nao)/(ki+0.12*nai));
eca = (R*T/(2*F))*log(cao/casub);
 
nai_tot = ihna+ina_ttxr+ina_ttxs+3.0*inak+3.0*inaca+ist+ibna;
ki_tot = ihk+iks+ikr+ik1+ibk-2.0*inak+isus+ito;

nai = nai - ddt*(nai_tot)/(F*vi);
ki = ki - ddt*(ki_tot)/(F*vi);

