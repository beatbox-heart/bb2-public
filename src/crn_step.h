/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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
  This code is based on the M.Courtemanche's `error-free' electronic 
  version of equations and parameters, with the following copyright notice: 
*/

/*********************************************************************/
/* This is a distribution version of the code for the CRN
   human atrial cell model as described in:
   Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying
   human atrial action potential properties: insights from a mathematical
   model. Am J Physiol (1998);275:H301-H321. */

/* This code is copyrighted by Marc Courtemanche, 1998. */

/* It may be used, modified and distributed with my permission only, 
   with proper credit given and reference to the original manuscript 
   provided. */
/*********************************************************************/

/* computed quantities */
  double Ek,Ena,Eca,ina,icaL,ik,ibna,ibk,ibca,inak,ik1;
  double ito,gkur,ikur,ikr,iks;
  double inaca,naidot,kidot;
  double caidot,irel;
  double itr,iup,iupleak;
  double fnak,caflux,icap,sigma;
  /* utility variables */
  double a,b,tau,inf,vshift;

  /* compute Ek */
  Ek = 26.71*log(kc/ki);
  
  /* compute Ena */
  Ena = 26.71*log(nac/nai);
  
  /* compute Eca */
  Eca = 13.35*log(cac/cai);
  
  /* compute sodium current */
  ina = Cm*gna*m*m*m*h*j*(V-Ena);
  
  /* compute transient outward current */
  ito = Cm*gto*oa*oa*oa*oi*(V-Ek); 
  
  /* compute ultra-rapid potassium current */
  gkur = gkur_scale*(0.005+0.05/(1+exp((V-15)/-13))); 
  ikur = Cm*gkur*ua*ua*ua*ui*(V-Ek);
  
  /* compute the rapid delayed outward rectifier K current */
  ikr = Cm*gkr*xr*(V-Ek)/(1+exp((V+15)/22.4)); 
  
  /* compute the slow delayed outward rectifier K current */
  iks = Cm*gks*xs*xs*(V-Ek);
  
  /* compute calcium current */
  icaL = Cm*gcaL*d*f*fca*(V-ErL);

  /* update the fca gate immediately */
  inf = 1/(1+pow(cai/0.00035,1.0)); 
  tau = 2.0; 
  fca = inf + (fca-inf)*exp(-dt/tau); 

  /* compute time independent potassium current */
  ik1 = Cm*gk1*(V-Ek)/(1+exp(0.07*(V+80))); 
  
  /* compute ibna background current */
  ibna = Cm*gbna*(V-Ena);
  
  /* compute potassium background current */
  ibk = Cm*gbk*(V-Ek);
	  
  /* compute ibca background current */
  ibca = Cm*gbca*(V-Eca);
	  
  /* compute inak sodium-potassium pump current, LR-style */
  sigma = (exp(nac/67.3)-1)/7.0;
  fnak = 1/(1+0.1245*exp(-0.1*V*F/(R*T))+0.0365*sigma*exp(-V*F/(R*T)));
  inak = Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));
	  
  /* compute icap calcium pump current LR-style */
  icap = Cm*icapbar*cai/(cai+kmcap); 
  
  /* compute inaca exchanger current LR-style */
  inaca = Cm*knacalr/(pow(kmnalr,3.0)+pow(nac,3.0))/(kmcalr+cac)/
    (1+ksatlr*exp((gammalr-1)*V*F/(R*T)))*
      (nai*nai*nai*cac*exp(V*gammalr*F/(R*T))-nac*nac*nac*cai*
       exp(V*(gammalr-1)*F/(R*T)));  
  
  /* compute naidot sodium concentration derivative */
  naidot = (-3*inak-3*inaca-ibna-ina)/(F*vi);
  
  /* compute kidot potassium concentration derivative */
  kidot = (2*inak-ik1-ito-ikur-ikr-iks-ibk)/(F*vi);
  
  /* calcium buffer dynamics */
  cmdn = cmdnbar*cai/(cai+kmcmdn);
  trpn = trpnbar*cai/(cai+kmtrpn);
  csqn = csqnbar*carel/(carel+kmcsqn);
  
  /* SR calcium handling */
  irel = grelbar*uu*uu*vv*ww*(carel-cai);
  iup = iupbar/(1+kmup/cai);
  iupleak = kupleak*caup;
  itr = (caup-carel)/tautr;
  
  /* compute caidot calcium concentration derivative */
  /* using steady-state buffer approximation */
  caidot = ((2*inaca-icap-icaL-ibca)/(2*F*vi)+
	    (iupleak-iup)*vup/vi+irel*vrel/vi)/
	      (1+trpnbar*kmtrpn/(cai*cai+2*cai*kmtrpn+kmtrpn*kmtrpn)+
	       cmdnbar*kmcmdn/(cai*cai+2*cai*kmcmdn+kmcmdn*kmcmdn)); 
  
  /* update caup calcium in uptake compartment */
  caup = caup + dt*(iup-itr*vrel/vup-iupleak);
  
  /* update carel calcium in release compartment */
  carel = carel + dt*((itr-irel)/(1+csqnbar*kmcsqn/
				  (carel*carel+2*carel*kmcsqn+kmcsqn*kmcsqn)));
  
  /* update all concentrations */
  nai = nai + dt*naidot;
  ki = ki + dt*kidot;
  cai = cai + dt*caidot;
  
  /* update ina m gate */
  a = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
  if (fabs(V+47.13) < 1e-10) a = 3.2; /* denominator = 0 */
  b = 0.08*exp(-V/11);
  tau = 1/(a+b); inf = a*tau;
  m = inf + (m-inf)*exp(-dt/tau);
  
  /* update ina h gate */
  if (V >= -40.0)
    {
      a  = 0.0;
      b = 1/(0.13*(1+exp((V+10.66)/-11.1)));
    }
  else
    {
      a = 0.135*exp((V+80)/-6.8);
      b = 3.56*exp(0.079*V)+3.1e5*exp(0.35*V);
    }
  tau = 1/(a+b); inf = a*tau;
  h = inf + (h-inf)*exp(-dt/tau);
  
  /* update ina j gate */
  if (V >= -40.0)
    {
      a  = 0.0;
      b = 0.3*exp(-2.535e-7*V)/(1+exp(-0.1*(V+32)));
    }
  else
    {
      a = (-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V))*(V+37.78)/
	(1+exp(0.311*(V+79.23)));
      b = 0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)));
    }
  tau = 1/(a+b); inf = a*tau;
  j = inf + (j-inf)*exp(-dt/tau);
  
  /* update oa ito gate */
  /* corrected for 37 deg */
  /* define an voltage shift due to junctional potential
     and effect of Cd++ */
  vshift = -10;
  a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));
  b = 0.65/(2.5+exp((V-vshift+72.0)/17.0)); 
  tau = 1/(a+b); 
  inf = 1/(1+exp((V-vshift+10.47)/-17.54));
  oa = inf + (oa-inf)*exp(-Tfac*dt/tau);

  /* update oi and ois ito gate */
  /* corrected for 37 deg */
  /* define an voltage shift due to junctional potential
     and effect of Cd++ */
  vshift = -10;
  a = 1/(18.53+exp((V-vshift+103.7)/10.95));
  b = 1/(35.56+exp((V-vshift-8.74)/-7.44));
  tau = 1/(a+b); 
  inf = 1/(1+exp((V-vshift+33.1)/5.3));
  oi = inf + (oi-inf)*exp(-Tfac*dt/tau);

  /* update ua ikur gate */
  /* corrected for 37 deg */
  /* define an voltage shift due to junctional potential
     and effect of Cd++ */
  vshift = -10;
  a = 0.65/(exp((V-vshift+0.0)/-8.5)+exp((V-vshift-40.0)/-59.0));
  b = 0.65/(2.5+exp((V-vshift+72.0)/17.0)); 
  tau = 1/(a+b); 
  inf = 1/(1+exp((V-vshift+20.3)/-9.6));
  ua = inf + (ua-inf)*exp(-Tfac*dt/tau);

  /* update ui ikur gate */
  /* corrected for 37 deg */
  /* define an voltage shift due to junctional potential
     and effect of Cd++ */
  vshift = -10;
  a = 1/(exp((V-vshift-195)/-28)+21);
  b = 1/(exp((V-vshift-168)/-16));
  tau = 1/(a+b); 
  inf = 1/(1+exp((V-vshift-109.45)/27.48));
  ui = inf + (ui-inf)*exp(-Tfac*dt/tau);

  /* update the xr ikr gate */
  vshift = 0;
  a = 0.0003*(V-vshift+14.1)/(1-exp((V-vshift+14.1)/-5));
  b = 0.000073898*(V-vshift-3.3328)/(exp((V-vshift-3.3328)/5.1237)-1);
  if (fabs(V-vshift+14.1) < 1e-10) a = 0.0015; /* denominator = 0 */ 
  if (fabs(V-vshift-3.3328) < 1e-10) b = 3.7836118e-4; /* denominator = 0 */ 
  tau = 1/(a+b); 
  inf = 1/(1+exp((V-vshift+14.1)/-6.5));
  xr = inf + (xr-inf)*exp(-dt/tau);
  
  /* update the xs ikr gate */
  vshift = 0;
  a = 0.00004*(V-vshift-19.9)/(1-exp((V-vshift-19.9)/-17));
  b = 0.000035*(V-vshift-19.9)/(exp((V-vshift-19.9)/9)-1);
  if (fabs(V-vshift-19.9) < 1e-10) /* denominator = 0 */
    {
      a = 0.00068; 
      b = 0.000315; 
    }
  /* tau reduced by 50% as described in manuscript */
  tau = 0.5/(a+b); 
  inf = sqrt(1/(1+exp((V-vshift-19.9)/-12.7)));
  xs = inf + (xs-inf)*exp(-dt/tau);

  /* update icaL d gate */
  vshift = 0; 
  a = 1/(1+exp((V-vshift+10)/-6.24));
  tau = a*(1-exp((V-vshift+10)/-6.24))/(0.035*(V-vshift+10));
  if (fabs(V-vshift+10) < 1e-10) tau = a*4.579; /* denominator = 0 */
  inf = 1/(1+exp((V-vshift+10)/-8));
  d = inf + (d-inf)*exp(-dt/tau);
	  
  /* update icaL f gate */
  vshift = 0;
  inf = exp(-(V-vshift+28)/6.9)/(1+exp(-(V-vshift+28)/6.9));
  tau = 1.5*2*3/(0.0197*exp(-0.0337*0.0337*(V-vshift+10)*
			    (V-vshift+10))+0.02);
  f = inf + (f-inf)*exp(-dt/tau); 

  /* update the SR gating variables */
  /* caflux is expected in umoles/ms, hence the factor of 1000 */
  /* 1e-15 is used to scale the volumes! */ 
  vshift=0;
  caflux = 1e3*(1e-15*vrel*irel-1e-15*(0.5*icaL-0.1*2*inaca)/(2*F));
  inf = 1/(1+exp(-(caflux-3.4175e-13-vshift)/13.67e-16));
  tau = 8.0;
  uu = inf + (uu-inf)*exp(-dt/tau); 
  inf = 1-1/(1+exp(-(caflux-6.835e-14-vshift)/13.67e-16));
  tau = 1.91+2.09/(1+exp(-(caflux-3.4175e-13-vshift)/13.67e-16));
  vv = inf + (vv-inf)*exp(-dt/tau); 
  inf = 1-1/(1+exp(-(V-40)/17.0));
  tau = 6.0*(1-exp(-(V-7.9)/5.0))/(1+0.3*exp(-(V-7.9)/5.0))/(V-7.9);
  if (fabs(V-7.9) < 1e-10) tau = 6*0.2/1.3;
  ww = inf + (ww-inf)*exp(-dt/tau);
	  
	  
  /* update membrane voltage */
  V = V - dt*(ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca)/Cm;
