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

#define TI2AB(g) real alp_s##g = g##_INF/TAU_##g; real bet_s##g = (1.0-g##_INF)/TAU_##g;

#define svolt V

real rec_iNaK=(1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
real rec_ipK=1./(1.+exp((25-svolt)/5.98));

/* compute steady state values and time constants  */
real AM=1./(1.+exp((-60.-svolt)/5.));
real BM=0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
real TAU_M=AM*BM;
real M_INF=1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
TI2AB(M);

real AH;
real BH;
if (svolt>=-40.)
  {
    AH=0.; 
    BH=(0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
  }
 else
   {
     AH=(0.057*exp(-(svolt+80.)/6.8));
     BH=(2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
   }
real TAU_H=1.0/(AH+BH);
real H_INF=1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
TI2AB(H);

real AJ;
real BJ;
if(svolt>=-40.)
  {
    AJ=0.;      
    BJ=(0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
  }
 else
   {
     AJ=(((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*
	  exp(-0.04391*svolt))*(svolt+37.78)/
	 (1.+exp(0.311*(svolt+79.23))));    
     BJ=(0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
   }
real TAU_J= 1.0/(AJ+BJ);
real J_INF=H_INF;
TI2AB(J);

real Xr1_INF=1./(1.+exp((-26.-svolt)/7.));
real axr1=450./(1.+exp((-45.-svolt)/10.));
real bxr1=6./(1.+exp((svolt-(-30.))/11.5));
real TAU_Xr1=axr1*bxr1;
TI2AB(Xr1);

real Xr2_INF=1./(1.+exp((svolt-(-88.))/24.));
real axr2=3./(1.+exp((-60.-svolt)/20.));
real bxr2=1.12/(1.+exp((svolt-60.)/20.));
real TAU_Xr2=axr2*bxr2;
TI2AB(Xr2);

real Xs_INF=1./(1.+exp((-5.-svolt)/14.));
real Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
real Bxs=(1./(1.+exp((svolt-35.)/15.)));
real TAU_Xs=Axs*Bxs+80;
TI2AB(Xs);

#ifdef EPI
    real R_INF=1./(1.+exp((20-svolt)/6.));
    real S_INF=1./(1.+exp((svolt+20)/5.));
    real TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    real TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif
#ifdef ENDO
    real R_INF=1./(1.+exp((20-svolt)/6.));
    real S_INF=1./(1.+exp((svolt+28)/5.));
    real TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    real TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
#endif
#ifdef MCELL
    real R_INF=1./(1.+exp((20-svolt)/6.));
    real S_INF=1./(1.+exp((svolt+20)/5.));
    real TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    real TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif
TI2AB(R);
TI2AB(S);
      
real D_INF=1./(1.+exp((-8-svolt)/7.5));
real Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
real Bd=1.4/(1.+exp((svolt+5)/5));
real Cd=1./(1.+exp((50-svolt)/20));
real TAU_D=Ad*Bd+Cd;
TI2AB(D);

real F_INF=1./(1.+exp((svolt+20)/7));
real Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
real Bf=200./(1+exp((13-svolt)/10.));
real Cf=(180./(1+exp((svolt+30)/10)))+20;
real TAU_F=Af+Bf+Cf;
TI2AB(F);

real F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
real Af2=600*exp(-(svolt+25)*(svolt+25)/170);
real Bf2=31/(1.+exp((25-svolt)/10));
real Cf2=16/(1.+exp((svolt+30)/10));
real TAU_F2=Af2+Bf2+Cf2;
TI2AB(F2);

#undef TI2AB
#undef svolt
