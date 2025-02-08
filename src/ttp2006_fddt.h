/**
 * Copyright (C) (2010-2025) Vadim Biktashev, Irina Biktasheva et al. 
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

#define TI2AB(g) real alp_##g = g##_INF/TAU_##g; real bet_##g = (1.0-g##_INF)/TAU_##g;

#define sm sM
#define sh sH
#define sj sJ
#define sxr1 sXr1
#define sxr2 sXr2
#define sxs sXs
#define sr sR
#define ss sS
#define sd sD
#define sf sF
#define sf2 sF2
#define sfcass FCaSS

#define Istim 0
#define dsvolt dV
#define svolt V

/* These could be calculated once at startup */
/* if the parameters can be assumed constant */
double inverseVcF2=1/(2*Vc*F);
double inverseVcF=1./(Vc*F);
double inversevssF2=1/(2*Vss*F);

/* update free calcium concentrations */
/* CaSR */
real bjsr=Bufsr-CaSRtot+Kbufsr;
real CaSR=(sqrt(bjsr*bjsr+4*Kbufsr*CaSRtot)-bjsr)/2;
/* CaSS */
real bcss=Bufss-CaSStot+Kbufss;
real CaSS=(sqrt(bcss*bcss+4*Kbufss*CaSStot)-bcss)/2;
/* Cai */
real bc=Bufc-Caitot+Kbufc;
real Cai=(sqrt(bc*bc+4*Kbufc*Caitot)-bc)/2;
    
/* Needed to compute currents */
real Ek=RTONF*(log((Ko/Ki)));
real Ena=RTONF*(log((Nao/Nai)));
real Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
real Eca=0.5*RTONF*(log((Cao/Cai)));
real Ak1=0.1/(1.+exp(0.06*(svolt-Ek-200)));
real Bk1=(3.*exp(0.0002*(svolt-Ek+100))+
	  exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
real rec_iK1=Ak1/(Ak1+Bk1);

/* Compute currents */
real INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
real ICaL=GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
  (0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.);
real Ito=Gto*sr*ss*(svolt-Ek);
real IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
real IKs=Gks*sxs*sxs*(svolt-Eks);
real IK1=GK1*rec_iK1*(svolt-Ek);
real INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
  (1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
  (exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
   exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
real INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
real IpCa=GpCa*Cai/(KpCa+Cai);
real IpK=GpK*rec_ipK*(svolt-Ek);
real IbNa=GbNa*(svolt-Ena);
real IbCa=GbCa*(svolt-Eca);

/* Determine total current */
real (sItot) = IKr    +
  IKs   +
  IK1   +
  Ito   +
  INa   +
  IbNa  +
  ICaL  +
  IbCa  +
  INaK  +
  INaCa +
  IpCa  +
  IpK   +
  Istim;

/* RyR */
real kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR)));
real k1=k1_/kCaSR;
real k2=k2_*kCaSR;
real alp_sRR=k4;
real bet_sRR=k2*CaSS;
real sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);

/* Ca currents in SR */
real Irel=Vrel*sOO*(CaSR-CaSS);
real Ileak=Vleak*(CaSR-Cai);
real Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
real Ixfer=Vxfer*(CaSS-Cai);

/* update concentrations     */
/* calcium in SR */
real dCaSRtot=(Iup-Irel-Ileak);
/* calcium in sub-space */
real dCaSStot=(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
/* bulk calcium concentration */
real dCaitot=((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
/* sodium concentration in Intracellular space */
real dNai=-(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
/* potassium concentration in Intracellular space */    
real dKi=-(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;

/* non tabulated gates definition */
real FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
real TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;
TI2AB(FCaSS);

/* update voltage */
real dsvolt= (-sItot);

#undef TI2AB

#undef sm          
#undef sh
#undef sj
#undef sxr1
#undef sxr2
#undef sxs
#undef sr
#undef ss
#undef sd
#undef sf
#undef sf2
#undef sfcass
#undef Istim
#undef dsvolt
#undef svolt

