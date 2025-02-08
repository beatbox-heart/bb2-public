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


/* Na+ current I_{Na} */
_(ENa	, RTF*log(Nac/Nai))
_(_m	, 1./(1+exp(-(V+27.12)/8.21)))
_(tm	, 0.000042*exp(-sqr((V+25.57)/28.8)+0.000024))
_(_h1	, 1./(1+exp((V+63.6)/5.3)))
_(_h2	, _h1)
_(th1	, 0.03/(1+exp(V+35.1)/3.2)+0.0003)
_(th2	, 0.12/(1+exp(V+35.1)/3.2)+0.003)
_(INa   , PNa*cub(m)*(0.9*h1+0.1*h2)*Nac*F*(exp((V-ENa)*FRT)-1)*hhx(-V*FRT))

/* Ca2+ current I_{Ca,L} */
_(_dL	, 1./(1+exp(-(V+9.0)/5.8)))
_(tdL	, 0.0027*exp(-sqr((V+35.0)/30.0))+0.002)
_(_fL1	, 1./(1+exp((V+27.4)/7.1)))
_(_fL2	, _fL1)
_(tfL1	, 0.161*exp(-sqr((V+40.0)/14.4))+0.01)
_(tfL2	, WHICH(1.3323,0.3323)*exp(-sqr((V+40.0)/14.4))+0.0626)
_(fCa	, Cad/(Cad+kCa))
_(ICaL	, _gCaL*dL*(fCa*fL1+(1-fCa)*fL2)*(V-ECaapp))

/* Transient outward K+ current I_t */
_(EK	, RTF*log(Kc/Ki))
_(_r	, 1./(1+exp(-(V-1.0)/11.0)))
_(tr	, 0.0035*exp(-sqr(V/30.0))+0.0015)
_(_s	, 1./(1+exp((V+40.5)/11.5)))
_(ts	, 0.4812*exp(-sqr((V+52.45)/14.97))+0.01414)
_(It	, _gt*r*s*(V-EK))

/* Sustained outward K+ current I_{sus} */
_(_rsus	, 1./(1+exp(-(V+4.3)/8.0)))
_(trsus	, 0.009/(1+exp((V+5.0)/12.0))+0.0005)
_(_ssus	, 0.4/(1+exp((V+20.0)/10.0))+0.6)
_(tssus	, 0.047/(1+exp((V+60)/10))+0.3)
_(Isus	, _gsus*rsus*ssus*(V-EK))

/* Delayed rectifier K+ currents I_{K,s} and I_{K,r} */
_(_n	, 1./(1+exp(-(V-19.9)/12.7)))
_(tn	, 0.7+0.4*exp(-sqr((V-20)/20)))
_(IKs	, _gKs*n*(V-EK))
_(_pa	, 1./(1+exp(-(V+15)/6)))
_(tpa	, 0.03118+0.21718*exp(-sqr((V+20.1376)/22.1996)))
_(pi	, 1./(1+exp((V+55)/24)))
_(IKr	, _gKr*pa*pi*(V-EK))

/* Inward rectifier K+ current: I_{K1} */
_(IK1	, _gK1*pow(Kc,0.4457)*(V-EK)/(1+exp(1.5*(V-EK+3.6)*FRT)))

/* Background inward currents I_{B,Na} and I_{B,Ca} */
_(ECa	, 0.5*RTF*log(Cac/Cai))
_(IBNa	, _gBNa*(V-ENa))
_(IBCa	, _gBCa*(V-ECa))

/* Pump and exchanger currents: I_{NaK}, I_{CaP}, I_{NaCa} */
_(INaK	, _INaK/(1+kNaKK/Kc)/(1+pow(kNaKNa/Nai,1.5))*(V+150)/(V+200))
_(ICaP	, _ICaP/(1+kCaP/Cai))
_(INaCa	, kNaCa*(cub(Nai)*Cac*exp(gam*V*FRT)-cub(Nac)*Cai*exp((gam-1)*V*FRT))/(1+dNaCa*(cub(Nac)*Cai+cub(Nai)*Cac)))

/* ??? description not found in the paper */
_(Idi	, (Cad-Cai)*2*F*Vold/tdi)

/* Ca2+ handling by the sarcoplasmic reticulum */
_(Iup 	, _Iup*(Cai/kcyca - kxcs*kxcs*Caup/ksrca)/((Cai+kcyca)/kcyca + kxcs*(Caup+ksrca)/ksrca))
_(Itr	, (Caup-Carel)*2*F*Volrel/ttr)
_(Irel	, alprel*(Carel-Cai)/sqr(1+0.25/F2))


_(ract	, 203.8*(1/qrt(1+kreli/Cai) + 1/qrt(1+kreld*Cad)))
_(rinact, 33.96 + 339.6/qrt(1+kreli/Cai))
