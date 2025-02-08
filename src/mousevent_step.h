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
All variables.
*/
/*
double vm,cai,cass,cajsr,cansr,ltrpnca,htrpnca,o,c1,c2,c3,c4,i1,i2,i3,pc1,pc2,po1,po2;
double pryr,cna3,cna2,cna1,ona,ifna,i1na,i2na,icna2,icna3;
double nai,ki;
double atoff,itof,nks,atos,itos,aur,iur,akss,ikss,cko,ck1,ck2,ok,ik;

The equations must be the same as printed in the paper, i.e. although it may be a dynamical
variable, it must be calculated by an algebriac formula (conservation of probability) etc.
*/

/*
currents
*/

double ipca = 0.0, ibna = 0.0, oclca = 0.0, iclca = 0.0, ik1 = 0.0, sigma = 0.0, fnak = 0.0, inak = 0.0;
double alphaa = 0.0, betaa = 0.0, alphai = 0.0, betai = 0.0, iktof = 0.0;
double ass = 0.0, iss = 0.0, tautas = 0.0, tautis = 0.0, iktos = 0.0; 
double alphan = 0.0, betan = 0.0, iks = 0.0 /*, nks = 0.0*/;
double icab = 0.0;
double inaca = 0.0;
double currentikss = 0.0, taukss = 0.0;
double ena = 0.0, ek = 0.0, ecan  = 0.0;
double tauaur = 0.0, tauiur = 0.0, ikur = 0.0;
double alphaa0  = 0.0, betaa0  = 0.0, alphaa1  = 0.0, betaa1  = 0.0, ikr  = 0.0, alphai_ikr = 0.0, betai_ikr = 0.0;
double alpha = 0.0, beta = 0.0, gammagamma = 0.0, kpcf = 0.0, ical = 0.0;

double alphana11 = 0.0,  alphana12 = 0.0,  alphana13 = 0.0, betana11 = 0.0,  betana12 = 0.0;
double betana13 = 0.0,  alphana3 = 0.0,  betana3 = 0.0,  alphana2 = 0.0,  betana2 =0.0;
double alphana4 = 0.0,  betana4 = 0.0,  alphana5 = 0.0,  betana5 = 0.0,ina =0.0;
double istim=0.0;

double jleak = 0.0,jxfer = 0.0,jup = 0.0,jtrpn = 0.0,jtr = 0.0,jrel = 0.0;
double bi = 0.0, bss = 0.0, bjsr = 0.0;

double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, temp5 = 0.0, temp6 = 0.0, temp7 = 0.0, temp8 = 0.0, temp9 = 0.0, temp10 = 0.0;

double soicr_ss;

/*
Calcium concentration
*/

bi = pow(1.0+cmdntot*kmcmdn/((kmcmdn+cai)*(kmcmdn+cai)),-1);
bss = pow(1.0+cmdntot*kmcmdn/((kmcmdn+cass)*(kmcmdn+cass)),-1);
bjsr = pow(1.0+csqntot*kmcsqn/((kmcsqn+cajsr)*(kmcsqn+cajsr)),-1);

/*
Calcium fluxes
*/

// jrel = nu1*(po1+po2)*(cajsr-cass)*pryr; // cajsr = carel
// Sung et al. 2011
soicr_ss = soicr_max/(1.0 + exp((soicr_thresh - cajsr)/soicr_slope));

 jrel = nu1*(po1+po2)*(cajsr-cass)*(pryr+soicr_switch*soicr); // cajsr = carel
 jtr = (cansr - cajsr)/tautr;
 jxfer = (cass - cai)/tauxfer;
 jleak = nu2*(cansr-cai); // cansr = caup
 jup = nu3*cai*cai/(kmup*kmup+cai*cai);
 jtrpn = kplushtrpn*cai*(htrpntot-htrpnca)-kminushtrpn*htrpnca+kplusltrpn*cai*(ltrpntot-ltrpnca)-kminusltrpn*ltrpnca;

/*
RyR receptors
*/

pc1 = 1.0-(pc2+po1+po2);

/*
ICaL
*/

c1 = 1.0-(o+c2+c3+c4+i1+i2+i3);

temp1 = 1.0/(1.0+0.12*exp((vm+12.0)/10.0));
temp2 = 0.4*exp((vm+12.0)/10.0);
temp3 = 1.0+0.7*exp(-(vm+40.0)*(vm+40.0)/10.0); // it should be a 100, not 10.
temp4 = 0.75*exp(-(vm+20.0)*(vm+20.0)/400.0);

alpha = temp2*(temp3 - temp4)*temp1;
beta = 0.05*exp(-(vm+12.0)/13.0);
gammagamma = kpcmax*cass/(kpchalf+cass);
kpcf = 13.0*(1.0-exp(-(vm+14.5)*(vm+14.5)/100.0));
ical = gcal*o*(vm-ecal);

/*
INa
*/

 alphana11 = 3.802/(0.1027*exp(-(vm+2.5)/17.0)+0.2*exp(-(vm+2.5)/150.0));
 alphana12 = 3.802/(0.1027*exp(-(vm+2.5)/15.0)+0.23*exp(-(vm+2.5)/150.0));
 alphana13 = 3.802/(0.1027*exp(-(vm+2.5)/12.0)+0.25*exp(-(vm+2.5)/150.0));
 betana11 = 0.1917*exp(-(vm+2.5)/20.3);
 betana12 = 0.2*exp(-(vm-2.5)/20.3);
 betana13 = 0.22*exp(-(vm-7.5)/20.3);
 alphana3 = exp(-(vm+7.0)/7.7)*7.0e-7;
 betana3 = 0.0084+0.00002*(vm+7.0);
 alphana2 = 1.0/(0.393956+0.188495*exp(-(vm+7.0)/16.6));
 betana2 = alphana13*alphana2*alphana3/(betana13*betana3);
 alphana4 = alphana2/1000.0;
 betana4 = alphana3;
 alphana5 = alphana2/95000.0;
 betana5 = alphana3/50.0;
 cna3 = 1.0-(ona+cna1+cna2+ifna+i1na+i2na+icna2+icna3);
 ena = (R*temp/F)*log((0.9*nao+0.1*ko)/(0.9*nai+0.1*ki));
 ina = gna*ona*(vm-ena);

/*
Background calcium current
*/

ecan = (R*temp/(2*F))*log(cao/cai);
icab = gcab*(vm - ecan);

/*
Fast transient outward current
*/

alphaa = 0.18064*exp(0.03577*(vm+30.0));
betaa = 0.3956*exp(-0.06237*(vm+30.0));
alphai = 0.000152*exp(-(vm+13.5)/7.0)/(1.0+0.067083*exp((vm+33.5)/7.0));
betai = 0.00095*exp((vm+33.5)/7.0)/(1.0+0.051335*exp((vm+33.5)/7.0));
ek = (R*temp/F)*log(ko/ki);
iktof = gktof*atoff*atoff*atoff*itof*(vm-ek);

/*
Slow transient outward current
*/

ass = 1.0/(1.0+exp(-(vm+22.5)/7.7));
iss = 1.0/(1.0+exp((vm+45.2)/5.7));
tautas = 2.058+0.493*exp(-0.0629*vm);
tautis = 270.0 + 1050.0/(1.0+exp((vm+45.2)/5.7));
iktos = gktos*atos*itos*(vm-ek);

/*
IK1
*/

ik1 = 0.2938*ko*(vm - ek)/((ko + 210.0)*(1.0 + exp(0.0896*(vm  - ek))));

/*
slow delayed rectifier current IKs
*/

alphan = 0.00000481333*(vm+26.5)/(1.0 - exp(-0.128*(vm+26.5))); 
betan = 0.0000953333*exp(-0.038*(vm+26.5));
iks = gks*nks*nks*(vm-ek);

/*
Ultrarapidly activating delayed rectifier K+ current
*/

// tauaur = 2.058 + 0.493*exp(-0.0629*vm);
tauaur = tautas; // the same as for Slow transient outward current.
tauiur = 1200.0 - 170.0/(1.0+exp((vm+45.2)/5.7));
ikur = gkur*aur*iur*(vm-ek);

// by now, ass and iss are also already computed.

/*
Noninactivating steady-state K+ current
*/

taukss = 39.3*exp(-0.0862*vm)+13.17;
currentikss = gkss*akss*ikss*(vm-ek);

/*
Calcium pump
*/

ipca = ipcamax*cai*cai/(kmpca*kmpca+cai*cai);

/*
Sodium calcium exchanger
*/

temp1 = nai*nai*nai*cao*exp(eta*vm*F/(R*temp));
temp2 = nao*nao*nao*cai*exp((eta-1)*vm*F/(R*temp));
temp3 = 1.0/(kmna*kmna*kmna+nao*nao*nao);
temp4 = 1.0/(kmca+cao);
temp5 = 1.0/(1.0+ksat*exp((eta-1)*vm*F/(R*temp)));

inaca = knaca*(temp1 - temp2)*temp3*temp4*temp5;

/*
Background sodium current
*/

ibna = gnab*(vm - ena);

/*
Rapid delayed rectifier K+ current (mERG)
*/

alphaa0 = 0.022348*exp(0.01176*vm);
betaa0 = 0.047002*exp(-0.0631*vm);
alphaa1 = 0.013733*exp(0.038198*vm);
betaa1 = 0.0000689*exp(-0.04178*vm);
alphai_ikr = 0.090821*exp(0.023391*(vm+5.0));
betai_ikr = 0.006497*exp(-0.03268*(vm+5.0));
cko = 1.0-(ck1+ck2+ok+ik);
ikr = ok*gkr*(vm-(R*temp/F)*log((0.98*ko+0.02*nao)/(0.98*ki+0.02*nai)));

/*
sodium-potassium pump
*/

sigma = (exp(nao/67300.0)-1.0)/7.0;
fnak = 1.0+0.1245*exp(-0.1*vm*F/(R*temp))+0.0365*sigma*exp(-vm*F/(R*temp));
fnak = 1.0/fnak;
inak = imaxnak*fnak*ko/((ko+kmko)*(1.0+pow(kmnai/nai,3./2.)));

/*
Calcium activated chloride current
*/

oclca = 0.2/(1.0 + exp(-(vm - 46.7)/7.8));

iclca = gclca*oclca*cai*(vm - ecl)/(cai + kmcl);

/*
solve differential equations using first order Euler.
*/

/*
Calcium concentration
*/

 cai = cai + dt*(bi*(jleak+jxfer-jup-jtrpn-(icab-2*inaca+ipca)*acap*cm/(2*vmyo*F)));
 cass = cass + dt*(bss*(jrel*vjsr/vss-jxfer*vmyo/vss-ical*acap*cm/(2*vss*F)));
 cajsr = cajsr + dt*(bjsr*(jtr-jrel));
 cansr = cansr + dt*((jup-jleak)*vmyo/vnsr-jtr*vjsr/vnsr);

/*
Calcium fluxes
*/

// the PRyR takes care of the ical based activation of RyR
pryr = pryr + dt*(-0.04*pryr-0.1*ical/icalmax*exp(-(vm-5.0)*(vm-5.0)/648.0));

/*
Calcium buffering
*/

ltrpnca = ltrpnca + dt*(kplusltrpn*cai*(ltrpntot-ltrpnca)-kminusltrpn*ltrpnca);
htrpnca  = htrpnca + dt*(kplushtrpn*cai*(htrpntot-htrpnca)-kminushtrpn*htrpnca);

/*
RyR receptors
*/

 po1 = po1 + dt*(kplusa*pow(cass,n)*pc1-kminusa*po1-kplusb*pow(cass,m)*po1+kminusb*po2 - kplusc*po1+kminusc*pc2);
 po2 = po2 + dt*(kplusb*pow(cass,m)*po1-kminusb*po2);
 pc2 = pc2 + dt*(kplusc*po1-kminusc*pc2);
 pc1 = 1.0-(pc2+po1+po2);

// Sung et al. 2011
soicr = soicr_ss - (soicr_ss - soicr)*exp(-dt/soicr_tau);

/*
ICaL
*/

o =  o + dt*(alpha*c4-4*beta*o+kpcb*i1-gammagamma*o+0.001*(alpha*i2-kpcf*o));
c2 = c2 + dt*(4*alpha*c1-beta*c2+2*beta*c3-3*alpha*c2);
c3 = c3 + dt*(3*alpha*c2-2*beta*c3+3*beta*c4-2*alpha*c3);
c4 = c4 + dt*(2*alpha*c3-3*beta*c4+4*beta*o-alpha*c4+0.01*(4*kpcb*beta*i1-alpha*gammagamma*c4)+0.002*(4*beta*i2-kpcf*c4)+4*beta*kpcb*i3-gammagamma*kpcf*c4);
i1 = i1 + dt*(gammagamma*o-kpcb*i1+0.001*(alpha*i3-kpcf*i1)+0.01*(alpha*gammagamma*c4-4.0*beta*kpcf*i1));
i2 = i2 + dt*(0.001*(kpcf*o-alpha*i2)+kpcb*i3-gammagamma*i2+0.002*(kpcf*c4-4*beta*i2));
i3 = i3 + dt*(0.001*(kpcf*i1-alpha*i3)+gammagamma*i2-kpcb*i3+gammagamma*kpcf*c4-4*beta*kpcb*i3);


/*sodium dynamics*/
 nai =nai - dt*(ina + ibna + 3.0*inaca + 3.0*inak)*acap*cm/(vmyo*F);

/* potassium dynamics*/

ki= ki - dt*(iktof+iktos+ik1+iks+currentikss+ikur+ikr-2*inak+istim)*acap*cm/(vmyo*F);

/*
INa ODEs
*/

cna2  = cna2 + dt*(alphana11*cna3-betana11*cna2+betana12*cna1-alphana12*cna2+alphana3*icna2-betana3*cna2);
cna1  = cna1 + dt*(alphana12*cna2-betana12*cna1+betana13*ona-alphana13*cna1+alphana3*ifna-betana3*cna1);
ona   = ona + dt*(alphana13*cna1-betana13*ona+betana2*ifna-alphana2*ona);
ifna  = ifna + dt*(alphana2*ona-betana2*ifna+betana3*cna1-alphana3*ifna+betana4*i1na-alphana4*ifna+alphana12*icna2-betana12*ifna);
i1na  = i1na + dt*(alphana4*ifna-betana4*i1na+betana5*i2na-alphana5*i1na);
i2na  = i2na + dt*(alphana5*i1na-betana5*i2na);
icna2 = icna2 + dt*(alphana11*icna3-betana11*icna2+betana12*ifna-alphana12*icna2+betana3*cna2-alphana3*icna2);
icna3 = icna3 + dt*(betana11*icna2-alphana11*icna3+betana3*cna3-alphana3*icna3);


/*
Fast transient outward current
*/


atoff = atoff + dt*(alphaa*(1-atoff)-betaa*atoff);
itof  = itof  + dt*(alphai*(1-itof)-betai*itof);

/*
Slow transient outward current
*/

atos = atos + dt*(ass-atos)/tautas;
itos = itos + dt*(iss-itos)/tautis;

/*
slow delayed rectifier current IKur
*/

nks = nks + dt*(alphan*(1.0-nks)-betan*nks);

/*
Noninactivating steady-state K+ current
*/

akss = akss + dt*(ass-akss)/taukss;
ikss = ikss;

/*
Ultrarapidly activating delayed rectifier K+ current
*/

aur = aur + dt*(ass-aur)/tauaur;
iur = iur + dt*(iss-iur)/tauiur;

/*
IKr
*/

 ck1 = ck1 + dt*(alphaa0*cko-betaa0*ck1+kb*ck2-kf*ck1);
 ck2 = ck2 + dt*(kf*ck1-kb*ck2+betaa1*ok-alphaa1*ck2);
 ok = ok + dt*(alphaa1*ck2-betaa1*ok+betai_ikr*ik-alphai_ikr*ok);
 ik = ik + dt*(alphai_ikr*ok-betai_ikr*ik);

/* the currents on RHS are ALREADY normalised to capacitance thru' normalised maximum conductances. */

vm = vm - dt*(ical+ipca+inaca+icab+ina+ibna+inak+iktof+iktos+ik1+iks+ikur+currentikss+ikr+iclca+istim);


