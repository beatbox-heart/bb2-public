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

/* ODE right-hand sides for CRN98, except V-controlled gates */

#define small (1e-10)
#define TI2AB(g) real alp_##g = inf_##g/tau_##g; real bet_##g = (1.0-inf_##g)/tau_##g;
// #define TI2D(g) real d_##g = (inf_##g-g)/tau_##g;
real Ek = 26.71*log(kc/ki);		/* potssium reversal potential */
real Ena = 26.71*log(nac/nai);	/* sodium reversal potential */
real Eca = 13.35*log(cac/cai);	/* calcium reversal potential */
real ina = Cm*gna*m*m*m*h*j*(V-Ena);	/* fast sodium current */
real ito = Cm*gto*oa*oa*oa*oi*(V-Ek);/* transient outward current */
real ikur = Cm*gkur_scale*gkur_rel*ua*ua*ua*ui*(V-Ek);  /* ultra-rapid potassium current */
real ikr = Cm*gkr*gkr_rel*(V-Ek)*xr;	/* rapid delayed outward rectifier K current */
real iks = Cm*gks*xs*xs*(V-Ek);	/* slow delayed outward rectifier K current */
real icaL = Cm*gcaL*d*f*fca*(V-ErL);	/* calcium current */
real inf_fca = 1.0/(1.0+(cai/0.00035)); 
real tau_fca = 2.0; 
TI2AB(fca);
real ik1 = Cm*gk1*(V-Ek)*gk1_rel;	/* time independent potassium current */
real ibna = Cm*gbna*(V-Ena);		/* sodium background current */
real ibk = Cm*gbk*(V-Ek);		/* potassium background current */
real ibca = Cm*gbca*(V-Eca);		/* calcium background current */
real sigma = (exp(nac/67.3)-1)/7.0;  /* sodium-potassium pump current, LR-style */
real fnak = 1.0/(1.0+fnak1+fnak2*sigma);
real inak = Cm*inakbar*fnak*kc/(kc+kmko)/(1+pow(kmnai/nai,1.5));
real icap = Cm*icapbar*cai/(cai+kmcap); /* calcium pump current, LR-style */
real inaca = 
  Cm*knacalr*(
    cub(nai)*cac*expVgammalrF_RT 
    -
    cub(nac)*cai*expVgammalr1F_RT
  )/(
     (cub(kmnalr)+cub(nac))
     *(kmcalr+cac)
     *(1+ksatlr*expVgammalr1F_RT)
  );				/* Na-Ca exchanger current LR-style */ 
real d_nai = (-3*inak-3*inaca-ibna-ina)/(F*vi); /* sodium concentration rate */
real d_ki = (2*inak-ik1-ito-ikur-ikr-iks-ibk)/(F*vi); /* potassium concentration rate */
real cmdn = cmdnbar*cai/(cai+kmcmdn);  /* calcium buffer: calmodulin */
real trpn = trpnbar*cai/(cai+kmtrpn);  /* calcium buffer: troponin */
real csqn = csqnbar*carel/(carel+kmcsqn); /* calcium buffer: calsequestrin */
real irel = grelbar*uu*uu*vv*ww*(carel-cai); /* SR calcium handling: relocation */
real iup = iupbar/(1+kmup/cai); /* SR calcium handling: uptake */
real iupleak = kupleak*caup; /* SR calcium handling: leak */
real itr = (caup-carel)/tautr; /* SR calcium handling: translocation */
real d_cai =
  (
    (2*inaca-icap-icaL-ibca)/(2*F*vi)
    + (iupleak-iup)*vup/vi
    + irel*vrel/vi
  )/(
    1
    + trpnbar*kmtrpn/(cai*cai+2*cai*kmtrpn+kmtrpn*kmtrpn)
    + cmdnbar*kmcmdn/(cai*cai+2*cai*kmcmdn+kmcmdn*kmcmdn)
  ); /* calcium concentration rate in steady-state buffer approximation */
real d_caup = (iup-itr*vrel/vup-iupleak); /* rate of calcium in uptake compartment */
real d_carel = (itr-irel) /* rate of carel calcium in release compartment */
  /(1+csqnbar*kmcsqn/(carel*carel+2*carel*kmcsqn+kmcsqn*kmcsqn));
/* SR gating variables */
/* caflux is expected in umoles/ms, hence the factor of 1000 */
/* 1e-15 is used to scale the volumes! */ 
real caflux = 1e3*(1e-15*vrel*irel-1e-15*(0.5*icaL-0.1*2*inaca)/(2*F));
real inf_uu = 1/(1+exp(-(caflux-3.4175e-13)/13.67e-16));
real tau_uu = 8.0;
TI2AB(uu);
real inf_vv = 1-1/(1+exp(-(caflux-6.835e-14)/13.67e-16));
real tau_vv = 1.91+2.09/(1+exp(-(caflux-3.4175e-13)/13.67e-16));
TI2AB(vv);
/* membrane voltage rate */
real d_V = - (ina+icaL+icap+ik1+ito+ikur+ikr+iks+ibna+ibk+ibca+inak+inaca)/Cm;
#undef TI2AB
#undef small
