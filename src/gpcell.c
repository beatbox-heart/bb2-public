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

/* Noble'91 guinea pig ventricle model */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

#define neqn 60
typedef real gpvector[neqn];

void gpstep(real *y, real *dy);

typedef int STR;
#define _(name,value) static real name;
#include "gpcell.h"
#undef _

static gpvector _Y;

RHS_HEAD(gpcell,neqn) {
  gpstep(u,du);
} RHS_TAIL(gpcell)


#define delta4          0.00001   /* smaller value */
#define Faraday         965.0e2   /* Faraday constant */

#define _(name,num,val) gpc_##name=num,
enum {
  #include "gpc_var.h"
  gpc_total
};
#undef _

RHS_CREATE_HEAD(gpcell) {
  gpvector _F;
  free(S);

  #undef ACCEPTP
  #define _(name,value) if (!acceptp(#name"=",&(name),value,0,RNONE,w,var,&ivar,v0)) return(0);
  #include "gpcell.h"
  #undef _

  #define _(name,num,val) _Y[num]=val;
  #include "gpc_var.h"
  #undef _

  gpstep(_Y, _F);

  MALLOC(*u,sizeof(gpvector));
  memcpy (*u, _Y, sizeof(gpvector));

} RHS_CREATE_TAIL(noble,neqn)

static real sqr(real a) {return a*a;}
/* static real cub(real a) {return a*a*a;} */

#define delta3 (1.e-4)
real hhexp(real E) {
  if (fabs(E)<delta3) return 1;
  else return (E/(1.-exp(-E)));
}

static void Yset(short n, real *y, real *dy, real alpha, real beta) {
  y[n] = alpha / (alpha + beta);
  dy[n] = 0.0;
}

static void RateSet(short n, real *y, real *dy, real alpha, real beta) {
  dy[n] = (alpha - y[n] * (alpha + beta));
}


void gpstep(real *y, real *dy) {

  /* extract variable values from the array */

  #define _(name,num,val) real name=y[num];
  #include "gpc_var.h"
  #undef _
  #define D(var) dy[gpc_##var]

  real ENa,Emh,EK,ECa;
  real IbNa,INaCa,IbCa;
  real Ito,IsiK,INaK,IbK;
  real IsiNa;
  #define IsiCa iCa
  real Itr,Iup,Irel;
  real alpha,beta,CaoNai3,CaiNao3,VFRT,eVS,eV,df,dfhh;

  /***** Gates ********/

  RateSet(gpc_x,y,dy,0.5*exp(0.0826*(V+50))/(1+exp(0.057*(V+50))),
	1.3*exp(-0.06*(V+20))/(1+exp(-0.04*(V+20))));
  RateSet(gpc_d,y,dy,360*hhexp(0.25*(V+19)),360*hhexp(-0.1*(V+19)));
  RateSet(gpc_f,y,dy,12.5*hhexp(-0.25*(V+34)),25/(1+exp(-0.25*(V+34))));
	alpha = 20.0 * exp(-0.125*(V+75));
	beta = 2000.0 / (1+320.0*exp(-0.1*(V+75)));
	alpha += beta*gB;	/* this is LQT-3 model */
	beta -= beta*gB;
  RateSet(gpc_h,y,dy,alpha,beta);
  Yset   (gpc_m,y,dy,2000.0*hhexp(0.1*(V+41)),8000.0*exp(-0.056*(V+66)));
	m=y[gpc_m];
  RateSet(gpc_r,y,dy,0.033*exp(-V/17.0),33.0/(1+exp(-(V+10)/8.0)));
  D(q) = 333 * (1 / (1 + exp(-0.2*(V+4))) - q);

  /***** Reversal potentials ********/

  if (Cai < 0.00000001)
	Cai = 0.00000001;
  if (Nai < delta4)
	Nai = delta4;
  ENa = RT_F * log(Nao / Nai);
  Emh = RT_F * log( (Nao+PNaK*Ko) / (Nai+PNaK*Ki) );
  ECa = RT_F * log(Cao / Cai) / 2;
  EK  = RT_F * log(Ko / Ki);
  VFRT = V/RT_F;

  /***** Channel currents ********/

  IK = x * iKm * (Ki - Ko * exp(-VFRT)) / 140;
  IK1 = gK1*(Ko/(KmK1+Ko))*(V-EK)/(1+exp(steepK1*(V-EK+10-ShiftK1)/RT_F));
  Ito = Gto * (V - EK) * q * r;
  IbK = gbK * (V - EK);
  INa = gNa * (V - Emh) * h*m*m*m;
  IbNa = gbNa * (V - ENa);
	eV=exp(-VFRT);
	eVS=exp(VSFRT);
	df=d*f*PCa*eVS;
	dfhh=df*PCaK*hhexp(VFRT-VSFRT);
  IsiK  = dfhh*(Ki-Ko*eV);
  IsiNa = dfhh*(Nai-Nao*eV);
  IsiCa = df*2*hhexp(2*(VFRT-VSFRT))*eVS*(Cai-Cao*eV*eV);
  IbCa = gbCa * (V - ECa);

  /***** Pump currents ********/

  INaK = Pump * (Nai/(KmNa+Nai)) * (Ko/(Km+Ko));
	CaoNai3=Cao*Nai*Nai*Nai;
	CaiNao3=Cai*Nao3;
  INaCa = kNaCa
	* (CaoNai3*exp(yNaCa*VFRT)-CaiNao3*exp(-(1-yNaCa)*VFRT))
	/ (1 + DNaCa*(CaiNao3 + CaoNai3));
  INaCa=(INaCa>iNaCam)?(iNaCam):(INaCa<-iNaCam)?(-iNaCam):(INaCa);

  /* Ca sequestration flows */

  Iup = (0.4*Cai - 0.03*KKK*Caup)/(Cai + KKK*Caup + KK_K);
  Itr = 50 * (Caup - Carel);
  Irel = sqr( fact / (fact + 0.25)) * KmCa2 * Carel;

  /* Differentials */

  D(V) = -(IK+IK1+Ito+IsiK+IbK+INaK+INa+IbNa+IsiNa+INaCa+IsiCa+IbCa+Iext)/C+IV;
	alpha = 500.0 * sqr(Cai / (Cai + KmCa*((IsiCa<-0.5)?0.1:1)));
	beta = alpha + 60.0;
  D(fact)  = (1-fprod-fact)*alpha - beta*fact;
  D(fprod) = beta*fact - fprod;
  D(Ki)    = -(IK + IK1 + IsiK + IbK + Ito + INaK/(1-nNaK))/ViF;
  D(Nai)   = -(INa + IbNa*Nao/140 + IsiNa + INaK*nNaK/(nNaK-1) + 3*INaCa)/ViF;
  D(Catrp) = 1.e5*Cai*(CTrop - Catrp) - 200*Catrp;
  D(Cacmd) = 1.e5*Cai*(MTrop - Cacmd) - 50*Cacmd;
  D(Caup)  = kup*Iup-Itr;
  D(Carel) = ktr*Itr-Irel;
  D(Cai)   = krel*Irel-Iup-(IsiCa+IbCa-2*INaCa)/(2*ViF)-D(Catrp)-D(Cacmd);

  /* Save newly computed values to the array */

  #define _(name,num,val) y[num]=name;
  #include "gpc_var.h"
  #undef _
}
