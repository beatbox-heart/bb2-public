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

/*
 * Human atrium model from: 
 * A. Nygren, C. Fiset, L. Firek, J.W. Clark, D.S. Lindblad, R.B. Clark and W.R. Giles, 
 * ``Mathematical Model of an Adult Human Atrial Cell. The Role of K^+ Currents in Repolarization'',
 * Circ. Res. 82:63-81 (1998)
 */

/* Versions*.
 * 0 coded myself
 * 1 after comparison with HZ's
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
			    #ifdef gamma
			    #undef gamma
			    #endif
#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"


/* To account for disrepancies in the Circ Res (CR) and Henggui's (HZ) versions */
#define CR(a,b) (a)
#define HZ(a,b) (b)
#define WHICH(a,b) (b)

static real sqr(real a) {return a*a;}
static real cub(real a) {return a*a*a;}
static real qrt(real a) {return sqr(sqr(a));}
/* H&H exponent */
static real hhx(real x){ if(-1.e-4<x && x<1.e-4) return 1+x/2.+x*x/12.; else return x/(1-exp(-x));}


					/* UNITS */
#define _(name,m) static real name=1.;
#include "ha_un.h"
#undef _				
					/* INDICES OF DYNAMIC VARIABLES */
enum {
  #define _(name,i) ha_##name,
  #include "ha_var.h"
  #undef _
  N		    /* total number of variables */
};

typedef struct {			/* PARAMETERS' STRUCTURE */
  #define _(name,v) real name;
  #include "ha_par.h"
  #undef _
} STR;


RHS_HEAD(ha,N) {
					/* DECLARATIONS */
                                        /* parameters */
  #define _(name,v) DEVICE_CONST(real,name)
  #include "ha_par.h"
  #undef _
                                        /* dynamic variables */
  #define _(name,v) real name=u[__LINE__-1];
  #include "ha_var.h"
  #undef _
					/* rates */
  #define _(name,v) real *d##name=du+__LINE__-1;
  #include "ha_var.h"
  #undef _
					/* currents and other aux vars */
  #define _(name,v) real name;
  #include "ha_cur.h"
  #undef _

					/* CALCULATIONS */

/* Currents */

  #define _(name,v) name=v; 
  #include "ha_cur.h"
  #undef _

/* Intracellular Ca buffering */
  *dOC	= 2e5 * Cai*(1-OC) - 476*OC;
  *dOTC	= 78400 * Cai*(1-OTC) - 392*OTC;
  *dOTMgC = 2e5 * Cai*(1-OTMgC-OTMgMg) - 6.6*OTMgC;
  *dOTMgMg = 2e3 * Mgi*(1-OTMgC-OTMgMg) - 666*OTMgMg;		

/* Intracellular ion concentrations */
  *dNai	= -(INa+IBNa+3*INaK+3*INaCa-PhiNaen)/(Voli*F);		/* NB: PhiNaen sign inverted */
  *dKi	= -(It+Isus+IK1+IKs+IKr+WHICH(-2*INaK,+2*INaK))/(Voli*F);
  *dCai	= -(WHICH(-Idi,ICaL) +IBCa+ICaP-2*INaCa+Iup-Irel)/(2*Voli*F) - (
    0.08*(*dOTC) + 0.16*(*dOTMgC) + 0.045*(*dOC)
  );
  *dCad	= -(ICaL+WHICH(-Idi,+Idi))/(2*Vold*F);

/* Ca2+ handling by the sarcoplasmic reticulum */
  *dOCalse = 480*Carel*(1-OCalse) - 400*OCalse;
  *dCarel = (Itr-Irel)/(2*Volrel*F) - 31*(*dOCalse);
  *dCaup =(Iup-Itr)/(2*Volup*F);
  *dF1	= rrecov*(1-F1-F2) - ract*F1;
  *dF2	= ract*F1 - rinact*F2;
  
/* Cleft space concentrations */
  *dNac	= (Nab-Nac)/tNa + (INa+IBNa+3*INaK+3*INaCa-PhiNaen)/(Volc*F);	/* NB: PhiNaen sign inverted */
  *dKc	= (Kb-Kc)/tK + (It+Isus+IK1+IKs+IKr+WHICH(-2*INaK,+2*INaK))/(Volc*F);
  *dCac = (Cab-Cac)/tCa + (ICaL+IBCa+ICaP-2*INaCa)/(2*Volc*F);

/* Gates */
  #define gate(name) u[ha_##name] = (_##name)+(name-(_##name))*exp(-ht/(t##name)); du[ha_##name]=0;
  gate(m); 
  gate(h1); 
  gate(h2); 
  gate(dL);
  gate(fL1);
  gate(fL2);
  gate(r)
  gate(s)
  gate(rsus)
  gate(ssus)
  gate(n)
  gate(pa)

/* Voltage */
  *dV = -1./Cm*( INa+ICaL+It+Isus+IK1+IBNa+IBCa+INaK+ICaP+INaCa+WHICH(0,IKr+IKs)-Iext );

} RHS_TAIL(ha )

RHS_CREATE_HEAD(ha ) {
  #define _(name,v) ACCEPTP(name,v,RNONE,RNONE);
  #include "ha_par.h"
  #undef _
  CALLOC(*u,N,sizeof(real));
  #define _(name,v) (*u)[ha_##name]=v;
  #include "ha_var.h"
  #undef _
} RHS_CREATE_TAIL(ha,1L*N)
