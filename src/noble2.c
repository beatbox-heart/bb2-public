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

/* D.Noble's collection of heart models as in Oxsoft Heart
 * Two-compartment modification for simulation external electric current
 */
#include <assert.h>
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

#include "noble.h"
#define OWN
#include "noble2.h"

typedef struct {
  real IE;
  char out[MAXPATH];
} STR;

RHS_HEAD(noble2,neqn2) {
  DEVICE_CONST(real,IE)
  static Arneq Yp, Fp, Alphap, Betap;
  static Arneq Yn, Fn, Alphan, Betan;
  int v;
  
  for (v=0;v<neqn;v++) {Yp[v]=u[h2p[v]];Alphap[v]=alpha[v];Betap[v]=beta[v];}
  Yp[0]+=0.5*u[1];
  FTN(Yp,Fp,Alphap,Betap);

  for (v=0;v<neqn;v++) {Yn[v]=u[h2n[v]];Alphan[v]=alpha[v];Betan[v]=beta[v];}
  Yn[0]-=0.5*u[1];
  FTN(Yn,Fn,Alphan,Betan);
  
  for (v=0;v<6;v++) {
    u [h2p[v]] = 0.5*(Yp[v] + Yn[v]);
    du[h2p[v]] = 0.5*(Fp[v] + Fn[v]);
  }
  for (v=6;v<10;v++) {
    u[ h2p[v]] = Yp[v]; u [h2n[v]] = Yn[v];  
    du[h2p[v]] = Fp[v]; du[h2n[v]] = Fn[v];  
  }
  for (v=10;v<neqn;v++) {
    u [h2p[v]] = 0.5*(Yp[v] + Yn[v]);
    du[h2p[v]] = 0.5*(Fp[v] + Fn[v]);
  }
  
  du[0]+=IE;
} RHS_TAIL(noble2)

RHS_CREATE_HEAD(noble2) {
  char temp[MAXPATH];

  int v;
   char temp3[] = "tmp/bbtp.XXXXXX";
  int errorbbtp = 0;
  
  ACCEPTP(IE,0,RNONE,RNONE);
  ACCEPTS(out,NULLFILE);
  if NOT(Filout3=fopen(S->out,"wt")) ABORT("\ncannot open %s",S->out);

/*  tmpnam(temp); */

  system("mkdir -p tmp");
  errorbbtp = 0;
  errorbbtp = mkstemp(temp3);

  if NOT(Filin1=fopen(temp3,"w+t")) ABORT("\ncannot open temporary file line 75 of noble2.c");
  fputs(w,Filin1);
  fputc('$',Filin1);
  rewind(Filin1);
  START(); 
  FIRST();
  FTN(/*&Tstart, */_Y, _F, alpha, beta);
  /* NB!! Fast=0;*/ /* Originally was set automatically by DESOLVE/ADAMS.. */
  fclose(Filin1);
  unlink(temp3);
  fprintf(Filout3,"Initial values of dependent variables (Pascal/C):\n");
  for (v=0;v<10;v++) fprintf(Filout3,"%10d/%1d ",v+1,v);
  for (v=0;v<neqn;v++) {
    if (v%10==0) fprintf(Filout3,"\n");
    fprintf(Filout3,"%12.7lf ",_Y[v]);
  }
  fprintf(Filout3,"\n");
  fclose(Filout3);
  MALLOC(*u,sizeof(real)*neqn2);
  for(v=0;v<neqn;v++) (*u)[h2n[v]]=(*u)[h2p[v]]=_Y[v];
} RHS_CREATE_TAIL(noble2,neqn2)
