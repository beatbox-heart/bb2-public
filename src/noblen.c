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
 * Multi-compartment modification for simulation external electric current
 * (n+1) compartments, each has separate set of fast gating vars
 * and characterised  by resistivity R[k] and partial surface C[k], k=0..n
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
/* 
 * Fastest variables:
 *  7 [6] --- iCa activation
 *  8 [7] --- iCa inactivation
 *  9 [8] --- iNa inactivation
 *  10 [9] -- iNa activation
 * - multiplicated to simulate different parts of the
 * cell membrane, and new variable I added to simulate electric current
 * determining potential differences between the parts.

 * Allocation of variables in Beatbox v-indices:
 * Index		Contents, in terms of HEART indices
 * 0..neqn-1		0..neqn, including fast vars of compartment #0
 * neqn			external current
 * For k-th comparment, k=1..n
 * neqn+4*(k-1)+1=neqn+4*k-3	(6) of k-th part
 * neqn+4*(k-1)+2=neqn+4*k-2	(7) of k-th part
 * neqn+4*(k-1)+3=neqn+4*k-1	(8) of k-th part
 * neqn+4*(k-1)+4=neqn+4*k	(9) of k-th part
 * 
 * Ergo: neqnn=N=neqn+4*n+1
 */

typedef int indexmap[neqn];
typedef struct {
  real IE;
  char out[MAXPATH];
  int n, neqnn;
  indexmap *h2k;    /* [n+1] */
  real *R;	    /* [n+1] */
  real *C;	    /* [n+1] */
  Arneq *Y;	    /* [n+1] */
  Arneq *F;	    /* [n+1] */
} STR;

RHS_HEAD(noblen,S->neqnn) {
  DEVICE_CONST(real,IE)
  DEVICE_CONST(int,n)
  DEVICE_ARRAY(indexmap,h2k)
  DEVICE_ARRAY(real,R)
  DEVICE_ARRAY(real,C)
  DEVICE_ARRAY(Arneq,Y)
  DEVICE_ARRAY(Arneq,F)
  static Arneq Alpha, Beta;
  int k,v;

  for (k=0;k<=n;k++) {
    memcpy (Alpha,alpha,sizeof(Arneq));
    memcpy (Beta,beta,sizeof(Arneq));
    for (v=0;v<neqn;v++) Y[k][v]=u[h2k[k][v]];
    Y[k][0]+=R[k]*u[neqn];
    FTN(Y[k],F[k],Alpha,Beta);
  }

  for(v=0;v<neqn;v++) if(v<6 || v>9) {
    u[v]=du[v]=0;
    for(k=0;k<=n;k++) {
      u [v]+=C[k]*Y[k][v];
      du[v]+=C[k]*F[k][v];
     }
  }

  for(v=6;v<=9;v++) for(k=0;k<=n;k++) {
      u [h2k[k][v]]=Y[k][v];
      du[h2k[k][v]]=F[k][v];
  }

  du[0]+=IE;
} RHS_TAIL(noblen)

RHS_CREATE_HEAD(noblen) {
  char temp[MAXPATH];
  real sum;
  int k, v;

  char temp5[] = "tmp/bbtp.XXXXXX";
  int errorbbtp = 0;

  ACCEPTP(IE,0,RNONE,RNONE);
  ACCEPTS(out,NULLFILE);
  if NOT(Filout3=fopen(S->out,"wt")) ABORT("\ncannot open %s",S->out);

/*  tmpnam(temp); */
  system("mkdir -p tmp");
  errorbbtp = 0;
  errorbbtp = mkstemp(temp5);

  if NOT(Filin1=fopen(temp5,"w+t")) ABORT("\ncannot open temporary file");
  fputs(w,Filin1);
  fputc('$',Filin1);
  rewind(Filin1);
  START(); 
  FIRST();
  FTN(/*&Tstart, */_Y, _F, alpha, beta);
  /* NB!! Fast=0;*/ /* Originally was set automatically by DESOLVE/ADAMS.. */
  fclose(Filin1);
  unlink(temp5);
  fprintf(Filout3,"Initial values of dependent variables (Pascal/C):\n");
  for (v=0;v<10;v++) fprintf(Filout3,"%10d/%1d ",v+1,v);
  for (v=0;v<neqn;v++) {
    if (v%10==0) fprintf(Filout3,"\n");
    fprintf(Filout3,"%12.7lf ",_Y[v]);
  }
  fprintf(Filout3,"\n");
  fclose(Filout3);

  ACCEPTI(n,INONE,0,INONE);
  S->neqnn=neqn+4*(S->n)+1;
  CALLOC(S->h2k,(S->n)+1,sizeof(indexmap));
  CALLOC(S->R,(S->n)+1,sizeof(real));
  CALLOC(S->C,(S->n)+1,sizeof(real));
  CALLOC(S->Y,(S->n)+1,sizeof(Arneq));
  CALLOC(S->F,(S->n)+1,sizeof(Arneq));
  for(k=0;k<=(S->n);k++) {
     for(v=0;v<neqn;v++) (S->h2k)[k][v]=v;
     for(v=6;v<=9;v++) (S->h2k)[k][v]=k?(neqn+4*k+v-9):v;
     snprintf(temp,MAXPATH,"R%d=",k);
     if (!acceptr(temp,&(S->R[k]),RNONE,RNONE,RNONE,w)) return(0);
     snprintf(temp,MAXPATH,"C%d=",k);
     if (!acceptr(temp,&(S->C[k]),RNONE,RSUCC(0.),1.0,w)) return(0);
  }
  
  for(sum=0,k=0;k<=(S->n);k++) sum+=S->C[k]; 
  assert(sum!=0);
  if (sum!=1.0) {
    MESSAGE("\n/*Sk make in total %lg not exactly 1.0 and rescaled*/\n",sum);
    for(k=0;k<=(S->n);k++) S->C[k]/=sum;
  }

  MALLOC(*u,sizeof(real)*(S->neqnn));
  for(v=0;v<neqn;v++) for(k=0;k<=(S->n);k++) (*u)[(S->h2k)[k][v]]=_Y[v];
} RHS_CREATE_TAIL(noblen,S->neqnn)
