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
 * D.Noble's collection of heart models as in Oxsoft Heart.
 * Modification: gbK as a variable parameter.
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
extern double gbK;

typedef struct {
  real IE;
  real gbk;
  char out[MAXPATH];
} STR;

RHS_HEAD(noble_gbk,neqn) {
  DEVICE_CONST(real,IE)
  DEVICE_CONST(real,gbk)
  Arneq Alpha, Beta;
  gbK=gbk;
  memcpy (Alpha, alpha, sizeof(Arneq));
  memcpy (Beta, beta, sizeof(Arneq));
  FTN(/*0,*/u,du,Alpha,Beta);
  du[0]+=IE;
} RHS_TAIL(noble_gbk)

RHS_CREATE_HEAD(noble_gbk) {
  char temp[MAXPATH];
  int v;

  char temp4[] = "tmp/bbtp.XXXXXX";
  int errorbbtp = 0;
  
  ACCEPTP(IE,0,RNONE,RNONE);
  ACCEPTS(out,NULLFILE);
  if NOT(Filout3=fopen(S->out,"wt")) ABORT("\ncannot open %s",S->out);

/*  tmpnam(temp); */

  system("mkdir -p tmp");
  errorbbtp = 0;
  errorbbtp = mkstemp(temp4);

  if NOT(Filin1=fopen(temp4,"w+t")) ABORT("\ncannot open temporary file line 57 of noble_gbk");
  fputs(w,Filin1);
  fputc('$',Filin1);
  rewind(Filin1);
  START(); 
  FIRST();
  FTN(/*&Tstart, */_Y, _F, alpha, beta);
  /* NB!! Fast=0;*/ /* Originally was set automatically by DESOLVE/ADAMS.. */
  fclose(Filin1);
  unlink(temp4);
  fprintf(Filout3,"Initial values of dependent variables (Pascal/C):\n");
  for (v=0;v<10;v++) fprintf(Filout3,"%10d/%1d ",v+1,v);
  for (v=0;v<neqn;v++) {
    if (v%10==0) fprintf(Filout3,"\n");
    fprintf(Filout3,"%12.7lf ",_Y[v]);
  }
  fprintf(Filout3,"\n");
  fclose(Filout3);
  MALLOC(*u,sizeof(Arneq));
  memcpy (*u, _Y, sizeof(Arneq));
  ACCEPTP(gbk,gbK,RNONE,RNONE);
} RHS_CREATE_TAIL(noble_gbk,neqn)
