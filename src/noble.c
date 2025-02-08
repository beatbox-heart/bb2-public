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

/* D.Noble's collection of heart models as in Oxsoft Heart */
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

typedef struct {
  real IE;
  char out[MAXPATH];
} STR;

RHS_HEAD(noble,neqn) {
  DEVICE_CONST(real,IE)
  Arneq Alpha, Beta;
  memcpy (Alpha, alpha, sizeof(Arneq));
  memcpy (Beta, beta, sizeof(Arneq));
  FTN(/*0,*/u,du,Alpha,Beta);
  du[0]+=IE;
} RHS_TAIL(noble)

RHS_CREATE_HEAD(noble) {
  char temp[MAXPATH];
/* there should be a /tmp here, but I will see these files being created before accessing the /tmp */
  char temp2[] = "tmp/bbtp.XXXXXX";
  int v;
  int errorbbtp = 0;
  
  ACCEPTP(IE,0,RNONE,RNONE);
  ACCEPTS(out,NULLFILE);
  if NOT(Filout3=fopen(S->out,"wt")) ABORT("\ncannot open %s",S->out);

/*  tmpnam(temp); */
  system("mkdir -p tmp");
  errorbbtp = 0;
  errorbbtp = mkstemp(temp2);
  if NOT(Filin1=fopen(temp2,"w+t")) ABORT("\ncannot open temporary file");
  if(errorbbtp==-1) ABORT("\ncannot open temporary file. line 48 in noble.c");
  fputs(w,Filin1);
  fputc('$',Filin1);
  rewind(Filin1);
  START(); 
  FIRST();
  FTN(/*&Tstart, */_Y, _F, alpha, beta);
  /* NB!! Fast=0;*/ /* Originally was set automatically by DESOLVE/ADAMS.. */
  fclose(Filin1);
  unlink(temp2);
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
} RHS_CREATE_TAIL(noble,neqn)
