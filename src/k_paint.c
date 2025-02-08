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

/* Show 2D image according to given k codes */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "k_.h"

#if defined(NOX11)
  #include "nograph.h"
  NOGRAPH_DUMMY(k_paint)
#else

#include "windraw.h"

extern int Verbose;            /* defined in main */

typedef struct {
  BGIWindow wnd;
  #include "k_code.h"
  INT nabs;
  REAL absmin, absmax;
  INT nord;
  REAL ordmin, ordmax;
  /* INT appmin, appmax; */
  double abs, ord, /* app, */ col;
  colortype *screen;
  FILE *file;
  char filename[MAXPATH];
  int ver, append;
} STR;

RUN_HEAD(k_paint)
{
  DEVICE_CONST(BGIWindow,wnd);
  #include "k_def.h"
  DEVICE_CONST(REAL,absmin) DEVICE_CONST(REAL,absmax) DEVICE_VAR(double,abs) DEVICE_CONST(INT,nabs);
  DEVICE_CONST(REAL,ordmin) DEVICE_CONST(REAL,ordmax) DEVICE_VAR(double,ord) DEVICE_CONST(INT,nord); 
  /* DEVICE_CONST(INT,appmin) DEVICE_CONST(INT,appmax) DEVICE_VAR(double,app); */
  DEVICE_VAR(REAL,col)
  DEVICE_ARRAY(colortype,screen)
  DEVICE_CONST(FILE *,file) DEVICE_ARRAY(char,filename) DEVICE_CONST(int,ver) DEVICE_CONST(int,append)
  int c, i, /* iapp, */ iord, iabs;

  SetWindow(wnd);
  if NOT(SetLimits(0,nabs,0,nord)) return 0;
  k_on();
  for(iabs=0;iabs<nabs;iabs++) { *abs=absmin+iabs*(absmax-absmin)/(nabs-1);
  for(iord=0;iord<nord;iord++) { *ord=ordmin+iord*(ordmax-ordmin)/(nord-1);
    *col=0;
    /* for(*app=iapp=appmin;iapp<=appmax;*app=++iapp) { */
      #include "k_exec.h"
    /* } */
    i=iord*nabs+iabs; c=((int)*col)%(MAXCOLORS+1);
    if (c!=screen[i]) {Bar1(iabs,iord,c);screen[i]=c;}
  }}
  k_off();

  if (file) {
    if (ver) {snprintf(buf,MAXSTRLEN,"version %s\n", filename); system(buf);}
    if (!append) rewind(file);
    fwrite(screen,nabs*nord,1,file);
  }
}
RUN_TAIL(k_paint)

DESTROY_HEAD(k_paint)
{
  #include "k_free.h"
  FREE(S->screen);
}
DESTROY_TAIL(k_paint)

CREATE_HEAD(k_paint)
{
  DEVICE_IS_SPACELESS;
  ACCEPT_WINDOW(wnd);
  
  k_on();									CHK(NULL);

  memcpy(loctb,deftb,sizeof(*deftb));
  tb_insert_real(loctb,"abs",&(S->abs));    CHK("abs");
  tb_insert_real(loctb,"ord",&(S->ord));    CHK("ord");
  /* tb_insert_real(loctb,"app",&(S->app));    CHK("app"); */
  tb_insert_real(loctb,"col",&(S->col));    CHK("col");
  #include "k_comp.h"
  if(!used(S->data,S->ncode,&(S->col)))EXPECTED_ERROR("/*WARNING: variable \'col\' never assigned!!*/");

  k_off();

  ACCEPTL(nabs,INONE,2,INONE);
  ACCEPTR(absmin,0,INONE,INONE);
  ACCEPTR(absmax,S->nabs-1,INONE,INONE);
  ACCEPTL(nord,INONE,2,INONE);
  ACCEPTR(ordmin,0,INONE,INONE);
  ACCEPTR(ordmax,S->nord-1,INONE,INONE);

  /* ACCEPTI(appmin,0,INONE,INONE); */
  /* ACCEPTI(appmax,0,S->appmin,INONE); */

  ACCEPTI(append,1,0,1);
  ACCEPTI(ver,0,0,1);
  ACCEPTF(file,S->append?"ab":"wb",null);

  CALLOC(S->screen,S->nabs*S->nord,1);
  
}
CREATE_TAIL(k_paint,0)

#endif
