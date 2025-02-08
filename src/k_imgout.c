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

/* Output an image file where r,g,b components are defined by arbitrary k-expressions. */
/* These k-expressions are to be in the [0,1] range */

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
#include "pipe.h"

extern int Verbose;            /* defined in main */

typedef struct {
  int append;
  char filter[1024];		/* generated image will be piped to this unix command */
  char filtercode[1024];        /* calculate a number to make part of the filter */
  pp_fn compiled;		/* k_code for the filter number calculation */
  #include "k_code.h"
  INT i, j;
  int width, height;
  REAL r, g, b;
} STR;

#define MAXCHAR 255
static unsigned int byteof(real u) {
  int i=u*MAXCHAR;
  if (i<0) return 0;
  if (i>MAXCHAR) return MAXCHAR;
  return i;
}

/***************/
RUN_HEAD(k_imgout) {
  DEVICE_ARRAY(char, filter)
  DEVICE_CONST(pp_fn, compiled)
  #include "k_def.h"
  DEVICE_VAR(INT,i) DEVICE_CONST(int,width) 
  DEVICE_VAR(INT,j) DEVICE_CONST(int,height) 
  DEVICE_VAR(double,r) DEVICE_VAR(double,g) DEVICE_VAR(double,b)
  char l[2048];
  PIPE *p;

  k_on();
  snprintf(l,2048,filter,*(real *)execute(compiled));
  k_off();
  if NOT(p=pipeto(l)) return 1;

  fprintf (p->f,"P6\n%d %d\n%d\n",width,height,MAXCHAR);

  k_on();
  for ((*j)=0;(*j)<height;(*j)++) {
    for ((*i)=0;(*i)<width;(*i)++) {
      *r=*g=*b=0;
      #include "k_exec.h"
      putc(byteof(*r),p->f);
      putc(byteof(*g),p->f);
      putc(byteof(*b),p->f);
      /* fprintf(stdout,"%ld %ld %g %g %g\n",*i,*j,*r,*g,*b); */
    }
  }
  k_off();
  pipeclose(p);
} RUN_TAIL(k_imgout)

/*******************/
DESTROY_HEAD(k_imgout)
  #include "k_free.h"
  FREE(S->compiled);
DESTROY_TAIL(k_imgout)

/*******************/
CREATE_HEAD(k_imgout) {
  DEVICE_IS_SPACELESS;

  k_on();									CHK(NULL);

  ACCEPTS(filter,"cat - > t=%06.0f.ppm");
  ACCEPTS(filtercode,"t");
  S->compiled=compile(S->filtercode,deftb,t_real); CHK(S->filtercode);

  ACCEPTI(width,INONE,2,INONE);
  ACCEPTI(height,INONE,2,INONE);

  memcpy(loctb,deftb,sizeof(*deftb));
  tb_insert_abstract(loctb,"i",t_int,&(S->i),0,f_ro); CHK("i");
  tb_insert_abstract(loctb,"j",t_int,&(S->j),0,f_ro); CHK("j");
  tb_insert_real(loctb,"r",&(S->r));    CHK("r");
  tb_insert_real(loctb,"g",&(S->g));    CHK("g");
  tb_insert_real(loctb,"b",&(S->b));    CHK("b");
  #include "k_comp.h"
  if(!used(S->data,S->ncode,&(S->r)))EXPECTED_ERROR("/*WARNING: variable \'r\' never assigned!!*/");
  if(!used(S->data,S->ncode,&(S->g)))EXPECTED_ERROR("/*WARNING: variable \'g\' never assigned!!*/");
  if(!used(S->data,S->ncode,&(S->b)))EXPECTED_ERROR("/*WARNING: variable \'b\' never assigned!!*/");
  k_off();

  /**
   *	Background colour components.
   *	Used only when geometry is active, to show void points.
   *	If provided when geometry is off, the values will be ignored.
   **/
  /* if (GEOMETRY_ON) { */
  /*   ACCEPTI(bgr, 255, 0, 255); */
  /*   ACCEPTI(bgg, 255, 0, 255); */
  /*   ACCEPTI(bgb, 255, 0, 255); */
  /* } else if ( find_key("bgr=",w) || find_key("bgg=",w) || find_key("bgb=",w) ) { */
  /*   MESSAGE("Background colour parameters (bgr, bgg, bgb) are only used when geometry is active.\n\tThe value(s) provided will be ignored."); */
  /* } */

} CREATE_TAIL(k_imgout,0)



