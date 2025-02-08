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

/* Produce an image or image stream (one image per z value) via ppm and an arbitrary unix filter */

#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "k_.h"
#include "pipe.h"
#include "byte.h"
#include "mpi_io_choice.h"

#define CODELENGTH 1024

typedef struct {
  int append;
  char filter[CODELENGTH];	/* generated image will be piped to this unix command */
  char code[CODELENGTH];	/* calculate a number to make part of the filter */
  char filtercode[CODELENGTH];	/* .., alias of */
  pp_fn compiled;		/* k_code for the filter number calculation */
  int r, g, b;			/* layers wherefrom to extract the colour components */
  int autoscale;		/* will scale the components automatically */
  real r0,r1;			/* range of values to transfrom node to a byte for the red comp. */
  real g0,g1;			/* -||- red comp. */
  real b0,b1;			/* -||- blue comp. */
  int bgr,bgg,bgb;		/* background (no-tissue) colour */
  char debugname[MAXPATH];
  FILE *debug;
} STR;

#define MAXCHAR 255

/*#define makename(name,stem,no) (snprintf((name),MAXPATH,"%s%04d",(stem),(no)))*/
/* #define makename(name,stem,no) (strchr((stem),'%')?snprintf((name),MAXPATH,(stem),(no)):snprintf((name),MAXPATH,"%s%04d",(stem),(no))) */

static void getlimits(Space s,int c,real *c0,real *c1);

RUN_HEAD(imgout)
{
  DEVICE_ARRAY(char, filter);
  DEVICE_CONST(pp_fn, compiled);
  DEVICE_CONST(int,autoscale);
  DEVICE_CONST(int,r) DEVICE_CONST(real,r0) DEVICE_CONST(real,r1);
  DEVICE_CONST(int,g) DEVICE_CONST(real,g0) DEVICE_CONST(real,g1);
  DEVICE_CONST(int,b) DEVICE_CONST(real,b0) DEVICE_CONST(real,b1);
  DEVICE_CONST(int,bgr);
  DEVICE_CONST(int,bgg);
  DEVICE_CONST(int,bgb);
  DEVICE_CONST(FILE *,debug);
  int x0=s.x0, x1=s.x1, nx=x1-x0+1;
  int y0=s.y0, y1=s.y1, ny=y1-y0+1;
  int z0=s.z0, z1=s.z1;
  int x, y, z;
  char l[2048];
  PIPE *p;
  /* real *u; */

  k_on();
  snprintf(l,2048,filter,*(real *)execute(compiled));
  k_off();
  if NOT(p=pipeto(l)) return 1;

  if (autoscale) {
    getlimits(s,r,&r0,&r1);
    getlimits(s,g,&g0,&g1);
    getlimits(s,b,&b0,&b1);
    if (debug) fprintf(debug,"imgout t=%ld r:%g..%g g:%g..%g b:%g..%g\n",t,r0,r1,g0,g1,b0,b1);
  }

  for(z=z0;z<=z1;z++) {
    /* a separate image for each z-section */
    fprintf (p->f,"P6\n%d %d\n%d\n",nx,ny,MAXCHAR);
    for(y=y0;y<=y1;y++) {
      for(x=x0;x<=x1;x++) {
	if (isTissue(x,y,z)) {
	  putc(Byte(x,y,z,r,r0,r1),p->f);
	  putc(Byte(x,y,z,g,g0,g1),p->f);
	  putc(Byte(x,y,z,b,b0,b1),p->f);
	} else { /*  Void. Use background colour. */
	  putc((unsigned) bgr,p->f);
	  putc((unsigned) bgg,p->f);
	  putc((unsigned) bgb,p->f);
	} /*  else */
      } /* for x */
    } /* for y */
  } /* for z */
  pipeclose(p);
}
RUN_TAIL(imgout)

DESTROY_HEAD(imgout)
{
  FREE(S->compiled);
}
DESTROY_TAIL(imgout)

CREATE_HEAD(imgout)
{
  int code_accepted=0;
  
  DEVICE_IS_RECTANGULAR;

  ACCEPTS(filter,"cat - > t=%06.0f.ppm");
  ACCEPTS(code,"t");
  ACCEPTS(filtercode,S->code);
  if (find_key("code=",w) && find_key("filtercode=",w))
      MESSAGE("/* both code= and filtercode= are present; the latter takes effect */\n");
  memcpy(S->code,S->filtercode,CODELENGTH);
  k_on();
  S->compiled=compile(S->code,deftb,t_real); CHK(S->code);
  k_off();
  
  ACCEPTI(r,INONE,-1,(int)vmax-1); 
  ACCEPTI(g,INONE,-1,(int)vmax-1);
  ACCEPTI(b,INONE,-1,(int)vmax-1);
  ACCEPTI(autoscale,0,0,1);		/* will scale the components automatically */
  if (S->autoscale) {
    #define CHECK(p) if (find_key(#p"=",w)) \
      MESSAGE("warning: parameter "#p" will be ignored since autoscale=%d is not zero\n",S->autoscale)
    CHECK(r0);
    CHECK(r1);
    CHECK(g0);
    CHECK(b1);
    CHECK(b0);
    CHECK(b1);
  } else {
    ACCEPTR(r0,RNONE,RNONE,RNONE);  
    ACCEPTR(r1,RNONE,RNONE,RNONE); 
    ASSERT(S->r0!=S->r1);
    ACCEPTR(g0,RNONE,RNONE,RNONE);
    ACCEPTR(g1,RNONE,RNONE,RNONE);
    ASSERT(S->g0!=S->g1);
    ACCEPTR(b0,RNONE,RNONE,RNONE);
    ACCEPTR(b1,RNONE,RNONE,RNONE);
    ASSERT(S->b0!=S->b1);
  }
  /**
   *	Background colour components.
   *	Used only when geometry is active, to show void points.
   *	If provided when geometry is off, the values will be ignored.
   **/
  if (GEOMETRY_ON) {
    ACCEPTI(bgr, 255, 0, 255);
    ACCEPTI(bgg, 255, 0, 255);
    ACCEPTI(bgb, 255, 0, 255);
  } else if ( find_key("bgr=",w) || find_key("bgg=",w) || find_key("bgb=",w) ) {
    MESSAGE("Background colour parameters (bgr, bgg, bgb) are only used when geometry is active.\n\tThe value(s) provided will be ignored.");
  }
  ACCEPTF(debug,"wt","");
}
CREATE_TAIL(imgout,0)


static void getlimits(Space s,int c,real *c0,real *c1) {
  real u;
  int x,y,z;
  if (0<=c && c<vmax) {
    *c0=MAXREAL; *c1=-MAXREAL;
    for(z=s.z0;z<=s.z1;z++) {
      for(y=s.y0;y<=s.y1;y++) {
	for(x=s.x0;x<=s.x1;x++) {
	  if (isTissue(x,y,z)) {
	    u=New[ind(x,y,z,c)];
	    if (u<*c0) *c0=u;
	    if (u>*c1) *c1=u;
	  }
	}
      }
    }
  } else {
    *c0=0; *c1=1;
  }
}
