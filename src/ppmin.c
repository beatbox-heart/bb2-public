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

/* "Inverse" of ppmout. */
/* Sequential only.     */

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "sequence.h"
#include "byte.h"

#define R 0
#define G 1
#define B 2
#define BUF(z,y,x,c) (Buf[(c)+3*((x)+local_xsize*((y)+local_ysize*(z)))])

typedef struct {
  char file[MAXPATH];		/* input file name or mask */
  char filecode[MAXPATH];	/* k_code producing file number   */
  pp_fn compiled;		/* .., compiled */
  int r, g, b;	/* layers whereto to import the colour components */
  real r0,r1;	/* range of values to transfrom node to a byte for the red comp. */
  real g0,g1;	/* -||- red comp. */
  real b0,b1;	/* -||- blue comp. */
} STR;

static inline real unbyte(unsigned b,unsigned bmax,real u0,real u1)
{
  return u0+(b*(u1-u0))/bmax;
}

#define MAXCHAR 255
RUN_HEAD(ppmin) {
  DEVICE_ARRAY(char, file)
  DEVICE_CONST(pp_fn, compiled)
  DEVICE_CONST(int,r); DEVICE_CONST(real,r0); DEVICE_CONST(real,r1);
  DEVICE_CONST(int,g); DEVICE_CONST(real,g0); DEVICE_CONST(real,g1);
  DEVICE_CONST(int,b); DEVICE_CONST(real,b0); DEVICE_CONST(real,b1);
  int nx, ny, nz, max;
  int x, y, z;
  long i;
  char tag[16], nl;
  unsigned char inbuf[3];
  size_t start, finish, total;
  
  INT filenum;			/* integer k-value */
  char filename[MAXPATH];	/* different if file is a mask */
  FILE *f;			/* file handle */

  if (compiled) {
    k_on();
    filenum=*(INT *)execute(compiled); CHK(S->filecode);
    snprintf(filename,MAXPATH,file,(long)filenum);
    k_off();
  } else {
    strncpy(filename,file,MAXPATH);
  }
  if (NOT(f=fopen(filename,"r")))
    EXPECTED_ERROR("could not open file %s for reading: %s\n",filename,strerror(errno));

  fscanf(f,"%s ",tag);
  ASSERT(0==strcmp(tag,"P6"));	   
  fscanf(f,"%u ",&nx); 
  fscanf(f,"%u ",&ny);
  fscanf(f,"%u",&max);
  fscanf(f,"%c",&nl);
  start=ftell(f);
  fseek(f,0L,SEEK_END);
  finish=ftell(f);
  fseek(f,start,SEEK_SET);
  total=finish-start;
  ASSERT(total % (nx*ny*3) == 0);
  nz=total/(nx*ny*3);
  if (Verbose) MESSAGE("File %s sizes nx=%d ny=%d nz=%d max=%d\n",filename,nx,ny,nz,max);
  ASSERT(nx==(s.x1-s.x0) + 1);
  ASSERT(ny==(s.y1-s.y0) + 1);
  ASSERT(nz==(s.z1-s.z0) + 1);
  i=0;
  for (z=s.z0;z<=s.z1;z++) {
    for (y=s.y0;y<=s.y1;y++) { 
      for (x=s.x0;x<=s.x1;x++) {
	if (3!=fread(inbuf,sizeof(*inbuf),3,f)) EXPECTED_ERROR("bad read triple number %ld from %s\n",i,filename);
	New[ind(x,y,z,r)]=unbyte(inbuf[R],max,r0,r1);
	New[ind(x,y,z,g)]=unbyte(inbuf[G],max,g0,g1);
	New[ind(x,y,z,b)]=unbyte(inbuf[B],max,b0,b1);
	i++;
      } /*  for x */
    } /*  for y */
  } /*  for z */
  MESSAGE("ppmin from '%s' at t=%ld\n",filename,t);
  fclose(f);
  return 1;
} RUN_TAIL(ppmin)

DESTROY_HEAD(ppmin)
DESTROY_TAIL(ppmin)
  
CREATE_HEAD(ppmin) {
  ACCEPTS(file,NULL);
  if (file[0]=='\0') EXPECTED_ERROR("input file name or mask is required\n");
  ACCEPTS(filecode,"");
  int formatted=(NULL != strchr(file,'%'));
  int coded=('\0' != filecode[0]);
  int decision=formatted*1+coded*2;
  switch (decision) {
  case 0*1+0*2: /* not formatted, not coded */
    S->compiled=NULL;
    break;
  case 1*1+0*2: /* formatted, not coded */
    MESSAGE("/* WARNING: filename '%s' contains format symbol '%' but there is no filecode parameter */",file);
    S->compiled=NULL;
    break;
  case 0*1+1*2: /* not formatted, coded */
    MESSAGE("/* WARNING: filecode='%s' is provided but filename '%s' contains no format symbol; filecode will be ignored */",filecode,file);
    S->compiled=NULL;
    break;
  case 1*1+1*2: /* formatted, coded */
    S->compiled=compile(S->filecode,deftb,t_int); CHK(S->filecode);
    break;
  default:
    EXPECTED_ERROR("this cannot be: formatted=%d coded=%d decision=%d\n",formatted,coded,decision);
  }

  ACCEPTI(r,INONE,-1,(int)vmax-1);
  ACCEPTR(r0,RNONE,RNONE,RNONE);
  ACCEPTR(r1,RNONE,RNONE,RNONE);
  
  ACCEPTI(g,INONE,-1,(int)vmax-1);
  ACCEPTR(g0,RNONE,RNONE,RNONE);
  ACCEPTR(g1,RNONE,RNONE,RNONE);
  
  ACCEPTI(b,INONE,-1,(int)vmax-1);
  ACCEPTR(b0,RNONE,RNONE,RNONE);
  ACCEPTR(b1,RNONE,RNONE,RNONE);
} CREATE_TAIL(ppmin,0)

#undef R
#undef G
#undef B
