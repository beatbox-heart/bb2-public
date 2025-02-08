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

/* "Inverse" of byteout. */
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

typedef struct {
  char file[MAXPATH];		/* input file name or mask */
  char filecode[MAXPATH];	/* k_code producing file number   */
  pp_fn compiled;		/* .., compiled */
  int v0,v1;			/* range of layers */
  real u0,u1;			/* Range of values corresponding to a byte */
  long local_size;	/* Size of the output buffer */
  unsigned char *Buf;   /* Output buffer. Maybe too big for stack so use heap */
  char debugname[MAXPATH];
  FILE *debug;
  int debugWriter;
} STR;

static inline real unbyte(unsigned b,unsigned bmax,real u0,real u1)
{
  return u0+(b*(u1-u0))/bmax;
}

#define MAXCHAR 255
RUN_HEAD(bytein)
{
  DEVICE_ARRAY(char, file)
  DEVICE_CONST(pp_fn, compiled)
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  DEVICE_CONST(real,u0);
  DEVICE_CONST(real,u1);
  DEVICE_CONST(long,local_size);
  DEVICE_ARRAY(unsigned char,Buf);
  DEVICE_CONST(int, debugWriter)
  DEVICE_CONST(FILE *,debug)

  long v, x, y, z;
  long nv = (v1-v0) + 1;
  long nx = (s.x1-s.x0) + 1;
  long ny = (s.y1-s.y0) + 1;
  long nz = (s.z1-s.z0) + 1;
  /* The order of bytes in the file:                        */
  /* v is the fastest, then x, ... and z is the slowest     */
#define BUF(x,y,z,v) (Buf[(v)+nv*((x)+nx*((y)+ny*(z)))])

  
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

  if (local_size!=fread(Buf,sizeof(*Buf),local_size,f)) EXPECTED_ERROR("bad read %d bytes from %s\n",local_size,filename);

  for(z=0;z<nz;z++){
    for(y=0;y<ny;y++) { 
      for(x=0;x<nx;x++) {
	if (!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)) {
	  for (v=0;v<nv;v++)
	    New[ind(x,y,z,v)]=unbyte(BUF(x,y,z,v),MAXCHAR,u0,u1);
	} /*  if isTissue */
      } /*  for x */
    } /*  for y */
  } /*  for z */
  
  if (debug && debugWriter) 
    fprintf(debug,"bytein [%ld:%ld]x[%ld:%ld]x[%ld:%ld]x[%ld:%ld] from %s at t=%ld\n",
	    (long)s.x0, (long)s.x1, (long)s.y0, (long)s.y1, (long)s.z0, (long)s.z1, (long)v0, (long)v1, filename, t);
  fclose(f);
  return 1;
}
RUN_TAIL(bytein)

DESTROY_HEAD(bytein)
{
  FREE(S->Buf);
  SAFE_CLOSE(S->debug);
}
DESTROY_TAIL(bytein)

CREATE_HEAD(bytein)
{
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
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT (v1 != v0);
  ACCEPTR(u0,RNONE,RNONE,RNONE);
  ACCEPTR(u1,RNONE,RNONE,RNONE);
  ASSERT(S->u1!=S->u0);

  /*  Dimensions of the local space. */
  int local_vsize = v1 - v0 + 1;
  int local_xsize = (dev->s.x1 - dev->s.x0) + 1;
  int local_ysize = (dev->s.y1 - dev->s.y0) + 1;
  int local_zsize = (dev->s.z1 - dev->s.z0) + 1;
  unsigned char *p;
  S->local_size = local_vsize*local_xsize*local_ysize*local_zsize;
  CALLOC(S->Buf, S->local_size, sizeof(unsigned char));
  
  S->debugWriter = 1;
  ACCEPTF(debug,"wt","");
}
CREATE_TAIL(bytein,0)

#undef R
#undef G
#undef B
