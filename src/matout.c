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
 * Dump contents of a 4D subset at selected times to a Matlab data file as a 5D array.
 * Options: output in original floating point, or discretized to desired precision. 
 */

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "mpi_io_choice.h"

/* MAT-File Data Types */
#define miINT8   (1)  /* 8-bit signed */
#define miUINT8  (2)  /* 8-bit unsigned */
#define miINT16  (3)  /* 32-bit signed */
#define miUINT16 (4)  /* 32-bit, unsigned */
#define miINT32  (5)  /* 32-bit signed */
#define miUINT32 (6)  /* 32-bit, unsigned */
#define miSINGLE (7)  /* IEEE 754 single format */
#define miDOUBLE (9)  /* IEEE 754 double format */
#define miINT64  (12) /* 32-bit signed */
#define miUINT64 (13) /* 32-bit, unsigned */
#define miMATRIX (14) /* MATLAB array */

#define FWRITE(ptr,size,nitems,...) \
  if (fwrite((ptr),(size),(nitems),file)<(nitems)) EXPECTED_ERROR("writing %d x %d bytes to %s\n",(nitems),(size),filename)

#define HEADLEN (116)
typedef struct {
  char head [HEADLEN];		/* Text at the heading of the file */
  char filename[MAXPATH];	/* Name of the output file */
  FILE *file;			/* File handle for sequential I/O */
  char array_name[MAXPATH];	/* Matlab name of the array written */
  char array_type[32];		/* Matlab name of the type of array written */
  int v0, v1;			/* Range of layers */
  real u0, u1;			/* Discretization limits */
  int dx;			/* Output step in x-direction */
  int dy;			/* .. in y-direction */
  int dz;			/* .. in z-direction */
  int dv;			/* .. in v-direction */
  char debugname[MAXPATH];	/* the file to which we print the debugging info, */
  FILE *debug;			/*    may or may not coincide with common debug file */
  int intact;			/* If we are to output the whole grid 1:1 */
  long numel;			/* Number of elements per write */
  int ArrayType;		/* Type of MAT-file array */
  int TypeLength;		/* Its length */
  int nt;			/* Total number of time points written so far */
  long int datasize_pos;		/* Position of the data element size field */
  long int datastart_pos;		/* Measure data size from here */
  long int nt_pos;		/* Position of the field containing the 5th dimension (time) */
  long int arraysize_pos;		/* Position of the array subelement size field */
  long int arraystart_pos;	/* Measure array size from here */
  long int current_pos;		/* Current end of array: maybe shorter than end of file by padding */
} STR;

/* Message to a local or to the common debug file */
/* #define DEBUG(...) if (debug) {fprintf(debug, __VA_ARGS__); FFLUSH(debug);} */

static int measure_size(STR *S, int nt);

/*  ----------------------------------------------------------------// */
RUN_HEAD(matout)
{
  DEVICE_CONST(FILE *,file);		/* File should be direct ("random") access */
  DEVICE_ARRAY(char,filename);		/* its name */
  DEVICE_VAR(int,nt);			/* Total number of time points written so far */
  DEVICE_CONST(int,intact);		/* If we are to output the whole grid 1:1 */
  DEVICE_CONST(long,numel);		/* Number of elements per write; computed at start time if needed */
  DEVICE_CONST(int,ArrayType);		/* Type of MAT-file array */
  DEVICE_CONST(int,TypeLength);		/* Type of MAT-file array */
  DEVICE_CONST(int,dx);			/* Output step in x-direction */
  DEVICE_CONST(int,dy);			/* .. in y-direction */
  DEVICE_CONST(int,dz);			/* .. in z-direction */
  DEVICE_CONST(int,dv);			/* .. in v-direction */
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  DEVICE_CONST(real,u0);		/* lower end of the discretization range */
  DEVICE_CONST(real,u1);		/* upper end of the discretization range */
  int x, y, z, v; 			/* Loop indices */
  real u;				/* Grid value */
  char buffer[16];			/* Buffer for output values. Must be longer than any supported datatype */
  
  if (!file) return 1;
  if (intact) {
    FWRITE(New,sizeof(real),numel);
  } else { /* if intact */
    /* "Chinese code" for the type cases. */
    /* NB the order of indices must agree with that */
    /* defined by index() macro to allow "intact" option. */
    #define DISCLOOP(TYPE,MINVALUE,MAXVALUE)				\
      for (x=s.x0; x<=s.x1; x+=dx) {					\
	for (y=s.y0; y<=s.y1; y+=dy) {					\
	  for (z=s.z0; z<=s.z1; z+=dz) {				\
	    for (v=v0; v<=v1; v+=dv) {				\
	      u=New[ind(x,y,z,v)];					\
	      *(TYPE *)buffer=(u<u0)?MINVALUE:(u>u1)?MAXVALUE:(TYPE)((MINVALUE*(u1-u)+MAXVALUE*(u-u0))/(u1-u0)); \
	      FWRITE(buffer,TypeLength,1); \
	    } /* for v */						\
	  } /* for z */							\
	} /* for y */							\
      } /* for x */
    #define CONTLOOP(TYPE) 						\
      for (x=s.x0; x<=s.x1; x+=dx) {					\
	for (y=s.y0; y<=s.y1; y+=dy) {					\
	  for (z=s.z0; z<=s.z1; z+=dz) {				\
	    for (v=v0; v<=v1; v+=dv) {				\
	      u=New[ind(x,y,z,v)];					\
	      *(TYPE *)buffer=(TYPE)u;					\
	      FWRITE(buffer,TypeLength,1);				\
	    } /* for v */						\
	  } /* for z */							\
	} /* for y */							\
      } /* for x */
    switch (ArrayType) {
    case miINT8:   DISCLOOP(char,SCHAR_MIN,SCHAR_MAX); break;
    case miUINT8:  DISCLOOP(unsigned char,0,UCHAR_MAX); break;
    case miINT16:  DISCLOOP(short,SHRT_MIN,SHRT_MAX); break;
    case miUINT16: DISCLOOP(unsigned short,0,USHRT_MAX); break;
    case miINT32:  DISCLOOP(int,INT_MIN,INT_MAX); break;
    case miUINT32: DISCLOOP(unsigned int,0,UINT_MAX); break;
    case miINT64:  DISCLOOP(long,LONG_MIN,LONG_MAX); break;
    case miUINT64: DISCLOOP(unsigned long,0,ULONG_MAX); break;
    case miSINGLE: CONTLOOP(float); break;
    case miDOUBLE: CONTLOOP(double); break;
    default: EXPECTED_ERROR("illegal ArrayType=%d\n");
    } /* switch ArrayType */
  } /* if intact else */
  (*nt)++;
  if NOT(measure_size(S,*nt)) return 0;
  FFLUSH(file);
  MESSAGE("matout [%d:%d:%d]x[%d:%d:%d]x[%d:%d:%d]x[%d:%d:%d]x[0:%d] to %s at t=%ld\n",
	  s.x0, dx, s.x1, s.y0, dy, s.y1, s.z0, dz, s.z1, v0, dv, v1, (*nt)-1, filename,t);
}
RUN_TAIL(matout)

/* ---------------------------------------------------------------- */

DESTROY_HEAD(matout) {
  SAFE_CLOSE(S->file);
  SAFE_CLOSE(S->debug);
} DESTROY_TAIL(matout)

/*  ----------------------------------------------------------------// */

/* Matlab array types (classes) */
#define mxDOUBLE_CLASS (6)  /* Double precision array  	*/
#define mxSINGLE_CLASS (7)  /* Single precision array  */
#define mxINT8_CLASS   (8)  /* 8-bit, signed integer   */
#define mxUINT8_CLASS  (9)  /* 8-bit, unsigned integer  */
#define mxINT16_CLASS  (10) /* 16-bit, signed integer   */
#define mxUINT16_CLASS (11) /* 16-bit, unsigned integer */
#define mxINT32_CLASS  (12) /* 32-bit, signed integer   */
#define mxUINT32_CLASS (13) /* 32-bit, unsigned integer */
#define mxINT64_CLASS  (14) /* 64-bit, signed integer   */
#define mxUINT64_CLASS (15) /* 64-bit, unsigned integer */

/* Flags - not used now */
#define flagLOGICAL (1<<1)
#define flagGLOBAL  (1<<2)
#define flagCOMPLEX (1<<3)

/* Writing fields of various types */
#define TEXT(len,val) {char *buf=calloc(len+1,1); strncpy(buf,(char *)(val),len); FWRITE(buf,1,len); free(buf);}
#define INT1(val)     {unsigned char buf=(unsigned char)(val); FWRITE(&buf,1,1);}
#define INT2(val)     {short buf=(short)(val); FWRITE(&buf,2,1);}
#define INT4(val)     {int buf=(int)(val); FWRITE(&buf,4,1);}
#define REAL8(val)    {double buf=(double)(val); FWRITE(&buf,8,1);}
#define DUMMY(len)    if (len) {char *dummy=(char *)calloc(len,1); FWRITE(dummy,1,len); free(dummy);}
#define PADDING(len)  DUMMY(len)

/* Remember position in the file */
#define MARK(p) S->p##_pos = ftell(S->file)
#define SEEK(p) fseek(S->file, S->p##_pos, SEEK_SET)
#define STUB(len,p)  MARK(p); DUMMY(len)

CREATE_HEAD(matout)
{  
  DEVICE_IS_RECTANGULAR;

  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT (v1 >= v0);
  Space s=dev->s;
  char defaulthead[HEADLEN+1]={'\0'};
  char safehead[HEADLEN+1]={'\0'};
  unsigned char ArrayClass;
  unsigned char flags;
  int x, y, z, v, xlen, ylen, zlen, vlen;
  char ArrayName[MAXPATH]={'\0'};
  int ArrayNameLength;
  /* int ArrayType; */
  time_t now;

  /* Get the user defined parameters */
  ACCEPTF(file,"wb",NULL);
  ACCEPTS(array_name,"beatbox");
  ACCEPTS(array_type,"native");
  if (strcasestr(array_type,"int")) {
    ACCEPTR(u0,RNONE,RNONE,RNONE);
    ACCEPTR(u1,RNONE,RNONE,RNONE);
    ASSERT(u0 < u1);
  }
  ACCEPTI(dx,1,1,xmax);
  ACCEPTI(dy,1,1,ymax);
  ACCEPTI(dz,1,1,zmax);
  ACCEPTI(dv,1,1,vmax);
  ACCEPTF(debug,"wt","");

  /* Prepare information for writing */
  time(&now);
  snprintf(defaulthead,HEADLEN,"Written at %24.24s by %s (%s) compiled %s %s",
	  ctime(&now),VERSTRING,VARIATION,__DATE__,__TIME__);
  ACCEPTS(head,defaulthead);
  /* Futile: either accepts is protected against too long user input, */
  /* or it isn't .. */
  strncpy(safehead,head,HEADLEN);
  strncpy(head,safehead,HEADLEN);

  if (0==strcasecmp(array_type,"native")) {
    /* Definitely will do conversion */
    switch (sizeof(real)) {
    case 4: strcpy(S->array_type,"single"); break;
    case 8: strcpy(S->array_type,"double"); break;
    default: EXPECTED_ERROR("what sort of strange reals of length % are we dealing with??\n", (int)sizeof(real));
    }
    MESSAGE("/* output type will be %s */\n",array_type);
  } else { /* if native */
    /* Perhaps will output whole */
    if (s.x0==0 && s.x1==xmax-1 &&
	s.y0==0 && s.y1==ymax-1 &&
	s.z0==0 && s.z1==zmax-1 &&
	v0==0 && v1==vmax-1 &&
	dx==1 && dy==1 && dz==1 && dv==1) {
      S->intact=1;
    }
  } /* if native else */
  
  if (0) {}
  #define CASE(a) else if (0==strcasecmp(array_type,a))
  CASE("int8")   {ArrayClass=mxINT8_CLASS;   S->ArrayType=miINT8;   S->TypeLength=1;}
  CASE("uint8")  {ArrayClass=mxUINT8_CLASS;  S->ArrayType=miUINT8;  S->TypeLength=1;}
  CASE("int16")  {ArrayClass=mxINT16_CLASS;  S->ArrayType=miINT16;  S->TypeLength=2;}
  CASE("uint16") {ArrayClass=mxUINT16_CLASS; S->ArrayType=miUINT16; S->TypeLength=2;}
  CASE("int32")  {ArrayClass=mxINT32_CLASS;  S->ArrayType=miINT32;  S->TypeLength=4;}
  CASE("uint32") {ArrayClass=mxUINT32_CLASS; S->ArrayType=miUINT32; S->TypeLength=4;}
  CASE("int64")  {ArrayClass=mxINT64_CLASS;  S->ArrayType=miINT64;  S->TypeLength=8;}
  CASE("uint64") {ArrayClass=mxUINT64_CLASS; S->ArrayType=miUINT64; S->TypeLength=8;}
  CASE("single") {ArrayClass=mxSINGLE_CLASS; S->ArrayType=miSINGLE; S->TypeLength=4;}
  CASE("double") {ArrayClass=mxDOUBLE_CLASS; S->ArrayType=miDOUBLE; S->TypeLength=8;}
  #undef CASE

  #define COUNT(q) q##len=0; for (q=dev->s.q##0;q<=dev->s.q##1;q+=d##q) q##len++
  COUNT(x);
  COUNT(y);
  COUNT(z);
  #undef COUNT
  vlen=0; for (v=v0;v<=v1;v+=dv) vlen++;
  S->numel=xlen*ylen*zlen*vlen;
  MESSAGE("/* xlen=%d ylen=%d zlen=%d vlen=%d numel=%ld */\n",xlen,ylen,zlen,vlen,S->numel);
  
  flags=0;
  strncpy(ArrayName,array_name,MAXPATH);
  if (0!=strcmp(array_name,ArrayName)) MESSAGE("/* Warning: array name will be '%s' */\n",ArrayName);
  ArrayNameLength=strlen(ArrayName);
  if (ArrayNameLength%8) ArrayNameLength=8*((ArrayNameLength/8)+1);

  /* Header */
  /* Descriptive text */ 	/* sybsys data offset */
  TEXT(116,head); 		DUMMY(4);
  /* sybsys data offset */      /* version */ /* endian indicator */
  DUMMY(4); 			INT2(0x0100); TEXT(2,"IM");

  /* Data element */
  /* data type */ /* number of bytes */
  INT4(miMATRIX); STUB(4,datasize); MARK(datastart);	 	/* tag */
  /* .. subelements */
  INT4(miUINT32); INT4(8); 					/* flags */
  INT1(ArrayClass); INT1(flags); INT2(0);  INT4(0);    		/* flags */

  INT4(miINT32); INT4(5*4);	 				/* dimensions */
  INT4(vlen); INT4(zlen);					/* dimensions */
  INT4(ylen); INT4(xlen);					/* dimensions */
  STUB(4,nt); PADDING(4);					/* dimensions */

  INT4(miINT8); INT4(ArrayNameLength);				/* name */
  TEXT(ArrayNameLength,ArrayName);				/* name */
  INT4(S->ArrayType); STUB(4,arraysize); MARK(arraystart);	/* array */

  DEBUG("datasize->%ld=%lx, datastart->%ld=%lx nt->%ld=%lx arraysize->%ld=%lx arraystart->%ld=%lx  current->%ld=%lx\n",
	S->datasize_pos,S->datasize_pos,S->datastart_pos,S->datastart_pos,S->nt_pos,S->nt_pos,
	S->arraysize_pos,S->arraysize_pos,S->arraystart_pos,S->arraystart_pos,S->current_pos,S->current_pos);

  if NOT(measure_size(S,0)) return 0;
}
CREATE_TAIL(matout,0)

/* Measure sizes, fill corresponding fields and pad up if necessary, */
/* assuming we are at the end of the array data. */
static int measure_size(STR *S, int nt)
{
  DEVICE_CONST(FILE *,file);
  DEVICE_ARRAY(char,filename);
  DEVICE_CONST(FILE *,debug);
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  long datasize, arraysize;
  int remainder;
  
  /* Take the measurements */
  MARK(current);
  datasize = (int) (S->current_pos - (S->datastart_pos));
  arraysize = (int) (S->current_pos - (S->arraystart_pos));

  /* Pad the end if necessary */
  remainder = arraysize % 8;
  if (remainder) DUMMY(8-remainder);

  /* Fill in the measurements' fields */
  SEEK(datasize);  INT4(datasize);
  SEEK(nt);        INT4(nt);
  SEEK(arraysize); INT4(arraysize);

  /* Return to the writing point */
  SEEK(current);
  DEBUG("t=%ld current->%ld=%lx nt=%d datasize=%ld arraysize=%ld padding=%d\n",
	t,S->current_pos,S->current_pos,nt,datasize,arraysize,remainder?(8-remainder):0);
  return 1;
}

