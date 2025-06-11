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

/* Record values in the given 4D cube to a text file */

#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "sequence.h"
#include "mpi_io_choice.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RECORD_SEPARATOR_LENGTH 2
typedef struct {
  int v0, v1;
  char filehead[1024];
  int timestamp;
  char timestampformat[1024];
  int timestampwidth;
  sequence filesequence;
  FILE *file;
  char filename[MAXPATH];
  int append;
  char format[1024];
  int formatwidth;
  char vsep[RECORD_SEPARATOR_LENGTH +1];
  char zsep[RECORD_SEPARATOR_LENGTH +1];
  char ysep[RECORD_SEPARATOR_LENGTH +1];
  char xsep[RECORD_SEPARATOR_LENGTH +1];
  char recordsep[RECORD_SEPARATOR_LENGTH +1];
} STR;

RUN_HEAD(record)
{
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
  DEVICE_ARRAY(char,filehead);
  DEVICE_CONST(int,timestamp);
  DEVICE_ARRAY(char,timestampformat);
  DEVICE_CONST(int,timestampwidth);
  DEVICE_VAR(sequence,filesequence);  
  DEVICE_CONST(FILE *,file);
  DEVICE_ARRAY(char,format);
  DEVICE_CONST(int,formatwidth);
  DEVICE_ARRAY(char,vsep);
  DEVICE_ARRAY(char,zsep);
  DEVICE_ARRAY(char,ysep);
  DEVICE_ARRAY(char,xsep);
  DEVICE_ARRAY(char,recordsep);
  int v, z, y, x;
  char *nextSep;
  FILE *f=NULL;
  if (filesequence->mask[0]) {
    thisq(filesequence);
    f=filesequence->f;
    if (*filehead) fputs(filehead,f);
  } else if (file) {
    f=file;
  }
  if (!f) return 1;

  if (timestamp) {
    char *tbuf;
    CALLOC(tbuf,timestampwidth+1,1);
    snprintf(tbuf,timestampwidth,timestampformat, (long)t);
    tbuf[timestampwidth]='\0';
    fputs(tbuf,f);
    FREE(tbuf);
  }
  
/*  Write data, selecting hierarchical separator */
/*  Record->X->Y->Z->V (slowest->fastest). */

  for (x=s.x0; x<=s.x1; x++) {
    for (y=s.y0; y<=s.y1; y++) {
      for (z=s.z0; z<=s.z1; z++) {
	for (v=v0; v<=v1; v++) {
	  snprintf(buf, formatwidth+1, format, (double) New[ind(x,y,z,v)]);
	  strcat(buf,
		 (v<v1)? vsep :
		 (z<s.z1)? zsep :
		 (y<s.y1)? ysep :
		 (x<s.x1)? xsep :
		 recordsep
		 );
	  fputs(buf,f);
	} /* for v */
      } /* for z */
    } /* for y */
  } /* for x */
  if (filesequence->mask[0]) 
    nextq(filesequence);
  else if (file)
    fflush(f);
}
RUN_TAIL(record)

DESTROY_HEAD(record)
{
  SAFE_CLOSE(S->file);
  SAFE_CLOSE(S->filesequence.f);
}
DESTROY_TAIL(record)


/* ----------------------------------------------------------------- */
CREATE_HEAD(record)
{
  DEVICE_IS_RECTANGULAR;
  
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT(v1 >= v0);
  ACCEPTS(filehead,"");
  /* Add newline to filehead if it is nonempty and does not have one already */
  if (filehead[0]!='\0' && filehead[strlen(filehead)-1]!='\n') {
    strcat(filehead,"\n");
  }	
  ACCEPTI(timestamp,0,0,1);
  ACCEPTS(timestampformat,"%-6ld\t");
  ACCEPTI(timestampwidth,8,1,INONE);
  
#define CHECK_SEPARATOR(s)						\
  if (strlen(S->s) > RECORD_SEPARATOR_LENGTH) {				\
    MESSAGE("WARNING from %s: %s is longer than the required separator size of %d and will be truncated.\n", \
	    dev->n, #s, RECORD_SEPARATOR_LENGTH);			\
  } else if (strlen(S->s) < RECORD_SEPARATOR_LENGTH) {			\
    int i;								\
    for (i=strlen(S->s); i<RECORD_SEPARATOR_LENGTH; i++) 		\
      S->s[i]=' ';							\
  }									\
  S->s[RECORD_SEPARATOR_LENGTH] = '\0';
  
  ACCEPTS(vsep,"  "); CHECK_SEPARATOR(vsep);
  ACCEPTS(zsep,", "); CHECK_SEPARATOR(zsep);
  ACCEPTS(ysep,"; "); CHECK_SEPARATOR(ysep);
  ACCEPTS(xsep," \n"); CHECK_SEPARATOR(xsep);
  ACCEPTS(recordsep,"\n\n"); CHECK_SEPARATOR(recordsep);
#undef CHECK_SEPARATOR
  
  ACCEPTS(format,"%38f");
  ACCEPTI(formatwidth,38,2,INONE);
  /* check the consistency of the two on some arbitrary chosen numbers */
  {
    int i, l;
    char s[1024];
    real testvalue[6]={0, 1, M_PI, -1, 0.1, -1.e6};
    for (i=0;i<6;i++) {
      snprintf(s, formatwidth+1, format, testvalue[i]);
      l=strlen(s);
      if (formatwidth!=l) {
	MESSAGE("/* Warning: test value "REALF" printed by format '%s' produces string '%s' of length %d not %d */\n",
		testvalue[i],format,s,l,formatwidth);
	/* break; */
      }
      /* else { */
      /* 	MESSAGE("/\* "REALF" [%s] -> '%s' [%ld] *\/\n",testvalue[i],format,s,l); */
      /* } */
    }
  }
  
  /*  Sequential version */
  ACCEPTI(append,1,0,1);
  ACCEPTF(file,append?"at":"wt","");
  ACCEPTQ(filesequence,"wt","");
  if (filesequence->mask[0]=='\0' && file==NULL)
    EXPECTED_ERROR("Must specify 'file' or 'filesequence' in '%s' device.\n",dev->n);
  if (filesequence->mask[0]!='\0' && file!=NULL)
    EXPECTED_ERROR("Must specify only one of 'file' or 'filesequence' in '%s' device.\n",dev->n);
  if (filesequence->mask[0]!='\0' && find_key("append=",w)) {
    MESSAGE("Parameter 'append' is not used with 'filesequence' in '%s' device. This parameter will be ignored.\n",dev->n);
  }
  if (file && *filehead) {
    fputs(filehead,file);
    FFLUSH(file);	
  }
}
CREATE_TAIL(record,0)

#undef RECORD_SEPARATOR_LENGTH
