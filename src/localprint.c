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

/* For selected gridpoints, write values of specified local k-expressions to file. */
/* Sequential-only for the moment. */

#include <assert.h>
#include <stdarg.h>
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
#include "mpi_io_choice.h"

typedef struct {
  NFILE(file);			/* output file  		*/
  int append;			/* 1 if append to the file	*/
  char filehead[4096];		/* beginning of the file	*/
  char headformat[1024];	/* beginning of the record	*/
  char headcode[1024];		/* ..may contain numerical part */
  pp_fn headcompiled;		/* precompiled k-code for that  */
  int ncode;			/* num of k-values per field	*/
  pp_fn *code; /*[ncode]*/	/* precompiled k-codes for that */
  p_real *data; /*[ncode]*/	/* where to allocate the results before printing */
  p_tb loctb;			/* local k_table 		*/
  real *u;			/* [vmax] k-vars		*/
  real *geom;			/* [geom_vmax] k-vars		*/
  real x, y, z;			/* k-vars promised to be real in the manual	*/
  int hasq, hasu, hasg;		/* flags of what k-code depend on 		*/
  char format[1024];		/* format for all outputs per field, if any 	*/
  char valuesep[80];		/* between values in a field 	*/
  char fieldsep[80];		/* between fields in a record	*/
  char recordsep[80];		/* between records 		*/
} STR;

#undef SEPARATORS
#define SEPARATORS ";"
#define FORMATCHAR '%'

#define output(...) {fprintf(file,__VA_ARGS__);}
#define outputs(s) {fputs(s,file);}

/***************/
RUN_HEAD(localprint)
{
  DEVICE_CONST(FILE *,file);
  DEVICE_ARRAY(char,filehead);
  DEVICE_ARRAY(char,headformat);
  DEVICE_CONST(pp_fn,headcompiled);
  DEVICE_CONST(int,ncode);
  DEVICE_ARRAY(pp_fn,code);
  DEVICE_ARRAY(p_real,data);
  DEVICE_CONST(int,hasq);
  DEVICE_CONST(int,hasu);
  DEVICE_CONST(int,hasg);
  DEVICE_ARRAY(char,format);  
  DEVICE_ARRAY(char,valuesep);
  DEVICE_ARRAY(char,fieldsep);
  DEVICE_ARRAY(char,recordsep);
  int icode, size;
  char *p;
  int iout=0;

  k_on();
  if (headformat) {
    if (headcompiled) {
      output(headformat,(real)(*(REAL *)execute(headcompiled)));
    } else {
      outputs(headformat);
    }
  }

#define COMMANDS							\
  {									\
    int icode;								\
    p_vb *result;							\
    pp_fn the_code;							\
    if ((iout)>0) output("%s", fieldsep);				\
    if (hasq) {S->x=*x; S->y=*y; S->z=*z;}				\
    if (hasu) memcpy(S->u,u,vmax*sizeof(real));				\
    if (hasg) memcpy(S->geom,Geom+geom_ind((*x),(*y),(*z),0),geom_vmax*sizeof(real)); \
    for (icode=0;icode<ncode;icode++) {					\
      the_code = code[icode];						\
      if (!the_code) break;						\
      result = execute(the_code);					\
      CHK(NULL);							\
      if (*format) {							\
	memcpy(data[icode],result,sizetable[t_real]);			\
      } else {								\
	if (icode) output("%s", valuesep);				\
	output("%s",prt(execute(code[icode]),res_type(code[icode])));	\
      } /* if format else */						\
    } /* for icode */							\
    iout++;								\
  } /* COMMANDS */							

  DO_FOR_ALL_POINTS; 
#undef COMMANDS
  
  if (*format) {
    /* alas C cannot do it nicely, apparently */
    switch (ncode) {
    case  1: output(format,data[0]); break;
    case  2: output(format,data[0],data[1]); break;
    case  3: output(format,data[0],data[1],data[2]); break;
    case  4: output(format,data[0],data[1],data[2],data[3]); break;
    case  5: output(format,data[0],data[1],data[2],data[3],data[4]); break;
    case  6: output(format,data[0],data[1],data[2],data[3],data[4],data[5]); break;
    case  7: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6]); break;
    case  8: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7]); break;
    case  9: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8]); break;
    case 10: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9]); break;
    case 11: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],
		    data[10]); break;
    case 12: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],
		    data[10],data[11]); break;
    case 13: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],
		    data[10],data[11],data[12]); break;
    case 14: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],
		    data[10],data[11],data[12],data[13]); break;
    case 15: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],
		    data[10],data[11],data[12],data[13],data[14]); break;
    case 16: output(format,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8],data[9],
		    data[10],data[11],data[12],data[13],data[14],data[15]); break;
    default:
      EXPECTED_ERROR("ncode=%d is unexpected\n");
    } /* switch ncode */
  } /* if format */
  
  k_off();
  output("%s", recordsep);
}
RUN_TAIL(localprint)

/**********************/
DESTROY_HEAD(localprint)
{
  int icode;
  for(icode=0;icode<S->ncode;icode++) FREE(S->code[icode]);
  FREE(S->code);
  FREE(S->headcompiled);
  SAFE_CLOSE(S->file);
}
DESTROY_TAIL(localprint)

/*********************/
CREATE_HEAD(localprint)
{
  CALLOC(S->u,vmax,sizeof(real));
  if (geom_vmax) CALLOC(S->geom,geom_vmax,sizeof(real));
  
  ACCEPTF(file,"wt","");
  ACCEPTS(filehead,"");

  if (file && *filehead) {
    output("%s",filehead);
    FFLUSH(file);	
  }
 
  ACCEPTS(headformat,"");
  if (S->headformat[0]) { /* need curly bracket as ACCEPTS is a funny macro */
    ACCEPTS(headcode,"");
  } else {
    S->headcode[0]='\0';
  }
  if (S->headcode[0]) {
    S->headcompiled=compile(S->headcode,deftb,t_real); CHK(S->headcode);
  } else {
    S->headcompiled=NULL;
  }

  BEGINBLOCK("list=",buf); {
    int iv;
    char name[maxname];
    int icode;
    char *pcode;
    char *s1=strdup(buf);
    int totalsize=0;

    k_on();				CHK(NULL);
    S->loctb = tb_new();
    memcpy(S->loctb,deftb,sizeof(*deftb));
    S->hasq=0;
    if (k_expr_depends (buf,"x")) {
      S->hasq=1;
      tb_insert_real_ro(S->loctb,"x",&(S->x));	CHK("x");
    }
    if (k_expr_depends (buf,"y")) {
      S->hasq=1;
      tb_insert_real_ro(S->loctb,"y",&(S->y));	CHK("y");
    }
    if (k_expr_depends (buf,"z")) {
      S->hasq=1;
      tb_insert_real_ro(S->loctb,"z",&(S->z));	CHK("z");
    }
    S->hasu=0;
    for (iv=0;iv<vmax;iv++) {
      snprintf(name,maxname,"u%d",iv);
      if (k_expr_depends (buf,name)) {
	S->hasu=1;
	tb_insert_real_ro(S->loctb,name,&(S->u[iv]));  CHK(name);
      }
    } /* for iv */
    if (geom_vmax) {
      S->hasg=0;
      for (iv=0;iv<geom_vmax;iv++) {
	snprintf(name,maxname,"geom%d",iv);
	if (k_expr_depends (buf,name)) {
	  S->hasg=1;
	  tb_insert_real_ro(S->loctb,name,&(S->geom[iv]));  CHK(name);
	}
      } /* for iv */
    } /* if geom_vmax */
  
    for(S->ncode=0;NULL!=(pcode=strtok(S->ncode?NULL:s1,SEPARATORS));S->ncode+=(*pcode!=0));
    if (!S->ncode) MESSAGE("/*WARNING: no expressions in \"%s\"*/",buf);
    FREE(s1);
    
    if (S->ncode) {  
      if NOT(S->code=calloc(S->ncode,sizeof(pp_fn)))
	      ABORT("not enough memory for code array of %d",S->ncode);
      
      for(icode=0;icode<S->ncode;icode+=(*pcode!=0)) {
	if NOT(pcode=strtok(icode?NULL:buf,SEPARATORS)) EXPECTED_ERROR("internal error");
	if (!*pcode) continue;
	/* here we presume ALL expression shall return real incl int() k-function */
	S->code[icode] = compile(pcode,S->loctb,t_real); CHK(pcode);
	MESSAGE("\x01""\n\t%s%c",pcode,SEPARATORS[0]);
	totalsize+=sizetable[res_type(S->code[icode])];
      } /* for icode */
      if (totalsize>MAXSTRLEN) EXPECTED_ERROR("MAXSTRLEN not sufficient");
    } /* if S->ncode */
  } ENDBLOCK;

  ACCEPTS(format,"");
  if (*format) {
    char *p;
    int countformat=0;
    for (p=format;*p;p++) if (*p==FORMATCHAR) countformat++;
    if (countformat<S->ncode) {
      MESSAGE("WARNING: format='%s' has too few '%c' signs (%d) for %d code results\n",
	      format,FORMATCHAR,countformat,S->ncode);
    }
    if (find_key("valuesep=",w)) {
      MESSAGE("Parameter 'format' is specified hence 'valuesep' is spurious and will be ignored\n");
    }
  } else {
    ACCEPTS(valuesep,",");
  }
  ACCEPTS(fieldsep," ");
  ACCEPTS(recordsep,"\n");
  
  k_off();
}
CREATE_TAIL(localprint,0)

