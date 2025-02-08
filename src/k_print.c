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

/* Write values of specified global k-expressions to file */

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
  #if MPI
  MPI_Comm comm;		/* communicator for this device instance */
  #endif
  int root; 			/* rank of the reporting thread */
  NFILE(file);			/* output file  		*/
  int append;			/* 1 if append to the file	*/
  char filehead[80];		/* beginning of the file	*/
  char headformat[80];		/* beginning of the record	*/
  char headcode[80];		/* ..may contain numerical part */
  pp_fn headcompiled;		/* precompiled k-code for that  */
  INT i, N;			/* k-vars counting fields	*/
  int ncode;			/* num of k-values per field	*/
  pp_fn *code;			/* precompiled k-codes for that */
  char valuesep[80];		/* between values in a field 	*/
  char fieldsep[80];		/* between fields in a record	*/
  char recordsep[80];		/* between records 		*/
} STR;

#undef SEPARATORS
#define SEPARATORS ";"

/***************/
RUN_HEAD(k_print) {
#define output(...) {if(mpi_rank==S->root) fprintf(file,__VA_ARGS__);}
#define outputs(s) {if(mpi_rank==S->root) fputs(s,file);}
  DEVICE_CONST(FILE *,file)
  DEVICE_ARRAY(char,headformat)
  DEVICE_CONST(pp_fn,headcompiled)
  DEVICE_CONST(INT,N) DEVICE_VAR(INT,i)
  DEVICE_CONST(int,ncode)
  DEVICE_ARRAY(pp_fn,code)
  DEVICE_ARRAY(char,valuesep)
  DEVICE_ARRAY(char,fieldsep)
  DEVICE_ARRAY(char,recordsep)
  int icode, size;
  char *p;
  
  k_on();
  if (headformat) {
    if (headcompiled) {
      output(headformat,(real)(*(REAL *)execute(headcompiled)));
    } else {
      outputs(headformat);
    }
  }
  if (ncode) {
    for((*i)=0;(*i)<N;(*i)++) {
      if ((*i)>0) output("%s", fieldsep);
  	for(icode=0;icode<ncode;icode++) {
  	  if (icode) output("%s", valuesep);
	  output("%s",prt(execute(code[icode]),res_type(code[icode])));
  	} /* for */
      } /* for *i */
    } /* if ncode */
    k_off();
    output("%s", recordsep);
    if(mpi_rank==S->root) FFLUSH(file);

} RUN_TAIL(k_print)

/*********************/
DESTROY_HEAD(k_print) {
  int icode;
  for(icode=0;icode<S->ncode;icode++) FREE(S->code[icode]);
  FREE(S->code);
  FREE(S->headcompiled);
  SAFE_CLOSE(S->file);
} DESTROY_TAIL(k_print)

/*********************/
CREATE_HEAD(k_print)
{
  DEVICE_IS_SPACELESS;
  /* DEVICE_MUST_BE_NOWHERE */

  ACCEPTI(append,1,0,1);
#if MPI
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    EXPECTED_ERROR("Could not create communicator.\n");
#else
  S->root = mpi_rank;
#endif

  ACCEPTS(filehead,"");

  if (mpi_rank==S->root) {
    ACCEPTF(file,S->append?"at":"wt","stdout");
    if(*S->filehead) {fprintf(S->file,"%s\n",S->filehead); FFLUSH(S->file);}
  } else {
    accepts("file=",&(S->filename[0]), "", w);
    S->file=NULL;
  }
  
  k_on();				
  CHK(NULL);
  memcpy(loctb,deftb,sizeof(*deftb));	
  tb_insert_int(loctb,  "i",&(S->i));   CHK("i");
  tb_insert_int(loctb,  "N",&(S->N));   CHK("N");
  /* #if MPI */
  /* tb_insert_fun(loctb,	"u",_u,4);      CHK("u"); */
  /* #endif */

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
  
  ACCEPTL(N,1,1,LNONE);  
  BEGINBLOCK("list=",buf); {
    int icode;
    char *pcode;
    char *s1=strdup(buf);
    int totalsize=0;
    
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
	S->code[icode] = compile(pcode,loctb,t_real); CHK(pcode);
	MESSAGE("\x01""\n\t%s%c",pcode,SEPARATORS[0]);
	totalsize+=sizetable[res_type(S->code[icode])];
      } /* for icode */
      if (totalsize>MAXSTRLEN) EXPECTED_ERROR("MAXSTRLEN not sufficient");
    } /* if S->ncode */
  } ENDBLOCK;

  ACCEPTS(valuesep,",");
  ACCEPTS(fieldsep," ");
  ACCEPTS(recordsep,"\n");
  
  k_off();

} CREATE_TAIL(k_print,0)

