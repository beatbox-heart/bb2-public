/**
 * Copyright (C) (2010-2026) Vadim Biktashev, Irina Biktasheva et al.
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
/* 									           */
/* Works both sequentially and in MPI: each process prints its points into	   */
/* its own text buffer, every point's text tagged by its global grid position,	   */
/* and the root collects the buffers, sorts the points and outputs them, so the   */
/* file does not depend on the domain decomposition (see textbuf.h). The order of */
/* points is z-major (x fastest), in the sequential mode as well.	   	   */

#if MPI
#include <mpi.h>
#endif

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
#include "sequence.h"
#include "mpi_io_choice.h"
#include "textbuf.h"

/* runHere is only set in MPI, despite the comment in device.h */
#if MPI
#define RUNS_HERE(s) ((s).runHere)
#else
#define RUNS_HERE(s) 1
#endif

#define MAXCODE 16		/* max number of k-expressions if 'format' is used */

typedef struct {
  NFILE(file);			/* output file  		*/
  int append;			/* ..., appended if 1; rewritten if 0 */
  sequence filesequence;	/* sequence of files, as an alternative to one output file */
  char filehead[4096];		/* beginning of the file	*/
  char headformat[1024];	/* beginning of the record	*/
  char headcode[1024];		/* ..may contain numerical part */
  pp_fn headcompiled;		/* precompiled k-code for that  */
  int ncode;			/* num of k-values per field	*/
  pp_fn *code; /*[ncode]*/	/* precompiled k-codes for that */
  p_tb loctb;			/* local k_table 		*/
  real *u;			/* [vmax] k-vars		*/
  real *geom;			/* [geom_vmax] k-vars		*/
  real x, y, z;			/* k-vars promised to be real in the manual	*/
  int hasq, hasu, hasg;		/* flags of what k-code depend on 		*/
  char format[1024];		/* format for all outputs per field, if any 	*/
  char valuesep[80];		/* between values in a field 	*/
  char fieldsep[80];		/* between fields in a record	*/
  char recordsep[80];		/* between records 		*/
  textpieces mine;		/* texts of this process's points in the current record */
#if MPI
  textpieces all;		/* ..., of all processes (root only) */
#endif
  textbuf out;			/* the record as it goes to the file (root only) */
#if MPI
  MPI_Comm comm;		/* communicator for this device instance */
  int root;			/* rank of the reporting process in ALL_ACTIVE_PROCS */
  int commroot;			/* ..., and in comm */
#endif
} STR;

#undef SEPARATORS
#define SEPARATORS ";"
#define FORMATCHAR '%'

/* User's format applied to all values of a point */
static int txt_format(textbuf *b, const char *format, int ncode, REAL *v) {
  /* alas C cannot do it nicely, apparently */
  switch (ncode) {
  case  1: txt_printf(b,format,v[0]); break;
  case  2: txt_printf(b,format,v[0],v[1]); break;
  case  3: txt_printf(b,format,v[0],v[1],v[2]); break;
  case  4: txt_printf(b,format,v[0],v[1],v[2],v[3]); break;
  case  5: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4]); break;
  case  6: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5]); break;
  case  7: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6]); break;
  case  8: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7]); break;
  case  9: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8]); break;
  case 10: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9]); break;
  case 11: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
		      v[10]); break;
  case 12: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
		      v[10],v[11]); break;
  case 13: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
		      v[10],v[11],v[12]); break;
  case 14: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
		      v[10],v[11],v[12],v[13]); break;
  case 15: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
		      v[10],v[11],v[12],v[13],v[14]); break;
  case 16: txt_printf(b,format,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
		      v[10],v[11],v[12],v[13],v[14],v[15]); break;
  default:
    EXPECTED_ERROR("ncode=%d is unexpected\n",ncode);
  } /* switch ncode */
  return 1;
}

/* Sort the points of the record and assemble it in S->out, with all separators. */
/* Sequentially, the points are those of this process; in MPI, of all processes.  */
static void assemble_record(STR *S, const char *head, textpieces *rec) {
  size_t i;
  tp_sort(rec);
  S->out.n=0;
  if (S->filesequence.mask[0]) txt_puts(&(S->out),S->filehead);
  txt_puts(&(S->out),head);
  for (i=0;i<rec->n;i++) {
    if (i) txt_puts(&(S->out),S->fieldsep);
    txt_write(&(S->out),rec->txt.s+rec->pc[i].off,(size_t)rec->pc[i].len);
  }
  txt_puts(&(S->out),S->recordsep);
}

/***************/
RUN_HEAD(localprint)
{
  DEVICE_CONST(FILE *,file);
  DEVICE_VAR(sequence,filesequence);
  DEVICE_ARRAY(char,headformat);
  DEVICE_CONST(pp_fn,headcompiled);
  DEVICE_CONST(int,ncode);
  DEVICE_ARRAY(pp_fn,code);
  DEVICE_CONST(int,hasq);
  DEVICE_CONST(int,hasu);
  DEVICE_CONST(int,hasg);
  DEVICE_ARRAY(char,format);
  DEVICE_ARRAY(char,valuesep);
  char head[MAXSTRLEN]="";
#if MPI
  int imroot=(mpi_rank==S->root);
#else
#define imroot 1
#endif

  if (!filesequence->mask[0] && !S->filename[0]) return 1;

  k_on();
  if (imroot && *headformat) {
    if (headcompiled) {
      snprintf(head,MAXSTRLEN,headformat,(real)(*(REAL *)execute(headcompiled)));
    } else {
      snprintf(head,MAXSTRLEN,"%s",headformat);
    }
  }

  /* Texts of all points of this process, each tagged by its position in the grid */
  tp_reset(&(S->mine));
#define COMMANDS							\
  {									\
    int icode;								\
    p_vb *result;							\
    REAL val[MAXCODE];							\
    tp_begin(&(S->mine),tp_key(*x,*y,*z));				\
    if (hasq) {S->x=*x; S->y=*y; S->z=*z;}				\
    if (hasu) memcpy(S->u,u,vmax*sizeof(real));				\
    if (hasg) memcpy(S->geom,Geom+geom_ind((*x),(*y),(*z),0),geom_vmax*sizeof(real)); \
    for (icode=0;icode<ncode;icode++) {					\
      result = execute(code[icode]);					\
      CHK(NULL);							\
      if (*format) {							\
	val[icode]=*(REAL *)result;					\
      } else {								\
	if (icode) txt_puts(&(S->mine.txt),valuesep);			\
	txt_puts(&(S->mine.txt),prt(result,res_type(code[icode])));	\
      } /* if format else */						\
    } /* for icode */							\
    if (*format && ncode) txt_format(&(S->mine.txt),format,ncode,val);	\
    tp_end(&(S->mine));							\
  } /* COMMANDS */

  DO_FOR_ALL_POINTS;
#undef COMMANDS
  k_off();

#if MPI
  if (!tp_gather(&(S->mine),&(S->all),S->commroot,S->comm)) return 0;
  if (imroot) assemble_record(S,head,&(S->all));
#else
  assemble_record(S,head,&(S->mine));
#endif

  /* Now the root has the whole record in S->out and can write it out */
  if (filesequence->mask[0]) {
    thisq(filesequence);
#if MPI
    /* An MPI file is not truncated on opening; this is collective. */
    MPIDO(MPI_File_set_size(filesequence->f,0),"localprint: could not truncate file");
    if (imroot)
      MPIDO(MPI_File_write_at(filesequence->f,0,S->out.s,(int)S->out.n,MPI_CHAR,MPI_STATUS_IGNORE),
	    "localprint: could not write to file");
#else
    fwrite(S->out.s,1,S->out.n,filesequence->f);
#endif
    nextq(filesequence);
  } else if (file) {			/* only the root has it */
    fwrite(S->out.s,1,S->out.n,file);
    fflush(file);
  }
#undef imroot
}
RUN_TAIL(localprint)

/**********************/
DESTROY_HEAD(localprint)
{
  int icode;
  for(icode=0;icode<S->ncode;icode++) FREE(S->code[icode]);
  FREE(S->code);
  FREE(S->headcompiled);
  FREE(S->u);
  FREE(S->geom);
  tp_free(&(S->mine));
#if MPI
  tp_free(&(S->all));
#endif
  txt_free(&(S->out));
  SAFE_CLOSE(S->file);
#if !MPI
  SAFE_CLOSE(S->filesequence.f);
#endif
}
DESTROY_TAIL(localprint)

/*********************/
CREATE_HEAD(localprint)
{
  CALLOC(S->u,vmax,sizeof(real));
  if (geom_vmax) CALLOC(S->geom,geom_vmax,sizeof(real));

  ACCEPTI(append,1,0,1);
#if MPI
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    EXPECTED_ERROR("Could not create communicator.\n");
  if (dev->s.runHere) {
    /* Where is the root in the device's communicator */
    int mine=-1, me;
    MPI_Comm_rank(S->comm,&me);
    if (mpi_rank==S->root) mine=me;
    MPIDO(MPI_Allreduce(&mine,&(S->commroot),1,MPI_INT,MPI_MAX,S->comm),
	  "localprint: could not identify the root");
  }
  /* Only the root opens the file; others only need to know its name */
  if (mpi_rank==S->root) {
    if (!acceptf("file=",append?"at":"wt","",S->filename,&(S->file),w)) return 0;
  } else {
    if (!accepts("file=",S->filename,"",w)) return 0;
    S->file=NULL;
  }
#else
  if (!acceptf("file=",append?"at":"wt","",S->filename,&(S->file),w)) return 0;
#endif
  ACCEPTQ(filesequence,"wt","");
  /* In MPI, filesequence is only read by the instances that run */
  if (RUNS_HERE(dev->s)) {
    if (filesequence->mask[0]=='\0' && S->filename[0]=='\0')
      EXPECTED_ERROR("Must specify 'file' or 'filesequence' in '%s' device.\n",dev->n);
    if (filesequence->mask[0]!='\0' && S->filename[0]!='\0')
      EXPECTED_ERROR("Must specify only one of 'file' or 'filesequence' in '%s' device.\n",dev->n);
    if (filesequence->mask[0]!='\0' && find_key("append=",w)) {
      MESSAGE("Parameter 'append' is not used with 'filesequence' in '%s' device. This parameter will be ignored.\n",dev->n);
    }
  }
  ACCEPTS(filehead,"");
  /* Add newline to filehead if it is nonempty and does not have one already */
  if (filehead[0]!='\0' && filehead[strlen(filehead)-1]!='\n') {
    strcat(filehead,"\n");
  }
  if (S->file && *filehead) {
    fputs(filehead,S->file);
    FFLUSH(S->file);
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
    S->hasg=0;
    if (geom_vmax) {
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
      } /* for icode */
    } /* if S->ncode */
  } ENDBLOCK;

  ACCEPTS(format,"");
  if (*format) {
    char *p;
    int countformat=0;
    if (S->ncode>MAXCODE)
      EXPECTED_ERROR("'format' can only be used with up to %d expressions, not %d\n",MAXCODE,S->ncode);
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

