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

/* Generic device executing k-codes for grid and global var data, 
  including "phase" inital conditions */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"

#include "k_.h"

/*++++++++++++++++++++++++++++++++++++++++*/ extern int idev;

#undef SEPARATORS
#define SEPARATORS ";"
#define BLANK " \n\t\r"

typedef struct {
  int advance;
  #include "k_code.h"
  real *u;		/*[vmax]*/
  real *geom;		/*[geom_vmax]*/
  /* int restricted; */
  real x, y, z;
  INT np, nv;
  real *a;		/*[np*nv]*/
  double *p;		/*[nv]*/
  double phaseu, phasep;
  NFILE(file);
  NFILE(debug);
  p_tb loctb;
  int debugWriter;
} STR;

RUN_HEAD(k_func)
{
  DEVICE_CONST(INT, np)
  DEVICE_CONST(INT, nv)
  DEVICE_CONST(int, debugWriter)
  DEVICE_CONST(FILE *,debug)
  DEVICE_CONST(p_tb, loctb)
  SPACE_CONST(int,restricted)
  int ix, iy, iz, iv, ip, ip1, icode, n, touchPoint;
  double p, q;

  /* DEBUG("enter\n"); */

  if (debug && debugWriter) {
    fprintf(debug, " t=%ld", (long)t);
    fflush(debug);
  }
  k_on();
  assert(DV==1);
  for (S->z=iz=s.z0;iz<=s.z1;S->z=++iz) {
    for (S->y=iy=s.y0;iy<=s.y1;S->y=++iy) {
      for (S->x=ix=s.x0;ix<=s.x1;S->x=++ix) {
	if (!s.nowhere && debug) {
	  fprintf(debug," (%d %d %d)",(int)ix,(int)iy,(int)iz);
	  fflush(debug);
	}
		
	touchPoint = !s.nowhere && (restricted==0 || isTissue(ix,iy,iz));
	/* No need to maintain local values if running nowhere. */
	if(touchPoint) memcpy(S->u,New+ind(ix,iy,iz,0),vmax*sizeof(real));
	if(geom_vmax) memcpy(S->geom,Geom+geom_ind(ix,iy,iz,0),geom_vmax*sizeof(real));
		
	for(icode=0;icode<S->ncode;icode++) {
	  void *code, *result;
	  if NOT(code=(S->code)[icode]) continue;
	  result = execute(code); CHK(NULL);
	  if ((!s.nowhere || debugWriter) && debug) {
	    fprintf(debug," %s=%s",var_name(loctb,(S->data)[icode],n), prt(result, res_type(code)));
	    fflush(debug);
	  }
	  memcpy((S->data)[icode],result,sizetable[res_type(code)]);
		
	  if (((S->data)[icode])==&(S->phasep)) {
	    if (np) {
	      S->phasep = S->phasep*(np)*M_1_PI*0.5;
	      while(S->phasep<0) S->phasep=S->phasep+(np);
	      while(S->phasep>=(np)) S->phasep=S->phasep-(np);
	      ip=S->phasep; ip1=(ip+1)%np; p=S->phasep-ip; q=1-p;
	      assert(ip<np);
	      for (iv=0;iv<nv;iv++){
		S->p[iv]=S->a[iv+nv*ip]*q+S->a[iv+nv*ip1]*p;
	      }
	    } else {
	      for (iv=0;iv<nv;iv++){
		S->p[iv]=S->a[iv+nv*ip];
	      }
	    }
	  }
	  
	  if (((S->data)[icode])==&(S->phaseu)) {
	    if (np) {
	      S->phaseu = S->phaseu*(np)*M_1_PI*0.5;
	      while(S->phaseu<      0) S->phaseu=S->phaseu+(np);
	      while(S->phaseu>=(np)) S->phaseu=S->phaseu-(np);
	      ip=S->phaseu; ip1=(ip+1)%np; p=S->phaseu-ip; q=1-p;
	      assert(ip<np);
	      for (iv=0;iv<nv;iv++){
		S->u[iv]=S->a[iv+nv*ip]*q+S->a[iv+nv*ip1]*p;
	      }
	    } else {
	      for (iv=0;iv<nv;iv++) S->u[iv]=S->a[iv+nv*ip];
	    }
	  }
	}
	
	if(touchPoint) memcpy(New+ind(ix,iy,iz,0),S->u,vmax*sizeof(real));
	if(geom_vmax) memcpy(Geom+geom_ind(ix,iy,iz,0),S->geom,geom_vmax*sizeof(real));
      } /* for ix */
    } /* for iy */
  } /* for iz */
  k_off();
  if (debug && debugWriter) {
    fprintf(debug,"\n");
    fflush(debug);
  }
  /* DEBUG("exit\n"); */
}
RUN_TAIL(k_func)

DESTROY_HEAD(k_func)
{
  #include "k_free.h"
  FREE(S->u);
  if (S->geom) FREE(S->geom);
  FREE(S->loctb);
  /* This is the first instance of S->file in k_func.c. */
  /* We can check that it was opened (or at least the pointer defined) somewhere else  */
  /* by a simple hack using ftell. */
  SAFE_CLOSE(S->file);
  SAFE_CLOSE(S->debug);
} 
DESTROY_TAIL(k_func)

CREATE_HEAD(k_func)
{
  #define DELIM "\\/,; \r\t\n"
  char *p, name[80];
  int iv, ip, nv, np;
  real *array, *pv;

  ACCEPTI(advance,0,0,1);
  /* ACCEPTI(restricted,1,0,1); */

  CALLOC(S->u,vmax,sizeof(real));
  if (geom_vmax) CALLOC(S->geom,geom_vmax,sizeof(real));
	
#if MPI
  /* Identifying the root process allows only one process to write debug output. */
  if (dev->s.nowhere) {
    S->debugWriter = (mpi_rank == 0);
  } else {
    S->debugWriter = (mpi_rank == getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0));
  }
#else
  S->debugWriter = 1;
#endif
	
  ACCEPTF(file,"rt","");
  if(S->file) {
    fgets(buf,MAXSTRLEN,S->file);	/* Count words in the first line */
    if NOT(p=strtok(buf,DELIM)) EXPECTED_ERROR("first line of %s seems empty, cannot recognize",S->filename);
    for(nv=1;NULL!=(p=strtok(NULL,DELIM));nv++);
    rewind(S->file);				/* Count lines in the file */
    for(np=0;fgets(buf,MAXSTRLEN,S->file);np++);
    switch(np) {
    case 0: EXPECTED_ERROR("no lines in %s",S->filename);
    case 1: MESSAGE("/* WARNING: one line found in %s, phase values will be ingored",S->filename);
    }
    rewind(S->file);					/* read the array from file */
    MESSAGE("\n/* Table of %d x %d is assumed in %s */\n",nv,np,S->filename);
    CALLOC(array,nv*np,sizeof(real));
    CALLOC(pv,nv,sizeof(real));
    for (ip=0;ip<np;ip++) {
      if (ip>0){
	for(iv=0;iv<nv;iv++){
	  array[iv+ip*nv]=array[iv+(ip-1)*nv];
	}
      }
      fgets(buf,MAXSTRLEN,S->file);
      p = strtok(buf,DELIM);
      iv=0;
      while (p!=NULL) {
	sscanf(p,REALF,array+iv+ip*nv);
	iv++;
	p = strtok(NULL,DELIM);
      } /*  while */
    } /*  for */
    /* fclose(S->file); */ // this is commented out, as the file is closed elsewhere.
    S->a = array;
    S->np = np;
    S->nv = nv;
    S->p  = pv;
  } else {
    S->a=S->p=NULL;
    S->np=S->nv=0;
  }
  
  ACCEPTF(debug,"wt","");
  k_on();				CHK(NULL);
  S->loctb = tb_new();
  memcpy(S->loctb,deftb,sizeof(*deftb));
  tb_insert_real_ro(S->loctb,"x",&(S->x));	CHK("x");
  tb_insert_real_ro(S->loctb,"y",&(S->y));	CHK("y");
  tb_insert_real_ro(S->loctb,"z",&(S->z));	CHK("z");

  /* Read-only k-variables np and nv count dimensions of the array */
  /* contained in the file; if there is no file, there is no np nor nv */
  if (S->file) {
    tb_insert_int_ro(S->loctb,"np",&(S->np)); CHK("np");
    tb_insert_int_ro(S->loctb,"nv",&(S->nv)); CHK("nv");
  }
	
	/*  Random function is safe for local use only. */
	/* tb_insert_fun(S->loctb,"rnd",_rnd,2);	CHK("rnd"); */
	
  for(iv=0;iv<vmax;iv++) {
    snprintf(name,80,"u%d",iv);
    tb_insert_real(S->loctb,name,&(S->u[iv]));  CHK(name);
  }
  if (geom_vmax) {
    for(iv=0;iv<geom_vmax;iv++) {
      snprintf(name,80,"geom%d",iv);
      tb_insert_real(S->loctb,name,&(S->geom[iv]));  CHK(name);
    }
  }
	
  if (S->file) {
    tb_insert_real(S->loctb,"phasep",&(S->phasep));	CHK("phasep");
    tb_insert_real(S->loctb,"phaseu",&(S->phaseu));	CHK("phaseu");
    for(iv=0;iv<nv;iv++) {
      snprintf(name,80,"p%d",iv);
      tb_insert_real(S->loctb,name, &(S->p[iv]));	CHK(name);
    }
  }
  #define loctb (S->loctb)
  #include "k_comp.h"
  #undef loctb
  if (!S->ncode)
    EXPECTED_ERROR("no statements in \"%s\"",buf);
  if (S->file) {
    if(!used(S->data,S->ncode,&(S->phaseu)) && !used(S->data,S->ncode,&(S->phasep))){
      MESSAGE("/*WARNING: file was specified but \'phase[up]\' not used in the expression!!*/");
    }
  }
  k_off();

  #if MPI
  if (advance) run_k_func (dev->s,S,dev->sync,dev->alwaysRun);
  #else
  if (advance) run_k_func (dev->s,S);
  #endif

  
}
CREATE_TAIL(k_func,1)

