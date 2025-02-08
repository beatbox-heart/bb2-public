/**
 * Copyright (C) (2010-2018) Vadim Biktashev, Irina Biktasheva et al. 
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

/* Execute code produced by Karpov's compiler */

#include <math.h>

#ifdef __APPLE__
#include <limits.h>
#include <float.h>
#define MINDOUBLE __DBL_MIN__
#define MAXDOUBLE __DBL_MAX__
#else
#include <values.h>
#endif

#include <signal.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "system.h"
#include "dynamic.h"
#include "k_.h"


int err_code;
char *err_ptr;

/* Error messages */
char *err_msg[maxerr]={
  #define _(tag,msg) msg,
  #include "k_errors.h"
  #undef _
};

extern p_tb sys_tab;

extern jmp_buf _env;

static char *vp;

#define stk_size 4096

static  char lostk[stk_size];
static  char *sp;
#define CLEAR()  sp=lostk+stk_size
#define PUSH(T,E)  {sp-=sizeof(*(T)0); *((T)sp)=E;}
#define SPUSH(T,E) {sp-=sizeof(*(T)0); *((T)sp)=E; if(sp<=lostk) return(STKOVR);}
#define POP(T,A) if(sp>=lostk+stk_size) return(STKUDR); else {A=*((T)sp); sp+=sizeof(*(T)0);}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  push_ic(void){
  err_ptr="PUSH";
  SPUSH(p_int, *((INT *)vp));
  vp+=sizeof(INT);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  push_rc(void) {
  err_ptr="PUSH";
  SPUSH(p_real, *((REAL *)vp));
  vp+=sizeof(REAL);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  push_int(void){
  err_ptr="PUSH";
  {
    size_t len=sizeof(long int);  
    long int **VP;
    long int *p;
    long int *SP;
    sp-=len; 
    SP=(long int *)sp;
    VP=(long int **)vp;
    p=*VP;
    *SP = *p;
    if(sp<=lostk) return(STKOVR);
  }
  vp+=sizeof(p_int);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  push_real(void){
  err_ptr="PUSH";
  SPUSH(p_real, *(*((p_real *)vp)));
  vp+=sizeof(p_real);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  push_adr(void){
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  add_int(void){
  INT t1,t2;
  err_ptr = "ADD";
  POP(p_int,t2);
  POP(p_int,t1);
  PUSH(p_int, t1+t2);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  add_real(void){
  REAL t1,t2;
  err_ptr = "ADD";
  POP(p_real,t2);
  POP(p_real,t1);
  PUSH(p_real, t1+t2);
  return(err_code);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  neg_int(void){
  INT t;
  err_ptr = "NEG";
  POP(p_int,t);
  PUSH(p_int, -t);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  neg_real(void){
  REAL t;
  err_ptr = "NEG";
  POP(p_real,t);
  PUSH(p_real, -t);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  sub_int(void){
  INT t1, t2;
  err_ptr = "SUB";
  POP(p_int,t2);
  POP(p_int,t1);
  PUSH(p_int, t1-t2);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  sub_real(void){
  REAL t1, t2;
  err_ptr = "SUB";
  POP(p_real,t2);
  POP(p_real,t1);
  PUSH(p_real, t1-t2);
  return(err_code);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  mult_int(void){
  INT t1,t2;
  err_ptr = "MULT";
  POP(p_int,t2);
  POP(p_int,t1);
  PUSH(p_int, t1*t2);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  mult_real(void){
  REAL t1,t2;
  err_ptr = "MULT";
  POP(p_real,t2);
  POP(p_real,t1);
  PUSH(p_real, t1*t2);
  return(err_code);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  mod_int(void){
  INT t1, t2;
  err_ptr = "MOD";
  POP(p_int,t2);
  if (!t2) return(DIVERR);
  POP(p_int,t1);
  PUSH(p_int, t1%t2);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  div_int(void){
  INT t1, t2;
  err_ptr = "DIV";
  POP(p_int,t2);
  if (!t2) return(DIVERR);
  POP(p_int,t1);
  PUSH(p_int, t1/t2);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  div_real(void){
  REAL t1, t2;
  err_ptr = "DIV";
  POP(p_real,t2);
  if (!t2)  return(DIVERR);
  POP(p_real,t1);
  PUSH(p_real, t1/t2);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  powr_real(void){
  REAL t1, t2;
  err_ptr = "POWR";
  POP(p_real,t2);
  POP(p_real,t1);
  PUSH(p_real,pow(t1, t2));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  int_real(void){
  INT t;
 
  /* err_ptr = "TYPECNF"; */
  POP(p_int,t);
  err_ptr="PUSH";
  SPUSH(p_real, ((REAL) t));
  
  return(0);
}


/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  int_real_p(void){
  REAL t1;
  INT t2;
  err_ptr = "TYPECNV";
  POP(p_real,t1);
  POP(p_int,t2);
  PUSH(p_real, ((REAL) t2));
  SPUSH(p_real, t1);
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  real_int(void){
  REAL t;
  err_ptr = "TYPECNV";
  POP(p_real,t);
  SPUSH(p_int, ((INT) t));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  _stop(void){
  return(-1);
}

static inline int  call_0(void){
  typedef REAL (* p_f0)(void);
  p_f0 func;
  err_ptr = "CALL";
  /*func = *((p_f0 *) vp) ++;*/
  func = *((p_f0 *) vp);
  vp+=sizeof(p_f0);
  SPUSH(p_real,(REAL) func());
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_1(void){
  REAL a1;
  typedef REAL (* p_f1)(REAL);
  p_f1 func;
  err_ptr = "CALL";
  POP(p_real,a1);
  /*func = *((p_f1 *) vp) ++;*/
  func = *((p_f1 *) vp);
  vp+=sizeof(p_f1);
  PUSH(p_real, (REAL) func(a1));
  return(err_code);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_2(void){
  REAL a1, a2;
  typedef REAL (* p_f2)(REAL,REAL);
  p_f2 func;
  err_ptr = "CALL";
  POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f2 *) vp) ++;*/
  func = *((p_f2 *) vp);
  vp+=sizeof(p_f2);
  PUSH(p_real, (REAL) func(a1,a2));
  return(err_code);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_3(void){
  REAL a1, a2, a3;
  typedef REAL (* p_f3)(REAL,REAL,REAL);
  p_f3 func;
  err_ptr = "CALL";
  POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f3 *) vp) ++;*/
  func = *((p_f3 *) vp);
  vp+=sizeof(p_f3);
  PUSH(p_real, (REAL) func(a1,a2,a3));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_4(void){
  REAL a1, a2, a3, a4;
  typedef REAL (* p_f4)(REAL,REAL,REAL,REAL);
  p_f4 func;
  err_ptr = "CALL";
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f4 *) vp) ++;*/
  func = *((p_f4 *) vp);
  vp+=sizeof(p_f4);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_5(void){
  REAL a1,a2,a3,a4,a5;
  typedef REAL (* p_f5)(REAL,REAL,REAL,REAL,REAL);
  p_f5 func;
  err_ptr = "CALL";
  POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f5 *) vp) ++;*/
  func = *((p_f5 *) vp);
  vp+=sizeof(p_f5);
  PUSH(p_real, (REAL)  func(a1,a2,a3,a4,a5));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_6(void){
  REAL a1,a2,a3,a4,a5,a6;
  typedef REAL (* p_f6)(REAL,REAL,REAL,REAL,REAL,REAL);
  p_f6 func;
  err_ptr = "CALL";
  POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f6 *) vp) ++;*/
  func = *((p_f6 *) vp);
  vp+=sizeof(p_f6);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_7(void){
  REAL a1,a2,a3,a4,a5,a6,a7;
  typedef REAL (* p_f7)(REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f7 func;
  err_ptr = "CALL";
  POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f7 *) vp) ++;*/
  func = *((p_f7 *) vp);
  vp+=sizeof(p_f7);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_8(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8;
  typedef REAL (* p_f8)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f8 func;
  err_ptr = "CALL";
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f8 *) vp) ++;*/
  func = *((p_f8 *) vp);
  vp+=sizeof(p_f8);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_9(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9;
  typedef REAL (* p_f9)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f9 func;
  err_ptr = "CALL";
  POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f9 *) vp) ++;*/
  func = *((p_f9 *) vp);
  vp+=sizeof(p_f9);
  PUSH(p_real, (REAL)  func(a1,a2,a3,a4,a5,a6,a7,a8,a9));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_10(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9,a10;
  typedef REAL (* p_f10)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f10 func;
  err_ptr = "CALL";
  POP(p_real,a10);POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f10 *) vp) ++;*/
  func = *((p_f10 *) vp);
  vp+=sizeof(p_f10);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_11(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11;
  typedef REAL (* p_f11)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f11 func;
  err_ptr = "CALL";
  POP(p_real,a11);POP(p_real,a10);POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f11 *) vp) ++;*/
  func = *((p_f11 *) vp);
  vp+=sizeof(p_f11);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_12(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12;
  typedef REAL (* p_f12)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f12 func;
  err_ptr = "CALL";
  POP(p_real,a12);POP(p_real,a11);POP(p_real,a10);POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f12 *) vp) ++;*/
  func = *((p_f12 *) vp);
  vp+=sizeof(p_f12);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_13(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13;
  typedef REAL (* p_f13)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f13 func;
  err_ptr = "CALL";
  POP(p_real,a13);
  POP(p_real,a12);POP(p_real,a11);POP(p_real,a10);POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f13 *) vp) ++;*/
  func = *((p_f13 *) vp);
  vp+=sizeof(p_f13);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_14(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14;
  typedef REAL (* p_f14)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f14 func;
  err_ptr = "CALL";
  POP(p_real,a14);POP(p_real,a13);
  POP(p_real,a12);POP(p_real,a11);POP(p_real,a10);POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f14 *) vp) ++;*/
  func = *((p_f14 *) vp);
  vp+=sizeof(p_f14);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14));
  return(0);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
static inline int  call_15(void){
  REAL a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
  typedef REAL (* p_f15)(REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL,REAL);
  p_f15 func;
  err_ptr = "CALL";
  POP(p_real,a15);POP(p_real,a14);POP(p_real,a13);
  POP(p_real,a12);POP(p_real,a11);POP(p_real,a10);POP(p_real,a9);
  POP(p_real,a8);POP(p_real,a7);POP(p_real,a6);POP(p_real,a5);
  POP(p_real,a4);POP(p_real,a3);POP(p_real,a2);POP(p_real,a1);
  /*func = *((p_f15 *) vp) ++;*/
  func = *((p_f15 *) vp);
  vp+=sizeof(p_f15);
  PUSH(p_real, (REAL) func(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15));
  return(0);
}

/* for explicit type cast */
/* INT Int(REAL a) { return (INT)rint((REAL)a);} */
/* quick fix: linux does not like the orthodoxal version */
REAL Int(REAL a) {return rint(a);}
/* REAL Real(REAL a) {return (REAL) a;} - silly function; bring it back when sure it's needed */

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
int k_init(void) {
  if (NULL==(sys_tab=(p_tb)tb_new())) return 0;
  #define VI(name) if (!tb_insert_abstract(sys_tab,#name,t_undf,(p_vd)name,0,f_vi)) return 0
  VI(push_ic);VI(push_rc);
  VI(push_int);VI(push_real);VI(push_adr);
  VI(add_int);VI(add_real);
  VI(neg_int);VI(neg_real);
  VI(sub_int);VI(sub_real);
  VI(mult_int);VI(mult_real);
  VI(mod_int);VI(div_int);VI(div_real);
  VI(powr_real);VI(int_real);
  VI(int_real_p);VI(real_int);
  VI(_stop);
  VI(call_0);VI(call_1);VI(call_2);VI(call_3);
  VI(call_4);VI(call_5);VI(call_6);VI(call_7);
  VI(call_8);VI(call_9);VI(call_10);VI(call_11);
  VI(call_12);VI(call_13);VI(call_14);VI(call_15);
  #undef VI
  #define FU(name,args) if (!tb_insert_fun(sys_tab,#name,name,args)) return 0
  FU(fabs,1);FU(ceil,1);FU(floor,1);FU(fmod,2);
  FU(sqrt,1);FU(hypot,2);FU(pow,2);
  FU(exp,1);FU(log,1);FU(log10,1);
  FU(sin,1);FU(cos,1);FU(tan,1);
  FU(asin,1);FU(acos,1);FU(atan,1);FU(atan2,2);
  #undef FU
  /* if (!tb_insert_abstract(sys_tab,"int",t_int,(p_vd)Int,1,f_fn)) return 0;   */
  if (!tb_insert_abstract(sys_tab,"int",t_real,(p_vd)Int,1,f_fn)) return 0;  
   /* if (!tb_insert_fun(sys_tab,"real",Real,1)) return 0; - weird, bring it back when sure  */
  return 1;
}

void k_term(void) {
  FREE(sys_tab); 
}

static void fp_trap(int dummy) {
  err_ptr = "EXECUTE";
  longjmp(_env, FPTRAP);
  return;
}

void k_on(void) {
  signal(SIGFPE, fp_trap);
  return;
}

void k_off(void) {
  return;
}


/* ==== THE INTERPRETER ======= */

void *execute(void *v_ptr) {
  p_vi   opc;
  if((err_code=setjmp(_env))!=0) return(NULL);/* for ftp_trap jump */
  vp = (char *)((p_int)v_ptr + 1);            /* skip the result type code */
  CLEAR();				      /* clear the stack */
  while(!err_code) {
    opc = *((p_vi *)vp);
    vp+=sizeof(p_vi);			      /* get next instruction */
    if ((err_code=opc()) >0) return(NULL);    /* execute it */
  }
  err_code = 0;
  return(sp);
}

/*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,*/
int res_type(void *v_ptr){
  return(*(int *) v_ptr);
}


