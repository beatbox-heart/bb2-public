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

/* Generic definitions related to Karpov's compiler */

#ifndef __K_H
#define __K_H

/* Pointers to the interpreter's data types */
#ifndef INT
#define INT long
#define DFMT "%ld"
#define INTFMT DFMT
#endif


#ifndef REAL
#define REAL double
#define FFMT "%lf"
#define EFMT "%le"
#define REALFMT "%lg"
#endif

typedef INT  *p_int;
typedef REAL *p_real;
typedef int  (*p_vi)(void);
typedef REAL (*p_fn)();
typedef p_fn *pp_fn;

/* define p_vb as a pointer to variable of maximal length */
#ifdef BORLANDC
  #if ( sizeof(p_int) > sizeof(p_real) )
    typedef p_int p_vb;
  #else
    typedef p_real p_vb;
  #endif
#else
  typedef void *p_vb;
#endif

/* define p_vd as a pointer to anything of maximal length */
#ifdef BORLANDC
  #if sizeof(p_fn) >= sizeof(p_vb)
    typedef p_fn p_vd;
  #else
    typedef p_vb p_vd;
  #endif
#else
  typedef p_vb p_vd;
#endif

typedef char *STRING;

/*========= error messages' codes =============*/
enum {
  #define _(tag,msg) tag,
  #include "k_errors.h"
  #undef _
  maxerr
};

/*------------------------------*/

#define execerr ILGVI       /* last execution message  */
#define compierr INCOMP     /* last compilation message  */
/* #define maxerr 30          total # of messages      */

extern   int    err_code;  	/*      error code          */
extern   char *err_ptr;  	/*      error pointer       */
extern   char *err_msg[maxerr];	/*     error messages  */

/*------------------------------*/

/* %%%%%%%%%%%%%% S Y M B O L   T A B L E S  %%%%%%%%%%%%%%%% */

#define maxarg   15        /*  max # of function arguments     */
#define maxname  32        /*  max identifier lentgh */
#define maxtab 1024        /*  max symbol table length */

/* ========= Codes of the compiler's objects types ================== */

typedef enum {
  t_undf, 		/* undefined */
  t_int, 		/* integer */
  t_real, 		/* floating point */
  t_str, 		/* string */
  numtypes
} k_type;

extern char *defaultfmt[numtypes];
extern char *typenametable[numtypes];


/*    Record of a symbol table */

typedef struct{
   char nm[maxname];     /* name of the fun/var */
   unsigned char  tp;    /* type code */
   unsigned char  np;    /* # of parameters, for functions */
   p_vd           ad;    /* address */
   int  rf;    /* # of the next record for a name starting */
			 /* 	with the same letter */
} t_ln;

typedef struct{
   int index['z'-'@'+1];
   /* index[c] = # of the first record for a name */
   /* starting with "c," or 0 if none */
   /* index[0] = # of the first free record in the table */
   t_ln array[maxtab];  /*  the array of records */
} k_table;

typedef k_table *p_tb;

/*==============================*/

/* ====== create a new symbol table ========= */
/* return pointer to the created table or NULL on failure */
p_tb    tb_new(void);

/* ========== find name in the table ========= */
/* return # of record containing name, or 0 if not found */
int tb_find(p_tb table, char name[]);

/* ========== find address in the table ========= */
/* return # of record containing address, or 0 if not found */
int tb_findaddr(p_tb table, p_vd addr);

/* ----------------- */

/* insert the object in the table and return its # or 0 on failure */
int tb_insert_int  (p_tb table,char name[],p_int  addr);
int tb_insert_int_ro  (p_tb table,char name[],p_int  addr);
int tb_insert_real (p_tb table,char name[],p_real addr);
int tb_insert_real_ro (p_tb table,char name[],p_real addr);
int tb_insert_str  (p_tb table,char name[],char * addr);
int tb_insert_fun  (p_tb table,char name[],p_fn   addr,int npar);
int tb_insert_abstract(p_tb table,char name[],int type,p_vd addr,int npar,int flagg);

/* ----------------- */
/* */
int tb_delete(p_tb table, char name[]);

#define tb_type(tb,n) ((tb)->array[(n)-1].tp & ~_mask)
#define tb_flag(tb,n) ((tb)->array[(n)-1].tp &  _mask)
#define tb_npar(tb,n) ((tb)->array[(n)-1].np)
#define tb_addr(tb,n) ((tb)->array[(n)-1].ad)
#define tb_name(tb,n) ((tb)->array[(n)-1].nm)
#define var_name(tb,a,n) ((n=tb_findaddr(tb,a))?tb_name(tb,n):"(unknown)")

/* ========== Initialization of the system ============== */
int  k_init(void);
void k_term(void);

void k_on(void);
void k_off(void);


/* ----------------- */
void  *compile(                 /*  Compile the expression */
	      STRING line,     /*  The expression */
	      p_tb table,      /*  User's symbol table */
	      int plant        /*  Planned expression type */
);
/* Compiles the expression into a sequence of virtual instructions
   and returns the pointer to the beginning of this sequence
   (dynamically allocated, to be freed by the user after use).
   The type of the result is inserted into the 1st two bytes
   of the virtual code. It is identical to 'plant' if the latter
   is not tp_undf.
   In case of error in the expression, err_ptr points to the
   erroneous place in the expression, and err_code is set to a nonzero
   value.
*/

/* -------- Calculate the expression ----------------- */

void *execute(
              void *v_ptr   /* pointer to the virtual code */
);
/* Interprets the sequence of virtual instructions, generated by 'compile'.
   'v_ptr' points to the beginning of this sequence.
   This function returns pointer to the result of calculations
   (of whatever type it is). 
   In case of error, 'err_code' is set to a non-zero, and err_ptr
   points to some additional information about the error
*/

/* ----------------- */
int res_type(
              void *v_ptr   /* pointer to the virtual code */
);
/* Returns the code of the type of the expression,
   contained in the first two bytes of the virtual code */

/* Form the error message
 based on current err_code and err_ptr; s is assumed to point to
 the beginning of the erroneous expression.
 Returns pointer to the beginning of the error message (static object). */

char *prterr(char *);

/* Form the standard representation of the variable var of type type.
   Returns pointer to the beginning of the message (static object). */
char *prt(void *var, int type);

#include "extern.h"
#define CHK(s) if (err_code) EXPECTED_ERROR(prterr(s))

/* -------- Attributes' bits  -------------------- */
#if 1
/* 8< mpi_attribute_bits */
  EXTERN int f_vb, f_fn, f_vi, f_rs, f_ro, _mask;
#ifdef OWN
  int f_vb=0x10;       /* usual variable       */
  int f_fn=0x20;       /* function             */
  int f_vi=0x40;       /* virtual instruction  */
  int f_rs=0x80;       /* reserved attribute   */

  int f_ro=0xF0;  	/* read-only */
  int _mask=0x1F0;
#endif
/* >8 mpi_attribute_bits */
#endif 
EXTERN p_tb sys_tab;
EXTERN p_tb deftb;
EXTERN p_tb loctb;
EXTERN size_t sizetable[numtypes]
#ifdef OWN
  ={0,sizeof(INT),sizeof(REAL)}
#endif
    ;

/* check array of pointers **data of length n for	*/
/* whether it contains address *addr			*/
int used (double **data, int n, double *addr);

/* Check of given expression is just one variable of the needed type */
/* from given symtable, and return its address */
p_vb is_variable (char *expr, p_tb tab, k_type type);

/* check if given (uncompiled) expression depends 	*/
/* on variable varname 					*/
int k_expr_depends (char *expr, char *varname);

#endif

