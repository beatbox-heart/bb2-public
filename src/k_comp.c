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

/* Karpov's run-time compiler of arithmetic expressions */

#include      <stdio.h>
#include      <setjmp.h>
#include      <string.h>
#include      <stdarg.h>

#include "dynamic.h"
#define OWN
#include      "k_.h"
#undef OWN

/* %%%%%%%%% SYMBOL TABLE %%%%%%%%%%%%%% */

/*  M A C R O     number of the first free record in the table  */
#define LFREE     (table->index[0])

#define LINDEX(c)     (table->index[(int)(c-'@')])
#define ARR(n)        (table->array[n-1])
#define ERROR(code,ptr) {err_code=code;err_ptr=ptr;return 0;}
#define RETERR(retc)  return(err_code = retc)
#define ERRJMP(code,ptr) {err_ptr = ptr; longjmp(_env,code);}

/* ---------------------*/
static  int getline2(p_tb table){
	int n = LFREE;
	if (n == 0) {
		ERROR(TABOVR, "insert");
	}
	LFREE = ARR(n).rf;
	ARR(n).rf = 0;
	return(n);
}

/* ---------------------*/
static  int rlsline(p_tb table, int n){
  ARR(n).rf = LFREE;
  LFREE = n;
  return(0);
}

/* ---------------------*/
/* Create new symbol table and return the pointer to it */
p_tb tb_new(void) {
  int i;
  p_tb table = (p_tb) Calloc(1, sizeof(k_table));
  if (table == NULL) ERROR(LACKDM,"tb_new");
  LFREE = 1;
  for(i=1;i<=maxtab;i++) ARR(i).rf=(i==maxtab)?0:(i+1);
  return(table);
}

/* ---------------------*/
/* Find symbol 'name' in 'table' and return number of the record (1,...)
or 0 if not found */
int tb_find(p_tb table, char name[]) {
  int l = (int)LINDEX(*name);
  if (!table) ERROR(UTABAB,"tb_find");
  if (l) do {
    if (strcmp(name, ARR(l).nm) == 0) return(l);
  } while ((l=(int) ARR(l).rf) != 0);
  return(0);
}

/* ========== find address in the table ========= */
/* return # of 1st record containing address, or 0 if not found */
int tb_findaddr(p_tb table, p_vd addr) {
  int i;
  if (!table) ERROR(UTABAB,"tb_findaddr");
  for (i=1;i<=maxtab;i++){
	if(ARR(i).ad==addr) return i;
	}
  return 0;
}

/* ---------------------*/
/* Insert symbol with specified attributes into the symbol table.
  Returns number of the record (1,..) or 0 in the case of error */
static int tb_insert(
	p_tb	table,		/* the table */
	char	name[],		/* the name */
	int		type,		/* its type tp_undf..tp_dbl */
	p_vd	addr,		/* its address */
	int		npar,		/* # of pars for fun-s or 0 */
	int		flagg		/* attribute (fun, virt.instr, var or RO) */
) {
  int new, last;
  if (!table) ERROR(UTABAB,"tb_insert");
  if (!addr) ERROR(ZEROAD,"tb_insert");
  if ((new = getline2(table)) == 0) ERROR(TABOVR,"tb_insert");
  if ((last = LINDEX(*name)) == 0) {
    LINDEX(*name) = new;
  } else {
    for(;;) {
      if (strcmp(name, ARR(last).nm) == 0) ERROR(DUPNAM, "tb_insert");
      if (!ARR(last).rf) break;
      last = ARR(last).rf;
    }
    ARR(last).rf = new;
  }
  strcpy(ARR(new).nm, name);
  if (flagg & ~_mask) ERROR(INCORF,"tb_insert");
  if ((type>4)||(type<0)) ERROR(INCORT,"tb_insert");
  ARR(new).tp = ((unsigned char)type) | ((unsigned char)(flagg&_mask));
  ARR(new).ad = addr;
  ARR(new).np = (unsigned char)(npar & maxarg);
  return(new);
}

int tb_insert_abstract(p_tb table,char name[],int type,p_vd addr,int npar,int flagg) {
  return tb_insert(table,name,type,addr,npar,flagg);
}

int tb_insert_int  (p_tb table,char name[],p_int  addr) {
  return tb_insert(table,name,t_int,(p_vd)addr,0,f_vb);
}

int tb_insert_int_ro  (p_tb table,char name[],p_int  addr) {
  return tb_insert(table,name,t_int,(p_vd)addr,0,f_ro);
}

int tb_insert_real (p_tb table,char name[],p_real addr) {
  return tb_insert(table,name,t_real,(p_vd)addr,0,f_vb);
}

int tb_insert_real_ro (p_tb table,char name[],p_real addr) {
  return tb_insert(table,name,t_real,(p_vd)addr,0,f_ro);
}

int numerator=0; int denominator=0;
int tb_insert_str  (p_tb table,char name[],char * addr) {
  fprintf(stderr,"tb_insert_str(%s)\n",name);
  numerator=numerator/denominator;
  return tb_insert(table,name,t_str,(p_vd)addr,0,f_vb);
}

int tb_insert_fun  (p_tb table,char name[],p_fn   addr,int npar) {
  return tb_insert(table,name,t_real,(p_vd)addr,npar,f_fn);   
}

/* ---------------------*/
int tb_delete(p_tb table, char name[]) {
  int n, pred;
  if (!table) ERROR(UTABAB,"tb_delete");
  if ((n = LINDEX(*name)) == 0) ERROR(INFCAL,"tb_delete");
  if (strcmp(name, ARR(n).nm) == 0) {
     LINDEX(*name) = ARR(n).rf;
   } else {
        do {
        pred = n;
        if ((n = ARR(n).rf) == 0) ERROR(INFCAL, "tb_delete");
       }
       while (strcmp(name, ARR(n).nm));
       ARR(pred).rf = ARR(n).rf;
     }
  rlsline(table, n);
  return(n);
}

/* ====================== */
/* Globals for the compiler */

#define seg_size (1024*1024)    /* Maximal segmemt size for a virtual code */
jmp_buf _env;              /* ???? */
p_tb    sys_tab = NULL;    /* System (standard) symbol table */
static  p_tb    user_tab;  /* Current user symbol table */
static  char *line_ptr; /* ???? */
static char *lp, *hp;  /* pointers to begin and end of the code */


/* --------------------------------------------- */
/* ============== SCANER ======================= */
/* --------------------------------------------- */

/* ---------- Lexemes' classes ----------------- */
#define cl_eop    1
#define cl_right  2
#define cl_com    3
#define cl_dot    4
#define cl_let    5
#define cl_add    6
#define cl_uadd   7
#define cl_mult   8
#define cl_powr   9
#define cl_left   10
#define idbrk     11
#define ident     12
#define iconst    13
#define fconst    14
#define econst    15
#define user_fn   16
#define sys_fn    17

static  int     cl_lex=cl_eop;

/* ---- lexeme's number in a class  ------------ */
static  int     lex;

/* ----------- classes of characters -------------------- */
#define EOLN    cl_eop
#define LIT     ident
#define DIG     iconst
#define IGNOR   '\040'
#define ERR     0

static  char    litera;
static  char    cl_lit;
static  char    clit[128] ={
EOLN,  IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR,
IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR,
IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR, IGNOR,
ERR, ERR, ERR, LIT, cl_mult, ERR, ERR, cl_left, cl_right, cl_mult, cl_add,
cl_com, cl_add, cl_dot, cl_mult, DIG, DIG, DIG, DIG, DIG, DIG, DIG, DIG, DIG,
DIG, ERR, ERR, ERR, ERR, ERR, ERR, ERR, LIT, LIT, LIT, LIT, LIT, LIT, LIT,
LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT,
LIT, LIT, LIT, LIT, ERR, ERR, ERR, cl_powr, LIT, ERR, LIT, LIT, LIT, LIT, LIT,
LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT, LIT,
LIT, LIT, LIT, LIT, LIT, LIT, ERR, ERR, ERR, ERR, ERR
};

/* ---- scaner's buffer -------------------------- */
#define sb_size 64

static  char    scan_buf[sb_size];
static  char    *sb_ptr;

#define back_space      --line_ptr;

/* ----------------------------*/
static void sp_ignore(void){
  while(clit[*line_ptr & 0x7F] == IGNOR) ++ line_ptr;
  return;
}

/* ----------------------------- */
static void get_ch(void) {
  cl_lit = clit[litera = *line_ptr++ & 0x7F];
  return;
}

static int scaner(void) {

/* ======== START  OF  SCANER ================== */
  sp_ignore();
  err_ptr = (char *)line_ptr;
  get_ch();
  switch(cl_lit) {
  /* ''''''''''''''''''''''''''''' */
  case LIT:
    sb_ptr = scan_buf;
    *sb_ptr++ = litera;
    get_ch();
    while((cl_lit == LIT) || (cl_lit == DIG)) {
      *sb_ptr++ = litera;
      get_ch();
    }
    *sb_ptr++ = '\0';
    if ((sb_ptr - scan_buf) > maxname) * (scan_buf + maxname) ='\0';
      back_space
      sp_ignore();
      get_ch();
    if (cl_lit == cl_left) {
      cl_lex = idbrk;
    } else {
      cl_lex = ident;
      back_space
    }
    lex = tb_find(user_tab, scan_buf);
    break;
  /* ''''''''''''''''''''''''''''' */
  case DIG:
    sb_ptr = scan_buf;
    *sb_ptr++ = litera;
    get_ch();
    while(cl_lit == DIG) {
      *sb_ptr++ = litera;
      get_ch();
    }
    cl_lex = iconst;
    if (cl_lit == cl_dot) {
      cl_lex = fconst;
      *sb_ptr++ = litera;
      get_ch();
      while(cl_lit == DIG) {
        *sb_ptr++ = litera;
        get_ch();
      }
    } /* if cl_dot */
    if (litera == 'e') {
      cl_lex = econst;
      *sb_ptr++ = litera;
      get_ch();
      if (cl_lit == cl_add) {
        *sb_ptr++ = litera;
        get_ch();
      }
      if (cl_lit != DIG) RETERR(INVCON);
      while(cl_lit == DIG) {
        *sb_ptr++ = litera;
        get_ch();
      }
    }   /*    if exp   */
    *sb_ptr++ = '\0';
    back_space
    break;
  /* ''''''''''''''''''''''''''''' */
  case cl_add:
  case cl_mult:
    cl_lex = (int) cl_lit;
    lex = (int) litera;
    if ((litera == '*') && (*line_ptr == '*')){
      get_ch();
      cl_lex = cl_powr;
    }
    break;
  /* ''''''''''''''''''''''''''''' */
  case EOLN:
  case cl_right:
  case cl_com:
  case cl_dot:
  case cl_let:
  case cl_left:
  case cl_powr:
    cl_lex = (int) cl_lit;
    break;
  /* ''''''''''''''''''''''''''''' */
  case ERR:
    RETERR(INVSYM);         /*    invalid symbol  */
  /* ''''''''''''''''''''''''''''' */
  default:
    RETERR(INVSYM);         /*    invalid symbol  */
  /* ''''''''''''''''''''''''''''' */
  } /* of switch */
  return(0);
} /* of scaner */
/* ========================================================= */

/* make global array ? */

static char *t_name[numtypes] = {
   #define D(type,name) name
   D(t_undf, "undef"),
   D(t_int, "int"),
   D(t_real, "real"),
   D(t_str,  "str")
   #undef D
};

static char *tp_nm(int tc) {
  return(t_name[tc]);
}


#define MAXSTRLEN 8192 		/* also defined in device.h but is different scope */
static char s[MAXSTRLEN];
static char *Sprintf(char *fmt, ...) {
  va_list argptr;
  va_start(argptr, fmt);
  vsnprintf(s, MAXSTRLEN, fmt, argptr);
  va_end(argptr);
  return(s);
}

/* ------------ Code generators  -------------------- */

/*#define CODE(type,tab,item) {*((type *)hp) ++ = (type)addr(tab,item);}*/
/*#define CODE(type,item) {*((type *)hp) ++ = (type)(item);}*/
#define CODE(type,item) {*((type *)hp) = (type)(item);  hp+=sizeof(type);}


static  void gnr_fn(p_tb tab, int fnum) {
  int n;
  Sprintf("call_%d",tb_npar(tab, fnum));
  n = tb_find(sys_tab,s);
  if (n==0) ERRJMP(ILGVI,s);
  CODE(p_vi,tb_addr(sys_tab,n));
  CODE(p_fn,tb_addr(tab,fnum));
  return;
}

static void gnr0(char *name) {
  int n = tb_find(user_tab, name);
  if ((n != 0) && (f_vi&tb_flag(user_tab,n))!=0) {
    CODE(p_vi,tb_addr(user_tab,n));
  } else {
    n = tb_find(sys_tab, name);
    if ((n==0) || (f_vi&tb_flag(sys_tab,n))==0) ERRJMP(ILGVI,name);
    CODE(p_vi,tb_addr(sys_tab,n));
  }  /*    of else  */
  return;
}

static void gnr1 (char *name, int par) {
  int n = tb_find(user_tab, name);
  if ((n != 0) && f_vi&tb_flag(user_tab, n)) {
    CODE(p_vi,tb_addr(user_tab,n));
    CODE(p_real,tb_addr(user_tab,par));
  } else {
    n=tb_find(sys_tab, name);
    if ((n == 0) || (f_vi&tb_flag(sys_tab, n))==0) ERRJMP(ILGVI,name);
    CODE(p_vi,tb_addr(sys_tab,n));
    CODE(p_real,tb_addr(user_tab, par));
  }
  return;
}

static void gnr1_ic(INT ic_val) {
  int n = tb_find(user_tab, "push_ic");
  if ((n != 0) && f_vi&tb_flag(user_tab, n)) {
    CODE(p_vi,tb_addr(user_tab,n));
    CODE(INT,ic_val);
  } else {
    n=tb_find(sys_tab,"push_ic");
    if ((n == 0) || (f_vi&tb_flag(sys_tab, n))==0) ERRJMP(ILGVI,"push_ic");
    CODE(p_vi,tb_addr(sys_tab,n));
    CODE(INT,ic_val);
  }
  return;
}

static void gnr1_rc(REAL rc_val) {
  int n = tb_find(user_tab, "push_rc");
  if ((n != 0) && f_vi&tb_flag(user_tab, n)) {
    CODE(p_vi,tb_addr(user_tab, n));
    CODE(REAL,rc_val);
  } else {
    n=tb_find(sys_tab,"push_rc");
    if ((n == 0) || (f_vi&tb_flag(sys_tab, n))==0) ERRJMP(ILGVI,"push_rc");
    CODE(p_vi,tb_addr(sys_tab,n));
    CODE(REAL,rc_val);
  }
  return;
}

/* =========== STACK ==================*/

#define maxstk 50

typedef struct {
  int *ptr;
  int body[maxstk];
} t_stack;

static  t_stack stk_opc;
static  t_stack stk_type;

/* ---------------------------- */
void clear(t_stack *stack) {
  stack->ptr=stack->body+maxstk;
  return;
}

/* ---------------------------- */
int first(t_stack *stack) {
  if (stack->ptr == stack->body+maxstk) ERRJMP(STKUDR,"COMPILE/first");
  return(*stack->ptr);
}

/* ---------------------------- */
int pop(t_stack *stack) {
  if (stack->ptr == stack->body+maxstk) ERRJMP(STKUDR,"COMPILE/pop");
  return(*stack->ptr ++);
}

/* ---------------------------- */
void push(int e, t_stack *stack) {
  if (stack->ptr == stack->body) ERRJMP(STKOVR,"COMPILE");
  *(-- stack->ptr) = e;
  return;
}

/* ------ Stack priorities of operators -------- */
static unsigned char prior[sys_fn]={
  cl_eop, iconst, cl_dot, iconst, cl_let, cl_uadd, cl_uadd,
  cl_powr, cl_powr, cl_right, cl_right, iconst, iconst,
  iconst, iconst, cl_right, cl_right
};

static int st_prior(int opc) {
  return((int)prior[opc-1]);
}

/*------- Table of adjacency ----------------- */
static int adj[econst]={
  077111, 0677, 077110, 077777, 077110, 077010, 077010,
  077010, 077010, 077110, 077110, 0677, 0677, 0677, 0677
};

static int adjac(int *old, int *new) {
  int n=*new;
  int mask=(1<<(n-1));
  int res;
  if (((res=(*old & mask))==0) && (n==cl_add)){
    mask<<=1;
    if ((res=(*old & mask)) !=0) *new=++n;
  }
  *old=adj[--n];
  return(res);
}

static int pr_usfn(int id, int n) {
  /* int tc; */
  if (n > maxarg) RETERR(TOOMAN);
  if (n != tb_npar(user_tab, id)) RETERR(INFCAL); /* Invalid function call */
  /*
  while(n --) {
    tc = pop(& stk_type);
    if (tc != tp_dbl) return(INFCAL);
  }
  */
  push(tb_type(user_tab,id), & stk_type);
  gnr_fn(user_tab, id);
  return(0);
}

static int pr_stfn(int id, int n) {
  /* int tc; */
  if (n > maxarg) RETERR(TOOMAN);
  if (n != tb_npar(sys_tab, id)) RETERR(INFCAL);
  /* while(n --) tc = pop(& stk_type); */
  push(tb_type(sys_tab,id), & stk_type);
  gnr_fn(sys_tab, id);
  return(0);
}

static int pr_opc(int opc) {
  int n, id, tc, tc1, tc2;
  switch(opc) {
  /* --------------------------------- */
  case iconst:
    if (cl_lex==cl_let) RETERR(ASSERR); /* assignment to a constant */
    break;
  /* --------------------------------- */
  case ident:
    id=pop(& stk_opc);
    if (cl_lex == cl_let) gnr1("push_adr", id);
    else gnr1(Sprintf("push_%s",tp_nm(first(&stk_type))),id);
    break;
  /* --------------------------------- */
  case user_fn:
  case sys_fn:
    n=pop(& stk_opc);
    n=pop(& stk_opc);
    /* break not needed ? */
  /* --------------------------------- */
  case cl_left:
    RETERR(MISBRT);                     /* missing right bracket */
  /* --------------------------------- */
  case cl_right:
    n = pop(& stk_opc);
    if (n == user_fn) {         /* call of a user function */
      n = pop(& stk_opc);
      id = pop(& stk_opc);
      return(pr_usfn(id, n));
    }
    if (n == sys_fn) {
      n = pop(& stk_opc);       /* call of a system function */
      id = pop(& stk_opc);
      return(pr_stfn(id, n));
    }
    if (n != cl_left) RETERR(UNSBRT); /* excessive right bracket */
    break;
  /* --------------------------------- */
  case cl_com:
    n=pop(& stk_opc);
    push(1+ pop(& stk_opc), & stk_opc);
    push(n, & stk_opc);
    #if 0
    if (first(& stk_type) != t_real) {
      gnr0( Sprintf( "%s_real",tp_nm(pop(&stk_type)) ) );
      /* push(tp_dbl, & stk_type); */
    }
    #endif
    if ((tc=pop(& stk_type)) != t_real) {
      gnr0( Sprintf( "%s_real",tp_nm(tc) ) );
    }
    break;
  /* --------------------------------- */
  case cl_let:
    tc2=pop(& stk_type);
    tc1=pop(& stk_type);
    if (tc1 != tc2) gnr0(Sprintf("%s_%s",tp_nm(tc1),tp_nm(tc2)));
    gnr0(Sprintf("let_%s",tp_nm(tc1)));
    break;
  /* --------------------------------- */
  case cl_add:
    n=pop(& stk_opc);
    tc2=pop(& stk_type);
    tc1=pop(& stk_type);
    if (tc1 < tc2) {
      tc=tc2;
      gnr0(Sprintf("%s_%s_p",tp_nm(tc1),tp_nm(tc2)));
    } else if (tc1 > tc2) {
      tc=tc1;
      gnr0(Sprintf("%s_%s",tp_nm(tc2),tp_nm(tc1)));
    } else
      tc=tc1;
    gnr0(Sprintf("%s_%s", ((n == '+') ? "add" : "sub"),tp_nm(tc)));
    push(tc, & stk_type);
    break;
  /* --------------------------------- */
  case cl_uadd:
    n=pop(& stk_opc);
    if (n == '-') gnr0(Sprintf("neg_%s",tp_nm(first(& stk_type))));
    break;
  /* --------------------------------- */
  case cl_mult:
    n=pop(& stk_opc);
    tc2=pop(& stk_type);
    tc1=pop(& stk_type);
    if (n == '%') {
      if ((tc1!=t_int)||(tc2!=t_int)) RETERR(INCOMP); /* type mismatch */
      tc=tc1;
      gnr0("mod_int");
    } else {
      if (tc1 < tc2){
        tc=tc2;
        gnr0(Sprintf("%s_%s_p",tp_nm(tc1),tp_nm(tc2)));
      } else if (tc1 > tc2) {
        tc=tc1;
        gnr0(Sprintf("%s_%s",tp_nm(tc2),tp_nm(tc1)));
      } else
        tc=tc1;
      gnr0(Sprintf("%s_%s",((n == '*') ? "mult" : "div"), tp_nm(tc)));
    }
    push(tc, & stk_type);
    break;
  /* --------------------------------- */
   case cl_powr:
     tc2=pop(& stk_type);
     tc1=pop(& stk_type);
     if (tc2 != t_real) gnr0(Sprintf("%s_real",tp_nm(tc2)));
     if (tc1 != t_real) gnr0(Sprintf("%s_real_p",tp_nm(tc1)));
     gnr0("powr_real");
     push(t_real, & stk_type);
     break;
  /* --------------------------------- */
   default:
     RETERR(SYNTAX);
   } /* of switch */
   return(0);
} /* of pr_opc */

/* ============================================= */

void *compile(STRING line, p_tb table, int plant) {
		    void *retval;
  #define bs 1
  int n,tc;
  INT ic_val;
  REAL rc_val;
  /*unsigned!!*/ static int old_elm=0X7E49;
  /* int scaner(); */

  /* --------------------------------------------- */
  /* =========START  OF  COMPILE================== */
  /* --------------------------------------------- */

  if ((err_code = setjmp(_env)) != 0) return(NULL);
  if (sys_tab == NULL) k_init();
  if ((line == (STRING) 0) || (*line == '\0')) ERRJMP(EXPRAB,"COMPILE");
  line_ptr = line;
  err_ptr = (char *)line;
  if (table == (p_tb) 0) ERRJMP(UTABAB,"COMPILE");
  user_tab=table;
  clear(& stk_opc);
  push(bs, & stk_opc);
  push(bs, & stk_opc);
  clear(& stk_type);
  if (NULL==(lp=hp=Calloc(seg_size,1))) ERRJMP(LACKDM,"COMPILE");
  CODE(INT, t_undf);
  err_code=0;
  do {
    if (scaner()) break;
                                        /* invalid order */
    if (!adjac(&old_elm,&cl_lex)) {err_code = SYNTAX;break;}
    switch(cl_lex) {
    case  cl_eop:
      while(st_prior(first(& stk_opc)) > cl_lex)
        if (pr_opc(pop(& stk_opc))) break;
      if (err_code)
        return(NULL);
      if ((plant != t_undf) && (first(& stk_type) != plant)) {
        gnr0(Sprintf("%s_%s",tp_nm(pop(& stk_type)), tp_nm(plant)));
        push(plant, & stk_type);
      }
      *((int *)lp) = pop(& stk_type);
      gnr0("_stop");
      /*return(Realloc((void *)lp,(size_t)((STRING)hp - (STRING)lp)));*/
      retval=(void *)Realloc((void *)lp,(size_t)((STRING)hp - (STRING)lp));
      return(retval);
    case  cl_right:
    case  cl_com:
    case  cl_let:
    case  cl_powr:
    case  cl_left:
      while(st_prior(first(& stk_opc)) > cl_lex)
        if (pr_opc(pop(& stk_opc))) break;
      push(cl_lex,& stk_opc);
      break;
    case  cl_add:
    case  cl_uadd:
    case  cl_mult:
      while(st_prior(first(& stk_opc)) > cl_lex)
        if (pr_opc(pop(& stk_opc))) break;
      push(lex,& stk_opc);
      push(cl_lex,& stk_opc);
      break;
    case  ident:
      if (lex == 0) {
        err_code = ILGVB;         /*  unknown variable */
      } else if ((tc=tb_type(user_tab, lex)) == t_undf) {
        err_code=UNDFTP;          /*  unknown variable's type */
      } else {
        push(tc,& stk_type);
        push(lex,& stk_opc);
        push(ident,& stk_opc);
      }
      break;
    case  idbrk:
      if (lex && f_fn&tb_flag(user_tab, lex)) {
        if(tb_type(user_tab, lex) == t_undf) {
          err_code =UNDFTP;         /* unknown function's type */
        } else {
          push(lex,& stk_opc);
          push(0,& stk_opc);
          push(user_fn,& stk_opc);
          push(cl_com,& stk_opc);
        }
      } else {
        n=tb_find(sys_tab, scan_buf);
	if (!n || (f_fn&tb_flag(sys_tab, n))==0) {
          err_code =ILGFN;         /*  unknown function */
        } else {
          push(n,& stk_opc);
          push(0,& stk_opc);
          push(sys_fn,& stk_opc);
          push(cl_com,& stk_opc);
        }
      }
      break;
    case  iconst:
      sscanf(scan_buf, DFMT , &ic_val);
      gnr1_ic(ic_val);
      push(iconst,& stk_opc);
      push(t_int,& stk_type);
      break;
    case  fconst:
      sscanf(scan_buf, FFMT, & rc_val);
      gnr1_rc(rc_val);
      push(iconst,& stk_opc);
      push(t_real,& stk_type);
      break;
    case  econst:
      sscanf(scan_buf, EFMT, & rc_val);
      gnr1_rc(rc_val);
      push(iconst,& stk_opc);
      push(t_real,& stk_type);
      break;
    default:
      err_code = SYNTAX; /* syntax error */
    } /* of switch */
  } while(! err_code);
  FREE(lp);
  return(NULL);
} /* of compile */

/* -!-!-!-!-!-!-!-!-!-!-!-!-!-!- */

static char format[80], msg[1024];
char *prterr(char *s) {    /* form the error message */
  if ((err_code <= execerr) || (err_code > compierr)){
    snprintf(msg,1024,"%s:\n ", err_ptr);
	if (s) {strcat(msg,"\""); strcat(msg,s); strcat(msg, "\" - ");}
  } else if (s) {
    snprintf(format,80,"\n%%s\n%%%ldc\n-", (err_ptr - s + 1));
    snprintf(msg,1024,format,s,'^');
  }
  strcat(msg,err_msg[err_code]);
  strcat(msg,"\n");
  return(&(msg[0]));
}

char *defaultfmt[numtypes]={"",INTFMT,REALFMT,"%s"};
char *typenametable[numtypes]={"none","int","real","string"};

extern FILE *debug;

/* "Print" a k variable, i.e. convert its value to a text string .	   */
/* NB for numerical variables, this returns the address of a static array, */
/* so should not be called again * until the previous result is utilized!  */
char *prt(void *var, int type) {  /* write the variable of type type */
  static char s[80];
				/*if (debug) fprintf(debug,"\nprt %p %d->",var,type);*/
  if (!var) {
    snprintf(s,80,"(null)");    
  } else switch (type) {
    case t_int: snprintf(s,80,INTFMT,*(INT *)var); break;
    case t_real: snprintf(s,80,REALFMT,*(REAL *)var); break;
    case t_str: return ((char *)var);
    default: *s=0; break;
  }
				/*if (debug) fprintf(debug,"%s\n",s);*/
  return &(s[0]);
}



