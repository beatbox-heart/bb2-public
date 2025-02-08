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

/* ROUTINES NEEDED DURING INITIALIZATION */

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "k_.h"
#include "rhs.h"
#include "bgi.h"
#define OWN
#include "qpp.h"
#undef OWN
#include "mpi_io_choice.h"

/*++++++++++++++++++++++++++++++++++++++++*/ extern int idev;


                        /***** Initial values of exported variables *****/
int depth=0;                   /* current depth of include  */

                        /***** Imported variables *****/
extern int narg;               /* defined in main */
extern char **arg;             /* defined in main */
extern int Verbose;            /* defined in main */
extern FILE *res;              /* defined in main */
extern FILE *debug;            /* defined in main */

extern INT Graph;              /* defined in beatbox.c */
extern INT graphon;            /* defined in bgi.c */
extern INT online;             /* defined in bgi.c */
extern INT XMAX;               /* defined in bgi.c */
extern INT YMAX;               /* defined in bgi.c */
extern INT WINX;               /* defined in bgi.c */
extern INT WINY;               /* defined in bgi.c */


/* Variable versions of CGA colours, defined in bgi.h as constants */
#define _(color) static INT var##color=color
_(BLACK);          /* dark colors */
_(BLUE);
_(GREEN);
_(CYAN);
_(RED);
_(MAGENTA);
_(BROWN);
_(LIGHTGRAY);
_(DARKGRAY);           /* light colors */
_(LIGHTBLUE);
_(LIGHTGREEN);
_(LIGHTCYAN);
_(LIGHTRED);
_(LIGHTMAGENTA);
_(YELLOW);
_(WHITE);
#undef _

                        /***** "Local" variables *****/
static union{FILE *f;char *c;}in[MAXDEPTH]; /* input file(s)/string macros */
static long int    inf    = LONG_MAX;       /* value of "infinity" */
static double      pi     = M_PI;           /* pi number */
static double      always = 1.0;            /* use a non-zero zero */
static double      never  = 0.0;            /* this one must be zero */

#undef MARK_ERROR
#define MARK_ERROR return FAILURE;

#define ERR (-2)

/* get char from current input file/string macro, count line and pos */
static int inpgetc(void) {
  int c;
  if (*inname[depth-1]==PASTEBEGIN) { /* Beginning of string macro */
    c = *(in[depth-1].c+(size_t)inppos[depth-1]);
    inppos[depth-1]++;
    if (!c) c=EOF;
  } else {
    c = fgetc(in[depth-1].f);
    if(c=='\n') {inppos[depth-1]=0;inpline[depth-1]++;}
    else {inppos[depth-1]++;}
  }
  return(c);
}

/* Remove C and C++ style comments; nested comments not understood */
static int key=0;               /* finite automaton switch */
static int prevch=0;            /* to remember the character after '/' */
static int uncommgetc(void) {
  int ch, c;
  if (prevch) {c=prevch; prevch=0; return c;}
  c=0;
  while (!c) {
    ch=inpgetc();
    switch(key) {
    case 0:                     /* plain text */
      if (ch=='/') key=1;       /* which comment ? */
      else c=ch;                /* get it */
      break;
    case 1:                     /* which comment ? */
      if (ch=='*') key=2;       /* skip C comment */
      else if (ch=='/') key=4;  /* skip C++ comment */
      else {prevch=ch;c='/';key=0;} /* not a comment */
      break;
    case 2:                     /* skip C comment */
      if (ch=='*') key=3;       /* maybe, end of comment */
      break;
    case 3:                     /* maybe, end of comment */
      if (ch=='/') key=0;       /* yes it is */
      else if (ch=='*') key=3;  /* again, may be */
      else key=2;               /* continue skipping */
      break;
    case 4:                     /* skip C++ comment */
      if (ch=='\n') {key=0;c='\n';}
      break;
    default:
      MESSAGE("ERROR IN UNCOMMGETC!");
      exit(1);
    } /* switch key */
  } /* while !c */
  return c;
}

/*
 * Opens for reading an included <file> or a string [macro]
 */
int qopen (char *name)
{
  if (depth+1>=MAXDEPTH)
    EXPECTED_ERROR("depth %c exceeded in call of %s",MAXDEPTH,name);
  if (*name!=PASTEBEGIN) {
    if NOT(in[depth].f = fopen(name,"rt"))
      EXPECTED_ERROR("cannot open %s",name);
  } else {
    int n=tb_find(deftb,name);
    if (!n) { 
      int argn;
      if (1==sscanf(name,"[%d]",&argn) && argn>=0) {
	MESSAGE("\nError: this script should be given at least %d command-line parameter%s\n",argn,argn>1?"s":""); 
      } else {
	MESSAGE("\nUnknown string macro %s",name); 
      }
      return 0; 
    }
    if (tb_flag(deftb,n)!=f_rs) { 
        MESSAGE("\n%s is not a string macro",name); 
        return 0; 
    }
    in[depth].c = (char *) tb_addr(deftb,n);
  }
  strcpy(inname[depth],name);
  inppos[depth]=0; 
  inpline[depth]=1;
  depth++;
  return 1;
}

int qclose (void)
{
  depth--;
  if (!depth) return 0;
  if(*inname[depth]!=PASTEBEGIN) fclose(in[depth].f);
  return 1;
}

/* Perform includes, syntax <file> or [string] or `syscmd` */
/* 1998/01/26: escape character, e.g. #[ does not open string macro  */
static int qgetc (void)
{
  char *p;
  char name[MAXPATH]; //  = {'\0'};

  char temp1[MAXPATH] = "tmp/bbtp.XXXXXX";
  int errorbbtp_qpp = 0;

  int c=uncommgetc();
  if (strchr(ESCAPE,c)) {
    c=uncommgetc();
    if (c==EOF) {
      if (!qclose()) c=EOF; else c = qgetc();
    }
  } else if (c==INCLUDEBEGIN) {
    p=name;
    while(EOF!=(c=qgetc())) { /* allow include filename to be a string macro */
     if (c==INCLUDEEND) break;
     *(p++)=c;
    }
    *p=0;
    if (!qopen(name)) 
      c=ERR; 
    else 
      c = qgetc();
  } else if (c==CATCHBEGIN) {
    char cmd[MAXSTRLEN];
    p=cmd; // put a break here.
    while(EOF!=(c=uncommgetc())) {
     if (c==CATCHEND) break;
     *(p++)=c;
    }
    *p=0;
    /*
      The strcat may give rise to spurious files with similar names to
      simulation files. This is a simple solution to the use of
      tmpnam, a function that is not valid on all systems,
      e.g. Hector. A simple hack? Another way of dealing with this is
      to use a different char string altogether, if this name is local
      to qgetc. Check before doing.
     */
/* tmpnam(name); */ // original code.

/*  mkstemp(name);
    strcat(name,"t"); */

    system("mkdir -p tmp");
    errorbbtp_qpp = 0;
    errorbbtp_qpp = mkstemp(temp1);


    strcat(cmd,">"); strcat(cmd,temp1);
    system(cmd);

    if (!qopen(temp1))
        c=ERR; 
    else 
        c = qgetc();
    }
    else if (c==PASTEBEGIN){
    // p=name; 
    p = temp1;
    *p++=PASTEBEGIN;
    while(EOF!=(c=uncommgetc())) {
     if (p>=temp1+MAXPATH-1) {
         MESSAGE("too long name %s..",temp1);
         return EOF;
     }
     *(p++)=c;
     if (c==PASTEEND) break;
    }

    *p=0;

//    if (!qopen(name)) c=ERR; else c = qgetc();
      if (!qopen(temp1)) c=ERR; else c = qgetc();
  } else if (c==EOF) {
    if (!qclose()) c=EOF; else c = qgetc();
  }
  return c;
}

/*
 * Extracts a command from the script.
 */
/* Read to next terminator or eof (error). */
/* Return num of bytes read, or 0 if eof (error) is the first char */
int read_command (char *s, int strlen)
{
  int c;
  char *p, *e;
  int blocks=0;
  for(p=s,e=p+strlen-1;p<e;p++) {
    *p=c=qgetc(); // this command reads in the complete line.
    if (c==ERR) {
       MESSAGE("\nError occurred in input file \"%s\" at line %ld pos %ld",
               inname[depth-1], inpline[depth-1], inppos[depth-1]); 
               return 0;
    }
    if (c==EOF) { 
        MESSAGE("unexpected end of file"); 
        return 0;
    }
    if (c==BLOCKBEGIN) blocks++;
    if (c==BLOCKEND) blocks--;
    if (!blocks && strchr(TERMINATOR,*p)) break;
  }
  *p=0;
  if (blocks) EXPECTED_ERROR("block not closed"); /* Can this ever happen? */
  return((int)(p-s));
}

int calc_with_table (void *resaddr, int restype, char *expr, p_tb table)
{
  void *p;
  pp_fn code;
  if (restype==0)                  {err_code=UNDFTP; return(0);}
  if (restype<0 || restype>t_real) {err_code=INCORT; return(0);}
  code = compile (expr, table, restype); CHK(expr);
  p = execute(code);                     CHK(NULL);
  FREE (code);
  memcpy(resaddr,p,sizetable[restype]);
  return(1);
}

/*
 * Evaluates a script expression
 */
int calc (void *resaddr, int restype, char *expr)
{
  return calc_with_table (resaddr, restype, expr, deftb);
}

/*
 * Check if a given name is a valid name for a k_variable
 */
int valid_name (char *s)
{
  char *p;
  for (p=s;*p;p++) {
    if (p==s) {
      if ((0==isalpha(*p)) && (*p!='_')) return 0;
    } else {
      if ((0==isalnum(*p)) && (*p!='_')) return 0;
    }
  }
  return 1;
}


/*
 * Defines a k_variable
 */
/* TODO: this function allocates memory (newvar)             */
/* which is never freed, so technically it is a memory leak. */
int def (char *s) {
  char *typename, *parname, *expr, tabname[MAXPATH+2];
  p_vb newvar=NULL;
  int restype, flag;
  typename = strtok(s, SEPARATORS);
  parname  = strtok(NULL, SEPARATORS"=");
  expr     = strtok(NULL, "" /*SEPARATORS*/);
  flag=f_vb|f_rs;
  if      (0==strcmp(typename,"int"))    restype=t_int;
  else if (0==strcmp(typename,"long"))   restype=t_int;
  else if (0==strcmp(typename,"float"))  restype=t_real;
  else if (0==strcmp(typename,"double")) restype=t_real;
  else if (0==strcmp(typename,"real"))   restype=t_real;
  else if (0==strcmp(typename,"str"))    {restype=t_undf; flag=f_rs;}
  else EXPECTED_ERROR("unknown type %s",typename);

  if (flag==f_rs)
    snprintf(tabname,MAXPATH+2,"%c%s%c",PASTEBEGIN,parname,PASTEEND);
  else {
    if (0==valid_name(parname))		/* is it a valid identifier name? */
      EXPECTED_ERROR("%s is not a valid name\n",parname);
    strcpy(tabname,parname);
  }   	

  if(tb_find(sys_tab,tabname))          /* is it in system table ? */
    EXPECTED_ERROR("Name %s is reserved",tabname);

  if (tb_find(deftb,tabname))           /* is it in user table ? */
    EXPECTED_ERROR("%s already declared\n",tabname);

  if(flag&f_vb) {                       /* arithmetic variable */
    /* CALLOC(newvar,1,sizetable[restype]); */
    if (0==(newvar=Calloc(1,sizetable[restype])))
      ABORT("could not allocate %ld bytes for %s",((long)(sizetable[restype])),"newvar");
    if (!tb_insert_abstract(deftb, tabname, restype, newvar, 0, flag))
      EXPECTED_ERROR("cannot declare %s:\n%s",tabname,prterr(NULL));
    if (expr==NULL && restype != t_str) 
      expr="0";
    if (!calc(newvar,restype,expr))
      EXPECTED_ERROR("cannot interpret %s",expr);
    if (Verbose)
      MESSAGE("\x01""\ndef %s %s %s$",typename,parname,prt(newvar,restype));
  } else if (flag==f_rs) {              /* string for pasting */
    if (NULL==expr) newvar=strdup("");
    else if (NULL==(newvar=(p_vb)strdup(expr)))
      EXPECTED_ERROR("not enough memory for %s",tabname);
    if (!tb_insert_abstract(deftb, tabname, restype, newvar, 0, flag))
      EXPECTED_ERROR("cannot declare %s:\n%s",tabname,prterr(NULL));
    if (Verbose)
      MESSAGE("\x01""\ndef %s %s %s$",typename,parname,newvar);
  } else {
    EXPECTED_ERROR("invalid flag %d",flag);
  }
  return 1;
}

/* TODO: this function allocates memory (newvar)             */
/* which is never freed, so technically it is a memory leak. */
int def_local(char *s, p_tb table){     
  char *typename, *parname, *expr, tabname[MAXPATH+2];
  p_vb newvar=NULL;
  int restype, flag;
  typename = strtok(s, SEPARATORS);
  parname  = strtok(NULL, SEPARATORS"=");
  expr     = strtok(NULL, "" /*SEPARATORS*/);
  flag=f_vb|f_rs;
  if      (0==strcmp(typename,"int"))    restype=t_int;
  else if (0==strcmp(typename,"long"))   restype=t_int;
  else if (0==strcmp(typename,"float"))  restype=t_real;
  else if (0==strcmp(typename,"double")) restype=t_real;
  else if (0==strcmp(typename,"real"))   restype=t_real;
  else if (0==strcmp(typename,"str"))   EXPECTED_ERROR("String macros cannot be defined as local variables.\n")
  else EXPECTED_ERROR("unknown type %s",typename)
  
  strcpy(tabname,parname);
  
  if(tb_find(sys_tab,tabname))          /* is it in system table ? */
    EXPECTED_ERROR("Name %s is reserved",tabname);
  
  if (tb_find(table,tabname))		/* is it in given user table ? */
    EXPECTED_ERROR("%s already declared\n",tabname);
  
  if(flag&f_vb) {			/* arithmetic variable */
    /* CALLOC(newvar,1,sizetable[restype]); */
    if (0==(newvar=Calloc(1,sizetable[restype])))
      ABORT("could not allocate %ld bytes for %s",((long)(sizetable[restype])),"newvar");
    if (!tb_insert_abstract(table, tabname, restype, newvar, 0, flag))
      EXPECTED_ERROR("cannot declare %s:\n%s",tabname,prterr(NULL));
    if (expr==NULL && restype != t_str) 
      expr="0";
    if (!calc_with_table(newvar,restype,expr,table))
      EXPECTED_ERROR("cannot interpret %s",expr);
    if (Verbose)
      MESSAGE("\x01""\n\tdef %s %s=%s;",typename,parname,expr);
  } else {
    EXPECTED_ERROR("invalid flag %d",flag);
  }
  return 1;
}

/*
 * Defines a device's name as a k_variable.
 */
int def_dev (Device *d) {
  char tabname[maxname+2];
  snprintf(tabname,maxname+2,"%c%s%c",PASTEBEGIN,d->n,PASTEEND);
  if(tb_find(sys_tab,tabname))          /* is it in system table ? */
    EXPECTED_ERROR("%s cannot be used as device name: it is reserved",tabname);
  if (tb_find(deftb,tabname))           /* is it in user table ? */
    EXPECTED_ERROR("name %s already in use\n",tabname);
  if (!tb_insert_abstract(deftb, tabname, t_undf, (p_vd)d, 0, f_rs))
    EXPECTED_ERROR("cannot declare %s:\n%s",tabname,prterr(NULL));
  return 1;
}

/*
 * Obtains a reference to a device from its name.
 */
int get_dev (Name n, Device **d) {
  char tabname[maxname+2];
  int i;
  snprintf(tabname,maxname+2,"%c%s%c",PASTEBEGIN,n,PASTEEND);
  i = tb_find(deftb,tabname);
  if (!i)
    EXPECTED_ERROR("Device name %s was never declared\n",tabname);
  if (tb_type(deftb,i)!=t_undf || tb_flag(deftb,i)!=f_rs)
    EXPECTED_ERROR("%s is not a device name",tabname);
  *d = (Device *)tb_addr(deftb,i);
  return 1;
}

/* Find first token in src */
static int token (char **dest, char *src) {
  strcpy(buf,src);
  *dest = strtok(buf,SEPARATORS);
  if (*dest) return 1;
  *dest = src;
  return 0;
}

/*
 * Extracts a keyword from a command.
 *
 * Point *w to first word in s, and return pointer to next after.
 */

char *first_word (char *s, char **w,char *delim) {
  char delimdflt[]=" \n\r\t";
  char *rest;
  char *d;
  if (delim) d = delim;
  else d = delimdflt;
  for(*w=s;**w&&strchr(d,**w);(*w)++);  /* find first good char within s */
  for(rest=*w;;rest++) {                /* find first bad char after */
    if(!*rest) break;
    if(strchr(d,*rest)) break;
  };
  if (*rest) {
    *(rest++)='\0';                     /* mark end of found word */
    return(rest);                       /* rest begins next */
  } else
    return(NULL);                       /* no rest */
}

/* Find key in buffer, return to the next char after key or NULL if not */
char *find_key (const char *key, char *s) {
  char *p=s;
  if (s==NULL) return NULL;
  for(p=s;p<s+strlen(s);p+=strlen(key)) {
    if NOT(p=strstr(p,key)) return NULL;
    if (p==s || strchr(SEPARATORS, p[-1])) return(p+strlen(key));
  }
  return NULL;
}


/**************************************************
Accept a variable:
  ACCEPTI: int
  ACCEPTL: long
  ACCEPTR: real
  ACCEPTS: string
Parameters:
  const char *name,    name of var - key to be found
  int *var,            addr of var
  int deflt,           default value; may be NONE      }  Not
  int minv,            allowed min value; may be NONE  }  for
  int maxv,            allowed max value; may be NONE  }  strings
  const char *w        input string
  where NONE::={INONE,LNONE,RNONE}
Also:
  ACCEPTB: block BEGIN ... END { }
  ACCEPT_CONDITION: the Time structure
  ACCEPT_SPACE: the Space structure
  ACCEPT_WINDOW: the BGIWindow structure
*****************************************************/

/*
 * Accepts an integer parameter.
 */
int accepti (
  const char *name,int *var,int deflt, int minv, int maxv, char *w
) {
  char *p = find_key(name,w);
  INT dummy;
  if (!p && deflt!=INONE) {
    *var = deflt;
    if (Verbose) MESSAGE4("\x01%s%d%s%c",name,*var,DFLT,SEPARATORS[0]);
    return 1;
  } else if (!p && deflt==INONE) {
    MESSAGE2("ERROR: cannot find word \"%s\" in the string:\n%s\n",name,w);
    return 0;

  } else if (!token(&p,p)) {
    MESSAGE2("ERROR parsing %s%s\n",name,p);
    return(0);

  } else if (!calc(&dummy,t_int,p) || (*var=(int)dummy,*var!=dummy)) {
    MESSAGE2("ERROR reading %s%s\n",name,p);
    return(0);
  } else if (minv!=INONE && *var<minv) {
    MESSAGE3("ERROR read value %s%d < allowed min=%d\n",name,*var,minv);
    return(0);
  } else if (maxv!=INONE && *var>maxv) {
    MESSAGE3("ERROR read value %s%d > allowed max=%d\n",name,*var,maxv);
    return(0);
  } else {
    MESSAGE3("\x01%s%d%c",name,*var,SEPARATORS[0]);
    return 1;
  }
}

/*
 * Accepts an entry of an integer parameter array.
 */
int acceptie (
  const char *mask,int *arr,int i,int deflt,int minv,int maxv,char *w
) {
  char name[MAXSTRLEN], *p;
  int *var=arr+i;
  snprintf(name,MAXSTRLEN,mask,i);
  if (0==strcmp(name,mask)) 
    EXPECTED_ERROR("Bad mask '%s' for an enumerated parameter name: substituted %d, obtained '%s'",
		   mask,i,name);
  return accepti(name,var,deflt,minv,maxv,w);
}

/*
 * Accepts a long integer parameter.
 */

int acceptl(
  const char *name,long *var,long deflt,long minv,long maxv,char *w
) {
  INT dummy;
  char *dflt="";
  char *p = find_key(name,w);
  if (!p && deflt!=LNONE) {
    *var = deflt;
    if (!Verbose) return 1;
    dflt = DFLT;
  } else if (!p && deflt==LNONE) {
    MESSAGE2("ERROR: cannot find word \"%s\" in the string:\n%s\n",name,w);
    return 0;
  } else if (!token(&p,p)) {
    MESSAGE2("ERROR parsing %s%s\n",name,p);
    return(0);
  } else if (!calc(&dummy,t_int,p) || (*var=(long)dummy,*var!=dummy)) {
    MESSAGE2("ERROR reading %s%s\n",name,p);
    return(0);
  } else if ((minv!=LNONE) && ((*var)<minv)) {
    MESSAGE3("ERROR read value %s%ld < allowed min=%ld\n",name,*var,minv);
    return(0);
  } else if (maxv!=LNONE && *var>maxv) {
    MESSAGE("ERROR read value %s%ld > allowed max=%ld (LNONE=%ld)\n",name,*var,maxv,LNONE);
    return(0);
  }
  if (*var==inf) {
    MESSAGE3("\x01%s*%s%c",name,dflt,SEPARATORS[0]);
    return 1;
  } else {
    MESSAGE4("\x01%s%ld%s%c",name,*var,dflt,SEPARATORS[0]);
    return 1;
  }
}

/*
 * Accepts a real parameter.
 */
int acceptr (
  const char *name,real *var,real deflt,real minv,real maxv,char *w
) {
  char *p = find_key(name,w);
  if (!p && deflt!=RNONE) {
    *var = deflt;
    if (!Verbose) return 1;
    MESSAGE4("\x01%s" REALF "%s%c",name,*var,DFLT,SEPARATORS[0]);
    return 1;
  } else if (!p && deflt==RNONE) {
    MESSAGE2("ERROR: cannot find word \"%s\" in the string:\n%s\n",name,w);
    return 0;
  } else if (!token(&p,p)) {
    MESSAGE2("ERROR parsing %s%s\n",name,p);
    return(0);
  } else if (!calc(var,t_real,p)) {
    MESSAGE2("ERROR reading %s%s\n",name,p);
    return(0);
  } else if (minv!=RNONE && *var<minv) {
    MESSAGE3("ERROR read value %s%f < allowed min=%f\n",name,*var,minv);
    return(0);
  } else if (maxv!=RNONE && *var>maxv) {
    MESSAGE3("ERROR read value %s%f > allowed max=%f\n",name,*var,maxv);
    return(0);
  } else {
    MESSAGE3("\x01%s" REALF "%c",name,*var,SEPARATORS[0]);
    return 1;
  }
}

/*
 * Accepts an entry of a  real parameter array.
 */
int acceptre (
  const char *mask,real *arr,int i,real deflt,real minv,real maxv,char *w
) {
  char name[MAXSTRLEN], *p;
  real *var=arr+i;
  snprintf(name,MAXSTRLEN,mask,i);
  if (0==strcmp(name,mask)) 
    EXPECTED_ERROR("Bad mask '%s' for an enumerated parameter name: substituted %d, obtained '%s'",
		   mask,i,name);
  return acceptr(name,var,deflt,minv,maxv,w);
}

/*
 * RHS only. 
 * Accepts a 'dependent parameter' which may take reference to a layer
 * instead of a specific value. 
 */
int acceptp(
  const char *name,
  real *var,
  real deflt,
  real minv,
  real maxv,
  char *w,
  Var *v,
  int *iv,
  int v0
) {
  #define SRC (v->src[*iv])
  #define DST (v->dst[*iv])
  INT dummy;
  char *p = find_key(name,w);
  if (p==NULL) {
    MESSAGE("\x01""\n\t");
    return(acceptr(name,var,deflt,minv,maxv,w));
  } else if (!token(&p,p)) {
    MESSAGE("ERROR parsing %s%s\n",name,p);
    return(0);
  } else if NOT(strchr(p,AT)) {
    MESSAGE("\x01""\n\t");
    return(acceptr(name,var,deflt,minv,maxv,w));
  } else if (*iv>=v->n) {
    MESSAGE("# of variable parameters exeeded",v->n);
    return(0);
  } else if (!calc(&dummy,t_int,++p) || (SRC=(int)dummy,SRC!=dummy)) {
    MESSAGE("ERROR reading %s%s\n",name,p);
    return(0);
  } else if (SRC>=vmax) {
    MESSAGE("ERROR read value %s%c%d > allowed max=%c%d\n",name,AT,SRC,AT,vmax-1);
    return(0);
  } else {
    *var=deflt;
    /* DST=var; */
    (v->dst[*iv])=var;
    MESSAGE("\x01""\n\t%s@%d%c",name,SRC,SEPARATORS[0]);
    SRC-=v0;	/* source layers are counted relative to the rhs base layer */
    (*iv)++;
    return 1;
  }
  #undef SRC
  #undef DST
}

/*
 * Accepts a real 'linked parameter' which may follow current value of a k-variable
 * or k-code instead of a fixed value. 
 */
int acceptrk (const char *name, 
	      REAL **varptr, REAL *varval, char **varsrc, pp_fn *varexe, p_tb vartab,
	      real deflt, real minv, real maxv, char *w)
{
  char *p = find_key(name,w);
  if (p==NULL) {
    *varptr=varval;
    *varsrc=NULL;
    *varexe=NULL;
    /* MESSAGE("\x01""\n\t"); */
    return(acceptr(name,varval,deflt,minv,maxv,w));
  } else if (!token(&p,p)) {
    MESSAGE("ERROR parsing %s%s\n",name,p);
    return(0);
  } else if (*p!=LINK) {
    /* MESSAGE("\x01""\n\t"); */
    *varptr=varval;
    *varsrc=NULL;
    *varexe=NULL;
    return(acceptr(name,varval,deflt,minv,maxv,w));
  } else if (*(++p)=='\0') {
    EXPECTED_ERROR("empty k-code for '%s' in '%s'\n", name, w);
  } else if (NULL!=(*varptr=is_variable(p,vartab,t_real))) {
    *varval=**varptr;
    *varexe=NULL;
    *varsrc=NULL;
    MESSAGE("\x01%s%c%s=" REALF "..%c",name,LINK,p,*varval,SEPARATORS[0]);
    return(1);
  } else {
    /* MESSAGE("\x01""\n\t"); */
    *varptr=NULL;
    *varsrc=strdup(p);
    *varexe=compile (p, vartab, t_real); CHK(p);
    *varval=*(REAL *)execute(*varexe); CHK(NULL);
    MESSAGE("\x01%s%c%s=" REALF "..%c",name,LINK,*varsrc,*varval,SEPARATORS[0]);
    return(1);
  }
}

/*
 * Accepts an integer 'linked parameter' which may follow current value of a k-variable
 * or k-code instead of a fixed value. 
 */
int acceptik (const char *name, 
	      INT **varptr, INT *varval, char **varsrc, pp_fn *varexe, p_tb vartab,
	      INT deflt, INT minv, INT maxv, char *w)
{
  char *p = find_key(name,w);
  if (p==NULL) {
    *varptr=varval;
    *varsrc=NULL;
    *varexe=NULL;
    /* MESSAGE("\x01""\n\t"); */
    return(acceptl(name,varval,deflt,minv,maxv,w));
  } else if (!token(&p,p)) {
    MESSAGE("ERROR parsing %s%s\n",name,p);
    return(0);
  } else if (*p!=LINK) {
    /* MESSAGE("\x01""\n\t"); */
    *varptr=varval;
    *varsrc=NULL;
    *varexe=NULL;
    return(acceptl(name,varval,deflt,minv,maxv,w));
  } else if (*(++p)=='\0') {
    EXPECTED_ERROR("empty k-code for '%s' in '%s'\n", name, w);
  } else {
    /* MESSAGE("\x01""\n\t"); */
    *varsrc=strdup(p);
    *varexe=compile (p, vartab, t_int); CHK(p);
    *varval=*(INT *)execute(*varexe); CHK(NULL);
    *varptr=NULL;
    MESSAGE("\x01%s%c%s=%ld..%c",name,LINK,*varsrc,(long)*varval,SEPARATORS[0]);
    return(1);
  }
}



/* Todo: this is unsafe if memory allocated at *var is shorter         */
/* than strings offered by *deflt or *p - need to protect against that */
int accepts (const char *name,char *var,const char *deflt, char *w) {
  char *p = find_key(name,w), *p1;
  if (!p && deflt) {
    strcpy(var,deflt);
    if (!Verbose) return 1;
    if (*deflt) MESSAGE4("\x01%s\"%s\"%s%c",name,var,DFLT,SEPARATORS[0]);
    return 1;
  } else if (!p && !deflt) {
    MESSAGE2("ERROR: cannot find word \"%s\" in the string:\n\"%s\"\n",name,w);
    return 0;
  } else if (*p==STRBEGIN) {
    for (p1=var,p++;*p&&*p!=STREND;p++,p1++) {
      if(*p=='\\') {
        switch(*++p) {
          case 'n': *p1='\n'; break;case 't': *p1='\t'; break;
          case 'b': *p1='\b'; break;case 'r': *p1='\r'; break;
          case 'f': *p1='\f'; break;case 'a': *p1='\a'; break;
          case '0': *p1='\0'; break;default:  *p1=*p;   break;
        };
      } else {
        *p1=*p;
      }
    }
    if (*p!=STREND)
      MESSAGE("/* WARNING: non-terminated string %c..%c in %s */",STRBEGIN,STREND,w);
    *p1='\0';
  } else if (!token(&p,p)) {
    MESSAGE2("ERROR parsing %s%s\n",name,p);
    return(0);
  } else {
    strcpy(var,p);
  }
  MESSAGE3("\x01%s\"%s\"%c",name,var,SEPARATORS[0]);
  return 1;
}

/* This is an attempt of a safer version of accepts */
int acceptsn (const char *name,char *var,int len,const char *deflt, char *w) {
  char *p = find_key(name,w), *p1;
  if (!p && deflt) {
    strncpy(var,deflt,len);
    if (!Verbose) return 1;
    if (*deflt) MESSAGE4("\x01%s\"%s\"%s%c",name,var,DFLT,SEPARATORS[0]);
    return 1;
  } else if (!p && !deflt) {
    MESSAGE2("ERROR: cannot find word \"%s\" in the string:\n\"%s\"\n",name,w);
    return 0;
  } else if (*p==STRBEGIN) {
    for (p1=var,p++;(*p)&&(*p!=STREND)&&(p1<var+len);p++,p1++) {
      if(*p=='\\') {
        switch(*++p) {
          case 'n': *p1='\n'; break;case 't': *p1='\t'; break;
          case 'b': *p1='\b'; break;case 'r': *p1='\r'; break;
          case 'f': *p1='\f'; break;case 'a': *p1='\a'; break;
          case '0': *p1='\0'; break;default:  *p1=*p;   break;
        };
      } else {
        *p1=*p;
      }
    }
    if (*p!=STREND)
      MESSAGE("/* WARNING: non-terminated string %c..%c in %s */",STRBEGIN,STREND,w);
    *p1='\0';
  } else if (!token(&p,p)) {
    MESSAGE2("ERROR parsing %s%s\n",name,p);
    return(0);
  } else {
    strncpy(var,p,len);
  }
  MESSAGE3("\x01%s\"%s\"%c",name,var,SEPARATORS[0]);
  return 1;
}

/*
 * Accepts a file parameter.
 */
int acceptf (const char *name,const char *mode,const char *deflt,char *fname,FILE **f, char *w) {
  /* DEBUG("acceptf enter\n"); */
  if NOT(accepts(name,fname, deflt, w)) {/* DEBUG("not accepts\n");  */return 0;}
  else if (0==strcmp(fname,"")) *f=NULL;
  else if (0==strcmp(fname,"stdout")) *f=stdout;
  else if (0==strcmp(fname,"stderr")) *f=stderr;
  else if (0==strcmp(fname,"log"))    *f=res;
  else if (0==strcmp(fname,"debug"))  *f=debug;
  else {
    if (0==stricmp(fname,null)) strcpy(fname,NULLFILE);
    if NOT(*f=fopen(fname,mode)) EXPECTED_ERROR("cannot open file %s in mode %s", fname, mode);
  }
  /* DEBUG("acceptf exit\n"); */
  return 1;
}

/*
 * Accepts a codestring parameter.
 */
int acceptc (const char *name,const k_type type,const char *deflt,char *codestring,void **code, char *w) {  
  if NOT(accepts(name,codestring, deflt, w)) return 0;
  else if (0==strcmp(codestring,"")) *code=NULL;
  else {*code = compile (codestring, loctb, type); CHK(codestring);}
  return 1;
}

/*
 * Accepts a block parameter.
 */
int acceptb (const char *name,char *var, char *w) {
  char *b, *e;
  if NOT(b = find_key(name,w)) EXPECTED_ERROR("cannot find block \"%s%c..%c\" in:\n %s",name,BLOCKBEGIN,BLOCKEND,w);
  while (*b && strchr(SEPARATORS,*b)) b++;
  if (*b != BLOCKBEGIN) EXPECTED_ERROR("no block \"%c..%c\" after \"%s\" in:\n %s",BLOCKBEGIN,BLOCKEND,name,w);
  e=++b;
  while (*e && *e!=BLOCKEND) e++;
  if (*e != BLOCKEND) {
    MESSAGE("cannot find end of block \"%c\" in %s",BLOCKEND,w);
    return 0;
  }
  strncpy(var,b,(size_t)(e-b));
  var[(size_t)(e-b)]=0;
  MESSAGE("\x01%s%c",name,BLOCKBEGIN);
  /* if(Verbose) MESSAGE("\x01\n"); */
  return 1;
}


void close_block(void) {
  if(Verbose) MESSAGE("\x01""\n");
  MESSAGE2("\x01%c%c",BLOCKEND,SEPARATORS[0]);
}


/* int accept_real_variable (REAL **v, char **vname, char *name, char *w) */
/* { */
/*   typedef struct {char var[80];} STR; */
/*   STR s; */
/*   STR *S=&s; */
/*   int i; */
/*   char key[80]; */
/*   snprintf(key,80,"%s=",name); */
/*   accepts(key,&(S->var[0]),"",w); */
/*   if (!strcmp(S->var,"")) { */
/*      EXPECTED_ERROR("cannot find word \"%s\" in the string:\n\"%s\"\n",name,w); */
/*   } */
/*   /\* k_on(); *\/ */
/*   if NOT(i=tb_find(deftb,S->var))  EXPECTED_ERROR("unknown symbol %s",S->var); */
/*   if (tb_type(deftb,i)!=t_real) EXPECTED_ERROR("symbol %s is not double",S->var); */
/*   if ((tb_flag(deftb,i)&f_ro)==f_ro) EXPECTED_ERROR("symbol '%s' is read-only",S->var); */
/*   if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("symbol %s is not a variable",S->var); */
/*   *v = (REAL *)tb_addr(deftb,i); */
/*   *vname = tb_name(deftb,i); */
/*   return 1; */
/* } */

/* int accept_real_variable_ro (REAL **v, char **vname, char *name, char *w) */
/* { */
/*   typedef struct {char var[80];} STR; */
/*   STR s; */
/*   STR *S=&s; */
/*   int i; */
/*   char key[80]; */
/*   snprintf(key,80,"%s=",name); */
/*   accepts(key,&(S->var[0]),"",w); */
/*   if (!strcmp(S->var,"")) { */
/*      EXPECTED_ERROR("cannot find word \"%s\" in the string:\n\"%s\"\n",name,w); */
/*   } */
/*   /\* k_on(); *\/ */
/*   if NOT(i=tb_find(deftb,S->var))  EXPECTED_ERROR("unknown symbol %s",S->var); */
/*   if (tb_type(deftb,i)!=t_real) EXPECTED_ERROR("symbol %s is not double",S->var); */
/*   if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("symbol %s is not a variable",S->var); */
/*   *v = (REAL *)tb_addr(deftb,i); */
/*   *vname = tb_name(deftb,i); */
/*   return 1; */
/* } */

int acceptkr (const char *key,const int writable,const char *deflt,REAL **varaddr,char *varname,char *w)
{
  int i;
  if NOT(accepts(key,varname,deflt,w)) return FAILURE;
  if (0==strcmp(varname,"")) {*varaddr=NULL;return SUCCESS;}
  if NOT(i=tb_find(deftb,varname))  EXPECTED_ERROR("unknown symbol %s",varname);
  if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("symbol %s is not a variable",varname);
  if (tb_type(deftb,i)!=t_real) EXPECTED_ERROR("symbol %s is not double",varname);
  if (writable && ((tb_flag(deftb,i)&f_ro) )==f_ro) EXPECTED_ERROR("symbol '%s' should be writable but is read-only",varname);
  *varaddr = (REAL *)tb_addr(deftb,i);
  return 1;
}

int acceptki (const char *key,const int writable,const char *deflt,INT **varaddr,char *varname,char *w)
{
  int i;
  if NOT(accepts(key,varname,deflt,w)) return FAILURE;
  if (0==strcmp(varname,"")) {*varaddr=NULL;return SUCCESS;}
  if NOT(i=tb_find(deftb,varname))  EXPECTED_ERROR("unknown symbol %s",varname);
  if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("symbol %s is not a variable",varname);
  if (tb_type(deftb,i)!=t_int) EXPECTED_ERROR("symbol %s is not integer",varname);
  if (writable && ((tb_flag(deftb,i)&f_ro) )==f_ro) EXPECTED_ERROR("symbol '%s' should be writable but is read-only",varname);
  *varaddr = (INT *)tb_addr(deftb,i);
  return 1;
}

int stateDimensionsExist (char *parameterString)
{
#define CHECK(par) if (find_key(#par"=",parameterString)) return 1
  CHECK(xmax);
  CHECK(ymax);
  CHECK(zmax);
#undef CHECK
  return 0;
}

int spaceParametersExist (char *parameterString)
{
/* #define CHECK(par) if (find_key(#par"=",parameterString)) return 1 */
#define CHECK(par) if (find_key(#par"=",parameterString)) {MESSAGE("!!! found "#par" !!!\n"); return 1;}
  CHECK(x0);
  CHECK(x1);
  CHECK(y0);
  CHECK(y1);
  CHECK(z0);
  CHECK(z1);
  CHECK(listpoints);
  CHECK(relist);
  CHECK(where);
#undef CHECK
  return 0;
}

int spaceIsBox (char *parameterString)
{
#define CHECK(par) if (find_key(#par"=",parameterString)) return 0
  CHECK(listpoints);
  CHECK(where);
  CHECK(relist);
#undef CHECK
  return 1;
}


/*
 * Accepts device's space parameters.
 */
int accept_space (Space *S, char *w)
{
  int iv;
  char name[maxname];
  int hasq;
  /* DEBUG("accept_space enter\n"); */
  ACCEPTI(nowhere,0,0,1);

  /* Exclude all space associations for a "nowhere" device */
  if (nowhere) {
    if (spaceParametersExist(w))
      EXPECTED_ERROR(" Space parameters are provided when 'nowhere' is set to 1. "
		     "Either set 'nowhere' to 0 (default), or remove the space parameters.\n");
    
    /* The box is restricted to first corner boundary/halo point, which is ordinarily unused. */
#if MPI 
    S->x0 = S->x1 = local_xmin - 1;
    S->y0 = S->y1 = local_ymin - 1;
    S->z0 = S->z1 = local_zmin - 1;
#else /* not MPI */
    S->x0 = S->x1 = 0;
    S->y0 = S->y1 = 0;
    S->z0 = S->z1 = 0;
#endif  /* not MPI */

    /* assign default "nowhere" values without further ado */
    S->restricted=1; /* the value is in fact immaterial? */
    S->relist=&never;
    S->listpoints=0; /* DEBUG("S->listpoints=%d\n",S->listpoints); */
    S->where[0]='\0';
    S->wherecode=NULL;
    S->hasu=0;
    S->hasg=0;
    S->u=NULL;
    S->geom=NULL;
    S->np=0;
    S->p=NULL;

    /* "nowhere" means no grid points, but every thread */
    S->runHere = 1;

    return SUCCESS;
  } /* if nowhere */
  /* DEBUG("somewhere\n"); */
  
  /* "nowhere" is false, so accept the space parameters. */
  ACCEPTI(restricted,1,0,1);
  ACCEPTSN(where,MAXSTRLEN,"");
  if (*where) {
    DEBUG("*where\n");
    memcpy(loctb,deftb,sizeof(*deftb));
    /* Check if 'where' expression depends on space */
    S->hasu=0;
    for (iv=0;iv<vmax;iv++) {
      snprintf(name,maxname,"u%d",iv);
      if (k_expr_depends (where,name)) S->hasu=1;
    }
    S->hasg=0;
    for (iv=0;iv<geom_vmax;iv++) {
      snprintf(name,maxname,"geom%d",iv);
      if (k_expr_depends (where,name)) S->hasg=1;
    }
    hasq=0;
    if (k_expr_depends (where,"x")) hasq=1;
    if (k_expr_depends (where,"y")) hasq=1;
    if (k_expr_depends (where,"z")) hasq=1;
    if (!hasq && !(S->hasu) && !(S->hasg))
      MESSAGE("/* Warning: the 'where' expression has no space dependencies! */");
    /* Allocate the space-dependence variables in the local symbol table */
    if (S->hasu) {
      CALLOC(S->u,vmax,sizeof(real));
      for(iv=0;iv<vmax;iv++) {
	snprintf(name,maxname,"u%d",iv);
	tb_insert_real(loctb,name,&(S->u[iv]));  CHK(name);
      }
    }
    if (S->hasg) {
      CALLOC(S->geom,geom_vmax,sizeof(real));
      if (geom_vmax) {
	for(iv=0;iv<geom_vmax;iv++) {
	  snprintf(name,maxname,"geom%d",iv);
	  tb_insert_real(loctb,name,&(S->geom[iv]));  CHK(name);
	}
      }
    }
    tb_insert_int (loctb,"x",&(S->x));	    CHK("x");
    tb_insert_int (loctb,"y",&(S->y));	    CHK("y");
    tb_insert_int (loctb,"z",&(S->z));	    CHK("z");
    /* Now the expression can be compiled */
    S->wherecode=compile(where, loctb, t_real); CHK(where);
    /* For the moment, 'where' feature demands listing the points */
    ACCEPTI(listpoints,1,1,1);  DEBUG("S->listpoints=%d\n",S->listpoints);
    /* And we need to know when if ever the point list is to be updated */
    ACCEPTKR(relist,0,"never");
  } else { /* not *where */
    /* DEBUG("not *where\n"); */
    S->wherecode=NULL;
    S->hasu=0;
    S->hasg=0;
    S->u=NULL;
    S->geom=NULL;
    /* Without 'where', we may or may not list the points */
    ACCEPTI(listpoints,(GEOMETRY_ON?1:0),0,1);  /* DEBUG("S->listpoints=%d\n",S->listpoints); */
    /* But there is no point in repeating it ever */
    /* DEBUG("listpoints=%d\n",listpoints); */
    ACCEPTKR(relist,0,"never");
    if (0!=strcmp(S->relistname,"never"))
      MESSAGE("/* Warning: the points list without 'where=' expression does not need "
	      "to be updated, so relist= parameter is superfluous! */");
  } /* not *where */

  /* Enclosing box */
  ACCEPTL(x0, SPACE_DEFAULT_X0, 0,     (int)xmax-1);
  ACCEPTL(x1, SPACE_DEFAULT_X1, S->x0, (int)xmax-1);
  ACCEPTL(y0, SPACE_DEFAULT_Y0, 0,     (int)ymax-1);
  ACCEPTL(y1, SPACE_DEFAULT_Y1, S->y0, (int)ymax-1);
  ACCEPTL(z0, SPACE_DEFAULT_Z0, 0,     (int)zmax-1);
  ACCEPTL(z1, SPACE_DEFAULT_Z1, S->z0, (int)zmax-1);
  /* For seq mode, local==global. */
  /* For par mode, locals will be constrained by the subdomain. */
  S->global_x0 = S->x0;
  S->global_x1 = S->x1;
  S->global_y0 = S->y0;
  S->global_y1 = S->y1;
  S->global_z0 = S->z0;
  S->global_z1 = S->z1;
#if MPI
  /* For parallel mode, local may be narrower than global */
                
  /*      Does the device box intersect with my subdomain?	*/
  if ((S->global_x0 < local_xmin   && S->global_x1 < local_xmin   && mpi_ix != 0 )        || /* < X */
      (S->global_x0 > local_xmax-1 && S->global_x1 > local_xmax-1 && mpi_ix != mpi_nx-1 ) || /* > X */
      (S->global_y0 < local_ymin   && S->global_y1 < local_ymin   && mpi_iy != 0 )        || /* < Y */
      (S->global_y0 > local_ymax-1 && S->global_y1 > local_ymax-1 && mpi_iy != mpi_ny-1 ) || /* > Y */
      (S->global_z0 < local_zmin   && S->global_z1 < local_zmin   && mpi_iz != 0 )        || /* < Z */
      (S->global_z0 > local_zmax-1 && S->global_z1 > local_zmax-1 && mpi_iz != mpi_nz-1 ) ){ /* > Z */
                        
    /*      Global space doesn't intersect with local subdomain. */
    /*      Disable device and set space to local minima         */
    /*      to prevent CREATE function problems.                 */
    S->runHere = 0;
    S->x0 = local_xmin - 1;
    S->x1 = local_xmin - 1;
    S->y0 = local_ymin - 1;
    S->y1 = local_ymin - 1;
    S->z0 = local_zmin - 1;
    S->z1 = local_zmin - 1;
    S->listpoints = 0;  DEBUG("S->listpoints=%d\n",S->listpoints);
    S->where[0] = '\0';
    S->wherecode = NULL;
    S->np = 0;
    S->p = NULL;
      
  } else { /* not "global outside local" */
      
    /* 	    Global box does intersect with the subdomain.   */
    /* 	    Now we constrain the local box by the subdomain */
    /*      Allow spaces to include absolute boundaries.    */      
    S->runHere = 1;

#define CONSTRAINT(q) \
    if (S->q##0 < local_##q##min) {			    \
      if (S->global_##q##0 == 0 && local_##q##min == 1) {   \
	S->q##0 = 0;					    \
      } else {						    \
	S->q##0 = local_##q##min;			    \
      }							    \
    }							    \
    if (S->q##1 > local_##q##max-1) {					\
      if (S->global_##q##1 == q##max-1 && local_##q##max == q##max-1) { \
	S->q##1 = q##max-1;						\
      } else {								\
	S->q##1 = local_##q##max-1;					\
      }									\
    }
    CONSTRAINT(x);
    CONSTRAINT(y);
    CONSTRAINT(z);
#undef CONSTRAINT
    
  } /* not "global outside local" */
#endif /* MPI */
    
  /* DEBUG("S->listpoints=%d\n",S->listpoints); */
  if (S->listpoints) {
    S->p=NULL; /* signal that it has not be alloc'd yet */
    /* DEBUG("about to list_space\n"); */
    if (SUCCESS != list_space(S))
      EXPECTED_ERROR("could not make list of points");
  }
  
  return SUCCESS;
} /* accept_space */

/* Make the list of points is a separate routine as it may be called on again */
int list_space (Space *S)
{
  DEVICE_CONST(int,listpoints);
  /* DEBUG("list_space: listpoints=%d\n",listpoints); */
  if (!listpoints) return SUCCESS;
  DEVICE_CONST(int,x0);
  DEVICE_CONST(int,x1);
  DEVICE_CONST(int,y0);
  DEVICE_CONST(int,y1);
  DEVICE_CONST(int,z0);
  DEVICE_CONST(int,z1);
  DEVICE_CONST(int,global_x0);
  DEVICE_CONST(int,global_x1);
  DEVICE_CONST(int,global_y0);
  DEVICE_CONST(int,global_y1);
  DEVICE_CONST(int,global_z0);
  DEVICE_CONST(int,global_z1);
  DEVICE_CONST(int,restricted);
  DEVICE_ARRAY(char,where);
  DEVICE_CONST(pp_fn,wherecode);
  DEVICE_CONST(int,hasu);
  DEVICE_CONST(int,hasg);
  DEVICE_VAR(INT,x);
  DEVICE_VAR(INT,y);
  DEVICE_VAR(INT,z);
  DEVICE_VAR(size_t,np);
  DEVICE_VAR(devicePoint *,p);
  S->relisted=1;

  /* If device's space is default then its point list is standard */
  /* DEBUG("numTissuePoints=%ld, TissuePoints=%p\n",numTissuePoints,TissuePoints); */
  if (TissuePoints) {
    DEBUG("this (%d:%d)x(%d:%d)x(%d:%d) vs dflt (%d:%d)x(%d:%d)x(%d:%d)\n",
    	  global_x0,global_x1,global_y0,global_y1,global_z0,global_z1,
    	  SPACE_DEFAULT_X0,SPACE_DEFAULT_X1,
    	  SPACE_DEFAULT_Y0,SPACE_DEFAULT_Y1,
    	  SPACE_DEFAULT_Z0,SPACE_DEFAULT_Z1);
    /* DEBUG("wherecode=%p\n",wherecode); */
    if (restricted==1 &&
	global_x0==SPACE_DEFAULT_X0 && global_x1==SPACE_DEFAULT_X1 &&
	global_y0==SPACE_DEFAULT_Y0 && global_y1==SPACE_DEFAULT_Y1 &&
	global_z0==SPACE_DEFAULT_Z0 && global_z1==SPACE_DEFAULT_Z1 &&
	wherecode==NULL) {
      *np=numTissuePoints;
      CALLOC((*p),numTissuePoints,sizeof(devicePoint));
      memcpy((*p),TissuePoints,numTissuePoints*sizeof(devicePoint));
      /* DEBUG("TissuePoints(%p) -> *p(%p) [%ld]\n",TissuePoints,*p,numTissuePoints); */
      /* DEBUG("TissuePoints=%p [0]->(%ld,%ld,%ld)\n", */
      /* 	    TissuePoints,TissuePoints[0].x,TissuePoints[0].y,TissuePoints[0].z); */
      /* DEBUG("p=%p [0]->(%ld,%ld,%ld)\n",(*p),(*p)[0].x,(*p)[0].y,(*p)[0].z); */
      
      /* *p=TissuePoints; */
      Debug("/* default list of %ld points */",(long)*np);
      return SUCCESS;
    }
    DEBUG("*np=%ld\n",(*np));
  }

  FREE(*p);
  CALLOC(*p,(x1-x0+1)*(y1-y0+1)*(z1-z0+1),sizeof(devicePoint));
  *np=0;
  /* DEBUG("where='%s' wherecode=%p\n",where,wherecode); */
  for ((*x)=x0;(*x)<=x1;(*x)++) {
    for ((*y)=y0;(*y)<=y1;(*y)++) {
      for ((*z)=z0;(*z)<=z1;(*z)++) {
	if ((restricted==1) && (isTissue(*x,*y,*z)==0)) continue;
	if (wherecode != NULL) {
	  if (hasu) memcpy(S->u,New+ind(*x,*y,*z,0),vmax*sizeof(real));
	  if (hasg) memcpy(S->geom,Geom+geom_ind(*x,*y,*z,0),geom_vmax*sizeof(real));
	  if (0.0==*(real *)execute(wherecode)) continue;
	}
	(*p)[*np].x=(*x);
	(*p)[*np].y=(*y);
	(*p)[*np].z=(*z);
	(*p)[*np].u=New+ind((*x),(*y),(*z),0);
	(*p)[*np].X=NULL;
	(*np)++;
      } /* for *z */
    } /* for *y */
  } /* for *x */
  /* DEBUG("*np=%ld = s->np=%ld\n",(*np),S->np); */
  REALLOC(*p,(*np)*sizeof(devicePoint));
  return SUCCESS;
} /* list_space */

int accept_window (BGIWindow *S, char *w)
{
  ACCEPTI(row0, 0, 0,INONE);
  ACCEPTI(row1, 0, 0,INONE);
  ACCEPTI(col0, 0, 0,INONE);
  ACCEPTI(col1, 0, 0,INONE);
  ACCEPTI(color,15,0,  255);
  return 1;
}

/* Typical constants for accept functions */

static real NEXTTO1=0.0; /* these are initialized */
static real PREDTO1=0.0; /*    in init_const() */

real RSUCC (real val)
{
  if      (val>0.0) return (val*NEXTTO1);
  else if (val<0.0) return (val*PREDTO1);
  else              return (FLT_MIN);
}

real RPRED (real val)
{
  if      (val>0.0) return (val*PREDTO1);
  else if (val<0.0) return (val*NEXTTO1);
  else		    return (-FLT_MIN);
}

static double _if   (double a, double b, double c) {if (a   ) return(b); else return(c);}
static double _ifne0(double a, double b, double c) {if (a!=0) return(b); else return(c);}
static double _ifeq0(double a, double b, double c) {if (a==0) return(b); else return(c);}
static double _ifgt0(double a, double b, double c) {if (a >0) return(b); else return(c);}
static double _ifge0(double a, double b, double c) {if (a>=0) return(b); else return(c);}
static double _iflt0(double a, double b, double c) {if (a <0) return(b); else return(c);}
static double _ifle0(double a, double b, double c) {if (a<=0) return(b); else return(c);}
static double _ifsign(double a, double b, double c, double d) {if (a<0) return b;else if (a>0) return d;else return c;}
static double _gt(double a, double b) {return (double)(a> b);}
static double _ge(double a, double b) {return (double)(a>=b);}
static double _lt(double a, double b) {return (double)(a< b);}
static double _le(double a, double b) {return (double)(a<=b);}
static double _eq(double a, double b) {return (double)(a==b);}
static double _ne(double a, double b) {return (double)(a!=b);}
static double _mod(double a, double b) {return fmod(a,b);}
static double _max(double a, double b) {return (a>b)?a:b;}
static double _min(double a, double b) {return (a<b)?a:b;}
static double _crop(double a, double b, double c) {return (a<b)?b:(a>c)?c:a;}
static double _crop01(double a, double b, double c) {return (a<b)?0:(a>c)?1:(a-b)/(c-b);}
static double j2(double x) {return jn(2,x);}
static double j3(double x) {return jn(3,x);}
/* A long way to define a function giving value at a certain point */
static double _U(int x,int y,int z,int v) {
  #if MPI
    double result;
    int the_rank;
    the_rank = getRankContainingPoint(x,y,z);
    if (the_rank<0) return RNONE;
    if (mpi_rank==the_rank) {
      result = (double) New[ind(x,y,z,v)];
    }
    MPIDO(MPI_Bcast(&result, 1, MPI_DOUBLE, the_rank, ALL_ACTIVE_PROCS),"Could not broadcast the result.");
    /* if (mpi_rank==0) fprintf(stdout,"_U(%d,%d,%d,%d)=%lg\n",x,y,z,v,result); */
    return result;
  #else
    return(double)New[ind(x,y,z,v)];
  #endif
}
/* Similar function, to access the geometry data */
static double _Geom(int x,int y,int z,int v) {
  if (!GEOMETRY_ON) return RNONE;
  #if MPI
    double result;
    int the_rank;
    the_rank = getRankContainingPoint(x,y,z);
    if (the_rank<0) return RNONE;
    if (mpi_rank==the_rank) {
      result = (double) Geom[geom_ind(x,y,z,v)];
    }
    MPIDO(MPI_Bcast(&result, 1, MPI_DOUBLE, the_rank, ALL_ACTIVE_PROCS),"Could not broadcast the result.");
    return result;
  #else
    return (double)Geom[geom_ind(x,y,z,v)];
  #endif
}
static double intp(
  double (f)(int,int,int,int), 
  double _x, double _y, double _z, double _v,
  int _vmax
) {
  int v=floor(_crop(_v,0,(double)_vmax-1));
  int x=floor(_crop(_x,SPACE_DEFAULT_X0,SPACE_DEFAULT_X1));double dx=_x-x;double px=1-dx;
  int y=floor(_crop(_y,SPACE_DEFAULT_Y0,SPACE_DEFAULT_Y1));double dy=_y-y;double py=1-dy;
  int z=floor(_crop(_z,SPACE_DEFAULT_Z0,SPACE_DEFAULT_Z1));double dz=_z-z;double pz=1-dz;
  double        result=         px*py*pz*f(x  ,y  ,z  ,v); if(dim==0) return result;
  if(dx)        result+=        dx*py*pz*f(x+1,y  ,z  ,v); if(dim==1) return result;
  if(dy)        result+=        px*dy*pz*f(x  ,y+1,z  ,v);
  if(dx&&dy)    result+=        dx*dy*pz*f(x+1,y+1,z  ,v); if(dim==2 ||dz==0) return result;  
                result+=        px*py*dz*f(x  ,y  ,z+1,v);
  if(dx)        result+=        dx*py*dz*f(x+1,y  ,z+1,v);
  if(dy)        result+=        px*dy*dz*f(x  ,y+1,z+1,v);
  if(dx&&dy)    result+=        dx*dy*dz*f(x+1,y+1,z+1,v); return result;
}
double _u(double _x, double _y, double _z, double _v) {
  return intp(_U,_x,_y,_z,_v,vmax);
}
static double _geom(double _x, double _y, double _z, double _v) {
  if (!GEOMETRY_ON) return RNONE;
  return intp(_Geom,_x,_y,_z,_v,geom_vmax);
}
#if MPI
#else
  /* Uniformly distributed in the [a,b] interval. Is "MT" generator better than this one? */
  static double _rnd(double a,double b) {return (a+(b-a)*drand48());}
  static double _rndseed(double a) {unsigned long seedval=floor(a); srand48(seedval); return 0;}
  /* Gaussian with mean a and stddev b, using Box-Muller 1958 transformation */
  static double _gauss(double a,double b) {return a+b*(sqrt(-2*log(_rnd(0,1)))*cos(2*pi*_rnd(0,1)));}
#endif

/* initiation of constants used in checking limits etc */
int init_const (void)
{
  real eps=1.0, eps2;           /* machine epsilon */
  volatile real nextto1;
  #define STRLEN MAXPATH
  static char s1[STRLEN];
  char *p;
  int iarg;

  do {
    eps2=eps;
    eps = eps*0.5;
    nextto1 = 1.0 + eps;
  } while (nextto1 > 1.0);
  macheps = eps2;
  NEXTTO1 = 1.0 + eps2;
  PREDTO1 = 1.0 - eps2;
  INONE = (-INT_MAX+1);
  LNONE = (-LONG_MAX+1L);
  RNONE = (-MAXREAL)*PREDTO1;
  real_inf = MAXREAL;
  /* predefined constants for def's */
  k_init();
  deftb = tb_new();                                             CHK("deftb");
  loctb = tb_new();                                             CHK("loctb");

  tb_insert_fun(deftb,"atan2",atan2,2);                         CHK("atan2");
  tb_insert_fun(deftb,"erf",erf,1);                             CHK("erf");
  tb_insert_fun(deftb,"hypot",hypot,2);                         CHK("hypot");
  tb_insert_fun(deftb,"J0",j0,1);                               CHK("j0");
  tb_insert_fun(deftb,"J1",j1,1);                               CHK("j1");
  tb_insert_fun(deftb,"J2",j2,1);                               CHK("j2");
  tb_insert_fun(deftb,"J3",j3,1);                               CHK("j3");
  tb_insert_fun(deftb,"tanh",tanh,1);                           CHK("tanh");
#define FU(name, args) tb_insert_fun(deftb,#name,_##name,args);       CHK(#name);
  FU(ifsign,4); FU(if,3); 
  FU(ifeq0,3); FU(ifne0,3); FU(ifgt0,3); FU(ifge0,3); FU(iflt0,3); FU(ifle0,3);
  FU(eq,2); FU(ne,2); FU(gt,2); FU(ge,2); FU(lt,2); FU(le,2);
  FU(mod,2); FU(max,2); FU(min,2); FU(crop,3); FU(crop01,3); FU(u,4); FU(geom,4);

  #if MPI
  #else
  FU(rnd,2); FU(rndseed,1); FU(gauss,2); /* these functions are not "safe" in parallel runs */
  #endif

  #undef FU

  /* Are there any read/write variables? */
  /* #define VB(name) tb_insert_int(deftb,#name,&name); CHK(#name); */
  /* #undef VB */

  /* Read-only integer variables */
  #define RO(name) tb_insert_abstract(deftb,#name,t_int,(p_vd)&name,0,f_ro); CHK(#name);
  RO(t);
  RO(inf);
  RO(xmax); RO(ymax); RO(zmax);
  RO(Graph); RO(graphon); RO(online); RO(XMAX); RO(YMAX); RO(WINX); RO(WINY); 
  #undef RO

  /* C-variable names of colours are different from their k-names */
  #define CLR(name) tb_insert_abstract(deftb,#name,t_int,(p_vd)&var##name,0,f_ro); CHK(#name);
  CLR(BLACK);
  CLR(BLUE);
  CLR(GREEN);
  CLR(CYAN);
  CLR(RED);
  CLR(MAGENTA);
  CLR(BROWN);
  CLR(LIGHTGRAY);
  CLR(DARKGRAY);
  CLR(LIGHTBLUE);
  CLR(LIGHTGREEN);
  CLR(LIGHTCYAN);
  CLR(LIGHTRED);
  CLR(LIGHTMAGENTA);
  CLR(YELLOW);
  CLR(WHITE);
  #undef CLR

  /* Read-only real variables */
  #define RO(name) tb_insert_abstract(deftb,#name,t_real,(p_vd)&name,0,f_ro); CHK(#name);
  RO(pi);
  RO(always); 
  RO(never); 
  RO(real_inf);
  #undef RO

  snprintf(buf,MAXSTRLEN,"int narg=%d",narg);
  if (!def(buf)) EXPECTED_ERROR("cannot declare \"%s\"",buf);

  for(iarg=0;iarg<=narg;iarg++) {
    if (iarg==0) {
      strncpy(s1,inname[0],STRLEN);
      /* if (NULL!=(p=strtok(s1,"."))) strcpy(s1,p); */
      p=strtok(s1,".");
      if (p!=NULL)
	memmove(s1,p,strlen(p)+1);
      snprintf(buf,MAXSTRLEN,"str %d %s",iarg,s1);
    } else {
      snprintf(buf,MAXSTRLEN,"str %d %s",iarg,arg[iarg-1]);
    }
    if (!def(buf)) EXPECTED_ERROR("cannot declare \"%s\"",buf);
  }

  return 1;
}

int term_const(void) { /* clear the symbol table and deallocate variables */
  k_term();
  return 1;
}
