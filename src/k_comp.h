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

/* Compile a k-program contained in a pgm={} block */

#include <regex.h>
#include "error.h"

#undef SEPARATORS
#define SEPARATORS ";"
#define BLANK " \n\t\r"
#define MAX_LOCAL_VARS 32

#define FIND_UNINITIALISED(RHS)						\
for (v=0; v<MAX_LOCAL_VARS; v++) { 					\
  if (uninitialised_vars[v] != NULL) {					\
    snprintf(pattern, 80, varFormat, uninitialised_vars[v]);		\
    if (regcomp(&compiled, pattern, REG_EXTENDED) != 0) {		\
      EXPECTED_ERROR("Could not compile regular expression required to evaluate k_func code.\n"); \
    }									\
    int nsub = compiled.re_nsub;					\
    regmatch_t matchptr[nsub];						\
    if (0!=(err = regexec(&compiled, RHS, nsub, matchptr, 0))) {	\
      if (err == REG_ESPACE) {						\
	EXPECTED_ERROR("Ran out of memory when evaluating regular expression.\n"); \
      }									\
    } else {								\
      EXPECTED_ERROR("The expression \'%s\' contains an uninitialised local variable, '%s'. \n" \
                     "Ensure that a value is	assigned to '%s' before use.\n"			\
		      ,RHS, uninitialised_vars[v], uninitialised_vars[v]); \
    }									\
    regfree(&compiled);							\
  }									\
}

BEGINBLOCK("pgm=",buf); {
  int idata, icode, v;
  char *pcode, *pdata, *initCode;

  /*  RegEx Necessities: needed to look for uninitialized variables */
  char *varFormat = "((^|[^A-Za-z0-9_])(%s)([^A-Za-z0-9_]|$))";
  char pattern[80];
  regex_t compiled;
  int err;

  /*  List of uninitialised variables */
  char *typename, *varname, *expr;
  char *uninitialised_vars[MAX_LOCAL_VARS];
  int num_uninitialised = 0;
  for (v=0;v<MAX_LOCAL_VARS;v++) {uninitialised_vars[v]=NULL;}
	
  /* Count the number of pieces in the content that are separated by ';' */
  char *s1=strdup(buf);
  for (S->ncode=0;strtok(S->ncode?NULL:s1,SEPARATORS);S->ncode++);
  if (!S->ncode) EXPECTED_ERROR("no statements in \"%s\"",buf);
  FREE(s1);


  CALLOC(S->code,S->ncode,sizeof(pp_fn))
  CALLOC(S->data,S->ncode,sizeof(p_real))
  char *lines[S->ncode];
  
  /* Split the block content to "lines" over ';' separators */
  for (icode=0;icode<S->ncode;icode++) {
    /* each line to start from a non-blank character */
    for (pdata=strtok(icode?NULL:buf,SEPARATORS);*pdata&&strchr(BLANK,*pdata);pdata++);
    lines[icode] = pdata;
  }

  for (icode=0;icode<S->ncode;icode++) {
    pdata = lines[icode];
    
    /* If there is nothing in this line, skip it. */
    if (!*pdata) {
      S->code[icode]=NULL;
      continue;
    }
    
    /*  Find LHS symbol: from start of the line to the first blank, if any */
    for (pcode=pdata;(*pcode)&&(*pcode!='=');pcode++) {
      if (strchr(BLANK,*pcode)) {
	*pcode='\0';
	pcode++;
	break;
      }
    }
    
    /*  If the LHS is 'def', define a local variable. */
    if (strcmp(pdata,"def")==0) {
      
      /*  Duplicate pcode because def_local() has side effects. */
      char *codeToCheck = strdup(pcode);
      
      /* MESSAGE("def %s\n", codeToCheck); */
      
      /*  Check statement to see if it's initialised. */
      typename = strtok(codeToCheck," ;\t\r\n!$");
      varname  = strtok(NULL," ;\t\r\n!$=");
      expr     = strtok(NULL,"");
      
      /*  Check value expression for uninitialised variables. */
      if (expr !=NULL) {
	FIND_UNINITIALISED(expr);
      }
      
      if (!def_local(pcode, loctb)) {
	EXPECTED_ERROR("Couldn't define local variable with \"%s\".\nDefinition syntax is 'def <type> <name> <value>'.\n", pcode);
      }
      
      /*  If variable is being initialised: */
      if (expr!=NULL) {
	/*  Compile a new statement with the new local var on the LHS and its initial value on the RHS. */
	idata=tb_find(loctb,varname);
	S->code[icode] = compile(expr,loctb,tb_type(loctb,idata)); CHK(expr);
	S->data[icode] = (p_vb) tb_addr(loctb,idata); /* !!!! this was missing !!!! */
      } else {
	/*  The variable is uninitialised, so add it to the list of variables to check for. */
	if(num_uninitialised>=MAX_LOCAL_VARS)
	  EXPECTED_ERROR("Exceeded maximum number of uninitialised local variables.\n");
	
	for (v=0; v<MAX_LOCAL_VARS; v++) {
	  if (!uninitialised_vars[v]) {
	    uninitialised_vars[v]=strdup(varname);
	    num_uninitialised++;
	    break;
	  }
	}
      }
     
      free(codeToCheck);

      /*  Don't compile this instruction. */
      pcode++;
      
    } else { /*  LHS is not 'def' */
      
      /*  Find the assignment operator ... */
      for (;(*pcode)&&(*pcode!='=');pcode++)
	if(strchr(BLANK,*pcode)) *pcode='\0';
      if (!*pcode)
	EXPECTED_ERROR(": no \'=\' sign in expression #%d: \'%s\'",icode+1,pdata); 
      /*  ... and remove it */
      *(pcode++)='\0';
      
      /*  By this point, the = sign will be wiped out and pcode will be pointing to  */
      /*  the first char of the RHS. */
      
      if NOT(idata=tb_find(loctb,pdata)) EXPECTED_ERROR("unknown symbol '%s'",pdata);
      if ((tb_flag(loctb,idata)&f_ro)==f_ro) EXPECTED_ERROR("symbol '%s' is read-only",pdata);
      if ((tb_flag(loctb,idata)&f_vb)==0) EXPECTED_ERROR("symbol '%s' is not a variable",pdata);
      
      FIND_UNINITIALISED(pcode)
	
      S->data[icode] = (p_vb) tb_addr(loctb,idata);
			
      /*
	Check if LHS exists in deftb (i.e is it global?).
	If so, check RHS for any local variables. 
	If present, report an error.
      */
      
      /*  Is the LHS a global variable? */
      if (tb_find(deftb,pdata)) {
	if (!dev->s.nowhere) {
	  EXPECTED_ERROR("`nowhere` must be true (1) in order to make assignments to global variables.\n");
	} else {
	  /* compile the expression using local symbol table */
	  S->code[icode] = compile(pcode,loctb,tb_type(loctb,idata)); CHK(pcode);
	  if (Verbose) MESSAGE("\x01""\n\t%s=%s;",pdata,pcode);
	  
	  /*  If the expression can also compile using deftb, the RHS does not contain local variables. */
	  pp_fn tmp;
	  tmp = compile(pcode,deftb,tb_type(loctb,idata));
	  free(tmp);
	  if (err_code) { /* it does not compile with global vars only, so contains local vars */
	    EXPECTED_ERROR(
	      "The line \"%s=%s%c\" attempts to assign a local value to a global variable (%s), which isn't allowed.\n",
	      pdata,pcode,SEPARATORS[0],pdata
	    );
	  }
	}
      } else { /*  LHS is a local variable. */				
	if (dev->s.nowhere) {
	  EXPECTED_ERROR("`nowhere` must be false (0) to make assignments to local variables.\n");
	} else {
	  for (v=0; v<MAX_LOCAL_VARS; v++) {
	    if (uninitialised_vars[v]!=NULL && strcmp(uninitialised_vars[v], pdata)==0) {
	      if (num_uninitialised<=0) {
		EXPECTED_ERROR("Cannot remove initialised variable \'%s\'---list empty.\n",pdata);
	      }
	      free(uninitialised_vars[v]);
	      uninitialised_vars[v] = NULL;
	      num_uninitialised--;
	      break;
	    }
	  }	
	  S->code[icode] = compile(pcode,loctb,tb_type(loctb,idata)); CHK(pcode);
	  if (Verbose) MESSAGE("\x01""\n\t%s=%s;",pdata,pcode);
	}
      } /*  else LHS is local  */
    } /* else (LHS != 'def') */
  } /*  for icode */

  if (num_uninitialised) {
    for (v=0; v<MAX_LOCAL_VARS; v++) {
      if (uninitialised_vars[v] != NULL) {
	MESSAGE("/* WARNING: Local variable \'%s\' is unused. */\n",uninitialised_vars[v]);
	free(uninitialised_vars[v]);
      }
    }
  }
} ENDBLOCK;
#undef FIND_UNINITIALISED
#undef SEPARATORS
#undef BLANK
