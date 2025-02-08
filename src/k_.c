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

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include "beatbox.h"
#include "system.h"
#include "error.h"
#include "k_.h"

/* check array of pointers **data of length n for	*/
/* whether it contains address *addr			*/
int used (double **data,int n,double *addr)
{
  int i;
  for (i=0;i<n;i++) if (data[i]==addr) return 1;
  return 0;
}

/* Check of given expression is just one variable of the needed type */
/* from given symtable, and return its address: */
/* 'public' version of tb_find local k_comp.c. */
/* Caveat: assumes all entries in the table are meaningful; */
/* - follow the access method as in tb_find ?? */
p_vb is_variable (char *expr, p_tb tab, k_type type)
{
  int i;
  t_ln entry;
//     strncmp(const char *s1, const char *s2, size_t n);
  for (i=0; i<maxtab; i++) {
    entry=tab->array[i];
    if (entry.np != 0) continue;			/* not a variable */
    if ( (entry.tp & ~_mask) != type) continue;		/* wrong type */
    if ( (entry.tp & _mask) & ~(f_vb|f_ro)) continue;	/* must be r/o or plain var */
    if (strncmp(expr,entry.nm,maxname) != 0) continue;	/* wrong name */
    return entry.ad;
  }
  return NULL;
}

#define MAXSTRLEN 8192

/* check if given (uncompiled) expression depends 	*/
/* on variable varname 					*/
#ifdef MARK_ERROR
#undef MARK_ERROR
#endif
#define MARK_ERROR beatbox_abort("");
int k_expr_depends (char *expr, char *varname)
{
  char pattern[MAXSTRLEN];
  char varFormat[] = "((^|[^A-Za-z0-9_])(%s)([^A-Za-z0-9_]|$))";
  regex_t compiled;
  int nsub;
  regmatch_t *matchptr;
  int err;

  snprintf(pattern, MAXSTRLEN, varFormat, varname);
  if (regcomp(&compiled, pattern, REG_EXTENDED) != 0) {
    EXPECTED_ERROR("Could not compile regular expression required to evaluate k_func code '%s'.\n",pattern);
  }
  nsub = compiled.re_nsub;
  CALLOC(matchptr, nsub, sizeof (regmatch_t));
  err = regexec(&compiled, expr, nsub, matchptr, 0);
  if (err == REG_ESPACE) {
    EXPECTED_ERROR("Ran out of memory when evaluating regular expression '%s'.\n",pattern);
  }
  return (err != 0);
}
