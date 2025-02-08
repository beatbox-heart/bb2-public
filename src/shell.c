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

/* 
 * Execute a command in a subshell
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "error.h"
/* #include "bikt.h" */
/* #include "windraw.h" */
/* #include "k_.h" */
/* #include "mpi_io_choice.h" */

typedef struct {
  int advance;
  char *echo;
  int critical;
  char *cmd;
  char code[1024];
  pp_fn compiled;
} STR;

int doit(char *cmd, pp_fn code, char *echo) {
  char localbuf[2048];
  int rc;
  if (cmd==NULL) return 0;
  if (*cmd=='\0') return 0;
  k_on();
  snprintf(localbuf,2048,cmd,*(real *)execute(code));
  k_off();
  if (echo) {if (*echo) printf(echo,localbuf);}
  if (0!=(rc=system(localbuf)))
    MESSAGE("t=%d shell returned code %d\n",t,rc);
  return rc;
}

RUN_HEAD(shell)
  DEVICE_ARRAY(char, cmd)
  DEVICE_CONST(pp_fn, compiled)
  DEVICE_ARRAY(char, echo)
  DEVICE_CONST(int, critical)
  int rc=doit(cmd,compiled,echo);
  /* We will only be here if mpi_rank==0  --- see below */
  if (rc!=0 && critical) {
    #if MPI
      ABORT("Critical error executing a shell command");
    #else
      return 0;
    #endif
  }
RUN_TAIL(shell)

DESTROY_HEAD(shell)
  FREE(S->cmd);
  FREE(S->echo);
  FREE(S->compiled);
DESTROY_TAIL(shell)

CREATE_HEAD(shell) {
  DEVICE_IS_SPACELESS
  DEVICE_MUST_BE_NOWHERE
  char l[2048];
  int rc;
  CALLOC(S->cmd,MAXSTRLEN,sizeof(char));
  CALLOC(S->echo,MAXSTRLEN,sizeof(char));
  ACCEPTI(advance,0,0,1);
  ACCEPTI(critical,0,0,1);
  ACCEPTS(echo,"beatbox> %s\n");
  if (strlen(S->echo)) REALLOC(S->echo,strlen(S->echo)+1);
  ACCEPTS(cmd,"");
  ACCEPTS(code,"t");
  k_on();
  S->compiled=compile(S->code,deftb,t_real); CHK(S->code);
  k_off();
  if (strlen(S->cmd)) {
    REALLOC(S->cmd,strlen(S->cmd)+1);
  }
#if MPI
  /* To be on the safe side, we only run this device on one process. */
  if (mpi_rank == 0) {
    dev->s.runHere = 1;
    if (advance) {
      rc=doit(S->cmd,S->compiled,S->echo);
      if (rc!=0 && S->critical) return 0;
    }
  } else {
    dev->s.runHere = 0;
  }
#else
  if (advance) {
    rc=doit(S->cmd,S->compiled,S->echo);
    if (rc!=0 && S->critical) return 0;
  }
#endif
} CREATE_TAIL(shell,0)
