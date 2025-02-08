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

#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include "pipe.h"

PIPE *pipeto(char *cmd) {
  PIPE *p=calloc(1,sizeof(PIPE));
  int pid;
   p->n=tempnam(NULL,NULL); // original code

  if (-1==mkfifo(p->n,0600)) {perror("pipeto could not make a fifo"); return NULL;}
  switch (pid=fork()) {
  case -1: 
    perror("pipeto could not fork"); 
    return NULL;
  case 0: {
    int len=strlen(cmd)+strlen(p->n)+10;
    char *s=calloc(len,1);
    snprintf(s,len,"cat %s | %s",p->n,cmd); 
    system(s);
    /* printf("%s:%d %ld ben %ld command '%s' completed\n",__FILE__,__LINE__,(long int)getpid(),(long int)getppid(),s); */
    free(s);
    /* _exit(0);					      */
    exit(0);					     
  }
  default: 
    if ((p->f=fopen(p->n,"w"))==NULL) {perror("pipe could not write to fifo"); return NULL;}
    p->child=pid;
    return p; 
 }
}

int pipeclose(PIPE *p) {
  int fcloseret=0;
  int killret=0;
  int status=0;
  if (0!=(fcloseret=fclose(p->f))) perror("pipeclose fclose fifo input");
  /* waitpid(p->child,&status,WNOHANG); /\* if child already finished, there is no one there to report the status *\/ */
  waitpid(p->child,&status,0); /* VNB 2024/09/10: with WNOHANG, child procs remained <defunct>, fork resource run out  */
  if (0==(WIFEXITED(status))) fprintf(stderr,"pipeclose child status %08x\n",status);
  unlink(p->n);			      
  free(p->n);
  /* printf("%s:%d %ld finished %ld\n",__FILE__,__LINE__,(long int)getpid(),(long int)p->child); */
  return fcloseret | killret;
}
