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

/* Growable text buffers, and texts of grid points collected in global order. */
/* See textbuf.h for the idea. */

#if MPI
#include <mpi.h>
#endif

#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "textbuf.h"

#define PIECE_LL 3		/* number of long long's in a piece, for MPI transfers */

/****************/
/* Text buffers */

void txt_write(textbuf *b, const char *s, size_t n) {
  if (b->n+n+1 > b->max) {
    b->max = 2*(b->n+n+1);
    REALLOC(b->s,b->max);
  }
  memcpy(b->s+b->n,s,n);
  b->n += n;
  b->s[b->n]='\0';
}

void txt_puts(textbuf *b, const char *s) {
  txt_write(b,s,strlen(s));
}

void txt_printf(textbuf *b, const char *fmt, ...) {
  va_list ap;
  int n;
  va_start(ap,fmt);
  n=vsnprintf(b->s+b->n,b->max-b->n,fmt,ap);
  va_end(ap);
  if (n<0) ABORT("vsnprintf failed on format '%s'\n",fmt);
  if (b->n+n+1 > b->max) {	/* did not fit: enlarge and try again */
    b->max = 2*(b->n+n+1);
    REALLOC(b->s,b->max);
    va_start(ap,fmt);
    vsnprintf(b->s+b->n,b->max-b->n,fmt,ap);
    va_end(ap);
  }
  b->n += n;
}

void txt_free(textbuf *b) {
  FREE(b->s);
  b->n=b->max=0;
}

/*******************/
/* Pieces of text  */

void tp_reset(textpieces *p) {
  p->txt.n=0;
  if (p->txt.s) p->txt.s[0]='\0';
  p->n=0;
}

void tp_begin(textpieces *p, long long key) {
  if (p->n>=p->max) {
    p->max=2*p->max+16;
    REALLOC(p->pc,p->max*sizeof(piece));
  }
  p->pc[p->n].key=key;
  p->pc[p->n].off=p->txt.n;
  p->pc[p->n].len=0;
}

void tp_end(textpieces *p) {
  p->pc[p->n].len=p->txt.n-p->pc[p->n].off;
  p->n++;
}

static int piece_cmp(const void *a, const void *b) {
  long long ka=((const piece *)a)->key;
  long long kb=((const piece *)b)->key;
  return (ka>kb)-(ka<kb);
}

void tp_sort(textpieces *p) {
  if (p->n) qsort(p->pc,p->n,sizeof(piece),piece_cmp);
}

void tp_free(textpieces *p) {
  txt_free(&(p->txt));
  FREE(p->pc);
  p->n=p->max=0;
}

long long tp_key(long x, long y, long z) {
  return ((long long)z*ymax+y)*xmax+x;
}

long tp_key_z(long long key) {
  return (long)(key/((long long)xmax*ymax));
}

#if MPI
int tp_gather(textpieces *mine, textpieces *all, int root, MPI_Comm comm) {
  int nproc, me, i;
  long long sizes[2], *allsizes=NULL;
  int *pccount=NULL, *pcdispl=NULL, *txtcount=NULL, *txtdispl=NULL;
  long long totpc=0, tottxt=0;

  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&me);
  if (mine->n*PIECE_LL>INT_MAX || mine->txt.n>INT_MAX)
    ABORT("too much text output from process %d\n",mpi_rank);
  sizes[0]=mine->n;
  sizes[1]=mine->txt.n;
  if (me==root) CALLOC(allsizes,2*nproc,sizeof(long long));
  MPIDO(MPI_Gather(sizes,2,MPI_LONG_LONG,allsizes,2,MPI_LONG_LONG,root,comm),
	"could not gather sizes of text output");

  if (me==root) {
    CALLOC(pccount,nproc,sizeof(int));
    CALLOC(pcdispl,nproc,sizeof(int));
    CALLOC(txtcount,nproc,sizeof(int));
    CALLOC(txtdispl,nproc,sizeof(int));
    for (i=0;i<nproc;i++) {
      pccount[i]=allsizes[2*i]*PIECE_LL;
      pcdispl[i]=totpc*PIECE_LL;
      txtcount[i]=allsizes[2*i+1];
      txtdispl[i]=tottxt;
      totpc+=allsizes[2*i];
      tottxt+=allsizes[2*i+1];
      if (totpc*PIECE_LL>INT_MAX || tottxt>INT_MAX)
	ABORT("too much text output from all processes\n");
    }
    /* Make room for everything */
    tp_reset(all);
    if (all->max<totpc) {
      all->max=totpc;
      REALLOC(all->pc,all->max*sizeof(piece));
    }
    if (all->txt.max<tottxt+1) {
      all->txt.max=tottxt+1;
      REALLOC(all->txt.s,all->txt.max);
    }
  }
  MPIDO(MPI_Gatherv(mine->pc,(int)(mine->n*PIECE_LL),MPI_LONG_LONG,
		    me==root?all->pc:NULL,pccount,pcdispl,MPI_LONG_LONG,root,comm),
	"could not gather keys of text output");
  MPIDO(MPI_Gatherv(mine->txt.s,(int)mine->txt.n,MPI_CHAR,
		    me==root?all->txt.s:NULL,txtcount,txtdispl,MPI_CHAR,root,comm),
	"could not gather text output");

  if (me==root) {
    long long j;
    all->n=totpc;
    all->txt.n=tottxt;
    all->txt.s[tottxt]='\0';
    /* Offsets of the pieces were local, make them refer to all->txt */
    for (i=0;i<nproc;i++)
      for (j=pcdispl[i]/PIECE_LL;j<(pcdispl[i]+pccount[i])/PIECE_LL;j++)
	all->pc[j].off+=txtdispl[i];
    FREE(allsizes); FREE(pccount); FREE(pcdispl); FREE(txtcount); FREE(txtdispl);
  }
  return 1;
}
#endif
