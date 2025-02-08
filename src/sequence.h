/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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

typedef struct {
  char mask[MAXPATH];	/* file mask for sequence of files */
  char name[MAXPATH];	/* name of the current file */
  char postproc[10*MAXPATH]; /* post-processing command template */
#if MPI
  MPI_File f;		/* file descriptor of the current file */
  int mpiMode;		/* its MPI mode */
  MPI_Comm comm;	/* parallel i/o communicator */
#else
  FILE *f;		/* file descriptor of the current file */
#endif
  char mode[8];		/* its (desired in MPI) mode */
  int startfrom;	/* starting number for the sequence */
  int autonumber;       /* automatically start first available, if this is nonzero */
  int ignoreempty;      /* ignore empty files when looking for automatic starting number */
  int n;		/* serial number of the next file */
  int mute;		/* open/close events will not be reported if this is nonzero */
  int trial; 		/* the current file is a trial, do not postprocess */
} sequence;

int thisq(sequence *q); /* prepare to write into the current file of the sequence */
int nextq(sequence *q); /* finish with the current file and proceed to the next */

#if MPI
/* MPI version includes runHere so that a new communicator can be created. */
#define ACCEPTQ(b,c,d)   if (!acceptq(#b"=",c,d,dev->s.runHere,&(S->b),w)) return(0); sequence *b=&(S->b)
int acceptq (const char *name,const char *mode,const char *deflt,int runHere,sequence *q, char *w);
#else
#define ACCEPTQ(b,c,d)   if (!acceptq(#b"=",c,d,&(S->b),w)) return(0); sequence *b=&(S->b)
int acceptq (const char *name,const char *mode,const char *deflt,sequence *q, char *w);
#endif

