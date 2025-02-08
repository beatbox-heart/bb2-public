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

#ifndef _BEATBOX_H_
#define _BEATBOX_H_ 

/* Most important common definitions */

#define VERLEN 16
#define VARLEN 16

#ifdef MAIN
char VERSTRING[VERLEN]=PACKAGE_STRING;
#if MPI
	char VARIATION[VARLEN]="MPI";
#else
	char VARIATION[VARLEN]="sequential";
#endif
#else
extern char VERSTRING[];
extern char VARIATION[];
extern int Verbose;
extern int idev;
#endif

#define MAXDEPTH 10

#define real double
/* #if MPI */
/* #ifdef MPI_REAL */
/* #undef MPI_REAL */
/* #endif */
/* #define MPI_REAL MPI_DOUBLE */
/* #endif */

#define REALF_IO "%38f" // Specifying precision to ensure exact length for I/O types.
#define REALF_IO_WIDTH 38

#define REALF "%lg"

#define MAXREAL DBL_MAX

/* return codes */
#define SUCCESS 1
#define FAILURE 0
#define MARK_ERROR return(FAILURE);

void ANY_MESSAGE(int urgent,char *fmt,...);
#define MESSAGE(...) ANY_MESSAGE(0,__VA_ARGS__)
#define MESSAGE0(s) MESSAGE(s)
#define MESSAGE1(s,a1) MESSAGE(s,a1)
#define MESSAGE2(s,a1,a2) MESSAGE(s,a1,a2)
#define MESSAGE3(s,a1,a2,a3) MESSAGE(s,a1,a2,a3)
#define MESSAGE4(s,a1,a2,a3,a4) MESSAGE(s,a1,a2,a3,a4)
#define URGENT_MESSAGE(...) ANY_MESSAGE(1,__VA_ARGS__)

void Debug(char *fmt, ...);
#if MPI
#define DEBUG(...)			\
  if (debug) {				\
    fprintf(debug,"#%d %s:%d t=%ld idev=%d ",mpi_rank,__FILE__,__LINE__,t,idev);\
    fprintf(debug, __VA_ARGS__);	\
    fflush(debug);			\
  }
#else
#define DEBUG(...)			\
  if (debug) {				\
    fprintf(debug,"%s:%d t=%ld idev=%d ",__FILE__,__LINE__,t,idev);	\
    fprintf(debug, __VA_ARGS__);	\
    fflush(debug);			\
  }
#endif

int nofflush(void *f);

#include "dynamic.h"

#endif

