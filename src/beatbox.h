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

#define REALF_IO "%38f" // Specifying precision to ensure exact length for I/O types.
#define REALF_IO_WIDTH 38

#define REALF "%lg"

#define MAXREAL DBL_MAX

/* return codes */
#define SUCCESS 1
#define FAILURE 0
#define MARK_ERROR return(FAILURE);

int nofflush(void *f);

#include "dynamic.h"

#endif

