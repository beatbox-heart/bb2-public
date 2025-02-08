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

#include "extern.h"

#define REALLOC(p,a) if(0==((p)=Realloc((p),(a)))) ABORT("could not allocate %ld bytes for %s",((long)(a)),#p);
/* #define CALLOC(p,a,b) if(0==(p=Calloc(a,b))) ABORT("could not allocate %ld bytes for %s",((long)(a))*((long)(b)),#p); */
/* TODO: consider replacing this with the following, and then phasing out the third argument */
/* #define CALLOC(p,a,b) if(0==(p=Calloc(a,sizeof(*p)))) ABORT("could not allocate %ld bytes for %s",((long)(a))*((long)(sizeof(*p))),#p); */
/* or even */
#define CALLOC(p,a,b) {							\
  if (sizeof(*(p))!=(b)) ABORT("wrong size of *%s in calloc: expected %d, given %d",#p,sizeof(*(p)),(b)); \
  if (0==((p)=Calloc((a),sizeof(*(p))))) ABORT("could not allocate %ld bytes for %s",((long)(a))*((long)(sizeof(*(p)))),#p); \
}
#define MALLOC(p,a) if(0==((p)=Malloc(a))) ABORT("could not allocate %ld bytes for %s",((long)(a)),#p);

#include <stdlib.h>
#define FREE(addr) if(addr) {free(addr);(addr)=NULL;}

void * Calloc(size_t nitems,size_t size);
void * Malloc(size_t size);
void * Realloc(void *block, size_t size);
