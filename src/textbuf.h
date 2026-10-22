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

/* Growable text buffers, and texts of grid points collected in global order.	*/
/* 										*/
/* Used by devices that print text about many grid points, e.g. singz and	*/
/* localprint. Each process prints its points into a 'textpieces' set, one	*/
/* piece of text per point, tagged by a key: the point's global position.	*/
/* The root collects the sets of all processes, sorts the pieces by their keys	*/
/* and outputs them in that order, so the output does not depend on the domain	*/
/* decomposition and is the same as in the sequential mode.			*/

#ifndef _TEXTBUF_H_
#define _TEXTBUF_H_

#include <stddef.h>
#if MPI
#include <mpi.h>
#endif

/* Growable text buffer */
typedef struct {
  char *s;			/* the text, zero-terminated 	*/
  size_t n;			/* its length 			*/
  size_t max;			/* allocated size 		*/
} textbuf;

void txt_write(textbuf *b, const char *s, size_t n);	/* add n characters */
void txt_puts(textbuf *b, const char *s);		/* add a string */
void txt_printf(textbuf *b, const char *fmt, ...);	/* add formatted output */
void txt_free(textbuf *b);

/* A piece of text in a textbuf, for one grid point */
typedef struct {
  long long key;		/* global position of the point, defines the order of output */
  long long off;		/* offset of the text in the buffer */
  long long len;		/* its length */
} piece;

/* Pieces of text of several grid points */
typedef struct {
  textbuf txt;			/* their texts, concatenated */
  piece *pc; /*[max]*/		/* their keys and positions in txt */
  size_t n, max;		/* how many of them, and space allocated */
} textpieces;

void tp_reset(textpieces *p);			/* make empty, keeping the memory */
void tp_begin(textpieces *p, long long key);	/* start a new piece; then print it into p->txt */
void tp_end(textpieces *p);			/* the text of the current piece is complete */
void tp_sort(textpieces *p);			/* sort the pieces by their keys */
void tp_free(textpieces *p);

/* Key of grid point (x,y,z): pieces sorted by keys go z-major, x fastest */
long long tp_key(long x, long y, long z);
long tp_key_z(long long key);			/* z of the grid point with this key */

#if MPI
/* Concatenate the sets 'mine' of all processes of comm into 'all' at its rank 'root' */
int tp_gather(textpieces *mine, textpieces *all, int root, MPI_Comm comm);
#endif

#endif /* _TEXTBUF_H_ */
