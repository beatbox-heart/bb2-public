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

/* diffweights_str.h: the fields of a device's STR used by the connection weights
 * code of diffweights_code.h; include inside the STR typedef. */

  p_tb loctb;				/* Local symtable; NULL if none of PARs are k-expressions */
  real *u;				/* Copy of local u vector; NULL if none of PARs depends on it */
  real *geom;				/* ..., geom vector, ... */
  INT x, y, z;				/* Copy of current coords for (non-NULL only) loctb */
  PAR(real,Dl);				/* Parallel diffusion coefficient (ANISO only) */
  PAR(real,Dt);				/* Transverse diffusion coefficient (ANISO only) */
  PAR(real,D);				/* Scalar diffusivity for isotropic or singular anisotropic */
  real hx;				/* Space step */
  int manypoint; 			/* == NINEPOINT, NINETEENPT in ezdiff < ezspiral, ezscroll */
  int conservative;			/* Whether to use conservative scheme */
  int harmonic;				/* .. conservative scheme is harmonic rather than arithmetic */
  int nnb;				/* Num of neighbours in the stencil */
  int w0, w1;				/* Layers range to keep the stencil weights */
  int nD;				/* Num of diff tensor components */
  int d0, d1;				/* Layers range to keep the diff tensor; must for var diff in MPI */
  int check_positivity;			/* Flag to test for "bad" points in anisotropy */
  int isolated;				/* Whether border of device space is no-flux, or connects to tissue beyond */
  size_t nadjacent;			/* Num of points of adjacent tissue, see diff_adjacent() */
  real *Xpool;				/* Weights of all points if not in layers w0..w1 */
  real *Ypool;				/* Diff tensors of all points if not in layers d0..d1 */
