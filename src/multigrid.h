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

/* multigrid.h: the types of the multigrid solver of the elliptic device,
 * shared with other devices that use it. See multigrid_code.h for how to
 * use it. Include after diffweights.h. */

#ifndef MULTIGRID_H
#define MULTIGRID_H

/*==============================================*/
/* Important defined/derived types              */

/* Choice of the smoother */
typedef enum {
  Jacobi,
  GS
} smoother_t;

/* Choice of inter-grid transfer method */
typedef enum {
  correct,
  replace
} transfer_t;

/* Choice of the exact analytical solution.  */
/* Obviously, this is only for debugging !!! */
typedef enum {
  zero,
  sinv,
  ones,
  zfk
} exact_t;

/* Enumerated neighbours.    */
/* 1D precede 2D precede 3D. */
typedef enum {
  ooo,			/* 0D */
  poo, moo,		/* 1D 3-point stencil */
  opo, omo,		/* 2D 5-point stencil */
  ppo, pmo, mpo, mmo,	/* 2D 9-point stencil */
  oop, oom,		/* 3D 7-point stencil */
  pop, pom, mop, mom,	/* 3D 19-point stencil */
  opp, opm, omp, omm,	/* 3D 19-point stencil */
  ppp, ppm, pmp, pmm,	/* 3D 27-point stencil */
  mpp, mpm, mmp, mmm,	/* 3D 27-point stencil */
  maxnb
} neighbour_t;

/* "Systems of coordinates":									*/
typedef size_t coord_t;		/* 3D local (netto) x,y,z coordinates: x=0 <=> local_xmin etc   */
typedef ssize_t scoord_t;	/* Signed version of the same: for displacements		*/
typedef size_t index_t;		/* Flattened 1D array index: 0..xmax_ymax_zmax_vmax             */
typedef ssize_t sindex_t;	/* Signed version of the same: for displacements	        */
typedef size_t list_t;		/* Number of a point in the point list: numPoints         */

/* Conversions between coordinates:								*/
/* gridIndex(level,x,y,z):        coord_t^3 |-> index_t						*/
/*: (z+S->zmin[level]+S->nz[level]*(y+S->ymin[level]+S->ny[level]*(x+S->xmin[level]))) - same only longer  */
/* S->gridPoints[level][i].x etc: list_t    |-> coord_t					        */
/* pointIndex[level][i]:          list_t    |-> index_t						*/
/* gridX(level,ind) etc:          index_t   |-> coord_t				        	*/
/* pointNumber(level,ind):        index_t   |-> list_t						*/

/* Enumerated contributors.  */
/* 1D precede 2D precede 3D. */
typedef enum {
  OOO,			/* 0D */
  POO, MOO,		/* 1D */
  OPO, OMO,		
  PPO, PMO, MPO, MMO,	/* 2D */
  OOP, OOM,		
  POP, POM, MOP, MOM,	
  OPP, OPM, OMP, OMM,	
  PPP, PPM, PMP, PMM,
  MPP, MPM, MMP, MMM,	/* 3D */
  maxcn
} contr_t;

/* Enumerator of colors (as in multicolor GS) */
typedef int color_t; /* 0..7 actually */

/* An entry of the list of points, RMcF style */
typedef struct {
  real *K;				/* Neibours connection strengths (entries of A matrices) */
  coord_t x,y,z;			/* Coordinates */
} AnisoPoint;      

/* A link to value outside space */
typedef struct {
  index_t to;				/* Target point in the grid  */
  real *from;				/* Pointer to outside value  */
  real K;				/* Connection strength       */
} bclink;

#endif /* MULTIGRID_H */
