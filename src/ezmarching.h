/**
 * Copyright (C) (2010-2018) Vadim Biktashev, Irina Biktasheva et al. 
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

/* Modifications for ezview: 2010-18, Vadim Biktashev v.n.biktashev@exeter.ac.uk */

#ifndef _EZMARCHING_
#define _EZMARCHING_

/* ------------------------------------------------------------------------- 
 * DB: All floating point numbers in ezmarching are of type GLReal. Defined here
 * to be GLfloat. If the precision of the numerical data is double then
 * implicit conversion will take place when calculating the contour intersect
 * and finding the cube index. You may wish to experiment and try setting
 * this to double or to Real.  The functions glVertex* and glColor* functions
 * have many forms so I use macros to make changes easier. 
 *
 * Note that we often need to define points (and vectors) in 3 space
 * dimensions. These are used both in calculations and as data to pass to
 * OpenGL for rendering.  All such variables are defined as arrays of size
 * 3. In cases where we need arrays of coordinates, the index for each
 * dimension should be the last one. This means that the X, Y and Z
 * components are contiguous in memory and enables us to use glVertex3fv
 * instead of glVertex3f. 
 * ------------------------------------------------------------------------- */

/* 
 * VNB: some typedefs and macros moved to ezview.h as they are needed to describe
 * "quasi-global" variables, i.e. fields of the STR structure 
*/

/* typedef GLfloat       GLReal; -> ezview.h */

#define GLVERTEX3V(v)     glVertex3fv(v)
#define GLVERTEX3(x,y,z)  glVertex3f((x),(y),(z))
#define GLCOLOR3(r,g,b)   glColor3f((r),(g),(b))
#define GLCOLOR4(r,g,b,a) glColor4f((r),(g),(b),(a))
#define GLNORMAL3(x,y,z)  glNormal3f((x),(y),(z))
#define GLNORMAL3V(v)     glNormal3fv(v)

/* ------------------------------------------------------------------------- 
 *
 *            Probably you do not want to change anything below              
 *
 * ------------------------------------------------------------------------- */

/* typedef unsigned char CubeEdge; -> ezview.h */
typedef unsigned int  CubeVertex;     /* Numer between 1 and 8. */
typedef unsigned int  CubeIndex;      /* Number between 0 and 255. Could use
				       * unsigned char, but storage is not an
				       * issue and the machine is probably
				       * better at multiplying ints. Some
				       * machines are better at multiplying
				       * unsigned ints than signed ints. */
typedef unsigned int  Axis;           /* Either X, Y or Z. */


/* #define MAX_TRIANGLE_LIST_LENGTH -> ezview.h */
#define NUM_GENERATORS            21  /* Number of basic marching cubes cases
				       * which can be used to build the rest
				       * of the table. 1987 paper says 15, but
				       * with 1995 modification this goes up
				       * to 21 to deal with complimentary
				       * cases. */
#define END                        0  /* End of edge list marker. */

#define EQUAL_ZERO(n)  ((n) == 0.0)   /* In this code, this should be
				       * ok because we always clamp
				       * values afterwards and we're
				       * not too * concerned about
				       * underflow and overflow. */

#define X 0                           /* Indices for vector components */
#define Y 1
#define Z 2

#define A 0                           /* Params for implicit form of plane. */
#define B 1
#define C 2
#define D 3


/* Vector operations used exclusively by the functions filament_in_cube() and
 * find_two_pts () which finds intersections of u and v triangles. */

#define SUB(ans,v1,v2) \
  (ans)[X] = (v1)[X] - (v2)[X]; \
  (ans)[Y] = (v1)[Y] - (v2)[Y]; \
  (ans)[Z] = (v1)[Z] - (v2)[Z]

#define DOT(v1, v2) ((v1)[X]*(v2)[X] + (v1)[Y]*(v2)[Y] + (v1)[Z]*(v2)[Z])

#define POINT_ALONG_LINE(ans, start, length, direction) \
  (ans)[X] = (start)[X] + (length)*(direction)[X]; \
  (ans)[Y] = (start)[Y] + (length)*(direction)[Y]; \
  (ans)[Z] = (start)[Z] + (length)*(direction)[Z]

#endif  /* _EZMARCHING_ */
