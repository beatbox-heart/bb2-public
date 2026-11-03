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

/* diffweights.h: the stencils of the "universal" diffusion device diff, shared
 * with other devices that need the same connection weights (e.g. elliptic).
 * Include before the device's STR typedef; then include diffweights_str.h
 * into the STR, and diffweights_code.h after it. */

#ifndef DIFFWEIGHTS_H
#define DIFFWEIGHTS_H

/* Handle tokens for CPP prescan */
#define JOIN2(a,b) a##b
#define JOIN3(a,b,c) a##b##c
#define JOIN4(a,b,c,d) a##b##c##d
#define APPLY(a,b) a b

/* Constants in stencils named as in popular formulas */
static const real sixth      = 1./6.;
static const real quarter    = 1./4.;
static const real twosixths  = 2./6.;
static const real half       = 1./2.;
static const real foursixths = 4./6.;

#ifndef UNIT_VECTOR_TOLERANCE
#define UNIT_VECTOR_TOLERANCE 0.01
#endif

/* Neighbourhoods and diff arrays descriptions 		*/
/* for different cases.					*/
/* ..nnb = # of neighbours,				*/
/* ..nD  = # of diff tensor components.			*/
/* NB: unique components listed lexically,		*/
/* then others added as synonyms, 			*/
/* for code convenience and clarity.			*/
/* STD?D??? macros will be used in diff3.h.		*/
/* 3-bit list of stencils is excessive in some cases.	*/

/* 1+3; 1 : 1D iso (can't be otherwise in 1D!) */
enum i1_nbrs {i10, i1p, i1m, i1nnb};
enum i1_D    {i1D11, i1nD};
#define STC1D000 i1
#define STC1D001 i1  
#define STC1D010 i1  
#define STC1D011 i1
#define STC1D100 i1  
#define STC1D101 i1  
#define STC1D110 i1  
#define STC1D111 i1
#define STC1D(a,m,v) JOIN4(STC1D,a,m,v)

/* 1+5; 1 : 2D iso */
enum i2_nbrs {i200, i2p0, i2m0, i20p, i20m, i2nnb};
enum i2_D    {i2D11, i2nD, i2D22=i2D11};
/* 1+9; 1 : 2D manypoint */
enum m2_nbrs {m200, m2p0, m2m0, m20p, m20m,
	      m2pp, m2pm, m2mp, m2mm, m2nnb};
enum m2_D    {m2D11, m2nD, m2D22=m2D11};
/* 1+9; 3 : 2D aniso */
enum a2_nbrs {a200, 
              a2p0, a2m0, a20p, a20m,
	      a2pp, a2pm, a2mp, a2mm, a2nnb};
enum a2_D    {a2D11, a2D12, a2D22, a2nD, a2D21=a2D12};
#define STC2D000 i2
#define STC2D001 i2  
#define STC2D010 m2  
#define STC2D011 m2
#define STC2D100 a2  
#define STC2D101 a2  
#define STC2D110 a2  
#define STC2D111 a2
#define STC2D(a,m,v) JOIN4(STC2D,a,m,v)

/* 1+7; 1 : 3D iso */
enum i3_nbrs {i3000, i3p00, i3m00, i30p0, i30m0, i300p, i300m, i3nnb};
enum i3_D    {i3D11, i3nD, i3D22=i3D11, i3D33=i3D11};
/* 1+19, 1 : 3D manypoint */
enum m3_nbrs {m3000, m3p00, m3m00, m30p0, m30m0, m300p, m300m,
	      m3pp0, m3pm0, m3mp0, m3mm0, m3p0p, m3p0m,
	      m3m0p, m3m0m, m30pp, m30pm, m30mp, m30mm, m3nnb};
enum m3_D    {m3D11, m3nD, m3D22=m3D11, m3D33=m3D11};
/* 1+19, 6 : 3D aniso */
enum a3_nbrs {a3000, a3p00, a3m00, a30p0, a30m0, a300p, a300m,
              a3pp0, a3pm0, a3mp0, a3mm0, a3p0p, a3p0m,
	      a3m0p, a3m0m, a30pp, a30pm, a30mp, a30mm, 
	      a3nnb};
enum a3_D    {a3D11, a3D12, a3D13, a3D22, a3D23, a3D33, a3nD, a3D21=a3D12, a3D31=a3D13, a3D32=a3D23};
#define STC3D000 i3
#define STC3D001 i3  
#define STC3D010 m3  
#define STC3D011 m3
#define STC3D100 a3  
#define STC3D101 a3  
#define STC3D110 a3  
#define STC3D111 a3
#define STC3D(a,m,v) JOIN4(STC3D,a,m,v)

/* Displacements (dx,dy,dz) of the neighbours of the stencil for the given
   anisotropy flag, manypoint flag and dimension, in the order of its enum
   (e.g. a3_nbrs); returns the number of neighbours nnb. */
static inline int diff_stencil_offsets(int aniso, int many, int dimension, int off[a3nnb][3])
{
  /* the largest stencils of each dimension list the smaller ones first */
  static const int o3[a3nnb][3]={
    { 0, 0, 0},
    {+1, 0, 0},{-1, 0, 0},{ 0,+1, 0},{ 0,-1, 0},{ 0, 0,+1},{ 0, 0,-1},
    {+1,+1, 0},{+1,-1, 0},{-1,+1, 0},{-1,-1, 0},
    {+1, 0,+1},{+1, 0,-1},{-1, 0,+1},{-1, 0,-1},
    { 0,+1,+1},{ 0,+1,-1},{ 0,-1,+1},{ 0,-1,-1}};
  static const int o2[a2nnb][3]={
    { 0, 0, 0},
    {+1, 0, 0},{-1, 0, 0},{ 0,+1, 0},{ 0,-1, 0},
    {+1,+1, 0},{+1,-1, 0},{-1,+1, 0},{-1,-1, 0}};
  const int (*o)[3];
  int nnb, nb, k;
  switch (dimension) {
  case 1: o=o3; nnb=i1nnb; break;		/* the first 3 of o3 */
  case 2: o=o2; nnb=(aniso||many)?a2nnb:i2nnb; break;
  case 3: o=o3; nnb=(aniso||many)?a3nnb:i3nnb; break;
  default: return 0;
  }
  for (nb=0;nb<nnb;nb++) for (k=0;k<3;k++) off[nb][k]=o[nb][k];
  return nnb;
} /* diff_stencil_offsets */

#endif /* DIFFWEIGHTS_H */
