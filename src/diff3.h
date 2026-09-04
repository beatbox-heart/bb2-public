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

switch(dim) {
 case 1: {
#define _(x) (x),(0),(0)
#define STC STC1D(ANISO,MANY,VARIA)
   DB("1D version with GEOM=%d ANISO=%d MANY=%d VARIA=%d\n",GEOM,ANISO,MANY,VARIA);   
#include "diff1d.h"
#undef STC
#undef _
 } break; /* case 1 */
 case 2: {
#define _(x,y) (x),(y),(0)
#define STC STC2D(ANISO,MANY,VARIA)
   DB("2D version with GEOM=%d ANISO=%d MANY=%d VARIA=%d\n",GEOM,ANISO,MANY,VARIA);   
#include "diff2d.h"
#undef STC
#undef _
 } break; /* case 2 */
 case 3: {
#define _(x,y,z) (x),(y),(z)
#define STC STC3D(ANISO,MANY,VARIA)
   DB("3D version with GEOM=%d ANISO=%d MANY=%d VARIA=%d\n",GEOM,ANISO,MANY,VARIA);   
#include "diff3d.h"
#undef STC
#undef _
 } break; /* case 3 */
 default:
   EXPECTED_ERROR("unexpected dimension %d",dim);
   break;
 } /* switch dim */
