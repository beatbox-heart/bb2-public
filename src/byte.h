/**
 * Copyright (C) (2010-2023) Vadim Biktashev, Irina Biktasheva et al. 
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

/* Functional version. Safer programming. */
unsigned Byte(int x, int y, int z, int v, real u0, real u1);

/* Macro version of the same. Better optimizing? (marginal speed-up) */
#define BYTE(x,y,z,v,u0,u1) ( \
  ((v)<0)?MAXCHAR \
  : (BYTEU=New[ind((x),(y),(z),(v))])<(u0)?0 \
  : BYTEU>(u1)? MAXCHAR \
  : MAXCHAR*(BYTEU-(u0))/((u1)-(u0)) \
)
