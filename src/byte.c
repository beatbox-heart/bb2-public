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

/* Scale value at given grid point to unsigned char value */

#include "beatbox.h"
#include "state.h"
#include "byte.h"

#define MAXCHAR 255

unsigned Byte(int x, int y, int z, int v, real u0, real u1) {
  real u;
  if (v<0) return MAXCHAR;
  u=New[ind(x,y,z,v)];
  if (u<u0) return 0;
  if (u>u1) return MAXCHAR;
  return MAXCHAR*(u-u0)/(u1-u0);
}
