/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
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

#include "bgi.h"

typedef unsigned char colortype;

void SetWindow(BGIWindow W);
int SetLimits(real Absmin, real Absmax, real Ordmin, real Ordmax);
void Frame (Space s);
void Clean(void);
void Bar (int x, int y, colortype c);
void Bar1 (int x, int y, colortype c);
void Pixel (real x, real y, colortype color);
void Mark (real x, real y, colortype color, int size);
void Line(real x0, real y0, real x1, real y1, colortype color);
void MoveTo(real x, real y);
void LineTo(real x, real y, colortype color);
