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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h> /* needed to print error messages */
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#define OWN
#include "windraw.h"

			    /* static int debug=0; */
			    #define DB if (debug) fprintf(debug,
			    #define BD );  
					

/* extern int update_now; */
static BGIWindow w;
static real absmin, absmax, ordmin, ordmax;

static real _x, _y, _z;
#define _X(x) (((w.col1-1)*((x)-absmin)+(w.col0+1)*(absmax-(x)))/(absmax-absmin))
#define _Y(y) (((w.row0+1)*((y)-ordmin)+(w.row1-1)*(ordmax-(y)))/(ordmax-ordmin))
/*#define _Y(y) (((w.row1-1)*((y)-ordmin)+(w.row0+1)*(ordmax-(y)))/(ordmax-ordmin))*/
#define _Z(z) (((w.row0+1)*((z)-(int)zmin)+(w.row1-1)*((int)zmax-(z)))/((int)zmax-(int)zmin))
/* #define d 1 */
static real cropX(real x,BGIWindow w) {
  real _x=_X(x);
  if (_x<(w.col0+1)) return (w.col0+1);
  if (_x>(w.col1-1)) return (w.col1-1);
  return _x;
}
#define X(x) cropX(x,w)
static real cropY(real y,BGIWindow w) {
  real _y=_Y(y);
  if (_y<(w.row0+1)) return (w.row0+1);
  if (_y>(w.row1-1)) return (w.row1-1);
  return _y;
}
#define Y(y) cropY(y,w)
static real cropZ(real z,BGIWindow w,int zmin,int zmax) {
  real _z=_Z(z);
  if (_z<(w.row0+1)) return (w.row0+1);
  if (_z>(w.row1-1)) return (w.row1-1);
  return _z;
}
#define Z(z) cropZ(z,w,zmin,zmax)

extern int ndev;		/* num of devices in current run */
extern Device dev[MAXDEV];	/* the array of devices */

void SetWindow(BGIWindow W) {w=W;}

int SetLimits(real Absmin, real Absmax, real Ordmin, real Ordmax){
  ASSERT(Absmin!=Absmax);
  ASSERT(Ordmin!=Ordmax);
  absmin=Absmin;absmax=Absmax;
  ordmin=Ordmin;ordmax=Ordmax;
  return 1;
}

/*********************************************************/
void Frame(Space s) {
  int absmin=0, ordmin=0, zmin=0;
  if (!graphon) return;
  setwritemode(0);
  setcolor(w.color/16);
  /*rectangle(X(absmin),Y(ordmin),X(absmax),Y(ordmax));*/
  rectangle(w.col0,w.row0,w.col1,w.row1);
  setcolor(w.color%16);
  switch(w.area) {
  case 0: break;
  case 1: rectangle(X(s.x0),Y(s.y0),X(s.x1+1)+1,Y(s.y1+1)-1); break;
  case 2: rectangle(X(s.x0),Z(s.z0),X(s.x1+1)+1,Z(s.z1+1)-1); break;
  }
  /* update_now=1; */
}

void Clean(void) {
  if (!graphon) return;
  setfillstyle(SOLID_FILL,BLACK);
  bar(w.col0+1,w.row0+1,w.col1-1,w.row1-1);
  /* update_now=1; */
}

void Bar (int x, int y, colortype c) {
  int absmin=0, ordmin=0, zmin=0;
  if (!graphon) return;
  setfillstyle(SOLID_FILL,c);
  switch(w.area) {
  case 0: break;
  case 1: bar(X(x)+1,Y(y)-1,X(x+1),Y(y+1)); break;
  case 2: bar(X(x)+1,Z(y)-1,X(x+1),Z(y+1)); break;
  }
  /* update_now=1; */
}

void Bar1 (int x,int y, colortype c) {
  if (!graphon) return;
  setfillstyle(SOLID_FILL,c);
  bar(X(x),Y(y+1),X(x+1),Y(y));
  /* update_now=1; */
}

void  _putpixel(int x, int y, int color) {
  putpixel(x,y,color);
}

void Pixel (real x, real y, colortype color) {
  if (!graphon) return;
   _putpixel(X(x),Y(y),color);
}

				#include <stdio.h>

void Mark (real x, real y, colortype color, int size) {
  if (!graphon) return;
				DB "Mark %f %f %d %d\n",x,y,color,size BD
  setfillstyle(SOLID_FILL,color);
  setcolor(color);
  fillellipse(X(x),Y(y),size,size);
}

void _line(int x1, int y1, int x2, int y2) {
  line(x1,y1,x2,y2);
}

void Line(real x0, real y0, real x1, real y1, colortype color) {
  if (!graphon) return;
  setcolor(color);
  _line(X(x0),Y(y0),X(x1),Y(y1));
}

void MoveTo(real x, real y) {
	if (!graphon) return;
	DB "MoveTo %f %f \n",x,y BD
	moveto(X(x),Y(y));
}

void LineTo(real x, real y, colortype color) {
  if (!graphon) return;
				DB "LineTo %f %f %d\n",x,y,color BD
  setcolor(color);
  lineto(X(x),Y(y));
}
