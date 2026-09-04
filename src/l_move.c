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

/* Shift layer in space by a given (dx,dy,dz,dv) vector, */
/* similar to Mover from Barkley EZSpiral                */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "qpp.h"

typedef struct {
  int v0;				/* bottom source layer */
  int v1;				/* top source layer    */
  int dv;				/* layers' shift       */
  PAR(long,dx);				/* components          */
  PAR(long,dy);				/*   of shift          */
  PAR(long,dz);				/*   vector            */
} STR;

static void Mover(STR *S,Space s,int m_x,int m_y,int m_z,int m_v);

/****************/
RUN_HEAD(l_move)
{
  DEVICE_PAR(long,dx);
  DEVICE_PAR(long,dy);
  DEVICE_PAR(long,dz);
  DEVICE_CONST(int,dv);
  int r_x=dx;				/* movement         */
  int r_y=dy;				/*   that remains   */
  int r_z=dz;				/*   to be done     */
  int m_x, m_y, m_z;			/* current move     */
  while (r_x!=0 || r_y!=0 || r_z!=0) {
    m_x=(r_x>0)?1:(r_x<0)?(-1):0;	/* move at most     */
    m_y=(r_y>0)?1:(r_y<0)?(-1):0;	/*   1 step in      */
    m_z=(r_z>0)?1:(r_z<0)?(-1):0;	/*   space 	    */
    Mover(S,s,m_x,m_y,m_z,dv);		/*   at a time      */
    r_x-=m_x;				/* what 	    */
    r_y-=m_y;				/*   remains        */
    r_z-=m_z; 				/*   to do          */
  }
}
RUN_TAIL(l_move)

/* Variation upon Barkley's theme. */
/* Todo: incorporate and optimize as this is called only once */
static void Mover(STR *S,Space s,int m_x,int m_y,int m_z,int m_v)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  int     i0, i1, i_start, i_stop, i_step;
  int     j0, j1, j_start, j_stop, j_step;
  int     k0, k1, k_start, k_stop, k_step;
  int     l0, l1, l_start, l_stop, l_step;

  if (m_x==0 && m_y==0 && m_z==0) return;
  /* Debug("m_x=%d m_y=%d m_z=%d\n",m_x,m_y,m_z); */

  /* these checks are redundant, correct usage ensured by hand */
  /* ASSERT(m_x<=1);  ASSERT(m_x>=-1); */
  /* ASSERT(m_y<=1);  ASSERT(m_y>=-1); */
  /* ASSERT(m_z<=1);  ASSERT(m_z>=-1); */

  if (m_x > 0) { i_start = s.x1; i_stop = s.x0; i_step = -1;}
          else { i_start = s.x0; i_stop = s.x1; i_step = 1; }
  if (m_y > 0) { j_start = s.y1; j_stop = s.y0; j_step = -1;}
          else { j_start = s.y0; j_stop = s.y1; j_step = 1; }
  if (m_z > 0) { k_start = s.z1; k_stop = s.z0; k_step = -1;}
          else { k_start = s.z0; k_stop = s.z1; k_step = 1; }
  if (m_v > 0) { l_start =   v1; l_stop =   v0; l_step = -1;}
          else { l_start =   v0; l_stop =   v1; l_step = 1; }

#if MPI
  haloSwap();
#endif
  
  /* Debug("i:%d:%d:%d j:%d:%d:%d k:%d:%d:%d\n", */
  /* 	i_start,i_stop,i_step, */
  /* 	j_start,j_stop,j_step, */
  /* 	k_start,k_stop,k_step); */

  for (i0=i_start; i0!=(i_stop+i_step); i0+=i_step) {
    for (j0=j_start; j0!=(j_stop+j_step); j0+=j_step) {
      for (k0=k_start; k0!=(k_stop+k_step); k0+=k_step) {
	for (l0=l_start; l0!=(l_stop+l_step); l0+=l_step) {
	  i1=i0+m_x;
	  j1=j0+m_y;
	  k1=k0+m_z;
	  l1=l0+m_v;
	  if ( isTissue(i0,j0,k0) && isTissue(i1,j1,k1))
	    New[ind(i1,j1,k1,l1)]=New[ind(i0,j0,k0,l0)];
	} /* for l0 */
      } /* for k0 */
    } /* for j0 */
  } /* for i0 */
} /* Mover */


DESTROY_HEAD(l_move)
DESTROY_TAIL(l_move)

CREATE_HEAD(l_move)
{
  DEVICE_IS_RECTANGULAR;
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,v0,v0,(int)vmax-1);
  ACCEPTI(dv,0,-v0,vmax-1-v1);
  if (ONE) {
    ACCEPTIK(dx,INONE,-xmax+1,xmax-1);
  } else {
    ACCEPTInoK(dx,0,0,0);
  }
  if (TWO) {
    ACCEPTIK(dy,INONE,-ymax+1,ymax-1);
  } else {
    ACCEPTInoK(dy,0,0,0);
  }
  if (TRI) {
    ACCEPTIK(dz,INONE,-zmax+1,zmax-1);
  } else {
    ACCEPTInoK(dz,0,0,0);
  }
}
CREATE_TAIL(l_move,1)

