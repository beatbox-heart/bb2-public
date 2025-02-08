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

#define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#define G(u,v)          ( (u)-(v) )

/* ------------------------------------- 
 * Defining the kinetics macros:         
 * U_KINETICS(u,uth) and V_KINETICS(u,v) 
 * ------------------------------------- */

#if EXPLICIT   
     /* Explicit u-kinetics */

  #define U_KINETICS(u,uth) ( (u)+dt_o_eps*(u)*(one-(u))*((u)-(uth)) )

#else          
     /* DB: Implicit u-kinetics */
     /* The original (Physica 49D) scheme can be obtained by defining
      * F2m and F2p to be one */

  #define F1m(u,uth) ( dt_o_eps*(one-(u))*((u)-(uth)) )
  #define F1p(u,uth) ( dt_o_eps*(u)*((u)-(uth)) )
  #define F2m(u,uth) ( one + dt_o_2eps*((uth)-(u)*(u)) )
  #define F2p(u,uth) ( one + dt_o_2eps*(two*(u) -(uth)-(u)*(u)) )

  #define U_KINETICS(u,uth) (                                      \
          ( (u) < (uth) ) ?                                        \
          (u)/(one-F1m(u,uth)*F2m(u,uth) ) :                       \
          ((u)+F1p(u,uth)*F2p(u,uth))/(one+F1p(u,uth)*F2p(u,uth)) )

#endif

#define V_KINETICS(u,v) ( (v)+dt*G(u,v) )

/* ------------------------------------- 
   * Main Loop (almost all work done here) 
   * ------------------------------------- */
  for (z=s.z0;z<=s.z1;z++) {
    for (y=s.y0;y<=s.y1;y++) {
      index = ind(s.x0,y,z,0); 
      for (x=s.x0;x<=s.x1;x++) {
	if (u[index]<delta) {  
	  u[index] = zero;
	  v[index] = V_KINETICS(zero,v[index]);
	}
	else { 
	  u_thresh = U_THRESHOLD(v[index]);
	  v[index] = V_KINETICS(u[index],v[index]);
	  u[index] = U_KINETICS(u[index],u_thresh);
	}
	index+=DX;
      }
    }
  }
/*-------------------------------*/

#undef V_KINETICS
#undef F2p
#undef F2m
#undef F1p
#undef F1m
#undef U_KINETICS
#undef G
#undef U_THRESHOLD
