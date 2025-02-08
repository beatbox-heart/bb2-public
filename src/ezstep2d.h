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

/* DB: I rely on the compiler to take care of most optimization.  The only place
 * where I give the compiler help is with the Laplacian formulas because
 * none of the compiles do the obvious thing.  */

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

/* ---------------------------------------- 
 * Defining the diffusion macros:           
 * U_DIFFUSION, V_DIFFUSION, ZERO_USED_SUMS 
 * ---------------------------------------- */

#if V_DIFF_ON        
     /* v is diffusing */
  #define U_DIFFUSION(index_1)     (dt_o_wh2   * sigma_u[index_1]) 
  #define V_DIFFUSION(index_1)     (dtDv_o_wh2 * sigma_v[index_1]) 
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = sigma_v[index_1] = zero;
#else
     /* v is not diffusing */
  #define U_DIFFUSION(index_1)     (dt_o_wh2 * sigma_u[index_1]) 
  #define V_DIFFUSION(index_1)     zero
  #define ZERO_USED_SUMS(index_1)  sigma_u[index_1] = zero;
#endif

/* -------------------------------- 
 * Defining the spatial-sum macros: 
 * ADD_TO_U_SUM, ADD_TO_V_SUM     
 * -------------------------------- */

/* 2D 5-point Laplacian formula */
#define ADD_TO_U_SUM5(index,index_2)		\
{ 						\
  real stupid_cc = u[index];			\
  sigma_u[(index_2)]  -= four*stupid_cc;	\
  sigma_u[(index_2)+J_INC] += stupid_cc;	\
  sigma_u[(index_2)-J_INC] += stupid_cc;	\
  sigma_u[(index_2)+I_INC] += stupid_cc;	\
  sigma_u[(index_2)-I_INC] += stupid_cc;	\
}
#define ADD_TO_V_SUM5(index,index_2) 		\
{						\
  real stupid_cc = v[index];		     	\
  sigma_v[(index_2)]  -= four*stupid_cc;	\
  sigma_v[(index_2)+J_INC] += stupid_cc;	\
  sigma_v[(index_2)-J_INC] += stupid_cc;	\
  sigma_v[(index_2)+I_INC] += stupid_cc;	\
  sigma_v[(index_2)-I_INC] += stupid_cc;	\
}
/* 2D 9-point Laplacian formula */
#define ADD_TO_U_SUM9(index,index_2)		\
{ 						\
  real stupid_cc = u[index];		      	\
  sigma_u[(index_2)]  -= twenty*stupid_cc;	\
  sigma_u[(index_2)+J_INC]  += four*stupid_cc;	\
  sigma_u[(index_2)-J_INC]  += four*stupid_cc;	\
  sigma_u[(index_2)+I_INC]  += four*stupid_cc;	\
  sigma_u[(index_2)-I_INC]  += four*stupid_cc;	\
  sigma_u[(index_2)+I_INC+J_INC] += stupid_cc;	\
  sigma_u[(index_2)+I_INC-J_INC] += stupid_cc;	\
  sigma_u[(index_2)-I_INC+J_INC] += stupid_cc;	\
  sigma_u[(index_2)-I_INC-J_INC] += stupid_cc;	\
}
#define ADD_TO_V_SUM9(index,index_2) 		\
{						\
  real stupid_cc = v[index];			\
  sigma_v[(index_2)]  -= twenty*stupid_cc;	\
  sigma_v[(index_2)+J_INC]  += four*stupid_cc;	\
  sigma_v[(index_2)-J_INC]  += four*stupid_cc;	\
  sigma_v[(index_2)+I_INC]  += four*stupid_cc;	\
  sigma_v[(index_2)-I_INC]  += four*stupid_cc;	\
  sigma_v[(index_2)+I_INC+J_INC] += stupid_cc;	\
  sigma_v[(index_2)+I_INC-J_INC] += stupid_cc;	\
  sigma_v[(index_2)-I_INC+J_INC] += stupid_cc;	\
  sigma_v[(index_2)-I_INC-J_INC] += stupid_cc;	\
}

#if MANYPOINT 	/* 9-point Laplacian formula */
  #define ADD_TO_U_SUM(index,index_2) ADD_TO_U_SUM9(index,index_2)
#else		/* 5-point Laplacian formula */
  #define ADD_TO_U_SUM(index,index_2) ADD_TO_U_SUM5(index,index_2)
#endif         

#if !V_DIFF_ON	/* v is not diffusing */
  #define ADD_TO_V_SUM(index,index_2) 
#else		/* v is diffusing */
#if MANYPOINT 	/* 9-point Laplacian formula */
  #define ADD_TO_V_SUM(index,index_2) ADD_TO_V_SUM9(index,index_2)
#else		/* 5-point Laplacian formula */
  #define ADD_TO_V_SUM(index,index_2) ADD_TO_V_SUM5(index,index_2) 
#endif 
#endif

/* ------------------------------------- 
   * Main Loop (almost all work done here) 
   * ------------------------------------- */
  for (y=s.y0;y<=s.y1;y++) {
    index = ind(s.x0,y,0,0); 
    index_1 = index + (*s1)*DV;
    index_2 = index + (*s2)*DV;
    for (x=s.x0;x<=s.x1;x++) {
#if SPLIT        /* (split=1 is best for large time steps, else use split=0) */
      u[index] += U_DIFFUSION(index_1);
      v[index] += V_DIFFUSION(index_1);
      ZERO_USED_SUMS(index_1);
#endif

      if(u[index]<delta) {  
	u[index] = zero;
	v[index] = V_KINETICS(zero,v[index]);
#if !SPLIT
	u[index] += U_DIFFUSION(index_1);
	v[index] += V_DIFFUSION(index_1);
#endif
      }
      else { 
	u_thresh = U_THRESHOLD(v[index]);
	v[index] = V_KINETICS(u[index],v[index]);
	u[index] = U_KINETICS(u[index],u_thresh);
#if !SPLIT
	u[index] += U_DIFFUSION(index_1);
	v[index] += V_DIFFUSION(index_1);
#endif
	ADD_TO_U_SUM(index,index_2);
      }
      ADD_TO_V_SUM(index,index_2);
      ZERO_USED_SUMS(index_1);
      
      index+=DX; index_1+=DX; index_2+=DX;
    }
  }
/*-------------------------------*/

#undef ADD_TO_V_SUM
#undef ADD_TO_U_SUM
#undef ZERO_USED_SUMS
#undef V_DIFFUSION
#undef U_DIFFUSION
#undef V_KINETICS
#undef F2p
#undef F2m
#undef F1p
#undef F1m
#undef U_KINETICS
#undef G
#undef U_THRESHOLD
