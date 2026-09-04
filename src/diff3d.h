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

/***************************************************************/
/* 3D diffusion                                                */
/***************************************************************/

real Gam=1.0/(hx*hx);			//DB("Gam=%lg\n",(double)Gam); /* prescale */
real f1, f2, f3;			/* fibre cosines */
real Dlh, Dth, Ddh, Dsh, dW;		/* convenience intermediates */
real dD11, dD12, dD13, dD21, dD22, dD23, dD31, dD32, dD33; /* diff tensor gradient */
real c1, c2, c3;			/* .., contraction: convection vector */
/* */					//DB("loctb=%p\n",S->loctb);

DB("3D version: GEOM=%d ANISO=%d MANY=%d VARIA=%d\n", GEOM, ANISO, MANY, VARIA);

#if !VARIA
  k_on();
#if ANISO
  DEVICE_PAR(real,Dl);			/* longitudinal */	//DB("Dl=%g\n",Dl);
  DEVICE_PAR(real,Dt);			/* transversal */	//DB("Dt=%g\n",Dt);
  real D=(S->D)=			/* scalar for exceptional points */
    (S->Dptr)?(*(real *)(S->Dptr)):
    (S->Dexe)?(*(real *)execute(S->Dexe)):
    ( pow(Dl*Dt*Dt,1.0/3.0) );		//DB("D=%g\n",D);
  Dlh=Dl*Gam;				/* scaled longitudinal */
  Dth=Dt*Gam;				/* scaled transversal */
  Ddh=Dlh-Dth;				/* scaled difference */
  Dsh=D*Gam;				/* scaled scalar */
  /* */					//DB("Dl=%g D1t=%g D=%g Dth=%g Ddh=%g Dsh=%g\n",Dl,Dt,D,Dth,Ddh,Dsh);
#else /* not ANISO */
  DEVICE_PAR(real,D);			/* scalar diffusivity only */
  Dsh=D*Gam;				/* scaled scalar */
#endif /* not ANISO */
  k_off();
#endif /* if !VARIA */

/*---------------------------------------------------*/
/* First round: the diff tensor and ipgrid signature */
/* */					DB("first round\n");
for (ip=0;ip<np;ip++) { 
  P=s->p[ip];
  x=P.x;
  y=P.y;
  z=P.z;
  CLEAR_WEIGHTS;
#if !MPI
  IPNBR(x,y,z)=ip; 			/* count this point in the ipgrid */
#endif
  /* */					//DB("ip=%ld xyz=%ld,%ld,%ld Gam=%g\n",ip,x,y,z,Gam);
#if VARIA
  if (S->loctb) {			/* yes local dependencies */
    k_on();
    S->x=x;				/* copy point coords */
    S->y=y;				/*        into local */
    S->z=z;				/*          symtable */
    if (S->u) memcpy(S->u,P.u,vmax*sizeof(real)); /* ditto, vector of u, if relevant */
#if GEOM
    if (S->geom) memcpy(S->geom,&(G(0,x,y,z)),geom_vmax*sizeof(real)); /* ditto, vector of geom */
#endif
#if ANISO
    DEVICE_PAR(real,Dl);		/* longitudinal */	//DB("Dl=%g\n",Dl);
    DEVICE_PAR(real,Dt);		/* transversal */	//DB("Dt=%g\n",Dt);
    real D=S->D=			/* scalar for exceptional points */
      (S->Dptr)?(*(real *)(S->Dptr)):
      (S->Dexe)?(*(real *)execute(S->Dexe)):
      ( pow(Dl*Dt*Dt,1.0/3.0) );	/* default if not given in input file */ //DB("D=%g\n",D);
    Dlh=Dl*Gam;				/* scaled longitudinal */
    Dth=Dt*Gam;				/* scaled transversal */
    Ddh=Dlh-Dth;			/* scaled difference */
    Dsh=D*Gam;				/* scaled scalar */
    /* */				//DB("Dl=%g D1t=%g D=%g Dth=%g Ddh=%g Dsh=%g\n",Dl,Dt,D,Dth,Ddh,Dsh);
#else /* not ANISO */
    DEVICE_PAR(real,D);			/* scalar diffusivity only */
    Dsh=D*Gam;				/* scaled scalar */
#endif /* not ANISO */
    k_off();
  } /* if loctb */
#endif /* VARIA */
  /* Now compute the tensor components */
#if ANISO
  f1=G(1,x,y,z);			/* fiber     */
  f2=G(2,x,y,z);			/* direction */
  f3=G(3,x,y,z);			/*   cosines */
  /* */					//DB("f=%g,%g,%g\n",f1,f2,f3);
  if (isUnitVector(f1,f2,f3,UNIT_VECTOR_TOLERANCE)) {
    D(11) = Ddh*f1*f1 + Dth;
    D(22) = Ddh*f2*f2 + Dth;
    D(33) = Ddh*f3*f3 + Dth;
    D(12) = Ddh*f1*f2;
    D(13) = Ddh*f1*f3;
    D(23) = Ddh*f2*f3;
  } else { /* not isUnitVector */
    D(11) = Dsh;
    D(22) = Dsh;
    D(33) = Dsh;
    D(12) = 0.0;
    D(13) = 0.0;
    D(23) = 0.0;
  } //DB("D[][]=\n%g\t%g\t%g\n%g\t%g\t%g\n%g\t%g\t%g\n\n",D(11),D(12),D(13),D(12),D(22),D(23),D(13),D(23),D(33));
#else /* not ANISO */
  D(11) = Dsh;				/* ==D(22)==D(33) */
#endif /* not ANISO */  
 } /* for ip */

#if (VARIA || ANISO) && MPI
/* Variable diffusivity needs neighbours' tensors */
/* TODO: make sure this is done only if needed    */
/* TODO: make haloSwap subarrays distinct for every device that might need it, */
/*                             and include into them only the relevant layers! */
/* TODO: consider making exchange buffers smaller by selecting only needed points ??? */
ASSERT(d0>=0);			/* here diff tensor has to be on the grid */
haloSwap();			/* make it available for neighbouring subdomains */
#endif /* (VARIA || ANISO) && MPI */

/*---------------------------------*/
/* Second round: periphery weights */
/* */					// DB("second round\n");
for (ip=0;ip<np;ip++) {
  P=s->p[ip];
  x=P.x;
  y=P.y;
  z=P.z;
  
  /*=============================================*/
#if (!MANY && !ANISO && !VARIA)		/* 000 */
                                        // DB("000\n");
  dW=D(11);
  if (G(0,x+1,y  ,z  )) W(p00) += dW;
  if (G(0,x-1,y  ,z  )) W(m00) += dW;
  dW=D(22);
  if (G(0,x  ,y+1,z  )) W(0p0) += dW;
  if (G(0,x  ,y-1,z  )) W(0m0) += dW;
  dW=D(33);
  if (G(0,x  ,y  ,z+1)) W(00p) += dW;
  if (G(0,x  ,y  ,z-1)) W(00m) += dW;
  
#elif (!MANY && !ANISO && VARIA)	/* 001 */
                                        // DB("001\n");
  
  if (conservative) {
    /* Newer: conservative ("mimetic"?) scheme, recap of diff2dv variable scalar diffusion module */
    if (G(0,x+1,y,z)) W(p00) += MEAN( D(11) , DNB(11,x+1,y,z) );
    if (G(0,x-1,y,z)) W(m00) += MEAN( D(11) , DNB(11,x-1,y,z) );
    if (G(0,x,y+1,z)) W(0p0) += MEAN( D(22) , DNB(22,x,y+1,z) );
    if (G(0,x,y-1,z)) W(0m0) += MEAN( D(22) , DNB(22,x,y-1,z) );
    if (G(0,x,y,z+1)) W(00p) += MEAN( D(33) , DNB(33,x,y,z+1) );
    if (G(0,x,y,z-1)) W(00m) += MEAN( D(33) , DNB(33,x,y,z-1) );
  } else {
    /* Old: 2D reduction of RMcF's 3D formula */
    /* Use the tensor */
    dW=D(11);
    if (G(0,x+1,y  ,z  )) W(p00) += dW;
    if (G(0,x-1,y  ,z  )) W(m00) += dW;
    dW=D(22);
    if (G(0,x  ,y+1,z  )) W(0p0) += dW;
    if (G(0,x  ,y-1,z  )) W(0m0) += dW;
    dW=D(33);
    if (G(0,x  ,y  ,z+1)) W(00p) += dW;
    if (G(0,x  ,y  ,z-1)) W(00m) += dW;
    
    if (G(0,x+1,y,z) && G(0,x-1,y,z)) {
      dD11 = DNB(11,x+1,y,z) - DNB(11,x-1,y,z);
    } else {
      dD11 = 0.0;
    }
    if (G(0,x,y+1,z) && G(0,x,y-1,z)) {
      dD22 = DNB(22,x,y+1,z) - DNB(22,x,y-1,z);
    } else {
      dD22 = 0.0;
    }
    if (G(0,x,y,z+1) && G(0,x,y,z-1)) {
      dD33 = DNB(33,x,y,z+1) - DNB(33,x,y,z-1);
    } else {
      dD33 = 0.0;
    }
    /* Use the diff coeff's derivative */
    c1 = dD11;
    c2 = dD22;
    c3 = dD33;
    if (G(0,x+1,y  ,z  )) W(p00)+=quarter*c1;
    if (G(0,x-1,y  ,z  )) W(m00)-=quarter*c1;
    if (G(0,x  ,y+1,z  )) W(0p0)+=quarter*c2;
    if (G(0,x  ,y-1,z  )) W(0m0)-=quarter*c2;
    if (G(0,x  ,y  ,z+1)) W(00p)+=quarter*c3;
    if (G(0,x  ,y  ,z-1)) W(00m)-=quarter*c3;
  }
  
#elif (!MANY && ANISO)	/* 010 + 011*/
  // DB("010+011\n");
  
  if (conservative) {
    /* Newer: conservative ("mimetic"?) scheme, recap of diff2dv variable scalar diffusion module */
    if (G(0,x+1,y,z)) W(p00) += MEAN( D(11) , DNB(11,x+1,y,z) );
    if (G(0,x-1,y,z)) W(m00) += MEAN( D(11) , DNB(11,x-1,y,z) );
    if (G(0,x,y+1,z)) W(0p0) += MEAN( D(22) , DNB(22,x,y+1,z) );
    if (G(0,x,y-1,z)) W(0m0) += MEAN( D(22) , DNB(22,x,y-1,z) );
    if (G(0,x,y,z+1)) W(00p) += MEAN( D(33) , DNB(33,x,y,z+1) );
    if (G(0,x,y,z-1)) W(00m) += MEAN( D(33) , DNB(33,x,y,z-1) );
    
    if (G(0,x+1,y+1,z)) W(pp0) += half*MEAN( D(12) , DNB(12,x+1,y+1,z) );
    if (G(0,x-1,y+1,z)) W(mp0) -= half*MEAN( D(12) , DNB(12,x-1,y+1,z) );
    if (G(0,x+1,y-1,z)) W(pm0) -= half*MEAN( D(12) , DNB(12,x+1,y-1,z) );
    if (G(0,x-1,y-1,z)) W(mm0) += half*MEAN( D(12) , DNB(12,x-1,y-1,z) );
    
    if (G(0,x+1,y,z+1)) W(p0p) += half*MEAN( D(13) , DNB(13,x+1,y,z+1) );
    if (G(0,x-1,y,z+1)) W(m0p) -= half*MEAN( D(13) , DNB(13,x-1,y,z+1) );
    if (G(0,x+1,y,z-1)) W(p0m) -= half*MEAN( D(13) , DNB(13,x+1,y,z-1) );
    if (G(0,x-1,y,z-1)) W(m0m) += half*MEAN( D(13) , DNB(13,x-1,y,z-1) );
    
    if (G(0,x,y+1,z+1)) W(0pp) += half*MEAN( D(23) , DNB(23,x,y+1,z+1) );
    if (G(0,x,y-1,z+1)) W(0mp) -= half*MEAN( D(23) , DNB(23,x,y-1,z+1) );
    if (G(0,x,y+1,z-1)) W(0pm) -= half*MEAN( D(23) , DNB(23,x,y+1,z-1) );
    if (G(0,x,y-1,z-1)) W(0mm) += half*MEAN( D(23) , DNB(23,x,y-1,z-1) );
  } else {
    /* Old: 1D reduction of RMcF's 3D formula */
    /* Use the tensor */
    dW=D(11);
    if (G(0,x+1,y  ,z  )) W(p00) += dW;
    if (G(0,x-1,y  ,z  )) W(m00) += dW;
    dW=D(22);
    if (G(0,x  ,y+1,z  )) W(0p0) += dW;
    if (G(0,x  ,y-1,z  )) W(0m0) += dW;
    dW=D(33);
    if (G(0,x  ,y  ,z+1)) W(00p) += dW;
    if (G(0,x  ,y  ,z-1)) W(00m) += dW;
    
    dW=half*D(12);
    if (G(0,x+1,y+1,z  )) W(pp0) += dW;
    if (G(0,x+1,y-1,z  )) W(pm0) -= dW;
    if (G(0,x-1,y-1,z  )) W(mm0) += dW;
    if (G(0,x-1,y+1,z  )) W(mp0) -= dW;
    dW=half*D(13);
    if (G(0,x+1,y  ,z+1)) W(p0p) += dW;
    if (G(0,x+1,y  ,z-1)) W(p0m) -= dW;
    if (G(0,x-1,y  ,z-1)) W(m0m) += dW;
    if (G(0,x-1,y  ,z+1)) W(m0p) -= dW;
    dW=half*D(23);
    if (G(0,x  ,y+1,z+1)) W(0pp) += dW;
    if (G(0,x  ,y+1,z-1)) W(0pm) -= dW;
    if (G(0,x  ,y-1,z-1)) W(0mm) += dW;
    if (G(0,x  ,y-1,z+1)) W(0mp) -= dW;
    
    /* Use the tensor's derivatives */
    if (G(0,x+1,y,z) && G(0,x-1,y,z)) {
      dD11 = DNB(11,x+1,y,z) - DNB(11,x-1,y,z);
      dD12 = DNB(12,x+1,y,z) - DNB(12,x-1,y,z);
      dD13 = DNB(13,x+1,y,z) - DNB(13,x-1,y,z);
    } else {
      dD11 = 0.0;
      dD12 = 0.0;
      dD13 = 0.0;
    }
    if (G(0,x,y+1,z) && G(0,x,y-1,z)) {
      dD21 = DNB(21,x,y+1,z) - DNB(21,x,y-1,z);
      dD22 = DNB(22,x,y+1,z) - DNB(22,x,y-1,z);
      dD23 = DNB(23,x,y+1,z) - DNB(23,x,y-1,z);
    } else {
      dD21 = 0.0;
      dD22 = 0.0;
      dD23 = 0.0;
    }
    if (G(0,x,y,z+1) && G(0,x,y,z-1)) {
      dD31 = DNB(31,x,y,z+1) - DNB(31,x,y,z-1);
      dD32 = DNB(32,x,y,z+1) - DNB(32,x,y,z-1);
      dD33 = DNB(33,x,y,z+1) - DNB(33,x,y,z-1);
    } else {
      dD31 = 0.0;
      dD32 = 0.0;
      dD33 = 0.0;
    }
    c1 = dD11 + dD21 + dD31;
    c2 = dD12 + dD22 + dD32;
    c3 = dD13 + dD23 + dD33;
    
    if (G(0,x+1,y  ,z  )) W(p00)+=quarter*c1;
    if (G(0,x-1,y  ,z  )) W(m00)-=quarter*c1;
    if (G(0,x  ,y+1,z  )) W(0p0)+=quarter*c2;
    if (G(0,x  ,y-1,z  )) W(0m0)-=quarter*c2;
    if (G(0,x  ,y  ,z+1)) W(00p)+=quarter*c3;
    if (G(0,x  ,y  ,z-1)) W(00m)-=quarter*c3;
  }
  
#elif (MANY && !ANISO && !VARIA)	/* 100 */
                                        // DB("100\n");
  
  /* Correct the weights to make stencil more spherical */
  /* as in EZSCROLL */
  dW=twosixths*D(11); /* D(11)==D(22)==D(33) in this case */
  if (G(0,x+1,y  ,z  )) W(p00) += dW;
  if (G(0,x-1,y  ,z  )) W(m00) += dW;
  if (G(0,x  ,y+1,z  )) W(0p0) += dW;
  if (G(0,x  ,y-1,z  )) W(0m0) += dW;
  if (G(0,x  ,y  ,z+1)) W(00p) += dW;
  if (G(0,x  ,y  ,z-1)) W(00m) += dW;
  dW=sixth*D(11);
  if (G(0,x+1,y+1,z  )) W(pp0) += dW;
  if (G(0,x+1,y-1,z  )) W(pm0) += dW;
  if (G(0,x-1,y-1,z  )) W(mm0) += dW;
  if (G(0,x-1,y+1,z  )) W(mp0) += dW;
  if (G(0,x+1,y  ,z+1)) W(p0p) += dW;
  if (G(0,x+1,y  ,z-1)) W(p0m) += dW;
  if (G(0,x-1,y  ,z-1)) W(m0m) += dW;
  if (G(0,x-1,y  ,z+1)) W(m0p) += dW;
  if (G(0,x  ,y+1,z+1)) W(0pp) += dW;
  if (G(0,x  ,y+1,z-1)) W(0pm) += dW;
  if (G(0,x  ,y-1,z-1)) W(0mm) += dW;
  if (G(0,x  ,y-1,z+1)) W(0mp) += dW;
  
#elif (MANY && !ANISO && VARIA)		/* 101 */
  
  EXPECTED_ERROR("This should not happen. MANY is incompatiable with VARIA");
  
#elif (MANY && ANISO && !VARIA)		/* 110 */
  
  EXPECTED_ERROR("This should not happen. MANY is incompatiable with ANISO");
  
#elif (MANY && ANISO && VARIA)		/* 111 */
  
  EXPECTED_ERROR("This should not happen. MANY is incompatiable with ANISO and VARIA");
  
#endif /* MANY? ANISO? VARIA? */
  /*=============================================*/
  
 } /* for ip */

/*----------------------------------------------*/
/* Third round: central weight for conservation */
/* */						// DB("third round\n");
for (ip=0;ip<np;ip++) {
  P=s->p[ip];
  XP[0]=0;
  for (nb=1;nb<nnb;nb++) XP[0]-=XP[nb];
  /* */ 					/* WGTS3; */
 } /* for ip */
