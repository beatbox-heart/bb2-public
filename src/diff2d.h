/***************************************************************/
/* 2D diffusion                                                */
/***************************************************************/

real Gam=1.0/(hx*hx);			//DB("Gam=%lg\n",(double)Gam); /* prescale */
real f1, f2;				/* fibre cosines */
real Dlh, Dth, Ddh, Dsh, dW;		/* convenience intermediates */
real dD11, dD12, dD21, dD22;		/* diff tensor gradient */
real c1, c2;				/* contraction of ... : convection vector */
/* */					//DB("loctb=%p\n",S->loctb);

DB("2D version: GEOM=%d ANISO=%d MANY=%d VARIA=%d\n", GEOM, ANISO, MANY, VARIA);

#if !VARIA
  k_on();
#if ANISO
  DEVICE_PAR(real,Dl);			/* longitudinal */	//DB("Dl=%g\n",Dl);
  DEVICE_PAR(real,Dt);			/* transversal */	//DB("Dt=%g\n",Dt);
  real D=(S->D)=			/* scalar for exceptional points */
    (S->Dptr)?(*(real *)(S->Dptr)):
    (S->Dexe)?(*(real *)execute(S->Dexe)):
    ( sqrt(Dl*Dt) );			//DB("D=%g\n",D);
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
  assert(P.z==0);
  CLEAR_WEIGHTS;
#if !MPI
  IPNBR(x,y,0)=ip; 			/* count this point in the ipgrid */
#endif
  /* */					//DB("ip=%ld xy=%ld,%ld Gam=%g\n",ip,x,y,Gam);
#if VARIA
  if (S->loctb) {			/* yes local dependencies */
    k_on();
    S->x=x;				/* copy point coords into */
    S->y=y;				/*         local symtable */
    if (S->u) memcpy(S->u,P.u,vmax*sizeof(real)); /* ditto, vector of u, if relevant */
#if GEOM
    if (S->geom) memcpy(S->geom,&(G(0,x,y)),geom_vmax*sizeof(real)); /* ditto, vector of geom */
#endif
#if ANISO
    DEVICE_PAR(real,Dl);		/* longitudinal */	//DB("Dl=%g\n",Dl);
    DEVICE_PAR(real,Dt);		/* transversal */	//DB("Dt=%g\n",Dt);
    real D=S->D=			/* scalar for exceptional points */
      (S->Dptr)?(*(real *)(S->Dptr)):
      (S->Dexe)?(*(real *)execute(S->Dexe)):
      ( sqrt(Dl*Dt) );			/* default if not given in input file */ //DB("D=%g\n",D);
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
  f1=G(1,x,y);				/* fiber direction */
  f2=G(2,x,y);				/*         cosines */
  /* */					//DB("f=%g,%g\n",f1,f2);
  if (isUnitVector(f1,f2,0,UNIT_VECTOR_TOLERANCE)) {
    D(11) = Ddh*f1*f1 + Dth;
    D(22) = Ddh*f2*f2 + Dth;
    D(12) = Ddh*f1*f2;
  } else { /* not isUnitVector */
    D(11) = Dsh;
    D(22) = Dsh;
    D(12) = 0.0;
  } 					//DB("D[][]=\n%g\t%g\t%g\n\n",D(11),D(12),D(22));
#else /* not ANISO */
  D(11) = Dsh;				/* ==D(22) */
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
/* */					DB("second round\n");
for (ip=0;ip<np;ip++) {
  P=s->p[ip];
  x=P.x;
  y=P.y;
  
  /* Use the tensor */
#if MANY
  /* Correct the weights to make stencil more spherical */
  /* This is  Patra-Karttunen variant, as in EZSPIRAL */
  dW=foursixths*D(11); /* D(11)==D(22) in this case */ 
  if (G(0,x+1,y  )) W(p0) += dW;
  if (G(0,x-1,y  )) W(m0) += dW;
  if (G(0,x  ,y+1)) W(0p) += dW;
  if (G(0,x  ,y-1)) W(0m) += dW;
  dW=sixth*D(11);
  if (G(0,x+1,y+1)) W(pp) += dW;
  if (G(0,x+1,y-1)) W(pm) += dW;
  if (G(0,x-1,y-1)) W(mm) += dW;
  if (G(0,x-1,y+1)) W(mp) += dW;
#else /* not MANY */
  dW=D(11);
  if (G(0,x+1,y  )) W(p0) += dW;
  if (G(0,x-1,y  )) W(m0) += dW;
  dW=D(22);
  if (G(0,x  ,y+1)) W(0p) += dW;
  if (G(0,x  ,y-1)) W(0m) += dW;
#endif  /* not MANY */

#if ANISO
  dW=half*D(12);
  if (G(0,x+1,y+1)) W(pp) += dW;
  if (G(0,x+1,y-1)) W(pm) -= dW;
  if (G(0,x-1,y-1)) W(mm) += dW;
  if (G(0,x-1,y+1)) W(mp) -= dW;
#endif /* ANISO */
  
#if (VARIA || ANISO)
  /* Use the tensor's derivatives */
  /* TODO: make them right into cj, w/o intermediate dDij ? */
  if (G(0,x+1,y) && G(0,x-1,y)) {
    dD11 = DNB(11,x+1,y) - DNB(11,x-1,y);
#if ANISO
    dD12 = DNB(12,x+1,y) - DNB(12,x-1,y);
#endif
  } else {
    dD11 = 0.0;
#if ANISO
    dD12 = 0.0;
#endif
  }
  if (G(0,x,y+1) && G(0,x,y-1)) {
#if ANISO
    dD21 = DNB(21,x,y+1) - DNB(21,x,y-1);
#endif
    dD22 = DNB(22,x,y+1) - DNB(22,x,y-1);
  } else {
#if ANISO
    dD21 = 0.0;
#endif
    dD22 = 0.0;
  }
#if ANISO
  c1 = dD11 + dD21;
  c2 = dD12 + dD22;
#else
  c1 = dD11;
  c2 = dD22;
#endif
  
  if (G(0,x+1,y  )) W(p0)+=quarter*c1;
  if (G(0,x-1,y  )) W(m0)-=quarter*c1;
  if (G(0,x  ,y+1)) W(0p)+=quarter*c2;
  if (G(0,x  ,y-1)) W(0m)-=quarter*c2;
#endif /* VARIA || ANISO */
 } /* for ip */

/*----------------------------------------------*/
/* Third round: central weight for conservation */
/* */						DB("third round\n");
for (ip=0;ip<np;ip++) {
  P=s->p[ip];
  XP[0]=0;
  for (nb=1;nb<nnb;nb++) XP[0]-=XP[nb];
  /* */ 					WGTS2;
 } /* for ip */
