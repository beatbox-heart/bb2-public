/***************************************************************/
/* 1D diffusion                                                */
/***************************************************************/

real Gam=1.0/(hx*hx);			//DB("Gam=%lg\n",(double)Gam); /* prescale */
real Dsh;				/* convenience intermediate */
real dD11;				/* diff tensor gradient */
real c1;				/* .., contraction: convection vector */
/**/					//DB("loctb=%p\n",S->loctb);

DB("1D version: GEOM=%d ANISO=%d MANY=%d VARIA=%d\n", GEOM, ANISO, MANY, VARIA);

#if ANISO
EXPECTED_ERROR("this should not happen");
#endif /* ANISO */

#if MANY
EXPECTED_ERROR("this should not happen");
#endif /* MANY */

#if !VARIA
if ((S->loctb) == NULL) {		/* variable diff but no local dependencies */
  k_on();
  DEVICE_PAR(real,D);			/* scalar diffusivity only */
  Dsh=D*Gam;				/* scaled scalar */
  k_off();				DB("D=%g Gam=%g Dsh=%g\n",D,Gam,Dsh);
 } /* if !loctb */
#endif /* VARIA */

/*---------------------------------------------------*/
/* First round: the diff tensor and ipgrid signature */
/* */					DB("first round\n");
for (ip=0;ip<np;ip++) { 
  P=s->p[ip];
  x=P.x;
  assert(P.y==0);
  assert(P.z==0);
  CLEAR_WEIGHTS;
#if !MPI
  IPNBR(x,0,0)=ip; 			/* count this point in the ipgrid */
#endif
  /* */					//DB("ip=%ld xy=%ld,%ld Gam=%g\n",ip,x,y,Gam);
#if VARIA
  if (S->loctb) {			/* yes local dependencies */
    k_on();
    S->x=x;				/* copy point coord */
    if (S->u) memcpy(S->u,P.u,vmax*sizeof(real)); /* ditto, vector of u, if relevant */
#if GEOM
    if (S->geom) memcpy(S->geom,&(G(0,x)),geom_vmax*sizeof(real)); /* ditto, vector of geom */
#endif
    DEVICE_PAR(real,D);			/* scalar diffusivity only */
    Dsh=D*Gam;				/* scaled scalar */
    k_off();				
  } /* if loctb */
#endif /* VARIA */
  D(11) = Dsh;				DB("x=%ld D11=%g\n",x,D(11));
 } /* for ip */

#if (VARIA || ANISO) && MPI
/* Variable diffusivity needs neighbours' diff coeff */
ASSERT(d0>=0);			/* here diff coeff has to be on the grid */
haloSwap();			/* make it available for neighbouring subdomains */
#endif /* (VARIA || ANISO) && MPI */

/*---------------------------------*/
/* Second round: periphery weights */
/* */					DB("second round\n");
for (ip=0;ip<np;ip++) {
  P=s->p[ip];
  x=P.x;				DB("x=%ld D11=%g\n",x,D(11));
  
  /* Use the tensor */
  if (G(0,x+1)) W(p) += D(11);
  if (G(0,x-1)) W(m) += D(11);
  
#if VARIA /* no ANISO in 1D */
  /* Use the diff coeff's derivative */
  /* TODO: make them right into c1, w/o intermediate dD11 ? */
  if (G(0,x+1) && G(0,x-1)) {
    dD11 = DNB(11,x+1) - DNB(11,x-1);
  } else {
    dD11 = 0.0;
  }
  c1 = dD11;
  
  if (G(0,x+1)) W(p)+=0.25*c1;
  if (G(0,x-1)) W(m)-=0.25*c1;
#endif /* VARIA */
 } /* for ip */

/*----------------------------------------------*/
/* Third round: central weight for conservation */
/* */						DB("third round\n");
for (ip=0;ip<np;ip++) {
  P=s->p[ip];
  XP[0]=0;
  for (nb=1;nb<nnb;nb++) XP[0]-=XP[nb];
  /* */ 					WGTS1;
  /* if (debug) { */
  /*   printf("debug=%p\n",debug); */
  /*   fprintf(debug,"mpi_rank=%d\n",mpi_rank); */
  /*   fprintf(debug,"t=%ld\n",t); */
  /*   fprintf(debug,"idev=%d\n",idev); */
  /*   fprintf(debug,"ip=%ld\n",ip); */
  /*   fprintf(debug,"&P=%p\n",&P); */
  /*   fprintf(debug,"P.x=%d\n",(int)P.x); */
  /*   fprintf(debug,"P.y=%d\n",(int)P.y); */
  /*   fprintf(debug,"XP=%p\n",XP); */
  /*   fprintf(debug,"i1m=%d\n",i1m); */
  /*   fprintf(debug,"i10=%d\n",i10); */
  /*   fprintf(debug,"i1p=%d\n",i1p); */
  /*   fprintf(debug, 						 */
  /* 	    "#%d %s:%d t=%ld idev=%d "				 */
  /* 	    "%ld:(%ld,%ld) W=n"				 */
  /* 	    "%gt%gt%gn",					 */
  /* 	    mpi_rank,__FILE__,__LINE__,t,idev,			 */
  /* 	    ip,P.x,P.y,						 */
  /* 	    XP[i1m],XP[i10],XP[i1p]);				 */
  /*   fflush(debug); */
  /* } */
 } /* for ip */
