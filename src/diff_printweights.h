/* Macros to print out the weights of the current point, */
/* for debug purposes only. Remove when not needed.      */
#define WGTS1 							\
  if (debug) {							\
    fprintf(debug, 						\
	    "#%d %s:%d t=%ld idev=%d "				\
	    "%ld:(%ld,%ld) W=\n"				\
	    "%g\t%g\t%g\n",					\
	    mpi_rank,__FILE__,__LINE__,t,idev,			\
	    ip,P.x,P.y,						\
	    XP[i1m],XP[i10],XP[i1p]);				\
    fflush(debug);						\
  }
#define WGTS2							\
  if (debug) {							\
    switch (nnb) {						\
    case i2nnb:							\
      fprintf(debug,						\
	      "#%d %s:%d t=%ld idev=%d "			\
	      "%ld:(%ld,%ld) W=\n"				\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n",					\
	      mpi_rank,__FILE__,__LINE__,t,idev,		\
	      ip,P.x,P.y,					\
	      0.0,     XP[i20p],0.0,				\
	      XP[i2m0],XP[i200],XP[i2p0],			\
	      0.0,     XP[i20m],0.0);				\
      break;							\
    case m2nnb:							\
      fprintf(debug,						\
	      "#%d %s:%d t=%ld idev=%d "			\
	      "%ld:(%ld,%ld) W=\n"				\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n",					\
	      mpi_rank,__FILE__,__LINE__,t,idev,		\
	      ip,P.x,P.y,					\
	      XP[m2mp],XP[m20p],XP[m2pp],			\
	      XP[m2m0],XP[m200],XP[m2p0],			\
	      XP[m2mm],XP[m20m],XP[m2pm]);			\
      break;							\
    default: EXPECTED_ERROR("unknown stencil size %d\n",nnb);	\
    } /* switch nnb */						\
    fflush(debug);						\
  }
#define WGTS3							\
  if (debug) {							\
    switch (nnb) {						\
    case i3nnb:							\
      fprintf(debug,						\
	      "#%d %s:%d t=%ld idev=%d "			\
	      "%ld(%ld,%ld,%ld) W=\n"				\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "-----------\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "-----------\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "===========\n\n",				\
	      mpi_rank,__FILE__,__LINE__,t,idev,		\
	      ip,P.x,P.y,P.z,					\
	      0.0      ,0.0      ,0.0,				\
	      0.0      ,XP[i300p],0.0,				\
	      0.0      ,0.0      ,0.0,				\
	      0.0      ,XP[i30p0],0.0,				\
	      XP[i3m00],XP[i3000],XP[i3p00],			\
	      0.0      ,XP[i30m0],0.0,				\
	      0.0      ,0.0      ,0.0,				\
	      0.0      ,XP[i300m],0.0,				\
	      0.0      ,0.0	,0.0				\
	      );						\
      break;							\
    case m3nnb:							\
      fprintf(debug,						\
	      "#%d %s:%d t=%ld idev=%d "			\
	      "%ld(%ld,%ld,%ld) W=\n"				\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "-----------\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "-----------\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "%g\t%g\t%g\n"					\
	      "===========\n\n",				\
	      mpi_rank,__FILE__,__LINE__,t,idev,		\
	      ip,P.x,P.y,P.z,					\
	      0.0      ,XP[m30pp],0.0,				\
	      XP[m3m0p],XP[m300p],XP[m3p0p],			\
	      0.0      ,XP[m30mp],0.0,				\
	      XP[m3mp0],XP[m30p0],XP[m3pp0],			\
	      XP[m3m00],XP[m3000],XP[m3p00],			\
	      XP[m3mm0],XP[m30m0],XP[m3pm0],			\
	      0.0      ,XP[m30pm],0.0,				\
	      XP[m3m0m],XP[m300m],XP[m3p0m],			\
	      0.0      ,XP[m30mm],0.0				\
	      );						\
      break;							\
    default: EXPECTED_ERROR("unknown stencil size %d\n",nnb);	\
    }; /* switch nnb */						\
    fflush(debug);						\
  }
