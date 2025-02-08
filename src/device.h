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

#ifndef _DEVICE_H_
#define _DEVICE_H_
#ifdef _OWN
#undef _OWN
#define OWN
#endif

#if MPI
#include <mpi.h>
#endif

#include "beatbox.h"
#include "extern.h"
#include "error.h"
#include "k_.h"

extern FILE *res;			/* log file */
extern FILE *debug;			/* debug file */
#define MAXDEV 1024			/* Maximum number of Devices that can be used */
#define MAXSTRLEN 8192			/* length of various buffers, including that for reading input file */
EXTERN char buf[MAXSTRLEN];		/* Common general-purpose buffer */
typedef REAL *Condition;		/* Condition when the Device will operate */

/* typedef struct { ... } devicePoint;  - defined in state.h */

/*! Description of the set of points on which the Device runs */
typedef struct {
  ssize_t x0, y0, z0,			/* set the brutto limits on which the Device will */
          x1, y1, z1;			/* operate in this thread */
  ssize_t global_x0, global_x1,		/* ..., in the overall grid; */
          global_y0, global_y1,		/* local == global in seq mode */
          global_z0, global_z1;
  int restricted;			/* flag to restrict device space to tissue points */
  int listpoints;			/* flag to catalogue all Device points for speed */
  REAL *relist;				/* address of k-var flagging that relisting is due */
  char relistname[maxname];		/* .., name of */
  int relisted;				/* flag that relisting was done by Xs, if any, not renewed yet */
  char where[MAXSTRLEN];		/* the source k-code defining grid points the device runs */
  pp_fn wherecode;			/* .., pointer to compiled k-code */
  int hasu, hasg;			/* whether the 'where' expression uses u[] and geom[] variables */
  INT x,y,z;	      			/* local k-vars accessible to the 'where' expression */
  real *u;				/* [vmax], .. */
  real *geom;				/* [geom_vmax], .. */
  size_t np;				/* number of listed grid points in the device space */
  devicePoint *p;			/* p[np] list of those points */
  int runHere;				/* Is ==1 in sequential mode, otherwise what is says on the tin */
  int nowhere;				/* All the rest is irrelevant if this one is nonzero? */
} Space;


/*! Run over the all the points of the Device */
#define DO_FOR_ALL_POINTS				\
  real *u;						\
  if (s.listpoints) {					\
    long ip;						\
    for (ip=0;ip<s.np;ip++) {				\
      devicePoint P=s.p[ip];				\
      u=P.u;						\
      INT *x=&(P.x);					\
      INT *y=&(P.y);					\
      INT *z=&(P.z);					\
      COMMANDS;						\
    }							\
  } else if (s.restricted) { /* not listpoints */	\
    SPACE_VAR(INT,x);					\
    SPACE_VAR(INT,y);					\
    SPACE_VAR(INT,z);					\
    for ((*z)=s.z0;(*z)<=s.z1;(*z)++) {			\
      for ((*y)=s.y0;(*y)<=s.y1;(*y)++) {		\
	for ((*x)=s.x0;(*x)<=s.x1;(*x)++) {		\
	  if (isTissue(*x,*y,*z)) {			\
	    u = New + ind(*x,*y,*z,0);			\
	    COMMANDS;					\
	  } /*  if isTissue */				\
	} /*  for *x */					\
      } /*  for *y */					\
    } /*  for *z */					\
  } else { /* not listpoints nor restricted */		\
    SPACE_VAR(INT,x);					\
    SPACE_VAR(INT,y);					\
    SPACE_VAR(INT,z);					\
    for ((*z)=s.z0;(*z)<=s.z1;(*z)++) {			\
      for ((*y)=s.y0;(*y)<=s.y1;(*y)++) {		\
	for ((*x)=s.x0;(*x)<=s.x1;(*x)++) {		\
	  u = New + ind(*x,*y,*z,0);			\
	  COMMANDS;					\
	} /*  for *x */					\
      } /*  for *y */					\
    } /*  for *z */					\
  } /* if listpoints or restricted */

/* Graphics display information for some Devices */
typedef struct {int row0,col0,row1,col1,color,area;} BGIWindow;

/* Internal persistent storage used by the Device */
typedef void *Par;

#if MPI
/**********************/
/* PARALLEL */
#define PROC(name) int name(Space s,Par par,int sync,int alwaysRun)
#else
#define PROC(name) int name(Space s,Par par)
#endif
typedef PROC(Proc);

/* Device name */
typedef char Name[32];

/* Device structure */
typedef struct {	
    Condition c; /* When the device is to be run */
    Space s;     /* Description of the 3D grid subset on which the device operates */
    Par par;     /* Store for persistent data within the device */
    Proc *p;     /* Pointer to the device's run function */
    Proc *d;     /* Pointer to the device's destroy function */
    Name n;      /* Device name */
#if MPI
    int sync;
    int alwaysRun;
#endif
} Device;

/* Function that creates the device given its description in the input file */
typedef int Create (Device *dev, char *w);

/*==========================================*/
/* Shortcuts for std pieces of device codes */

/*! Declaration of a Device's constant parameter: value defined at the time of device creation */
#define DEVICE_CONST(type,name) type name=S->name;

/*! Declaration of a Device's variable: may change during run, preserved between runs */
#define DEVICE_VAR(type,name)   type *name=&(S->name);

/*! Declaration of a Device's array of variables */
#define DEVICE_ARRAY(type,name) type *name=&(S->name[0]);

/*! Declaration of a parameter that may be recalculated before each run */
/* If pointer is NULL then use code, assign any result to name in S and to local */
/* #define DEVICE_PAR(type,name)   type name=(S->name##ptr)?(S->name=*(type *)(S->name##ptr)):(memcpy(&(S->name),execute(S->name##exe),sizeof(type)),S->name); */

/*! Declaration of a parameter that may be recalculated at any point. */
/* Exec k-code if given, in any case assign to the local the value from the pointer */
#define DEVICE_PAR(type,name) type name=(S->name=*(type *)((S->name##exe)?((type *)execute(S->name##exe)):(S->name##ptr)));

/*! Entry for recalculable parameter in the Device's structure */
/* If pointer is NULL then use code, and name in S to contain any result */
#define PAR(type,name) type name; type *name##ptr; char *name##src; pp_fn name##exe;

/*! Get value of a recalculable parameter at any point and moment of device work */
/* If pointer is NULL then use code, assign any result to name  and return */
#define VAL(name) ((S->name##ptr)?(name=S->name=*(S->name##ptr)):(memcpy(&(S->name),execute(S->name##exe),sizeof(name)),name=S->name))

/*! Set value of a recalculable parameter at any point and moment of device work */
/* */
#define SETVAL(name,value) {S->name##ptr=&(S->name); S->name=value;} /* leave code untouched */
/* #define SETEXPR(name,expo) {S->name##ptr=NULL; *(S->code)=compile...} - TODO if and when needed */

/*! Copy all three hypostases of a recalculable parameter within S */
#define CPYPAR(dst,src) { \
  S->dst=S->src; \
  S->dst##ptr=S->src##ptr; \
  FREE(S->name##code); \
  S->dst##ptr=strap(S->src##ptr); \
}



/*! Declaration of a Space's constant parameter: value defined at the time of device creation */
#define SPACE_CONST(type,name) type name=s.name;

/*! Declaration of a Space's variable: may change during run, preserved between runs */
#define SPACE_VAR(type,name)   type *name=&(s.name);

/*! Declaration of a Space's array of variables */
#define SPACE_ARRAY(type,name) type *name=&(s.name[0]);


#if MPI
#define RUN_HEAD(name)		 \
PROC (run_##name) {		 \
   if (sync) {haloSwap();}	 \
   if (s.runHere || alwaysRun) { \
      STR *S = (STR *) par;

#define RUN_TAIL(name)		 \
  } /* if s.runHere */		 \
  return 1;			 \
}

#define DEVICE_REQUIRES_SYNC dev->sync = 1;

#define DEVICE_ALWAYS_RUNS dev->alwaysRun = 1;

/* Sync is set to 0 to prevent haloSwaps in delegated devices
 * This isn't just for efficiency, but also avoids deadlocks if
 * any device's runHere is false.
 */
#define DELEGATE_TO_DEVICE(name) run_##name(s,w,par,0,alwaysRun);

#else /* not MPI */
/**********************/
/* SEQUENTIAL */

#define RUN_HEAD(name)	\
PROC (run_##name) {	\
  STR *S = (STR *) par;

#define RUN_TAIL(name)	\
  return 1;		\
}

#define DEVICE_REQUIRES_SYNC
#define DEVICE_ALWAYS_RUNS

#define DELEGATE_TO_DEVICE(name) run_##name(s,w,par);

#endif /* not MPI */
/**********************/

#define DESTROY_HEAD(name)	\
PROC (destroy_##name) {		\
  STR *S = (STR *) par;

#define DESTROY_TAIL(name)	\
  FREE(par);			\
  return 1;			\
}

#define SAFE_CLOSE(file)        \
  if ((file) != NULL && 	\
      (file) != res && 		\
      (file) != debug && 	\
      (file) != stdout && 	\
      (file) != stderr && 	\
      ftell(file) !=-1 		\
      ) {			\
    fseek((file),0,SEEK_END);	\
    fflush(file);		\
    fclose(file);		\
  }				\
  (file)=NULL;

#define CREATE_HEAD(name)			\
int create_##name (Device *dev, char *w)  {	\
  STR *S = (STR *) Calloc(1,sizeof(STR));	\
  if (!S) ABORT("cannot create %s",#name);

#define CREATE_TAIL(name,areaval)	\
  dev->p = (Proc *) run_##name;		\
  dev->d = (Proc *) destroy_##name;	\
  dev->par  = S;			\
  return(1);				\
}

#define DEVICE_IS_SPACELESS	    \
if (spaceParametersExist(w)) {	    \
  MESSAGE("/* WARNING: %s does not use space parameters (x0,x1,y0,y1,z0,z1,where,list,when_list).\nThe parameters provided will be ignored. */",dev->n); \
}

#define DEVICE_HAS_DEFAULT_SPACE	  \
if (spaceParametersExist(w)) {		  \
  EXPECTED_ERROR("%s insists on using the default space.\nPlease try again without space parameters (x0,x1,y0,y1,z0,z1,where,list,when_list).\n",dev->n); \
}

#define DEVICE_IS_RECTANGULAR	  \
if (!spaceIsBox(w)) {		  \
  EXPECTED_ERROR("%s insists on using a box space.\nPlease try again without point list parameters (where,list,when_list).\n",dev->n); \
}


#define DEVICE_MUST_BE_NOWHERE		\
if (dev->s.nowhere != 1) {		\
  MESSAGE("/* WARNING: %s requires that nowhere is equal to 1.\nAdd 'nowhere=1' to remove this warning.*/\n",dev->n); \
  dev->s.nowhere=1; \
}

#define DEVICE_OPERATES_ON_A_SINGLE_POINT					\
  if ((!spaceIsBox(w))||(find_key("x1=",w))||(find_key("y1=",w))||(find_key("z1=",w))) {	\
    EXPECTED_ERROR("/* %s does not use the x1,y1,z1,where,list,when_list space parameters"	\
		   " -- it only operates on a single point.\n"			\
		   "Please try again without these parameters.*/\n",dev->n);	\
  }

/*  Some handy utilities for devices running on MPI. */
#if MPI

/*  Error handling variables. */
EXTERN MPI_Status 	status;
EXTERN int 		mpi_errno;
EXTERN char 		error_string[MPI_MAX_ERROR_STRING];
EXTERN int 		error_length;

#define CHECK_MPI_SUCCESS(error_message)                      \
  if(mpi_errno != MPI_SUCCESS) {                              \
    MPI_Error_string(mpi_errno, error_string, &error_length); \
    ABORT("Process %d: %s\n\tMPI Message: %s\n", mpi_rank, error_message, error_string); \
  }

/* this defines its own buf in case there is another one outside */
#define MPIDO(cmd,...)							\
  mpi_errno = cmd;							\
  if (mpi_errno != MPI_SUCCESS) {					\
    MPI_Error_string(mpi_errno, error_string, &error_length);		\
    char buf[MAXSTRLEN];						\
    snprintf(buf,MAXSTRLEN,__VA_ARGS__);				\
    ABORT("%s\nMPI Message: %s\n",buf,error_string);			\
  }

/* 	Define new group and communicator for I/O, 
 *	including only instances with runHere = true.	*/
int deviceCommunicator (int runHere, MPI_Comm *new_comm);
int deviceCommunicatorWithFirstRank (int runHere, MPI_Comm *new_comm, int *first);

/*  	Exchange internal boundaries. */
int haloSwap(void);

#endif /*  MPI */
#undef EXTERN
#endif /* end of include guard: _DEVICE_H_ */
