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

/* INITIALIZATION: READING PARAMETERS AND FORMING THE COMPUTATION LOOP */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "screen.h"
#include "init.h"
#include "qpp.h"
#include "bikt.h"
#include "k_.h"

/* Declare list of all possible devices */
#define D(name) Create create_##name;
#define S(name) Create create_##name;
#include "devlist.h"
#undef D
#undef S

/* Imported variables */
extern FILE     *debug;
extern int      Verbose;
extern INT	Graph;
extern int      idev;
extern int      ndev;
extern Device   dev[];
extern char *device_name;
extern long int inf;

static int create (Name name, Create c, Device *d, char *rest)
{
  char when_name[MAXSTRLEN];
  if (!accepts("name=",&(d->n[0]),name,rest)) return FAILURE;
  if (!acceptkr("when=",0,"always",&(d->c),when_name,rest)) return FAILURE;  
  if (!accept_space(&(d->s),rest)) return FAILURE;
  if (strcmp(d->n,name)) if(!def_dev(d)) return FAILURE;
#if MPI
  d->sync      = 0;
  d->alwaysRun = 0;
#endif
  return(c(d,rest));
}

#define ERR_CREATE EXPECTED_ERROR(\
    "\nError occurred in input file \"%s\" command \"%s\" ending at line %ld pos %ld",\
    inname[depth-1], w, inpline[depth-1], inppos[depth-1])

int init (void)
{
  static char s[MAXSTRLEN];
  char *rest;
  char *w;
  int begun;
  int state_called;
  
  /* this function defines tolerances, and also several functions. */
  if (!init_const()) EXPECTED_ERROR("initializing constants"); 

  idev=0; begun=0; state_called=0;
  for(;;) {
    /* this is what reads in the text line. */
    if NOT(read_command(s,MAXSTRLEN)) EXPECTED_ERROR("reading command");
    if (s[0]=='\0') continue;		/* null command - OK */
    if (s[0]==(char)EOF) break;		/* top level EOF */
    rest=first_word(s,&w," \t\n\r;$");	/* extract the command */
    device_name=w;
    Debug("\n#%d %s:",mpi_rank,w);
    #define CASE(a) else if (0==stricmp(w,a)) 
    if NOT(*w) continue;		/* empty cmd, skip */
    CASE("rem") continue;		/* "remark" comment, skip */
    CASE("if") {			/* conditional command */
      double condition;
      rest=first_word(rest,&w," \t\n\r;$");
      if (!calc(&condition,t_real,w)) ERR_CREATE;
      if (!condition) continue;
      rest=first_word(rest,&w," \t\n\r;$"); /* parse the rest of the command */
    }
    CASE("def") {			/* define k-variable or string macro  */
      if(!def(rest)) ERR_CREATE; 
      begun=1;
      continue;
    }
    CASE("state") {			/* create the computational grid */
      if (begun) MESSAGE("\n");
      MESSAGE("state ");
      if(!state(rest)) ERR_CREATE; 
      MESSAGE("$");
      begun=0;
      state_called = 1;
      continue;
    }
    CASE("screen") {			/* create the BGI window */
#if MPI
      MESSAGE("\n/* The screen command is disabled when using MPI. Your simulation will continue without it. */");
#else
      if (Graph) {
	MESSAGE("\nscreen ");
	// this initialises the screen.
	if (makescreen(rest)) {
	  MESSAGE("$");
	} else {
	  MESSAGE("\n/* Will proceed without onscreen graphics */ $");
	}
      } else {
	MESSAGE("\n/* With nograph option, the 'screen' command is ignored */");
      }
      begun=0;
#endif
      continue;
    }    
    CASE("end") break;			/* end of input stream */

    if (idev>=MAXDEV) {
      MESSAGE("\n The number of devices in the input file exceeded the maximum of %d allowed in this compile.\n",MAXDEV);
      ERR_CREATE;
    }
    
#if MPI
    /*  Don't create devices on idle processes. */
    #define D(name)                                                     \
    CASE(#name) {                                                       \
      if (!I_AM_IDLE) {							\
        if (!state_called) {						\
          MESSAGE("\nDevices cannot be created before state is called."); \
          ERR_CREATE;                                              	\
        }                                                               \
        MESSAGE("\n%s ",#name);						\
        if(!create(#name,create_##name,dev+idev,rest)) ERR_CREATE;	\
        MESSAGE("$");                                                	\
        begun=0;                                                        \
      }                                                                 \
    }
    
    /*  Disable Sequential-only devices. */
    #define S(name)                                                     \
    CASE(#name) {                                                       \
      if (!I_AM_IDLE) {							\
        if (!state_called) {						\
          MESSAGE("\nDevices cannot be created before state is called."); \
          ERR_CREATE;							\
        }                                                               \
        MESSAGE("\n/* The %s device is disabled when using MPI. Your simulation will continue without it. */",#name); \
        begun=0;                                                        \
        continue;                                                       \
      }                                                                 \
    }
#else
    #define D(name)                                                     \
    CASE(#name) {                                                       \
      if (!state_called) {						\
        MESSAGE("\nDevices cannot be created before state is called."); \
	ERR_CREATE;							\
      }                                                                 \
      MESSAGE("\n%s ",#name);						\
      if (!create(#name,create_##name,dev+idev,rest)) ERR_CREATE;	\
      MESSAGE("$");                                                  	\
      begun=0;                                                          \
    }
      
    #define S(name)                                                     \
    CASE(#name) {                                                       \
      if(!state_called) {                                               \
        MESSAGE("\nDevices cannot be created before state is called."); \
	ERR_CREATE;							\
      }                                                                 \
      MESSAGE("\n%s ",#name);                                           \
      if (!create(#name,create_##name,dev+idev,rest)) ERR_CREATE;	\
      MESSAGE("$");                                                  	\
      begun=0;                                                          \
    }
#endif

    #include "devlist.h"
    #undef D
    #undef S

    else { MESSAGE("\nUnknown device name '%s'\n",w); ERR_CREATE; }
    idev++;
  } /*  for(;;) */
  Debug("\n#%d end of input file\n",mpi_rank);
    
  MESSAGE("\nend of input file $");
  MESSAGE("\nLoop of %d devices created:",ndev=idev);
  for (idev=0;idev<ndev;idev++) {
    if (idev%5==0) MESSAGE("\x01""\n  ");
        MESSAGE(" (%d)%s",idev,dev[idev].n);
  }
  
  return(SUCCESS);
} /* init */

/* free all dynamical objects */
void term (void)
{
  Device d;
  if (ndev) for (idev=ndev-1;idev>=0;idev--) {
    d=dev[idev];
#if MPI
    d.d(d.s,d.par,d.sync,d.alwaysRun);
#else
    d.d(d.s,d.par);			/* terminate each device in turn */
#endif
  }
  state_free ();			/* free the state array(s) */
  term_const ();			/* clear the symbol table and free variables */
} /* term */

