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

/*
26/9/2012
SK
First attempt at making the ecg device.
Firstly, get it working for constant D.
Formula for ECG taken from:
Simulation of cardiac activity and the ECG using a
heart model with a reaction-diffusion
action potential
Omer Berenfeld, and Shimon Abboudt
Med. Eng. Phys. Vol. 18. No. 8. pp 615-625, 1996

Test would be:
1. Construct heterogenous TTNP strand and compute ECG (serial and parallel).
Compare that to the ECG in CINC 2012 paper.
need to activate the when option.
*/
#include <math.h>
#include <assert.h>
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
  real D;		/* Isotropic diffusion coefficient. */
  real hx;		/* Space step */
  int xl, yl, zl;	/* lead position */
  int v0;		/* layer for the voltage field */
#if MPI
  MPI_Comm comm;	/* Communicator for reduction operation. */
  int root;		/* Rank of device root in ALL_ACTIVE_PROCS */
  MPI_Op global_reduce;	/* Global reduction operation to use. */
#endif
} STR;

/****************/
RUN_HEAD(ecg)
{
  int x, y, z;
  real *u;
  DEVICE_CONST(real,hx);
  DEVICE_CONST(real,D);
  DEVICE_CONST(int, xl);
  DEVICE_CONST(int, yl);
  DEVICE_CONST(int, zl);
  DEVICE_CONST(int, v0);
  //#if MPI // I am using usual MPI functions for these simple calculations. No need to invoke DEVICE_CONST macros
  //	DEVICE_CONST(MPI_Comm, comm)
  //	DEVICE_CONST(int, root) // this definition of root is never used. could rid of it at test time.
  //	DEVICE_CONST(MPI_Op, global_reduce)
  //#endif
  real gam=D/(2.0*hx);
  real ecg = 0.0;
  real total_ecg; // this one is on the root processor only.
  real r, r3;
  real rx, ry, rz;
  for(x=s.x0;x<=s.x1;x++) {
    for(y=s.y0;y<=s.y1;y++){
      for(z=s.z0;z<=s.z1;z++){
	if(isTissue(x,y,z)){
	  rx = (x-xl);
	  ry = (y-yl);
	  rz = (z-zl);
	  r = sqrt(rx*rx + ry*ry + rz*rz);
	  if(fabs(r)<0.0001) r = 0.0001;
	  r3 = r*r*r;
	  u=New+ind(x,y,z,v0);
	  if(dim>=1){
	    if(isTissue(x+1,y,z)&&isTissue(x-1,y,z))
	      ecg+=(u[+DX] - u[-DX])*rx/r3;
	  }                              
	  if(dim>=2){
	    if(isTissue(x,y+1,z)&&isTissue(x,y-1,z))
	      ecg+=(u[+DY] - u[-DY])*ry/r3;
	  }
	  if(dim>=3){
	    if(isTissue(x,y,z+1)&&isTissue(x,y,z-1))
	      ecg+=(u[+DZ] - u[-DZ])*rz/r3;
	  }
	} /* if isTissue() */
      } /* for z */
    } /* for y */
  } /* for x */
  
  ecg = gam*ecg;
  
  // now you need to reduce it, and put it into a file.
#if MPI
  /* Using default communicator for the moment since ecg is for the whole model. SK notes: http://stackoverflow.com/questions/7915891/mpi-communicator-error */
  MPI_Reduce(&ecg, &total_ecg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  if(mpi_rank==0){
    FILE *ecgfile;
    ecgfile = fopen("ecg.dat","a+");
    fprintf(ecgfile,"%ld\t%g\n",t,total_ecg);
    fclose(ecgfile);
  }
#else
  FILE *ecgfile;
  ecgfile = fopen("ecg.dat","a+");
  fprintf(ecgfile,"%ld\t%g\n",t,ecg);
  fclose(ecgfile);
#endif
}
RUN_TAIL(ecg)

DESTROY_HEAD(ecg)
DESTROY_TAIL(ecg)

CREATE_HEAD(ecg)
{
  DEVICE_REQUIRES_SYNC;
  DEVICE_HAS_DEFAULT_SPACE;

  ACCEPTR(D,RNONE,0.,RNONE);
  ACCEPTR(hx,RNONE,0.,RNONE);
  ACCEPTI(xl,INONE,0,+10000);
  ACCEPTI(yl,INONE,0,+10000);
  ACCEPTI(zl,INONE,0,+10000);
  ACCEPTI(v0,INONE,0,(int)vmax-1);
}
CREATE_TAIL(ecg,1)

