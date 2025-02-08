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

/* Output a 4D box to a file, each mode value reduced to unsigned char	*/
/* File format: pure information, without markers 			*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "sequence.h"
#include "byte.h"
#include "mpi_io_choice.h"

typedef struct {
  int append;
  sequence file;
  int v0, v1;		/* Range of layers */
  real u0,u1;		/* Range of values to transform node to a byte */
  int bgbyte;		/* Byte representing the void */
  long local_size;	/* Size of the output buffer */
  unsigned char *Buf;   /* Output buffer. Maybe too big for stack so use heap */
#if MPI
  MPI_Comm comm;	/* Communicator: the group of all processes active in this device */
  int root;		/* One processor chosen to write the header and to pre-fill */
  MPI_Datatype filetype;
  long global_size;	/*  Size of the pre-filling buffer */
  unsigned char *Prefill; /* Pre-filling buffer. Maybe too big for stack so use heap */
#endif
  char debugname[MAXPATH];
  FILE *debug;
  int debugWriter;
} STR;

RUN_HEAD(byteout)
{
  DEVICE_CONST(int,append);
  DEVICE_VAR(sequence,file);
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  DEVICE_CONST(real,u0);
  DEVICE_CONST(real,u1);
  DEVICE_CONST(int,bgbyte);
  DEVICE_CONST(long,local_size);
  DEVICE_ARRAY(unsigned char,Buf);
#if MPI
  DEVICE_CONST(int, root);
  DEVICE_CONST(MPI_Comm, comm);
  DEVICE_CONST(MPI_Datatype,filetype);
  DEVICE_CONST(long,global_size);
  DEVICE_ARRAY(unsigned char,Prefill);
#endif
  DEVICE_CONST(int, debugWriter)
  DEVICE_CONST(FILE *,debug)
  long v, x, y, z;
  long nv = (v1-v0) + 1;
  long nx = (s.x1-s.x0) + 1;
  long ny = (s.y1-s.y0) + 1;
  long nz = (s.z1-s.z0) + 1;
  /* The order of bytes in the file:                        */
  /* v is the fastest, then x, ... and z is the slowest     */
#define BUF(x,y,z,v) (Buf[(v)+nv*((x)+nx*((y)+ny*(z)))])
  
  if (!append) thisq(file);
  if (!file->f) return 1;
	
#if MPI
  /* Preallocate file size */
  MPIDO(MPI_File_preallocate(file->f,global_size),"Couldn't preallocate file");

  /* Only one process shall prefill the raster, if needed */
  if (mpi_rank==root && num_empty_subdoms>0) {
    MPIDO(MPI_File_set_view(file->f,0,MPI_CHAR,filetype,"native",MPI_INFO_NULL),"Couldn't set view for prefill");
    MPIDO(MPI_File_write(file->f,Prefill,global_size,MPI_CHAR,MPI_STATUS_IGNORE),"Couldn't prefill the raster\n");
  }	

  /*  Set view for collective writing (again) */
  MPIDO(MPI_File_set_view(file->f,0,MPI_CHAR,filetype,"native",MPI_INFO_NULL),"Couldn't set view for writing");
#endif

  for(z=0;z<nz;z++){
    for(y=0;y<ny;y++) { 
      for(x=0;x<nx;x++) {
	if (!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)) {
	  for (v=0;v<nv;v++)
	    BUF(x,y,z,v) = (unsigned char) Byte(s.x0+x,s.y0+y,s.z0+z,v0+v,u0,u1);
	} else { /*  Void. Use background colour. */
	  for (v=0;v<nv;v++)
	    BUF(x,y,z,v) = (unsigned char) bgbyte;
	} /*  if isTissue else */
      } /*  for x */
    } /*  for y */
  } /*  for z */

#if MPI
  MPIDO(MPI_FILE_WRITE(file->f,Buf,local_size,MPI_CHAR,MPI_STATUS_IGNORE),"Couldn't write raster to file\n");
  if (debug && debugWriter) 
    fprintf(debug,"byteout [%ld:%ld]x[%ld:%ld]x[%ld:%ld]x[%ld:%ld] to %s at t=%ld\n",
	    (long)s.global_x0, (long)s.global_x1,
	    (long)s.global_y0, (long)s.global_y1,
	    (long)s.global_z0, (long)s.global_z1,
	    (long)v0, (long)v1, file->name, t);
#else
  fwrite(Buf,sizeof(unsigned char),local_size,file->f);
  if (append) FFLUSH(file->f);
  if (debug && debugWriter) 
    fprintf(debug,"byteout [%ld:%ld]x[%ld:%ld]x[%ld:%ld]x[%ld:%ld] to %s at t=%ld\n",
	    (long)s.x0, (long)s.x1, (long)s.y0, (long)s.y1, (long)s.z0, (long)s.z1, (long)v0, (long)v1, file->name, t);
#endif
  
  if (!append) nextq(file);
}
RUN_TAIL(byteout)

/*******************/
DESTROY_HEAD(byteout)
{
  FREE(S->Buf);
#if MPI
  if (S->file.f) MPI_File_close(&(S->file.f)); 
  if (S->Prefill) FREE(S->Prefill);
#else
  SAFE_CLOSE(S->file.f);
#endif
  SAFE_CLOSE(S->debug);
  /* Erase the last file if empty */
  /* Process 0 should be ok for that even if it is not 'runHere' */
  if (mpi_rank==0) {
    long filesize;
    FILE *f;
    /* Using sequential approach */
    f=fopen(S->file.name,"r");
    fseek(f,0,SEEK_END);
    filesize=ftell(f);
    fclose(f);
    if (filesize==0) unlink(S->file.name);
  }
  S->file.f=NULL;
} DESTROY_TAIL(byteout)

/********************/
CREATE_HEAD(byteout)
{
  ACCEPTI(append,0,0,1);
  ACCEPTQ(file,S->append?"ab":"wb",NULL);
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT(v1>=v0);
  ACCEPTR(u0,RNONE,RNONE,RNONE);
  ACCEPTR(u1,RNONE,RNONE,RNONE);
  ASSERT(u1!=u0);
  if (GEOMETRY_ON) {
    ACCEPTI(bgbyte, 255, 0, 255);
  } else if ( find_key("bgbyte=",w) ) {
    MESSAGE("Background byte parameter bgbyte is only used when geometry is active.\n"
	    "\tThe value provided will be ignored.");
    S->bgbyte=0;
  }

  
  /*  Dimensions of the local space. */
  int local_vsize = (v1 - v0) + 1;
  int local_xsize = (dev->s.x1 - dev->s.x0) + 1;
  int local_ysize = (dev->s.y1 - dev->s.y0) + 1;
  int local_zsize = (dev->s.z1 - dev->s.z0) + 1;
  unsigned char *p;
  S->local_size = local_vsize*local_xsize*local_ysize*local_zsize;
  CALLOC(S->Buf, S->local_size, sizeof(unsigned char));
  for (p=S->Buf; p<S->Buf+S->local_size; p++) *p=S->bgbyte;
 
#if MPI
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    ABORT("Could not create communicator.\n");

  if (dev->s.runHere) {
    #define NDIMS 4
    int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
    MPI_Datatype filetype;

    /* NB: the legacy sequential version stores bytes in "Fortran" order, */
    /* the MPI implementation uses MPI_ORDER_C, */
    /* hence we need to reverse the coordinate order */
    /* similar to how it is done in ppmout. */
    
    /*  Dimensions of the global space. */
    int vsize = (v1 - v0) + 1;
    int xsize = (dev->s.global_x1 - dev->s.global_x0) + 1;
    int ysize = (dev->s.global_y1 - dev->s.global_y0) + 1;
    int zsize = (dev->s.global_z1 - dev->s.global_z0) + 1;
    S->global_size = vsize*xsize*ysize*zsize;
    sizes[0] = zsize;
    sizes[1] = ysize;
    sizes[2] = xsize;
    sizes[3] = vsize;
    /* Prefilling will be done only by root */
    if (num_empty_subdoms>0 && mpi_rank==S->root) {
      MALLOC(S->Prefill, S->global_size);
      for (p=S->Prefill; p<S->Prefill+S->global_size; p++) *p=S->bgbyte;
    } else {
      /* Zero pointer signals that this process does not prefill */
      S->Prefill=NULL;
    }
    
    /*  Dimensions of the local space. */
    subsizes[0] = local_zsize;
    subsizes[1] = local_ysize;
    subsizes[2] = local_xsize;
    subsizes[3] = local_vsize;

    /* Position of the local space in the global space */
    starts[0] = dev->s.z0 - dev->s.global_z0;
    starts[1] = dev->s.y0 - dev->s.global_y0;
    starts[2] = dev->s.x0 - dev->s.global_x0;
    starts[3] = 0;

    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_CHAR,&filetype),"Couldn't define filetype.");
    MPIDO(MPI_Type_commit(&filetype),"Couldn't commit filetype.");

    S->filetype = filetype;
    
    #undef NDIMS
  }/*  if runHere */
#endif

#if MPI
  /* Identifying the root process allows only one process to write debug output. */
  if (dev->s.nowhere) {
    S->debugWriter = (mpi_rank == 0);
  } else {
    S->debugWriter = (mpi_rank == getRankContainingPoint(dev->s.global_x0,dev->s.global_y0,dev->s.global_z0));
  }
#else
  S->debugWriter = 1;
#endif
  ACCEPTF(debug,"wt","");
}
CREATE_TAIL(byteout,0)
