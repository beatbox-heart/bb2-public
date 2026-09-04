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

/* Produce a PGM file */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "sequence.h"
#include "byte.h"
#include "mpi_io_choice.h"

#define BUF(z,y,x) (Buf[(x)+local_xsize*((y)+local_ysize*(z))])

typedef struct {
  int append;
  sequence file;
  int echo;
  int g;			/* layer wherefrom to extract the grayscale value */
  real g0,g1;			/* range of values to transfrom node to a byte for the gray value. */
  int bgg;			/* gray value representing void background */
#if MPI
  MPI_Comm comm;	/* Communicator: the group of all processes active in this device */
  int root;		/* One processor chosen to write the header and to pre-fill */
  MPI_Datatype pointType;
  MPI_Datatype filetype;
  int count;		/*  Number of pointTypes in local space. */
  int xsize;		/*  x-dimension of global space */
  int ysize;		/*  y-dimension of global space */
  int zsize;		/*  z-dimension of global space */
  int local_xsize;	/*  x-dimension of local space */
  int local_ysize;	/*  y-dimension of local space */
  int local_zsize;	/*  z-dimension of local space */
  unsigned char *Buf;   /*  output buffer */
  unsigned char *prefill;/*  buffer for pre-filling the files */
  int raster_size;	/*  byte size of the pre-filling buffer */
#endif
} STR;

#define MAXCHAR 255
RUN_HEAD(pgmout)
{
  DEVICE_CONST(int,append);
  DEVICE_VAR(sequence,file);
  DEVICE_CONST(int,echo);
  DEVICE_CONST(int,g) DEVICE_CONST(real,g0) DEVICE_CONST(real,g1);
  DEVICE_CONST(int,bgg);

  int x, y, z;
  int nx = (s.x1-s.x0) + 1;
  int ny = (s.y1-s.y0) + 1;
  int nz = (s.z1-s.z0) + 1;
#if MPI
  DEVICE_CONST(int, root);
  DEVICE_CONST(MPI_Comm, comm);
  DEVICE_CONST(MPI_Datatype,filetype);
  DEVICE_CONST(MPI_Datatype,pointType);
  DEVICE_CONST(int,count);
  DEVICE_CONST(int,xsize);
  DEVICE_CONST(int,ysize);
  DEVICE_CONST(int,zsize);
  DEVICE_CONST(int,local_xsize);
  DEVICE_CONST(int,local_ysize);
  DEVICE_CONST(int,local_zsize);
  DEVICE_ARRAY(unsigned char,prefill);
  DEVICE_CONST(int,raster_size);
  /* The order of bytes in the file:                        */
  /* x is the fastest, then y, and z is the slowest */
  /* unsigned char buf[local_zsize][local_ysize][local_xsize] */ 
  /* The buffer may be too big for the stack, hence go for the heap  instead */
  DEVICE_ARRAY(unsigned char,Buf);
  char header[256]; /* it is hoped that 256 is sufficient but no check is done! */
  int headerLength;
  
  /*  For displacement */
  MPI_Aint lb;
  MPI_Aint char_extent;
  MPI_Offset headerSkip;
  MPI_Type_get_extent(MPI_CHAR,&lb, &char_extent);
  MPI_Offset prealloc_size;
#endif
  if (!append) thisq(file);
  if (!file->f) return 1;
  
#if MPI
  snprintf (header,256,"P5\n%d %d\n%d\n",xsize,ysize,MAXCHAR);
  headerLength = strlen(header);
  
  /* Preallocate file size */
  prealloc_size = headerLength + raster_size;
  MPIDO(MPI_File_preallocate(file->f,prealloc_size),"Couldn't preallocate file");
  
  /* Header goes with no displacement, from top of file */
  MPIDO(MPI_File_set_view(file->f,0,MPI_CHAR,MPI_CHAR,"native", MPI_INFO_NULL),"Couldn't set view for header.");

  /* Only one process shall write the header */
  if (mpi_rank==root) {
    MPIDO(MPI_File_write(file->f,header,headerLength,MPI_CHAR,MPI_STATUS_IGNORE),"Couldn't write header to file.");
    /* and prefill the raster, if needed */
    if (num_empty_subdoms>0) { 
      MPIDO(MPI_File_write(file->f,prefill,raster_size,MPI_CHAR,MPI_STATUS_IGNORE),"Couldn't prefill the raster\n");
    }
  }
  
  /*  Set view for collective writing. Displace by length of header. */
  headerSkip = char_extent * headerLength;
  MPIDO(MPI_File_set_view(file->f,headerSkip,pointType,filetype,"native",MPI_INFO_NULL),"Couldn't set view the raster.");

  for(z=0;z<local_zsize;z++) {
    for(y=0;y<local_ysize;y++) { /* RMF code was back-to-front here, for no apparent reason */
      for(x=0;x<local_xsize;x++) {
	if (!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)) {
	  BUF(z,y,x) = (unsigned char) Byte(s.x0+x,s.y0+y,s.z0+z,g,g0,g1);
	} else { /*  Void. Use background colour. */
	  BUF(z,y,x) = (unsigned char) bgg;
	} /*  else */
      } /*  for x */
    } /*  for y */
  } /*  for z */
  MPIDO(MPI_FILE_WRITE(file->f,Buf,count,pointType,MPI_STATUS_IGNORE),"Couldn't write raster to file\n");
#else
  /*  Sequential version */
  fprintf (file->f,"P5\n%d %d\n%d\n",nx,ny,MAXCHAR);
  for(z=0;z<nz;z++){
    for(y=0;y<ny;y++) { 
      for(x=0;x<nx;x++) {
	if (!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)) {
	  putc(Byte(s.x0+x,s.y0+y,s.z0+z,g,g0,g1),file->f);
	} else { /*  Void. Use background colour. */
	  putc((unsigned)bgg,file->f);
	}
      } /*  for x */
    } /*  for y */
    if (append) FFLUSH(file->f);
  } /*  for z */
#endif
  if (echo) MESSAGE("pgmout to '%s' at t=%ld\n",file->name,t);
  if (!append) nextq(file);
}
RUN_TAIL(pgmout)

DESTROY_HEAD(pgmout) {
#if MPI
  if (S->file.f) MPI_File_close(&(S->file.f)); 
  if (S->prefill) FREE(S->prefill);
  FREE(S->Buf);
#else
  SAFE_CLOSE(S->file.f);
#endif
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
} DESTROY_TAIL(pgmout)
  
CREATE_HEAD(pgmout)
{
  DEVICE_IS_RECTANGULAR;

  ACCEPTI(append,0,0,1);
  ACCEPTQ(file,S->append?"at":"wt",NULL);
  ACCEPTI(echo,1,0,1);

  ACCEPTI(g,INONE,-1,(int)vmax-1);
  ACCEPTR(g0,RNONE,RNONE,RNONE);
  ACCEPTR(g1,RNONE,RNONE,RNONE);
  ASSERT(S->g0!=S->g1);
  
  /**
   *	Background colour components.
   *	Used only when geometry is active, to show void points.
   *	If provided when geometry is off, the values will be ignored.
   **/
  if (GEOMETRY_ON) {
    ACCEPTI(bgg, 255, 0, 255);
  } else if ( find_key("bgr=",w) || find_key("bgb=",w) ) {
    MESSAGE("Background colour parameters (bgr, bgb) are only used when geometry is active.\n"
	    "\tThe value(s) provided will be ignored.");
  }

#if MPI
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
    ABORT("Could not create communicator.\n");

  if (dev->s.runHere) {

    /*  Point type specific to pgmout (i.e. 1 byte per pixel). */
    MPI_Datatype pointType;
    MPIDO(MPI_Type_contiguous(1,MPI_CHAR,&pointType),"Couldn't create pixel type (pointType).");
    MPIDO(MPI_Type_commit(&pointType),"Couldn't commit pixel type (pointType).");
      
    /*  Generic 3D filetype based on space. */
    /*  Filetype */
    MPI_Datatype filetype;
    int count;
    
    #define NDIMS 3
    int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
    
    /*  Image dimensions. */
    int xsize = (dev->s.global_x1 - dev->s.global_x0) + 1;
    int ysize = (dev->s.global_y1 - dev->s.global_y0) + 1;
    int zsize = (dev->s.global_z1 - dev->s.global_z0) + 1;
    int raster_size = xsize*ysize*zsize*3;;
    unsigned char *prefill;
    
    /*  Local buffer dimensions */
    int local_xsize = (dev->s.x1 - dev->s.x0) + 1;
    int local_ysize = (dev->s.y1 - dev->s.y0) + 1;
    int local_zsize = (dev->s.z1 - dev->s.z0) + 1;
    
    sizes[0] = zsize;
    sizes[1] = ysize;
    sizes[2] = xsize;
    
    subsizes[0] = local_zsize;
    subsizes[1] = local_ysize;
    subsizes[2] = local_xsize;
    count = local_xsize*local_ysize*local_zsize;
    
    starts[0] = dev->s.z0 - dev->s.global_z0;
    starts[1] = dev->s.y0 - dev->s.global_y0;
    starts[2] = dev->s.x0 - dev->s.global_x0;

    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,pointType,&filetype),"Couldn't define filetype.");
    MPIDO(MPI_Type_commit(&filetype),"Couldn't commit filetype.");
      
    S->pointType = pointType;
    S->filetype = filetype;
    S->count = count;
    S->xsize = xsize;
    S->ysize = ysize;
    S->zsize = zsize;
    S->local_xsize = local_xsize;
    S->local_ysize = local_ysize;
    S->local_zsize = local_zsize;
    #undef NDIMS

    /* Prefilling will be done only by root */
    if (num_empty_subdoms>0 && mpi_rank==S->root) {
      unsigned char *p;
      MALLOC(prefill,raster_size);
      for (p=prefill;p<prefill+raster_size;p++) {
	p[0]=(unsigned char) (S->bgg);
      }
    } else {
      /* Zero pointer signals that this process does not prefill */
      prefill=NULL;
    }
    S->prefill=prefill;
    S->raster_size=raster_size;

    /* Byte buffer allocated on heap once rather than on stack every time */
    unsigned char *Buf;
    CALLOC(Buf,1L*local_zsize*local_ysize*local_xsize,sizeof(unsigned char));
    S->Buf=Buf;
    
  } /*  if s.runHere */
#endif
  
}
CREATE_TAIL(pgmout,0)

