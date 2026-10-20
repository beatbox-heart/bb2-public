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
SK, Sept. 2012
Produce a VTK file with Bytes this time.
3 October 2012: 
It now works for the 2D FHN with 50 x 50 grid.
The 50 x 50 grid was hard coded to allow testing
and finishing vtkout2.c

To do:
1. Generalise the parallel 2D device for FHN
2. Generalise it to 3D
*/


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

/* extern FILE *debug; */

/* #define DEBUG if (0) MESSAGE */

#define R 0
#define G 1
#define B 2

#define BUF(z,y,x) (Buf[(x)+local_xsize*((y)+local_ysize*(z))])

typedef struct {
	int append;
	sequence file;
	int r, g, b;		/* layers wherefrom to extract the colour components */
	int bgr,bgg,bgb;
	real r0,r1;			/* range of values to transfrom node to a byte for the red comp. */
	real g0,g1;			/* -||- red comp. */
	real b0,b1;			/* -||- blue comp. */
	int local_xsize;	/*  x-dimension of local space */
	int local_ysize;	/*  y-dimension of local space */
	int local_zsize;	/*  z-dimension of local space */
	unsigned char *Buf;	/*  output buffer */
#if MPI
	MPI_Comm comm;		/* Communicator: the group of all processes active in this device */
	int root;		/* One processor chosen to write the header and to pre-fill */
	MPI_Datatype pointType;
	MPI_Datatype filetype;
	int count;		/*  Number of pointTypes in local space. */
	unsigned char *prefill;	/*  buffer for pre-filling the files */
	int raster_size;	/*  byte size of the pre-filling buffer */
	MPI_Offset offset;	/*  where the next image starts in the file */
#endif
} STR;

/* Only the r layer is output, as a scalar field; g and b are not used */
#define MAXCHAR 255
RUN_HEAD(vtkout2)
	DEVICE_CONST(int,append)
	DEVICE_VAR(sequence,file)
	DEVICE_CONST(int,r) DEVICE_CONST(real,r0) DEVICE_CONST(real,r1)
	DEVICE_CONST(int,bgr)
	DEVICE_CONST(int,local_xsize)
	DEVICE_CONST(int,local_ysize)
	DEVICE_CONST(int,local_zsize)
	/* The buffer may be too big for the stack, hence go for the heap instead */
	DEVICE_ARRAY(unsigned char,Buf)
	int x, y, z;
	int xsize = (s.global_x1-s.global_x0) + 1;
	int ysize = (s.global_y1-s.global_y0) + 1;
	int zsize = (s.global_z1-s.global_z0) + 1;
	char header[1024];
	int headerLength;
#if MPI
	DEVICE_CONST(int,root)
	DEVICE_CONST(MPI_Datatype,filetype)
	DEVICE_CONST(MPI_Datatype,pointType)
	DEVICE_CONST(int,count)
	DEVICE_ARRAY(unsigned char,prefill)
	DEVICE_CONST(int,raster_size)
	MPI_Offset headerSkip;
#endif
	if (!append) thisq(file);
	if (!file->f) return 1;

	snprintf(header,1024,
		"# vtk DataFile Version 2.0\n"
		"vtk output\n"
		"BINARY\n"
		"DATASET STRUCTURED_POINTS\n"
		"DIMENSIONS %d %d %d\n"
		"SPACING 1 1 1\n"
		"ORIGIN 0 0 0\n"
		"POINT_DATA %d\n"
		"SCALARS ImageFile unsigned_char 1\n"
		"LOOKUP_TABLE default\n",
		xsize,ysize,zsize,xsize*ysize*zsize);
	headerLength = strlen(header);

	/* The order of bytes in the file: x is the fastest, then y, and z is the slowest */
	for(z=0;z<local_zsize;z++) {
		for(y=0;y<local_ysize;y++) {
			for(x=0;x<local_xsize;x++) {
				if(!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)){
					BUF(z,y,x) = (unsigned char) Byte(s.x0+x,s.y0+y,s.z0+z,r,r0,r1);
				}else{ /*  Void. Use background colour. */
					BUF(z,y,x) = (unsigned char) bgr;
				} /*  else */
			} /*  for x */
		} /*  for y */
	} /*  for z */

#if MPI
	/* Set file size; MPI_File_preallocate hangs intermittently in Open MPI 5.0 OMPIO */
	MPIDO(MPI_File_set_size(file->f,S->offset+headerLength+raster_size),"Couldn't set file size");

	/* Header goes from the start of the image */
	MPIDO(MPI_File_set_view(file->f,S->offset,MPI_CHAR,MPI_CHAR,"native",MPI_INFO_NULL),"Couldn't set view for header.");

	/* Only one process shall write the header */
	if (mpi_rank==root) {
		MPIDO(MPI_File_write(file->f,header,headerLength,MPI_CHAR,MPI_STATUS_IGNORE),"Couldn't write header to file.");
		/* and prefill the raster, if needed */
		if (prefill) {
			MPIDO(MPI_File_write(file->f,prefill,raster_size,MPI_CHAR,MPI_STATUS_IGNORE),"Couldn't prefill the raster\n");
		}
	}

	/*  Set view for collective writing. Displace by length of header. */
	headerSkip = S->offset + headerLength;
	MPIDO(MPI_File_set_view(file->f,headerSkip,pointType,filetype,"native",MPI_INFO_NULL),"Couldn't set view for writing.");
	MPIDO(MPI_FILE_WRITE(file->f,Buf,count,pointType,MPI_STATUS_IGNORE),"Couldn't write raster to file\n");
	/* In append mode, the images follow one another in the same file */
	if (append) S->offset += headerLength+raster_size;
#else
	fputs(header,file->f);
	fwrite(Buf,sizeof(unsigned char),(size_t)local_xsize*local_ysize*local_zsize,file->f);
	if (append) FFLUSH(file->f);
#endif
	if (!append) nextq(file);
RUN_TAIL(vtkout2)

DESTROY_HEAD(vtkout2)
#if MPI
	if (S->file.f) nextq(&(S->file)); /* append mode */
	if (S->prefill) FREE(S->prefill);
#else
	if (S->file.f) nextq(&(S->file)); /* append mode */
#endif
	if (S->Buf) FREE(S->Buf);
	S->file.f=NULL;
DESTROY_TAIL(vtkout2)

CREATE_HEAD(vtkout2)
{
	DEVICE_IS_RECTANGULAR;
	
	ACCEPTI(append,0,0,1);
	ACCEPTQ(file,S->append?"at":"wt",NULL);
	ACCEPTI(r,INONE,-1,(int)vmax-1); ACCEPTR(r0,RNONE,RNONE,RNONE); ACCEPTR(r1,RNONE,RNONE,RNONE); ASSERT(S->r0!=S->r1);
	/* g and b are not used; they are optional, for syntax compatibility with ppmout */
	ACCEPTI(g,-1,-1,(int)vmax-1); ACCEPTR(g0,0,RNONE,RNONE); ACCEPTR(g1,1,RNONE,RNONE);
	ACCEPTI(b,-1,-1,(int)vmax-1); ACCEPTR(b0,0,RNONE,RNONE); ACCEPTR(b1,1,RNONE,RNONE);
	/**
	 *	Background colour components.
	 *	Used only when geometry is active, to show void points.
	 *	If provided when geometry is off, the values will be ignored.
	**/
	if(GEOMETRY_ON){
		ACCEPTI(bgr, 255, 0, 255); ACCEPTI(bgg, 255, 0, 255); ACCEPTI(bgb, 255, 0, 255);
	}else if ( find_key("bgr=",w) || find_key("bgg=",w) || find_key("bgb=",w) ){
			MESSAGE("Background colour parameters (bgr, bgg, bgb) are only used when geometry is active.\n\tThe value(s) provided will be ignored.");
	}

	/*  Local buffer dimensions */
	S->local_xsize = (dev->s.x1 - dev->s.x0) + 1;
	S->local_ysize = (dev->s.y1 - dev->s.y0) + 1;
	S->local_zsize = (dev->s.z1 - dev->s.z0) + 1;
	S->Buf = NULL;

#if MPI
	S->prefill = NULL;
	if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &(S->comm), &(S->root)))
		ABORT("Could not create communicator.\n");

 if(dev->s.runHere){
	S->offset=0;
	if (S->append && S->file.f) MPIDO(MPI_File_get_size(S->file.f,&(S->offset)),"Couldn't get file size");
	
	/*  Point type: 1 byte per point. */
	MPI_Datatype pointType;
	MPIDO(MPI_Type_contiguous(1,MPI_CHAR,&pointType),"Couldn't create pixel type (pointType).");
	MPIDO(MPI_Type_commit(&pointType),"Couldn't commit pixel type (pointType).");
	
	/*  Generic 3D filetype based on space. */
	MPI_Datatype filetype;
	
	#define NDIMS 3
	int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
	
	sizes[0] = (dev->s.global_z1 - dev->s.global_z0) + 1;
	sizes[1] = (dev->s.global_y1 - dev->s.global_y0) + 1;
	sizes[2] = (dev->s.global_x1 - dev->s.global_x0) + 1;

	subsizes[0] = S->local_zsize;
	subsizes[1] = S->local_ysize;
	subsizes[2] = S->local_xsize;
	
	starts[0] = dev->s.z0 - dev->s.global_z0;
	starts[1] = dev->s.y0 - dev->s.global_y0;
	starts[2] = dev->s.x0 - dev->s.global_x0;
	
	MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,pointType,&filetype),"Couldn't define filetype.");
	MPIDO(MPI_Type_commit(&filetype),"Couldn't commit filetype.");
	
	S->pointType = pointType;
	S->filetype = filetype;
	S->count = S->local_xsize*S->local_ysize*S->local_zsize;
	S->raster_size = sizes[0]*sizes[1]*sizes[2];
	#undef NDIMS

	/* Prefilling will be done only by root */
	if (num_empty_subdoms>0 && mpi_rank==S->root) {
		unsigned char *p;
		MALLOC(S->prefill,S->raster_size);
		for (p=S->prefill;p<S->prefill+S->raster_size;p++) *p=(unsigned char) (S->bgr);
	}

	CALLOC(S->Buf,1L*S->local_zsize*S->local_ysize*S->local_xsize,sizeof(unsigned char));
 }/*  if s.runHere */
#else
	CALLOC(S->Buf,1L*S->local_zsize*S->local_ysize*S->local_xsize,sizeof(unsigned char));
#endif
	
} CREATE_TAIL(vtkout2,0)

#undef R
#undef G
#undef B
