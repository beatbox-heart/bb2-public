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

typedef struct {
	int append;
	sequence file;
	int r, g, b;		/* layers wherefrom to extract the colour components */
	int bgr,bgg,bgb;
	real r0,r1;			/* range of values to transfrom node to a byte for the red comp. */
	real g0,g1;			/* -||- red comp. */
	real b0,b1;			/* -||- blue comp. */
#if MPI
	MPI_Datatype pointType;
	MPI_Datatype filetype;
	int count;		/*  Number of pointTypes in local space. */
	int xsize;		/*  x-dimension of global space */
	int ysize;		/*  y-dimension of global space */
	int zsize;		/*  z-dimension of global space */
	int local_xsize;	/*  x-dimension of local space */
	int local_ysize;	/*  y-dimension of local space */
	int local_zsize;	/*  z-dimension of local space */
#endif
} STR;

#define MAXCHAR 255
RUN_HEAD(vtkout2)
	DEVICE_CONST(int,append)
	DEVICE_VAR(sequence,file)
	DEVICE_CONST(int,r) DEVICE_CONST(real,r0) DEVICE_CONST(real,r1)
	DEVICE_CONST(int,g) DEVICE_CONST(real,g0) DEVICE_CONST(real,g1)
	DEVICE_CONST(int,b) DEVICE_CONST(real,b0) DEVICE_CONST(real,b1)
	DEVICE_CONST(int,bgr)
	DEVICE_CONST(int,bgg)
	DEVICE_CONST(int,bgb)
	
	int x, y, z;
	int nx = (s.x1-s.x0) + 1; // the value of this during MPI run is not good. do not use this in mpi
	int ny = (s.y1-s.y0) + 1;
	int nz = (s.z1-s.z0) + 1;
	char buf_seq[nz][ny][nx];
 #if MPI
	DEVICE_CONST(MPI_Datatype,filetype)
	DEVICE_CONST(MPI_Datatype,pointType)
	DEVICE_CONST(int,count)
	DEVICE_CONST(int,xsize)
	DEVICE_CONST(int,ysize)
	DEVICE_CONST(int,zsize)
	DEVICE_CONST(int,local_xsize)
	DEVICE_CONST(int,local_ysize)
	DEVICE_CONST(int,local_zsize)	
	char buf[local_zsize][local_ysize][local_xsize]; /*  3 is for RGB. */
	char header[1024];
	int headerLength;
	
	/*  For displacement */
	MPI_Aint lb;
	MPI_Aint char_extent;
	MPI_Offset headerSkip;
	MPI_Type_get_extent(MPI_CHAR,&lb, &char_extent);
	MPI_Offset prealloc_size;
#endif
	if (!file->f) return 1;
        DEBUG("");
	if (!append) nextq(file); DEBUG("next finished ok\n");
	
#if MPI
snprintf(header,1024,"# vtk DataFile Version 2.0\n\
vtk output\n\
BINARY\n\
DATASET STRUCTURED_POINTS\n\
DIMENSIONS %d %d %d\n\
SPACING 1 1 1\n\
ORIGIN 0 0 0\n\
POINT_DATA %d\n\
SCALARS ImageFile char 1\n\
LOOKUP_TABLE default\n", xsize,ysize,zsize,xsize*ysize*zsize);
headerLength = strlen(header);
	
	/* Preallocate file size */
	prealloc_size = headerLength + xsize*ysize*zsize;
	MPI_File_preallocate(file->f,prealloc_size);
	
	/*  Set view with no displacement (top of file). */
	mpi_errno = MPI_File_set_view(file->f, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
	CHECK_MPI_SUCCESS("Couldn't set view for header.")
	
	/*  First process in the space writes the header. */
	if(s.x0 == s.global_x0 && s.y0 == s.global_y0 && s.z0 == s.global_z0){
		mpi_errno = MPI_File_write(file->f, header, headerLength, MPI_CHAR, MPI_STATUS_IGNORE);
		CHECK_MPI_SUCCESS("Couldn't write to file.")
	}
	
	/*  Set view for collective writing. Displace by length of header. */
	headerSkip = char_extent * headerLength;
	mpi_errno = MPI_File_set_view(file->f, headerSkip, pointType, filetype, "native", MPI_INFO_NULL);
	CHECK_MPI_SUCCESS("Couldn't set view for writing.")
#else
//	fprintf (file->f,"P6\n%d %d\n%d\n",nx,ny,MAXCHAR);
		fprintf(file->f, "# vtk DataFile Version 2.0\n");
		fprintf(file->f, "vtk output\n");
		fprintf(file->f, "BINARY\n");
		fprintf(file->f, "DATASET STRUCTURED_POINTS\n");
		fprintf(file->f, "DIMENSIONS %d %d %d\n",nx,ny,nz); // N+1 columns, M+1+header rows
		fprintf(file->f, "SPACING 1 1 1\n");
		fprintf(file->f, "ORIGIN 0 0 0\n");
		fprintf(file->f, "POINT_DATA %d\n",nx*ny*nz); // total number of points on structured data.
		fprintf(file->f, "SCALARS ImageFile char 1\n");
		fprintf(file->f, "LOOKUP_TABLE default\n");
#endif

#if MPI
	for(z=0;z<local_zsize;z++) {
		for(y=local_ysize-1;y>=0;y--) {
			for(x=0;x<local_xsize;x++) {
				if(!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)){
					buf[z][y][x] = Byte(s.x0+x,s.y0+y,s.z0+z,r,r0,r1);
				}else{ /*  Void. Use background colour. */
					buf[z][y][x] = (unsigned) bgr;
				} /*  else */
			} /*  for x */
		} /*  for y */
	} /*  for z */

/*
MPI File Write is not the same as fwrite in serial.
*/
	if(blocking_or_contiguous==0)
	mpi_errno = MPI_File_write(file->f, buf, count, pointType, MPI_STATUS_IGNORE);
	else
	mpi_errno = MPI_File_write_all(file->f, buf, count, pointType, MPI_STATUS_IGNORE);
	CHECK_MPI_SUCCESS("Could not write to file.")
#else
	/*  Sequential version */
	for(z=0;z<nz;z++){
		for(y=0;y<ny;y++) { 
			for(x=0;x<nx;x++) {
				if(!GEOMETRY_ON || isTissue(s.x0+x,s.y0+y,s.z0+z)){
					buf_seq[z][y][x] = Byte(s.x0+x,s.y0+y,s.z0+z,r,r0,r1);
				}else{ /*  Void. Use background colour. */
					buf_seq[z][y][x] = (unsigned)bgr;
				}
			} /*  for x */
		} /*  for y */
	} /*  for z */

fwrite(buf_seq,sizeof(char),nx*ny*nz,file->f);

#endif
RUN_TAIL(vtkout2)

DESTROY_HEAD(vtkout2)
#if MPI
  if (S->file.f) MPI_File_close(&(S->file.f)); S->file.f=NULL;
#else
  SAFE_CLOSE(S->file.f);
#endif
DESTROY_TAIL(vtkout2)

CREATE_HEAD(vtkout2) {
	
	ACCEPTI(append,0,0,1);
	ACCEPTQ(file,S->append?"at":"wt",NULL);
	ACCEPTI(r,INONE,-1,(int)vmax-1); ACCEPTR(r0,RNONE,RNONE,RNONE); ACCEPTR(r1,RNONE,RNONE,RNONE); ASSERT(S->r0!=S->r1);
	ACCEPTI(g,INONE,-1,(int)vmax-1); ACCEPTR(g0,RNONE,RNONE,RNONE); ACCEPTR(g1,RNONE,RNONE,RNONE); ASSERT(S->g0!=S->g1);
	ACCEPTI(b,INONE,-1,(int)vmax-1); ACCEPTR(b0,RNONE,RNONE,RNONE); ACCEPTR(b1,RNONE,RNONE,RNONE); ASSERT(S->b0!=S->b1);
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

#if MPI
 if(dev->s.runHere){
	
	/*  Point type specific to ppmout (i.e. 3 chars per pixel). */
	MPI_Datatype pointType;
	mpi_errno = MPI_Type_contiguous(1,MPI_CHAR,&pointType);
	CHECK_MPI_SUCCESS("Couldn't create pixel type (pointType).")
	mpi_errno = MPI_Type_commit(&pointType);
	CHECK_MPI_SUCCESS("Couldn't commit pixel type (pointType).")
	
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
	
	mpi_errno = MPI_Type_create_subarray(NDIMS, sizes, subsizes, starts, MPI_ORDER_C, pointType, &filetype);
	CHECK_MPI_SUCCESS("Couldn't define filetype.")
	mpi_errno = MPI_Type_commit(&filetype);
	CHECK_MPI_SUCCESS("Couldn't commit filetype.")
	
	S->pointType = pointType;
	S->filetype = filetype;
	S->count = count;
	S->xsize = xsize;
	S->ysize = ysize;
	S->zsize = zsize;
	S->local_xsize = local_xsize;
	S->local_ysize = local_ysize;
	S->local_zsize = local_zsize;
 }/*  if s.runHere */
#endif
	
} CREATE_TAIL(vtkout2,0)

#undef R
#undef G
#undef B
