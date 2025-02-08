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
SK
5/4/2013
A morph of record.c as in the SVN directory.
This code replaces the x!=s.x1, y!=s.y1, ... code with
x!=xmax-1, y=!ymax-1, ...
Doing so facilitates visualising 1D results using gnuplot.
At each output time, and when the output is appended to a file,
it appears as a row entry from x=0 to x=xmax-1. That allows the
use of gnuplot pm3d.
*/

/* Record values in the given 4D cube to a text file */
#if MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "system.h"
#include <limits.h>

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "mpi_io_choice.h"

/* 8< mpi_record_parameter_structure */
#define RECORD_SEPARATOR_LENGTH 2
typedef struct {
  int v0, v1;
	char filename[MAXPATH];
#if MPI
	MPI_File file;
	MPI_Datatype V_Type;
	int value_length;
#else
	FILE *file;
#endif
	int append;
	char filehead[80];
	char vsep[RECORD_SEPARATOR_LENGTH +1];
	char zsep[RECORD_SEPARATOR_LENGTH +1];
	char ysep[RECORD_SEPARATOR_LENGTH +1];
	char xsep[RECORD_SEPARATOR_LENGTH +1];
	char recordsep[RECORD_SEPARATOR_LENGTH +1];
} STR;
/* >8 mpi_record_parameter_structure */

RUN_HEAD(skrecord)
  DEVICE_CONST(int, v0);
  DEVICE_CONST(int, v1);
	DEVICE_ARRAY(char,vsep)
	DEVICE_ARRAY(char,zsep)
	DEVICE_ARRAY(char,ysep)
	DEVICE_ARRAY(char,xsep)
	DEVICE_ARRAY(char,recordsep)
	int v, z, y, x;
	char *nextSep;

#if MPI
	DEVICE_CONST(MPI_File, file)
	DEVICE_CONST(MPI_Datatype, V_Type)
	int buf_length = S->value_length + 1;
	char buf[buf_length];
#else
	DEVICE_CONST(FILE *,file)
	if (!file) return 1;
        int buf_length = MAXSTRLEN;
#endif
/* 8< mpi_record_writing_to_file */
/*  Write data, selecting hierarchical separator */
/*  Record->X->Y->Z->V (most->least significant). */

	for(x=s.x0; x<=s.x1; x++){
		for(y=s.y0; y<=s.y1; y++){
			for(z=s.z0; z<=s.z1; z++){
				for(v=v0; v<=v1; v++){
					if(v!=v1){
						nextSep = vsep;
					}else{
						if(z!=s.z1){
							nextSep = zsep;
						}else{
							if(y!=s.y1){ /* This sets up the nextSep for the corner */
								nextSep = ysep;
							}else{
								if(x!=(xmax-1)){ /* To facilitate use with gnuplot. use record for usual record */
									nextSep = xsep;
								}else{
									nextSep = recordsep; /* the default recordsep is a newline char */
								}
							}
						}
					}
					snprintf(buf, buf_length, REALF_IO, New[ind(x,y,z,v)]);
					strcat(buf, nextSep);
#if MPI
					mpi_errno = MPI_File_write(file, buf, 1, V_Type, MPI_STATUS_IGNORE);
					CHECK_MPI_SUCCESS("Couldn't write to file.")
#else
					fprintf(file, "%s", buf);
#endif
				}
			}
		}
	}
#if MPI==0
	if(fflush_or_not==1)
	fflush(file);
#endif
RUN_TAIL(skrecord)
/* >8 mpi_record_writing_to_file */

DESTROY_HEAD(skrecord)
#if MPI
	if(s.runHere){
		MPI_File_close(&(S->file));
	}
#else
  SAFE_CLOSE(S->file);
#endif
DESTROY_TAIL(skrecord)


/* ----------------------------------------------------------------- */
CREATE_HEAD(skrecord)
  ACCEPTI(v0,INONE,0,vmax-1);
  ACCEPTI(v1,INONE,0,vmax-1);
  ASSERT(v1 != v0);
	ACCEPTI(append,1,0,1);
	ACCEPTS(filehead,"");
	 /* Add newline to filehead, only if one exists,
	  * otherwise we'll have an empty line at the top of the file. */
	if(strcmp(S->filehead,"")!=0){
		strcat(S->filehead, "\n");
	}
	MESSAGE("Filehead = %s\n",S->filehead);
	
	#define CHECK_SEPARATOR(s)																															\
	if(strlen(S->s) > RECORD_SEPARATOR_LENGTH){ 																										\
		MESSAGE("WARNING from %s: %s is longer than the required separator size of %d and will be truncated.\n", dev->n, #s, RECORD_SEPARATOR_LENGTH);	\
	}else if(strlen(S->s) < RECORD_SEPARATOR_LENGTH){																									\
		int i;																																			\
		for(i=strlen(S->s); i < RECORD_SEPARATOR_LENGTH; i++){																						\
			S->s[i] = ' ';																																\
		}																																				\
	}																																					\
	S->s[RECORD_SEPARATOR_LENGTH] = '\0';

/* changed default xsep from \n to empty space char. */

	  ACCEPTS(vsep,"  "); CHECK_SEPARATOR(vsep)
	  ACCEPTS(zsep,", "); CHECK_SEPARATOR(zsep)
	  ACCEPTS(ysep,"; "); CHECK_SEPARATOR(ysep)
	  ACCEPTS(xsep," "); CHECK_SEPARATOR(xsep)
	  ACCEPTS(recordsep," \n"); CHECK_SEPARATOR(recordsep)
	
	#undef CHECK_SEPARATOR

#if MPI
	if (!accepts("file=",&(S->filename[0]),"",w)) return(0);

	/*  Temporary Datatypes for file I/O */
	MPI_Datatype V_Type, Point_Type, Z_Type, Y_Type, X_Type;

	/*  For defining MPI datatypes */
	int *blocklengths = calloc(3, sizeof(int));
	MPI_Aint *disp = calloc(3, sizeof(MPI_Aint));
	MPI_Datatype *old_types = calloc(3, sizeof(MPI_Datatype));

	/*  Extents */
	MPI_Aint char_extent,v_extent,point_extent,z_extent,y_extent,x_extent,lb;
	int num_vars = (dev->v1 - dev->v0) + 1;
	int count;

	/*  Lengths */
	int filehead_length = strlen(S->filehead);
	int value_length = REALF_IO_WIDTH + RECORD_SEPARATOR_LENGTH;
	MPI_Type_get_extent(MPI_CHAR,&lb,&char_extent);

	/*  Represents a single value with separator. */
	count = value_length;
	mpi_errno = MPI_Type_contiguous(count,MPI_CHAR,&V_Type);
	CHECK_MPI_SUCCESS("Couldn't create V_Type.")
	mpi_errno = MPI_Type_commit(&V_Type);
	CHECK_MPI_SUCCESS("Couldn't commit V_Type.")
	MPI_Type_get_extent(V_Type, &lb, &v_extent);

	/*  Represents one point, i.e. num_vars values with separators. */
	count = num_vars;
	mpi_errno = MPI_Type_contiguous(count,V_Type,&Point_Type);
	CHECK_MPI_SUCCESS("Couldn't create Point_Type.")
	mpi_errno = MPI_Type_commit(&Point_Type);
	CHECK_MPI_SUCCESS("Couldn't commit Point_Type.")
	MPI_Type_get_extent(Point_Type, &lb, &point_extent);

	/*  Struct representing local contribution to ZV plane */
	blocklengths[0] = 1;
	blocklengths[1] = (dev->s.z1 - dev->s.z0) + 1;
	blocklengths[2] = 1;

	disp[0] = 0;
	disp[1] = (dev->s.z0 - dev->s.global_z0) * point_extent;
	disp[2] = ((dev->s.global_z1 - dev->s.global_z0) + 1) * point_extent;

	old_types[0] = MPI_LB;
	old_types[1] = Point_Type;
	old_types[2] = MPI_UB;

	mpi_errno = MPI_Type_create_struct(3, blocklengths, disp, old_types, &Z_Type);
	CHECK_MPI_SUCCESS("Couldn't create Z_Type.")
	mpi_errno = MPI_Type_commit(&Z_Type);
	CHECK_MPI_SUCCESS("Couldn't commit Z_Type.")
	MPI_Type_get_extent(Z_Type, &lb, &z_extent);

	/*  Struct representing local contribution to YZV hyperplane */
	blocklengths[0] = 1;
	blocklengths[1] = (dev->s.y1 - dev->s.y0) + 1;
	blocklengths[2] = 1;

	disp[0] = 0;
	disp[1] = (dev->s.y0 - dev->s.global_y0) * z_extent;
	disp[2] = ((dev->s.global_y1 - dev->s.global_y0) + 1) * z_extent;

	old_types[0] = MPI_LB;
	old_types[1] = Z_Type;
	old_types[2] = MPI_UB;

	mpi_errno = MPI_Type_create_struct(3, blocklengths, disp, old_types, &Y_Type);
	CHECK_MPI_SUCCESS("Couldn't create Y_Type.")
	mpi_errno = MPI_Type_commit(&Y_Type);
	CHECK_MPI_SUCCESS("Couldn't commit Y_Type.")
	MPI_Type_get_extent(Y_Type, &lb, &y_extent);

	/*  Struct representing local contribution to XYZV hypercube */
	blocklengths[0] = 1;
	blocklengths[1] = (dev->s.x1 - dev->s.x0) + 1;
	blocklengths[2] = 1;

	disp[0] = 0;
	disp[1] = (dev->s.x0 - dev->s.global_x0) * y_extent;
	disp[2] = ((dev->s.global_x1 - dev->s.global_x0) + 1) * y_extent;

	old_types[0] = MPI_LB;
	old_types[1] = Y_Type;
	old_types[2] = MPI_UB;

	mpi_errno = MPI_Type_create_struct(3, blocklengths, disp, old_types, &X_Type);
	CHECK_MPI_SUCCESS("Couldn't create X_Type.")
	mpi_errno = MPI_Type_commit(&X_Type);
	CHECK_MPI_SUCCESS("Couldn't commit X_Type.")
	MPI_Type_get_extent(X_Type, &lb, &x_extent);
	/*---------------------------------------------------------------------*/

/* 8< record_initialising_file */
	/*  Communicator for active dump instances. */
	MPI_Comm comm_record;
	int firstRank; /*  Rank of first processor in communicator, for writing headers. */

	/*  File */
	MPI_File file;
	int mode;
	int displacement=0;
	MPI_Offset filesize; /*  existing file size, for appending. */

	if(!deviceCommunicatorWithFirstRank(dev->s.runHere, &comm_record, &firstRank)){
		ABORT("Could not create communicator.");
	}
	
	if(dev->s.runHere){ /*  Active devices only!	 */
		/*  Set access mode */
		if(S->append){
			mode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND;
		}else{
			mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
		}
	
		/*  Open file and set view */
		mpi_errno = MPI_File_open(comm_record, S->filename, mode, MPI_INFO_NULL, &file);		
	
		/*  If requested by user, append to EOF */
		if(S->append){
			MPI_File_get_size(file, &filesize);
			displacement = (int)filesize;
		}
	
		/*  Write file header, if we have one. */
		if(*S->filehead){
			mpi_errno = MPI_File_set_view(file, displacement, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

			if(mpi_rank == firstRank){
				MPI_File_write(file, S->filehead, filehead_length, MPI_CHAR, MPI_STATUS_IGNORE);
			}
			displacement += filehead_length;
		}
		
		/*  Have to reset view for X_Type. */
		mpi_errno = MPI_File_set_view(file, displacement, MPI_CHAR, X_Type, "native", MPI_INFO_NULL);
	
		/*  Assign the file to the parameter structure */
		S->file = file;
		S->value_length = value_length;
		S->V_Type = V_Type;
	}

#else
	/*  Sequential version */
	ACCEPTF(file,S->append?"at":"wt",NULL);
	if(*S->filehead) {
		fprintf(S->file,"%s",S->filehead);
		if(fflush_or_not==1)
		fflush(S->file);	
	}
#endif
/* >8 record_initialising_file */
CREATE_TAIL(skrecord,0)

#undef RECORD_SEPARATOR_LENGTH
