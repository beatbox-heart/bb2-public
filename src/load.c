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
 * Load contents of a 4D subset from a binary file
 */

#include <stdlib.h>
#include <stdio.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "mpi_io_choice.h"

typedef struct {
  int v0, v1;			/* Range of layers */
  char filename[MAXPATH];
#if MPI
  MPI_File file;		/*  File handle for collective I/O. */
  MPI_Datatype sourceType;	/*  Datatype describing subarray of New to be read in. */
#else
  FILE *file;			/*  File handle for sequential I/O. */
#endif
  int rewind;
} STR;

/* "Respectful" alias, to tell apart from eponymous parameter */
void (*Rewind) (FILE *stream) = rewind;

#if MPI
#else
#define error(size) { \
    if (feof(file)) { \
      EXPECTED_ERROR("end of file in '%s' while reading %d number(s); aborting\n",S->filename,size); \
    } else { \
      EXPECTED_ERROR("error %d in '%s' while reading %d number(s); aborting\n",ferror(file),S->filename,size); \
    } \
  }
#endif

RUN_HEAD(load)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  /********************** MPI VERSION ********************************/
#if MPI
  DEVICE_CONST(MPI_File, file);
  DEVICE_CONST(MPI_Datatype, sourceType);
  DEVICE_CONST(int, rewind);
  if (!file) return 1;
  
  /*  Can I get a rewind? */
  if(rewind) MPI_File_seek(file,0,MPI_SEEK_SET);
  
  if (blocking_or_contiguous==0) {
    MPIDO(MPI_File_read(file, New, 1, sourceType, MPI_STATUS_IGNORE),"Couldn't read from file.");
  } else { /* in the hope that the read works the same as the write! */
    MPIDO(MPI_File_read_all(file, New, 1, sourceType, MPI_STATUS_IGNORE),"Couldn't read from file.");
  }
	
#else /* not MPI */
  /********************* SEQUENTIAL VERSION **************************/
  DEVICE_CONST(FILE *,file);
  DEVICE_CONST(int, rewind);
  int x, y, z, v; /*  Loop indices */
  if (!file) return 1;
  if (rewind) Rewind(file);
  if (s.x0==0 && s.x1==xmax-1 && s.y0==0 && s.y1==ymax-1 && s.z0==0 && s.z1==zmax-1 && v0==0 && v1==vmax-1) {
    size_t size=xmax*ymax*zmax*vmax;
    if (size!=fread(New,sizeof(real),size,file)) error(size);
  } else {
    for (x=s.x0;x<=s.x1;x++) {
      for (y=s.y0;y<=s.y1;y++) {
	for (z=s.z0;z<=s.z1;z++) {
	  for(v=v0;v<=v1;v++) {
	    if (1!=fread(New+ind(x,y,z,v),sizeof(real),1,file)) error(1);
	  }
	}
      }
    }
  }
#endif
}
RUN_TAIL(load)

DESTROY_HEAD(load)
#if MPI
	if(s.runHere){ /*  No file to close if the device isn't active. */
		if(S->file){
			MPI_File_close(&(S->file));
			S->file = NULL;
		}
	}
#else
  SAFE_CLOSE(S->file);
#endif
DESTROY_TAIL(load)

long load_file_size,  field_size;

CREATE_HEAD(load)
{
  DEVICE_IS_RECTANGULAR;

  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 >= v0 );

  ACCEPTI(rewind,0,0,1);
  
#if MPI /*  Prepare for collective reading. */
  
	/*  File */
  MPI_File file;
  MPI_Offset filesize;
  
  /*  Communicator for active dump instances. */
  MPI_Comm comm_load;
  
  if (!deviceCommunicator(dev->s.runHere, &comm_load))
    ABORT("Could not create communicator.\n");
  
  /*  Set up the file	 */
  if (dev->s.runHere) { /*  Active devices only!	 */
    
    /*	Only define types on active instances,
     *	as inactive instances will give invalid
     *	subarray dimensions. */
#include "dump_types.h"
    
    /*  Get filename */
    if (!accepts("file=",&(S->filename[0]),NULL,w)) return(0);
    
    field_size = (v1-v0+1)
      * (dev->s.global_x1-dev->s.global_x0+1)
      * (dev->s.global_y1-dev->s.global_y0+1)
      * (dev->s.global_z1-dev->s.global_z0+1);
    
    /*  Open file and set view */
    MPIDO(MPI_File_open(comm_load, S->filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file),"Couldn't open file");
    MPIDO(MPI_File_set_view(file, 0, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL),"Couldn't set view");
    
    if(MPI_File_seek(file, 0, MPI_SEEK_END) != MPI_SUCCESS){
      ABORT("Process %d could not find end of file %s",mpi_rank, S->filename);
    }
    
    MPI_File_seek(file,0,MPI_SEEK_SET); /*  Equivalent to rewind() */
    
    MPI_File_get_size(file, &filesize);
    if (filesize == 0) MESSAGE("/* WARNING empty file %s */",S->filename);
    if (filesize < field_size) EXPECTED_ERROR("\n file size %ld < field size %ld\n",filesize,field_size);
    if (filesize%field_size!=0) MESSAGE(
					"/* WARNING field size %ld is not divisor of file size %ld */",
					field_size, 
					filesize
					);
    
    /*  Assign file to the parameter structure. */
    S->file = file;
    S->sourceType = sourceType;
  }
  
#else /*  Sequential Version */
  ACCEPTF(file,"rb",NULL);
  
  field_size = (v1-v0+1)
    * (dev->s.x1-dev->s.x0+1)
    * (dev->s.y1-dev->s.y0+1)
    * (dev->s.z1-dev->s.z0+1);
  ASSERT(field_size>0);
  
  if (0!=fseek(S->file,0,SEEK_END)) ABORT("could not find end of file %s",S->filename);
  load_file_size=ftell(S->file);
  MESSAGE("/* size of %s is %ld */",S->filename,load_file_size);
  Rewind(S->file);
  if (load_file_size==0) MESSAGE("/* WARNING empty file %s */",S->filename);
  if (load_file_size < field_size) EXPECTED_ERROR("\n file size %ld < field size %ld\n",load_file_size,field_size);
  if (load_file_size%field_size!=0) MESSAGE(
					    "/* WARNING field size %ld is not divisor of file size %ld */",
					    field_size, 
					    load_file_size
					    );
#endif
}
CREATE_TAIL(load,0)
    
