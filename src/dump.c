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
 * Dump contents of a 4D subset to a binary file
 */

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

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
  MPI_Datatype sourceType;	/*  Datatype describing subarray of New to be written. */
#else
  FILE *file;			/*  File handle for sequential I/O. */
#endif
  int append;
} STR;

/*  ----------------------------------------------------------------// */
RUN_HEAD(dump)
{
  DEVICE_CONST(int,v0);
  DEVICE_CONST(int,v1);
  /********************** MPI VERSION ********************************/
#if MPI
  DEVICE_CONST(int, append);
  DEVICE_CONST(MPI_File, file);
  DEVICE_CONST(MPI_Datatype, sourceType);

  /*  Rewind file if not appending */
  if(!append) MPIDO(MPI_File_seek(file,0,MPI_SEEK_SET),"could not rewind file");
  
  /* This writes to a part of file allocated to this process, according to sourceType */
  MPIDO(MPI_FILE_WRITE(file, New, 1, sourceType, MPI_STATUS_IGNORE),"Couldn't write to file.");
  
#else /* not MPI */
  /********************* SEQUENTIAL VERSION **************************/
  DEVICE_CONST(FILE *,file)
  DEVICE_CONST(int, append)
  int x, y, z, v; /*  Loop indices */
  if (!file) return 1;
  if (!append) rewind(file);
  if (s.x0==0 && s.x1==xmax-1 && s.y0==0 && s.y1==ymax-1 && s.z0==0 && s.z1==zmax-1 && v0==0 && v1==vmax-1) {
    fwrite(New,sizeof(real),xmax*ymax*zmax*vmax,file);
  } else {
    for(x=s.x0;x<=s.x1;x++){
      for(y=s.y0;y<=s.y1;y++){
	for(z=s.z0;z<=s.z1;z++){
	  for(v=v0;v<=v1;v++){
	    fwrite(New+ind(x,y,z,v),sizeof(real),1,file);
	  }
	}
      }
    }
  }
  FFLUSH(file);
  /********************************************************************/
#endif
  MESSAGE("dumped [%d:%d]x[%d:%d]x[%d:%d]x[%d:%d] to %s at t=%ld\n",
	  s.x0, s.x1, s.y0, s.y1, s.z0, s.z1, v0, v1, S->filename,t);
}
RUN_TAIL(dump)

/* ---------------------------------------------------------------- */

DESTROY_HEAD(dump)
#if MPI
 if (s.runHere) { /*  No file to close if the device isn't active. */
   MPI_File_close(&(S->file));
 }
#else
  SAFE_CLOSE(S->file);
#endif
DESTROY_TAIL(dump)

#if MPI
#define sizereport(msg) {MPIDO(MPI_File_get_size(file,&currentsize),"Could not get size at " msg);printf("Process %d %s: filesize=%d\n",mpi_rank,msg,(int)currentsize);}
static MPI_Offset currentsize;
#endif


/*  ----------------------------------------------------------------// */

CREATE_HEAD(dump)
{
  DEVICE_IS_RECTANGULAR;
  
  ACCEPTI(v0,INONE,0,(int)vmax-1);
  ACCEPTI(v1,INONE,0,(int)vmax-1);
  ASSERT( v1 >= v0 );

  ACCEPTI(append,0,0,1);

#if MPI /*  Prepare for collective writing. */
  if (!accepts("file=",&(S->filename[0]),NULL,w)) return(0);

  /* Number of data elements to be written */
  int xsize = (dev->s.global_x1 - dev->s.global_x0) + 1;
  int ysize = (dev->s.global_y1 - dev->s.global_y0) + 1;
  int zsize = (dev->s.global_z1 - dev->s.global_z0) + 1;
  int vsize = (v1 - v0) + 1;
  int boxsize = xsize*ysize*zsize*vsize;
	
  MPI_Comm comm_dump;	/*  Communicator for active dump instances. */
  int root;		/* One processor chosen to pre-fill */
  
  MPI_File file;	/* The dump file written to */
  int mode;		/* mode: write or append */
  MPI_Aint lb;          /* Type's lower bound */
  MPI_Offset displacement; /* how much to skip from the beginning of file */
  MPI_Aint typesize;	/* size of data to be written according to defined MPI type */

  /* We need one process guaranteed to be engaged in this device */
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &comm_dump, &root))
    ABORT("Could not create communicator.\n");
  
  if (dev->s.runHere) { /*  Active devices only!	 */
    
    /* Only define types on active instances,  */
    /* as inactive instances will give invalid */
    /* subarray dimensions.                    */
    #include "dump_types.h"
    
    /*  Set access mode */
    if (S->append) {
      mode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND;
    } else {
      mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    }

    /* It would appear MPI cannot overwrite an existing file, so take precaution */
    if (mpi_rank==root) {
      if (fexist(S->filename)) {
	if (Verbose) URGENT_MESSAGE("File %s exists; deleting it\n",S->filename);
	if (0!=unlink(S->filename)) {
	  URGENT_MESSAGE("Could not delete %s: %s\n",S->filename,strerror(errno));
	}
      }
    }
    MPIDO(MPI_Barrier(comm_dump),"Could not stop at barrier before opening file");
    
    /*  Open and preallocate file */
    MPIDO(MPI_File_open(comm_dump, S->filename, mode, MPI_INFO_NULL, &file),"Couldn't open file");
    /*  If requested by user, append to EOF */
    if (S->append) {
      MPIDO(MPI_File_get_size(file, &displacement),"Couldn't get existing file size");
    } else {
      displacement = (MPI_Aint)0;
    }
    MPIDO(MPI_File_preallocate(file,displacement + boxsize*sizeof(real)),"Couldn't preallocate file");
 
    /* Prefill the file if there are empty subdomains.    */
    /* This is done by one process only, which needs then */
    /* to rewind back to the same position.               */
    if (num_empty_subdoms>0 && mpi_rank==root) {
      real *prefill;	 /* prefill data buffer */
      int i;		 /* prefill buffer counter */
      MPI_Offset offset; /* file pointer position before prefill */
      CALLOC(prefill,boxsize,sizeof(real));
      for (i=0;i<boxsize;i++) prefill[i]=0.0; /* or whatever value represents void - TBC */
      MPIDO(MPI_File_get_position(file,&offset ),"Could not remember position before prefill"); 
      MPIDO(MPI_File_write(file,prefill,boxsize,MPI_DOUBLE,MPI_STATUS_IGNORE),"Couldn't prefill the file");
      MPIDO(MPI_File_seek(file,offset,MPI_SEEK_SET),"could not rewind to before prefill");
      FREE(prefill);
    } 

    /* Set view for collective write from now on */
    MPIDO(MPI_File_set_view(file,displacement,MPI_DOUBLE,filetype,"native",MPI_INFO_NULL),"Couldn't set start of file");

    /* Check there was no misunderstanding about amount of data to be written */
    MPIDO(MPI_Type_get_extent(filetype,&lb,&typesize),"Couldn't get filetype extent");	
    ASSERT(typesize == boxsize*sizeof(real));
   
    /*  Assign the file to the parameter structure */
    S->file = file;
    S->sourceType = sourceType;
  } /* if runHere */

#else /*  Sequential version */
  ACCEPTF(file,S->append?"ab":"wb",NULL);
#endif
}
CREATE_TAIL(dump,0)
