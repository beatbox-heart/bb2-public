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

/* Control point device. Saves and restores global variables and the
 * grid.  Optionally, keeps control point in a directory, with global
 * data and grid subsets for individual processes in separate files.
 * This makes it impossible to restore to a different partition scheme
 * than it was saved from, but has the advantage for very big grids in
 * that it avoids limitations of MPI parallel writing.
 */

#include <dirent.h>	
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <stdarg.h>
#include "system.h"
#include "beatbox.h"
#include "k_.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "mpi_io_choice.h"

typedef struct {
  char file[MAXPATH];		/* this will be a directory name if bythread!=0 */
  int bythread;			/* save grid subsets separately */
  int enrich;			/* k-variables stored in ctl point are created if not existed in this run */
  long last;			/* when this device was called last time */
  char debugname[MAXPATH];	/* the file to which we print the debugging info, */
  FILE *debug;			/*    may or may not coincide with common debug file */
#if MPI
  MPI_Datatype k_var_Type;	/* MPI        */
  MPI_Datatype filetype;	/*   parallel */
  MPI_Datatype sourceType;	/*   i/o      */
  MPI_Aint     filetype_extent;	/*   stuff    */
#endif
} STR;

#define maxstr 256
typedef struct {
  char          nm[maxname];
  unsigned char tp;
  unsigned char val[maxstr];
} k_var;

static k_var empty_var={"", 0, ""};

/* Return values of some static functions */
#define FAILURE 0
#define SUCCESS 1

#define NOT(a) (0==(a))

/* Message to a local or to the common debug file */
/* #define DEBUG(...) if (debug) {fprintf(debug, __VA_ARGS__); FFLUSH(debug);} */

static void fputext (const char *oldfid, char *newfid, const char *newext)
{
  char *p;
  p = strrchr(oldfid,'.');
  if (p) {
    strncpy(newfid,oldfid,p-oldfid);
    p = newfid + (p-oldfid);
    *(p++)='.';
    strcpy(p,newext);
  } else {
    snprintf(newfid,MAXPATH,"%s.%s",oldfid,newext);
  }
}

static int fexist (const char *fid, FILE *debug)
{
  FILE *f=fopen(fid,"r");
  if (f) fclose(f);
  return (f!=NULL);
}

/* Zero return means it does not exist */
static int dexist (const char *dirname, FILE *debug)
{
  DIR *d=opendir(dirname);
  if (d) closedir(d);
  else DEBUG("Could not open directory %s: %s\n",dirname,strerror(errno));
  return (d!=NULL);
}

/* Zero return means failure */
static int deldir (const char *dirname)
{
  char cmd[MAXPATH];
  snprintf(cmd,MAXPATH,"rm -rf %s",dirname);
  return (0==system(cmd)) ? SUCCESS : FAILURE;
}

/* Zero return means failure */
static int frename (const char *oldfid, const char *newfid,FILE *debug)
{
  int retcode;
  if (fexist(newfid,debug)) {
    char bakname[MAXPATH];
    fputext(newfid, bakname, "BAK");
    if (fexist(bakname,debug)) {
      if (0!=remove(bakname)) {
	DEBUG("Could not delete %s: %s\n",bakname,strerror(errno));
	return FAILURE;
      } /* if not remove */
    } /* if fexist backname */
    if (0!=rename(newfid,bakname)) {
      DEBUG("Could not rename %s->%s\n",newfid,bakname);
      return FAILURE;
    }
  }
  if (0!=rename(oldfid,newfid)) {
    DEBUG("Could not rename %s->%s\n",oldfid,newfid);
    return FAILURE;
  }
  return SUCCESS;
}

/* Zero return means failure */
static int drename(const char *olddirname, const char *newdirname,FILE *debug)
{
  int retcode;
  if (dexist(newdirname,debug)) {
    char bakname[MAXPATH];
    fputext(newdirname,bakname,"BAK");
    if (FAILURE==deldir(bakname)) {DEBUG("Could not delete %s\n",bakname);return FAILURE;}
    if (0!=rename(newdirname,bakname)) {DEBUG("Could not rename %s->%s\n",newdirname,bakname);return FAILURE;}
  }
  if (0!=rename(olddirname,newdirname)) {DEBUG("Could not rename %s->%s\n",olddirname,newdirname);return FAILURE;}
  return SUCCESS;
}

/* Macros for frequent things done differently within MPI and without */
#if MPI
  #define FREAD(addr, size, type, mpitype, ...)		  			\
  if (bythread) {							  	\
    if (size!=fread(addr,sizeof(type),size,globalfile)) ABORT(__VA_ARGS__);	\
  } else {									\
    MPIDO(MPI_File_read(f, addr, size, mpitype, MPI_STATUS_IGNORE),__VA_ARGS__);\
  }
  #define FWRITE(addr, size, type, mpitype, ...)				\
  if (bythread) {								\
   if (size != fwrite(addr,sizeof(type),size,globalfile)) ABORT(__VA_ARGS__);	\
  } else {									\
   MPIDO(MPI_File_write(f, addr, size, mpitype, MPI_STATUS_IGNORE),__VA_ARGS__);\
  }
  #define FCLOSE(...) MPIDO(MPI_File_close(&f),__VA_ARGS__); f = MPI_FILE_NULL;
#else /* if MPI */
  #define FREAD(addr, size, type, mpitype, ...)					\
  if (bythread) {							  	\
    if (size!=fread(addr,sizeof(type),size,globalfile)) ABORT(__VA_ARGS__);	\
  } else {									\
    if (size!=fread(addr,sizeof(type),size,f)) ABORT(__VA_ARGS__);		\
  }
  #define FWRITE(addr, size, type, mpitype, ...)				\
  if (bythread) {								\
   if (size != fwrite(addr,sizeof(type),size,globalfile)) ABORT(__VA_ARGS__);	\
  } else {									\
    if (size != fwrite(addr,sizeof(type),size,f)) ABORT(__VA_ARGS__);		\
  }
  #define FCLOSE(...) if (0!=fclose(f)) ABORT(__VA_ARGS__); f=NULL;				
#endif /* if MPI else */

RUN_HEAD(ctlpoint)
{
  DEVICE_ARRAY(char,file);
  DEVICE_CONST(int,bythread);
  DEVICE_CONST(int,enrich);
  DEVICE_VAR(long,last);
  DEVICE_CONST(FILE *,debug);

  #if MPI
    DEVICE_CONST(MPI_Datatype, k_var_Type);
    DEVICE_CONST(MPI_Datatype, filetype);
    DEVICE_CONST(MPI_Datatype, sourceType);
    DEVICE_CONST(MPI_Aint, filetype_extent);
    MPI_File f;
    MPI_Offset offset;
  #else
    FILE *f;
  #endif
    
  FILE *globalfile;			/* the file for global variables */
  char globalfilename[MAXPATH];		/* .., the name of */
  FILE *subsetfile;			/* the file for this processes's grid subset */
  char subsetfilename[MAXPATH];		/* .., the name of */
  size_t subsetnum=xlen*ylen*zlen*vmax;	/* .., the number of values in */

  /* New values of principal global vars for comparison */
  char V_[32];
  INT xm_,ym_,zm_,vm_,mpinx_,mpiny_,mpinz_;
  
  /* Buffers needed since prt() returns a static address */
  #define PRTBUFLEN 80
  char oldbuf[PRTBUFLEN], newbuf[PRTBUFLEN];

  /* Sundry */
  char *p;
  int i;
  k_var var;
  
  DEBUG("\n#%d ctlpoint called at %ld\n",mpi_rank,t);
    
  if (*last==LNONE && (bythread?dexist(file,debug):fexist(file,debug))) {
    /**********************************************/
    /*  Read from previously saved control point  */
    DEBUG("Reading:");
    if (bythread) {
      snprintf(globalfilename,MAXPATH,"%s/global",file);
      if (NULL==(globalfile=fopen(globalfilename,"r")))
	EXPECTED_ERROR("Directory %s exists but file %s not found or not readable\n",
		       file,globalfilename);
      snprintf(subsetfilename,MAXPATH,"%s/%ld_%ld_%ld",file,(long)mpi_ix,(long)mpi_iy,(long)mpi_iz);
      if (NULL==(subsetfile=fopen(subsetfilename,"r")))
	ABORT("Directory %s exists but file %s not found or not readable\n",
	      file,subsetfilename);
    } else { /* if bythread */
      #if MPI
        MPIDO(MPI_File_open(ALL_ACTIVE_PROCS, file, MPI_MODE_RDONLY, MPI_INFO_NULL, &f),
  	    "File %s exists but could not be opened for reading.", file);
        MPIDO(MPI_File_set_view(f, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL),"Could not set view.");
      #else
        if (NULL==(f=fopen(file,"rb")))
          ABORT("File %s exists but could not be opened for reading.", file);
      #endif
    } /* if bythread else */
    
    /*  Read first block of global data. */
#define _(name, new, old, size, type, mpitype, ktype, severity)			\
      FREAD(new, size, type, mpitype, "Error reading %s from %s.",name,file);	\
      if (new!=old && bcmp(new,old,sizeof(type)*size)) {			\
	memcpy(newbuf,prt(new,ktype),PRTBUFLEN);				\
	memcpy(oldbuf,prt(old,ktype),PRTBUFLEN);				\
	if (severity || bythread) {						\
	  ABORT("#%d: wrong %s read from %s: %s instead of %s. This is a critical error.\n", \
		mpi_rank,name,file,newbuf,oldbuf);				\
	} else {								\
	  URGENT_MESSAGE("#%d: wrong %s read from %s: %s instead of %s. Proceed in hope...\n", \
			 mpi_rank,name,file,newbuf,oldbuf);			\
	}									\
      } else									\
	DEBUG("\t%s=%s,",name,prt(new,ktype));
    #include "ctlpoint.h"
    #undef _

    /*  Read the big thing : this is very different between MPI and sequential */
    if (bythread) {
      if (subsetnum!=fread(New,sizeof(real),subsetnum,subsetfile))
        ABORT("Error reading the grid from %s: %s.",subsetfilename,strerror(errno));
      DEBUG("\tgrid[%ldx%ldx%ldx%ld],",(long)xlen,(long)ylen,(long)zlen,vm_);
    } else { /* if bythread */
    #if MPI
      MPIDO(MPI_File_get_position(f, &offset),"Could not get position of the file.");
      MPIDO(MPI_File_set_view(f, offset, sourceType, filetype, "native", MPI_INFO_NULL),"Could not set view.");
      MPIDO(MPI_File_read(f, New, 1, sourceType, MPI_STATUS_IGNORE),"Could not read grid from %s.",file);
      haloSwap();
      DEBUG("\t#%d subgrid[%ldx%ldx%ldx%ld],",mpi_rank,xlen,ylen,zlen,vm_);
      offset += filetype_extent;
      MPIDO(MPI_File_set_view(f, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL),"Could not set view.");
    #else
      if (xmax*ymax*zmax*vm_!=fread(New,sizeof(real),xmax*ymax*zmax*vm_,f))
        ABORT("Error reading the grid from %s: %s.",file,strerror(errno));
      DEBUG("\tgrid[%ldx%ldx%ldx%ld],",xmax,ymax,zmax,vm_);
    #endif
    } /* if bythread else */
    
    /* Read the k-table */
    for (;;) {
      FREAD(&var, 1, k_var, k_var_Type,"Error reading a k-var from %s.", file);
      if NOT(var.nm[0]) break;                /* eof mark */
      if NOT(i=tb_find(deftb,var.nm)) {
	  if (!enrich) continue;
	  if NOT(i=tb_insert_abstract(
				      deftb,var.nm,var.tp,(p_vd)Calloc(1,sizetable[var.tp]),0,f_vb|f_rs
				      )) EXPECTED_ERROR("Cannot define k-var %s.",var.nm);
	} /*  if not defined. */
      if (tb_type(deftb,i)!=var.tp) EXPECTED_ERROR("Variable %s from %s is of wrong type.",var.nm,file);
      if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("%s is not a variable.",var.nm);
      memcpy(tb_addr(deftb,i),&(var.val),sizetable[var.tp]);
      DEBUG("\t%s=%s",var.nm,prt(&(var.val),var.tp));
    } /*  for (;;) */

    /* Close the file(s) */
    if (bythread) {
      if (0!=fclose(globalfile))
	URGENT_MESSAGE("Could not close %s after reading: %s.\n",
		       globalfile,strerror(errno));
      globalfile=NULL;
      if (0!=fclose(subsetfile))
	URGENT_MESSAGE("Could not close %s after reading: %s.\n",
		       subsetfile,strerror(errno));
      subsetfile=NULL;
    } else { /* if bythread */
      FCLOSE("Could not close %s after reading: %s.\n",file,strerror(errno));
    } /* if bythread else */
    
    DEBUG("\treading completed\n");
    MESSAGE("Control point read from '%s' at t=%ld.\n",file,t);
    
  } else { /*  if (last or fexist(file)) */
    /***************************************************************/
    /* this device read or saved before, now save first time/again */
    DEBUG("Writing:");
    vm_=vmax;

    /* Temp filename as modified target filename */
    strcpy(buf,file);
    for (p=buf+strlen(buf)-1;p>=buf&&*p=='~';p--);
    if (p<buf) EXPECTED_ERROR("Could not generate tempname for %s",file);
    *p='~';

    if (bythread) {
      if (mpi_rank==0) {
	/* Erase the file or directory with such a name, if exists, and make a new one */
        if (fexist(buf,debug) || dexist(buf,debug)) {
	  if (FAILURE==deldir(buf)) {
	    URGENT_MESSAGE("#%d: Could not delete %s: %s\n",mpi_rank,buf,strerror(errno));
	  } 
	}
        if (0!=mkdir(buf,0700)) 
	  ABORT("#%d: Couldn't create directory %s: %s\n",mpi_rank,buf,strerror(errno));
	/* And open the file for writing the global data */
	snprintf(globalfilename,MAXPATH,"%s/global",buf);
	if (NULL==(globalfile=fopen(globalfilename,"w")))
	  EXPECTED_ERROR("Couldn't open %s for writing: %s\n",
			 globalfilename,strerror(errno));
      } /* if mpi_rank==0 */
    } else {
      #if MPI
        /* Apparently MPI cannot overwrite an existing file, so take precaution. */
        /* Here we do not expect the new file to be a directory - is that wise?  */
        if (mpi_rank==0) {
          if (fexist(buf,debug)) {
	    if (Verbose) MESSAGE("\nFile %s exists; deleting it.\n",buf);
	    if (0!=unlink(buf)) MESSAGE("Could not delete %s: %s.\n",buf,strerror(errno));
          }
        }
        MPIDO(MPI_Barrier(ALL_ACTIVE_PROCS),"Could not stop at barrier before opening file.");
        MPIDO(MPI_File_open(ALL_ACTIVE_PROCS, buf, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f),
    	  "Could not open %s for writing",buf);
        /*  Set view with no displacement (top of file). */
        MPIDO( MPI_File_set_view(f, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL),"Could not set view.");
      #else
	/* Sequential access should just overwrite any existing file            */
        /* Here we do not expect the new file to be a directory - is that wise? */
        if NOT(f=fopen(buf,"wb")) EXPECTED_ERROR("Could not open %s for writing.",buf);
      #endif
    }
      
    /*  Write first block of global data, only by process 0 */
    if (mpi_rank==0) {
#define _(name, new, old, size, type, mpitype, ktype, severity)			\
	FWRITE(old, size, type, mpitype,"Error writing %s to %s.",name,buf);	\
        DEBUG("\t%s=%s,",name,prt(new,ktype));
      #include "ctlpoint.h"
      #undef _
    }

    /* Save the big thing : this is very different between MPI and sequential */
    if (bythread) {
      /* Make sure the directory was created before opening files in it */
      #if MPI
      MPIDO(MPI_Barrier(ALL_ACTIVE_PROCS),"Could not stop at barrier before opening subset files");
      #endif
      /* Every process writes its subset of the grid to a separate file */
      snprintf(subsetfilename,MAXPATH,"%s/%ld_%ld_%ld",buf,mpi_ix,mpi_iy,mpi_iz);
      if NOT(subsetfile=fopen(subsetfilename,"w"))
	      ABORT("#%d could not open %s for writing: %s\n",
		    mpi_rank,subsetfilename,strerror(errno));
      if (subsetnum!=fwrite(New,sizeof(real),subsetnum,subsetfile))
	ABORT("#%d could not write subset of the grid to %s at t=%ld: %s\n",
	      mpi_rank,subsetfilename,t,strerror(errno));
    } else { /* if bythread */
      #if MPI
        /* Every process writes its subset of the grid to its part of the common file */
        MPIDO(MPI_File_get_position(f, &offset),"Could not get file position.");
        MPIDO(MPI_Bcast(&offset, 1, MPI_INT, 0, ALL_ACTIVE_PROCS),"Could not broadcast offset.");
        MPIDO(MPI_File_set_view(f,offset,sourceType,filetype,"native",MPI_INFO_NULL),"Could not set view.");
        MPIDO(MPI_File_write_all(f, New, 1, sourceType, MPI_STATUS_IGNORE),"Error writing subgrid to %s.",buf);
        DEBUG("\t#%d subgrid[%ldx%ldx%ldx%ld],",mpi_rank,xlen,ylen,zlen,vmax);
        offset += filetype_extent;
        MPIDO(MPI_File_set_view(f,offset,MPI_CHAR,MPI_CHAR,"native",MPI_INFO_NULL),"Could not set view.");
      #else
        if (xmax*ymax*zmax*vmax!=fwrite(New,sizeof(real),xmax*ymax*zmax*vmax,f))
          ABORT("Error reading the grid from %s\n",file);
        DEBUG("\t#%d grid[%ldx%ldx%ldx%ld],",mpi_rank,xmax,ymax,zmax,vmax);
      #endif
    } /* if bythread else */
    
    /*  Save the k-table, process 0 alone again */
    if (mpi_rank == 0) {
      for (i=1;i<=maxtab;i++) {
	if(!*tb_name(deftb,i)) continue;
	if((tb_flag(deftb,i)&f_vb)==0) continue;
	memcpy(var.nm,tb_name(deftb,i),maxname);
	var.tp=tb_type(deftb,i);
	memcpy(&(var.val),tb_addr(deftb,i),sizetable[var.tp]);
	FWRITE(&var, 1, k_var, k_var_Type, "Error writing %s to %s.",var.nm,buf);
	DEBUG("\t%s=%s,",var.nm,prt(var.val,var.tp));
      } /*  for i-->maxtab */
      FWRITE(&empty_var, 1, k_var, k_var_Type,"error writing empty k-var as eof to %s\n",buf);
    } /*  if(mpi_rank == 0) */

    /* Close the file(s) */
    if (bythread) {
      if (mpi_rank==0) {
	if (0!=fclose(globalfile))
	  URGENT_MESSAGE("Could not close %s after writing: %s.\n",
			 globalfilename, strerror(errno));
	globalfile=NULL;
      }
      if (0!=fclose(subsetfile))
	URGENT_MESSAGE("Could not close %s after writing: %s.\n",
		       subsetfilename, strerror(errno));
      subsetfile=NULL;								
    } else {
      FCLOSE("Could not close %s after writing: %s.\n",buf, strerror(errno));
    }

    /* Make sure all files were closed before renaming the directory */
    #if MPI
    MPIDO(MPI_Barrier(ALL_ACTIVE_PROCS),"Could not stop at barrier after closing subset files");
    #endif
    if (mpi_rank==0) {
      if (SUCCESS!=(bythread?drename(buf,file,debug):frename(buf,file,debug)))
	URGENT_MESSAGE("#%d: Could not rename %s to %s: %s\n",mpi_rank,buf,file,strerror(errno));
    }
    DEBUG("\twriting completed\n");
    MESSAGE("Control point written to '%s' at t=%ld\n",file,t);
    
  } /*  if (last or fexist(file)) else */
  *last = t;
}
RUN_TAIL(ctlpoint)

/********************/
DESTROY_HEAD(ctlpoint) {
  SAFE_CLOSE(S->debug);
} DESTROY_TAIL(ctlpoint)

/*******************/
CREATE_HEAD(ctlpoint)
{
  DEVICE_IS_SPACELESS;
        
  ACCEPTS(file,NULL);
  ACCEPTI(bythread,0,0,1);
  ACCEPTI(enrich,0,0,1);
  ACCEPTF(debug,"wt","");
  S->last=LNONE;
  int v0 = 0;
  int v1 = vmax-1;
  
#if MPI
  Space *s=&(dev->s);
  
  /* Define the k-variables type in MPI parlance */
  MPI_Datatype k_var_Type;
  k_var var;
  int lengths[3] = {maxname, 1, maxstr};
  MPI_Aint *displacements = calloc(3, sizeof(MPI_Aint));
  MPI_Datatype old_datatypes[3] = {MPI_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};
  int i;
  MPI_Get_address(var.nm,	&displacements[0]);
  MPI_Get_address(&var.tp,	&displacements[1]);
  MPI_Get_address(var.val,	&displacements[2]);
  for (i = 2; i>=0; i--) displacements[i] -= displacements[0];
  MPI_Type_create_struct(3, lengths, displacements, old_datatypes, &k_var_Type);
  MPI_Type_commit(&k_var_Type);

  /*  Global space must be whole medium (inc absolute boundaries.) */
  s->runHere = 1; /*  Must run on all processes. */
  s->global_x0 = 0;
  s->global_x1 = xmax-1;
  s->global_y0 = 0;
  s->global_y1 = ymax-1;
  s->global_z0 = 0;
  s->global_z1 = zmax-1;
  s->x0 = (mpi_ix == 0)        ? s->global_x0 : local_xmin;
  s->x1 = (mpi_ix == mpi_nx-1) ? s->global_x1 : local_xmax-1;
  s->y0 = (mpi_iy == 0)        ? s->global_y0 : local_ymin;
  s->y1 = (mpi_iy == mpi_ny-1) ? s->global_y1 : local_ymax-1;
  s->z0 = (mpi_iz == 0)        ? s->global_z0 : local_zmin;
  s->z1 = (mpi_iz == mpi_nz-1) ? s->global_z1 : local_zmax-1;

  /*  Define types for New. */
  #include "dump_types.h"
  MPI_Aint lb, filetype_extent;
  MPIDO(MPI_Type_get_extent(filetype,&lb,&filetype_extent),"Could not get filetype extent");

  S->k_var_Type = k_var_Type;
  S->filetype = filetype;
  S->sourceType = sourceType;
  S->filetype_extent = filetype_extent;
#endif
}
CREATE_TAIL(ctlpoint,0)

