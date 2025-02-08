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

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "k_.h"
#include "state.h"
#include "rhs.h"
#include "qpp.h"

#include "sequence.h"

#if MPI
int getMPIFileMode (const char *fm);
#endif

extern int idev;
extern char *device_name;
extern int Verbose;


static size_t fsize (const char *fid) {
  FILE *f=fopen(fid,"r");
  size_t filesize;
  fseek(f,0,SEEK_END);
  filesize=ftell(f);
  fclose(f);
  return filesize;
}


/****************************************************************************************/
/* This function opens the (new) current file in a sequence to write into it.		*/
/* If the sequence is unnumbered, then it backs up the previous version if any.		*/
/****************************************************************************************/
int thisq (sequence *q) 
{
  /* If a file is already open, it is an error */
  if (q->f) ABORT("file already open, q->name='%s' - this should not happen\n",q->name);

  if (strchr(q->mask,'%')) {
    /* A template is given rather than a specific filename */
    snprintf(q->name,MAXPATH,q->mask,q->n);
    if (0==strcmp(q->name,q->mask))
      EXPECTED_ERROR("bad mask: print %d in format \"%s\" gives \"%s\"",q->n,q->mask,q->name);
#if MPI
    MPIDO(MPI_File_open(q->comm, q->name, q->mpiMode, MPI_INFO_NULL, &(q->f)),"Could not open file.");
#else
    if NOT(q->f=fopen(q->name,q->mode))	EXPECTED_ERROR(
      "Cannot open file %s in mode %s. "
      "Check that the directory exists and that you have the necessary permissions.", 
      q->name, q->mode
    );
#endif
  } else {
    /* Otherwise it is a specific file name; we want to backup the original file if any */
    char bakname[MAXPATH];
    strcpy(q->name,q->mask);
    snprintf(bakname,MAXPATH,"%s.BAK",q->name);
    frename(q->name,bakname);
#if MPI
    MPIDO(MPI_File_open(q->comm, q->name, q->mpiMode, MPI_INFO_NULL, &(q->f)),"Could not open file.");
#else
    if NOT(q->f=fopen(q->name,q->mode)) EXPECTED_ERROR(	
      "Cannot open file %s in mode %s. "
      "Check that the directory exists and that you have the necessary permissions.",
      q->name, q->mode
    );
#endif
  }
  if (q->mute==0) MESSAGE("/* t=%ld device %d(%s) opened %s */\n", 
			  t, idev, device_name, q->name);
  return 1;
}

/****************************************************************************************/
/* This function closes the current file and increments the sequence counter,		*/
/* to prepare for the next.								*/
/* 2018/04/23: Do not leave zero-length files.                                          */
/****************************************************************************************/
int nextq (sequence *q) 
{
  long position;
  /* If a file is already open, need to close it first */
  if (!q->f) ABORT("file not open, q->name='%s' - this should not happen\n",q->name);

#if MPI
  MPI_Barrier(q->comm);
  MPIDO(MPI_File_close(&(q->f)),"could not close file");
#else
  position=ftell(q->f);
  if (0!=fclose(q->f)) {
    MESSAGE("/* Error %d:'%s' closing %s at t=%d, device %d(%s) */\n",
	    errno, strerror(errno),
	    q->name,t,idev,device_name);
    MESSAGE("/* This will be ignored and computations proceed futher. */\n");
  }
  if (position==0) unlink(q->name); /* Do not leave zero-length files */
#endif
  q->f=NULL;
  if (q->mute==0) MESSAGE("/* t=%ld device %d(%s) closed %s */\n",
			    t, idev, device_name, q->name);
  /* Post-processing of the just closed file, unless it was trial opening */
  if (mpi_rank==0 && q->trial==0 && *(q->postproc)!='\0') {
    char cmd[128*1024];
    int rc;
    /* assume the name of the file is mentioned no more than 8 times in the command */
    snprintf(cmd,MAXPATH,q->postproc,q->name,q->name,q->name,q->name,q->name,q->name,q->name,q->name);
    if (q->mute==0) MESSAGE("/* t=%ld device %d(%s) is about to execute command %s */\n", 
			    t, idev, device_name, cmd);
    rc=system(cmd);
    if (q->mute==0 || rc!=0) MESSAGE("/* t=%ld device %d(%s): command '%s' returned %d */\n",
				     t, idev, device_name, cmd, rc);
  }

  (q->n)++;
  return 1;
}




/*************************************************************************************/
/* This function accepts parameters of a file sequence,                              */
/* and initializes that sequence, including a trial opening of the first file in it. */
/* Experimental extension: a postprocessing command is forked (from thread 0 in MPI) */
/*************************************************************************************/
#if MPI
int acceptq (const char *name,const char *mode,const char *deflt,int runHere,sequence *q, char *w)
#else
int acceptq (const char *name,const char *mode,const char *deflt,sequence *q, char *w)
#endif
{
  char *pnum; /* pointer to the initial value of the counter */
  char onlyname[8*1024];
  char postname[8*1024];
#if MPI
  if (!deviceCommunicator(runHere, &(q->comm))) ABORT("Could not create communicator.\n");
  q->mpiMode = getMPIFileMode(mode);
  if (runHere)
#endif
{
  strcpy(q->mode,mode);
  q->f = NULL; /*  Prevents error on first close. We assume q->f non-NULL iff file is open. */

  if NOT(accepts(name,&(q->mask[0]),deflt,w)) return 0; /* Read in the sequence template. */
  if NOT(accepts("postproc=",&(q->postproc[0]),"",w)) return 0; /* Read in the postprocessing command. */
  if NOT(accepti("autonumber=",&(q->autonumber),1,0,1,w)) return 0; /* Read in the autonumbering flag. */
  if NOT(accepti("startfrom=",&(q->startfrom),0,INONE,INONE,w)) return 0; /* Read in the default initial number. */
  if NOT(accepti("ignoreempty=",&(q->ignoreempty),0,0,1,w)) return 0; /* Read in the flag to ignore empty file. */

  /* Empty or 'null' filename means output to nowhere. */
  /* Might be needed if device has other functions apart from output to the file sequence. */
  if (0==strcmp(q->mask,"") || 0==strcmp(q->mask,null)) {
    strcpy(q->name,null);
    q->f=NULL;
    MESSAGE("\x01%s\"%s\"%c",name,NULLFILE,SEPARATORS[0]);
    return 1;
  } 

  /* Asterisk in the beginning of the template means: */
  /* do report opening new files in the sequence  */
  if ((q->mask)[0]=='*') {
    q->mute=0;
    memmove((q->mask),(q->mask)+1,MAXPATH-1);
  } else {
    q->mute=1;
  }


  /*****************************************************************************/
  /* We check that the (first) file can be open, but then immediately close it */
  /* to leave the actual opening and sequencing order to the calling device.   */

  /*  If the template does not contain a format specifier, test the literal filename. */
  if (NULL==(strchr(q->mask,'%'))) {
    strcpy(q->name,q->mask);
#if MPI
    MPIDO(MPI_File_open(q->comm, q->name, q->mpiMode, MPI_INFO_NULL, &(q->f)),"Could not open file.");
#else
    if NOT(q->f=fopen(q->name,q->mode)) EXPECTED_ERROR(
      "Could not open %s in mode %s. "
      "Check that the directory exists and that you have the necessary permissions.",
    q->name,q->mode);
#endif
    /* It was open ok, so now can close it. */
#if MPI
    MPI_Barrier(q->comm);
    MPIDO(MPI_File_close(&(q->f)),"could not close file");
#else
    fclose(q->f);
#endif
    q->f=NULL;
    MESSAGE("\x01%s\"%s\"%c",name,q->name,SEPARATORS[0]);
    return 1;

  } 
  /*  The template contains a format specifier. */
  else {

    /* Check if the format specifier contains an initial value of the counter */
    strtok(q->mask,",");
    pnum=strtok(NULL,",");
    if (pnum==NULL)
      q->n=0;	/* No it does not so we assume zero. */
    else /* Check if it is a number */
      if (1!=sscanf(pnum,"%d",&(q->n))) 
      EXPECTED_ERROR("cannot recognise \"%s\" as starting serial number",pnum);

    /* Start from first non-existent file matching that template. */
    if (q->autonumber) {
      /* Scan through the filenames of the sequence, in order */
      for((q->n)=q->startfrom;;(q->n)++) {
	snprintf(q->name,MAXPATH,q->mask,q->n);
	/* Insure against deadlocks due to weird formatting characters */
	if (0==strcmp(q->name,q->mask)) 
	  EXPECTED_ERROR("Weird mask '%s' produces itself while printing %d through it\n",q->mask,q->n);
	/* We are after the first file of this sequence that does not yet exist, */
	/* i.e. can not be open for reading. */
#if MPI
	if (MPI_File_open(q->comm, q->name, MPI_MODE_RDONLY, MPI_INFO_NULL, &(q->f)) != MPI_SUCCESS) break;
	MPI_File_close(&(q->f));
#else
	if NOT((q->f)=fopen(q->name,"r")) break;
	fclose(q->f);
#endif
      }
      q->f=NULL; /* no file is expected to be open at this point */
      
      if (q->ignoreempty) { /* discard all empty files with biggest numbers */
	while ((q->n) > (q->startfrom)) {
	  (q->n)--;
	  snprintf(q->name,MAXPATH,q->mask,q->n);
	  if (fsize(q->name)>0) {
	    (q->n)++;
	    break;
	  }
	}
      }
    } else {
      q->n=q->startfrom;
    }
    /* Now q->n contains the initial value of the counter; */
    /* try to open the corresponding file name. */
    q->trial=1;
    if NOT(thisq(q)) ABORT("could not open \"%s\"=\"%s\",%d",q->name,q->mask,q->n);
    if NOT(nextq(q)) ABORT("could not close \"%s\"=\"%s\",%d",q->name,q->mask,q->n-1);
    q->n--;    /* Reset the counter back to the initial value. */
    q->trial=0;

    MESSAGE("\x01/* check: %s\"%s\",%d=\"%s\" */",name,q->mask,q->n,q->name);
  }
} /*  if(runHere) */
  return 1;
} /* acceptq */


#if MPI
/*  Returns an MPI file access mode roughly  */
/*  equivalent to the given fopen() mode. */
int getMPIFileMode (const char *fm) {
  int MPImode;
	
  if (		0==strcmp(fm,"r") || 0==strcmp(fm,"rb")){
    MPImode = MPI_MODE_RDONLY;
    /* MESSAGE("Mode = MPI_MODE_RDONLY");fflush(stdout); // DEBUG */
  } else if (	0==strcmp(fm,"w") || 0==strcmp(fm,"wb") || 
		0==strcmp(fm,"a") || 0==strcmp(fm,"ab") ||
		0==strcmp(fm,"a+")|| 0==strcmp(fm,"ab+")||
		0==strcmp(fm,"a+b")){
    MPImode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    /* MESSAGE("Mode = MPI_MODE_RDONLY | MPI_MODE_CREATE");fflush(stdout); // DEBUG */
  } else if (	0==strcmp(fm,"r+")|| 0==strcmp(fm,"rb+") || 
		0==strcmp(fm,"r+b")){
    MPImode = MPI_MODE_RDWR;
    /* MESSAGE("Mode = MPI_MODE_RDWR");fflush(stdout); // DEBUG */
  } else if (	0==strcmp(fm,"w+")|| 0==strcmp(fm,"wt") ||
		0==strcmp(fm,"w+b")){
    MPImode = MPI_MODE_RDWR | MPI_MODE_CREATE;
    /* MESSAGE("Mode = MPI_MODE_RDWR | MPI_MODE_CREATE");fflush(stdout); // DEBUG */
  } else{
    EXPECTED_ERROR("Could not create equivalent MPI file access mode given \"%s\"\n", fm);
  }
  
  return MPImode;
}
#endif
