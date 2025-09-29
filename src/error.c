#include <stdarg.h>
#include <stdio.h>

#include "beatbox.h"
#include "bikt.h"
#include "state.h"
#include "device.h"
#include "mpi_io_choice.h"
#include "error.h"

extern int Mute;                 /* defined in beatbox.c */

/* Generic message procedure, writes both to stdout and log file. */
/* Works from any thread only if expected that only one thread may generate it. */
/* Outputs to stdout and log file separately if in sequential mode, and all to stdout in MPI mode. */
void ANY_MESSAGE(int urgent, char *fmt, ...)
{
  char s[MAXSTRLEN], *p;
  va_list argptr;
  if (mpi_rank==0 || urgent) { /*  Avoid repeated messages on MPI unless urgent */
    va_start(argptr, fmt);
    vsnprintf(s, MAXSTRLEN, fmt, argptr);/* form the message */
    va_end(argptr);
    p=s;
    if (s[0]=='\x01') p++;                    /* \x01 - stdout mutes */
#if MPI
    fputs(p,stdout);
    FFLUSH(stdout);
#else
    if (res) {				    /* to log file - always */
      fputs(p,res);
      FFLUSH(res);
    } else {
      fputs(p,stdout);
      FFLUSH(stdout);
    }
    if (p==s && Mute==0) {		    /* to stdout */
      fputs(p,stdout); 
      FFLUSH(stdout);
    }
#endif
  }
}

void Debug(char *fmt, ...)
{
  va_list argptr;
  if NOT(debug) return;
  va_start(argptr, fmt);
  vfprintf(debug, fmt, argptr);             /* form the message */
  va_end(argptr);
  /* FFLUSH(debug); */ /* debug is about thoroughness, not performance! */
  fflush(debug);
}

int nofflush(void *f) {return 0;}

static char Sprintf_buffer[MAXSTRLEN];
char *Sprintf(char *fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  vsnprintf(Sprintf_buffer, MAXSTRLEN, fmt, argptr);
  va_end(argptr);
  return(&(Sprintf_buffer[0]));
}


