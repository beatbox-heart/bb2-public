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

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifdef __APPLE__
#ifndef _DARWIN_C_SOURCE
#define _DARWIN_C_SOURCE
#endif
#endif

#include <signal.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include "pipe.h"

static void remove_fifo(PIPE *p) {
  if (p->n && unlink(p->n) == -1 && errno != ENOENT)
    perror("pipeto could not remove fifo");
  if (p->dir && rmdir(p->dir) == -1 && errno != ENOENT)
    perror("pipeto could not remove temporary directory");
}

static void free_pipe(PIPE *p) {
  if (!p) return;
  free(p->n);
  free(p->dir);
  free(p);
}

static void discard_pipe(PIPE *p) {
  remove_fifo(p);
  free_pipe(p);
}

static int wait_for_child(pid_t child, int *status) {
  pid_t result;
  do {
    result = waitpid(child, status, 0);
  } while (result == -1 && errno == EINTR);
  if (result == -1) {
    perror("pipe waitpid");
    return -1;
  }
  return 0;
}

static int stop_child(pid_t child) {
  int status;
  if (kill(child,SIGTERM)==-1 && errno!=ESRCH) {
    perror("pipeto could not stop child");
    return -1;
  }
  return wait_for_child(child,&status);
}

PIPE *pipeto(char *cmd) {
  PIPE *p=calloc(1,sizeof(PIPE));
  const char *tmpdir=getenv("TMPDIR");
  const char *base=(tmpdir && *tmpdir)?tmpdir:"/tmp";
  size_t dirlen=strlen(base)+sizeof("/beatbox-XXXXXX");
  size_t fifolen;
  int fd, flags, pathlen;
  pid_t pid;

  if (!p) {
    perror("pipeto could not allocate state");
    return NULL;
  }
  if (!cmd) {
    errno=EINVAL;
    perror("pipeto received a null command");
    free_pipe(p);
    return NULL;
  }
  p->dir=malloc(dirlen);
  if (!p->dir) {
    perror("pipeto could not allocate temporary directory name");
    free_pipe(p);
    return NULL;
  }
  pathlen=snprintf(p->dir,dirlen,"%s/beatbox-XXXXXX",base);
  if (pathlen<0 || (size_t)pathlen>=dirlen) {
    errno=ENAMETOOLONG;
    perror("pipeto temporary directory name");
    free_pipe(p);
    return NULL;
  }
  if (!mkdtemp(p->dir)) {
    perror("pipeto could not make temporary directory");
    free_pipe(p);
    return NULL;
  }
  fifolen=strlen(p->dir)+sizeof("/output");
  p->n=malloc(fifolen);
  if (!p->n) {
    perror("pipeto could not allocate fifo name");
    discard_pipe(p);
    return NULL;
  }
  pathlen=snprintf(p->n,fifolen,"%s/output",p->dir);
  if (pathlen<0 || (size_t)pathlen>=fifolen) {
    errno=ENAMETOOLONG;
    perror("pipeto fifo name");
    discard_pipe(p);
    return NULL;
  }
  if (mkfifo(p->n,0600)==-1) {
    perror("pipeto could not make a fifo");
    discard_pipe(p);
    return NULL;
  }
  switch (pid=fork()) {
  case -1: 
    perror("pipeto could not fork"); 
    discard_pipe(p);
    return NULL;
  case 0: {
    int input=open(p->n,O_RDONLY);
    if (input==-1) {
      perror("pipeto child could not open fifo");
      _exit(126);
    }
    if (input!=STDIN_FILENO) {
      if (dup2(input,STDIN_FILENO)==-1) {
        perror("pipeto child could not connect fifo to stdin");
        close(input);
        _exit(126);
      }
      close(input);
    }
    execl("/bin/sh","sh","-c",cmd,(char *)NULL);
    perror("pipeto child could not execute command");
    _exit(127);
  }
  default: 
    p->child=pid;
    for (;;) {
      fd=open(p->n,O_WRONLY|O_NONBLOCK);
      if (fd!=-1) break;
      if (errno!=ENXIO && errno!=EINTR) {
        perror("pipeto could not open fifo for writing");
        stop_child(pid);
        discard_pipe(p);
        return NULL;
      }
      {
        int status;
        pid_t child_result=waitpid(pid,&status,WNOHANG);
        if (child_result==pid) {
          fprintf(stderr,"pipeto child exited before opening fifo (status %08x)\n",status);
          discard_pipe(p);
          return NULL;
        }
        if (child_result==-1 && errno!=EINTR) {
          perror("pipeto could not check child status");
          stop_child(pid);
          discard_pipe(p);
          return NULL;
        }
      }
      {
        struct timespec delay={0,10000000};
        while (nanosleep(&delay,&delay)==-1) {
          if (errno==EINTR) continue;
          perror("pipeto could not wait for fifo reader");
          stop_child(pid);
          discard_pipe(p);
          return NULL;
        }
      }
    }
    flags=fcntl(fd,F_GETFL);
    if (flags==-1 || fcntl(fd,F_SETFL,flags & ~O_NONBLOCK)==-1) {
      perror("pipeto could not configure fifo writer");
      close(fd);
      stop_child(pid);
      discard_pipe(p);
      return NULL;
    }
    p->f=fdopen(fd,"w");
    if (!p->f) {
      perror("pipeto could not open fifo stream");
      close(fd);
      stop_child(pid);
      discard_pipe(p);
      return NULL;
    }
    return p;
 }
}

int pipeclose(PIPE *p) {
  int fcloseret=0;
  int status=0;
  if (!p) {
    errno=EINVAL;
    perror("pipeclose received a null pipe");
    return -1;
  }
  if (0!=(fcloseret=fclose(p->f))) perror("pipeclose fclose fifo input");
  if (wait_for_child(p->child,&status)==-1) fcloseret=-1;
  else if (!WIFEXITED(status)) fprintf(stderr,"pipeclose child status %08x\n",status);
  if (p->n && unlink(p->n)==-1 && errno!=ENOENT) {
    perror("pipeclose could not remove fifo");
    fcloseret=-1;
  }
  if (p->dir && rmdir(p->dir)==-1 && errno!=ENOENT) {
    perror("pipeclose could not remove temporary directory");
    fcloseret=-1;
  }
  free_pipe(p);
  return fcloseret;
}
