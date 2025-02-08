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

/* Hold on for a number of seconds or until user presses a key. */
/* Has no effect in parallel mode. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"

#define nanosecond (1.0e-9)

typedef struct {
  real seconds;
  int echo;
} STR;

RUN_HEAD(pause)
#if MPI
#else
  DEVICE_CONST(real,seconds)
  DEVICE_CONST(int,echo)
  if (seconds>=0) {
    struct timespec tim;
    tim.tv_sec = (long) floor(seconds);
    tim.tv_nsec = (long) (fmod(seconds,1.0)/(nanosecond));
    if (echo) MESSAGE("waiting for %f seconds...\n",seconds);
    if (nanosleep(&tim, NULL) < 0) MESSAGE("Nano sleep system call failed \n");
  } else {
    printf("press Enter key in this window to continue\n");
    fgetc(stdin);
  }
#endif
RUN_TAIL(pause)

DESTROY_HEAD(pause)
DESTROY_TAIL(pause)

CREATE_HEAD(pause) {
  ACCEPTR(seconds,-1.0,RNONE,RNONE);
  ACCEPTI(echo,(seconds<0)?1:0,0,1);
} CREATE_TAIL(pause,0)
