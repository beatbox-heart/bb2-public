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

/* MODULE FOR STOPPING THE JOB */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"

extern int ReturnCode;

typedef int STR;	/* dummy declaration; for compatibility only */

RUN_HEAD(stop)
  if (mpi_rank==0) MESSAGE("\nSTOP AT %ld\n",t);
  ReturnCode=EXIT_SUCCESS;
  return(0);
RUN_TAIL(stop)

DESTROY_HEAD(stop)
DESTROY_TAIL(stop)

CREATE_HEAD(stop)
CREATE_TAIL(stop,0)

