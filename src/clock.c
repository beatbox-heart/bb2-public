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

/* SHOW/PRINT THE STATE OF A TIMER */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"

typedef struct {
  BGIWindow wnd;
} STR;

RUN_HEAD(clock)
{
#if MPI
  MESSAGE("\rt=%-6ld",t);
#elif defined(X11)
  DEVICE_CONST(BGIWindow,wnd);
  char l[80];
  snprintf (l,80,"t=%-6ld",t);
  crt_text(l,wnd.row0,wnd.col0,wnd.color);
#else
  MESSAGE("\rt=%-6ld",t);
#endif
}
RUN_TAIL(clock)

DESTROY_HEAD(clock)
DESTROY_TAIL(clock)

CREATE_HEAD(clock)
{
  DEVICE_IS_SPACELESS;
  ACCEPT_WINDOW(wnd);
}
CREATE_TAIL(clock,0)


