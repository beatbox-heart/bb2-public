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
#include "k_.h"
#include "bikt.h"

typedef struct {
  char format[80];
  char code[80];
  pp_fn compiled;
  BGIWindow wnd;
} STR;

RUN_HEAD(k_clock)
{
  DEVICE_CONST(BGIWindow,wnd);
  DEVICE_ARRAY(char, format)
  DEVICE_CONST(pp_fn, compiled)
  char l[80];
  k_on();
  if (strchr(format,'%')) snprintf(l,80,format,*(real *)execute(compiled));
  else snprintf(l,80,"%s%s",format,prt(execute(compiled),t_real));
  k_off();
  crt_text(l,wnd.row0,wnd.col0,wnd.color);
}
RUN_TAIL(k_clock)

DESTROY_HEAD(k_clock)
{
  FREE(S->compiled);
}
DESTROY_TAIL(k_clock)

CREATE_HEAD(k_clock)
{
  DEVICE_IS_SPACELESS;
  ACCEPT_WINDOW(wnd);
  k_on();
  ACCEPTS(format,"t=%-6lf");
  ACCEPTS(code,"t");
  S->compiled=compile(S->code,deftb,t_real); CHK(S->code);
  k_off();
}
CREATE_TAIL(k_clock,0)


