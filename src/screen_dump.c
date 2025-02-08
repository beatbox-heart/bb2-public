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

/* Produce an XWD file from current X screen */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "bikt.h"
#include "bgi.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "sequence.h"
#include "byte.h"

#ifdef DEBUG
#undef DEBUG
#endif
#define DEBUG if (0) MESSAGE

typedef struct {
  sequence file;
  char fmt[16];
} STR;

RUN_HEAD(screen_dump)
{
  DEVICE_VAR(sequence,file)	   
  DEVICE_ARRAY(char,fmt)	   
  thisq(file); 		   DEBUG("thisq file is %s\n",file->name);
                           DEBUG("file %s is in the sequence\n",file->name);
  fclose(file->f);	   DEBUG("file %s closed\n",file->name);
  file->f=NULL;		   DEBUG("file handle %s nulled\n",file->name);
  dump_window(file->name,fmt); DEBUG("dump_window done to %s\n",file->name);
  file->f=fopen(file->name,"a"); DEBUG("file %s reopen, to make subsequent nextq happy\n",file->name);
  nextq(file); 		   DEBUG("nextq file is %s\n",file->name);
}
RUN_TAIL(screen_dump)

DESTROY_HEAD(screen_dump)
S->fmt[0]='\0';
DESTROY_TAIL(screen_dump)

CREATE_HEAD(screen_dump) {
  ACCEPTQ(file,"wb",NULL);
  ACCEPTS(fmt,"png");
  if (stricmp(S->fmt,"jpg") &&
      stricmp(S->fmt,"jpeg") &&
      stricmp(S->fmt,"gif") &&
      stricmp(S->fmt,"ppm") &&
      stricmp(S->fmt,"png")
      ) EXPECTED_ERROR("unknown format '%s'\n",S->fmt);
} CREATE_TAIL(screen_dump,0)

