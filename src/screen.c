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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "k_.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "bgi.h"

#include "screen.h"

extern INT XMAX;
extern INT YMAX;
extern INT WINX;
extern INT WINY;
extern INT online;

int screen_already_created=0;

typedef struct {
  int XMAX, YMAX, WINX, WINY, online;
} STR;

/* makescreen is performed once just when reading from parameter file */
int makescreen (char *w) {
  STR *S=(STR *)Malloc(sizeof(STR));
  if(!S) EXPECTED_ERROR("not enough memory");
  if (screen_already_created) EXPECTED_ERROR("screen already made");

  /* ACCEPTI implicitly defines eponymous local variable, */
  {
    ACCEPTI(XMAX,639,1,INONE);
    ACCEPTI(YMAX,479,1,INONE);
    ACCEPTI(WINX,-1,INONE,INONE);
    ACCEPTI(WINY,-1,INONE,INONE);
    ACCEPTI(online,0,0,1);
  }
  /* so we have to assign to the global ones outside the block: */
  XMAX=S->XMAX;
  YMAX=S->YMAX;
  WINX=S->WINX;
  WINY=S->WINY;
  online=S->online;
  
  FREE(S);
  if (!opengraph()) {
    graphon=0;
    return 0;
  }
  update_graph();
  return (1);
}

