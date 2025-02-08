/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
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


/* SYSTEM-DEPENDENT PART */

#include "system.h"

#ifdef UNIX

#include <string.h>
#include <stdio.h>
#include <ctype.h>

void fcloseall(void) {;}
void nosound (void) {;}
void sound (int a) {;}
void delay (int a) {;}
int kbhit(void) {return(0);}

#include "bgi.h"
int getch(void) {
  if(graphon)update_graph(); 
  return getchar();
}

int Strcmp(const char *s1, const char *s2){
  if (!s1 &&!s2) return 0;
  else if (!s1) return 1;
  else if (!s2) return -1;
  else return strcmp(s1, s2);
}

int stricmp(const char *s1, const char *s2) {
  unsigned char c1, c2;
  if (!s1 &&!s2) return 0;
  else if (!s1) return 1;
  else if (!s2) return -1;
  else for(;;) {
    c1=*(s1++); c2=*(s2++); /* both at once! */
    if (!c1 || !c2) break;
    if (c1 != c2) break;
  } 
  return ((signed int)c1-(signed int)c2);
}

#endif

