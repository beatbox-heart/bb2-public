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


/* SYSTEM-DEPENDENT HEADERS */

/* Currently all non-unix variants eliminated     */
/* as they have not been supported for many years */
#ifndef __SYSTEM_H
#define __SYSTEM_H

#define UNIX

/******************************************************/
#ifdef UNIX
# include <unistd.h>

# define NULLFILE "/dev/null"

# ifndef max
#   define max(a,b) (((a)>(b))?(a):(b))
#   define min(a,b) (((a)<(b))?(a):(b))
# endif

# ifndef MAXPATH
#   ifdef FILENAME_MAX
#     define MAXPATH FILENAME_MAX
#   else  
#     define MAXPATH 4096
#   endif
# endif

  #define fabsl fabs
  void fcloseall(void);
  int kbhit(void);
  void sound (int a);
  void delay (int a);
  void nosound (void);
  int getch(void);
  int Strcmp (const char *s1, const char *s2); 
  int stricmp (const char *s1, const char *s2);
#endif

#endif
