/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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

/* Generic error reporting routines */

#ifndef __ERROR_H
#define __ERROR_H

/* Most urgent fail: */
/* Print msg on error and do a long jump to main */
/* This is required when the failure may occur not in all processes */
#if MPI
#define ABORT(...) {							\
  URGENT_MESSAGE("\nOne-process error in %s:%u: ",__FILE__,__LINE__);	\
  beatbox_abort(__VA_ARGS__);						\
  }
#else
#define ABORT(...) {							\
  URGENT_MESSAGE("\nError in %s:%u: ",__FILE__,__LINE__);		\
  beatbox_abort(__VA_ARGS__);						\
  }
#endif
void beatbox_abort(char const *fmt, ...);

/* Lower class of error. 
 * Asks Beatbox to fail gracefully, but can get by
 * with only process 0 producing a message.
 *
 * Useful for user errors in the script, etc.
 */
#if MPI
#define EXPECTED_ERROR(...) {					\
  MESSAGE("\nCommon error in %s:%u: ",__FILE__,__LINE__);       \
  MESSAGE(__VA_ARGS__);MARK_ERROR}
#else
#define EXPECTED_ERROR(...) {					\
  MESSAGE("\nError in %s:%u: ",__FILE__,__LINE__);       \
  MESSAGE(__VA_ARGS__);MARK_ERROR}
#endif

#define ASSERT(p) { if(0==(p)) EXPECTED_ERROR("Assertion failed:\n %s in %s:%d\n",#p,__FILE__,__LINE__); }

#endif
