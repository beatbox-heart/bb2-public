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

/* Choice of the type of MPI IO to use:
        0 for blocking, 
        1 for contiguous. 
   In our case, contiguous is the one that gives scaling.
*/
#define blocking_or_contiguous 1 

#if (blocking_or_contiguous==0) 
#define MPI_FILE_WRITE MPI_File_write
#else
#define MPI_FILE_WRITE MPI_File_write_all
#endif

/* Whether to flush buffers:
     1 for use fflush (local workstations) and 
     0 for not using fflush (HPC systems, clusters)
*/
#define fflush_or_not 1

/* NB here the positive answer to "if .. or" is presumed opposite     */
/* to that presumed in the previous case, but that's how MCF wrote it */
#if fflush_or_not
#define FFLUSH fflush
#else
#define FFLUSH nofflush
#endif
