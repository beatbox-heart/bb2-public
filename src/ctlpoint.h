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

_("version",	&V_ ,	&VERSTRING,	VERLEN, char, MPI_CHAR, t_str,	0)
_("xmax",   	&xm_,	&xmax,		1, long, MPI_LONG, t_int,	1)
_("ymax",   	&ym_,	&ymax,		1, long, MPI_LONG, t_int,	1)
_("zmax",   	&zm_,	&zmax,		1, long, MPI_LONG, t_int,	1)
_("vmax",   	&vm_,	&vmax,		1, long, MPI_LONG, t_int,	1)
_("mpi_nx",   	&mpinx_,&mpi_nx,	1, long, MPI_LONG, t_int,	0)
_("mpi_ny",   	&mpiny_,&mpi_ny,	1, long, MPI_LONG, t_int,	0)
_("mpi_nz",   	&mpinz_,&mpi_nz,	1, long, MPI_LONG, t_int,	0)
