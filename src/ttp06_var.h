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

/* Variables as listed in Variables.h */
/* Initial conditions from Variables.cc unless stated otherwise */

/* Non-gating variables */
_(Volt,-86.2)	/* from Main.cc */
_(Cai,0.00007)	/* from Main.cc */
_(CaSR,1.3)	/* from Main.cc */
_(CaSS,0.00007)	/* from Main.cc */
_(Nai,7.67)	/* from Main.cc */
_(Ki,138.3)	/* from Main.cc */


/* states of voltage and time dependent gates */
/* INa */
_(M, 0)
_(H, 0.75)
_(J, 0.75)
/* IKr */
/* IKr1 */
_(Xr1, 0.)
/* IKr2 */
_(Xr2, 1.)
/* IKs */
_(Xs, 0.)
/* Ito1 */
_(R, 0.)
_(S, 1.)
/* ICa */
_(D, 0.)
_(F, 1.)
_(F2,1.)
_(FCass,1.)
/* Irel */
_(RR,1.)
