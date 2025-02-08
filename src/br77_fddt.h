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

/* (non-tabulated) currents */
real es=-82.3-13.0287*log(ca);
real ik1=gk1*ik1std;
real ix1=gix1*x1;
real ina=(gna*m*m*m*h*j+gnac)*(V-ena);
real is=gs*d*f*(V-es);

/* The rates of the non-gate dynamic variables: the main result */
real d_V = -(ik1+ix1+ina+is);
real d_ca = -0.0000001*is + 0.07*(0.0000001-ca);
