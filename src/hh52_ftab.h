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

/* gate transition rates */
real alpha_n = 0.01 * (-V + 10.0) /
  (exp ((-V + 10.0) / 10.0) - 1.0);
real beta_n = 0.125 * exp (-V / 80.0);

real alpha_m = 0.1 * (-V + 25.0) /
  (exp ((-V + 25.0) / 10.0) - 1.0);
real beta_m = 4.0 * exp (-V / 18.0);

real alpha_h = 0.07 * exp (-V / 20.0);
real beta_h = 1.0 / (exp ((-V + 30.0) / 10.0) + 1.0);
