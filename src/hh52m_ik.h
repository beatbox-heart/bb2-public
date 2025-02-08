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

/* _(state_name, intial_condition) */
_(closed4,(1-n)*(1-n)*(1-n)*(1-n))
_(closed3,4*n*(1-n)*(1-n)*(1-n))
_(closed2,6*n*n*(1-n)*(1-n))
_(closed1,4*n*n*n*(1-n))
_(O_K,n*n*n*n)

/* _RATE(from_state, to_state, direct_TR, reverse_TR) */
_RATE(closed4,closed3,4*alpha_n,beta_n)
_RATE(closed3,closed2,3*alpha_n,2*beta_n)
_RATE(closed2,closed1,2*alpha_n,3*beta_n)
_RATE(closed1,O_K,alpha_n,4*beta_n)
