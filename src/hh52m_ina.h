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
_(C3,(1-m)*(1-m)*(1-m)*h)
_(C2,3*m*(1-m)*(1-m)*h)
_(C1,3*m*m*(1-m)*h)
_(O_Na,m*m*m*h)

_(I3,(1-m)*(1-m)*(1-m)*(1-h))
_(I2,3*m*(1-m)*(1-m)*(1-h))
_(I1,3*m*m*(1-m)*(1-h))
_(I0,m*m*m*(1-h))

/* _RATE(from_state, to_state, direct_TR, reverse_TR) */
_RATE(C3,C2,3*alpha_m,beta_m)
_RATE(C2,C1,2*alpha_m,2*beta_m)
_RATE(C1,O_Na,alpha_m,3*beta_m)

_RATE(I3,I2,3*alpha_m,beta_m)
_RATE(I2,I1,2*alpha_m,2*beta_m)
_RATE(I1,I0,alpha_m,3*beta_m)

_RATE(I3,C3,alpha_h,beta_h)
_RATE(I2,C2,alpha_h,beta_h)
_RATE(I1,C1,alpha_h,beta_h)
_RATE(I0,O_Na,alpha_h,beta_h)
