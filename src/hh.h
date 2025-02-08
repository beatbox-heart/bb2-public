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

/* List of pars  */
/* Name	   Value */
_(C,        1.   )
/* maximum conductance of current I_Na (mS/cm^2) */
_(gNa, 120);
/* maximum conductance of current I_K (mS/cm^2)*  */
_(gK, 36);
/* maximum conductance of leakage current: I_l (mS/cm^2) */
_(gl, 0.3);

#if AUTHENTIC
/* the sign from HH1952 is according new convention */
/* reversal potential of current I_Na (mV)  */
_(VNa, 115);
/* reversal potential of current I_K (mV)  */
_(VK, -12);
/* reversal potential of current I_l (mV)  */
_(Vl, 10.613);
#else
_(VNa,   -115.   )
_(VK,      12.   )
_(Vl,     -10.613)
/* artifical accelration/deceleration for research purposes */
_(km, 	     1.0 )
_(kn, 	     1.0 )
_(kh, 	     1.0 )
#endif /* AUTHENTIC */
