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

real alp_x1 = (0.0005*exp(0.083*((V)+50))/(exp(0.057*((V)+50))+1.0));
real bet_x1 = (0.0013*exp(-0.06*((V)+20))/(exp(-0.04*((V)+20))+1.0));
real alp_m =  (-((V)+47)/(exp(-0.1*((V)+47))-1.0));
real bet_m =  (40*exp(-0.056*((V)+72)));
real alp_h =  (0.126*exp(-0.25*((V)+77)));
real bet_h =  (1.7/(exp(-0.082*((V)+22.5))+1.0));
real alp_j =  (0.055*exp(-0.25*((V)+78))/(exp(-0.2*((V)+78))+1.0));
real bet_j =  (0.3/(exp(-0.1*((V)+32))+1.0));
real alp_d =  (0.095*exp(-0.01*((V)-5))/(exp(-0.072*((V)-5))+1.0));
real bet_d =  (0.07*exp(-0.017*((V)+44))/(exp(0.05*((V)+44))+1.0));
real alp_f =  (0.012*exp(-0.008*((V)+28))/(exp(0.15*((V)+28))+1.0));
real bet_f =  (0.0065*exp(-0.02*((V)+30))/(exp(-0.2*((V)+30))+1.0));
real ik1std = (V==-23.)
  ?(2.817939549)
  :(0.35*(
    4.0*(exp(0.04*(V+85))-1.0)/(exp(0.08*(V+53))+exp(0.04*(V+53)))
    +
    0.2*(V+23)/(1.0-exp(-0.04*(V+23)))
  ));
real gix1=0.8*(exp(0.04*(V+77))-1.0)/exp(0.04*(V+35));
