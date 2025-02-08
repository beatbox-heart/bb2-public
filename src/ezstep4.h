/**
 * Copyright (C) (2010-2018) Vadim Biktashev, Irina Biktasheva et al. 
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

switch (ezdim) {
          case 0:
            #include "ezstep0d.h"
            break;
          case 1:
            EXPECTED_ERROR("%s: 1D is not implemented\n",__FILE__);
            break;
          case 2:
            #include "ezstep2d.h"
            Impose_boundary_conditions_2D(s, S, u, sigma_u);
            if (v_diff_on)        
              Impose_boundary_conditions_2D(s, S, v, sigma_v);
            break;
          case 3:
            #include "ezstep3d.h"
            Impose_boundary_conditions_3D(s, S, u, sigma_u);
            if (v_diff_on)        
              Impose_boundary_conditions_3D(s, S, v, sigma_v);
            break;
          default:
            EXPECTED_ERROR("unexpected dimensionality %d in %s\n",dim,__FILE__);
            break;
          } /* switch dim */
  
