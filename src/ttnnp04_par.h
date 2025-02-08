/**
 * Copyright (C) (2010-2023) Vadim Biktashev, Irina Biktasheva et al. 
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

/* TenTussher-Noble-Noble-Panfilov 2004 model, from ten Tusscher's website */
/* From Main.cc */

/* ERRATA:
   We found a typo in the parameter values in the ms describing this model.
   Below we list the value the parameter should have and which is used
   in this implementation but that is wrong in the ms:
   GpCa , 0.825nS/pF
 */

/*-----------------------------------------------------------------------------
  ELECTROPHYSIOLOGICAL PARAMETERS:
-----------------------------------------------------------------------------*/

/* External concentrations */
_( Ko , 5.4)
_( Cao , 2.0)
_( Nao , 140.0)

/* Intracellular volumes */
_( Vc , 0.016404)
_( Vsr , 0.001094)

/* Calcium dynamics */
_( Bufc , 0.15)
_( Kbufc , 0.001)
_( Bufsr , 10.)
_( Kbufsr , 0.3)
_( taufca , 2.)
_( taug , 2.)
_( Vmaxup , 0.000425)
_( Kup , 0.00025)

/* Constants */
_( R , 8314.472)
_( F , 96485.3415)
_( T , 310.0)
_( RTONF , (R*T)/F)

/* Cellular capacitance          */
_( CAPACITANCE , 0.185)

/* Parameters for currents */
/* Parameters for IKr */
_( Gkr , 0.096)
/* Parameters for Iks */
_( pKNa , 0.03)
_( Gks , _ver(0.245,0.245,0.062))
/* Parameters for Ik1 */
_( GK1 , 5.405)
/* Parameters for Ito */
_( Gto , _ver(0.294,0.073,0.294))
/* Parameters for INa */
_( GNa , 14.838)
/* Parameters for IbNa */
_( GbNa , 0.00029)
/* Parameters for INaK */
_( KmK , 1.0)
_( KmNa , 40.0)
_( knak , 1.362)
/* Parameters for ICaL */
_( GCaL , 0.000175)
/* Parameters for IbCa */
_( GbCa , 0.000592)
/* Parameters for INaCa */
_( knaca , 1000)
_( KmNai , 87.5)
_( KmCa , 1.38)
_( ksat , 0.1)
_( n , 0.35)
/* Parameters for IpCa */
_( GpCa , 0.825)
_( KpCa , 0.0005)
/* Parameters for IpK) */
_( GpK , 0.0146)
