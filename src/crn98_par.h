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

_(  fscale1 ,	1.0) //				( 1 )
_(  fscale2 ,	1.0) //				( 1 )
_(  vcell ,	 20100.0) //			( um3 )
_(  vi ,	 vcell*0.68) //			( 1 )
_(  vup ,	 fscale1*0.0552*vcell) //	( 1 )
_(  vrel ,	 fscale2*0.0048*vcell) //	( 1 )
/* _(  T ,	 	 310) //			( Celcius ) */
/* _(  Tfac ,	 3) //				( 1 )  */
_(  Csp ,	 1e+6) //			( pF/cm2 )
_(  Cm ,	 100.0) //			( pF )
/* _(  F ,	 	 96.4867) //			( coul/mmol ) */
/* _(  R ,	 	 8.3143) //			( J K-1 mol-1 ) */
_(  kb ,	 5.4) //			( mM )
_(  nab ,	 140) //			( mM )
_(  cab ,	 1.8) //			( mM )
_(  nac ,	 nab) //			( mM )
_(  cac ,	 cab) //			( mM )
_(  kc ,	 kb) //				( mM )
_(  gna ,	 7.8) //			( nS/pF ) 
_(  gto ,	 0.1652) //			( nS/pF ) 
_(  gkr ,	 0.029411765) //		( nS/pF )
_(  gks ,	 0.12941176) //			( nS/pF )
_(  gcaL ,	 0.12375) //			( nS/pF ) 
_(  ErL ,	 65.0) //			( mV )
_(  gk1 ,	 0.09) //			( nS/pF ) 
_(  gbna ,	 0.0006744375) //		( nS/pF ) 
_(  gbk ,	 0.0) //			( 1 )
_(  gbca ,	 0.001131) //			( nS/pF ) 
_(  inakbar ,	 0.59933874) //			( pA/pF )
_(  kmnai ,	 10.0) //			( mM )
_(  kmko ,	 1.5) //			( mM )
_(  icapbar ,	 0.275) //			( pA/pF )
_(  kmcap ,	 0.0005) //			( mM )
_(  knacalr ,	 1600.0) //			( pA/pF )
_(  kmnalr ,	 87.5) //			( mM )
_(  kmcalr ,	 1.38) //			( mM )
_(  ksatlr ,	 0.1) //			( 1 )
_(  trpnbar ,	 0.070) //			( mM )
_(  cmdnbar ,	 0.050) //			( mM )
_(  csqnbar ,	 10) //				( mM )
_(  kmcsqn ,	 0.8) //			( mM )
_(  kmcmdn ,	 0.00238) //			( mM )
_(  kmtrpn ,	 0.0005) //			( mM )
_(  grelbar ,	 30.0) //			( ms-1 ) 
_(  kmup ,	 0.00092) //			( mM ) 
_(  iupbar ,	 0.005) //			( mM/ms )
_(  caupmax ,	 15.0) //			( mM )
_(  kupleak ,	 iupbar/caupmax) //		( ms-1 )
_(  tautr ,	 180.0) //			( ms )
_(  gkur_scale , 1.0) //			( 1 ) /* gkur_scale as in Marshall et al. Pflugers Arch. ??(2011):?? */
