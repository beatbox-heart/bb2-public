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

/* TenTussher-Panfilov 2006 model, from ten Tusscher's website */
/* From Main.cc */


/* We discovered some typo's in parameter values in the publication:
   Vrel=40.8; should be 0.102
   k4=0.000015; should be 0.005
   Vc=16.404; should be 16404
   Vsr=1.094; should be 1094
   Vss=0.05468; should be 54.68

*/   


/*-----------------------------------------------------------------------------
  ELECTROPHYSIOLOGICAL PARAMETERS:
-----------------------------------------------------------------------------*/

/* External concentrations */
_(Ko,5.4)
_(Cao,2.0)
_(Nao,140.0)

/* Intracellular volumes */
_(Vc,0.016404)
_(Vsr,0.001094)
_(Vss,0.00005468)

/* Calcium buffering dynamics */
_(Bufc,0.2)
_(Kbufc,0.001)
_(Bufsr,10.)
_(Kbufsr,0.3)
_(Bufss,0.4)
_(Kbufss,0.00025)

/* Intracellular calcium flux dynamics */
_(Vmaxup,0.006375)
_(Kup,0.00025)
_(Vrel,0.102) /* 40.8) */
_(k1_,0.15)
_(k2_,0.045)
_(k3,0.060)
_(k4,0.005) /* 0.000015) */
_(EC,1.5)
_(maxsr,2.5)
_(minsr,1.)
_(Vleak,0.00036)
_(Vxfer,0.0038)

/* Cellular capacitance          */
_(CAPACITANCE,0.185)

/* Parameters for currents */
/* Parameters for IKr */
_(Gkr,0.153)
/* Parameters for Iks */
_(pKNa,0.03)
#ifdef EPI
_(Gks,0.392)
#endif
#ifdef ENDO
_(Gks,0.392)
#endif
#ifdef MCELL
_(Gks,0.098)
#endif
/* Parameters for Ik1 */
_(GK1,5.405)
/* Parameters for Ito */
#ifdef EPI
_(Gto,0.294)
#endif
#ifdef ENDO
_(Gto,0.073)
#endif
#ifdef MCELL
_(Gto,0.294)
#endif
/* Parameters for INa */
_(GNa,14.838)
/* Parameters for IbNa */
_(GbNa,0.00029)
/* Parameters for INaK */
_(KmK,1.0)
_(KmNa,40.0)
_(knak,2.724)
/* Parameters for ICaL */
_(GCaL,0.00003980)  
/* Parameters for IbCa */
_(GbCa,0.000592)
/* Parameters for INaCa */
_(knaca,1000)
_(KmNai,87.5)
_(KmCa,1.38)
_(ksat,0.1)
_(n,0.35)
/* Parameters for IpCa */
_(GpCa,0.1238)
_(KpCa,0.0005)
/* Parameters for IpK */
_(GpK,0.0146)
