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


_(V,-86.5)
_(nai,10.0)
_(ki,145.0)
_(cai,caiinit)
_(cansr,cansr_max / (1.0 + (kmup / caiinit)))
_(cajsr,cansr_max / (1.0 + (kmup / caiinit)))
_(caiontotold,0.0)
_(dcaiontot,0.0)
_(dvdt,0.0)
_(dcai,0.0)
_(csqn,0.0)
_(tjsrol,0.0)
_(h, 0.0)
_(j, 0.0)
_(d, 0.0)
_(f, 0.0)
_(x, 0.0)
_(m, 0.0)
_(test, 0.0)
_(tdvdtmax, 0.0)
_(delta_cai2, 0.0)
_(tcicr,0.0)
_(grelcicr_max,0.0)
_(xi,0.0)


/* If CELLTEST is defined, the following variables are declared, thus allowing
   their time courses to be traced in the script.
  
   All of the following variables have been initialised to 0, as they first occur in the
   LHS of equations, hence their normal declaration in lrd_step.h.

   The corresponding declarations in lrd_step.h will be removed automatically.

   The variable declarations below are grouped by the figures in which they first appear
   in C.-H. Luo and Y.ÊRudy. A dynamic model of the cardiac ventricular action potential. 
   I. simulations of ionic currents and concentration changes. 
   Circulation Research, 74(6):1071Ð1096, Jun 1994. */

#ifdef CELLTEST
	/*  Fig 13 */
	_(ina,0)
	_(ica,0)
	_(ik,0)
	_(iv,0)
	_(inaca,0)
	
	/*  Fig 14 */
	_(trpn,0)
	_(cmdn,0)
	
	/*  Fig 15 */
	_(irelcicr,0)
	_(ireljsrol,0)
	_(iup,0)
	_(ileak,0)
	_(itr,0)
	
	/*  Fig 16 */
	_(fca,0)
	
	/*  Fig 17 */
	_(icat,0)
	_(icak,0)
	_(icana,0)
	
	/*  Fig 18 */
	_(ik1,0)
	_(ikp,0)
	_(inak,0)
	_(ipca,0)
	_(inab,0)
	_(icab,0)
#endif
