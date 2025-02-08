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


_(V,	    -74.2525*mV)	/* transmembr voltage */
_(Nac,	    130.011*mM)		/* cleft sodium */
_(Kc,	    5.3581*mM)		/* cleft potassium */
_(Cac,	    1.8147*mM)		/* cleft calcium */
_(Nai,	    8.5547*mM)		/* internal sodium */
_(Ki,	    129.435*mM)		/* internal potassium */
_(Cai,	    6.7290e-5*mM)	/* internal calcium */
_(Cad,	    7.2495e-5*mM)	/* ??? some calcium */
_(Caup,	    0.6646*mM)		/* uptake calcium */
_(Carel,    0.6465*mM)		/* release calcium */
_(m,	    3.2017e-3)		/* INa activation */
_(h1,	    0.8814)		/* INa inactivation 1 */
_(h2,	    0.8742)		/* INa inactivation 2 */
_(dL,	    WHICH(1.3005e-5,1.3e-5))	/* ICa activation */
_(fL1,	    0.9986)		/* ICa inactivation 1 */
_(fL2,	    0.9986)		/* ICa inactivation 2 */
_(r,	    1.0678e-3)		/* It activation  */
_(s,	    0.9490)		/* It inactivation  */
_(rsus,	    1.5949e-4)		/* Isus activation  */
_(ssus,	    0.9912)		/* Isus inactivation  */
_(n,	    4.8357e-3)		/* ???  */
_(pa,	    0.0001)		/* ???  */
_(F1,	    0.4284)		/* ???	*/
_(F2,	    0.0028)		/* ???	*/
_(OC,	    0.0275)		/* ???	*/
_(OTC,	    0.0133)		/* ???	*/
_(OTMgC,    0.1961)		/* ???	*/
_(OTMgMg,   0.7094)		/* ???	*/
_(OCalse,   WHICH(0.4369,0.4469)) /* ???	*/
