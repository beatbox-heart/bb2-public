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


_(V 	,0,	-93.0)		/* voltage */
_(Nai 	,1,	5.0)		/* internal sodium */
/*_(Ko 	,2,	4.0)            // external potassium - constant */
_(Cai 	,3,	0.00000763)     /* internal calcium */
/* 4=i_f activation */
_(x 	,5,	0.001)          /* iK activation */
_(d 	,6,	0.0)            /* iCa activation */
_(f	,7,	1.0)            /* iCa inactivation */
_(h	,8,	0.9951)         /* iNa inactivation */
_(m	,9,	0.0015)         /* iNa activation */
_(Ki	,10,	140.0)          /* internal potassium */
_(Caup	,11,	0.3013)         /* SR uptake store calcium */
_(Carel	,12,	0.2989)         /* SR release store calcium */
/* 13=repriming variable */
/* 14=slow iCa inactivation */
/* 15=activation of maintained iCa */
/* 16=inactivation of maintained iCa */
_(r	,17,	1.0)
_(q	,18,	0.0)
_(Cacmd	,19,	0.0003)
_(Catrp	,20,	0.0002)
/* 21=magnesium bound to troponin */
/* 22=calcium bound to M troponin */
/* _(Cao	,23,	2.0) */
/* 24=mitochondrial calcium */
/* 25=precursor fraction for calcium release - excluded */
_(fact	,26,	0.0023)		/* activator fraction for calcium release */
_(fprod	,27,	0.0291)		/* product fraction for calcium release */
/* 28=light chain conformation in contraction model  - decoupled */
/* 29=cross bridge formation in contraction model - decoupled */
/* 30=total cytosol calcium (including cytosolic buffered calcium,) - decoupled */
/* 31=T type calcium channel activation */
/* 32=T type calcium channel inactivation */
/* 33=calcium bound to calcium indicator */
/* 34=internal ATP concentration - now constant paramer */
_(INa	,56,	0)
_(IK	,57,	0)
_(IK1	,58,	0)
_(iCa	,59,	0)
