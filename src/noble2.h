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

/* 
 * Fastest variables:
 *  7 [6] --- iCa activation
 *  8 [7] --- iCa inactivation
 *  9 [8] --- iNa inactivation
 *  10 [9] -- iNa activation
 * - duplicated to simulate anodal and cathodal parts of the
 * cell membrane, and new variable I added to simulate potential difference 
 * between the parts due to external electric current. 
 * Allocation of variables in Beatbox v-indices:
 * Index		Contents, in terms of HEART indices
 * 0			0 ( average trmbrn voltage )
 * 1			exogenic potential difference
 * 2--5			6--9 of anodal part
 * 6--9			6--9 of cathodal part
 * 10--14		1--5 common of both parts
 * 15--neqn+4		10--N-1 common of both parts
 * 
 * Ergo: neqn2=N=neqn+5, where 5=exogenic potential+m+h+d+f
 */
 
#ifndef neqn
#  define neqn 60
#elif neqn != 60
#  error neqn != 60
#endif

extern int neqn2;
extern int h2p[neqn];
extern int h2n[neqn];

#ifdef OWN
int neqn2=neqn+5;
int h2p[neqn]={0,10,11,12,13,14,2,3,4,5,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64}; 
int h2n[neqn]={0,10,11,12,13,14,6,7,8,9,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64}; 
#endif
