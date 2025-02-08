/**
 * Copyright (C) (2010-2022) Vadim Biktashev, Irina Biktasheva et al. 
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

/* _(Type,Name,Dflt,miN,maX) */

/* Front-end, Beatbox-style parameters */
_(real,a,RNONE,RNONE,RNONE)	/* kinetics parameter a				*/
_(real,b,RNONE,RNONE,RNONE)	/* kinetics parameter b				*/
_(real,eps,RNONE,RNONE,RNONE)	/* kinetics parameter epsilon			*/
_(real,D,RNONE,0.,RNONE)        /* Isotropic diffusion coefficient		*/
_(real,Dv,0,0.,RNONE)           /* .., for the second variable			*/
_(real,hx,0,0.,RNONE)       	/* Space step					*/
_(real,ht,RNONE,0.,RNONE)       /* Time step					*/
_(real,delta,0.,RNONE,RNONE)	/* shortcut/round-up tolerance			*/

/* Compile-time macros made run-time flags */
_(int,explicit,0,0,1)		/* if 1 then explicit u-kinetics, else implicit	*/
_(int,split,1,0,1)		/* if 1 then split diffusion and kinetics steps	*/
_(int,manypoint,0,0,1)		/* equivalent of NINEPOINT and NINETEENPT flags	*/
_(int,pbc_x,0,0,1)		/* equivalent of PBC_x				*/
_(int,pbc_y,0,0,1)		/* equivalent of PBC_y				*/
_(int,pbc_z,0,0,1)		/* equivalent of PBC_z				*/
