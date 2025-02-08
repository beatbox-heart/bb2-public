/**
 * Copyright (C) (2010-2025) Vadim Biktashev, Irina Biktasheva et al. 
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
_(INT,winx,-1,LNONE,LNONE);	    /* position (bottom,left) of the window */
_(INT,winy,0,LNONE,LNONE);	    /*   on the display */
_(INT,width,640,LNONE,LNONE);	    /* width of the window */
_(INT,height,640,LNONE,LNONE);	    /* hight of the window */
_(REAL,theta0,20,RNONE,RNONE)       /* starting */
_(REAL,phi0,60,RNONE,RNONE)         /* position  */
_(REAL,psi0,0,RNONE,RNONE)          /*   of the  */
_(REAL,distance0,5.0,RNONE,RNONE)   /*   camera    */
_(INT,autostart,1,0,1) 		    /* autostart flag */
