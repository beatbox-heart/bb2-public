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

/* _(Type,Name,Init) */

_(int,nx,xmax-2)				/* netto sizes of the grid  */
_(int,ny,ymax-2)				/* consistent with EZSCROLL */
_(int,nz,zmax-2)				/* eponymous variables      */
_(int,nv,vmax+geom_vmax)			/* number of fields 	    */
_(real,taboo,1e100) 				/* effective NaN: "prohibited" value for a field */

_(Display *,theDisplay,NULL)			/* standard */
_(int,theDWidth,0)				/* X11 */
_(int,theDHeight,0)				/* stuff */
_(Window,theWindow,0)				/* ... */

_(int,horspace,0);				/* window is smaller than the display */
_(int,vertspace,0);				/*   by these lengths */

_(ezmode_type,ezmode,MODE_SIMULATING);		/* mode flag */
_(int,lastx,0);					/* last mouse position,  */
_(int,lasty,0);					/*   used by Rotate() and Mouse_down() */

_(int,NORM_LIST,0)				/* full resolution surface GL list */
_(int,ROT_LIST,0)				/* low resolution surface GL list */
_(int,FLM_LIST,0)				/* filament GL list */
#if FIBRES
_(int,FIB_LIST,0)				/* fibres */
#endif
/* To do .. */
#if MARKER
_(int,MARKER_LIST,0)				/* marker GL list */
#endif

_(int,draw_mode,OFF) /* draw_mode can be in one of six states. start * with it off */
_(GLReal,last_lwidth,0)				/* line width previously set */

