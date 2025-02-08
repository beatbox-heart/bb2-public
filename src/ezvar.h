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

/* _(Type,Name) */
_(GLXContext,theGLXContext)
_(XTextProperty,theWindowName)
_(XTextProperty,theIconName)
_(GLReal,marching_h) 		/* Size of marching cube in graphics coordinates */
_(GLReal,convert_to_phy) 	/* Convertion factor from graphics to physical coordinates */
_(unsigned int,ii)		/* Position of the 3rd vertex of */
_(unsigned int,jj)		/*     the current marching cube */
_(unsigned int,kk)		/*        on the simulation grid */
_(unsigned int,inc) 		/* The length of the marching (graphics) cube in  terms */
                                /*            of the simulation grid. The marching cube */
                                /*                         is defined by (ii,jj,kk) and */
                                /*               (ii-inc,jj-inc,kk-inc).  Must be >= 1. */

_(char *,busychar)		/* position of busy indicator in the window title */

_(char *,filament_buffer)
_(char *,filament_buffer_end)
_(char *,filament_buffer_current)
