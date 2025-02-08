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

/* Modifications for ezview: 2010-18, Vadim Biktashev v.n.biktashev@exeter.ac.uk */

#ifndef _EZGRAPH3D_
#define _EZGRAPH3D_

#define WINDOW_TITLE "EZView"

#define WM_CTRLS_POS 0      /* If WM_CTRL_POS is 0 then the window is placed
			     * at (WINX, WINY). If WM_CTRL_POS is 1 then WINX
			     * and WINY are ignored and the location is
			     * determined by the window manager. */
#define WINSIZE      640   /* Window is square of this size in pixels, initially. */
#define PLOT_SIZE    1.2    /* This controls the size of the simulation
			     * volume in the view port: >1.0 for larger
			     * size, <1.0 for smaller size. */

#define PERSPECTIVE  1      /* 1 is perspective projection, 0 orthographic.
			     * see myReshape(). */

#define START_PAUSED 0      /* If 1 then window is opened in paused mode
			     * showing initial condition. */

#define ROTATE_SCALE   0.5   /* Number of degrees of rotation corresponding to
			      * one pixel of mouse motion. */
#define ROTATE_STEP   1.0   /* Number of degrees of rotation corresponding to
			     * one click of an arrow button. */
#define deg (M_PI/180.0)       /* rad per angle degree */
#define DISTANCE_FAC  1.01   /* By how much distance is increased/decreased by
			     * one click of + or -  button. */

#define many (10)	     /* shift+arrow moves this many pixels/degrees */

#define BUSYRUNNING ' '
#define BUSYDRAWING '-'
#define BUSYWAITING '+'
#define BUSYTURNING '='

/* --------------------------------------------- 
 * You probably should not change anything below 
 * --------------------------------------------- */

#define TRUE             1   
#define FALSE            0

#define U_FIELD          0     /* These are used to determine which field */
#define V_FIELD          1     /* (if any) is being plotted */
#define NO_FIELD        -1

#endif /*  _EZGRAPH3D_  */
