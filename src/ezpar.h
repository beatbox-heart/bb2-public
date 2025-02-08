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

/* Parameters that are definable in bbs and changeable in dialogue */

/* _(Type,Name,Init,reRead,reView,reLight,reMake,reDraw) */
/* T  N    I            R V L M D */
_(INT,imin,s.x0,	0,0,0,1,1) /* 3D          */
_(INT,imax,s.x1,        0,0,0,1,1) /*   subset    */
_(INT,jmin,s.y0,	0,0,0,1,1) /*   of        */
_(INT,jmax,s.y1,        0,0,0,1,1) /*   the box   */
_(INT,kmin,s.z0,	0,0,0,1,1) /*   to be     */
_(INT,kmax,s.z1,        0,0,0,1,1) /*   displayed */
_(INT,norm_res,1,	0,0,0,1,1) /* view and run mode resolution */
_(INT,rot_res,1,	0,0,0,1,1) /* rotate mode resolution */

/* T  N    I            R V L M D */
_(INT,show_surface,1,	0,0,0,0,1) /* whether to draw the surface */
_(INT,ulayer,0,		0,0,0,1,1) /* isosurface channel in input files */
_(REAL,uc,0,		0,0,0,1,1) /* isosurface critical value */
_(INT,vlayer,1,		0,0,0,1,1) /* surface coloring layer in input files */
_(REAL,vc,0,		0,0,0,1,1) /* .. default value */
_(REAL,vmin,0,		0,0,0,1,1) /* coloring range lower end */
_(REAL,maxv,1,		0,0,0,1,1) /* coloring range lower end (vmax reserved) */
_(INT,rlayer,0,		0,0,0,1,1) /* red       channel in input files */
_(REAL,rmin,0,		0,0,0,1,1) /* red       range lower end */
_(REAL,rmax,1,		0,0,0,1,1) /* red       range upper end */
_(INT,glayer,0,		0,0,0,1,1) /* green     channel in input files */
_(REAL,gmin,0,		0,0,0,1,1) /* green     range lower end */
_(REAL,gmax,1,		0,0,0,1,1) /* green     range upper end */
_(INT,blayer,0,		0,0,0,1,1) /* blue      channel in input files */
_(REAL,bmin,0,		0,0,0,1,1) /* blue      range lower end */
_(REAL,bmax,1,		0,0,0,1,1) /* blue      range upper end */
_(INT,alayer,0,		0,0,0,1,1) /* alpha     channel in input files */
_(REAL,amin,0,		0,0,0,1,1) /* alpha     range lower end */
_(REAL,amax,1,		0,0,0,1,1) /* alpha     range upper end */

/* T  N    I            R V L M D */
_(INT,color_mode,1,	0,0,0,1,1) /* surface colouring algorithm */
_(REAL,alphamin,0,	0,0,0,1,1) /* minimal opaqueness */
_(REAL,alphamax,0,	0,0,0,1,1) /* maximal opaqueness */

/* T  N    I            R V L M D */
_(INT,show_filament,0,	0,0,0,1,1) /* whether to draw filaments */
_(INT,layer1,0,		0,0,0,1,1) /* sing layer 1 */
_(REAL,const1,0,	0,0,0,1,1) /* sing const 1 */
_(INT,layer2,1,		0,0,0,1,1) /* sing layer 2 */
_(REAL,const2,0,	0,0,0,1,1) /* sing const 2 */
_(REAL,flm_r,1,		0,0,0,1,1) /* filament          */
_(REAL,flm_g,1,		0,0,0,1,1) /*    colour         */
_(REAL,flm_b,0.8,	0,0,0,1,1) /*    components     */
_(REAL,flm_wt,3.0,	0,0,0,1,1) /*    and the linewidth */
_(INT,flm_balls,0,	0,0,0,1,1) /* "", plot it by spheres */
_(REAL,flm_em,0,	0,0,0,1,1) /* emission property */

/* T  N    I            R V L M D */
_(REAL,theta,0,		0,1,0,0,1) /* view angle: rotation about z axis */
_(REAL,phi,0,		0,1,0,0,1) /* view angle: elevation about z=0 plane */
_(REAL,psi,0,		0,1,0,0,1) /* view angle: rotation around view line */
_(REAL,distance,5.0,	0,1,0,0,1) /* former DISTANCE of ezgraph3d.h */

_(INT,remove_backs,0,	0,0,0,1,1) /* remove (positive/negative/none) part of surface */
_(INT,clipping,0,	0,0,0,0,1) /* use the extra clipping plane */
_(INT,ezdepth,1,	0,0,0,0,1) /* use the depth buffer */

_(INT,write_images,0,	0,0,0,0,0) /* save all images */
_(INT,write_filament,0,	0,0,0,0,0) /* write all filament data */

/* T  N    I            R V L M D */
_(REAL,bg_r,0.0,	0,0,0,0,1) /* background        */
_(REAL,bg_g,0.0,	0,0,0,0,1) /*    colour         */
_(REAL,bg_b,0.0,	0,0,0,0,1) /*    components     */
_(REAL,bbox_r,0.5,	0,0,0,0,1) /* bounding box      */
_(REAL,bbox_g,0.5,	0,0,0,0,1) /*    colour         */
_(REAL,bbox_b,0.5,	0,0,0,0,1) /*    components,    */
_(REAL,bbox_a,0.5,	0,0,0,0,1) /*    its transparency */
_(REAL,bbox_wt,2.0,	0,0,0,0,1) /*    and the linewidth */

_(INT,verbose,1,	0,0,0,0,0) /* verbosity level */

/* T  N    I            R V L M D */
_(REAL,lightx,-50.0,	0,0,1,0,1) /* position  */
_(REAL,lighty,+10.0,	0,0,1,0,1) /*   of the */
_(REAL,lightz,+10.0,	0,0,1,0,1) /*   directional */
_(REAL,lightw,1,	0,0,1,0,1) /*   light */
_(REAL,amb_r,0.5,	0,0,1,0,1) /* color */
_(REAL,amb_g,0.5,	0,0,1,0,1) /*   of the */
_(REAL,amb_b,0.5,	0,0,1,0,1) /*   ambient */
_(REAL,amb_a,0.5,	0,0,1,0,1) /*   light */

#if MARKER
/* T  N    I            R V L M D */
_(INT,show_marker,0,	0,0,0,1,1) /* whether to show the marker */
_(REAL,marker_size,0,	0,0,0,0,1) /* size and             */
_(REAL,marker_x,xmax/2-1,0,0,0,0,1)/*    coordinates       */
_(REAL,marker_y,ymax/2-1,0,0,0,0,1)/*    of the marker     */
_(REAL,marker_z,zmax/2-1,0,0,0,0,1)/*    position          */
_(REAL,marker_r,1.0,	0,0,0,0,1) /* marker               */
_(REAL,marker_g,1.0,	0,0,0,0,1) /*    colour            */
_(REAL,marker_b,0.0,	0,0,0,0,1) /*    components,       */
_(REAL,marker_a,1,	0,0,0,0,1) /*    its transparency  */
_(REAL,marker_wt,2.0,	0,0,0,0,1) /*    and the linewidth */
#endif

#if FIBRES
/* T  N    I            R V L M D */
_(INT,show_fibres,0,	0,0,0,1,1) /* whether to show the fibres */
_(INT,fib_dx,max(xmax/10,1),	0,0,0,1,1) /* # grid points between fibre roots */
_(INT,fib_dy,fib_dx,	0,0,0,1,1) /* .., along y axis */
_(INT,fib_dz,fib_dx,	0,0,0,1,1) /* .., along z axis */
_(INT,fib_x0,imin,	0,0,0,1,1) /* fibre root grid start along x */
_(INT,fib_y0,jmin,	0,0,0,1,1) /* fibre root grid start along y */
_(INT,fib_z0,kmin,	0,0,0,1,1) /* fibre root grid start along z */
_(REAL,fib_step,1.0,  	0,0,0,1,1) /* length of one segment */
_(INT,fib_len,10*xmax,  0,0,0,1,1) /* max num of segments */
_(INT,fib_rlayer,0,	0,0,0,1,1) /* red       channel in input files */
_(REAL,fib_rmin,0,	0,0,0,1,1) /* red       range lower end */
_(REAL,fib_rmax,1,	0,0,0,1,1) /* red       range upper end */
_(INT,fib_glayer,0,	0,0,0,1,1) /* green     channel in input files */
_(REAL,fib_gmin,0,	0,0,0,1,1) /* green     range lower end */
_(REAL,fib_gmax,1,	0,0,0,1,1) /* green     range upper end */
_(INT,fib_blayer,0,	0,0,0,1,1) /* blue      channel in input files */
_(REAL,fib_bmin,0,	0,0,0,1,1) /* blue      range lower end */
_(REAL,fib_bmax,1,	0,0,0,1,1) /* blue      range upper end */
_(REAL,fib_wt,1.0,	0,0,0,0,1) /* linewidth */
#endif
