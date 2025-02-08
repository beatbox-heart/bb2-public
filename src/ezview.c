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

/* Wrap around ezview utility version 1.7 based on Mantel and Barkley's EZScroll visualizer */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "k_.h"

#if defined(NOX11) || defined (NOGL)
  #include "nograph.h"
  NOGRAPH_DUMMY(ezview)
#else
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h> 
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <GL/gl.h>
#include <GL/glx.h>

#define OWN
#include "ezview.h" /* STR is defined in there! */


RUN_HEAD(ezview)
{
  int rc;
  #define _(Type,Name,Dflt,miN,maX) DEVICE_PAR(Type,Name);
  #include "ezfix.h"
  #undef _
  #define _(Type,Name,Init,reRead,reView,reLight,reMake,reDraw) DEVICE_PAR(Type,Name);
  #include "ezpar.h"
  #undef _
  
  glXMakeCurrent(S->theDisplay, S->theWindow, S->theGLXContext);
  Make_lists(S);
  /* setView (S,theta0,phi0,psi0,distance0); */
  setView (S,theta,phi,psi,distance);
  Draw(S);
  rc=Event_check(S);
  if (write_images) Save_image(S);
  if (write_filament) Save_filament(S);

  return(rc==0);
} RUN_TAIL(ezview);
/* ========================================================================= */

DESTROY_HEAD(ezview)
{
  QuitX(S);
  SAFE_CLOSE(S->filament);
} DESTROY_TAIL(ezview);
/* ========================================================================= */

CREATE_HEAD(ezview)
{
  Space s=dev->s;
  /* Initialize quasi-global computable variables that need it */
  #define _(Type,Name,Init) S->Name=Init;
  #include "ezini.h"
  #undef _

  /* Accept the user-defined parameters */
  #define noneREAL RNONE
  #define noneINT LNONE
  /* #define acceptREAL(b,c,d,e) if (!acceptrk(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); REAL b=S->b */
  #define acceptREAL(b,c,d,e) ACCEPTRK(b,c,d,e) REAL b=S->b
  /* #define acceptINT(b,c,d,e) if (!acceptik(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); INT b=S->b */
  #define acceptINT(b,c,d,e) ACCEPTIK(b,c,d,e) INT b=S->b
  #define _(Type,Name,Dflt,miN,maX) accept##Type(Name,Dflt,miN,maX);
  #include "ezfix.h"
  #undef _
  #define _(Type,Name,Init,reRead,reView,reLight,reMake,reDraw) accept##Type(Name,Init,none##Type,none##Type);
  #include "ezpar.h"
  #undef _
  #undef acceptreal
  #undef acceptint
  #undef noneint
  #undef nonereal

  DEVICE_IS_RECTANGULAR;
  
  /* Accept other paramerters */
  ACCEPTF(filament,"wt","");
  if (filament) {
    CALLOC(S->filament_buffer,FILAMENT_BUFFER_SIZE,1);
    S->filament_buffer_current = S->filament_buffer;
    S->filament_buffer_end = S->filament_buffer+FILAMENT_BUFFER_SIZE;
  } else {
    S->filament_buffer_end=S->filament_buffer_current=S->filament_buffer=NULL;
  }
  
  ACCEPTS(images,"");
  if (*images) {
    ACCEPTS(imagescode,"t");
    k_on();
    S->imagescompiled=compile(imagescode,deftb,t_real); CHK(imagescode);
    k_off();
  }
  
  ACCEPTS(title,"ezview t=%.0f");
  ACCEPTS(titlecode,"t");
  if (*titlecode) {
    S->titlecompiled=compile(S->titlecode,deftb,t_real); CHK(S->titlecode);
  } else {
    S->titlecompiled=NULL;
  }

  /* Some checks */
  ASSERT(0<imin); ASSERT(imin<imax); ASSERT(imax<xmax-1); /* marginal grid slices disallowed */
  ASSERT(0<jmin); ASSERT(jmin<jmax); ASSERT(jmax<ymax-1); /* because normal calculations */
  ASSERT(0<kmin); ASSERT(kmin<kmax); ASSERT(kmax<zmax-1); /* use central differences */

  switch (color_mode) {
  case 0: case 1: case 2: case 3:
    ASSERT(0<=vlayer); ASSERT(vlayer<vmax); ASSERT(maxv!=vmin);
    break;
  case 4: case 5: /* negative layer number signals colour is not used */
    ASSERT(-1<=rlayer); ASSERT(rlayer<vmax); ASSERT(rmax!=rmin);
    ASSERT(-1<=glayer); ASSERT(glayer<vmax); ASSERT(gmax!=gmin);
    ASSERT(-1<=blayer); ASSERT(blayer<vmax); ASSERT(bmax!=bmin);
    ASSERT(-1<=alayer); ASSERT(alayer<vmax); ASSERT(amax!=amin);
    break;
  }
  if (show_surface!=0 && amax<=0)
    MESSAGE("/* WARNING: show_surface=%d, alphamax=%lg, the surface will not be visible */\n"
	    ,(int)show_surface,(double)alphamax);
  if ( (show_filament!=0 || write_filament!=0) && layer1==layer2)
    MESSAGE("/* WARNING: show_filament=%d, write_filament=%d, but layer1=%d = layer2=%d, the filament will not be found */\n",
	    (int)show_filament,(int)write_filament,(int)layer1,(int)layer2);
  if (write_images && (0==(*images)))
    MESSAGE("/* WARNING: write_images=%d but not images template specified; no images will be written */\n",
	    (int)write_images);

  #if FIBRES
  if (ANISOTROPY_ON) {
    if (show_fibres) {
      ASSERT(fib_dx>0);
      ASSERT(fib_dy>0);
      ASSERT(fib_dz>0);
      ASSERT(-1<=fib_rlayer); ASSERT(fib_rlayer<vmax); ASSERT(fib_rmax!=fib_rmin);
      ASSERT(-1<=fib_glayer); ASSERT(fib_glayer<vmax); ASSERT(fib_gmax!=fib_gmin);
      ASSERT(-1<=fib_blayer); ASSERT(fib_blayer<vmax); ASSERT(fib_bmax!=fib_bmin);
    }
  } else {
    if (find_key("show_fibres=",w)) MESSAGE("/* Warning: anisotropy is on so show_fibres will be ignored */\n");
  }
  #endif
  
  /* Initialize drawing routines */
  if NOT(Draw_ini(S)) return FAILURE;
} CREATE_TAIL(ezview,0);
/* ========================================================================= */

#endif
