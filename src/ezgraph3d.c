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

/* ------------------------------------------------------------------------- */
/* ezgraph3d.c -- Graphics routines (except for marching cubes).
 *
 * Copyright (C) 1996 - 1998, 2006, 2007 Dwight Barkley and Matthew Dowle
 *
 * RCS Information
 * ---------------------------
 * $Revision: 1.5.1.1 $
 * $Date: 2007/05/07 10:07:03 $
 * ------------------------------------------------------------------------- */

/* Modifications for ezview: 2010-18, Vadim Biktashev v.n.biktashev@exeter.ac.uk */

/* #include <limits.h> */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "system.h"
#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "bikt.h"
#include "k_.h"
#include "pipe.h"

#if defined(NOX11) || defined (NOGL)
#else
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h> 
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "ezview.h"
#include "ezgraph3d.h"
#include "ezmarching.h"

/* -------------------------------------------------------------------------
 * DB: This file contains all graphics manipulation routines. The functions that
 * compute the iso-surfaces and filaments and render these are in ezmarching.
 *
 * The important things to know about this file are: 
 *
 * (1) X11 is used to open the graphics window and handle the events
 * (i.e. key presses and pointer motions within the window).  InitX() is a
 * long routine that opens the window; Event_check() contains the event loop
 * that looks for window events; QuitX() quits X obviously. You should be
 * able to switch over to a higher-level method with only minor
 * modifications except for these routines.
 *
 * (2) After the window is open, OpenGL is used to handle all the rendering.
 * myReshape() must be called before anything can be plotted.  To understand
 * this function see the OpenGL manual.
 *
 * The routines near the top of this file handle the interactive graphics
 * through OpenGL calls.  Note: to add a new interactive feature, add to the
 * event loop in Event_check() and add a corresponding routine to perform the
 * desired task.
 *
 * (3) Note on the graphics modes.
 *
 * There are 3 modes the program can be in: 
 *   MODE_SIMULATING  Simulation progresses through time (it cannot be 
 *                    rotated while in this mode)
 *   MODE_VIEWING     Simulation is paused.  While in this mode, graphics 
 *                    lists are created for use in rotations.  New list must
 *                    be created for each basic change in what is plotted,
 *                    e.g switching u-field to v-field, turning on clipping, 
 *                    etc.
 *   MODE_ROTATING    Mouse button 1 is down so the view responds to mouse 
 *                    motion.  Graphics lists are called and no new 
 *                    iso-surfaces are computed (i.e. fast rotations).
 * ------------------------------------------------------------------------- */

/* 
 * Private functions 
 * ----------------- */
static void  Pause        	      (STR *S);
static void  Restart                  (STR *S);
static void  Toggle_surface_drawing   (STR *S);
static void  Toggle_backs_removal     (STR *S);
static void  Toggle_clipping_plane    (STR *S);
static void  Toggle_depth_buffer      (STR *S);
static void  Toggle_filament_plotting (STR *S);
static void  Mouse_up                 (STR *S);
static void  Mouse_down               (STR *S, int x, int y);
/* Extension to Dowle-Barkley control */
static void  Rotate                   (STR *S, int x, int y, real new_theta, real new_phi);
static void  Tilt                     (STR *S,float new_psi);
/* "Iris Explorer" style control */
static void  Pitch		      (STR *S, float angle);
static void  Yaw                      (STR *S, float angle);
static void  Roll		      (STR *S, float angle);
static void  Turn                     (STR *S, int x, int y);
/* pseudo-Euler coords (turn to true Euler when find out what it is?) */
/* static void  setView                  (STR *S, float new_theta, float new_phi,float new_psi, float new_distance); */
static void  getView                  (float *current_theta, float *current_phi, float *current_psi);
static void  viewInit                 (STR *S);
				       
static void  Draw_title		      (STR *S);
static void  Draw_view                (STR *S);
static void  Proportions	      (STR *S);
static void  myReshape                (STR *S, int w, int h);
static void  MoveWindow               (STR *S, int dx, int dy);
static Bool  WaitForNotify            (Display *d, XEvent *e, char *arg);
static void  showbusy		      (STR *S, char c);
static void  Pop		      (STR *S);
static void  Tink		      (STR *S);

static void  Assign_par               (STR *S);
static void  Save_all                 (STR *S);

static int   InitX                    (STR *S);
/* ========================================================================= */

static char BUSY[]={BUSYRUNNING,BUSYWAITING,BUSYTURNING};

/* Highest level plotting routine.  Call glCallList() for surface
 * and filament precompiled by Marching_cubes via Make_lists()
 * function, and draw the simpler objects - box, marker, trace.
 */
void Draw (STR *S)
{
  DEVICE_CONST(int,show_surface);
  DEVICE_CONST(int,show_filament);
  DEVICE_CONST(real,theta);
  DEVICE_CONST(real,phi);
  DEVICE_CONST(real,psi);
  DEVICE_CONST(int,ezdepth);
  DEVICE_CONST(real,bg_r);
  DEVICE_CONST(real,bg_g);
  DEVICE_CONST(real,bg_b);
  DEVICE_CONST(real,bbox_wt);
  DEVICE_CONST(ezmode_type,ezmode);
  DEVICE_CONST(char *,busychar);
  DEVICE_CONST(int,NORM_LIST);
  DEVICE_CONST(int,ROT_LIST);
  DEVICE_CONST(int,FLM_LIST);
#if MARKER
  DEVICE_CONST(int,show_marker);
  DEVICE_CONST(int,marker_size);
#endif
#if FIBRES
  DEVICE_CONST(int,FIB_LIST);
  DEVICE_CONST(int,show_fibres);
#endif
  char oldbusy=*busychar;
  if (!GRAPHICS) return;
  showbusy(S,BUSYDRAWING);
  message(3,"Draw at t=%d view=(%g,%g,%g)\n",(int)t,theta,phi,psi);

  Proportions(S);			   /* Set current 3D window */
  Draw_title(S);			   /* Define the window title */
  
  glClearColor(bg_r,bg_g,bg_b,0.0); 	   /* Define the background */
  glClear(GL_COLOR_BUFFER_BIT | (ezdepth*GL_DEPTH_BUFFER_BIT));  /* Clear the color buffer.  
					    * Clear the depth buffer if on */

  /* First draw (almost) all the opaque objects */
#if MARKER
  /* Marker is always opaque */
  if (show_marker && marker_size) {
    Draw_marker(S);
    gdebug(3);
  }
#endif
#if FIBRES
  /* Fibres are always opaque */
  if (show_fibres && FIB_LIST) {
    glCallList(FIB_LIST);
    message(4,"FIB_LIST=%d called\n",FIB_LIST);    
    gdebug(3);
  }
#endif
  /* Filament is always opaque */
  if (show_filament && FLM_LIST) {
    glCallList(FLM_LIST);
    gdebug(3);
    message(4,"FLM_LIST=%d called\n",FLM_LIST);
  } 
  /* Bounding box may be semi-transparent */
  if (bbox_wt) {
    Draw_bounding_box(S);
    gdebug(3);
  }
  /* The surfaces are (potentially) transparent */
  if (show_surface) {
    switch (ezmode) {
    case MODE_SIMULATING :
    case MODE_VIEWING :
      if (NORM_LIST) glCallList(NORM_LIST);
      break;
    case MODE_ROTATING :
      if (ROT_LIST) glCallList(ROT_LIST);
      break;
    }
    gdebug(3);
  }
  glEnd();	/* to pair up with the most recent glBegin()s of any Draw_* function */

#if DOUBLEBUFFER
  glXSwapBuffers(S->theDisplay,S->theWindow);
#else
  glFinish();
#endif
  showbusy(S,oldbusy);
} /* Draw */
/* ========================================================================= */

/* Redraw the existing list from a different view. */
/* Always do it in rotating mode. */
void Draw_view (STR *S)
{
  DEVICE_VAR(ezmode_type,ezmode);
  ezmode_type ezmode_save=(*ezmode);
  (*ezmode)=MODE_ROTATING;
  Draw(S);
  (*ezmode)=ezmode_save;
} /* Draw_view */
/* ========================================================================= */

/* (Re)write the window title with flags and values */
/* static char name[256]=WINDOW_TITLE;  */
/* static char *pname=&(name[0]); */
/* static char *busychar; */

static void Draw_title (STR *S)
{
  DEVICE_ARRAY(char,title);
  DEVICE_CONST(pp_fn,titlecompiled);
  DEVICE_ARRAY(char,shorttitle);
  DEVICE_ARRAY(char,longtitle);
  DEVICE_CONST(int,show_surface);
  DEVICE_CONST(int,show_filament);
  #if FIBRES
  DEVICE_CONST(int,show_fibres);
  #endif
  #if MARKER
  DEVICE_CONST(int,show_marker);
  #endif
  DEVICE_CONST(real,theta);
  DEVICE_CONST(real,phi);
  DEVICE_CONST(real,psi);
  DEVICE_CONST(real,distance);
  DEVICE_CONST(int,remove_backs);
  DEVICE_CONST(int,clipping);
  DEVICE_CONST(int,ezdepth);
  DEVICE_CONST(Display *,theDisplay);
  DEVICE_CONST(Window,theWindow);
  DEVICE_CONST(GLXContext, theGLXContext);
  DEVICE_VAR(XTextProperty,theWindowName);
  DEVICE_VAR(XTextProperty, theIconName);
  DEVICE_CONST(ezmode_type,ezmode);
  DEVICE_VAR(char *,busychar);
  k_on();
  snprintf(shorttitle,MAXWINDOWTITLE,title,*(real *)execute(titlecompiled));
  k_off();  
  snprintf(longtitle,MAXWINDOWTITLE,"%s <%c> [%c%c"
	   #if FIBRES
	   "%c"
	   #endif
	   #if MARKER
	   "%c"
	   #endif
	   "%s%c%c] (%.1f,%.1f,%.1f,%.2f)",
	   shorttitle,
	   BUSY[ezmode],
	   show_surface?   'S':'s',
	   show_filament?  'F':'f',
	   #if FIBRES
	   show_fibres?    'B':'b',
	   #endif
	   #if MARKER
	   show_marker?    'M':'m',
	   #endif
	   remove_backs>0? "+V": remove_backs<0? "-V": " v",
	   clipping?       'C':'c',
	   ezdepth?        'D':'d',
	   theta, phi, psi, distance
	   );
   glXMakeCurrent(theDisplay, theWindow, theGLXContext);
   XStringListToTextProperty(&longtitle,1,theWindowName);
   XSetWMName(theDisplay, theWindow, theWindowName);
   XSetWMIconName(theDisplay, theWindow, theIconName);
   XFlush(theDisplay);
   *busychar=longtitle+strlen(shorttitle)+2;
} /* Draw_title */
/* ========================================================================= */

/* Create the GL lists of the surface (both resolutions) and filaments, */
/* so they can be plotted and replotted quickly afterwards */
void Make_lists (STR *S)
{
  DEVICE_CONST(int,norm_res);
  DEVICE_CONST(int,rot_res);
  DEVICE_VAR(int,NORM_LIST);
  DEVICE_VAR(int,ROT_LIST);
  DEVICE_VAR(int,FLM_LIST);
  DEVICE_CONST(char *,busychar);
  char oldbusy;
  oldbusy=*busychar;
  showbusy(S,BUSYDRAWING);

  /* Make new graphics list for both resolutions. */
  if ((*NORM_LIST)==0) (*NORM_LIST)=glGenLists(1);
  if ((*FLM_LIST)==0) (*FLM_LIST)=glGenLists(1);
  message(4,"FLM_LIST=%d allocated\n",(*FLM_LIST));

  /* Always make full resolution surface and filament */
  Marching_cubes(S,norm_res,(*NORM_LIST),(*FLM_LIST));

  /* Make rotation resolution surface if different */
  if (rot_res != norm_res) {
    if ((*ROT_LIST)==0) (*ROT_LIST)=glGenLists(1);
    Marching_cubes(S,rot_res,(*ROT_LIST),0);
  } else {
    (*ROT_LIST)=(*NORM_LIST);
  }

#if FIBRES  
  Make_fibres(S);
#endif
  showbusy(S,oldbusy);
} /* Make_lists */
/* ========================================================================= */

/* Return to MODE_SIMULATING. This has the effect of resuming the
   simulation via a return from Event_check(). */
static void Restart (STR *S)
{
  DEVICE_VAR(ezmode_type,ezmode);
  if ((*ezmode) == MODE_VIEWING) 
    (*ezmode) = MODE_SIMULATING;
  if ((*ezmode) != MODE_SIMULATING) 
    fprintf(stderr,"unexpected ezmode=%d in Restart(S)\n",(int)(*ezmode));
} /* Restart */
/* ========================================================================= */

/* The simulation has paused -- go to MODE_VIEWING. */
/* Call Make_lists() already done in the outer loop in ezview.c */
static void Pause (STR *S)
{
  DEVICE_VAR(ezmode_type,ezmode);

  if ((*ezmode) == MODE_SIMULATING)
    (*ezmode) = MODE_VIEWING;
  if ((*ezmode) != MODE_VIEWING)
    fprintf(stderr,"unexpected ezmode=%d in Pause(S)\n",(int)(*ezmode));
} /* Pause */
/* ========================================================================= */

static void Toggle_surface_drawing (STR *S)
{
  DEVICE_VAR(INT,show_surface);
  if (*show_surface) (*show_surface) = FALSE;
  else (*show_surface) = TRUE;
  Draw(S);
} /* Toggle_surface_drawing */
/* ========================================================================= */
  
/* This switch can be in three states: plot full surface (0), remove
   high-level half (1) and remove low-level half (-1) of it */
static void Toggle_backs_removal (STR *S)
{
  DEVICE_VAR(INT,remove_backs);
  (*remove_backs)++;
  if ((*remove_backs)>=2) (*remove_backs)=-1;
  Make_lists(S);
  Draw(S);
} /* Toggle_backs_removal */
/* ========================================================================= */

static void Toggle_clipping_plane (STR *S)
{
  DEVICE_VAR(INT,clipping);
  if (*clipping) {
    glDisable(GL_CLIP_PLANE0);
    (*clipping) = FALSE;
  } else {
    glEnable(GL_CLIP_PLANE0);
    (*clipping) = TRUE;
  }
  Make_lists(S); /* not sure why this is needed, but doesn't work without it */
  Draw(S);
} /* Toggle_clipping_plane */
/* ========================================================================= */

static void Toggle_depth_buffer (STR *S)
{
  DEVICE_VAR(INT,ezdepth);
  if (glIsEnabled(GL_DEPTH_TEST)) {
    glDisable(GL_DEPTH_TEST);
    *ezdepth = 0; 
  } else {
    glEnable(GL_DEPTH_TEST);
    *ezdepth = 1;
  }
  Draw(S);
} /* Toggle_depth_buffer */
/* ========================================================================= */

static void Toggle_filament_plotting (STR *S)
{
  DEVICE_VAR(INT,show_filament);
  if (*show_filament) (*show_filament) = FALSE;
  else (*show_filament) = TRUE;
  Draw(S);
} /* Toggle_filament_plotting */
/* ========================================================================= */

#if FIBRES
static void Toggle_fibres_plotting (STR *S)
{
  DEVICE_VAR(INT,show_fibres);
  if (*show_fibres) (*show_fibres) = FALSE;
  else (*show_fibres) = TRUE;
  Draw(S);
} /* Toggle_fibres_plotting */
/* ========================================================================= */
#endif

#if MARKER	
static void Toggle_marker_drawing (STR *S)
{
  DEVICE_VAR(INT,show_marker);
  if (*show_marker) (*show_marker) = FALSE;
  else (*show_marker) = TRUE;
  Draw(S);
} /* Toggle_marker_drawing */
/* ========================================================================= */
#endif

static void Mouse_up (STR *S)
{
  DEVICE_VAR(ezmode_type,ezmode);
  if ((*ezmode) == MODE_ROTATING)
    (*ezmode) = MODE_VIEWING;
  if ((*ezmode) != MODE_VIEWING)
    fprintf(stderr,"unexpected ezmode=%d in Mouse_up()\n",(int)(*ezmode));
  Draw(S);
} /* Mouse_up */
/* ========================================================================= */

static void Mouse_down (STR *S,int x, int y)
{
  DEVICE_CONST(int,norm_res);
  DEVICE_CONST(int,rot_res);
  DEVICE_CONST(int,width);
  DEVICE_CONST(int,height);
  DEVICE_VAR(ezmode_type,ezmode);
  DEVICE_VAR(int,lastx);
  DEVICE_VAR(int,lasty);
  DEVICE_ARRAY(GLfloat,mouse_down_mx);
  /* If the mode is MODE_VIEWING then on mouse down the state goes to
   * MODE_ROTATING.  x and y are the location of the mouse when the button
   * was pressed and these are saved as lastx and lasty for use in
   * rotating(). */

  if ((*ezmode) != MODE_VIEWING) return;

  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX, mouse_down_mx);

  (*lastx) = x;    /* Rotate() uses these values to calculate the distance  */ 
  (*lasty) = y;    /* moved by the mouse since the last Draw() call */

  (*ezmode) = MODE_ROTATING;
  if (rot_res!=norm_res) Draw(S);
  message(4,"Mouse down at (%d,%d) within (%d,%d)\n",(*lastx),(*lasty),width,height);
} /* Mouse_down */
/* ========================================================================= */

static void Rotate (STR *S, int x, int y, real new_theta, real new_phi)
{
  DEVICE_VAR(real,theta);
  DEVICE_VAR(real,phi);
  DEVICE_CONST(real,psi);
  DEVICE_CONST(real,distance);
  DEVICE_CONST(ezmode_type,ezmode);
  DEVICE_VAR(int,lastx);
  DEVICE_VAR(int,lasty);
  /* Rotates volume either in response to changes in mouse location or to
   * angles new_theta, new_phi.  More specifically:
   *
   * if MODE_ROTATING: change the viewing angles theta and phi by the
   *   difference of (x,y) - (lastx,lasty), i.e. the difference between the
   *   current mouse location and the previous location (set in call to
   *   Rotate() or to Mouse_down()).  
   *
   * else (not in MODE_ROTATING): set theta and phi to incoming values. */

  if (ezmode == MODE_ROTATING) {
    (*theta) += ROTATE_SCALE * (x-(*lastx));
    (*phi)   -= ROTATE_SCALE * (y-(*lasty));
    (*lastx) = x;
    (*lasty) = y;
  } else {
    (*theta) = new_theta;
    (*phi)   = new_phi;
  }

  /* Require both theta and phi to be in the interval [-180,180).  Extending
   * the range of theta to [-180,180) prevents 'flipping' at theta equal 90
   * and -90. */

  while ((*theta) > 180.0) (*theta) -= 360.0;
  while ((*theta) <= -180) (*theta) += 360.0;
  while ((*phi) > 180.0) (*phi) -= 360.0;
  while ((*phi) <= -180) (*phi) += 360.0;

  setView(S, (GLfloat)(*theta), (GLfloat)(*phi), (GLfloat)psi, (GLfloat)distance);

} /* Rotate */
/* ========================================================================= */


/* Turns the volume round y axis prior to rotating by theta and phi. */
static void Tilt (STR *S,float new_psi)
{
  DEVICE_CONST(real,theta);
  DEVICE_CONST(real,phi);
  DEVICE_VAR(real,psi);
  DEVICE_CONST(real,distance);
  (*psi)   = new_psi;

  /* Require psi to be in the interval [-180,180). */

  while ((*psi) >   180.0) (*psi)-=360.0;
  while ((*psi) <= -180.0) (*psi)+=360.0;

  setView(S, (GLfloat)theta, (GLfloat)phi, (GLfloat)(*psi), (GLfloat)distance);
} /* Tilt */
/* ========================================================================= */

static void Pitch (STR *S,float angle)
{
  DEVICE_VAR(real,theta);
  DEVICE_VAR(real,phi);
  DEVICE_VAR(real,psi);
  DEVICE_VAR(real,distance);
  GLfloat m[16];
  float current_theta, current_phi, current_psi;
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX, m);

  glLoadIdentity();
  glTranslatef(0.0,0.0,-(*distance)); 
  glRotatef(-angle,1,0,0);
  glTranslatef(0.0,0.0,+(*distance));

  glMultMatrixf(m);

  /* indirect, in case of types conflict */
  getView(&current_theta, &current_phi, &current_psi);
  (*theta)=current_theta; (*phi)=current_phi; (*psi)=current_psi;
  while ((*theta) > 180.0) (*theta) -= 360.0; while ((*theta) <= -180) (*theta) += 360.0;
  while ((*phi) > 180.0) (*phi) -= 360.0; while ((*phi) <= -180) (*phi) += 360.0;
  while ((*psi) >   180.0) (*psi)-=360.0; while ((*psi) <= -180.0) (*psi)+=360.0;
} /* Pitch */
/* ========================================================================= */

static void Yaw (STR *S,float angle)
{
  DEVICE_VAR(real,theta);
  DEVICE_VAR(real,phi);
  DEVICE_VAR(real,psi);
  DEVICE_CONST(real,distance);
  GLfloat m[16];
  float current_theta, current_phi, current_psi;
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX, m);

  glLoadIdentity();
  glTranslatef(0.0,0.0,-distance); 
  glRotatef(angle,0,1,0);
  glTranslatef(0.0,0.0,+distance); 

  glMultMatrixf(m);

  getView(&current_theta, &current_phi, &current_psi);
  (*theta)=current_theta; (*phi)=current_phi; (*psi)=current_psi;
  while ((*theta) > 180.0) (*theta) -= 360.0; while ((*theta) <= -180) (*theta) += 360.0;
  while ((*phi) >   180.0) (*phi)   -= 360.0; while ((*phi)   <= -180) (*phi)   += 360.0;
  while ((*psi) >   180.0) (*psi)   -= 360.0; while ((*psi)   <= -180) (*psi)   += 360.0;
} /* Yaw */
/* ========================================================================= */

static void Roll (STR *S,float angle)
{
  DEVICE_VAR(real,theta);
  DEVICE_VAR(real,phi);
  DEVICE_VAR(real,psi);
  DEVICE_CONST(real,distance);
  GLfloat m[16];
  float current_theta, current_phi, current_psi;
  glMatrixMode(GL_MODELVIEW);
  glGetFloatv(GL_MODELVIEW_MATRIX, m);

  glLoadIdentity();
  glTranslatef(0.0,0.0,-distance); 
  glRotatef(-angle,0,0,1);
  glTranslatef(0.0,0.0,+distance); 

  glMultMatrixf(m);

  getView(&current_theta, &current_phi, &current_psi);
  (*theta)=current_theta; (*phi)=current_phi; (*psi)=current_psi;
  while ((*theta) > 180.0) (*theta) -= 360.0; while ((*theta) <= -180) (*theta) += 360.0;
  while ((*phi) > 180.0) (*phi) -= 360.0; while ((*phi) <= -180) (*phi) += 360.0;
  while ((*psi) >   180.0) (*psi)-=360.0; while ((*psi) <= -180.0) (*psi)+=360.0;
} /* Roll */
/* ========================================================================= */

static void Turn (STR *S,int x, int y)
{
  DEVICE_VAR(real,theta);
  DEVICE_VAR(real,phi);
  DEVICE_VAR(real,psi);
  DEVICE_CONST(real,distance);
  DEVICE_VAR(ezmode_type,ezmode);
  DEVICE_CONST(int,lastx);
  DEVICE_CONST(int,lasty);
  DEVICE_ARRAY(GLfloat,mouse_down_mx);
  /* In MODE_ROTATING, turn the box round in response to changes in mouse location.
   *
   * The model: the inscribed disk of radius r
   * is considered as 2D projection of a sphere,
   * that represents the SO3 coords that define the rotation, 
   * corresp. to new vs old mouse position. 
   */
  float new_theta, new_phi, new_psi;
  /* Radius and centre of the disk */
  float r=0.5*WINSIZE, xc=0.5*WINSIZE, yc=0.5*WINSIZE;
  float X0, Y0, Z0, P0, X1, Y1, Z1, P1, Xa, Ya, Za, crossnorm, dotprod, angle;

  if ((*ezmode) != MODE_ROTATING) return; 

  message(4,"(%d,%d) -> (%d,%d)\n",lastx,lasty,x,y);

  message(4,"xc=%g yc=%g r=%g\n",xc,yc,r);

  /* Old mouse position on the "sphere", normalized and cropped */
  /* NB mouse position is in "raster coord" which is upside down so Y needs to be negated */
  X0=(lastx-xc)/r; Y0=-(lasty-yc)/r;
  P0=X0*X0+Y0*Y0;
  message(4,"Old: x=%d X0=%g y=%d Y0=%g P0=%g\n",lastx,X0,lasty,Y0,P0);
  if (P0>1.0) {
    X0/=sqrt(P0);
    Y0/=sqrt(P0);
    P0=1.0;
  }
  Z0=sqrt(1.0-P0);

  /* New position */
  /* NB mouse position is in "raster coord" which is upside down so Y needs to be negated */
  X1=(x-xc)/r; Y1=-(y-yc)/r;
  P1=X1*X1+Y1*Y1;
  message(4,"New: x=%d X1=%g y=%d Y1=%g P1=%g\n",x,X1,y,Y1,P1);
  if (P1>1.0) {
    X1/=sqrt(P1);
    Y1/=sqrt(P1);
    P1=1.0;
  }
  Z1=sqrt(1.0-P1);

  /* Cross-product */
  Xa=Y0*Z1-Y1*Z0;
  Ya=Z0*X1-Z1*X0;
  Za=X0*Y1-X1*Y0;

  /* Dot-product */
  dotprod=X0*X1+Y0*Y1+Z0*Z1;

  message(4,"(%g,%g,%g) x (%g,%g,%g) = (%g,%g,%g) . (%g)",X0,Y0,Z0,X1,Y1,Z1,Xa,Ya,Za,dotprod);

  /* Its length is sin of rotation angle */
  crossnorm=sqrt(Xa*Xa+Ya*Ya+Za*Za);
  angle=asin(crossnorm)/deg;
  if (dotprod<0) angle=180.0-angle;
  /* Normalized vector as rotation axis */
  Xa/=crossnorm;
  Ya/=crossnorm;
  Za/=crossnorm;

  message(4,"-> %g about (%g,%g,%g)\n",angle,Xa,Ya,Za);

  /* Now we do the rotation by prepending appropriate trans mx */
  /* to the one recorded at mouse-down time */
  glLoadIdentity();
  glTranslatef(0.0,0.0,-distance); 
  glRotatef(angle,Xa,Ya,Za);
  glTranslatef(0.0,0.0,+distance); 
  glMultMatrixf(mouse_down_mx);

  /* The resulting SO3 coords */
  getView(&new_theta, &new_phi, &new_psi);
  (*theta)=new_theta; (*phi)=new_phi; (*psi)=new_psi;
  while ((*theta) > 180.0) (*theta) -= 360.0; while ((*theta) <= -180) (*theta) += 360.0;
  while ((*phi) > 180.0) (*phi) -= 360.0; while ((*phi) <= -180) (*phi) += 360.0;
  while ((*psi) >   180.0) (*psi)-=360.0; while ((*psi) <= -180.0) (*psi)+=360.0;
} /* Turn */
/* ========================================================================= */


void setView (STR *S,float new_theta, float new_phi, float new_psi, float new_distance)
{
  DEVICE_VAR(real,theta);
  DEVICE_VAR(real,phi);
  DEVICE_VAR(real,psi);
  DEVICE_VAR(real,distance);
  DEVICE_CONST(real,lightx);
  DEVICE_CONST(real,lighty);
  DEVICE_CONST(real,lightz);
  DEVICE_CONST(real,lightw);
  DEVICE_CONST(real,amb_r);
  DEVICE_CONST(real,amb_g);
  DEVICE_CONST(real,amb_b);
  DEVICE_CONST(real,amb_a);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_ARRAY(real,plot_ini);
  DEVICE_ARRAY(GLfloat,light_position0);
  DEVICE_ARRAY(GLfloat,lmodel_ambient);

  /* DB: OpenGL calls to set rendering (transformation matrix). Note that the
   * last call (which is the first action applied to a triangle being
   * rendered) is to center the 3D volume about the origin. There are then
   * rotations in (VNB: psi,) theta and phi.  The second translation pulls the view away
   * from the cube.  I have found that distance=5 is about right -- the
   * distance (as measured in unit lengths as they appear on the screen) the
   * volume is from the eye under normal view conditions is about 5.  In
   * principle this varies with the size of the window on the screen.  For
   * an orthographic projection this is irrelevant so long as the cube is in
   * the viewed volume. */

  (*theta)=new_theta; (*phi)=new_phi; (*psi)=new_psi; (*distance)=new_distance;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0,0.0,-new_distance); 
  glRotatef(-new_phi,1.0,0.0,0.0);
  glRotatef(new_theta,0.0,0.0,1.0);
  glRotatef(new_psi,0.0,1.0,0.0);
  glTranslatef(
    -(plot_ini[0]+plot_length[0]/2),
    -(plot_ini[1]+plot_length[1]/2),
    -(plot_ini[2]+plot_length[2]/2)
  );

  /* VNB: I want position of the light to be fixed in eye coordinates,
     rather than in the model coordinates */
  glPushMatrix();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  light_position0[0]=lightx;
  light_position0[1]=lighty;
  light_position0[2]=lightz;
  light_position0[3]=lightw;
  glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
  glPopMatrix();

  /* .. and (re)define the global ambient light in case it changed */
  lmodel_ambient[0]=amb_r;
  lmodel_ambient[1]=amb_g;
  lmodel_ambient[2]=amb_b;
  lmodel_ambient[3]=amb_a;
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
} /* setView */
/* ========================================================================= */

/* Reconstruct the current SO3 coord from the current transformation matrix */
/* Arbitrarily restric theta to [-pi/2,pi/2] range */
static void getView (float *current_theta, float *current_phi, float *current_psi)
{
  GLfloat m[16];
  double eulertheta=0, costh=0; /* initial assignments to keep compiler happy */
  double eulerphi;
  double eulerpsi;
  double R11, R12, R13, R21, R22, R23, R31, R32, R33;
  glGetFloatv(GL_MODELVIEW_MATRIX, m);
  #define r(i,j) (double)(m[(j)-1+((i)-1)*4])
  R11=r(1,1); R12=r(1,3); R13=r(1,2);
  R21=r(3,1); R22=r(3,3); R23=r(3,2);
  R31=r(2,1); R32=r(2,3); R33=r(2,2);
  /* http://gregslabaugh.name/publications/euler.pdf with psi<->phi */
  if (R31==1.0) {
    /* sinth=-1 costh=0 */
    /* (0          , -cosph*sinps-sinph*cosps,  sinph*sinps-cosph*cosps); */
    /* (0          ,  cosph*cosps-sinph*sinps, -sinph*cosps-cosph*sinps); */
    /* (1          ,               0         ,        0                ); */
    /**/
    /* (0          , -sin(phi+psi),  -cos(phi+psi) */
    /* (0          ,  cos(phi+psi),  -sin(phi+psi); */
    /* (1          ,  0         ,        0       ); */
    eulerpsi=0.0;
    eulerphi=atan2(-R12,-R13);
  } else if (R31==-1.0) {
    /* sinth=1 costh=0 */
    /* (0       , -cosph*sinps+sinph*cosps,  sinph*sinps+cosph*cosps); */
    /* (0       ,  cosph*cosps+sinph*sinps, -sinph*cosps+cosph*sinps); */
    /* (-1      ,              0          ,        0                ); */
    /**/
    /* (0       ,  sin(phi-psi),  cos(phi-psi); */
    /* (0       ,  cos(phi-psi),  sin(psi-phi); */
    /* (-1      ,       0      ,      0      ); */
    eulerpsi=0.0;
    eulerphi=atan2(R12,R13);
  } else {
    /* (costh*cosps, -cosph*sinps+sinph*sinth*cosps,  sinph*sinps+cosph*sinth*cosps); */
    /* (costh*sinps,  cosph*cosps+sinph*sinth*sinps, -sinph*cosps+cosph*sinth*sinps); */
    /* (-sinth     ,              sinph*costh      ,        cosph*costh            ); */
    eulertheta=-asin(R31);
    costh=cos(eulertheta);
    eulerphi=atan2(R32/costh,R33/costh);
    eulerpsi=atan2(R21/costh,R11/costh);
  }
  *current_theta=eulertheta/deg;
  *current_phi=-eulerphi/deg;
  *current_psi=eulerpsi/deg;
} /* getView */
/* ========================================================================= */

/* Current camera will be the "initial" */
static void  viewInit (STR *S)
{
  if (S->theta0==S->theta && S->phi0==S->phi && S->psi0==S->psi && S->distance0==S->distance) return;
  S->theta0=S->theta;
  S->phi0=S->phi;
  S->psi0=S->psi;
  S->distance0=S->distance;
  MESSAGE("\x01/* ezview idev=%d new initial camera: theta="REALF" phi="REALF" psi="REALF" distance="REALF" */\n",
	  idev,S->theta0,S->phi0,S->psi0,S->distance0);
}

/* ========================================================================= */

static void myReshape (STR *S, int w, int h)
{
  DEVICE_CONST(real,distance);
  DEVICE_VAR(INT,width);
  DEVICE_VAR(INT,height);

  GLfloat half_width;

  /* DB: half_width is the half width of the volume viewed in GL.  In general if
   * half_width is large, then rendered simulation volume will appear small,
   * and vice versa.  PLOT_SIZE in ezgraph3d.h allows adjustment of this
   * without changing any of the code below.
   *
   * For orthographic projection setting half_width this is straightforward:
   * The simulation volume has maximum side of 1 in graphics coordinates
   * (set Draw_ini()).  If half_width is set to sqrt(3)/2 = 0.866..., then a
   * unit cube, and hence the simulation volume, will always lie entirely
   * within the viewport so I choose this for the default.  For ortho
   * projection, the near are far values are set to large values and are
   * irrelevant so long as the volume lies between near and far.
   * 
   * For perspective projection the situation is more complicated because
   * the half_width applies at the value of near (see the OpenGL manual).
   * Hence the value of near and half_width must be set together so that a
   * unit cube fills the viewport but does not cut the near plane.  The
   * values below accomplish this.  The far value if irrelevant and is set to
   * a large value. */

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
#if PERSPECTIVE   
  /* Perspective projection */
  half_width = 0.7/PLOT_SIZE;
  glFrustum(-half_width, half_width, -half_width, half_width, 
	    distance-1., 20.);
#else            
  /* Orthographic projection */
  half_width = 0.866/PLOT_SIZE;
  glOrtho (-half_width, half_width, -half_width, half_width, -20., 20.); 
#endif
  glMatrixMode (GL_MODELVIEW);
  glViewport (0, 0, w, h);

  (*width)=w; (*height)=h;
} /* myReshape */
/* ========================================================================= */

/* Initialize everything necessary for 3D plotting.  */
int Draw_ini (STR *S)
{
  DEVICE_VAR(INT,show_surface);
  DEVICE_CONST(real,theta0);
  DEVICE_CONST(real,phi0);
  DEVICE_CONST(real,psi0);
  DEVICE_CONST(real,distance0);
  DEVICE_VAR(INT,clipping);
  DEVICE_VAR(INT,ezdepth);
  DEVICE_CONST(real,bg_r);
  DEVICE_CONST(real,bg_g);
  DEVICE_CONST(real,bg_b);
  DEVICE_CONST(int,autostart);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_VAR(INT,winx);
  DEVICE_VAR(INT,winy);
  DEVICE_VAR(ezmode_type,ezmode);

  /* DB: The lengths of the simulation volume in graphics coordinates are set.  I
   * choose to have the largest plot_length=1.  Thus the volume lies inside
   * the unit cube in graphics coordinates. */
  /* VNB: here is the initial setup, to allow the clipping plane; redone every Draw */
  Proportions(S);

  /* Initialization for marching cubes. */
  Marching_ini(S);

  /* DB: At this point everything has been initialized for finding
   * filaments without graphics.  Can return after setting a few
   * things.  Clipping must be set but is not important. */

  if( !GRAPHICS ) {
    (*show_surface) = FALSE;
    (*ezmode) = MODE_SIMULATING;
    (*clipping) = FALSE;
    showbusy(S,BUSYRUNNING);
    return SUCCESS;
  }

  /* Create X window and prepares for use by OpenGL */
  if NOT(InitX(S)) return FAILURE;

  /* The default should be smooth anyway */
  glShadeModel(GL_SMOOTH); 

  /* Enable the depth buffer if required. 
   * Standard value of the flag is 1. */
  if (*ezdepth!=0) {
    glEnable(GL_DEPTH_TEST);
    *ezdepth = 1;
  } 

  /* Enable light sources and properties */
  /* glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); - checked, not good!  */

  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  mat_emission);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffusive);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  /* glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1); */
  /* glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0); */
  /* glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 1.e-3); */

  glEnable(GL_COLOR_MATERIAL); 
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

  glEnable (GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);

  /* Set the background color.  */
  glClearColor(bg_r,bg_g,bg_b,0.0); 

  /* DB: Set up the clipping plane.  Note that clipping is not on enabled
   * initially, but could be by inserting a call to Toggle_clipping_plane().
   * Here the clipping plane is set to pass half way through the volume.
   * plane[3] controls this.  I start with clipping off (clipping=FALSE) */
  {
    GLdouble plane[4]; 
    plane[0] = 0.0; 
    plane[1] = 0.0;
    plane[2] = -1.0;
    plane[3] = plot_length[2]/2;
    setView(S, 0., 0., 0., distance0);
    glClipPlane(GL_CLIP_PLANE0, plane);
  }

  /* Initialize the view direction. */
  /* The lights are reset in there to maintain their position in the
     eye coordinate system. */
  setView(S, theta0,phi0,psi0,distance0);
  /*------------------------------------------------*/ printf("Draw_ini: setView(%g,%g,%g,%g)\n",theta0,phi0,psi0,distance0);
 
  /* Form the title bar of the window. */
  /* This is actually redone with every Draw() ?? */
  Draw_title(S);

  /* Start paused unless autostart flag is on */
  if (autostart) {
    (*ezmode) = MODE_SIMULATING;
    showbusy(S,BUSYRUNNING);
  } else {
    (*ezmode) = MODE_VIEWING;
    showbusy(S,BUSYWAITING);
  }

  return SUCCESS;
} /* Draw_ini */
/* ========================================================================= */

/* Interactively reassign any of the parameters listed in ezpar.h, */
/* in the parent command shell (TODO: proper dialog pop-up window?) */
enum par_t { /* integer and real for now; may be add strings one day */
  par_INT,
  par_REAL
};
enum par_count {
#define _(t,n,i,r,v,l,m,d) par_##n,
  #include "ezpar.h"
  numpar
#undef _
};
typedef struct {
  char *name; 
  void *addr; 
  enum par_t type;
  int reread, review, relight, remake, redraw;
} pardef;
static void  Assign_par (STR *S)
{
  DEVICE_CONST(real,theta);
  DEVICE_CONST(real,phi);
  DEVICE_CONST(real,psi);
  DEVICE_CONST(real,distance);
  DEVICE_CONST(Display *,theDisplay);
  DEVICE_CONST(Window,theWindow);

  /* This is on the stack so no conflict between device instances - is that important? */
  pardef pars[] = {
  #define _(t,n,i,r,v,l,m,d) {#n,(void *)&(S->n),par_##t,r,v,l,m,d},
    #include "ezpar.h"
    {"",NULL,par_INT,0,0,0,0,0}
  #undef _
  };
  #define SEP " \t\r\b\n="
  char buf[BUFLEN], *name, *value, *lasts;
  int ip, np;
  int changed;
  INT oldint, newint;
  REAL oldreal, newreal;

  while (1) {
    printf("parameter=");
    fgets(buf,BUFLEN,stdin);
    name=strtok_r(buf,SEP,&lasts);
    if (!name) break;
    np=-1;
    for(ip=0;ip<numpar;ip++) if(strcasecmp(buf,pars[ip].name)==0) np=ip;
    if (np==-1) {
      printf("Known parameters:");
      for (ip=0;ip<numpar;ip++) printf("%s%c",pars[ip].name,ip<numpar-1?',':'\n');
      continue;
    }
    switch (pars[np].type) {
    case par_INT:
      oldint=*(INT *)pars[np].addr;
      printf("old value="INTFMT" new value=",oldint);
      fgets(buf,BUFLEN,stdin);
      value=strtok_r(buf,SEP,&lasts);
      if (!value) continue;
      sscanf(value,INTFMT,&newint);
      changed=(newint!=oldint);
      *(INT *)pars[np].addr=newint;
      printf("%s="INTFMT"\n",pars[np].name,*(INT *)(pars[np].addr));
      break;
    case par_REAL:
      oldreal=*(REAL *)pars[np].addr;
      printf("old value="FFMT" new value=",oldreal);
      fgets(buf,BUFLEN,stdin);
      value=strtok_r(buf,SEP,&lasts);
      if (!value) continue;
      sscanf(value,FFMT,&newreal);
      changed=(newreal!=oldreal);
      *(REAL *)(pars[np].addr)=newreal;
      printf("%s="FFMT"\n",pars[np].name,*(REAL *)(pars[np].addr));
      break;
    }
    if (changed) {
      if (pars[np].review) {
	message(2,"need to setView due to a changed parameter\n");
	setView(S,theta,phi,psi,distance);  /* some parameters affect camera position */
      }
      if (pars[np].relight) {
	message(2,"re-lighting due to a changed parameter\n");
	setView(S,theta,phi,psi,distance);  /* light position is fixed in eye coords */
      }
      if (pars[np].remake) {
	message(2,"need to re-make lists due to a changed parameter\n");
	Make_lists(S);
      }
      if (pars[np].redraw) {
	message(2,"need to re-draw due to a changed parameter\n");
	Draw(S);
      }
    } /* if changed */
  } /* while(1) */
  XRaiseWindow(theDisplay,theWindow);
  XSetInputFocus(theDisplay, theWindow, RevertToParent, CurrentTime); 
} /* Assign_par */
/* ========================================================================= */

void Save_all (STR *S)
{
  Draw(S);
  Save_image(S);
  Save_filament(S);
} /* Save_all */
/* ========================================================================= */

void MoveWindow (STR *S, int dx, int dy)
{
  DEVICE_CONST(Display *,theDisplay);
  DEVICE_CONST(Window,theWindow);
  DEVICE_VAR(INT,winx);
  DEVICE_VAR(INT,winy);
  DEVICE_CONST(int,horspace);
  DEVICE_CONST(int,vertspace);
  (*winx) += dx;
  if (horspace) {
    while ((*winx)<0) (*winx)+=horspace;
    while ((*winx)>=horspace) (*winx)-=horspace;
  }
  (*winy) += dy;
  if (vertspace) {
    while ((*winy)<0) (*winy)+=vertspace;	
    while ((*winy)>=vertspace) (*winy)-=vertspace;
  }
  XMoveWindow(theDisplay,theWindow,(*winx),(*winy)); 
} /* MoveWindow */
/* ========================================================================= */

int Event_check (STR *S)
{
  DEVICE_VAR(INT,color_mode); 
  DEVICE_VAR(real,theta);       	/* covertly changed */
  DEVICE_VAR(real,phi);         	/*   by Turn        */
  DEVICE_VAR(real,psi);         	/*   and Tilt       */
  DEVICE_VAR(real,distance);
  DEVICE_CONST(real,theta0);
  DEVICE_CONST(real,phi0);
  DEVICE_CONST(real,psi0);
  DEVICE_CONST(real,distance0);
  DEVICE_CONST(Display *,theDisplay);
  DEVICE_CONST(Window,theWindow);
  DEVICE_VAR(ezmode_type,ezmode); 	/* covertly changeable by Restart() */
  #if MARKER
  DEVICE_CONST(int,imin);
  DEVICE_CONST(int,imax);
  DEVICE_CONST(int,jmin);
  DEVICE_CONST(int,jmax);
  DEVICE_CONST(int,kmin);
  DEVICE_CONST(int,kmax);
  DEVICE_CONST(real,marker_size);
  DEVICE_VAR(real,marker_x);
  DEVICE_VAR(real,marker_y);
  DEVICE_VAR(real,marker_z);
  #endif

  XEvent         theEvent;
  KeySym         theKeySym;
  #define theKeyBufferMaxLen 64
  char           theKeyBuffer[theKeyBufferMaxLen+1];
  XComposeStatus theComposeStatus;
  /* Standard XLib coding for these keys */
  #define SHIFT ShiftMask
  #define CTRL ControlMask
  /* Could not find out proper XLib coding for this 
     so found it experimentally, may not work on other platforms */
  /* #define ALT Mod1Mask */
  #define ALT 0x2000
  unsigned int   state;

  if ( !GRAPHICS ) return(0);

  /* VB: remember camera position, to record any difference afterwards */
  real old_theta=*theta;
  real old_phi=*phi;
  real old_psi=*psi;
  real old_distance=*distance;

  int stop_next_time=0;

  /* X Event Loop
   *
   * If there is something in the queue then each event is processed in
   * turn. When the queue is empty, and the mode (which may have been changed
   * by the events just processed) is MODE_SIMULATING then control is
   * returned (presumably to the main loop in main()).  However, when the
   * queue is empty and the mode is either MODE_ROTATING or MODE_VIEWING then
   * XNextEvent is called which blocks until the next X event is received
   * (eg. the space bar being pressed which sets the mode back to
   * MODE_SIMULATING and so control returns to main()). */

  if (*ezmode != MODE_SIMULATING) showbusy(S,BUSYWAITING);

  while (XPending(theDisplay) || ((*ezmode) != MODE_SIMULATING)) {

    XNextEvent(theDisplay, &theEvent);
    
    switch(theEvent.type) {      /* switch according to X event */
      
    case KeyPress:               /* A KeyPress event has occurred. */
      
      XLookupString((XKeyEvent *)&theEvent, theKeyBuffer,
		    theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
      state= (theEvent.xkey.state) & ( CTRL | SHIFT | ALT );
      message(4,"state=%04x -> %04x = %s:%s:%s\n",
			     theEvent.xkey.state,
			     state,
			     state&CTRL?"C":"-",
			     state&ALT?"A":"-",
			     state&SHIFT?"S":"-");


      switch(theKeySym) {        /* switch according to the pressed key */

      /* case XK_Escape: exit(0);   /\* hard exit, nothing saved *\/ */
      case XK_P: case XK_p: Pause(S); showbusy(S,BUSYWAITING); break;
      case XK_Q: case XK_q: return(1);
      case XK_N: case XK_n: stop_next_time=1; Restart(S); showbusy(S,BUSYRUNNING); break;
      case XK_space:  Restart(S); showbusy(S,BUSYRUNNING); break;

	/*******************************/
#if FIBRES
      case XK_B: case XK_b: Toggle_fibres_plotting(S); break;
#endif	
      case XK_C: case XK_c: Toggle_clipping_plane(S); break;
      case XK_D: case XK_d: Toggle_depth_buffer(S); break;
      case XK_F: case XK_f: Toggle_filament_plotting(S); break;
#if MARKER	
      case XK_M: case XK_m: Toggle_marker_drawing(S); break;
#endif
      case XK_O: case XK_o:
	switch(state) {
	case (0): setView(S, theta0,phi0,psi0,distance0); Draw(S); break;
	case (CTRL): viewInit(S); break;
	default:           Pop(S); break;
	}
	break;
      case XK_R: case XK_r: Make_lists(S); Draw(S); break;
      case XK_S: case XK_s: Toggle_surface_drawing(S); break;
      case XK_V: case XK_v: Toggle_backs_removal(S); break;
      case XK_W: case XK_w: Save_all(S); break;
      case XK_X: setView(S, -90,-90,180,(*distance)); Draw_view(S); break;
      case XK_x: setView(S,  90,  0, 90,(*distance)); Draw_view(S); break;
	
      case XK_Y: setView(S,   0,-90,180,(*distance)); Draw_view(S); break;
      case XK_y: setView(S,   0, 90,  0,(*distance)); Draw_view(S); break;
      case XK_Z: setView(S,   0,  0,  0,(*distance)); Draw_view(S); break;
      case XK_z: setView(S,   0,  0,180,(*distance)); Draw_view(S); break;

      case XK_0: (*color_mode)=0; Make_lists(S); Draw(S); break;
      case XK_1: (*color_mode)=1; Make_lists(S); Draw(S); break;
      case XK_2: (*color_mode)=2; Make_lists(S); Draw(S); break;
      case XK_3: (*color_mode)=3; Make_lists(S); Draw(S); break;
      case XK_4: (*color_mode)=4; Make_lists(S); Draw(S); break;
      case XK_5: (*color_mode)=5; Make_lists(S); Draw(S); break;
      case XK_6: (*color_mode)=6; Make_lists(S); Draw(S); break;
      case XK_7: (*color_mode)=7; Make_lists(S); Draw(S); break;
      case XK_8: (*color_mode)=8; Make_lists(S); Draw(S); break;
      case XK_9: (*color_mode)=9; Make_lists(S); Draw(S); break;

      case XK_Left:  
	switch(state) {
	case (0):	   Rotate(S,0,0,(*theta)-ROTATE_STEP,(*phi)); Draw_view(S); break;
	case (SHIFT):      Rotate(S,0,0,(*theta)-many*ROTATE_STEP,(*phi)); Draw_view(S); break;
	case (ALT):	   Yaw(S,-ROTATE_STEP); Draw_view(S); break;
	case (ALT+SHIFT):  Yaw(S,-many*ROTATE_STEP); Draw_view(S); break;
	case (CTRL):       MoveWindow(S,-1,0); break;
	case (CTRL+SHIFT): MoveWindow(S,-many,0); break;
#if MARKER	
	case (CTRL+ALT): if (marker_size>0 && (*marker_x)>imin) {(*marker_x)--; Draw(S);} else Tink(S); break;
	case (CTRL+ALT+SHIFT): if (marker_size>0 && (*marker_x)>=imin+many) {(*marker_x)-=many; Draw(S);} else Tink(S); break;
#endif
	default:           Pop(S); break;
	}
	break;
      case XK_Right: 
	switch(state) {
	case (0):          Rotate(S,0,0,(*theta)+ROTATE_STEP,(*phi)); Draw_view(S); break;
	case (SHIFT):      Rotate(S,0,0,(*theta)+many*ROTATE_STEP,(*phi)); Draw_view(S); break;
	case (ALT):        Yaw(S,+ROTATE_STEP); Draw_view(S); break;
	case (ALT+SHIFT):  Yaw(S,+many*ROTATE_STEP); Draw_view(S); break;
	case (CTRL):       MoveWindow(S,+1,0); break;
	case (CTRL+SHIFT): MoveWindow(S,+many,0); break;
#if MARKER	
	case (CTRL+ALT): if (marker_size>0 && (*marker_x)<imax) {(*marker_x)++; Draw(S);} else Tink(S); break;
	case (CTRL+ALT+SHIFT): if (marker_size>0 && (*marker_x)<=imax-many) {(*marker_x)+=many; Draw(S);} else Tink(S); break;
#endif
	default:           Pop(S); break;
	}
	break;
      case XK_Up: 
	switch(state) {
	case (0):          Rotate(S,0,0,(*theta),(*phi)+ROTATE_STEP); Draw_view(S); break;
	case (SHIFT):      Rotate(S,0,0,(*theta),(*phi)+many*ROTATE_STEP); Draw_view(S); break;
	case (ALT):        Pitch(S,+ROTATE_STEP); Draw_view(S); break;
	case (ALT+SHIFT):  Pitch(S,+many*ROTATE_STEP); Draw_view(S); break;
	case (CTRL):       MoveWindow(S,0,-1); break;
	case (CTRL+SHIFT): MoveWindow(S,0,-many); break;
#if MARKER	
	case (CTRL+ALT): if (marker_size>0 && (*marker_y)<jmax) {(*marker_y)++; Draw(S);} else Tink(S); break;
	case (CTRL+ALT+SHIFT): if (marker_size>0 && (*marker_y)<=jmax-many) {(*marker_y)+=many; Draw(S);} else Tink(S); break;
#endif
	default:           Pop(S); break;
	}
	break;
      case XK_Down:
	switch(state) {
	case (0):          Rotate(S,0,0,(*theta),(*phi)-ROTATE_STEP); Draw_view(S); break;
	case (SHIFT):      Rotate(S,0,0,(*theta),(*phi)-many*ROTATE_STEP); Draw_view(S); break;
	case (ALT):        Pitch(S,-ROTATE_STEP); Draw_view(S); break;
	case (ALT+SHIFT):  Pitch(S,-many*ROTATE_STEP); Draw_view(S); break;
	case (CTRL):       MoveWindow(S,0,+1); break;
	case (CTRL+SHIFT): MoveWindow(S,0,+many); break;
#if MARKER	
	case (CTRL+ALT): if (marker_size>0 && (*marker_y)>jmin) {(*marker_y)--; Draw(S);} else Tink(S); break;
	case (CTRL+ALT+SHIFT): if (marker_size>0 && (*marker_y)>=jmin+many) {(*marker_y)-=many; Draw(S);} else Tink(S); break;
#endif
	default:           Pop(S); break;
	}
	break;
      case XK_Page_Down: 
	switch(state) {
	case (0):          Tilt(S,(*psi)+ROTATE_STEP); Draw_view(S); break;
	case (SHIFT):      Tilt(S,(*psi)+many*ROTATE_STEP); Draw_view(S); break;
	case (ALT):        Roll(S,+ROTATE_STEP); Draw_view(S); break;
	case (ALT+SHIFT):  Roll(S,+many*ROTATE_STEP); Draw_view(S); break;
	/* case (CTRL):       m+=mstep; if (Read_fmt(0,m)) {Make_lists(S); Draw(S);} else {m-=mstep;} break; */
	/* case (CTRL+SHIFT): if (m!=mmax) {m=mmax; Read_fmt(0,m);}; Make_lists(S); Draw(S); break; */
#if MARKER	
	case (CTRL+ALT): if (marker_size>0 && (*marker_z)>kmin) {(*marker_z)--; Draw(S);} else Tink(S); break;
	case (CTRL+ALT+SHIFT): if (marker_size>0 && (*marker_z)>=kmin+many) {(*marker_z)-=many; Draw(S);} else Tink(S); break;
#endif
	default:           Pop(S); break;
	}
	break;
      case XK_Page_Up: 
	switch(state) {
	case (0):          Tilt(S,(*psi)-ROTATE_STEP); Draw_view(S); break;
	case (SHIFT):      Tilt(S,(*psi)-many*ROTATE_STEP); Draw_view(S); break;
	case (ALT):        Roll(S,-ROTATE_STEP); Draw_view(S); break;
	case (ALT+SHIFT):  Roll(S,-many*ROTATE_STEP); Draw_view(S); break;
	/* case (CTRL):       m-=mstep; Read_fmt(0,m); Make_lists(S); Draw(S); break; */
	/* case (CTRL+SHIFT): if (m!=mmin) {m=mmin; Read_fmt(0,m);}; Make_lists(S); Draw(S); break; */
#if MARKER	
	case (CTRL+ALT): if (marker_size>0 && (*marker_z)<kmax) {(*marker_z)++; Draw(S);} else Tink(S); break;
	case (CTRL+ALT+SHIFT): if (marker_size>0 && (*marker_z)<=kmax-many) {(*marker_z)+=many; Draw(S);} else Tink(S); break;
#endif
	default:           Pop(S); break;
	}
	break;

      case XK_minus: 
	switch(state) {
	/* case (0): case (SHIFT): */
	/*   mstep*=-1; Draw_title(); break; */
	case (CTRL): case (CTRL+SHIFT):
	  (*distance)*=DISTANCE_FAC; setView(S, (*theta),(*phi),(*psi),(*distance)); Draw_view(S); break;
	default:           Pop(S); break;
	}
	break;

      case XK_equal:
	switch(state) {
	case (0): case (SHIFT):
	  Assign_par(S); break;
	case (CTRL): case (CTRL+SHIFT):
	  (*distance)/=DISTANCE_FAC; setView(S, (*theta),(*phi),(*psi),(*distance)); Draw_view(S); break;
	default:           Pop(S); break;
	}
	break;

      case XK_Shift_L: case XK_Shift_R: case XK_Control_L: case XK_Control_R: case XK_Caps_Lock: 
      case XK_Shift_Lock: case XK_Meta_L: case XK_Meta_R: case XK_Alt_L: case XK_Alt_R: case XK_Mode_switch:

	break;

      default: 
	message(4,"Unknown key %4lx pressed\n",theKeySym);
	Pop(S); break;
      }  /* switch(theKeySym) */
      break;

    case EnterNotify:
      /* This case is necessary for window managers which do not set keyboard
       * focus to theWindow when the pointer enters theWindow. */
      /* XSetInputFocus(theDisplay, theWindow, RevertToPointerRoot, CurrentTime);  */
      XSetInputFocus(theDisplay, theWindow, RevertToParent, CurrentTime); 
      break;
      
    case Expose:
      /* Window mapped for the first time and if you uncover some part of the
       * window. If you start paused and you see nothing in the window then
       * its possible that the problem is that the first Expose event is not
       * being caught for some reason. */
      Draw(S);
      break;

    case ConfigureNotify:
      /* Window size is changed by the user (or the window is initially
       * opened). Note that InitX contains code that constrains the window
       * to be square. */
      myReshape(S,theEvent.xconfigure.width, theEvent.xconfigure.height);
      break;

    case ButtonPress:
      if (theEvent.xbutton.button == 1) {
	Mouse_down(S,theEvent.xbutton.x, theEvent.xbutton.y);
      }
      break;

    case ButtonRelease:
      if (theEvent.xbutton.button == 1) {
	Mouse_up(S);
      }
      break;

    case MotionNotify:
      if (theEvent.xmotion.state & Button1Mask) {
	/* Remove all the MotionNotify events currently in the queue leaving
	 * the last one in theEvent. This ensures the image 'sticks' to the
	 * mouse. */
	char oldbusy=*(S->busychar);
	showbusy(S,BUSYTURNING);
	while (XCheckWindowEvent(theDisplay,theWindow,
				 PointerMotionMask, &theEvent));
	/* Rotate(S,theEvent.xmotion.x, theEvent.xmotion.y, 0., 0.); */
	Turn(S,theEvent.xmotion.x, theEvent.xmotion.y);
	Draw(S); 
	*(S->busychar)=oldbusy;
      }
      break;

    } /* switch (theEvent.type) */
  } /* while (XPending(theDisplay) || (ezmode != MODE_SIMULATING)) */

  if (*theta != old_theta || *phi != old_phi || *psi != old_psi || *distance != old_distance)
    MESSAGE("\x01/* ezview idev=%d camera change: theta="REALF" phi="REALF" psi="REALF" distance="REALF" */\n",
	    idev,*theta,*phi,*psi,*distance);

  if (stop_next_time) Pause(S);
 
  return(0);
} /* Event_check */
/* ========================================================================= */

int InitX (STR *S)
{
  /* Initializes X and opens a window. */
  DEVICE_VAR(Display *,theDisplay);
  DEVICE_VAR(int,theDWidth);
  DEVICE_VAR(int,theDHeight);
  DEVICE_VAR(Window,theWindow);
  DEVICE_VAR(INT,winx);
  DEVICE_VAR(INT,winy);
  DEVICE_VAR(INT,width);
  DEVICE_VAR(INT,height);
  DEVICE_VAR(int,horspace);
  DEVICE_VAR(int,vertspace);
  
  DEVICE_VAR(GLXContext,theGLXContext);
  DEVICE_VAR(XTextProperty,theWindowName);
  DEVICE_VAR(XTextProperty,theIconName);

  XVisualInfo           *theVisualInfo;
  Colormap              theColormap;
  int                   theScreen; 
  int                   theDepth;
  char                  *theDisplayName;
  XEvent                event;
  Atom                  del_atom;
  XSizeHints            theSizeHints;
  XSetWindowAttributes  theSWA;
  char                  *name = WINDOW_TITLE;
  int                   num1,num2;
  int list[] = {GLX_RGBA,
#if DOUBLEBUFFER
		       GLX_DOUBLEBUFFER,
#endif
		       GLX_RED_SIZE, 1,
		       GLX_GREEN_SIZE, 1,
		       GLX_BLUE_SIZE, 1,
		       GLX_ALPHA_SIZE, 1,
		       GLX_DEPTH_SIZE, 1,
		       None } ;
  /* DB: DEPTH and DOUBLEBUFFER are the important ones.  Perhaps need to add error
   * checking when DOUBLEBUFFER and DEPTH not available.  In the aux library
   * the first entry in list was GLX_LEVEL, 0, */

  /* Open the display */
  theDisplayName = XDisplayName(NULL);
  if NOT((*theDisplay) = XOpenDisplay(NULL))
    EXPECTED_ERROR("ERROR: Could not open a connection to X on display %s\n",theDisplayName);
  if NOT(glXQueryExtension((*theDisplay), &num1, &num2))
    EXPECTED_ERROR("ERROR: No glx extension on display %s\n",theDisplayName);

  theScreen     = DefaultScreen(*theDisplay);
  theDepth      = DefaultDepth ((*theDisplay), theScreen);
  *theDWidth     = DisplayWidth ((*theDisplay), theScreen);
  *theDHeight    = DisplayHeight((*theDisplay), theScreen);

  /* How much space for moving the window without goint beyond screen */
  (*horspace)=(*theDWidth)-(*width);
  (*vertspace)=(*theDHeight)-(*height);


  if NOT(theVisualInfo = glXChooseVisual((*theDisplay), theScreen, list))
	  EXPECTED_ERROR("ERROR: Couldn't find visual");
  if NOT((*theGLXContext) = glXCreateContext((*theDisplay),theVisualInfo,None,GL_TRUE))
  /* Last parameter indicates that, if possible, then render directly to
   * graphics hardware and bypass the X server. This should be faster. */
	  EXPECTED_ERROR("ERROR: Can not create a context");
  if NOT(theColormap = XCreateColormap(
				       (*theDisplay),
				       RootWindow((*theDisplay),theVisualInfo->screen),
				       theVisualInfo->visual, AllocNone))
	  /* AllocAll would generate a BadMatch.  */
	  EXPECTED_ERROR("Couldn't create Colormap");
  theSWA.colormap = theColormap;
  theSWA.border_pixel = 0;
  theSWA.event_mask = (EnterWindowMask | KeyPressMask | StructureNotifyMask | 
		       ButtonPressMask | ButtonReleaseMask | ExposureMask | 
		       PointerMotionMask);
  
  /* Move window to within the screen */
  if (*horspace) {
    while (*winx<0) (*winx)+=(*horspace);
    while (*winx>=(*horspace)) (*winx)-=(*horspace);
  }
  if (*vertspace) {
    while (*winy<0) (*winy)+=(*vertspace);	
    while (*winy>=(*vertspace)) (*winy)-=(*vertspace);
  }
  if NOT((*theWindow) = XCreateWindow((*theDisplay),
				      RootWindow((*theDisplay), theVisualInfo->screen),
				      (*winx), (*winy), (*width), (*height), 0,
				      theVisualInfo->depth, InputOutput,
				      theVisualInfo->visual,
				      CWBorderPixel|CWColormap|CWEventMask, &theSWA))
	  EXPECTED_ERROR("couldn't create X window");

  /* Set window properties theWindowName and theIconName 
   * to whatever variable name contains at the moment */

  XStringListToTextProperty(&name,1,theWindowName);
  XStringListToTextProperty(&name,1,theIconName);

  theSizeHints.base_width = (*width);
  theSizeHints.base_height = (*height);
  theSizeHints.min_aspect.x = (*width);   /* Maintain x:y ratio */
  theSizeHints.max_aspect.x = (*width);
  theSizeHints.min_aspect.y = (*height);
  theSizeHints.max_aspect.y = (*height);

  theSizeHints.flags = PSize|PAspect;

  if(!(WM_CTRLS_POS)) theSizeHints.flags |= USPosition;
  /* Setting USPosition here seems to make the WM honor the x and y
   * specified by XCreateWindow above.  Not setting this should give control
   * of position to the WM.  Note that the root window of an application is
   * special in that the WM has special privileges over its properties so
   * this may not work on all platforms.  */

  XSetWMProperties((*theDisplay), (*theWindow), theWindowName, theIconName,
		   NULL, 0, &theSizeHints, NULL, NULL);

  /* Express interest in WM killing this application  */
  if ((del_atom = XInternAtom((*theDisplay), "WM_DELETE_WINDOW", TRUE)) != None) {
    XSetWMProtocols((*theDisplay), (*theWindow), &del_atom, 1); 
  }

  XMapWindow((*theDisplay), (*theWindow));
  XIfEvent((*theDisplay), &event, WaitForNotify, (char *)(*theWindow));

  glXMakeCurrent((*theDisplay), (*theWindow), (*theGLXContext));
  XSetInputFocus((*theDisplay), (*theWindow), RevertToParent, CurrentTime); 

  /* Print useful information. I suggest printing this at least once */
  message(2,"\n/* %s version %d of the X Window System, X%d R%d\n"
	  "Color plane depth...........%d %s\n"
	  "Display Width...............%d \n"
	  "Display Height..............%d \n"
	  "The display: %s */\n",
	  ServerVendor(*theDisplay), VendorRelease(*theDisplay),ProtocolVersion(*theDisplay),ProtocolRevision(*theDisplay),
	  theDepth,(theDepth==1)?"(monochrome)":"",
	  (*theDWidth),
	  (*theDHeight),
	  theDisplayName
	  );
  return SUCCESS;
} /* InitX */
/* ========================================================================= */

static Bool WaitForNotify (Display *d, XEvent *e, char *arg) 
{
  /*  As seen in the Porting Guide. */
  return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
} /* WaitForNotify */
/* ========================================================================= */

void QuitX (STR *S)
{
  if (GRAPHICS) {
    XDestroyWindow(S->theDisplay,S->theWindow);
    XCloseDisplay(S->theDisplay);
  }
  return;
} /* QuitX */
/* ========================================================================= */

/* /\* Make an enumerated backup copy of the file if it already exists *\/ */
/* static void backup (STR *S, const char *fid) { */
/*   int n; */
/*   char bfid[BUFLEN]; */
/*   char cmd[BUFLEN]; */
/*   if (!fexist(fid)) return; */
/*   for (n=0;;n++) { */
/*     snprintf(bfid,BUFLEN,"%s.%d",fid,n); */
/*     if (!fexist(bfid)) break; */
/*   } */
/*   snprintf(cmd,BUFLEN,"mv %s %s",fid,bfid); */
/*   message(1,"%s\n",cmd); */
/*   system(cmd); */
/* } /\* backup *\/ */
/* ========================================================================= */

/* Save to a file a png image of what is currently on the screen using glReadPixels + netpbm converter. */
void Save_image (STR *S)
{
  DEVICE_CONST(Display*,theDisplay);
  DEVICE_CONST(Window,theWindow);
  DEVICE_CONST(int,width);
  DEVICE_CONST(int,height);
  DEVICE_ARRAY(char,images);
  DEVICE_ARRAY(char,imagescode);
  DEVICE_CONST(pp_fn,imagescompiled);

  PIPE *p;
  unsigned char *buf;
  char l[4*MAXPATH];

  int bufsize=width*height*3;
  /* #define DEBUG(...) {message(4,"Save image: "); message(4,__VA_ARGS__);} */

  if NOT(*images) {MESSAGE("Save_image not given images\n"); return;}
  k_on();
  snprintf(l,4*MAXPATH,images,*(real *)execute(imagescompiled));
  k_off();
  if NOT(p=pipeto(l)) {MESSAGE("Save_image could not pipe to '%s'\n",l); return;}
  fprintf(p->f,"P6\n%d %d\n%d\n",width,height,255);

  if NOT(theWindow) {MESSAGE("Save_image called while there is no visual\n"); return;}
  MALLOC(buf,bufsize);
  XRaiseWindow(theDisplay,theWindow);
  glFlush();
  glPixelStorei(GL_PACK_ALIGNMENT,1);
  glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,buf);

  if NOT(fwrite(buf,1,bufsize,p->f)) {MESSAGE("Save_image error while piping to '%s'\n",l); return;}

  if (0!=pipeclose(p)) {MESSAGE("Save_image could not close pipe to '%s'\n",l); return;}
  free(buf);
  
} /* Save_image */
/* ========================================================================= */

/* Save the filament buffer to actual filament file, if any */
void Save_filament (STR *S)
{
  DEVICE_CONST(FILE *,filament);
  DEVICE_ARRAY(char,filament_buffer);
  DEVICE_CONST(char *,filament_buffer_end);
  DEVICE_VAR(char *,filament_buffer_current);
  
  if (!filament) return;
  if (!filament_buffer) return;
  fputs(filament_buffer,filament);
  *filament_buffer_current=filament_buffer;
  fflush(filament);
  
} /* Save_filament */
/* ========================================================================= */

static void showbusy (STR *S, char c)
{
  if (!GRAPHICS) return;
  DEVICE_CONST(Display*,theDisplay);
  DEVICE_CONST(Window,theWindow);
  DEVICE_CONST(GLXContext,theGLXContext);  
  DEVICE_VAR(XTextProperty,theWindowName);
  DEVICE_ARRAY(char,longtitle);
  DEVICE_CONST(char *,busychar);
  if (busychar) *busychar=c;
  glXMakeCurrent(theDisplay, theWindow, theGLXContext);
  XStringListToTextProperty(&longtitle, 1, theWindowName);
  XSetWMName(theDisplay, theWindow, theWindowName);
  XFlush(theDisplay);
} /* showbusy */
/* ========================================================================= */

/* Calculate 3D window based on the current box shape */
static void Proportions (STR *S)
{
  DEVICE_CONST(int,imin);
  DEVICE_CONST(int,imax);
  DEVICE_CONST(int,jmin);
  DEVICE_CONST(int,jmax);
  DEVICE_CONST(int,kmin);
  DEVICE_CONST(int,kmax);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_ARRAY(real,plot_ini);
  real Lx=imax-imin+1;
  real Ly=jmax-jmin+1;
  real Lz=kmax-kmin+1;
  real Lmax=max(max(Lx,Ly),Lz);
  plot_length[0] = Lx/Lmax;
  plot_length[1] = Ly/Lmax;
  plot_length[2] = Lz/Lmax;
  plot_ini[0] = (imin-1)/Lmax;
  plot_ini[1] = (jmin-1)/Lmax;
  plot_ini[2] = (kmin-1)/Lmax;
} /* Proportions */
/* ========================================================================= */

/* Sound effects: relying on whatever XKeyboardControl can do */
void Pop (STR *S)
{
  XKeyboardControl kc;
  unsigned long vm;
  if (!GRAPHICS) return;
  vm = KBBellPitch;
  kc.bell_pitch = 220;
  XChangeKeyboardControl(S->theDisplay, vm, &kc);
  XBell(S->theDisplay,100);
} /* Pop */
/* ========================================================================= */

void Tink (STR *S)
{
  XKeyboardControl kc;
  unsigned long vm;
  if (!GRAPHICS) return;
  vm = KBBellPitch;
  kc.bell_pitch = 880;
  XChangeKeyboardControl(S->theDisplay, vm, &kc);
  XBell(S->theDisplay,100);
} /* Tink */
/* ========================================================================= */

#endif
