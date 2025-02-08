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

/**
 * This is obsoleted by ezpaint
 */

/* #define DBG(msg) printf("########### %s:%d %s\n",__FILE__,__LINE__,msg) */
#define DBG(msg) 

/* Paint on screen */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "system.h"
#include "beatbox.h"

#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "k_.h"
#include "pipe.h"

#if defined(NOX11) || defined (NOGL)
  #include "nograph.h"
  NOGRAPH_DUMMY(k_paintgl)
#else
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h> 
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <GL/gl.h>
#include <GL/glx.h>

#define TRUE             1   
#define MAXWINDOWTITLE 512
#define Xmin (-1.0)
#define Ymin (-1.0)
#define Xmax (1.0)
#define Ymax (1.0)


extern int Verbose;            /* defined in main */

typedef struct {
  #include "k_code.h"
  INT nabs;
  REAL absmin, absmax;
  INT nord;
  REAL ordmin, ordmax;
  INT appmin, appmax;
  REAL abs, ord, app, red, grn, blu;

  char filter[1024];		/* generated image will be piped to this unix command */
  char filtercode[1024];	/* calculate a number to make part of the filter */
  pp_fn filtercompiled;		/* k_code for the filter number calculation */
  
  char title[1024];		/* window title template */
  char titlecode[1024];         /* calculate a number to make part of the title */
  char fulltitle[MAXWINDOWTITLE]; /* window title filled */
  pp_fn titlecompiled;		/* k_code for the title number calculation */

  Display *theDisplay;
  Window theWindow;
  GLXContext theGLXContext;
  XTextProperty theWindowName, theIconName;
  int winx, winy, width, height;
  XEvent event;
} STR;

static Bool WaitForNotify(Display *d, XEvent *e, char *arg) 
{
  /*  DB: As seen in the Porting Guide. */
  return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}

RUN_HEAD(k_paintgl) {
  DEVICE_CONST(Display *,theDisplay);
  DEVICE_CONST(Window,theWindow);
  DEVICE_CONST(int,width);
  DEVICE_CONST(int,height);
  DEVICE_CONST(GLXContext,theGLXContext);
  DEVICE_CONST(XTextProperty,theWindowName);
  DEVICE_CONST(XTextProperty,theIconName);
  DEVICE_ARRAY(char, title);
  DEVICE_ARRAY(char, fulltitle);
  DEVICE_CONST(pp_fn, titlecompiled);
  #include "k_def.h"

  DEVICE_ARRAY(char, filter);
  DEVICE_ARRAY(char, filtercode);
  DEVICE_CONST(pp_fn, filtercompiled);

  DEVICE_CONST(REAL,absmin) DEVICE_CONST(REAL,absmax) DEVICE_VAR(double,abs) DEVICE_CONST(INT,nabs);
  DEVICE_CONST(REAL,ordmin) DEVICE_CONST(REAL,ordmax) DEVICE_VAR(double,ord) DEVICE_CONST(INT,nord); 
  DEVICE_CONST(INT,appmin) DEVICE_CONST(INT,appmax) DEVICE_VAR(double,app);
  DEVICE_VAR(double,red) DEVICE_VAR(double,grn) DEVICE_VAR(double,blu);

  int iapp, iord, iabs;
  GLfloat X1, X2, Y1, Y2;

  DBG("started");
  
  k_on();
  snprintf(fulltitle,MAXWINDOWTITLE,title,*(real *)execute(titlecompiled));
  k_off();
  char *titleref=&(S->fulltitle[0]);
  DBG("title made");
  XStringListToTextProperty(&titleref,1,&(S->theWindowName));
  XStringListToTextProperty(&titleref,1,&(S->theIconName));
  XSetWMName    (theDisplay, theWindow, &(S->theWindowName));
  XSetWMIconName(theDisplay, theWindow, &(S->theIconName));
  DBG("title displayed");

  glXMakeCurrent(theDisplay, theWindow, theGLXContext);
  DBG("window is current");
  k_on();

  for(iabs=0;iabs<nabs;iabs++) { 
    *abs=absmin+iabs*(absmax-absmin)/(nabs-1);
    X1=Xmin+((iabs+0)*(Xmax-Xmin))/nabs;
    X2=Xmin+((iabs+1)*(Xmax-Xmin))/nabs;
    for(iord=0;iord<nord;iord++) { 
      *ord=ordmin+iord*(ordmax-ordmin)/(nord-1);
      Y1=Ymin+((iord+0)*(Ymax-Ymin))/nord;
      Y2=Ymin+((iord+1)*(Ymax-Ymin))/nord;
      *red=*grn=*blu=0;
      for(*app=iapp=appmin;iapp<=appmax;*app=++iapp) {
        #include "k_exec.h"
      }
      glColor3f((GLfloat)(*red),(GLfloat)(*grn),(GLfloat)(*blu));
      glRectf(X1,Y1,X2,Y2);
    }
  }
  DBG("picture plotted");
  k_off();
  glXSwapBuffers(theDisplay,theWindow); 
  DBG("buffers swapped");
  glFlush();
  DBG("GL flushed");
  XFlush(theDisplay);
  DBG("X11 flushed");
  if (*filter) {
    PIPE *p;
    int bufsize=width*height*3;
    char l[4*MAXPATH];
    char *buf;
    MALLOC(buf,bufsize);
    DBG("buf allocated");

    snprintf(l,4*MAXPATH,filter,*(real *)execute(filtercompiled));
    DBG("filtercompiled executed");
    if NOT(p=pipeto(l)) return 1;
    DBG("pipe open");
    fprintf(p->f,"P6\n%d %d\n%d\n",width,height,255);
    DBG("header written");

    if (!theWindow) ABORT("no window???");
    XRaiseWindow(theDisplay,theWindow);
    DBG("window raised");
    glFinish();
    DBG("GL finished");
    glPixelStorei(GL_PACK_ALIGNMENT,1);
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,buf);
    DBG("pixels read");
    fwrite(buf,1,bufsize,p->f);
    DBG("buffer written");
    
    if (0!=pipeclose(p)) ABORT("Save_image could not close pipe");
    DBG("pipe closed");
    free(buf);
    DBG("buf freed");
  }
  DBG("window filtered");
} RUN_TAIL(k_paintgl)

DESTROY_HEAD(k_paintgl)
  #include "k_free.h"
DESTROY_TAIL(k_paintgl)

CREATE_HEAD(k_paintgl)
{
  static char			*theDisplayName = NULL;
  static XVisualInfo		*theVisualInfo;
  static Colormap		theColormap;
  static int			theScreen; 
  static int			theDepth;
  static int			theDWidth;
  static int			theDHeight;
  static Atom			del_atom;
  static XSizeHints		theSizeHints;
  static XSetWindowAttributes	theSWA;
  /* static XTextProperty		theWindowName, theIconName; */
  static int			num1,num2;
  static int list[] = {GLX_RGBA,
		       GLX_DOUBLEBUFFER,
		       GLX_RED_SIZE, 1,
		       GLX_GREEN_SIZE, 1,
		       GLX_BLUE_SIZE, 1,
		       GLX_DEPTH_SIZE, 1,
		       None } ;
  static int horspace, vertspace;

  DEVICE_IS_SPACELESS;
  
  k_on();									CHK(NULL);

  memcpy(loctb,deftb,sizeof(*deftb));
  tb_insert_real(loctb,"abs",&(S->abs));    CHK("abs");
  tb_insert_real(loctb,"ord",&(S->ord));    CHK("ord");
  tb_insert_real(loctb,"app",&(S->app));    CHK("app");
  tb_insert_real(loctb,"red",&(S->red));    CHK("red");
  tb_insert_real(loctb,"grn",&(S->grn));    CHK("grn");
  tb_insert_real(loctb,"blu",&(S->blu));    CHK("blu");
  #include "k_comp.h"
  if(!used(S->data,S->ncode,&(S->red)))EXPECTED_ERROR("/*WARNING: variable 'red' never assigned!!*/");
  if(!used(S->data,S->ncode,&(S->grn)))EXPECTED_ERROR("/*WARNING: variable 'grn' never assigned!!*/");
  if(!used(S->data,S->ncode,&(S->blu)))EXPECTED_ERROR("/*WARNING: variable 'blu' never assigned!!*/");

  k_off();

  ACCEPTL(nabs,INONE,2,INONE);
  ACCEPTR(absmin,0,INONE,INONE);
  ACCEPTR(absmax,S->nabs-1,INONE,INONE);
  ACCEPTL(nord,INONE,2,INONE);
  ACCEPTR(ordmin,0,INONE,INONE);
  ACCEPTR(ordmax,S->nord-1,INONE,INONE);

  ACCEPTL(appmin,0,INONE,INONE);
  ACCEPTL(appmax,S->appmin,S->appmin,INONE);

  ACCEPTS(filter,"");
  if (*(S->filter)) {
    ACCEPTS(filtercode,"t");
    k_on();
    S->filtercompiled=compile(S->filtercode,deftb,t_real); CHK(S->filtercode);
    k_off();
  }
  
  ACCEPTI(winx,-1,INONE,INONE); /* default: next to right edge of screen */
  ACCEPTI(winy,+1,INONE,INONE); /* default: next to top edge of screen */
  ACCEPTI(width,512,1,INONE);
  ACCEPTI(height,512,1,INONE);

  ACCEPTS(title,"k_paintgl t=%.0f");
  ACCEPTS(titlecode,"t");
  S->titlecompiled=compile(S->titlecode,deftb,t_real); CHK(S->titlecode);

  /* Open the display */
  if NOT(S->theDisplay = XOpenDisplay(NULL)) 
	  EXPECTED_ERROR("Could not open a connection to X on display %s\n",
			 XDisplayName(theDisplayName));
  if NOT(glXQueryExtension(S->theDisplay, &num1, &num2)) 
	  EXPECTED_ERROR("No glx extension on display %s\n",
			 XDisplayName(theDisplayName));
  theScreen     = DefaultScreen(S->theDisplay);
  theDepth      = DefaultDepth (S->theDisplay, theScreen);
  theDWidth     = DisplayWidth (S->theDisplay, theScreen);
  theDHeight    = DisplayHeight(S->theDisplay, theScreen);

  if NOT(theVisualInfo = glXChooseVisual(S->theDisplay, theScreen, list)) 
	  EXPECTED_ERROR("ERROR: Couldn't find visual");
  if NOT(S->theGLXContext = glXCreateContext(S->theDisplay, theVisualInfo,None,GL_TRUE)) 
	  EXPECTED_ERROR("Can not create a context");
  if NOT(theColormap = XCreateColormap(
				       S->theDisplay,
				       RootWindow(S->theDisplay, theVisualInfo->screen),
				       theVisualInfo->visual, AllocNone))
	  EXPECTED_ERROR("Couldn't create Colormap");
  theSWA.colormap = theColormap;
  theSWA.border_pixel = 0;
  theSWA.event_mask = (EnterWindowMask | KeyPressMask | StructureNotifyMask |
  		       ButtonPressMask | ButtonReleaseMask | ExposureMask |
  		       PointerMotionMask);

  /* Move window to within the screen */
  horspace=theDWidth-width;
  if (horspace>0) {
    while (winx<0) winx+=horspace;
    while (winx>=horspace) winx-=horspace;
  }
  vertspace=theDHeight-height;
  if (vertspace>0) {
    while (winy<0) winy+=vertspace;	
    while (winy>=vertspace) winy-=vertspace;
  }
  if NOT(S->theWindow=XCreateWindow(S->theDisplay,
				 RootWindow(S->theDisplay, theVisualInfo->screen),
				 winx, winy, width, height, 0,
				 theVisualInfo->depth, InputOutput,
				 theVisualInfo->visual,
				 CWBorderPixel|CWColormap|CWEventMask, &theSWA))
	  EXPECTED_ERROR("couldn't create X window");

  snprintf(S->fulltitle,MAXWINDOWTITLE,S->title,*(real *)execute(S->titlecompiled));
  char *titleref=&(S->fulltitle[0]);
  XStringListToTextProperty(&titleref,1,&(S->theWindowName));
  XStringListToTextProperty(&titleref,1,&(S->theIconName));
  XSetWMName    (S->theDisplay,S->theWindow,&(S->theWindowName));
  XSetWMIconName(S->theDisplay,S->theWindow,&(S->theIconName));

  theSizeHints.base_width = width;
  theSizeHints.base_height = height;
  theSizeHints.min_aspect.x = width;   /* Maintain x:y ratio */
  theSizeHints.max_aspect.x = width;
  theSizeHints.min_aspect.y = height;
  theSizeHints.max_aspect.y = height;

  theSizeHints.flags = PSize|PAspect;
  theSizeHints.flags |= USPosition;
  XSetWMProperties(S->theDisplay, S->theWindow, &(S->theWindowName), &(S->theIconName),
		   NULL, 0, &theSizeHints, NULL, NULL);
  if ((del_atom = XInternAtom(S->theDisplay, "WM_DELETE_WINDOW", TRUE)) != None) {
    XSetWMProtocols(S->theDisplay, S->theWindow, &del_atom, 1); 
  }
  XMapWindow(S->theDisplay, S->theWindow);
  XIfEvent(S->theDisplay, &S->event, WaitForNotify, (char *)S->theWindow);

  glXMakeCurrent(S->theDisplay, S->theWindow, S->theGLXContext);
  MESSAGE("/*\n");
  MESSAGE("%s version %d of the X Window System, X%d R%d\n",
	  ServerVendor    (S->theDisplay),
	  VendorRelease   (S->theDisplay),
	  ProtocolVersion (S->theDisplay),
	  ProtocolRevision(S->theDisplay));

  if(theDepth==1) {
    MESSAGE("Color plane depth...........%d (monochrome)\n", theDepth);
  } else {
    printf("Color plane depth...........%d \n",             theDepth);
  }

  MESSAGE("Display Width...............%d \n", theDWidth);
  MESSAGE("Display Height..............%d \n", theDHeight);
  MESSAGE("The display %s\n", XDisplayName(theDisplayName));
  MESSAGE("*/\n");

  /* S->Xmin=1; S->Xmax=S->width; */
  /* S->Ymin=1; S->Ymax=S->height; */
 
}
CREATE_TAIL(k_paintgl,0)

#endif
