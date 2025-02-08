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
 * This device emulates D. Barkley's EZSPIRAL (version 3.1) 2D graphical output, with minor extensions. 
 * Does plot the tip trajectory if produced by singz (or any other device), 
 * but does not engage in dialogue. 
 * Options: paint k-functions for flexibility, or grid values for speed. 
 * Extra facility: put one or two markers on the screen (as in ezride).
 */

#include <assert.h>
#include <regex.h>
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
  NOGRAPH_DUMMY(ezpaint)
#else
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h> 
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <GL/gl.h>
#include <GL/glx.h>

#define TRUE 1   
#define MAXWINDOWTITLE 512
#define Xmin (-1.0)
#define Ymin (-1.0)
#define Xmax (1.0)
#define Ymax (1.0)
#define PX(x) (Xmin+((x)*(Xmax-Xmin))/nabs)
#define PY(y) (Ymin+((y)*(Ymax-Ymin))/nord)

extern int Verbose;            /* defined in main */

typedef enum {xy, yx, xz, zx, yz, zy} aspect_t;

typedef struct {
  /* The raster sizes */
  int nabs;			/* horizontal */
  int nord;			/* vertical */
  int napp;			/* depth  */

  /* k-function painting option */
  #include "k_code.h"		/* k-program calculating the colours */
  real absmin, absmax;  	/* real ..  */
  real ordmin, ordmax;		/* .. raster .. */
  real appmin, appmax;		/* .. limits */
  REAL abs, ord, app, red, grn, blu; /* k-variables involved */

  /* grid values painting option */
  char aspect[2];		/* orientation of the 2D projection of the 3D box */
  aspect_t aspectchoice;	/* numerical code of the same */
  int abs0, abs1;		/* integer (grid) .. */
  int ord0, ord1;		/* .. raster .. */
  int app0, app1;		/* .. limits */
  int rlayer, glayer, blayer;	/* grid layers involved */
  real rmin, gmin, bmin;	/* grid values scaling minima */
  real rmax, gmax, bmax;	/* grid values scaling maxima */
  real bg_r, bg_g, bg_b;	/* background colour to represent voids */

  /* plotting "tip path" */
  char xtip[1024];	        /* k-code how to calculate the tip position */
  char ytip[1024];	        /* k-code how to calculate the tip position */
  pp_fn xtipcompiled;           /* compiled k-code of the same */
  pp_fn ytipcompiled;           /* compiled k-code of the same */
  int ntip;                     /* number of entries in the lists */
  double *tipx;                 /* list ot x positions */
  double *tipy;                 /* list of y positions */
  int ntipmax;                  /* maximal .. */

  /* plotting markers */
  #define _(Type,Name,Init,miN,maX) PAR(Type,Name)
  #include "ezpaint.h" /* dynamically linked parameters */
  #undef _
  
  /* saving copy of the screen to files */
  char filter[1024];		/* generated image will be piped to this unix command */
  char filtercode[1024];	/* calculate a number to make part of the filter */
  pp_fn filtercompiled;		/* k_code for the filter number calculation */
  
  /* window title */
  char title[1024];		/* window title template */
  char titlecode[1024];         /* calculate a number to make part of the title */
  char fulltitle[MAXWINDOWTITLE]; /* window title filled */
  pp_fn titlecompiled;		/* k_code for the title number calculation */

  /* X and GL stuff */
  int winx, winy, width, height; /* window coordinates and sizes */
  int doublebuffer;		/* flag of whether to use GL double-buffering */
  Display *theDisplay;		
  Window theWindow;
  GLXContext theGLXContext;
  XTextProperty theWindowName, theIconName;
  XEvent event;
} STR;

static real crop(real a, real b, real c)
{
  return (a<=b)?b:(a>=c)?c:a;
}

static void  Draw_markers (STR *S);

/*  DB: As seen in the Porting Guide. */
static Bool WaitForNotify(Display *d, XEvent *e, char *arg) 
{
  return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}


/* Check GL errors */
static GLenum err;
#include "glerr.h"
#define CGL(msg) 							\
  while ( (err = glGetError()) != GL_NO_ERROR ) {			\
    MESSAGE("%s:%d %s GL error %d '%s'\n",__FILE__,__LINE__,msg,err,glerr(err)); \
  }


RUN_HEAD(ezpaint)
{
  DEVICE_CONST(Display *, theDisplay);
  DEVICE_CONST(Window, theWindow);
  DEVICE_CONST(int, width);
  DEVICE_CONST(int, height);
  DEVICE_CONST(GLXContext, theGLXContext);
  DEVICE_VAR(XTextProperty, theWindowName);
  DEVICE_VAR(XTextProperty, theIconName);
  DEVICE_CONST(int, doublebuffer);
  DEVICE_ARRAY(char, title);
  DEVICE_ARRAY(char, fulltitle);
  DEVICE_CONST(pp_fn, titlecompiled);
  #include "k_def.h"

  DEVICE_ARRAY(char, filter);
  DEVICE_ARRAY(char, filtercode);
  DEVICE_CONST(pp_fn, filtercompiled);

  DEVICE_CONST(int, nabs);
  DEVICE_CONST(int, nord); 
  DEVICE_CONST(int, napp); 
  DEVICE_CONST(aspect_t,aspectchoice);

  DEVICE_CONST(pp_fn, xtipcompiled);
  DEVICE_CONST(pp_fn, ytipcompiled);
  DEVICE_VAR(int, ntip);
  DEVICE_ARRAY(double, tipx);
  DEVICE_ARRAY(double, tipy);
  DEVICE_CONST(int, ntipmax);

  GLfloat X1, X2, Y1, Y2;
  int itip;
  REAL result;

  REAL the_xtip, the_ytip;
  
  /* Update the window and icon title */
  k_on();
  snprintf(fulltitle,MAXWINDOWTITLE,title,*(real *)execute(titlecompiled));
  k_off();
  glXMakeCurrent(theDisplay, theWindow, theGLXContext); 	CGL("glXMakeCurrent");
  XStringListToTextProperty(&fulltitle, 1, theWindowName);	CGL("XStringListToTextProperty");
  XStringListToTextProperty(&fulltitle, 1, theIconName);	CGL("XStringListToTextProperty");
  XSetWMName    (theDisplay, theWindow, theWindowName);		CGL("XSetWMName");
  XSetWMIconName(theDisplay, theWindow, theIconName);		CGL("XSetWMIconName");
  XFlush(theDisplay);						CGL("XFlush");

  /* Paint the raster */
  if (ncode) {
    /* Arbitrary raster and colours defined by k-functions */
    DEVICE_CONST(real,absmin) DEVICE_CONST(real,absmax) DEVICE_VAR(REAL,abs);
    DEVICE_CONST(real,ordmin) DEVICE_CONST(real,ordmax) DEVICE_VAR(REAL,ord);
    DEVICE_CONST(real,appmin) DEVICE_CONST(real,appmax) DEVICE_VAR(REAL,app);
    DEVICE_VAR(REAL,red) DEVICE_VAR(REAL,grn) DEVICE_VAR(REAL,blu);
    int iapp, iord, iabs;

    k_on();
    for(iabs=0;iabs<nabs;iabs++) { 
      *abs=(REAL)(absmin+iabs*(absmax-absmin)/(nabs-1));
      X1=PX(iabs+0);
      X2=PX(iabs+1);
      for(iord=0;iord<nord;iord++) { 
	*ord=(REAL)(ordmin+iord*(ordmax-ordmin)/(nord-1));
	Y1=PY(iord+0);
	Y2=PY(iord+1);
	*red=*grn=*blu=0;
	for(iapp=0;iapp<napp;iapp++) { 
	  *app=(REAL)((napp>1)?(appmin+iapp*(appmax-appmin)/(napp-1)):appmin);
          #include "k_exec.h"
	} /* for iapp */
	glColor3f((GLfloat)(*red),(GLfloat)(*grn),(GLfloat)(*blu));
	glRectf(X1,Y1,X2,Y2);
      } /* for iord */
    } /* for iabs */
    k_off();
    
  } else { /* if ncode */
    /* Raster as subset of grid and colours linearly depending on grid values */
    DEVICE_CONST(int,abs0) DEVICE_CONST(int,ord0) DEVICE_CONST(int,app0);
    DEVICE_CONST(int,abs1) DEVICE_CONST(int,ord1) DEVICE_CONST(int,app1);
    DEVICE_CONST(int,rlayer) DEVICE_CONST(real,rmin) DEVICE_CONST(real,rmax);
    DEVICE_CONST(int,glayer) DEVICE_CONST(real,gmin) DEVICE_CONST(real,gmax);
    DEVICE_CONST(int,blayer) DEVICE_CONST(real,bmin) DEVICE_CONST(real,bmax);
    DEVICE_CONST(real,bg_r)  DEVICE_CONST(real,bg_g) DEVICE_CONST(real,bg_b);

    /* Todo: should this stuff be precomputed? */
    real rbase=(rlayer>=0)?(-rmin/(rmax-rmin)):0;
    real gbase=(glayer>=0)?(-gmin/(gmax-gmin)):0;
    real bbase=(blayer>=0)?(-bmin/(bmax-bmin)):0;
    real rcoef=(rlayer>=0)?(1.0/(rmax-rmin)):0;
    real gcoef=(glayer>=0)?(1.0/(gmax-gmin)):0;
    real bcoef=(blayer>=0)?(1.0/(bmax-bmin)):0;
    int r=(rlayer>=0)?rlayer:0;
    int g=(glayer>=0)?glayer:0;
    int b=(blayer>=0)?blayer:0;
    real red, grn, blu;
    real rsum, gsum, bsum;
    int iapp, iord, iabs;
    int x, y, z;
    int nval;

    for (iabs=abs0; iabs<=abs1; iabs++) {
      X1=PX(iabs-abs0+0);
      X2=PX(iabs-abs0+1);
      for (iord=ord0; iord<=ord1; iord++) {
	Y1=PY(iord-ord0+0);
	Y2=PY(iord-ord0+1);
	rsum=gsum=bsum=0.0;
	nval=0;
	for (iapp=app0; iapp<=app1; iapp++) {
	  switch (aspectchoice) {
	  case xy: x=iabs; y=iord; z=iapp; break;
	  case yx: y=iabs; x=iord; z=iapp; break;
	  case xz: x=iabs; z=iord; y=iapp; break;
	  case zx: z=iabs; x=iord; y=iapp; break;
	  case yz: y=iabs; z=iord; x=iapp; break;
	  case zy: z=iabs; y=iord; x=iapp; break;
	  default: ABORT("illegal aspectchoice=%d\n",(int)aspectchoice);
	  }
	  if (isTissue(x,y,z)) {
	    rsum+=New[ind(x,y,z,r)];
	    gsum+=New[ind(x,y,z,g)];
	    bsum+=New[ind(x,y,z,b)];
	    nval++;
	  } /* isTissue */
	} /* for iapp */
	if (nval) {
	  red=(rlayer>=0)?crop(rbase+rsum/nval*rcoef,0.0,1.0):0;
	  grn=(glayer>=0)?crop(gbase+gsum/nval*gcoef,0.0,1.0):0;
	  blu=(blayer>=0)?crop(bbase+bsum/nval*bcoef,0.0,1.0):0;
	} else { /*  Void. Use background colour. */
	  red=bg_r;
	  grn=bg_g;
	  blu=bg_b;
	} /* if nsums else */
	glColor3f((GLfloat)(red),(GLfloat)(grn),(GLfloat)(blu));
	glRectf(X1,Y1,X2,Y2);
      } /* for iord */
    } /* for iabs */
  } /* if ncode else */
  
  /* Draw the tip trace if needed */
  if ( (xtipcompiled!=NULL) && (ytipcompiled!=NULL) ) {
    /* if stack is full, discard the earliest record */
    if (*ntip==ntipmax) {
      for (itip=0;itip<(*ntip)-1;itip++) {
	tipx[itip]=tipx[itip+1];
	tipy[itip]=tipy[itip+1];
      }
      (*ntip)--;
    }
    
    /* x and y codes to be in the same units (grid coords) as iabs and iord */
    k_on();
    if (xtipcompiled!=NULL) 
      the_xtip=*(REAL *)execute(xtipcompiled);
    if (ytipcompiled!=NULL) 
      the_ytip=*(REAL *)execute(ytipcompiled);
    if (xtipcompiled!=NULL && ytipcompiled!=NULL && the_xtip!=real_inf && the_ytip !=real_inf) {
      tipx[*ntip]=the_xtip;
      tipy[*ntip]=the_ytip;
      (*ntip)++;
    }
    k_off();
    
    /* Now plot the tip trace */
    /* To do: make these device parameters one day? */
#define TIP_PLOT_TYPE GL_LINE_STRIP
#define TIP_WT  1.0
#define TIP_R   1.0
#define TIP_G   1.0
#define TIP_B   1.0
    glLineWidth (TIP_WT);
    glBegin (TIP_PLOT_TYPE);
    glColor3f (TIP_R, TIP_G, TIP_B);
    for (itip=0; itip<(*ntip); itip++) {
      GLfloat abs, ord, X, Y;
      abs=tipx[itip];
      ord=tipy[itip];
      X=PX(abs+0.5);
      Y=PY(ord+0.5);
      glVertex2f (X, Y);
    }
    glEnd ();
  }

  Draw_markers(S);
  
  /* Important: read pixels before swapping buffers */
  if (*filter) {
    PIPE *p;
    int bufsize=width*height*3;
    char l[4*MAXPATH];
    char *buf;
    MALLOC(buf,bufsize);
    k_on();
    snprintf(l,4*MAXPATH,filter,*(real *)execute(filtercompiled));
    k_off();
    if NOT(p=pipeto(l)) ABORT("could not open pipe to %s\n",l);
    fprintf(p->f,"P6\n%d %d\n%d\n",width,height,255);
    if (!theWindow) ABORT("no window???");
    XRaiseWindow(theDisplay,theWindow);				CGL("XRaiseWindow");
    glFinish();							CGL("glFinish");
    glPixelStorei(GL_PACK_ALIGNMENT,1); 			CGL("glPixelStorei");
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,buf); CGL("glReadPixels");
    if (bufsize!=fwrite(buf,1,bufsize,p->f)) ABORT("could not write to pipe\n");
    if (0!=pipeclose(p)) ABORT("could not close pipe\n");
    free(buf);
  }

  if (doublebuffer) {
    glXSwapBuffers(theDisplay,theWindow);    			CGL("glXSwapBuffers");
  } else {
    glFlush();							CGL("glFlush");
  }
}
RUN_TAIL(ezpaint)

DESTROY_HEAD(ezpaint)
  #include "k_free.h"
DESTROY_TAIL(ezpaint)

CREATE_HEAD(ezpaint)
{
  char			*theDisplayName = NULL;
  XVisualInfo		*theVisualInfo;
  Colormap		theColormap;
  int			theScreen; 
  int			theDepth;
  int			theDWidth;
  int			theDHeight;
  Atom			del_atom;
  XSizeHints		theSizeHints;
  XSetWindowAttributes	theSWA;
  int			num1,num2;
  int list[] = {GLX_RGBA,
		GLX_RED_SIZE, 1,
		GLX_GREEN_SIZE, 1,
		GLX_BLUE_SIZE, 1,
		GLX_DEPTH_SIZE, 1,
		None, /* placeholder for possible GLX_DOUBLEBUFFER */
		None} ;
  int doublebuffer_pos=9; /* position of the placeholder in the list[] above */
  int horspace, vertspace;
  Space s=dev->s;
  char *p;
  regex_t re_compiled;

  DEVICE_IS_RECTANGULAR;
  
  /* The regex check is to make sure that an empty block does not count for making this decition */
  if (regcomp(&re_compiled, "^[[:space:]]*\\{[[:space:]]*\\}", REG_EXTENDED|REG_NOSUB) != 0)
    EXPECTED_ERROR("could not compile the regex for an empty block\n");
  if (NULL!=(p=find_key("pgm=",w)) && 0!=regexec(&re_compiled, p,0,NULL,0)) {
    /* Now we know the pgm= block is given and it is not empty */
    MESSAGE("/* k-colours version understood */");
    k_on();									CHK(NULL);
    
    memcpy(loctb,deftb,sizeof(*deftb));
    tb_insert_real(loctb,"abs",&(S->abs));    CHK("abs");
    tb_insert_real(loctb,"ord",&(S->ord));    CHK("ord");
    tb_insert_real(loctb,"app",&(S->app));    CHK("app");
    tb_insert_real(loctb,"red",&(S->red));    CHK("red");
    tb_insert_real(loctb,"grn",&(S->grn));    CHK("grn");
    tb_insert_real(loctb,"blu",&(S->blu));    CHK("blu");
#include "k_comp.h"
    if(!used(S->data,S->ncode,&(S->red)))EXPECTED_ERROR("/* WARNING: variable 'red' never assigned! */");
    if(!used(S->data,S->ncode,&(S->grn)))EXPECTED_ERROR("/* WARNING: variable 'grn' never assigned! */");
    if(!used(S->data,S->ncode,&(S->blu)))EXPECTED_ERROR("/* WARNING: variable 'blu' never assigned! */");
    
    k_off();
    
    ACCEPTI(nabs,INONE,2,INONE);
    ACCEPTR(absmin,0,RNONE,RNONE);
    ACCEPTR(absmax,S->nabs-1,RNONE,RNONE);
    ACCEPTI(nord,INONE,2,INONE);
    ACCEPTR(ordmin,0,RNONE,RNONE);
    ACCEPTR(ordmax,S->nord-1,RNONE,RNONE);
    ACCEPTI(napp,1,1,INONE);
    ACCEPTR(appmin,0,RNONE,RNONE);
    ACCEPTR(appmax,0,S->appmin,RNONE);
    if (napp==1) ASSERT(appmin == appmax);

    #define IGNORE(x) if (find_key(#x "=",w)) MESSAGE("/* pgm={} is present, so the given " #x " will be ignored /*\n")
    IGNORE(rlayer); IGNORE(rmin); IGNORE(rmax);
    IGNORE(glayer); IGNORE(gmin); IGNORE(gmax);
    IGNORE(blayer); IGNORE(bmin); IGNORE(bmax);
    IGNORE(bg_r); IGNORE(bg_g); IGNORE(bg_b);
    #undef IGNORE
    S->rlayer=S->glayer=S->glayer=-1;
    S->rmin=S->gmin=S->bmin=0;
    S->rmax=S->gmax=S->bmax=1;
  } else { /* if findkey(pgm=) */
    /* Now ... the pgm= block is not given, or is given but is empty */
    MESSAGE("/* grid colours version understood */");
    S->ncode=0;
    #define IGNORE(x) if (find_key(#x "=",w)) MESSAGE("/* no pgm={} is present, so the given " #x " will be ignored /*\n");
    IGNORE(nabs); IGNORE(absmin); IGNORE(absmin);
    IGNORE(nord); IGNORE(ordmin); IGNORE(ordmin);
    IGNORE(napp); IGNORE(appmin); IGNORE(appmin);
    #undef IGNORE
    
    ACCEPTSN(aspect,2,"xy");
    if (0) {}
    #define CASE(t,a,o,p) else if (0==strcmp(aspect,#t)) { S->aspectchoice=t;	\
	  S->abs0=s.a##0; S->ord0=s.o##0; S->app0=s.p##0;			\
	  S->abs1=s.a##1; S->ord1=s.o##1; S->app1=s.p##1; }
    CASE(xy,x,y,z)
    CASE(yx,y,x,z)
    CASE(xz,x,z,y)
    CASE(zx,z,x,y)
    CASE(yz,y,z,x)
    CASE(zy,z,y,x)
    #undef CASE
    else EXPECTED_ERROR("invalid aspect '%s'\n",aspect);
    
    S->nabs=S->abs1-S->abs0+1;
    S->nord=S->ord1-S->ord0+1;
    S->napp=S->app1-S->app0+1;
    MESSAGE("/* abs=%d:%d ord=%d:%d app=%d:%d */\n", S->abs0, S->abs1, S->ord0, S->ord1, S->app0, S->app1);
    S->absmin=S->absmax=0;
    S->ordmin=S->ordmax=0;
    S->appmin=S->appmax=0;
    ACCEPTI(rlayer,-1,-1,vmax-1);
    ACCEPTR(rmin,0,RNONE,RNONE);
    ACCEPTR(rmax,1,RNONE,RNONE);
    ACCEPTI(glayer,-1,-1,vmax-1);
    ACCEPTR(gmin,0,RNONE,RNONE);
    ACCEPTR(gmax,1,RNONE,RNONE);
    ACCEPTI(blayer,-1,-1,vmax-1);
    ACCEPTR(bmin,0,RNONE,RNONE);
    ACCEPTR(bmax,1,RNONE,RNONE);
    ASSERT(rmin!=rmax);
    ASSERT(gmin!=gmax);
    ASSERT(bmin!=bmax);
    ACCEPTR(bg_r,0,0,1);
    ACCEPTR(bg_g,0,0,1);
    ACCEPTR(bg_b,0,0,1);
  } /* if findkey(pgm=) else */
    
  ACCEPTS(xtip,"");
  ACCEPTS(ytip,"");
  if ( ('\0'!=*(S->xtip)) && ('\0'!=*(S->ytip)) ) {
    S->xtipcompiled=compile(S->xtip,deftb,t_real); CHK(S->xtip);
    S->ytipcompiled=compile(S->ytip,deftb,t_real); CHK(S->ytip);
    ACCEPTI(ntipmax,100,1,INONE);
    CALLOC(S->tipx,S->ntipmax,sizeof(double));
    CALLOC(S->tipy,S->ntipmax,sizeof(double));
    S->ntip=0;
  } else { /* if xtip && ytip */
    S->ntipmax=0;
    S->xtipcompiled=NULL;
    S->ytipcompiled=NULL;
    S->tipx=NULL;
    S->tipy=NULL;
  } /* if xtip && ytip else */

  /* #define acceptREAL(b,c,d,e) if (!acceptrk(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); REAL b=S->b */
  #define acceptREAL(b,c,d,e) ACCEPTRK(b,c,d,e) REAL b=S->b
  /* #define acceptINT(b,c,d,e) if (!acceptik(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); INT b=S->b */
  #define acceptINT(b,c,d,e) ACCEPTIK(b,c,d,e) INT b=S->b
  #define _(Type,Name,Dflt,miN,maX) accept##Type(Name,Dflt,miN,maX);
  #include "ezpaint.h"
  #undef _
  
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
  
  ACCEPTI(doublebuffer,0,0,1);
  if (doublebuffer) list[doublebuffer_pos]=GLX_DOUBLEBUFFER;

  ACCEPTS(title,"ezpaint t=%.0f");
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

  k_on();
  snprintf(S->fulltitle,MAXWINDOWTITLE,S->title,*(real *)execute(S->titlecompiled));
  k_off();
  char *titleref=&(S->fulltitle[0]);
  XStringListToTextProperty(&titleref,1,&(S->theWindowName));	CGL("XStringListToTextProperty");
  XStringListToTextProperty(&titleref,1,&(S->theIconName));	CGL("XStringListToTextProperty");
  XSetWMName(S->theDisplay,S->theWindow,&(S->theWindowName));	CGL("XSetWMName");
  XSetWMIconName(S->theDisplay,S->theWindow,&(S->theIconName));	CGL("XSetWMIconName");

  theSizeHints.base_width = width;
  theSizeHints.base_height = height;
  theSizeHints.min_aspect.x = width;   /* Maintain x:y ratio */
  theSizeHints.max_aspect.x = width;
  theSizeHints.min_aspect.y = height;
  theSizeHints.max_aspect.y = height;

  theSizeHints.flags = PSize|PAspect;
  theSizeHints.flags |= USPosition;
  XSetWMProperties(S->theDisplay, S->theWindow, &(S->theWindowName), &(S->theIconName),	
		   NULL, 0, &theSizeHints, NULL, NULL);		CGL("XSetWMProperties");
  if ((del_atom = XInternAtom(S->theDisplay, "WM_DELETE_WINDOW", TRUE)) != None) {
    XSetWMProtocols(S->theDisplay, S->theWindow, &del_atom, 1);	CGL("XSetWMProtocols");
  }
  XMapWindow(S->theDisplay, S->theWindow);			CGL("XMapWindow");
  XIfEvent(S->theDisplay, &S->event, WaitForNotify, (char *)S->theWindow);	CGL("XIfEvent");

  glXMakeCurrent(S->theDisplay, S->theWindow, S->theGLXContext);		CGL("glXMakeCurrent");
  if (Verbose) {
    MESSAGE("/*\n");
    MESSAGE("%s version %d of the X Window System, X%d R%d\n",
	    ServerVendor    (S->theDisplay),
	    VendorRelease   (S->theDisplay),
	    ProtocolVersion (S->theDisplay),
	    ProtocolRevision(S->theDisplay));
    
    if(theDepth==1) {
      MESSAGE("Color plane depth...........%d (monochrome)\n", theDepth);
    } else {
      MESSAGE("Color plane depth...........%d \n",             theDepth);
    }
    
    MESSAGE("Display Width...............%d \n", theDWidth);
    MESSAGE("Display Height..............%d \n", theDHeight);
    MESSAGE("The display %s\n", XDisplayName(theDisplayName));
    MESSAGE("*/\n");
  } /* if Verbose */
}
CREATE_TAIL(ezpaint,0)
/* ========================================================================= */

#define GLRECT(x1,y1,x2,y2)  glRectf((x1),(y1),(x2),(y2))
#define GLCOLOR3(r,g,b)      glColor3f((r),(g),(b))
#define GLVERTEX2(x,y)       glVertex2f((x),(y))
#define DRAW_MARKER(m) \
  glLineWidth(m##_wt); \
  GLCOLOR3(m##_r, m##_g, m##_b); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x-m##_size),PY(m##_y));GLVERTEX2(PX(m##_x+m##_size),PY(m##_y));glEnd(); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x),PY(m##_y-m##_size));GLVERTEX2(PX(m##_x),PY(m##_y+m##_size));glEnd(); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x-m##_size),PY(m##_y));GLVERTEX2(PX(m##_x+m##_size),PY(m##_y));glEnd(); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x),PY(m##_y-m##_size));GLVERTEX2(PX(m##_x),PY(m##_y+m##_size));glEnd()

static void Draw_markers (STR *S) 
{
  DEVICE_CONST(int, nabs);
  DEVICE_CONST(int, nord); 
  #define _(Type,Name,Dflt,miN,maX) DEVICE_PAR(Type,Name);
  #include "ezpaint.h" /* dynamically linked parameters */
  #undef _
  if (marker1_size) {
    DRAW_MARKER(marker1);
  }
  if (marker2_size) {
    DRAW_MARKER(marker2);
  }
}
/* ========================================================================= */

#endif
