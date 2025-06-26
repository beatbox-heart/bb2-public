/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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
 * Extra facility: put one or two markers on the screen (as in ezride).
 * Differs from ezpaint in that
 * - it only does grid values option but not k-functions,
 * - it work in MPI mode. 
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
  NOGRAPH_DUMMY(mpipaint)
#else
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h> 
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <GL/gl.h>
#include <GL/glx.h>

/* Macros surviving from EZSPIRAL */
#define TRUE 1   
#define MAXWINDOWTITLE 512
#define Xmin (-1.0)
#define Ymin (-1.0)
#define Xmax (1.0)
#define Ymax (1.0)
#define PX(x) (Xmin+((x)*(Xmax-Xmin))/nabs)
#define PY(y) (Ymin+((y)*(Ymax-Ymin))/nord)

/* Imported global variables */
extern int Verbose;		/* defined in beatbox.c */
extern int idev;		/* defined in beatbox.c */
extern unsigned char *gpoints;	/* defined in geometry.c */

typedef enum {xy, yx, xz, zx, yz, zy} aspect_t;

typedef struct {
  /* The raster sizes */
  int nabs;			/* horizontal */
  int nord;			/* vertical */
  int napp;			/* depth  */

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
  /* precomputed combinations thereof */
  real rbase, gbase, bbase;
  real rcoef, gcoef, bcoef;
  int r, g, b;

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
  #include "mpipaint.h" /* dynamically linked parameters */
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

  /* X11 and GL stuff */
  int winx, winy, width, height; /* window coordinates and sizes */
  int doublebuffer;		/* flag of whether to use GL double-buffering */
  Display *theDisplay;		
  Window theWindow;
  GLXContext theGLXContext;
  XTextProperty theWindowName, theIconName;
  XEvent event;

  /* MPI stuff */
  #if MPI
  MPI_Comm comm;		/* communicator for this device instance */
  int tagbase;			/* tag base for messages of this device instance */
  MPI_Datatype *theirType;	/* array of receive buffers in the output thread */
  MPI_Datatype myType;		/* send buffer for a work thread */
  MPI_Datatype mygType;		/* .., for the tissue data */
  real *rfield;			/* copies                 */
  real *gfield;			/*     of relevant fields */
  real *bfield;			/*     in the root thread */
  real *tfield;			/* .., of the geometry ("tissue") */
  #endif
  int root;			/* the rank of the thread that does the output */
  /*					- used in seq mode for formal purpose, too */
} STR;

#if MPI
/* Tagging the messages, for extra care */
enum {
  tag_l, /* limits */
  tag_r, /* red */
  tag_g, /* green */
  tag_b, /* blue */
  tag_t, /* tissue */
  tag_n	 /* # of tag variants */
};
#define tag(a) (tagbase*tag_n+tag_##a)
#endif

#define CROP(a,b) a=(b); if(a<0) a=0; if (a>1) a=1;

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


RUN_HEAD(mpipaint)
{
#if MPI
  DEVICE_CONST(MPI_Comm,comm);
  DEVICE_CONST(int,tagbase);
  DEVICE_CONST(int,root);
  DEVICE_ARRAY(MPI_Datatype,theirType);
  DEVICE_CONST(MPI_Datatype,myType);
  DEVICE_CONST(MPI_Datatype,mygType);
  DEVICE_CONST(int,rlayer);
  DEVICE_CONST(int,glayer);
  DEVICE_CONST(int,blayer);
  DEVICE_CONST(real *,rfield);
  DEVICE_CONST(real *,gfield);
  DEVICE_CONST(real *,bfield);
  DEVICE_CONST(real *,tfield);
  int rank;  
  
  if (mpi_rank != root) {
    /* Send my data to root */
#define _(a) if (a##layer>=0) \
      MPIDO(MPI_Send(New+DV*a##layer,1,myType,root,tag(a),comm), \
	    "mpipaint "#a"layer send failure")
    _(r); _(g); _(b);
#undef _
    if (geom_vmax) 
      MPIDO(MPI_Send(Geom,1,mygType,root,tag(t),comm),"mpipaint geometry send failure");
    return SUCCESS;
  }
  
  /* Receive data from all other threads */
  for (rank=0; rank<num_active_procs; rank++) if (rank!=root) {
#define _(a) if (a##layer>=0) \
	MPIDO(MPI_Recv(a##field,1,theirType[rank],rank,tag(a),comm,MPI_STATUS_IGNORE), \
	      "mpipaint "#a"layer receive from #%d failure",rank)
    _(r); _(g); _(b);
#undef _
    if (geom_vmax) MPIDO(MPI_Recv(tfield,1,theirType[rank],rank,tag(t),comm,MPI_STATUS_IGNORE),
	      "mpipaint geometry receive from #%d failure",rank);
  } /* for rank */
  
  /* Send/receive data to/from myself */
#define _(a) if (a##layer>=0) 						\
    MPIDO(MPI_Sendrecv(New+DV*a##layer, 1, myType, root, tag(a), 	\
		       a##field, 1, theirType[root], root, tag(a),	\
		       comm, MPI_STATUS_IGNORE),			\
	  "mpipaint "#a"layer sendrecv failure")
  _(r); _(g); _(b);

  if (geom_vmax) 							\
    MPIDO(MPI_Sendrecv(Geom, 1, mygType, root, tag(t),			\
		       tfield, 1, theirType[root], root, tag(t),	\
		       comm, MPI_STATUS_IGNORE),			\
	  "mpipaint geometry sendrecv failure");
#undef _
  
  /* Redefine access method for the fields */
#define VALUE(x,y,z,v) (((v)>=0)?(v##field[((x)*ymax+(y))*zmax+(z)]):0)
#else /* not MPI */
  /* The standard (sequential) access method for the fields */
#define VALUE(x,y,z,v) (((v)>=0)?(New[ind(x,y,z,v)]):0)
#endif /* not MPI */

  /*************************************/
  /* The graph work proper starts here */
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
  
  DEVICE_ARRAY(char, filter);
  DEVICE_ARRAY(char, filtercode);
  DEVICE_CONST(pp_fn, filtercompiled);
  
  DEVICE_CONST(aspect_t,aspectchoice);
  DEVICE_CONST(int,nabs)   DEVICE_CONST(int,nord)   DEVICE_CONST(int,napp);
  DEVICE_CONST(int,abs0)   DEVICE_CONST(int,ord0)   DEVICE_CONST(int,app0);
  DEVICE_CONST(int,abs1)   DEVICE_CONST(int,ord1)   DEVICE_CONST(int,app1);
  DEVICE_CONST(real,rbase) DEVICE_CONST(real,rcoef) DEVICE_CONST(int,r);
  DEVICE_CONST(real,gbase) DEVICE_CONST(real,gcoef) DEVICE_CONST(int,g);
  DEVICE_CONST(real,bbase) DEVICE_CONST(real,bcoef) DEVICE_CONST(int,b);
  DEVICE_CONST(real,bg_r)  DEVICE_CONST(real,bg_g)  DEVICE_CONST(real,bg_b);
  
  DEVICE_CONST(pp_fn, xtipcompiled);
  DEVICE_CONST(pp_fn, ytipcompiled);
  DEVICE_VAR(int, ntip);
  DEVICE_ARRAY(double, tipx);
  DEVICE_ARRAY(double, tipy);
  DEVICE_CONST(int, ntipmax);
  
  GLfloat X1, X2, Y1, Y2;
  real red, grn, blu;
  real rsum, gsum, bsum;
  int iapp, iord, iabs;
  int x, y, z;
  int nval;
  int itip;
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
	} /* switch (aspectchoice) */
	/* NB in MPI we are dealing with copies of global data, and */
	/* isTIssue appeals to Geom which has access only to local geom info */
#if MPI
	if ( tfield?(tfield[((x)*ymax+(y))*zmax+(z)]):1 )
#else
	if (isTissue(x,y,z))
#endif
	  {
	  rsum+=VALUE(x,y,z,r);
	  gsum+=VALUE(x,y,z,g);
	  bsum+=VALUE(x,y,z,b);
	  nval++;
	} /* if gpoints/isTissue */
#undef VALUE
      } /* for iapp */
      if (nval) {
	CROP(red,rbase+rsum/nval*rcoef);
	CROP(grn,gbase+gsum/nval*gcoef);
	CROP(blu,bbase+bsum/nval*bcoef);
      } else { /*  not nval: void, use background colour. */
	red=bg_r;
	grn=bg_g;
	blu=bg_b;
      } /* not nval */
      glColor3f((GLfloat)(red),(GLfloat)(grn),(GLfloat)(blu));
      glRectf(X1,Y1,X2,Y2);
    } /* for iord */
  } /* for iabs */
  
    /* Construct the tip trace if needed */
  if ( (xtipcompiled!=NULL) && (ytipcompiled!=NULL) ) {
    /* if stack is full, discard the earliest record */
    if (*ntip==ntipmax) {
      for (itip=0;itip<(*ntip)-1;itip++) {
	tipx[itip]=tipx[itip+1];
	tipy[itip]=tipy[itip+1];
      } /* for itip */
      (*ntip)--;
    } /* if (*ntip==ntipmax) */
    
    /* x and y codes to be in the same units (grid coords) as iabs and iord */
    k_on();
    if (xtipcompiled!=NULL) 
      memcpy(&the_xtip,execute(xtipcompiled),sizeof(REAL));
    if (ytipcompiled!=NULL) 
      memcpy(&the_ytip,execute(ytipcompiled),sizeof(REAL));
    if (xtipcompiled!=NULL && ytipcompiled!=NULL && the_xtip!=real_inf && the_ytip !=real_inf) {
      tipx[*ntip]=the_xtip;
      tipy[*ntip]=the_ytip;
      (*ntip)++;
    } /* if xtipcompiled ... */
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
RUN_TAIL(mpipaint)

DESTROY_HEAD(mpipaint)
DESTROY_TAIL(mpipaint)

CREATE_HEAD(mpipaint)
{
#if MPI
  enum {X, Y, Z, V, NDIMS};
  enum {X0, X1, Y0, Y1, Z0, Z1, NLIMS};
  int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
  int buf[NLIMS];  
  int tagbase;
  int rank;
  real *dummy;  
  MPI_Comm comm;		/* communicator for this device instance */
#endif
  int root;
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
  
  DEVICE_IS_RECTANGULAR;
  
  ACCEPTSN(aspect,2,"xy");
  if (0) {}
#define CASE(t,a,o,p) else if (0==strcmp(aspect,#t)) { S->aspectchoice=t; \
    S->abs0=s.global_##a##0; S->ord0=s.global_##o##0; S->app0=s.global_##p##0; \
    S->abs1=s.global_##a##1; S->ord1=s.global_##o##1; S->app1=s.global_##p##1; }
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
  ASSERT(rlayer>=0 || glayer>=0 || blayer>=0);  
  ACCEPTR(bg_r,0,0,1);
  ACCEPTR(bg_g,0,0,1);
  ACCEPTR(bg_b,0,0,1);

  /* Precompute mapping coeffs based on user-def'd pars */
  S->rbase=(rlayer>=0)?(-rmin/(rmax-rmin)):0;
  S->gbase=(glayer>=0)?(-gmin/(gmax-gmin)):0;
  S->bbase=(blayer>=0)?(-bmin/(bmax-bmin)):0;
  S->rcoef=(rlayer>=0)?(1.0/(rmax-rmin)):0;
  S->gcoef=(glayer>=0)?(1.0/(gmax-gmin)):0;
  S->bcoef=(blayer>=0)?(1.0/(bmax-bmin)):0;
  S->r=(rlayer>=0)?rlayer:0; /* these will be */
  S->g=(glayer>=0)?glayer:0; /* used in seq mode; */
  S->b=(blayer>=0)?blayer:0; /* mpi mode has *dummy instead */
  
  ACCEPTS(xtip,"");
  ACCEPTS(ytip,"");
  if ( ('\0'!=*xtip) && ('\0'!=*ytip) ) {
    S->xtipcompiled=compile(S->xtip,deftb,t_real); CHK(S->xtip);
    S->ytipcompiled=compile(S->ytip,deftb,t_real); CHK(S->ytip);
    ACCEPTI(ntipmax,100,1,INONE);
    CALLOC(S->tipx,S->ntipmax,sizeof(double));
    CALLOC(S->tipy,S->ntipmax,sizeof(double));
    S->ntip=0;
  } else { /* not xtip && ytip */
    S->ntipmax=0;
    S->xtipcompiled=NULL;
    S->ytipcompiled=NULL;
    S->tipx=NULL;
    S->tipy=NULL;
  } /* not xtip && ytip */
  
/* #define acceptREAL(b,c,d,e) if (!acceptrk(#b"=",&(S->b##ptr),&(S->b),&(S->b##src),&(S->b##exe),deftb,c,d,e,w)) return(0); REAL b=S->b */
#define acceptREAL(b,c,d,e) ACCEPTRK(b,c,d,e) REAL b=S->b
/* #define acceptINT(b,c,d,e) if (!acceptik(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); INT b=S->b */
#define acceptINT(b,c,d,e) ACCEPTIK(b,c,d,e) INT b=S->b
#define _(Type,Name,Dflt,miN,maX) accept##Type(Name,Dflt,miN,maX);
#include "mpipaint.h"
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
  
  ACCEPTS(title,"mpipaint t=%.0f");
  ACCEPTS(titlecode,"t");
  S->titlecompiled=compile(S->titlecode,deftb,t_real); CHK(S->titlecode);
  
  /*********************************/
#if MPI
  /* Create private communicator for this instance of the device. */
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &comm, &root))
    EXPECTED_ERROR("Could not create communicator.\n");
  S->root=root;
  S->comm=comm;

  /* To make messages unique for this device */
  S->tagbase=tagbase=idev;
  
  /*******************************************/
  /* Define exchange buffers for ALL threads */

  /* Send my limits to the root to make receipt buffer types there */
  if (mpi_rank!=root) {
    buf[X0]=s.x0; buf[X1]=s.x1;
    buf[Y0]=s.y0; buf[Y1]=s.y1;
    buf[Z0]=s.z0; buf[Z1]=s.z1;
    MPIDO(MPI_Send(buf,NLIMS,MPI_INT,root,tag(l),comm),"Could not send subdomain sizes.");
  }

  /* Create my send buffer, whether I am root or not */
  sizes[X] = xlen; subsizes[X] = s.x1-s.x0+1; starts[X] = s.x0-local_xmin;
  sizes[Y] = ylen; subsizes[Y] = s.y1-s.y0+1; starts[Y] = s.y0-local_ymin;
  sizes[Z] = zlen; subsizes[Z] = s.z1-s.z0+1; starts[Z] = s.z0-local_zmin;
  sizes[V] = vmax; subsizes[V] = 1;    starts[V] = 0;
  MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&(S->myType)),
	"Couldn't define myType.");
  MPIDO(MPI_Type_commit(&(S->myType)),"Couldn't commit myType.");

  /* .., for the geometry */
  if (geom_vmax) {
    sizes[X] = xlen; subsizes[X] = s.x1-s.x0+1; starts[X] = s.x0-local_xmin;
    sizes[Y] = ylen; subsizes[Y] = s.y1-s.y0+1; starts[Y] = s.y0-local_ymin;
    sizes[Z] = zlen; subsizes[Z] = s.z1-s.z0+1; starts[Z] = s.z0-local_zmin;
    sizes[V] = geom_vmax; subsizes[V] = 1;    starts[V] = 0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&(S->mygType)),
	  "Couldn't define myType.");
    MPIDO(MPI_Type_commit(&(S->mygType)),"Couldn't commit myType.");
  }
  
  if (mpi_rank==root) {
    /* Array of receipt buffer types for ALL threads */
    CALLOC(S->theirType,num_active_procs,sizeof(MPI_Datatype));    

    /* Create receipt buffer types for OTHER threads */
    for (rank=0; rank<num_active_procs; rank++) if (rank!=root) {
      MPIDO(MPI_Recv(buf,NLIMS,MPI_INT,rank,tag(l),comm,MPI_STATUS_IGNORE),"Could not receive subdomain sizes.");
      sizes[X]=xmax; subsizes[X]=buf[X1]-buf[X0]+1; starts[X]=buf[X0];
      sizes[Y]=ymax; subsizes[Y]=buf[Y1]-buf[Y0]+1; starts[Y]=buf[Y0];
      sizes[Z]=zmax; subsizes[Z]=buf[Z1]-buf[Z0]+1; starts[Z]=buf[Z0];
      sizes[V]=1;    subsizes[V]=1;                 starts[V]=0;
      MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&(S->theirType[rank])),
	    "Couldn't define theirType.");
      MPIDO(MPI_Type_commit(&(S->theirType[rank])),"Couldn't commit theirType.");
    } /* for rank */
    
    /* And now the receipt buffer type for MYSELF, for Sendrecv */
    sizes[X]=xmax; subsizes[X]=s.x1-s.x0+1; starts[X]=s.x0;
    sizes[Y]=ymax; subsizes[Y]=s.y1-s.y0+1; starts[Y]=s.y0;
    sizes[Z]=zmax; subsizes[Z]=s.z1-s.z0+1; starts[Z]=s.z0;
    sizes[V]=1;    subsizes[V]=1;	    starts[V]=0;
    MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&(S->theirType[root])),
	  "Couldn't define theirType.");
    MPIDO(MPI_Type_commit(&(S->theirType[root])),"Couldn't commit theirType.");
    
    /* And the buffers themselves */
    if (rlayer>=0) CALLOC(S->rfield,xmax*ymax*zmax,sizeof(real));
    if (glayer>=0) CALLOC(S->gfield,xmax*ymax*zmax,sizeof(real));
    if (blayer>=0) CALLOC(S->bfield,xmax*ymax*zmax,sizeof(real));
    if (geom_vmax>0) CALLOC(S->tfield,xmax*ymax*zmax,sizeof(real));

    if (rlayer>=0) dummy=S->rfield;
    else if (glayer>=0) dummy=S->gfield;
    else if (blayer>=0) dummy=S->bfield;
    else ABORT("at least one layer is expected");
    
    if (rlayer<0) S->rfield=dummy;
    if (glayer<0) S->gfield=dummy;
    if (blayer<0) S->bfield=dummy;
  } /* if (mpi_rank==root) */
#else /* not MPI */
  root = mpi_rank; /* presumably = 0, defined in state.h */
#endif /* not MPI */
  
  if (mpi_rank==root) { /* in seq mode this is true by definition of root */
    /*********************************/
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
    if NOT(theColormap = XCreateColormap(S->theDisplay,
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
	    EXPECTED_ERROR("couldn't create X11 window");
    
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
      MESSAGE("%s version %d of the X11 Window System, X%d R%d\n",
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
  } /* if (mpi_rank==root) */
}
CREATE_TAIL(mpipaint,0);

/* ========================================================================= */

#define GLRECT(x1,y1,x2,y2)  glRectf((x1),(y1),(x2),(y2))
#define GLCOLOR3(r,g,b)      glColor3f((r),(g),(b))
#define GLVERTEX2(x,y)       glVertex2f((x),(y))
#define DRAW_MARKER(m)		 \
  glLineWidth(m##_wt);							\
  GLCOLOR3(m##_r, m##_g, m##_b);					\
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x-m##_size),PY(m##_y));GLVERTEX2(PX(m##_x+m##_size),PY(m##_y));glEnd(); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x),PY(m##_y-m##_size));GLVERTEX2(PX(m##_x),PY(m##_y+m##_size));glEnd(); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x-m##_size),PY(m##_y));GLVERTEX2(PX(m##_x+m##_size),PY(m##_y));glEnd(); \
  glBegin(TIP_PLOT_TYPE);GLVERTEX2(PX(m##_x),PY(m##_y-m##_size));GLVERTEX2(PX(m##_x),PY(m##_y+m##_size));glEnd()

static void Draw_markers (STR *S) 
{
  DEVICE_CONST(int, nabs);
  DEVICE_CONST(int, nord); 
#define _(Type,Name,Dflt,miN,maX) DEVICE_PAR(Type,Name);
#include "mpipaint.h" /* dynamically linked parameters */
#undef _
  if (marker1_size) {
    DRAW_MARKER(marker1);
  }
  if (marker2_size) {
    DRAW_MARKER(marker2);
  }
}
/* ========================================================================= */

#endif /* X11 and GL */
