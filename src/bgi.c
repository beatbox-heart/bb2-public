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

/* Definitions for legacy Borland Graphics Interface emulator */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef X11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#endif

#include "system.h"
#include "k_.h"
#include "bikt.h"
#include "bgi.h"
#include "beatbox.h"

INT graphon=0;	/* graphics is on */
INT online=0;	/* show every primitive as it is made */

INT XMAX=639;	/* horizontal windoe size */
INT YMAX=479;	/* vertica window size */
INT WINX=-20;	/* user defined window location, in pixels, from screen left */
INT WINY=20;	/* user defined window location, in pixels, from screen top  */
		/* if either WINX or WINY < 0 then window       */
		/* location is measured from screen right/bottom */

char WindowName[1024];
char IconName[1024];

#if MPI || defined(NOX11)
/* make stubs */
static char NULLSTR[]="SORRY: any graphics package is absent in this version";
int Opengraph(void) {return 0;}
int _Cdecl EGAVGA_driver_far[1];
char *grapherrormsg(int errorcode) {return(NULLSTR);}
int graphresult(void) {return(grNotDetected);}
int registerfarbgidriver(void *driver) {return(0);}
void closegraph(void) {;}
void initgraph(int *gd, int *gm,const char *path) {;}
int getmaxx(void) {return(0);}
int getmaxy(void) {return(0);}
int getmaxcolor(void) {return(0);}
void setcolor(int color) {;}
void setwritemode(int mode) {;}
void settextstyle(int font, int direction, int charsize) {;}
void settextjustify(int horiz, int vert) {;}
void gettextsettings(struct textsettingstype  *texttypeinfo) {;}
int textwidth(char *text) {return 0;}
int textheight(char  *text) {return 0;}
void outtextxy(int x, int y,char  *textstring) {printf("%s\n", textstring);}
void putpixel(int x, int y, int color) {;}
void line(int x0, int y0, int x1, int y1) {;}
void rectangle(int lft, int top, int rgt, int bot) {;}
void setfillstyle(int pattern, int color) {;}
void bar(int lft, int top, int rgt, int bot) {;}
void fillellipse(int x, int y, int xradius, int yradius) {;}
unsigned imagesize(int left, int top, int right, int bottom) {return(0);}
void getimage(int left, int top, int right, int bottom,void  *bitmap) {;}
void putimage(int left, int top, const void  *bitmap, int op) {;}
void update_graph(void){;}
void dump_window(char *file, char *fmt){;}
void moveto(int x, int y) {;}
void lineto(int x, int y) {;}

#else

extern int _argc;
extern char **_argv;
/*extern char *VERSION;*/



/* #define DEGREE (1.0) */
#define CMAX 15	/* maximal color number  - not needed ?*/
#define BORDER       0 /*10*/      /* border width in pixels                       */
#define WIDTH ((int)XMAX+3+2*BORDER)
#define HEIGHT ((int)YMAX+3+2*BORDER)

#define BACKGROUND   0 /* background color (0 BLACK or 1 WHITE) */

/*#  define EV_MASK (KeyPressMask | EnterWindowMask)*/
/*#  define EV_MASK (KeyPressMask)*/
#define EV_MASK (0)


static Display       *theDisplay; 
static int            theDWidth;
static int            theDHeight;
static Window         theWindow;
static Pixmap         thePixmap;
static GC             theGC;
static Colormap       theColormap;
static int            theScreen; 
static int            theDepth;
static XFontStruct   *theFont=NULL;
static XTextProperty theWindowName;
static XTextProperty theIconName;
static unsigned long  theColors[CMAX+1];
static char BGIERRSTR[1024]="";
static int GRAPHRESULT=grOk;
static struct textsettingstype CURRENTTEXT;
static int CURRENTDRAWCOLOR=WHITE;
static int CURRENTFILLCOLOR=WHITE;
static int CURRENTFILLPATTERN=FillSolid;
int _Cdecl EGAVGA_driver_far[1]; /* dummy */


static int winx, winy; /* actual window location */
static int horspace, vertspace; /* free space in horizontal/vertical direction */

			    static int debug=1;
			    #define DB if (debug) fprintf(stderr,
			    #define BD ); FFLUSH(stderr);



int Opengraph(void) {
  int gdriver = VGA, gmode=VGAHI, errorcode;
  if (graphon) return 1;
  errorcode = registerfarbgidriver(EGAVGA_driver_far);
  if (errorcode < 0) {
     printf("Graphics error: %s\n", grapherrormsg(errorcode));
     /* printf("Press any key:"); */
     /* getch(); */
     return(0);
  }
  initgraph(&gdriver, &gmode, "");
  errorcode = graphresult();
  if (errorcode != grOk) {
     printf("Graphics error: %s\n", grapherrormsg(errorcode));
     /* printf("Press any key:"); */
     /* getch(); */
     return(0);
  }
  graphon=1;
  return 1;
}

void initgraph(int *gd,int *gm,const char *__pathtodriver) {
  /*float xsize, ysize;*/
  char *theDisplayName = NULL;
  XSizeHints    theSizeHints;
  char *iconname=IconName;
  char *windowname=WindowName;
  unsigned long background;
  XColor theRGBColor, theHardwareColor;
  /*int    theStatus;*/
  /*int    cv,i;*/
  if( (theDisplay = XOpenDisplay(NULL)) == NULL) {
    GRAPHRESULT=grNoInitGraph;
    snprintf(BGIERRSTR,1024,"Could not open a connection to X on display %s",XDisplayName(theDisplayName));
    fputs(BGIERRSTR,stderr);
    return;
  }
  theScreen     = DefaultScreen(theDisplay);
  theDepth      = DefaultDepth (theDisplay, theScreen);
  theDWidth     = DisplayWidth (theDisplay, theScreen);
  theDHeight    = DisplayHeight(theDisplay, theScreen);
  horspace=theDWidth-WIDTH;
  vertspace=theDHeight-HEIGHT;

  theColormap   = DefaultColormap(theDisplay, theScreen);
  background = BACKGROUND ? WhitePixel(theDisplay,theScreen) : BlackPixel(theDisplay,theScreen);
				/*DB "pos=%d,%d size=%d,%d\n",WINX,WINY,WIDTH,HEIGHT BD*/
  winx=WINX;
  if (horspace) {
    while (winx<0) winx+=horspace;
    while (winx>=horspace) winx-=horspace;
  }
  winy=WINY;
  if (vertspace) {
    while (winy<0) winy+=vertspace;	
    while (winy>=vertspace) winy-=vertspace;
  }
  theWindow = XCreateSimpleWindow(
    theDisplay,
    RootWindow(theDisplay, theScreen),
    winx,winy,WIDTH,HEIGHT,BORDER,
    background, background
  );
  XStringListToTextProperty(&windowname,1,&theWindowName);
  XStringListToTextProperty(&iconname,1,&theIconName);
  /* theSizeHints.flags = PSize | ((WINX>0 && WINY>0)?USPosition:PPosition); */
  theSizeHints.flags = PSize|USPosition;
  XSetWMProperties(theDisplay,theWindow,&theWindowName,&theIconName,NULL,0,&theSizeHints,NULL,NULL);
  XSelectInput(theDisplay, theWindow, EV_MASK);
  theGC = XCreateGC(theDisplay, theWindow, (unsigned long) 0, NULL);
  if(theGC==0) {
    XDestroyWindow(theDisplay,theWindow);
    GRAPHRESULT=grNotDetected;
    snprintf(BGIERRSTR,1024,"Could not create the Graphics Context");
    fputs(BGIERRSTR,stderr);
    return;
  }
  XMapWindow(theDisplay, theWindow); 
  XFlush(theDisplay);
  if NOT(thePixmap=XCreatePixmap(theDisplay,theWindow,(int)XMAX+3,(int)YMAX+3,theDepth))  
    printf("couldn't open Pixmap\n");
    
  if(theDepth==1) {
    printf("Warning: X window B&W\n");
  } else {
    #define _(i,name) if(XAllocNamedColor(theDisplay,theColormap,name,&theRGBColor,&theHardwareColor))\
    theColors[i]=theHardwareColor.pixel;\
    else printf("WARNING: cannot allocate color %s\n",name);
    /* color names chosen with pgm /usr/sbin/colorview */    
    _(BLACK,"black")	    _(BLUE,"darkblue")		_(GREEN,"darkgreen")	_(CYAN,"darkcyan")
    _(RED,"darkred")	    _(MAGENTA,"darkmagenta")	_(BROWN,"brown")	_(DARKGRAY,"gray33")
    _(LIGHTGRAY,"gray67")   _(LIGHTBLUE,"blue")		_(LIGHTGREEN,"green")   _(LIGHTCYAN,"cyan")
    _(LIGHTRED,"red")	    _(LIGHTMAGENTA,"magenta")	_(YELLOW,"yellow")	_(WHITE,"white")
    #undef _
  }

  /* The X graphics is open; can do the bgi rites now */
  graphon=1;
  settextstyle(DEFAULT_FONT,HORIZ_DIR,1);
  settextjustify(LEFT_TEXT,TOP_TEXT);
  setfillstyle(SOLID_FILL,BLACK);
  bar(0,0,(int)XMAX,(int)YMAX);		    /* clean out garbage */
  setfillstyle(SOLID_FILL,WHITE);
  setcolor(WHITE); 
  rectangle(-1,-1,(int)XMAX+1,(int)YMAX+1);   /* out window border */
  *gd=VGA; *gm=VGAHI;
  if (online) update_graph();
}
char * grapherrormsg(int errorcode) {return(BGIERRSTR);}
int  graphresult(void) {return(GRAPHRESULT);}
int   _Cdecl registerfarbgidriver(void *__driver) {return(0);}
void closegraph(void) {  
  if (!graphon) return;
  update_graph();
  /*getch();*/
  /*sleep(10);*/
  XFreeGC(theDisplay,theGC);
  XFreePixmap(theDisplay,thePixmap);
  XCloseDisplay(theDisplay);
  XFree(theWindowName.value);
  XFree(theIconName.value);
}
int getmaxx(void) {return((int)XMAX);}
int getmaxy(void) {return((int)YMAX);}
int getmaxcolor(void) {return(CMAX);}
void setcolor(int color) {CURRENTDRAWCOLOR=color;}
void  setwritemode(int mode) {;}	        /* can't find anything appropriate */
static int __lastx, __lasty;
void line(int x0, int y0, int x1, int y1) {
  XPoint thePoints[2];
  if (!graphon) return;
  __lastx=x1; __lasty=y1;
  x0++; y0++; x1++; y1++;
  thePoints[0].x=x0; thePoints[1].x=x1;
  thePoints[0].y=y0; thePoints[1].y=y1;
  XSetForeground(theDisplay, theGC,theColors[CURRENTDRAWCOLOR]);
  XDrawLines(theDisplay,thePixmap,theGC,thePoints,2,CoordModeOrigin);
  if (online) update_graph();
}
void moveto(int x, int y) {
  if (!graphon) return;
				    /*DB "moveto %d %d\n",x,y BD*/
  __lastx=x; __lasty=y;
  if (online) update_graph();
}
void lineto(int x, int y) {
  if (!graphon) return;
				    /*DB "lineto %d %d\n",x,y BD*/
  line(__lastx,__lasty,x,y);
  if (online) update_graph();
}
void rectangle(int lft, int top, int rgt, int bot) {
  XPoint thePoints[5];
  if (!graphon) return;
  lft++; top++; rgt++; bot++;
  thePoints[0].x=lft; thePoints[1].x=lft; thePoints[2].x=rgt; thePoints[3].x=rgt; thePoints[4].x=lft;
  thePoints[0].y=top; thePoints[1].y=bot; thePoints[2].y=bot; thePoints[3].y=top; thePoints[4].y=top;
  XSetForeground(theDisplay,theGC,theColors[CURRENTDRAWCOLOR]);
  XDrawLines(theDisplay,thePixmap,theGC,thePoints,5 ,CoordModeOrigin);
  if (online) update_graph();
}
void setfillstyle(int pattern, int color) {
  if (!graphon) return;
  CURRENTFILLCOLOR=color;
  switch (pattern) {
  case EMPTY_FILL:	CURRENTFILLPATTERN=FillStippled; break;
  case SOLID_FILL:	CURRENTFILLPATTERN=FillSolid; break;
  case LINE_FILL:	CURRENTFILLPATTERN=FillStippled; break;
  case LTSLASH_FILL:	CURRENTFILLPATTERN=FillStippled; break;
  case SLASH_FILL:	CURRENTFILLPATTERN=FillStippled; break;
  case BKSLASH_FILL:	CURRENTFILLPATTERN=FillStippled; break;
  case LTBKSLASH_FILL:	CURRENTFILLPATTERN=FillStippled; break;
  case HATCH_FILL:	CURRENTFILLPATTERN=FillTiled; break;
  case XHATCH_FILL: 	CURRENTFILLPATTERN=FillTiled; break;
  case INTERLEAVE_FILL:	CURRENTFILLPATTERN=FillTiled; break;
  case WIDE_DOT_FILL:	CURRENTFILLPATTERN=FillTiled; break;
  case CLOSE_DOT_FILL: 	CURRENTFILLPATTERN=FillTiled; break;
  }
}
void bar(int lft,int top,int rgt,int bot) {
  if (!graphon) return;
  lft++; top++; rgt++; bot++;
  XSetForeground(theDisplay,theGC,theColors[CURRENTFILLCOLOR]);
  XSetFillStyle(theDisplay,theGC,CURRENTFILLPATTERN);
  XFillRectangle(theDisplay,thePixmap,theGC,lft,top,rgt-lft+1,bot-top+1);
  if (online) update_graph();
}
void fillellipse(int x, int y, int xradius, int yradius) {
  if (!graphon) return;
				    /*DB "fillellipse %d %d %d %d\n",x,y,xradius,yradius BD*/
  x++; y++;
  XSetForeground(theDisplay,theGC,theColors[CURRENTFILLCOLOR]);
  XSetFillStyle(theDisplay,theGC,CURRENTFILLPATTERN);
  XFillArc(theDisplay,thePixmap,theGC,x,y,2*xradius+1,2*yradius+1,0,360*64);
  if (online) update_graph();
}
void putpixel(int x, int y, int color) {
  if (!graphon) return;
  x++; y++;
  color=color%(CMAX+1);
  XSetForeground(theDisplay, theGC,theColors[color]);
  XSetFillStyle(theDisplay,theGC,FillSolid);
  XDrawPoint(theDisplay,thePixmap,theGC,x,y);
  if (online) update_graph();
}
void  settextstyle(int font, int direction, int charsize) {
  if (!graphon) return;
  char fontname[1024];
  if (!charsize) charsize=1;
  if (direction!=HORIZ_DIR) {
    GRAPHRESULT=grFontNotFound;
    snprintf(BGIERRSTR,1024,"font direction must be horizontal in this implementation\n");
    fputs(BGIERRSTR,stderr);
    return;
  }
  CURRENTTEXT.font=font;
  CURRENTTEXT.direction=direction;
  CURRENTTEXT.charsize=charsize;
  if (theFont) {XFreeFont(theDisplay,theFont);  theFont=NULL;}
  snprintf(fontname,1024,"-adobe-courier-medium-r-normal--%u-*",CURRENTTEXT.charsize*12);
  if NOT(theFont=XLoadQueryFont(theDisplay,fontname)){
    snprintf(fontname,1024,"*--%u-*",CURRENTTEXT.charsize*8); 
    if NOT(theFont=XLoadQueryFont(theDisplay,fontname)){
      GRAPHRESULT=grFontNotFound;
      snprintf(BGIERRSTR,1024,"font %s not found",fontname);
      fputs(BGIERRSTR,stderr);
      return;
    }
  }
  XSetFont(theDisplay,theGC,theFont->fid);
}
void  settextjustify(int horiz, int vert) {
  if (!graphon) return;
  CURRENTTEXT.horiz=horiz;	    
  CURRENTTEXT.vert=vert;
}
void gettextsettings(struct textsettingstype  *texttypeinfo) {
  if (!graphon) return;
  if (!theFont) {fprintf(stderr,"gettextsettings: font not defined\n");return;}
  memcpy(texttypeinfo, &CURRENTTEXT, sizeof(struct textsettingstype));
}
void outtextxy(int x, int y,char  *text) {
  int dir, asc, desc;
  XCharStruct cs;
  if (!graphon) return;
  x++; y++;
  if (!theFont) {fprintf(stderr,"outtextxy: font not defined\n");return;}
  XSetForeground(theDisplay,theGC,theColors[CURRENTDRAWCOLOR]);
  XTextExtents(theFont,text,(int)strlen(text),&dir,&asc,&desc,&cs);
  switch(CURRENTTEXT.horiz) {
    case LEFT_TEXT:   x-=cs.lbearing; break;
    case CENTER_TEXT: x-=(cs.lbearing+cs.rbearing)/2; break;
    case RIGHT_TEXT:  x-=cs.rbearing; break;
  }
  switch(CURRENTTEXT.vert) {
    case TOP_TEXT:    y+=theFont->ascent; break;
    case CENTER_TEXT: y+=(theFont->ascent-theFont->descent)/2; break;
    case BOTTOM_TEXT: y-=theFont->descent; break;
  }
  XDrawString(theDisplay,thePixmap,theGC,x,y,text,(int)strlen(text));
  if (online) update_graph();
}
int textwidth(char *text) {
  if (!graphon) return 0;
  if (!theFont) {fprintf(stderr,"textwidth: font not defined\n");return 0;}
  return XTextWidth(theFont,text,(int)strlen(text));
 }
int textheight(char  *text) {
  if (!graphon) return 0;
  if (!theFont) {fprintf(stderr,"textheight: font not defined\n");return 0;}
  return theFont->ascent+theFont->descent+1;
}
typedef struct {int w,h; Pixmap m;} bitmapstruct;
unsigned imagesize(int left, int top, int right, int bottom) {
  if (!graphon) return 0;
  return(sizeof(bitmapstruct));
}
void getimage(int lft, int top, int rgt, int bot,void  *bitmap) {
  if (!graphon) return;
  bitmapstruct s;
  lft++; top++; rgt++; bot++;
  s.w=rgt-lft+1;
  s.h=bot-top+1;
  s.m=XCreatePixmap(theDisplay,thePixmap,s.w,s.h,theDepth);
  XCopyArea(theDisplay,thePixmap,s.m,theGC,lft,top,s.w,s.h,0,0);
  memcpy(bitmap,&s,sizeof(bitmapstruct));
}
void putimage(int lft, int top, const void  *s, int op) {
  if (!graphon) return;
  bitmapstruct *S=(bitmapstruct *)s;
  lft++; top++;
  XCopyArea(theDisplay,S->m,thePixmap,theGC,0,0,S->w,S->h,lft,top);
  if (online) update_graph();
}
void cleardevice(void){
  if (!graphon) return;
  int fill=CURRENTFILLPATTERN;
  int col=CURRENTFILLCOLOR;
  setfillstyle(SOLID_FILL,BLACK);
				update_graph();
				sleep(1);
  bar(0,0,(int)XMAX,(int)YMAX);
  setfillstyle(fill,col);
  if (online) update_graph();
}
void update_graph(void){
  if (!graphon) return;
  XCopyArea(theDisplay,thePixmap,theWindow,theGC,0,0,(int)XMAX+3,(int)YMAX+3,0,0);
  XFlush(theDisplay);
}
void dump_window(char *filename,char *fmt){
  char cmd[4096];
  if (!graphon) return;
  if (!theWindow) return;
  XRaiseWindow(theDisplay,theWindow);
  XCopyArea(theDisplay,thePixmap,theWindow,theGC,0,0,(int)XMAX+3,(int)YMAX+3,0,0);
  XFlush(theDisplay);
  if (0==stricmp(fmt,"ppm"))
    snprintf(cmd,4096,"xwd -id %lu -silent | xwdtopnm -quiet > %s",theWindow,filename);
  else if (0==stricmp(fmt,"gif"))
    snprintf(cmd,4096,"xwd -id %lu -silent | xwdtopnm -quiet | ppmtogif > %s",theWindow,filename);
  else if (0==stricmp(fmt,"jpeg") || 0==stricmp(fmt,"jpg"))
    snprintf(cmd,4096,"xwd -id %lu -silent | xwdtopnm -quiet | pnmtojpeg --optimize > %s",theWindow,filename);
  else if (0==stricmp(fmt,"png"))
    snprintf(cmd,4096,"xwd -id %lu -silent | xwdtopnm -quiet | pnmtopng > %s",theWindow,filename);
  else { MESSAGE("dump_window: unknown format %s",fmt); return; }
  system(cmd);
}
void setlinewidth(int width){
  XGCValues theGCValues;
  if (!graphon) return;
  XGetGCValues(theDisplay,theGC,GCLineStyle|GCCapStyle|GCJoinStyle,&theGCValues);
  XSetLineAttributes(theDisplay,theGC,width,
		     theGCValues.line_style,theGCValues.cap_style,theGCValues.join_style);
}

/********************************************************************/
#endif
