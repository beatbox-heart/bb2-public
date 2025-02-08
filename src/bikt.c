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

/* SOME MY AUXILIARY PROGRAMS */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "system.h"
#include "bgi.h"
#include "beatbox.h"
#define OWN
#include "bikt.h"

#ifdef IBMPC
char *just_fid(char *path) {
  static char drive[MAXDRIVE];
  static char dir[MAXDIR];
  static char file[MAXFILE];
  static char ext[MAXEXT];
  static char fid[MAXFILE+MAXEXT+2];
  fnsplit(path,drive,dir,file,ext);
  snprintf(fid,MAXFILE+MAXEXT+2,"%s%s",file,ext);
  return(fid);
}
#endif

#define MAXKBSTACK 128
static char kbstack[MAXKBSTACK];
static int kbstacktop=0;
int kbhitmy(void)  {return(kbstacktop? 1                   :kbhit());}
int getchmy(void) {return(kbstacktop?kbstack[--kbstacktop]:getch());}
int ungetchmy(int ch) {
  if(kbstacktop>=MAXKBSTACK) return(EOF);
  return(kbstack[kbstacktop++]=ch);
}
void kbclearmy(void) {kbstacktop=0;kbclear();}

void beep(unsigned Hz, unsigned ms) {
  sound(Hz); delay(ms); nosound();
}

#ifdef IBMPC
#define CHK_GRAPH							\
  if (!graphon) {                                                       \
    MESSAGE("%-79s\r",text);                                             \
    return;                                                             \
  }
#else
#define CHK_GRAPH							\
  if (!graphon) {                                                       \
    MESSAGE("%-79s\n",text);                                             \
    return;                                                             \
  }
#endif

void crt_text(char *text, int row, int col, int color) {
  struct textsettingstype old;
  int left, top, bottom, right;
  CHK_GRAPH;
  gettextsettings(&old);
  setcolor(MAGENTA);
  settextstyle(DEFAULT_FONT, HORIZ_DIR, 1);
  settextjustify(LEFT_TEXT,TOP_TEXT);
  left = col*textwidth("M");
  top = row*textheight("Q");
  right = left+textwidth(text);
  bottom = top+textheight(text);
      /*	printf("width=%d height=%d\n",textwidth("M"),textheight("Q")); */
	/* printf("box %d %d %d %d\n",left,top,bottom,right); */
  setfillstyle (SOLID_FILL,color/16);
  bar (left,top,right,bottom);
  setcolor (color%16);
  outtextxy (left,top,text);
  settextstyle (old.font,old.direction,old.charsize);
  settextjustify (old.horiz,old.vert);
}

#define MAX_SIZE 4096
static char next_text_bm[MAX_SIZE];
static int old_left, old_top;
static int cur_row, cur_col, cur_color;
static char cur_text[256];

void next_text(char *text, int row, int col, int color) {
  struct textsettingstype old;
  int left, top, bottom, right;
  int size;

  CHK_GRAPH;
  cur_row=row; cur_col=col; cur_color=color;
  strncpy(cur_text,text,256);
  hide_text();
  if (text[0]=='\0') return;
  gettextsettings(&old);
  settextstyle(DEFAULT_FONT, HORIZ_DIR, 1);
  settextjustify(LEFT_TEXT,TOP_TEXT);
  left = col*textwidth("M");
  top = row*textheight("Q");
  right = left+textwidth(text);
  bottom = top+textheight(text);
  size = imagesize(left,top,right,bottom);
  if (size<MAX_SIZE) {
    getimage(left,top,right,bottom,next_text_bm);
    old_left=left;
    old_top=top;
  }
  setfillstyle (SOLID_FILL,color/16);
  bar (left,top,right,bottom);
  setcolor (color%16);
  outtextxy (left,top,text);

  settextstyle (old.font,old.direction,old.charsize);
  settextjustify (old.horiz,old.vert);
}

void hide_text(void) {
  if (!graphon) return;
  if (next_text_bm[0]) {
    putimage(old_left,old_top,next_text_bm,COPY_PUT);
    next_text_bm[0]=0;
  }
}

void restore_text(void) {
  next_text(cur_text,cur_row,cur_col,cur_color);
}

char stop_text(char *text, int row, int col, int color) {
  char ch;
  next_text(text,row,col,color);
  ch=getchmy();
  next_text("",0,0,0);
  return(ch);
}
#undef MAX_SIZE

int fexist (const char *fid) {
  FILE *f=fopen(fid,"r");
  if (f) fclose(f);
  return(f!=NULL);
}

#ifdef IBMPC
void fputext(const char *oldfid, char *newfid, const char *newext) {
  char drive[MAXDRIVE], dir[MAXDIR], name[MAXFILE], ext[MAXEXT];
  fnsplit(oldfid,drive,dir,name,ext);
  snprintf(newfid,MAXPATH,"%s%s%s.%-3s",drive,dir,name,newext);
}
#else
void fputext(const char *oldfid, char *newfid, const char *newext) {
  char *p;
  p = strrchr(oldfid,'.');
  if (p) {
    strncpy(newfid,oldfid,p-oldfid);
    p = newfid + (p-oldfid);
    *(p++)='.';
    strcpy(p,newext);
  } else {
    snprintf(newfid,MAXPATH,"%s.%s",oldfid,newext);
  }
}
#endif

int frename(const char *oldfid, const char *newfid) {
  int retcode;
  if (fexist(newfid)) {
    char bakname[MAXPATH];
    fputext(newfid,bakname,"BAK");
    remove(bakname);
    rename(newfid,bakname);
  }
  retcode = rename(oldfid,newfid);
  return retcode;
}


