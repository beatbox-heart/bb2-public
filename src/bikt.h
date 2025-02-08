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

#define kbclear() {while(kbhit()) getch();}
int kbhitmy(void);
int getchmy(void);
int ungetchmy(int ch);
void kbclearmy(void);
void beep(unsigned Hz, unsigned ms);

/* Using 0 in place of NULL avoids warnings from NOT(int foo = somethingThatReturnsPointer()) */
#define NOT(a) (0==(a))

/* Write text onscreen in "crt" style */
void crt_text(char *text, int row, int col, int color);

/* write text "upon" the screen, wait a key and then restore. Return key */
char stop_text(char *text, int row, int col, int color);

/* remove prev msg if any and write the new one */
void next_text(char *text, int row, int col, int color);

/* remove text written by next_text */
void hide_text(void);

/* restore text written by next_text and then hidden by hide_text */
void restore_text(void);

/* Print and return size of free space in heap */
#define memavail(m)  (m=(unsigned long)farcoreleft())

#include "error.h"

#define OUp   72
#define ODwn  80
#define ORgt  77
#define OLft  75
#define OPUp  73
#define OPDn  81
#define OHom  71
#define OEnd  79
#define OIns  82
#define ODel  83
#define OTab   9
#define OEnt  13
#define OEsc  27
#define OIns  82
#define CUp  141
#define CDwn 145
#define CRgt 116
#define CLft 115
#define CHom 119
#define CEnd 117
#define CPUp 132
#define CPDn 118
#define CIns 146
#define CDel 147
#define Ctl(a) (a+1-'A')

char *just_fid(char *path);
#define EXENAME just_fid(argv[0])

/* Check whether file "name" exist */
int fexist (const char *fid);
/* Make *.bak version of a filename name */
void fputext(const char *oldfid, char *newfid, const char *newext);
/* rename file oldname to newname with making *.bak of newname */
int frename(const char *oldfid, const char *newfid);
/* Some useful macros */
#define MakeLE(a,b) {if((a)>(b)) (a)=(b);}
#define MakeGE(a,b) {if((a)<(b)) (a)=(b);}

#define STRCMP strcmp
#define STRSWITCH(s) {char *switchstr=(s); if(switchstr!=(s)) {
#define STRCASE(value) } else if (0==STRCMP(switchstr,(value))) {
#define STRDEFAULT } else {
#define STRENDSW }}

