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

/* ------------------------------------------------------------------------- *
 * ezscroll.h -- header file for EZ-Scroll
 *
 * Copyright (C) 1996 - 1998, 2006, 2007 Dwight Barkley and Matthew Dowle
 * ------------------------------------------------------------------------- */

/* Modifications for ezview: 2010-18, Vadim Biktashev v.n.biktashev@exeter.ac.uk */

#ifndef _EZVIEW_
#define _EZVIEW_

#include "extern.h"

#define GRAPHICS   1      /* If 1 then run with interactive
			     graphics. Otherwise only write to
			     files.  */
#define DOUBLEBUFFER 0	/* 0/1: On Mac, double-buffering does not
			   offer noticeable speedup, but makes image
			   saving out-of-sync if done in the
			   SIMULATING mode - VNB */

/* -------------------------------------------------------------------------- */
/* Special extra                                                              */
#define MARKER     1      /* if 1 then display a marker with specified coords */
#define FIBRES     1      /* if 1 then visualize fibre directions by lines    */
/* -------------------------------------------------------------------------- */

#ifndef M_PI
#define M_PI	      3.14159265358979323846
#endif

#define message(level,...) {if (S->verbose>=level) MESSAGE(__VA_ARGS__);}


#define BUFLEN 4096 	/* size of all sorts of string buffers */
#define MAXWINDOWTITLE 512 /* window titles longer than that are not feasible */
#define FILAMENT_BUFFER_SIZE (1024L*1024L) /* 1M should be enough even for large but resonable filaments? */

#define SUCCESS 1
#define FAILURE 0

#define OFF       0                   /* Draw_modes */ 
#define LINES     1
#define TRIANGLES 2
#define CURVES    3
#define FAN  	  4
#define STRIP	  5


/* --------------------------------------------------- 
 * All global variables from EZSCROLL are now packed into
 * the STR structure, to allow multiple instances of the
 * ezview device. 
 * --------------------------------------------------- */
typedef GLfloat       GLReal;
typedef enum {
  MODE_SIMULATING,
  MODE_VIEWING,
  MODE_ROTATING
} ezmode_type;
typedef unsigned char CubeEdge;       /* Number between 0 and 12. 0 means "no
				       * edge". Large tables of CubeEdge are
				       * defined so we want it as small as
				       * possible. */
#define MAX_TRIANGLE_LIST_LENGTH  16  /* Maximum triangles in marching cube is
				       * 5 (each with 3 edges) plus 1 for the
				       * end marker. */
#define MAX_EDGE_LIST_LENGTH      13  /* As above but without
				       * repetitions. Indexes 90 and 165 (4
				       * non touching triangles) have an
				       * intersect on each of the 12
				       * edges. Plus one for the end
				       * marker. */

typedef struct {
  #define _(Type,Name,Init,reRead,reView,reLight,reMake,reDraw) PAR(Type,Name)
  #include "ezpar.h" /* used-def parameters, changeable in dialogue */
  #undef _
  #define _(Type,Name,Dflt,miN,maX) PAR(Type,Name);
  #include "ezfix.h" /* used-def parameters, not changeable in dialogue */
  #undef _
  #define _(Type,Name,Init) Type Name;
  #include "ezini.h" /* quasiglobal variables, to be initialized */
  #undef _
  #define _(Type,Name) Type Name;
  #include "ezvar.h" /* quasiglobal variables, not to be initalized */
  #undef _
  #include "ezarr.h" /* arrays and other awkward quasiglobal objects */
} STR;
  
/* ------------------------------------------- 
 * Array de-facto constants
 * These arrays set the light sources and properties. 
 * Used in Draw_ini and in drawing routines in ezmarching.c
 * ------------------------------------------- */
static const GLfloat mat_emission[]  = {0.3, 0.3, 0.3, 1};// {0.5, 0.5, 0.5, 1.0};
static const GLfloat mat_diffusive[] = {1, 1, 1, 0};
static const GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
static const GLfloat mat_shininess[] = {100.0};
static const GLfloat mat_noemission[] = {0,0,0,1};

/* ------------------------------------------------------------------------- 
 * Fields (f,i,j,k)  --  array of fields: f=ulayer for u etc
 * ------------------------------------------------------------------------- */
#define Fields(f,i,j,k)   ((f<vmax)?(New[ind(i,j,k,f)]):(Geom[geom_ind(i,j,k,(f)-vmax)]))
//#define Geomfields(f,i,j,k)  Geom[geom_ind(i,j,k,(f)-vmax)]
/* Protected version, to avoid crash when stepped out of bounds 
   - the only place where taboo is used; to be eliminated */
#define FIELDS(f,i,j,k) (\
  (f<0)?taboo:(f>=nv)?taboo:\
  (i<1)?taboo:(i>nx)?taboo:\
  (j<1)?taboo:(j>ny)?taboo:\
  (k<1)?taboo:(k>nz)?taboo:\
  Fields(f,i,j,k))


/* ------------------------------------------- 
 * Prototypes for public functions defined in: 
 * ------------------------------------------- */

/* ezgraph3d.c 
 * ----------- */
void Make_lists   (STR *S);
void Draw         (STR *S);
int  Draw_ini     (STR *S);
int  Event_check  (STR *S);
void QuitX        (STR *S);
void Save_image   (STR *S);
void Save_filament(STR *S);
void setView      (STR *S, float new_theta, float new_phi,float new_psi, float new_distance);

/* ezmarching.c 
 * ------------ */
void Marching_ini	(STR *S);
void Marching_cubes	(STR *S,unsigned int resolution, int SRF_LIST, int FLM_LIST);
void Draw_bounding_box	(STR *S);
void Draw_marker	(STR *S);
void Make_fibres 	(STR *S);


#if DOUBLEBUFFER
#define gdebug(level) {}
#else
#define gdebug(level) if (S->verbose>=level) glFinish()
#endif

#endif /*  _EZVIEW_  */
