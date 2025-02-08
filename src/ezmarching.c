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

/* ------------------------------------------------------------------------- */
/* ezmarching.c -- Efficient implementation of the Marching Cubes Algorithm
 * first published by Bill Lorensen (1987).
 *
 * Copyright (C) 1996 - 1998, 2006, 2007 Dwight Barkley and Matthew Dowle
 *
 * RCS Information
 * ---------------------------
 * $Revision: 1.5.1.1 $
 * $Date: 2007/05/07 10:07:03 $
 * ------------------------------------------------------------------------- */

/* Modifications for beatbox: 2010-18, Vadim Biktashev v.n.biktashev@exeter.ac.uk */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "beatbox.h"
#include "state.h"
#include "system.h"
#include "device.h"
#include "bikt.h"

#if defined(NOX11) || defined (NOGL)
#else
#include <GL/gl.h>
#include <GL/glx.h>

#include "ezview.h"    
#include "ezgraph3d.h"
#include "ezmarching.h"


/* -------------------------------------------------------------------------
 *
 * This file contains all graphics functions that compute the iso-surfaces
 * and filaments, and render these.
 *
 * The important things to know about this file are:
 * 
 * (1) All functions containing OpenGL calls start with Draw_.  These
 * functions and Marching_cubes() are probably the only ones you would
 * need to change to affect the look of the graphics (contour colors etc.)
 *
 * (2) The U and V contours values are defined in ezmarching.h.  The filament
 * is defined as the intersection of these iso-surfaces.
 *
 * (3) There are three states the draw_mode can be in:
 *   OFF:       glEnd() has been called and nothing being passed to OpenGL.
 *   LINES:     glBegin(GL_LINES) has been called and points on line 
 *              segments are being passed.
 *   TRIANGLES: glBegin(TRIANGLES) has been called and triangle vertices 
 *              are being passed.
 * (VNB: also: CURVES, FAN and STRIP)
 *
 * (4) The last approx 1/4 of the code is devoted to initialization,
 * i.e. setting up the lookup table.  You shouldn't touch this.
 *
 * May 2007: added computation of surface normals for nicer 3D visualization.
 *
 * ------------------------------------------------------------------------- */

/* 
 * Global constants for this file only 
 * ----------------------------------- */
static const Axis        axis[13]  = {0,Z,Y,Z,Y,Z,Y,Z,Y,X,X,X,X}; 
                            /* Given a cube edge, axis tells you which
			     * coordinate axis the edge is on. */

/* 
 * Private functions 
 * ----------------- */

/* Functions for computing iso-surfaces and filaments */
static void       Compute_Surface_and_Filament (STR *S, int SRF_LIST, int FLM_LIST);
static void       Filament_in_cube        (STR *S, int lu, real uconst, int lv, real vconst);

static void       Surface_in_cube 	  (STR *S,int do_surface,int do_filament,
					   int lu,real uconst,int lv,real vconst);
static Bool       Find_two_pts            (GLReal point[2][3], 
					   GLReal point_r[4], 
					   GLReal p[3], GLReal dir[3]);
static void       Set_relative_edge_pos   (STR *S, int inc);
static GLReal     Contour_normal_intersect(STR *S, int layer, real contour, CubeEdge edge, 
					   GLReal *normal);
static int	  Set_surface_color	  (STR *S);

static real	  value			  (STR *S, int l, int i, int j, int k, real dflt);

/* OpenGL-specific drawing functions */
static void	  Draw_triangles	  (STR *S);
static void	  Draw_sphere		  (STR *S, GLReal centre[3], real radius, GLReal emission);
static void       Set_draw_mode_off       (STR *S);
static void       Set_draw_mode_lines     (STR *S, GLReal lwidth);
static void       Set_draw_mode_curves    (STR *S, GLReal lwidth, GLReal emission);
static void       Set_draw_mode_triangles (STR *S);
static void       Set_draw_mode_fan       (STR *S, GLReal emission);
static void       Set_draw_mode_strip     (STR *S, GLReal emission);
static void	  Set_filament_emission	  (STR *S, GLReal emission);

/* Functions for linking filament points into continuous chains */
static void	  Register_segment	  (GLReal point[2][3]);
static void	  Register_chains	  (STR *S);

/* Functions for output of filaments */
static void	  Filament_start	  (STR *S); 
static void 	  Filament_next	  	  (STR *S, GLReal point[3]); 
static void 	  Filament_break	  (STR *S); 
static void 	  Filament_frame  	  (STR *S); 
static void 	  Filament_print	  (STR *S, char *fmt, ...);
static void 	  Filament_rewind	  (STR *S);

/* Functions used in initialization */
static CubeIndex  Rotate_edge_list        (STR *S, CubeIndex index, int direction);
static CubeIndex  Rotate_edge_list_x      (STR *S, CubeIndex index);
static CubeIndex  Rotate_edge_list_y      (STR *S, CubeIndex index);
static void       Rotate_diagonal         (STR *S, CubeIndex index);
static void       Do_all_orientations     (STR *S, CubeIndex index);
static void       Generate_edge_table     (STR *S);
static void       Print_list              (CubeEdge *p);


static real value(STR *S, int l, int i, int j, int k,real dflt)
{
  DEVICE_CONST(int,nx);
  DEVICE_CONST(int,ny);
  DEVICE_CONST(int,nz);
  DEVICE_CONST(int,nv);
  int nneib, rad, i1, j1, k1;
  real vsum;
  static int times_reported=0; /* TODO: should this be privatized, too */
  real v;
  if (l<0)
    return dflt;
  else if (l>=nv)
    return dflt;
  else if (l>=vmax)
    return Fields(l,i,j,k);
  else if (isTissue(i,j,k)) {
    return Fields(l,i,j,k);
  } else {
    for(rad=1;rad<nx;rad++) {
      nneib=0;
      vsum=0.0;
      for (i1=i-rad;i1<=i+rad;i1++) {
	for (j1=j-rad;j1<=j+rad;j1++) {
	  for (k1=k-rad;k1<=k+rad;k1++) {
	    if (!isTissue(i1,j1,k1)) continue;
	    vsum+=Fields(l,i1,j1,k1);
	    nneib++;
	  }
	}
      }
      if (nneib) {
	v=vsum/nneib;
	break;
      }
    }
    if (nneib==0) {
      if (times_reported<=10) {
	MESSAGE("could not find good %d-value around (%d,%d,%d)\n",l,i,j,k);
	if (times_reported==10)
	  MESSAGE("(further similar reports suppressed)\n");
	times_reported++;
      }
      v=dflt;
    }
  }
  return v;
} /* value */
/* ========================================================================= */

void Marching_cubes (STR *S, unsigned int resolution, int SRF_LIST, int FLM_LIST)
{
  DEVICE_CONST(int,imin);
  DEVICE_CONST(int,imax);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_VAR(GLReal,marching_h);
  DEVICE_VAR(GLReal,convert_to_phy);
  DEVICE_VAR(unsigned int,inc);

  (*inc) = resolution;
  (*marching_h) = (*inc)*plot_length[0]/(imax-imin+1);
  (*convert_to_phy) = (*inc)*1.0/(*marching_h);
  Set_relative_edge_pos (S, *inc);
  Compute_Surface_and_Filament (S, SRF_LIST, FLM_LIST);
} /* Marching_cubes */
/* ========================================================================= */

static void Draw_triangles (STR *S)
{
  DEVICE_CONST(GLReal,marching_h);

  int n;

  Set_draw_mode_triangles(S);

  /* DB: I only color at the marching cube level, i.e. all triangles
    within a given marching cube get the same color.  One could color
    at the individual triangle level (inside the while loop). */
  if (0==Set_surface_color(S)) return;
  
  /* DB: Send OpenGL the positions of each triangle vertex.  Every
    group of 3 vertices is interpreted as defining one triangle. The
    triangle table is organized so that all we need do now is pass all
    the vertices in the list to OpenGL. */
  n = 0;
  while (S->triangle_list[U_FIELD][n] != END) {
    GLNORMAL3V(S->normals[S->triangle_list[U_FIELD][n]]);
    GLVERTEX3V(S->vertex[U_FIELD][S->triangle_list[U_FIELD][n]]);
    n++;
  }

  Set_draw_mode_off(S);
} /* Draw_triangles */
/* ========================================================================= */

/* Set the iso-surface color. */
/* VNB: there are several different algorithms which may employ up to
   4 layers of data */
static int Set_surface_color (STR *S)
{
  DEVICE_CONST(unsigned int,ii);
  DEVICE_CONST(unsigned int,jj);
  DEVICE_CONST(unsigned int,kk);
  DEVICE_CONST(int,color_mode);
  DEVICE_CONST(real,alphamin);
  DEVICE_CONST(real,alphamax);
  DEVICE_CONST(int,vlayer);
  DEVICE_CONST(real,vmin);
  DEVICE_CONST(real,maxv);
  DEVICE_CONST(real,vc);
  DEVICE_CONST(int,remove_backs);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_ARRAY(GLReal,pos);
  GLReal red, green, blue, alpha;
  real vvalue, rvalue, gvalue, bvalue, avalue;

  vvalue=value(S,vlayer,ii,jj,kk,vc);

  if (remove_backs) {
    /* Do not draw triangles if v(ii,jj,kk) > V_CONT */
    if (vvalue*remove_backs > vc*remove_backs) return 0;
  }

  switch(color_mode) {
  case 0:  /* Color (yellow-green) according to v-field */
    {
    DEVICE_CONST(int,vlayer);
    DEVICE_CONST(real,vmin);
    DEVICE_CONST(real,maxv);
    DEVICE_CONST(real,vc);
    vvalue=value(S,vlayer,ii,jj,kk,vc);
    red   = (maxv-vvalue)/(maxv-vmin);
    green = 0.5;
    blue  = 0.;
    alpha = alphamax;
    }
    break;

  case 1: /* Color (red-blue) according to other field */
    {
    DEVICE_CONST(int,vlayer);
    DEVICE_CONST(real,vmin);
    DEVICE_CONST(real,maxv);
    DEVICE_CONST(real,vc);
    vvalue=value(S,vlayer,ii,jj,kk,vc);
    red   = (maxv-vvalue)/(maxv-vmin);
    blue  = (vvalue-vmin)/(maxv-vmin);
    green = red*blue*4;
    alpha = alphamin+(alphamax-alphamin)*red;
    }
    break;

  case 2: /* Fixed color (orange), transparency according to other field */
    {
    DEVICE_CONST(int,vlayer);
    DEVICE_CONST(real,vmin);
    DEVICE_CONST(real,maxv);
    DEVICE_CONST(real,vc);
    vvalue=value(S,vlayer,ii,jj,kk,vc);
    red   = 1.;
    green = 0.5;
    blue  = 0;
    alpha = alphamin+(alphamax-alphamin)*(maxv-vvalue)/(maxv-vmin);
    }
    break;

  case 3: /* Similar to 1, with total brightness constant */
    {
    DEVICE_CONST(int,vlayer);
    DEVICE_CONST(real,vmin);
    DEVICE_CONST(real,maxv);
    DEVICE_CONST(real,vc);
    vvalue=value(S,vlayer,ii,jj,kk,vc);
    red   = (maxv-vvalue)/(maxv-vmin);
    blue  = (vvalue-vmin)/(maxv-vmin);
    red  = (red <0)?0:(red >1)?1:(red * red);
    blue = (blue<0)?0:(blue>1)?1:(blue*blue);
    green = 1.0 - red - blue;
    alpha = alphamin+(alphamax-alphamin)*red;
    }
    break;

  /* Backward compatibility: algs 4 and 5 have been 'merged' */
  case 4: case 5: /* Color components specified by FOUR given fields */
    {
    DEVICE_CONST(real,vc);
    DEVICE_CONST(int,rlayer);
    DEVICE_CONST(real,rmin);
    DEVICE_CONST(real,rmax);
    DEVICE_CONST(int,glayer);
    DEVICE_CONST(real,gmin);
    DEVICE_CONST(real,gmax);
    DEVICE_CONST(int,blayer);
    DEVICE_CONST(real,bmin);
    DEVICE_CONST(real,bmax);
    DEVICE_CONST(int,alayer);
    DEVICE_CONST(real,amin);
    DEVICE_CONST(real,amax);    
    rvalue = (rlayer>=0)?value(S,rlayer,ii,jj,kk,vc):0;
    gvalue = (glayer>=0)?value(S,glayer,ii,jj,kk,vc):0;
    bvalue = (blayer>=0)?value(S,blayer,ii,jj,kk,vc):0;
    avalue = (alayer>=0)?value(S,alayer,ii,jj,kk,vc):0;
    red   = (rvalue-rmin)/(rmax-rmin);
    green = (gvalue-gmin)/(gmax-gmin);
    blue  = (bvalue-bmin)/(bmax-bmin);
    alpha = alphamin+(alphamax-alphamin)*(avalue-amin)/(amax-amin);
    }
    break;

  default: /* Color according to position in space, Useful if not computing S->normals */
    red   = pos[Y]/plot_length[1]; 
    green = pos[Z]/plot_length[2]; 
    blue  = pos[X]/plot_length[0];
    alpha = alphamax;
    break;

  }

  GLCOLOR4(red, green, blue, alpha);
  return 1;
} /* Set_surface_color */

/*********************************************************/
/* Collecting filament segments into continuous "curves" */
/* TODO: the global objects here ought to be privatized! */

/* The Marching Cube algorithm generates points in duos. */
typedef struct listnode {
  GLReal piece[2][3];
  struct listnode *prev;
  struct listnode *next;
} listitemtype;

/* This is transformed into chains consisting of solo points */
typedef struct chainnode {
  GLReal point[3];
  struct chainnode *prev;
  struct chainnode *next;
} chainitemtype;

/* The ends of the linked list of points build by Register_segment() */
/* and consumed by Register_chains() */
static listitemtype *listfirst;
static listitemtype *listlast;

static int numlistitems=0;
static listitemtype *newlistitem (void)
{
  listitemtype *p=calloc(1,sizeof(listitemtype));
  if (!p) {
    fprintf(stderr,"could not alloc memory for list item; abort\n");
    exit(1);
  }
  numlistitems++;
  return p;
}

static void freelistitem (listitemtype *p)
{
  free(p);
  numlistitems--;
}

static int numchainitems=0;
static chainitemtype *newchainitem (void)
{
  chainitemtype *p=calloc(1,sizeof(chainitemtype));
  if (!p) {
    fprintf(stderr,"could not alloc memory for chain item; abort\n");
    exit(1);
  }
  numchainitems++;
  return p;
}

static void freechainitem (chainitemtype *p)
{
  free(p);
  numchainitems--;
}

static void pointcopy (chainitemtype *dest, GLReal src[3])
{
  int i;
  if (!dest) return;
  if (!src) return;
  for(i=0;i<3;i++) dest->point[i]=src[i];
}

/* are two given points close enough to be considered the same? */
static int closepoints (chainitemtype *item, GLReal point[3])
{
  GLReal small=1.e-4; /* even this still breaks them up */
  int i;
  for (i=0;i<3;i++) if (fabs(item->point[i]-point[i])>small) return 0;
  return 1;
}

/* Add the given piece to the linked list of segments */
static void Register_segment (GLReal point[2][3])
{
  int j;
  listitemtype *listitem;
  listitemtype *lastbutone;

  listitem=newlistitem();
  for(j=0;j<3;j++){
    listitem->piece[0][j] = point[0][j];
    listitem->piece[1][j] = point[1][j]; 
  }
  if (listfirst && listlast) {
    lastbutone=listlast;
    listlast=listitem;
    listlast->prev=lastbutone;
    listlast->next=NULL;
    lastbutone->next=listlast;
  } else if (listfirst || listlast) {
    fprintf(stderr,"this cannot be in %s:%d\n",__FILE__,__LINE__);
    exit(1);
  } else {
    listitem->prev=listitem->next=NULL;
    listfirst=listlast=listitem;
  }
} /* Register_segment */
/* ========================================================================= */

/* Exhausts the current list of points and builds */
/* continuos chains out of those points */
static void Register_chains (STR *S)
{
  DEVICE_CONST(int,flm_r);
  DEVICE_CONST(int,flm_g);
  DEVICE_CONST(int,flm_b);
  DEVICE_CONST(int,flm_wt);
  DEVICE_CONST(int,flm_balls);
  DEVICE_CONST(int,flm_em);
  DEVICE_CONST(GLReal,convert_to_phy);
  listitemtype *listnext;
  listitemtype *listitem;
  chainitemtype *chainfirst;
  chainitemtype *chainlast;
  chainitemtype *addbefore;
  chainitemtype *addafter;
  chainitemtype *chainitem;
  chainitemtype *chainnext;
  int trymore, matched;

  Filament_rewind(S);
  
  /* The input list of segments is still not empty, so ... */
  while (NULL!=listfirst) {

    /* start a new continuous chain from the first available segment */
    chainfirst=newchainitem();
    chainlast=newchainitem();

    pointcopy(chainfirst,listfirst->piece[0]);
    chainfirst->next=chainlast;
    chainfirst->prev=NULL;

    pointcopy(chainlast,listfirst->piece[1]);
    chainlast->next=NULL;
    chainlast->prev=chainfirst;

    /* unlink that segment from the original list */
    listnext=listfirst->next;
    freelistitem(listfirst);
    if (listnext) {
      listnext->prev=NULL;
      listfirst=listnext;
    } else {
      listfirst=listlast=NULL;
    }

    /* Try to grow the chain by trying out all available segments one by one */
    do {
      /* Loop through the rest of the original list to find next matching segment */
      listitem=listfirst;
      trymore=0;
      while(listitem!=NULL) {
	matched=1;

	/* there are four ways a given segment may match the existing chain */
	if (closepoints(chainfirst,listitem->piece[0])) {
	  addbefore=newchainitem();
	  pointcopy(addbefore,listitem->piece[1]);
	  addbefore->next=chainfirst;
	  addbefore->prev=NULL;
	  chainfirst->prev=addbefore;
	  chainfirst=addbefore;
	} else if (closepoints(chainfirst,listitem->piece[1])) {
	  addbefore=newchainitem();
	  pointcopy(addbefore,listitem->piece[0]);
	  addbefore->next=chainfirst;
	  addbefore->prev=NULL;
	  chainfirst->prev=addbefore;
	  chainfirst=addbefore;
	} else if (closepoints(chainlast,listitem->piece[0])) {
	  addafter=newchainitem();
	  pointcopy(addafter,listitem->piece[1]);
	  addafter->prev=chainlast;
	  addafter->next=NULL;
	  chainlast->next=addafter;
	  chainlast=addafter;
	} else if (closepoints(chainlast,listitem->piece[1])) {
	  addafter=newchainitem();
	  pointcopy(addafter,listitem->piece[0]);
	  addafter->prev=chainlast;
	  addafter->next=NULL;
	  chainlast->next=addafter;
	  chainlast=addafter;
	} else {
	  matched=0;
	}
	listnext=listitem->next;
	if (matched) {
	  /* remove the matching segment from the original list */
	  /* taken care in case it is the first or the last one there */
	  if (listitem->prev!=NULL) {
	    listitem->prev->next=listitem->next;
	  } else {
	    listfirst=listitem->next;
	  }
	  if (listitem->next!=NULL) {
	    listitem->next->prev=listitem->prev;
	  } else {
	    listlast=listitem->prev;
	  }
	  freelistitem(listitem);
	  trymore=1;
	}
	listitem=listnext;
      }; /* while listitem != NULL */
    } while(trymore);

    /* Now we can output the continuous chain. */
    /* This is either by lines or by balls. */
    if (chainfirst) {
      Filament_start (S); /* start of a new continuous chain */
      chainitem=chainfirst;
      do { /* output all chain points */
	Filament_next (S, chainitem->point);
	chainnext=chainitem->next;
	if (chainnext) chainnext->prev=NULL;
	freechainitem(chainitem);
	chainitem=chainnext;
      } while (chainitem!=NULL);
      Filament_break (S); /* end of a continuous chain */
    } /* if chainfirst */
  } /* while listfirst */

  if (numlistitems!=0) 
    message(1,"Warning: unused list items in the end of Register_chains, %s:%d\n",__FILE__,__LINE__);
  Filament_frame(S); /* to distinguish temporal frames in the output file */
} /* Register_chains */
/* ========================================================================= */


static void Filament_start (STR *S)
{
  DEVICE_CONST(real,flm_r);
  DEVICE_CONST(real,flm_g);
  DEVICE_CONST(real,flm_b);
  DEVICE_CONST(real,flm_wt);
  DEVICE_CONST(real,flm_em);
  DEVICE_CONST(int,flm_balls);
  GLCOLOR3(flm_r, flm_g, flm_b);
  if (flm_balls==0) {
    Set_draw_mode_curves(S,flm_wt,flm_em);
  }
} /* Filament_start */
/* ========================================================================= */

static void Filament_next (STR *S, GLReal point[3])
{
  DEVICE_CONST(real,convert_to_phy);
  DEVICE_CONST(int,flm_balls);
  DEVICE_CONST(real,flm_wt);
  DEVICE_CONST(real,flm_em);
  float x=convert_to_phy*point[X];
  float y=convert_to_phy*point[Y];
  float z=convert_to_phy*point[Z];
  Filament_print(S,"%ld %.5f %.5f %.5f\n",(long)t, x, y, z);
  message(4,"%ld: (%.5f %.5f %.5f)\n",(long)t, x, y, z);
	
  if (flm_balls)
    Draw_sphere(S,point,flm_wt,flm_em);
  else
    GLVERTEX3V(point);
} /* Filament_next */
/* ========================================================================= */

/* All points of this chain have been plotted/written */
static void Filament_break (STR *S)
{
  DEVICE_CONST(int,flm_balls);
  Filament_print(S,"# end of chain\n\n");
  if (!flm_balls) Set_draw_mode_off(S);
} /* Filament_break */
/* ========================================================================= */


static void Filament_frame (STR *S)
{
  Filament_print(S,"# end of frame\n\n");
} /* Filament_frame */
/* ========================================================================= */

/* the box with axes labels */
void Draw_bounding_box (STR *S)
{
  DEVICE_CONST(real,bbox_r);
  DEVICE_CONST(real,bbox_g);
  DEVICE_CONST(real,bbox_b);
  DEVICE_CONST(real,bbox_a);
  DEVICE_CONST(real,bbox_wt);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_ARRAY(real,plot_ini);
  real x0=plot_ini[0]; 
  real xm=plot_ini[0]+plot_length[0]*0.5;
  real x1=plot_ini[0]+plot_length[0];
  real y0=plot_ini[1]; 
  real ym=plot_ini[1]+plot_length[1]*0.5;
  real y1=plot_ini[1]+plot_length[1];
  real z0=plot_ini[2]; 
  real zm=plot_ini[2]+plot_length[2]*0.5;
  real z1=plot_ini[2]+plot_length[2];

  real cw=0.02;         /* char width compared to plot size */
  real ch=sqrt(2.0)*cw; /* char height compared to plot size */
  real hh=0.5*ch;       /* half-height */
  real sq=sqrt(0.5);    /* to make 45deg slopes */
  real d=0.5*ch*sq;     /* dist from axis */

  /* the box */
  Set_draw_mode_lines(S,bbox_wt);

  GLCOLOR4(bbox_r, bbox_g, bbox_b, bbox_a);

  GLVERTEX3(x0,y0,z0);  
  GLVERTEX3(x1,y0,z0);
  GLVERTEX3(x1,y0,z0);  
  GLVERTEX3(x1,y1,z0);
  GLVERTEX3(x1,y1,z0);
  GLVERTEX3(x0,y1,z0);
  GLVERTEX3(x0,y1,z0);
  GLVERTEX3(x0,y0,z0);

  GLVERTEX3(x0,y0,z0);
  GLVERTEX3(x0,y0,z1);
  GLVERTEX3(x1,y0,z0);
  GLVERTEX3(x1,y0,z1);
  GLVERTEX3(x1,y1,z0);
  GLVERTEX3(x1,y1,z1);
  GLVERTEX3(x0,y1,z0);
  GLVERTEX3(x0,y1,z1);

  GLVERTEX3(x0,y0,z1);  
  GLVERTEX3(x1,y0,z1);
  GLVERTEX3(x1,y0,z1);  
  GLVERTEX3(x1,y1,z1);
  GLVERTEX3(x1,y1,z1);  
  GLVERTEX3(x0,y1,z1);
  GLVERTEX3(x0,y1,z1);  
  GLVERTEX3(x0,y0,z1);

  /* the axes labels */
  Set_draw_mode_lines(S,1.5*bbox_wt);
  /* X */
  GLVERTEX3(xm-0.5*cw ,y0      -d,z0      -d);
  GLVERTEX3(xm+0.5*cw ,y0-d-sq*ch,z0-d-sq*ch);
  GLVERTEX3(xm+0.5*cw ,y0      -d,z0      -d);
  GLVERTEX3(xm-0.5*cw ,y0-d-sq*ch,z0-d-sq*ch);
  /* Y */
  GLVERTEX3(x0-d-sq*hh,ym        ,z0-d-sq*hh);
  GLVERTEX3(x0-d      ,ym+0.5*cw ,z0-d      );
  GLVERTEX3(x0-d-sq*ch,ym+0.5*cw ,z0-d-sq*ch);
  GLVERTEX3(x0-d      ,ym-0.5*cw ,z0-d      );
  /* Z */
  GLVERTEX3(x0-d      ,y0-d      ,zm-0.5*cw );
  GLVERTEX3(x0-d      ,y0-d      ,zm+0.5*cw );
  GLVERTEX3(x0-d      ,y0-d      ,zm+0.5*cw );
  GLVERTEX3(x0-d-sq*ch,y0-d-sq*ch,zm-0.5*cw );
  GLVERTEX3(x0-d-sq*ch,y0-d-sq*ch,zm-0.5*cw );
  GLVERTEX3(x0-d-sq*ch,y0-d-sq*ch,zm+0.5*cw );

  Set_draw_mode_off(S);
} /* Draw_bounding_box */
/* ========================================================================= */

#if MARKER
/* coords start with 1 and end with NQ, and 1 is the min index for output in beatbox */
#define gridX(x) ((x-imin)*plot_length[0]/(imax-imin+1))
#define gridY(y) ((y-jmin)*plot_length[1]/(jmax-jmin+1))
#define gridZ(z) ((z-kmin)*plot_length[2]/(kmax-kmin+1))
#define gridvertex(x,y,z) GLVERTEX3(gridX(x),gridY(y),gridZ(z))
void Draw_marker (STR *S)
{
  DEVICE_CONST(int,nx);
  DEVICE_CONST(int,ny);
  DEVICE_CONST(int,nz);
  DEVICE_CONST(int,nv);
  DEVICE_CONST(real,taboo);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_CONST(int,imin);
  DEVICE_CONST(int,imax);
  DEVICE_CONST(int,jmin);
  DEVICE_CONST(int,jmax);
  DEVICE_CONST(int,kmin);
  DEVICE_CONST(int,kmax);
  DEVICE_CONST(real,marker_size);
  DEVICE_CONST(real,marker_x);
  DEVICE_CONST(real,marker_y);
  DEVICE_CONST(real,marker_z);
  DEVICE_CONST(real,marker_r);
  DEVICE_CONST(real,marker_g);
  DEVICE_CONST(real,marker_b);
  DEVICE_CONST(real,marker_a);
  DEVICE_CONST(real,marker_wt);
  DEVICE_CONST(int,verbose);
  
  if (!marker_size) return;

  Set_draw_mode_lines(S,marker_wt);

  GLCOLOR4(marker_r, marker_g, marker_b, marker_a);

  /* The main body of the probe */
  gridvertex(marker_x+marker_size,marker_y,marker_z);
  gridvertex(marker_x-marker_size,marker_y,marker_z);

  gridvertex(marker_x,marker_y+marker_size,marker_z);
  gridvertex(marker_x,marker_y-marker_size,marker_z);
 
  gridvertex(marker_x,marker_y,marker_z+marker_size);
  gridvertex(marker_x,marker_y,marker_z-marker_size);

  /* Its shadows on the box faces */
  GLCOLOR4(marker_r, marker_g, marker_b, 0.5*marker_a);

  gridvertex(imin,marker_y+marker_size,marker_z);
  gridvertex(imin,marker_y-marker_size,marker_z);
  gridvertex(imin,marker_y,marker_z+marker_size);
  gridvertex(imin,marker_y,marker_z-marker_size);
  gridvertex(imax,marker_y+marker_size,marker_z);
  gridvertex(imax,marker_y-marker_size,marker_z);
  gridvertex(imax,marker_y,marker_z+marker_size);
  gridvertex(imax,marker_y,marker_z-marker_size);

  gridvertex(marker_x+marker_size,jmin,marker_z);
  gridvertex(marker_x-marker_size,jmin,marker_z);
  gridvertex(marker_x,jmin,marker_z+marker_size);
  gridvertex(marker_x,jmin,marker_z-marker_size);
  gridvertex(marker_x+marker_size,jmax,marker_z);
  gridvertex(marker_x-marker_size,jmax,marker_z);
  gridvertex(marker_x,jmax,marker_z+marker_size);
  gridvertex(marker_x,jmax,marker_z-marker_size);

  gridvertex(marker_x+marker_size,marker_y,kmin);
  gridvertex(marker_x-marker_size,marker_y,kmin);
  gridvertex(marker_x,marker_y+marker_size,kmin);
  gridvertex(marker_x,marker_y-marker_size,kmin);
  gridvertex(marker_x+marker_size,marker_y,kmax);
  gridvertex(marker_x-marker_size,marker_y,kmax);
  gridvertex(marker_x,marker_y+marker_size,kmax);
  gridvertex(marker_x,marker_y-marker_size,kmax);
  
  Set_draw_mode_off(S);

  if (verbose>=3) {
    int l;
    int x=floor(marker_x+0.5);
    int y=floor(marker_y+0.5);
    int z=floor(marker_z+0.5);
    MESSAGE("(x,y,z)=(%d,%d,%d) u=(",x,y,z);
    for (l=0;l<nv;l++)
      /* this is the only place where macro FIELDS is still used? */
      MESSAGE(FFMT"%s",FIELDS(l,x,y,z),((l==nv-1)?")\n":","));
  }
} /* Draw_marker */
#undef gridvertex
#undef gridZ
#undef gridY
#undef gridX
/* ========================================================================= */
#endif

/* Draw a "sphere" with center at given plot coords, with filament colour */
#define plotvertex(x,y,z) GLVERTEX3(c[X]+x,c[Y]+y,c[Z]+z)
/* "physical" spherical coords: inclination theta, azimuth phi */
#define spvertex(phi,theta) \
  GLNORMAL3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));	\
  plotvertex(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))
#define spnormal(phi,theta) GLNORMAL3(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))
void Draw_sphere (STR *S,GLReal c[3], real radius, GLReal emission)
{
  DEVICE_CONST(int,nx);
  DEVICE_CONST(int,ny);
  DEVICE_CONST(int,nz);
  DEVICE_CONST(int,nv);
  DEVICE_CONST(int,flm_r);
  DEVICE_CONST(int,flm_g);
  DEVICE_CONST(int,flm_b);
  DEVICE_CONST(int,flm_wt);
  DEVICE_CONST(int,flm_balls);
  DEVICE_CONST(int,flm_em);
  int nt=5;
  int np=10;
  int t, p;
  double theta, theta0, theta1, phi;
  double dtheta=M_PI/nt;
  double dphi=2.0*M_PI/np;
  int nmax = max(max(nx,ny),nz);
  double r=radius/(nmax-1);

  if (r<=0) return;

  GLCOLOR3(flm_r, flm_g, flm_b);

  /* Northern pole fan */
  Set_draw_mode_fan(S,emission);
  plotvertex(0,0,r);
  theta=dtheta;
  for (p=0;p<=np;p++) {
    phi=(p==np)?0:p*dphi; /* to avoid gaps due to round-off */
    spvertex(phi,theta);
    spnormal(phi+0.5*dphi,0.5*dtheta);
  }

  /* Southern pole fan */
  Set_draw_mode_fan(S,emission);
  plotvertex(0,0,-r);
  theta=M_PI-dtheta;
  for (p=0;p<=np;p++) {
    phi=(p==np)?0:p*dphi; /* to avoid gaps due to round-off */
    spvertex(phi,theta);
    spnormal(phi+0.5*dphi,M_PI-0.5*dtheta);
  }
  
  /* Intermediate latitude strips */
  for (t=1;t<nt-1;t++) {
    theta0=t*dtheta;
    theta1=(t+1)*dtheta;
    Set_draw_mode_strip(S,emission);
    for (p=0;p<=np;p++) {
      phi=(p==np)?0:p*dphi; /* to avoid gaps due to round-off */
      spvertex(phi,theta0);
      spvertex(phi,theta1);
      spnormal(phi+0.5*dphi,0.5*(theta0+theta1));
    }
  }

  Set_draw_mode_off(S);
} /* Draw_sphere */
#undef spvertex
#undef plotvertex
/* ========================================================================= */

static void Compute_Surface_and_Filament (STR *S, int SRF_LIST, int FLM_LIST)
{
  DEVICE_CONST(int,imin);
  DEVICE_CONST(int,imax);
  DEVICE_CONST(int,jmin);
  DEVICE_CONST(int,jmax);
  DEVICE_CONST(int,kmin);
  DEVICE_CONST(int,kmax);
  DEVICE_CONST(int,ulayer);
  DEVICE_CONST(real,uc);
  DEVICE_CONST(int,layer1);
  DEVICE_CONST(real,const1);
  DEVICE_CONST(int,layer2);
  DEVICE_CONST(real,const2);
  DEVICE_CONST(int,verbose);
  DEVICE_CONST(GLReal,marching_h);
  DEVICE_ARRAY(GLReal,pos);
  DEVICE_VAR(unsigned int,ii);
  DEVICE_VAR(unsigned int,jj);
  DEVICE_VAR(unsigned int,kk);
  DEVICE_CONST(unsigned int,inc);

  /* DEVICE_ARRAY(GLReal[13][3],S->vertex); */
  /* cube counters */
  int i, j, k, idir, jdir, kdir, istart, jstart, kstart, istop, jstop, kstop;
  float small=0.1;
  int do_filament;

  message(3,"Compute_Surface_and_Filament (%d,%d)\n",SRF_LIST,FLM_LIST);

  do_filament = (FLM_LIST!=0);

  if (verbose>=4) {
    if (ulayer==layer1 && uc==const1) {
      MESSAGE("will be calling Surface_in_cube(S,%d,%d,%d,%g,%d,%g)\n",
	     1,do_filament,ulayer,uc,layer2,const2);
    } else {
      MESSAGE("will be calling Surface_in_cube(S,%d,%d,%d,%g,%d,%g)\n",
	       1,0,ulayer,uc,layer2,const2);
      if (do_filament) {
	MESSAGE("will be calling Surface_in_cube(S,%d,%d,%d,%g,%d,%g)\n",
	       0,1,layer1,const1,layer2,const2);
      }
    }
  }
  /* The surface list is made within this cycle, the filament list is deferred */
  if (SRF_LIST) glNewList(SRF_LIST, GL_COMPILE);
  pos[X] = 0.0;
  for ((*ii)=imin+inc; (*ii)<=imax; (*ii)+=inc) {
    pos[X] += marching_h;

    /* Each triangle S->vertex has only one linearly interpolated component
     * (coordinate). The other 2 are vertices lie on a cube edge and are
     * known.  So for this Y-Z plane, we already know the X component of all
     * the edges not parallel to the X axis. The vertices on the edges
     * parallel to the X axis will be linearly interpolated as necessary
     * later. */

    S->vertex[U_FIELD][1][X] = 
    S->vertex[U_FIELD][2][X] = 
    S->vertex[U_FIELD][3][X] = 
    S->vertex[U_FIELD][4][X] = pos[X];
    S->vertex[U_FIELD][5][X] = 
    S->vertex[U_FIELD][6][X] = 
    S->vertex[U_FIELD][7][X] = 
    S->vertex[U_FIELD][8][X] = pos[X] - marching_h;

    pos[Y] = 0.0;
    for ((*jj)=jmin+inc; (*jj)<=jmax; (*jj)+=inc) {
      pos[Y] += marching_h;
    
      /* Similar to above, we can set the components of each edge not on the
	 Y axis. */

      S->vertex[U_FIELD][1][Y]  = 
      S->vertex[U_FIELD][5][Y]  = 
      S->vertex[U_FIELD][9][Y]  = 
      S->vertex[U_FIELD][10][Y] = pos[Y] - marching_h;
      S->vertex[U_FIELD][3][Y]  = 
      S->vertex[U_FIELD][7][Y]  = 
      S->vertex[U_FIELD][11][Y] = 
      S->vertex[U_FIELD][12][Y] = pos[Y];

      pos[Z] = 0.0;
      for ((*kk)=kmin+inc; (*kk)<=kmax; (*kk)+=inc) {
	pos[Z] += marching_h;

	if (ulayer==layer1 && uc==const1) {
	  message(5,"calling Surface_in_cube(S,%d,%d,%d,%g,%d,%g)\n",
		  1,do_filament,ulayer,uc,layer2,const2);
	  Surface_in_cube(S,1,do_filament,ulayer,uc,layer2,const2);
	  message(5,"Surface_in_cube returned\n");
	} else {
	  message(5,"calling Surface_in_cube(S,%d,%d,%d,%g,%d,%g)\n",
		  1,0,ulayer,uc,layer2,const2);
	  Surface_in_cube(S,1,0,ulayer,uc,layer2,const2);
	  message(5,"Surface_in_cube returned\n");
	  if (do_filament) {
	    message(5,"calling Surface_in_cube(S,%d,%d,%d,%g,%d,%g)\n",
		    0,1,layer1,const1,layer2,const2);
	    Surface_in_cube(S,0,1,layer1,const1,layer2,const2);
	    message(5,"Surface_in_cube returned\n");
	  } /* if do_filament */
	} /* if ulayer=layer1 ... else */
      } /* for k */
    } /* for j */
  } /* for i */
  if (SRF_LIST) glEndList();

  /* Now we can compile the filament list */
  /* using filament description stored during the above cycle. */
  /* Remember we decided we shall make the filament list in any case */
  if (FLM_LIST) glNewList(FLM_LIST, GL_COMPILE);
  Register_chains(S);
  if (FLM_LIST) glEndList();
} /* Compute_Surface_and_Filament */
/* ========================================================================= */

/* Calculate the index for the current cube for the given,
 * and draw the piece of surface and/or filament in it, as appropriate. 
 * Recall, the cube is defined by (ii,jj,kk) to (ii-inc,jj-inc,kk-inc)
 * (which are all static globals in this file). 
 * The index considers the cube S->vertex (ii,jj,kk) to be bit 3 of the
 * index.
 */
void Surface_in_cube (STR *S, int do_surface, int do_filament, int lu, real uconst, int lv, real vconst)
{
  if (!do_surface && !do_filament) {MESSAGE("Surface_in_cube without surface or filament"); return;}
  
  DEVICE_CONST(int,nx);
  DEVICE_CONST(int,ny);
  DEVICE_CONST(int,nz);
  DEVICE_CONST(int,nv);
  DEVICE_CONST(real,taboo);
  DEVICE_CONST(GLReal,marching_h);
  DEVICE_ARRAY(GLReal,pos);
  DEVICE_CONST(unsigned int,ii);
  DEVICE_CONST(unsigned int,jj);
  DEVICE_CONST(unsigned int,kk);
  DEVICE_CONST(unsigned int,inc);

  CubeIndex  index, other_index;
  CubeEdge   edge;
  int        n;
  int omit; 
  int di, dj, dk;
  int extra=1; /* consider playing with this? */

  real u100;
  real u101;
  real u111;
  real u110;
  real u000;
  real u001;
  real u011;
  real u010;

  if (lu<0 || lu>=nv) return;
  
  u100=Fields(lu,ii    ,jj-inc,kk-inc);
  u101=Fields(lu,ii    ,jj-inc,kk);
  u111=Fields(lu,ii    ,jj    ,kk);
  u110=Fields(lu,ii    ,jj    ,kk-inc);
  u000=Fields(lu,ii-inc,jj-inc,kk-inc);
  u001=Fields(lu,ii-inc,jj-inc,kk);
  u011=Fields(lu,ii-inc,jj    ,kk);
  u010=Fields(lu,ii-inc,jj    ,kk-inc);
  index =
    ((u100 >= uconst)     ) |
    ((u101 >= uconst) << 1) |
    ((u111 >= uconst) << 2) |
    ((u110 >= uconst) << 3) |
    ((u000 >= uconst) << 4) |
    ((u001 >= uconst) << 5) |	      
    ((u011 >= uconst) << 6) |
    ((u010 >= uconst) << 7);

 /* If the index is non-trivial then the contour passes through this
   * cube. Otherwise we are finished. */
  
  if ((index != 0) && (index != 255)) {
    
    /* Calculate the remaining components of any triangle vertices
     * known without interpolation.  We did the X and Y components
     * earlier. */
    
    S->vertex[U_FIELD][2][Z]  = 
    S->vertex[U_FIELD][6][Z]  = 
    S->vertex[U_FIELD][10][Z] = 
    S->vertex[U_FIELD][12][Z] = pos[Z];
    S->vertex[U_FIELD][4][Z]  = 
    S->vertex[U_FIELD][8][Z]  = 
    S->vertex[U_FIELD][9][Z]  = 
    S->vertex[U_FIELD][11][Z] = pos[Z] - marching_h;
    
    /* Linearly interpolate along the necessary cube edges to finish up
     * the vertices */
    
    n = 0;
    while ((edge = S->edge_table[index][n]) != END) {
      S->vertex[U_FIELD][edge][axis[edge]] =
	Contour_normal_intersect(S,lu,uconst,edge,S->normals[edge]);
      n++;
    }
    S->triangle_list[U_FIELD] = S->triangle_table[index];
    
    /* All the triangles have been found. Draw them if necessary */
    
    if (do_surface) Draw_triangles (S);
    
    if (do_filament) {   /* If filament is needed
			  * compute  it, otherwise 
			  * we are finished */
      
      /* Only draw filaments inside the body */
      omit=0;
      for (di=-inc-extra;di<=extra;di++) {
	for (dj=-inc-extra;dj<=extra;dj++) {
	  for (dk=-inc-extra;dk<=extra;dk++) {
	    if (isTissue(ii+di,jj+dj,kk+dk) == 0.0) {
	      omit=1;
	      break;
	    }
	  }
	  if (omit) break;
	}
	if (omit) break;
      }
      if (!omit) {
	  /* Calculate index for other lu */
	  other_index =
	    ((Fields(lv,ii    ,jj-inc,kk-inc) >= vconst)     ) |
	    ((Fields(lv,ii    ,jj-inc,kk)     >= vconst) << 1) |
	    ((Fields(lv,ii    ,jj    ,kk)     >= vconst) << 2) |
	    ((Fields(lv,ii    ,jj    ,kk-inc) >= vconst) << 3) |
	    ((Fields(lv,ii-inc,jj-inc,kk-inc) >= vconst) << 4) |
	    ((Fields(lv,ii-inc,jj-inc,kk)     >= vconst) << 5) |
	    ((Fields(lv,ii-inc,jj    ,kk)     >= vconst) << 6) |
	    ((Fields(lv,ii-inc,jj    ,kk-inc) >= vconst) << 7);
	  
	  if ((other_index != 0) && (other_index != 255)) {
	    
	    /* Both u and v contours pass through this cube, so there is
	     * potential for contour intersection (i.e. the filament) */
	    
	    S->triangle_list[V_FIELD] = S->triangle_table[other_index];
	    
	    /* Copy over the existing vertices. Each of which only one
	       component may differ (if an intersect exists). */
	    
	    for (n=1; n<=12; n++) {
	      S->vertex[V_FIELD][n][X] = S->vertex[U_FIELD][n][X];
	      S->vertex[V_FIELD][n][Y] = S->vertex[U_FIELD][n][Y];
	      S->vertex[V_FIELD][n][Z] = S->vertex[U_FIELD][n][Z];
	    }
	    
	    /* Do the interpolation */
	    
	    n = 0;
	    while ((edge = S->edge_table[other_index][n]) != END) {
	      S->vertex[V_FIELD][edge][axis[edge]] =
		Contour_normal_intersect(S,lv,vconst,edge,NULL);
	      n++;
	    }
	    
	    /* Now we have all the U and V triangles for this cube so we
	     * can look for the filament. Cube edges that do not contain a
	     * triangle S->vertex will contain junk values. */
	    
	    /* You may want to try organizing things differently when just
	     * filament is required. At the moment, the Y and Z components
	     * of the vertices are being set too early. You only need set
	     * them before the following call. However, they are only
	     * assignments: the X is ok because it would be unusual not to
	     * get a filament through some part of each Y-Z plane; the Z is
	     * only when the U is non-trivial. */
	    
	    Filament_in_cube (S,lu,uconst,lv,vconst);
	    
	  } /* if otherindex */
      } /* if !omit */
    } /* if do_filament */
  } /* if index */ 
} /* Surface_in_cube */
/* ========================================================================= */

static void Filament_in_cube (STR *S, int lu, real uconst, int lv, real vconst)
{
  /* At the current cube (ii,jj,kk) both U and V have a non-trivial index.
   * This function finds any triangle intersections (i.e. filament) in this
   * cube.
   *
   * We look at each pair of triangles and find their common line if it
   * exists. The filament thus consists of many small lines.  Experience has
   * shown that usually there are usually only 2 (or maybe 3) triangles in
   * each cube and for this small number of triangles the algorithm should be
   * efficient. Other methods were considered such as following the filament
   * through the cube but then issues such as what do you when there is more
   * than one filament in one cube complicate such a solution. */
  unsigned int  field;		  /* 0/1 field identifier */
  unsigned int  t_num[2];          /* Number of triangle (0-4) */
  unsigned int  t_pos[2];          /* Position in triangle list. */
  CubeEdge      t_cube_edge[2][3]; /* The cube edges of the current U
					   * and V triangle. */
  GLReal        t_edge_vect[2][5][3][3];
                                          /* All the triangle edge vectors. */
  GLReal        *te0, *te1, *te;   /* Convenience pointers. */
  GLReal        plane[2][5][4];    /* Plane coefficients for each
					   * triangle. */
  GLReal        *u_plane, *v_plane;/* Convenience pointers. */
  GLReal        p_t0[3], p_t1[3];  /* Vector going from p (the reference
					   * point on the common line) to
					   * the triangle S->vertex. */
  GLReal        point[2][3];       /* The line segment to be drawn. */
  GLReal        den,dc,db,ad;      /* Arithmetic storage. */
  unsigned int  tri_edge;          /* Loop variable. */
  GLReal        r,s;               /* Arithmetic storage. */
  GLReal        p[3], dir[3]; 
  GLReal        point_r[4];        /* Used to find the
					   * filament. Defines the intersects
					   * of the common plane line with
					   * each triangle edge. */
  unsigned int  points_found;      /* Intersection points found so far
					   * between the common plane line and
					   * triangle edges. */
  Axis          dim1,dim2;         /* Dimensions to perform line
					   * intersection in. */
  Axis          next_axis[3] = {Y,Z,X};
                                          /* Convenient array to cycle axes. */
  unsigned int  count;             /* How many times you've cycles the
					   * axes. */

  /* layer[U_FIELD]=lu; */
  /* layer[V_FIELD]=lv; */

  /* For each U triangle (then for each V triangle) determine the
   * coefficients of the implicit form of the plane through each
   * triangle. That is, find a,b,c,d in a*x + b*y + c*z + d = 0. This is done
   * _once only_ for each triangle. */

  for (field=0;field<=1;field++) { /* for either field */ 

    t_num[U_FIELD] = 0;  /* It doesn't matter in this loop, but use U_FIELD */
    t_pos[U_FIELD] = 0;

    while( (t_cube_edge[U_FIELD][0] = S->triangle_list[field][t_pos[U_FIELD]]) 
	   != END) { /* 1123 */

      /* t0, t1 and t2 are labels for the 3 triangle vertices.  
       * te0, te1, te2 are labels for the vectors along each triangle edge. */

      t_cube_edge[U_FIELD][1] = S->triangle_list[field][t_pos[U_FIELD]+1];
      t_cube_edge[U_FIELD][2] = S->triangle_list[field][t_pos[U_FIELD]+2];

      /* Calculate vectors along 2 of the edges of the triangle. Remember them
       * for later. */

      SUB(t_edge_vect[field][t_num[U_FIELD]][0],
	  S->vertex[field][t_cube_edge[U_FIELD][1]], 
	  S->vertex[field][t_cube_edge[U_FIELD][0]]);

      SUB(t_edge_vect[field][t_num[U_FIELD]][1],
	  S->vertex[field][t_cube_edge[U_FIELD][2]], 
	  S->vertex[field][t_cube_edge[U_FIELD][0]]);

      SUB(t_edge_vect[field][t_num[U_FIELD]][2],
	  S->vertex[field][t_cube_edge[U_FIELD][2]], 
	  S->vertex[field][t_cube_edge[U_FIELD][1]]);

      te0 = t_edge_vect[field][t_num[U_FIELD]][0];   /* te0 = t1 - t0 */
      te1 = t_edge_vect[field][t_num[U_FIELD]][2];   /* te1 = t2 - t0 */
                                               /* te2 = t2 - t1  which we may
						* need later */

      plane[field][t_num[U_FIELD]][A] = te0[Y]*te1[Z] - te0[Z]*te1[Y];
      plane[field][t_num[U_FIELD]][B] = te0[Z]*te1[X] - te0[X]*te1[Z];
      plane[field][t_num[U_FIELD]][C] = te0[X]*te1[Y] - te0[Y]*te1[X];    
      plane[field][t_num[U_FIELD]][D] =
	-(S->vertex[field][t_cube_edge[U_FIELD][1]][X]
                                             *plane[field][t_num[U_FIELD]][A] +
	  S->vertex[field][t_cube_edge[U_FIELD][1]][Y]
                                             *plane[field][t_num[U_FIELD]][B] +
	  S->vertex[field][t_cube_edge[U_FIELD][1]][Z]
                                             *plane[field][t_num[U_FIELD]][C]);

      t_num[U_FIELD] += 1;        
      t_pos[U_FIELD] += 3;
    } /* while t_cube_edge .. != END */
  } /* for field */

  /* Now look at each pair of triangles from the U and V cubes. */

  t_num[U_FIELD] = 0;           /* Triangle number (0-4) */
  t_pos[U_FIELD] = 0;

  while ( (t_cube_edge[U_FIELD][0] = S->triangle_list[U_FIELD][t_pos[U_FIELD]]) 
	  != END) { /* 1173 */

    /* For each triangle in the u cube */

    t_cube_edge[U_FIELD][1] = S->triangle_list[U_FIELD][t_pos[U_FIELD]+1];
    t_cube_edge[U_FIELD][2] = S->triangle_list[U_FIELD][t_pos[U_FIELD]+2];

    t_num[V_FIELD] = 0;
    t_pos[V_FIELD] = 0;

    while ( (t_cube_edge[V_FIELD][0] = S->triangle_list[V_FIELD][t_pos[V_FIELD]]) 
	    != END) { /* 1184 */
      /* For each triangle in the v cube */

      t_cube_edge[V_FIELD][1] = S->triangle_list[V_FIELD][t_pos[V_FIELD]+1];
      t_cube_edge[V_FIELD][2] = S->triangle_list[V_FIELD][t_pos[V_FIELD]+2];
      
      /* -----------------------------------------------------
       *  Start of code that find the common line of 2 planes. 
       * ----------------------------------------------------- */

      /* Find the common line :   x[] = p[] + t*dir[] */

      u_plane = plane[U_FIELD][t_num[U_FIELD]];
      v_plane = plane[V_FIELD][t_num[V_FIELD]];

      dir[X] = u_plane[B]*v_plane[C] - v_plane[B]*u_plane[C];
      dir[Y] = u_plane[C]*v_plane[A] - v_plane[C]*u_plane[A];
      dir[Z] = u_plane[A]*v_plane[B] - v_plane[A]*u_plane[B];

      den = DOT(dir,dir);

      if (EQUAL_ZERO(den)) {
	/* Planes are parallel so no intersect. 
	 * Overlapping co-planar triangles will be discarded here. */
      } else {

	dc = u_plane[D]*v_plane[C] - u_plane[C]*v_plane[D];
	db = u_plane[D]*v_plane[B] - u_plane[B]*v_plane[D];
	ad = u_plane[A]*v_plane[D] - v_plane[A]*u_plane[D];
	
	p[X] =  (dir[Y]*dc - dir[Z]*db) / den;
	p[Y] = -(dir[X]*dc + dir[Z]*ad) / den;
	p[Z] =  (dir[X]*db + dir[Y]*ad) / den;
	
	/* In fact p is the closest point on the line to the origin. But we
	 * do not use that information here. */
	
	/* --------------------------------------------------------- 
	 * Start of code that finds subset of line in both triangles 
	 * --------------------------------------------------------- */
	
	/* We find where the common line cuts the triangle edges. This will
	 * give us either 0, 2 or 4 points for sure. For the triangles to
	 * intersect, we must have 2 points from the U edges and 2 points
	 * from the V edges. We can stop looking as soon as these conditions
	 * are not met.
	 *
	 * By construction, the common line lies in the same plane as both
	 * triangles. Thus our arithmetic need only be in 2 dimension. We must
	 * be careful in the case where the triangle is normal to an
	 * axis. This will happen rarely under FPA, but is most likely at the
	 * initial condition where scalar values are set by the user.  You
	 * might want to see if there is a better way of handling such
	 * cases. I don't think the method I use here is the best. */

	points_found = 0;

	for (field=0;field<=1;field++) {  /* for either field */ 

	  SUB(p_t0, S->vertex[field][t_cube_edge[field][0]], p);

	  tri_edge = 0;
	  do /* for first 2 edges */ {

	    te = t_edge_vect[field][t_num[field]][tri_edge];
	    /* "te" <=> "triangle edge" */
				      
	    den = 0.0;
	    dim1 = X; dim2 = Y;
	    count = 0;
	    while (EQUAL_ZERO(den) && count < 3) {
	      dim1 = next_axis[dim1];
	      dim2 = next_axis[dim2];
	      den = te[dim2]*dir[dim1] - te[dim1]*dir[dim2];
	      count ++;
	    }
	    
	    if (EQUAL_ZERO(den)) {
	      /* Common line and triangle edge are parallel. No intersect. */
	      /* This happens quite a lot at the start of the simulation
	       * because scalar values are still close to their initial
	       * values. Holes can be seen. Later on, this doesn't happen
	       * at all. */
	    } else {
	      s = (dir[dim2]*p_t0[dim1] - dir[dim1]*p_t0[dim2]) / den;

	      if ((0.0 <= s) && (s <= 1.0)) {
		/* Intersection is on the triangle edge. If its actually on
		 * a S->vertex then it could still form the end of a filament
		 * segment (unlikely though). */

		r = (te[dim2]*p_t0[dim1] - te[dim1]*p_t0[dim2]) / den;
		/* This r is used later to give us an ordering on the points
		 * found. */

		point_r[points_found] = r;
		points_found++;
	      }
	    }
	    tri_edge += 1;
	  } while (tri_edge == 1);

	  /* We've done edges 0 and 1. If there have been no points found yet
	   * then the common line cannot pass through this triangle. If we've
	   * found 2 then we need to do no more and can look at the V
	   * triangle.  If we've found 1 then then the other points will be
	   * on the last edge. */

	  if ((points_found == 1) || (points_found == 3)) {
	  
	    /* This whole loop is repeated for the V field and is the reason
	     * this needs to be done for the 3rd edge of the V triangle. */

	    /* We know for sure that we're going to find another point on the
	       3rd edge. */

	    te = t_edge_vect[field][t_num[field]][2];

	    den = 0.0;
	    dim1 = X; dim2 = Y;
	    count = 0;
	    while (EQUAL_ZERO(den) && count < 3) {
	      dim1 = next_axis[dim1];
	      dim2 = next_axis[dim2];
	      den = te[dim2]*dir[dim1] - te[dim1]*dir[dim2];
	      count ++;
	    }

	    if (EQUAL_ZERO(den)) {
	      /* This should never happen. Test anyway. */
	    } else {
	      p_t1[dim1] = S->vertex[field][t_cube_edge[field][1]][dim1]- p[dim1];
	      p_t1[dim2] = S->vertex[field][t_cube_edge[field][1]][dim2]- p[dim2];
	      /* We don't need the remaining dimension for the calculation. */
	    
	      s = (dir[dim2]*p_t1[dim1] - dir[dim1]*p_t1[dim2]) / den;
	      if ((0.0 <= s) && (s <= 1.0)) {
		/* This _should_ always happen. Test anyway. */
		r = (te[dim2]*p_t1[dim1] - te[dim1]*p_t1[dim2]) / den;
	    
		point_r[points_found] = r;
		points_found++;
	      }
	    }
	  }

	  /* Found all the points we are going to on the U triangle. Now
	   * repeat for the V triangle as long as we found exactly 2
	   * points. */

	  if (points_found!=2) break;
	}

	if (points_found == 4) {
	  if (Find_two_pts(point, point_r, p, dir)) {
	    /* -----------------------------------------------------
	     * A piece of the filament has been found. Plot or write 
	     * it as necessary.  
	     * We know if this function has been called then at least
	     * "show_filament" is to be assumed nonzero.
	     * ----------------------------------------------------- */
	    Register_segment(point);
	  }
	}
      } /* end of dealing with common line */
      
      /* Try next V triangle. */
      t_num[V_FIELD] += 1;
      t_pos[V_FIELD] += 3;
    }
    
    /* Try next U triangle. */
    t_num[U_FIELD] += 1;
    t_pos[U_FIELD] += 3;
  } /* while t_cube_edge ... != END */
} /* Filament_in_cube */
/* ========================================================================= */

static Bool Find_two_pts (GLReal point[2][3], GLReal point_r[4], GLReal p[3], GLReal dir[3])
{
  /* This function takes the 4 r's in point_r[4] which give the intersections
   * of the common plane line with the triangle edges. It then draws the
   * subset of the common line which lies within both triangles.
   *
   * It is assumed that the r's have been correctly stored in point_r[4]. The
   * first 2 entries will be the U intersections, the 2nd two will be the
   * V. They probably will not be in order. 
   *
   * There must be a quicker way to do this, surely? I've tried Batcher's
   * optimal shuffle but the method below uses less 'ifs' I think because
   * this is a special case. We only need the middle 2 points which should
   * simplify things. */

  unsigned int u_min, v_min;
  unsigned int second_position, third_position;
  unsigned int higher_than_second, also_higher_than_second;

  /* Find order of the U pair */
  if (point_r[1] > point_r[0]) u_min = 0;
  else                         u_min = 1;

  /* Find order of the V pair */
  if (point_r[3] > point_r[2]) v_min = 0;
  else                         v_min = 1;

  /* Find minimum overall */
  if (point_r[v_min+2] > point_r[u_min]) {

    /* Lowest value is the lowest U. The next highest must be the lowest V for
     * there to be a common region. This will happen if the lowest V is lower
     * than the highest U. */

    second_position = v_min+2;
    higher_than_second = (1-u_min);
    also_higher_than_second = (1-v_min)+2;
  } 
  else {

    /* Lowest value is the lowest V. The next highest must be the lowest U for
     * there to be a common region. This will happen if the lowest U is lower
     * than the highest V. */

    second_position = u_min;
    higher_than_second = (1-v_min)+2;
    also_higher_than_second = (1-u_min);
  }

  if (point_r[second_position] < point_r[higher_than_second]) {

    /* There is a line but we don't know which is the third position. */

    if (point_r[higher_than_second] < point_r[also_higher_than_second]) {
      third_position = higher_than_second;
    }
    else {
      third_position = also_higher_than_second;
    }

    /* set the points */

    POINT_ALONG_LINE(point[0], p, point_r[second_position], dir);
    POINT_ALONG_LINE(point[1], p, point_r[third_position], dir);
    return TRUE;
  } else {
    return FALSE;
  }	
} /* Find_two_pts */
/* ========================================================================= */

/* Given a cube edge, relative_edge_pos tells you how to get from (ii,jj,kk)
 * to the end of the edge which is furthest along the edge's axis. The entry
 * for edge 1 tells you that you need to move down the Y axis from (ii,jj,kk)
 * and arrive at (ii,jj-1,kk). It also tells you that you do not need to move
 * along either the X or the Z axis. This information is used by
 * Contour_normal_intersect() which calculates the intersection of the
 * iso-surface with any given cube edge. The fact that the cube S->vertex "is
 * further along the edges axis" is important because then you know that the
 * intersect always lies behind you. This is why there is always at least one
 * '0' in each entry in the table.
 *
 * Note that in fact the 1's are replaced by "inc" in Set_relative_edge_pos()
 * This is because the graphics resolution is variable so you might want to
 * move along an axis by more than one numerical cube length. Doing this
 * saves some multiplications in Contour_normal_intersect(). */

static const unsigned int unit_relative_edge_pos[13][3] = {
  {0,0,0}, {0,1,0}, {0,0,0}, {0,0,0}, {0,0,1}, {1,1,0},
  {1,0,0}, {1,0,0}, {1,0,1}, {0,1,1}, {0,1,0}, {0,0,1}, {0,0,0} };

static void Set_relative_edge_pos (STR *S,int inc)
{
  int p, q;
  for (p=1; p<=12; p++) {
    for (q=0; q<=2; q++) {
      S->relative_edge_pos[p][q] = inc*unit_relative_edge_pos[p][q];
    }
  }
} /* Set_relative_edge_pos */
/* ========================================================================= */

static GLReal Contour_normal_intersect (STR *S, int layer, real contour, CubeEdge edge, GLReal *normal) 
{
  DEVICE_CONST(int,nx);
  DEVICE_CONST(int,ny);
  DEVICE_CONST(int,nz);
  DEVICE_CONST(int,nv);
  DEVICE_CONST(GLReal,marching_h);
  DEVICE_ARRAY(GLReal,pos);
  DEVICE_CONST(unsigned int,ii);
  DEVICE_CONST(unsigned int,jj);
  DEVICE_CONST(unsigned int,kk);
  DEVICE_CONST(unsigned int,inc);
  /* This function returns the position of the intersection of a contour with
   * the cube edge specified. It is assumed that the edge specified contains
   * a contour otherwise a divide by zero may occur. You can't go wrong as
   * long as you use the cube index to point to the edge table which will
   * tell you which edges have an intersect.
   *
   * The other field components at the intersect can be determined by the
   * calling routine because they are the same as the position of the cube
   * edge.
   *
   * May 2007, added computation of S->normals if normal is not NULL.  Normals
   * are found from the field gradient (which is normal to the field contour).
   * Gradients are computed by finite differences on the computational grid at
   * the two cube vertices on the current edge.  Then the same linear
   * interpolation used to find the contour intersect (triangle vextex) is
   * used to interpolate the gradient to the triangle S->vertex.  
   *
   */

  unsigned int  begin [3],
                end   [3];
  GLReal        ratio, norm;
  int k;

  /* Set end to be the position of the node in the numerical volume which is
   * furthest along the axis on the specified edge specified. See declaration
   * of relative_edge_pos. See Marching_cubes(). */

  end[X] = ii-S->relative_edge_pos[edge][X];
  end[Y] = jj-S->relative_edge_pos[edge][Y];
  end[Z] = kk-S->relative_edge_pos[edge][Z];

  /* The other end of the edge differs only along its axis by inc. */

  begin[X] = end[X];
  begin[Y] = end[Y];
  begin[Z] = end[Z];
  begin[axis[edge]] -= inc;

  /* Calculate the ratio that the intersect occurs along the cube edge from
   * the end. */

  ratio = (contour-Fields(layer,end[X],end[Y],end[Z]))
    / ( Fields(layer,begin[X],begin[Y],begin[Z])
       -Fields(layer,end[X],end[Y],end[Z]) );

  /* For the contour intersection (triangle S->vertex), all that remains is the
   * return statement at the end.
   *
   * If we are to compute normal also.  Use finite differences to compute
   * gradients at relevant cube vertices and interpolate to find the gradient
   * at triangle S->vertex. */

  if (normal != NULL) {

    normal[0] = 
      (  Fields(layer,begin[X]+1,begin[Y],begin[Z]) 
	 - Fields(layer,begin[X]-1,begin[Y],begin[Z]) ) * ratio
      + 
      (  Fields(layer,end[X]+1,end[Y],end[Z]) 
	 - Fields(layer,end[X]-1,end[Y],end[Z]) ) * (1. - ratio);

    normal[1] = 
      (  Fields(layer,begin[X],begin[Y]+1,begin[Z]) 
	 - Fields(layer,begin[X],begin[Y]-1,begin[Z]) ) * ratio
      + 
      (  Fields(layer,end[X],end[Y]+1,end[Z]) 
	 - Fields(layer,end[X],end[Y]-1,end[Z]) ) * (1. - ratio);

    normal[2] = 
      (  Fields(layer,begin[X],begin[Y],begin[Z]+1) 
	 - Fields(layer,begin[X],begin[Y],begin[Z]-1) ) * ratio
      + 
      (  Fields(layer,end[X],end[Y],end[Z]+1) 
	 - Fields(layer,end[X],end[Y],end[Z]-1) ) * (1. - ratio);

    /* Normalize.  The minus sign is because we want the normal in the
       direction of minus the gradient */

    norm = -sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    for(k=0;k<3;k++) normal[k] /= norm;
  }

  /* Now return the contour intersection (which gives the non-trivial
     coordinate of the current triangle S->vertex. */
  return (pos[axis[edge]] - marching_h*ratio);
} /* Contour_normal_intersect */
/* ========================================================================= */

static void Set_draw_mode_off (STR *S)
{
  DEVICE_VAR(int,draw_mode);
  if ((*draw_mode)==OFF) {
    MESSAGE("off again? this shouldn't be!\n");
  }
  glEnd();
  gdebug(4);
  (*draw_mode)=OFF;
} /* Set_draw_mode_off */
/* ========================================================================= */

static void Set_draw_mode_lines (STR *S,GLReal lwidth)
{
  DEVICE_CONST(int,clipping);
  DEVICE_VAR(int,draw_mode);
  DEVICE_VAR(GLReal,last_lwidth);
  if((*draw_mode)==LINES && lwidth==(*last_lwidth)) return;

  if((*draw_mode)!=OFF) {
    glEnd();
    gdebug(4);
  }
  if(clipping) glDisable(GL_CLIP_PLANE0);  /* we never clip lines */
  glDisable(GL_LIGHTING);                  /* we don't shine light on lines */
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  mat_noemission); /* nor let them emit ? */

  glLineWidth(lwidth);
  glBegin(GL_LINES);
  (*draw_mode)=LINES;
  (*last_lwidth)=lwidth;
} /* Set_draw_mode_lines */
/* ========================================================================= */

static void Set_draw_mode_curves (STR *S,GLReal lwidth,GLReal emission)
{
  DEVICE_CONST(int,clipping);
  DEVICE_VAR(int,draw_mode);
  /* if((*draw_mode)==CURVES) return; - continuous chains do not concatenate */

  if((*draw_mode)!=OFF) {
    glEnd();
    gdebug(4);
  }
  if(clipping) glDisable(GL_CLIP_PLANE0);  /* we never clip lines */
  glDisable(GL_LIGHTING);                  /* we don't shine light on lines */
  Set_filament_emission(S, emission);	   /* user says if filament emits */
  glLineWidth(lwidth);                     /* and how thick it is */
  glBegin(GL_LINE_STRIP);
  (*draw_mode)=CURVES;
} /* Set_draw_mode_curves */
/* ========================================================================= */

static void Set_filament_emission (STR *S, GLReal emission)
{
  DEVICE_CONST(real,flm_r);
  DEVICE_CONST(real,flm_g);
  DEVICE_CONST(real,flm_b);
  GLfloat flm_emission[4];
  flm_emission[0]=flm_r*emission;
  flm_emission[1]=flm_g*emission;
  flm_emission[2]=flm_b*emission;
  flm_emission[3]=1; /* apparenly this value is ignored anyway */
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  flm_emission);
} /* Set_filament_emission */

static void Set_draw_mode_triangles (STR *S)
{
  DEVICE_CONST(int,clipping);
  DEVICE_VAR(int,draw_mode);
  if((*draw_mode)==TRIANGLES) return;

  if((*draw_mode)!=OFF) {
    glEnd();
    gdebug(4);
  }

  if(clipping) {
    glEnable(GL_CLIP_PLANE0); 
    message(5,"clipping enabled in Set_draw_mode_triangles\n");
  } else {
    glDisable(GL_CLIP_PLANE0);
    message(5,"clipping disabled in Set_draw_mode_triangles\n");
  }
  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,  mat_emission);
  glBegin(GL_TRIANGLES);
  (*draw_mode)=TRIANGLES;
} /* Set_draw_mode_triangles */
/* ========================================================================= */

static void Set_draw_mode_fan (STR *S,GLReal emission)
{
  DEVICE_CONST(int,clipping);
  DEVICE_VAR(int,draw_mode);
  /* if((*draw_mode)==FAN) return; - fans do not concatenate */

  if((*draw_mode)!=OFF) {
    glEnd();
    gdebug(4);
  }

  if(clipping) glDisable(GL_CLIP_PLANE0);  /* we better not clip balls */
  glDisable(GL_LIGHTING);                  /* nor shine light on them */
  Set_filament_emission(S, emission);	   /* user says if filament emits */
  glEnable(GL_LIGHTING);
  glBegin(GL_TRIANGLE_FAN);
  (*draw_mode)=FAN;
} /* Set_draw_mode_fan */
/* ========================================================================= */

static void Set_draw_mode_strip (STR *S,GLReal emission)
{
  DEVICE_CONST(int,clipping);
  DEVICE_VAR(int,draw_mode);
  /* if((*draw_mode)==STRIP) return; - strips do not concatenate */
  if((*draw_mode)!=OFF) {
    glEnd();
    gdebug(4);
  }
  if(clipping) glDisable(GL_CLIP_PLANE0);  /* we better not clip balls */
  glDisable(GL_LIGHTING);                  /* nor shine light on them */
  Set_filament_emission(S, emission);	   /* user says if filament emits */
  glEnable(GL_LIGHTING);
  glBegin(GL_QUAD_STRIP);
  (*draw_mode)=STRIP;
} /* Set_draw_mode_strip */
/* ========================================================================= */

void Marching_ini (STR *S)
{
  /* DEVICE_ARRAY(CubeEdge[MAX_TRIANGLE_LIST_LENGTH],S->triangle_table); */
  /* DEVICE_ARRAY(CubeEdge[MAX_EDGE_LIST_LENGTH],S->edge_table); */
  /* DEVICE_ARRAY(CubeEdge *,S->triangle_list); */
  /* DEVICE_ARRAY(GLReal[13][3],S->vertex); */
  #include "ezgenerators.h"
  
  CubeIndex  index;
  int        gen;    

  /*  Initialize the triangle table to be empty  */
  for (index=0; index<=255; index++) {
    S->triangle_table[index][0] = END;
  }

  /* Copy all the basic 14 recipes into the triangle table and do all 24
   * orientations for each one. This is extremely inefficient. However, the
   * code is concise and easier to debug. Most entries are generated in more
   * than one way (and an error is generated if any inconsistencies
   * arise). Thus we can be confident our table is correct. This is only done
   * once, so efficiency is not an issue. */

  for (gen=0; gen<NUM_GENERATORS; gen++) {
    index = table_generators[gen][0];
    strcpy((char *)S->triangle_table[index], (char *)table_generators[gen]+1);
    Do_all_orientations(S, table_generators[gen][0]);
  }

  /* Generate the edge table from the triangle table by deleting all the
   * repeated edges. */

  Generate_edge_table(S);
  
  /*  Display the triangle table  */
  if(0) {
    MESSAGE ("Triangle table and corresponding edge table:\n");
    MESSAGE ("--------------------------------------------\n");
    for (index=0; index<=255; index++) { 
      MESSAGE ("\ntriangle_table[%3d] = ",index); 
      Print_list(S->triangle_table[index]);
      MESSAGE ("\n    S->edge_table[%3d] = ",index);
      Print_list(S->edge_table[index]);
    }
    MESSAGE ("\n"); 
  }
  /* I have checked the generated table against Lorensen's table which he
   * includes in the Visualization Toolkit Version 1.0. Randomly comparing
   * both tables and specifically looking at 'awkward' cases confirms that
   * they are essentially the same. He uses a different labeling scheme and
   * makes different topologically arbitrary choices to his own paper. He
   * uses a similar modification to the case table to deal with holes. I
   * follow the labeling scheme in his paper and the modifications suggested
   * by Montani et al.  */
  
} /* Marching_ini */
/* ========================================================================= */

static CubeIndex Rotate_edge_list (STR *S,CubeIndex index, int direction)
{
  /* This takes an entry in the S->triangle_table and generates the new list of
   * triangles obtained by rotating the cube around either the X or the Y
   * axis. */

  /* To generate the look up table, we simply rotate each of the above
 * generators through all 24 orientations. Under a rotation about the X axis,
 * x_mod_edge tells you where each edge is mapped. x_mod_S->vertex tells you
 * where each cube S->vertex is mapped. Similarly for a rotation about the Y
 * axis. */
  /* DEVICE_ARRAY(CubeEdge[MAX_TRIANGLE_LIST_LENGTH],S->triangle_table); */
  /* DEVICE_ARRAY(CubeEdge[MAX_EDGE_LIST_LENGTH],S->edge_table); */
  /* DEVICE_ARRAY(CubeEdge *,S->triangle_list); */
  /* DEVICE_ARRAY(GLReal[13][3],S->vertex); */
  

  CubeEdge   x_mod_edge[13]  = {0,3,12,7,11,1,10,5,9,4,2,8,6};
  CubeVertex x_mod_vertex[9] = {0,8,4,64,128,1,2,32,16};

  CubeEdge   y_mod_edge[13]  = {0,9,4,11,8,10,2,12,6,5,1,7,3};
  CubeVertex y_mod_vertex[9] = {0,16,1,8,128,32,2,4,64}; 

/* ------------------------------------------------------------------------- */

  CubeEdge    *mod_edge;
  CubeVertex  *mod_vertex;
  CubeIndex    new_index = 0;
  CubeEdge     new_edge_list[MAX_TRIANGLE_LIST_LENGTH];
  int          i;

  if (direction == 'x') {
    mod_edge   = x_mod_edge;
    mod_vertex = x_mod_vertex;
  }
  else {
    mod_edge   = y_mod_edge;
    mod_vertex = y_mod_vertex;
  }

  for (i=0; i<=7; i++) {
    if (index & (1<<i)) new_index |= mod_vertex[i+1];
  }

  i = 0;
  while (S->triangle_table[index][i] != END) {
    new_edge_list[i] = mod_edge[S->triangle_table[index][i]];
    i++;
  }
  new_edge_list[i] = END;

  if (S->triangle_table[new_index][0] == END) {
    /*  The first time this entry has been defined  */
    strcpy((char *)S->triangle_table[new_index], (char *)new_edge_list);
  }

  else {
    if (strcspn ((char *)S->triangle_table[new_index], (char *)new_edge_list) 
	!= 0) {
      /* Triangle lists are considered equivalent if they have the same
       * elements (even if they are in a different order).  I know this is a
       * weak test but its better than nothing. */
      fprintf (stderr,"Logic error generating lookup table!\n");
    }
  }

  return (new_index);
} /* Rotate_edge_list */
/* ========================================================================= */

static CubeIndex Rotate_edge_list_x (STR *S, CubeIndex index)
{
  return (Rotate_edge_list(S,index,'x'));
}
/* ========================================================================= */

static CubeIndex Rotate_edge_list_y (STR *S, CubeIndex index)
{
  return (Rotate_edge_list(S,index,'y'));
}
/* ========================================================================= */

static void Rotate_diagonal (STR *S,CubeIndex index)
{
  CubeIndex  lindex = index;
  int        i;

  for (i=0; i<=2; i++) lindex = Rotate_edge_list_x(S,Rotate_edge_list_y(S,lindex));
  if (lindex != index)
    fprintf(stderr,"Assertion failed in Rotate_diagonal\n");
}
/* ========================================================================= */

static void Do_all_orientations (STR *S, CubeIndex index)
{
  /* Given an entry in the S->triangle_table at location index, this function
   * rotates the cube through all 24 orientations and writes the resulting
   * edge lists to the edge table. */

  CubeIndex  lindex = index;   /* local copy of index */
  int        i;

  for (i=0; i<=3; i++) {
    Rotate_diagonal(S, lindex);
    lindex = Rotate_edge_list_x(S, lindex);
  }
  if (lindex != index)
    fprintf(stderr,"Assertion 1 failed in Do_all_orientations\n");
  
  lindex = Rotate_edge_list_y(S, Rotate_edge_list_y(S, lindex));
  for (i=0; i<=3; i++) {
    Rotate_diagonal(S, lindex);
    lindex = Rotate_edge_list_x(S, lindex);
  }
  lindex = Rotate_edge_list_y(S, Rotate_edge_list_y(S, lindex));
  if (lindex != index)
    fprintf(stderr,"Assertion 2 failed in Do_all_orientations\n");
} /* Do_all_orientations */
/* ========================================================================= */

static void Generate_edge_table (STR *S)
{
  /* This takes the triangle table and generates the S->edge_table by deleting
   * all the repeated vertices. */

  /* DEVICE_ARRAY(CubeEdge[MAX_TRIANGLE_LIST_LENGTH],S->triangle_table); */
  /* DEVICE_ARRAY(CubeEdge[MAX_EDGE_LIST_LENGTH],S->edge_table); */
  /* DEVICE_ARRAY(CubeEdge *,S->triangle_list); */
  /* DEVICE_ARRAY(GLReal[13][3],S->vertex); */

  CubeIndex  index;
  int        pos1, pos2;
  Bool       got_already;
  CubeEdge   new_edge;

  for (index = 0; index<=255; index++) {
    
    /* go through this entry and write any new edges to the edge table */
    pos1 = 0;
    S->edge_table[index][0] = END;

    while ((new_edge = S->triangle_table[index][pos1]) != END) {

      pos2 = 0;
      got_already = FALSE;

      while (S->edge_table[index][pos2] != END) {
	if (new_edge == S->edge_table[index][pos2]) got_already = TRUE;
	pos2++;
      }

      if (!got_already) {
	S->edge_table[index][pos2] = new_edge;
	S->edge_table[index][pos2+1] = END;
      }

      pos1++;
    }
  }
} /* Generate_edge_table */
/* ========================================================================= */

static void Print_list (CubeEdge *p)
{
  /* Simple routine which prints a list of edges as stored in the triangle
   * table and the edge table. */

  if (*p == END) return;                       /* Empty list */
  while (*(p+1) != END) MESSAGE ("%d, ",*p++);  /* Print all but the last */
  MESSAGE ("%d",*p);                            /* Print the last */
} /* Print_list */
/* ========================================================================= */

/* Writing to filament file is bufferized in case the decision is to be made later */
static void Filament_print (STR *S, char *fmt, ...)
{
  DEVICE_ARRAY(char,filament_buffer);
  DEVICE_CONST(char *,filament_buffer_end);
  DEVICE_VAR(char *,filament_buffer_current);
  size_t remains=filament_buffer_end-(*filament_buffer_current);
  size_t added;
  va_list argptr;
  if (!remains) return;
  va_start(argptr, fmt);
  added=vsnprintf((*filament_buffer_current), remains, fmt, argptr);
  va_end(argptr);
  if (added<remains)
    (*filament_buffer_current)+=added;
  else
    (*filament_buffer_current)+=remains;
} /* Filament_print */
/* ========================================================================= */

/* Reset filament buffer */
static void Filament_rewind (STR *S)
{
  DEVICE_ARRAY(char,filament_buffer);
  DEVICE_VAR(char *,filament_buffer_current);
  *filament_buffer_current=filament_buffer;
}
/* ========================================================================= */

#if FIBRES
/* Create the GL list of the fibres */
static int get_dir(STR *S,real x, real y, real z, real *kx, real *ky, real *kz, real *r, real *g, real *b);
void Make_fibres (STR *S)
{
  DEVICE_VAR(int,FIB_LIST);
  DEVICE_ARRAY(real,plot_length);
  DEVICE_CONST(int,fib_x0);
  DEVICE_CONST(int,fib_y0);
  DEVICE_CONST(int,fib_z0);
  DEVICE_CONST(int,fib_dx);
  DEVICE_CONST(int,fib_dy);
  DEVICE_CONST(int,fib_dz);
  DEVICE_CONST(int,imin);
  DEVICE_CONST(int,imax);
  DEVICE_CONST(int,jmax);
  DEVICE_CONST(int,kmax);
  DEVICE_CONST(int,fib_wt);
  DEVICE_CONST(int,flm_em);
  DEVICE_CONST(int,fib_len);
  DEVICE_CONST(real,fib_step);
  int dir, len;
  real x, startx;
  real y, starty;
  real z, startz;
  real kx, kx1, startkx;
  real ky, ky1, startky;
  real kz, kz1, startkz;
  real red, green, blue;
  real convert_to_graph=plot_length[0]/(imax-imin+1);

  if (!ANISOTROPY_ON || !fib_wt || !fib_len || !fib_step) {
    *FIB_LIST=0;
    return;
  }

  /* Make new graphics list for both resolutions. */
  if ((*FIB_LIST)==0) (*FIB_LIST)=glGenLists(1);
  message(4,"FIB_LIST=%d allocated\n",(*FIB_LIST));

  glNewList(*FIB_LIST, GL_COMPILE);

  for (startx=fib_x0; startx<=imax; startx+=fib_dx) {
    for (starty=fib_y0; starty<=jmax; starty+=fib_dy) {
      for (startz=fib_z0; startz<=kmax; startz+=fib_dz) {
	if (!get_dir(S,startx,starty,startz,&startkx,&startky,&startkz,&red,&green,&blue)) continue;
	for (dir=-1;dir<=1;dir+=2) {
	  x=startx; y=starty; z=startz;
	  kx=dir*startkx; ky=dir*startky; kz=dir*startkz;
	  Set_draw_mode_curves(S,fib_wt,flm_em);
	  GLCOLOR3(red,green,blue);
	  GLVERTEX3(convert_to_graph*x,convert_to_graph*y,convert_to_graph*z);
	  for (len=0;len<fib_len;len++) {
	    x+=kx*fib_step;
	    y+=ky*fib_step;
	    z+=kz*fib_step;
	    if (!get_dir(S,x,y,z,&kx1,&ky1,&kz1,&red,&green,&blue)) break;
	    GLCOLOR3(red,green,blue);
	    GLVERTEX3(convert_to_graph*x,convert_to_graph*y,convert_to_graph*z);
	    if ((kx*kx1+ky*ky1+kz*kz1)>0) {
	      kx=kx1; ky=ky1; kz=kz1;
	    } else {
	      kx=-kx1; ky=-ky1; kz=-kz1;
	    }
	  } /* for len */
	  Set_draw_mode_off(S);
	} /* for dir */
      } /* for startz */
    } /* for starty */
  } /* for startx */
      
  glEndList();
} /* Make_fibres */
/* ========================================================================= */

static int get_dir (STR *S,real x, real y, real z, real *kx, real *ky, real *kz, real *r, real *g, real *b)
{
  DEVICE_CONST(int,fib_rlayer);
  DEVICE_CONST(real,fib_rmin);
  DEVICE_CONST(real,fib_rmax);
  DEVICE_CONST(int,fib_glayer);
  DEVICE_CONST(real,fib_gmin);
  DEVICE_CONST(real,fib_gmax);
  DEVICE_CONST(int,fib_blayer);
  DEVICE_CONST(real,fib_bmin);
  DEVICE_CONST(real,fib_bmax);
  int ix, iy, iz, jx, jy, jz;
  real px, py, pz, pj;
  real rsum, gsum, bsum, kxsum, kysum, kzsum, sum;
  real rvalue, gvalue, bvalue;
  int status;

  /* Identify the relevant grid cell */
  ix=floor(x);
  iy=floor(y);
  iz=floor(z);
  /* Compute weighted average of required quantities at cell vertices */
  sum=rsum=gsum=bsum=kxsum=kysum=kzsum=0;
  for (jx=ix;jx<=ix+1;jx++) {
    px=1.0-fabs(x-jx);
    if (!px) continue;
    for (jy=iy;jy<=iy+1;jy++) {
      py=1.0-fabs(y-jy);
      if (!py) continue;
      for (jz=iz;jz<=iz+1;jz++) {
	pz=1.0-fabs(z-jz);
	if (!pz) continue;
	/* if (!isTissue(jx,jy,jz)) continues; /\* Ignore vertices in the void *\/ */
	status=isTissue(jx,jy,jz);
	/* printf ("%d %d %d -> %d\n",jx,jy,jz,status); */
	if (!status) continue; 		/* Ignore vertices in the void */
	pj=px*py*pz; 			/* The weight of this vertex */
	sum+=pj;
	if (fib_rlayer>=0) rsum+=pj*New[ind(jx,jy,jz,fib_rlayer)];
	if (fib_glayer>=0) gsum+=pj*New[ind(jx,jy,jz,fib_glayer)];
	if (fib_blayer>=0) bsum+=pj*New[ind(jx,jy,jz,fib_blayer)];
	kxsum+=pj*Geom[geom_ind(jx,jy,jz,GEOM_FIBRE_1)];
	kysum+=pj*Geom[geom_ind(jx,jy,jz,GEOM_FIBRE_2)];
	kzsum+=pj*Geom[geom_ind(jx,jy,jz,GEOM_FIBRE_3)];
      } /* for jz */
    } /* for jy */
  } /* for jz */
  if (!sum) return 0;			/* All vertices void */
  *kx=kxsum/sum;			/* Fibre direciton averages */
  *ky=kysum/sum;
  *kz=kzsum/sum;
  rvalue=rsum/sum;			/* Field values at the point */
  gvalue=gsum/sum;
  bvalue=bsum/sum;
  *r=(rvalue-fib_rmin)/(fib_rmax-fib_rmin);	/* Apply the colour scheme */
  *g=(gvalue-fib_gmin)/(fib_gmax-fib_gmin);
  *b=(bvalue-fib_bmin)/(fib_bmax-fib_bmin);
  return 1;
} /* get_dir */
/* ========================================================================= */

#endif

#endif
