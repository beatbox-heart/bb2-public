/**
 * Copyright (C) (2010-2020) Vadim Biktashev, Irina Biktasheva et al. 
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

/* bbxd: count non-empty subdomains for a set of decomposition of a given geometry */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define NOT(a) (0==(a))
#ifndef INT
#define INT long
#define DFMT "%ld"
#define INTFMT DFMT
#endif
#define GEOM_VOID   0
#define FFLUSH fflush
#define EXPECTED_ERROR ABORT
#define CALLOC(p,a,b) if(0==(p=Calloc(a,b))) ABORT("could not allocate %ld bytes for %s",1L*a*b,#p)
#define FREE(addr) if(addr) {free(addr);(addr)=NULL;}
#define ASSERT(p) { if(0==(p)) EXPECTED_ERROR("Assertion failed:\n %s\n",#p); }
#define gind(x,y,z) (((x)*ymax+(y))*zmax+(z))

typedef struct{
  int nx, ny, nz;
} Decomposition;

int Verbose=1;
int dim,ONE,TWO,TRI;
int x_offset,y_offset,z_offset;
INT xmax,ymax,zmax;
INT GEOMETRY_ON=1;
int *all_xmin, *all_xmax;	/* Local minima and maxima */
int *all_ymin, *all_ymax;	/*   of all coordinate     */
int *all_zmin, *all_zmax;	/*   axes' partitions      */
int num_subdoms;		/* Number of subdomains */
int num_active_procs;		/* Number of active (non-idle) processes. Could be < num_subdoms. */
unsigned char *gpoints;		/* "skeleton" of the full geometry   */

void ABORT(char *fmt, ...);
void MESSAGE(char *fmt, ...);
void *Calloc(size_t nitems,size_t size);
int getMeshDimensions (FILE *geomFile,INT *mesh_xmax,INT *mesh_ymax,INT *mesh_zmax, INT vd, INT vn);
int readGeomLine(FILE *geomFile,int *x,int *y,int *z, int *status,double *k1,double *k2,double *k3,double *depth,double *n);
void decomp_allocSubdomains(int NX,int NY,int NZ);
int decomp_checkSubdomains(int NX, int NY, int NZ);
void decomp_deallocSubdomains(void);

int main(int argc, char **argv) {
  int i;
  char *geometryname;
  FILE *geometry;
  int NX, NY, NZ;
  int NXmin, NXmax, NYmin, NYmax, NZmin, NZmax;
  int vd=-1;
  int vn=-1;

  if (argc!=5) ABORT("argv=%d should be 5\n",argc);

  geometryname=argv[1];
  if (1!=sscanf(argv[2],"%d",&Verbose)) ABORT("could not read '%s' as Verbose\n",argv[2]);
  if (Verbose) MESSAGE("/* Finding minimal enclosing box for geometry from %s ... */\n",
	  geometryname);
  if NOT(geometry=fopen(geometryname,"r")) ABORT("cannot open file %s for reading", geometryname);
  getMeshDimensions(geometry,&xmax,&ymax,&zmax,vd,vn);
  if (Verbose) MESSAGE("/* Grid size (%d x %d x %d)\n*/",xmax,ymax,zmax);
  if (xmax==2||ymax==2||zmax==2) {
    EXPECTED_ERROR("One or more dimensions (xmax,ymax,zmax) is equal to 2. Dimensions must be either 1 or >=3.");
  }
  TRI = zmax >=3 ? 1 : 0;
  TWO = ymax >=3 ? 1 : 0;
  ONE = xmax >=3 ? 1 : 0;
  dim = ONE + TWO + TRI;
  if(dim==1 && !ONE)		EXPECTED_ERROR("1-Dimensional simulations must be defined on the x axis.");
  if(dim==1 && GEOMETRY_ON)	EXPECTED_ERROR("Geometry requires a 2- or 3-dimensional simulation medium.");
  if(dim==2 && !(ONE && TWO))	EXPECTED_ERROR("2-Dimensional simulations must be defined on the x and y axes.");


  if (1!=sscanf(argv[3],"%d",&NXmin)) ABORT("could not read '%s' as NXmin\n",argv[3]);
  if (1!=sscanf(argv[4],"%d",&NXmax)) ABORT("could not read '%s' as NXmax\n",argv[4]);
  #define above 1.1
  #define below 0.9
  for (NX=NXmin;NX<=NXmax;NX++) {
    if (NX<1) continue;
    NYmax=NX*ymax/xmax*above;
    NYmin=NX*ymax/xmax*below;
    NZmax=NX*zmax/xmax*above;
    NZmin=NX*zmax/xmax*below;
    for (NY=NYmin;NY<=NYmax;NY++) {
      if (NY<1) continue;
      for (NZ=NZmin;NZ<=NZmax;NZ++) {
	if (NZ<1) continue;
	num_subdoms=NX*NY*NZ;
	if (Verbose>=2) MESSAGE("/* Domain decomposed as (%03d,%03d,%03d). */\n",NX,NY,NZ);
	
	decomp_allocSubdomains(NX,NY,NZ);
	num_active_procs = decomp_checkSubdomains(NX,NY,NZ);
	if (Verbose>=2) MESSAGE("/* %d processes, %d active, %d idle. */\n",
			     num_subdoms, num_active_procs, num_subdoms-num_active_procs
			     );
	decomp_deallocSubdomains();
	printf("%ld\t%ld\t%ld\t%ld\t%ld\n",(long)NX,(long)NY,(long)NZ,(long)num_subdoms,(long)num_active_procs);
	fflush(stdout);
      }
    }
  }
}
/*=================================================================================================*/
/*=================================================================================================*/


/************************/
void ABORT(char *fmt, ...)
{
  char s[4092], *p;
  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  exit(1);
}
/**************************/
void MESSAGE(char *fmt, ...)
{
  char s[4092], *p;
  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);
  return;
}
/**************************************/
void * Calloc(size_t nitems,size_t size)
{
  return(calloc(nitems, size));
}

/***************************************************************************************************/
int getMeshDimensions (FILE *geomFile,INT *mesh_xmax,INT *mesh_ymax,INT *mesh_zmax, INT vd, INT vn) {
  int first,
    x,y,z,
    x_in, y_in, z_in,
    min_x,max_x,min_y,max_y,min_z,max_z,
    geom_xlen, geom_ylen, geom_zlen;
  double k1, k2, k3;
  int status=GEOM_VOID;
  int *pstatus;
  double depth=0;
  double *pdepth;
  double n;
  double *pn;

  first = 1;
  rewind(geomFile);
  pstatus=(vd>=0)?NULL:(&status);
  pdepth=(vd>=0)?(&depth):NULL;
  pn=(vn>=0)?(&n):NULL;
  while (EOF!=readGeomLine(geomFile,&x,&y,&z,pstatus,&k1,&k2,&k3,pdepth,pn)) {

    if (depth>=0 || status != GEOM_VOID) {
      if (first) {
	min_x = max_x = x;
	min_y = max_y = y;
	min_z = max_z = z;
	first = 0;
      } else {
	if (x<min_x) min_x=x;
	if (x>max_x) max_x=x;
	if (y<min_y) min_y=y;
	if (y>max_y) max_y=y;
	if (z<min_z) min_z=z;
	if (z>max_z) max_z=z;
      }
    }
  }
  if (Verbose) MESSAGE("/* Geometry is contained in a box [%d:%d]x[%d:%d]x[%d:%d] */\n",min_x,max_x,min_y,max_y,min_z,max_z);
  if (vd>=0) {
    min_x--; max_x++;
    if (TWO) {min_y--; max_y++;}
    if (TRI) {min_z--; max_z++;}
    if (Verbose) MESSAGE("/* ... and extended to [%d:%d]x[%d:%d]x[%d:%d] */\n",min_x,max_x,min_y,max_y,min_z,max_z);
  }
  geom_xlen = ((max_x - min_x) + 1);
  geom_ylen = ((max_y - min_y) + 1);
  geom_zlen = ((max_z - min_z) + 1);
  
  x_offset = ((geom_xlen>1?1:0)-min_x);
  y_offset = ((geom_ylen>1?1:0)-min_y);
  z_offset = ((geom_zlen>1?1:0)-min_z);
  if (Verbose) MESSAGE("/* Geometry coordinates will be offset by (%d,%d,%d) */\n",x_offset,y_offset,z_offset);

  *mesh_xmax = geom_xlen + (geom_xlen>1?2:0); /* Num points + bounds. */
  *mesh_ymax = geom_ylen + (geom_ylen>1?2:0); /* Num points + bounds. */
  *mesh_zmax = geom_zlen + (geom_zlen>1?2:0); /* Num points + bounds. */
  
  FFLUSH(stdout);

  /* Now we need to create the gpoints array for subdomain allocation. */
  CALLOC(gpoints,(*mesh_xmax)*(*mesh_ymax)*(*mesh_zmax),sizeof(unsigned char));
  rewind(geomFile);
  pdepth=(vd>=0)?(&depth):NULL;
  pn=(vn>=0)?(&n):NULL;
  while (EOF!=readGeomLine(geomFile,&x_in,&y_in,&z_in,&status,&k1,&k2,&k3,pdepth,pn)) {
    x = x_in + x_offset;
    y = y_in + y_offset;
    z = z_in + z_offset;
    gpoints[gind(x,y,z)] = (status!=GEOM_VOID) ? 1 : 0;
  }
  return 1;
}


/*******************************************************************************************/
/* List of characters that can separate fields in a geometry line, other than simple blank */
/* Perhaps this should be extended to include all characters not used in numbers           */
char GEOMSEPS[]=",:;\t\r\b!$=";

/* Read buffer: should be distinct from qpp buffer as they are processed "in parallel" */
#define MAXSTRLEN 8192
static char buf[MAXSTRLEN];
int readGeomLine(FILE *geomFile,
			int *x,int *y,int *z, int *status,
			double *k1,double *k2,double *k3,
			double *depth, double *n) {
  char *p, *c;
  int rc;
  if (NULL==(fgets(buf,MAXSTRLEN,geomFile))) return EOF;
  /* Convert all separator characters to blanks */
  for (c=&(GEOMSEPS[0]);*c;c++) {
    while (NULL!=(p=strchr(buf,*c))) *p=' ';
  }
  if (depth==NULL) { /* old way: no depth value expected */
    rc=sscanf(buf,"%d %d %d %d %lf %lf %lf\n",x,y,z,status,k1,k2,k3);
  } else if (n==NULL) { /* newer way: expect depth value in place of status */
    rc=sscanf(buf,"%d %d %d %lf %lf %lf %lf\n",x,y,z,depth,k1,k2,k3);
    (*status)=((*depth)>0);
  } else { /* newest way: also the NBC nonhomogeneity */
    rc=sscanf(buf,"%d %d %d %lf %lf %lf %lf %lf\n",x,y,z,depth,k1,k2,k3,n);
    (*status)=((*depth)>0);
  }
  return rc;
}

/************************************************/
int decomp_checkSubdomains(int NX, int NY, int NZ)
{
  int ix,iy,iz; /* looping superindices */
  int x,y,z; /* indices within subdomain */
  int nonempty; /* flag for the current subdomain */
  int n;    /* number of nonempty subdomains -> new num_active_procs */
  ASSERT(num_subdoms == NX * NY * NZ);
  ASSERT(gpoints!=NULL);
  n=0;
  for (ix=0;ix<NX;ix++) {
    for (iy=0;iy<NY;iy++) {
      for (iz=0;iz<NZ;iz++) {
	nonempty=0;
	for (x=all_xmin[ix];x<=all_xmax[ix];x++) {
	  for (y=all_ymin[iy];y<=all_ymax[iy];y++) {
	    for (z=all_zmin[iz];z<=all_zmax[iz];z++) {
	      if (gpoints[gind(x,y,z)]) {nonempty=1; break;}
	    } /* for z */
	    if (nonempty) break;
	  } /* for y */
	  if (nonempty) break;
	} /* for x */
	n+=nonempty;
      } /* for iz */
    } /* for iy */
  } /* for ix */
  if (!n) EXPECTED_ERROR("No non-empty subdomains found in this geometry.");
  return n;
}

/************************************************/
void decomp_allocSubdomains(int NX,int NY,int NZ)
{
  int ix, iy, iz;
  /* Approximate length of subdomain along each axis.	*/
  int x_slice_size = (xmax-(ONE*2)) / NX;
  int y_slice_size = (ymax-(TWO*2)) / NY;
  int z_slice_size = (zmax-(TRI*2)) / NZ;
  
  /* Remainders from i[xyz] * [xyz]_slice_size */
  /* The remainder will be shared amongst the first processes on the axis */
  int x_remainder  = (xmax-(ONE*2)) % NX;
  int y_remainder  = (ymax-(TWO*2)) % NY;
  int z_remainder  = (zmax-(TRI*2)) % NZ;
  
  /* Calculate upper and lower bounds along x axis */
  CALLOC(all_xmin,NX,sizeof(int));
  CALLOC(all_xmax,NX,sizeof(int));
  for (ix=0;ix<NX;ix++) {
    all_xmin[ix] = ONE + ( ix    * x_slice_size) + (ix < x_remainder?  ix    :x_remainder);
    all_xmax[ix] = ONE + ((ix+1) * x_slice_size) + (ix < x_remainder? (ix+1) :x_remainder);
  }
  
  /* Calculate upper and lower bounds along y axis */
  CALLOC(all_ymin,NY,sizeof(int));
  CALLOC(all_ymax,NY,sizeof(int));
  for (iy=0;iy<NY;iy++) {
    all_ymin[iy] = TWO + ( iy    * y_slice_size) + (iy < y_remainder?  iy    :y_remainder);
    all_ymax[iy] = TWO + ((iy+1) * y_slice_size) + (iy < y_remainder? (iy+1) :y_remainder);
  }
  
  /* Calculate upper and lower bounds along z axis */
  CALLOC(all_zmin,NZ,sizeof(int));
  CALLOC(all_zmax,NZ,sizeof(int));
  for (iz=0;iz<NZ;iz++) {
    all_zmin[iz] = TRI + ( iz    * z_slice_size) + (iz < z_remainder?  iz    :z_remainder);
    all_zmax[iz] = TRI + ((iz+1) * z_slice_size) + (iz < z_remainder? (iz+1) :z_remainder);
  }
}

/************************************************/
void decomp_deallocSubdomains(void)
{
  FREE(all_zmax);
  FREE(all_zmin);
  FREE(all_ymax);
  FREE(all_ymin);
  FREE(all_xmax);
  FREE(all_xmin);
}
