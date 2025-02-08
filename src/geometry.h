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

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#ifdef _OWN
#undef _OWN
#define OWN
#endif
#include "extern.h"

#define UNIT_VECTOR_TOLERANCE 0.01

int isUnitVector (double v1, double v2, double v3, double tolerance);


/** \brief Gets medium dimensions from the geometry file.
	
    Gets dimensions of the containing box given coordinates of tissue points.
    Sets offsets for the reading of tissue points, as minima might not be 1.
    Calculated dimensions include a 1-point boundary wrapper.
    
    \author Ross McFarlane
    \date 2010-04-01
    \param geomFile is the pointer to geometry file.
    \param mesh_xmax the address where the recommended xmax should be stored.
    \param mesh_ymax the address where the recommended ymax should be stored.
    \param mesh_zmax the address where the recommended zmax should be stored.
    
    \return True (1) on success, 0 otherwise.
**/
/* output parameters mesh_xmax ... of type (INT *) are pointers to k-variables */
int getMeshDimensions(FILE *geomFile, INT *mesh_xmax,INT *mesh_ymax,INT *mesh_zmax,int padding);

/* Former side-effect of getMeshDimensions, now made separate */
void make_gpoints(FILE *geomFile);

/** \brief Loads geometry data.
    
    Loads geometry data from a comma-separated value file.
    Lines should be formatted as:
    x,y,z,status,fibre1,fibre2,fibre3
    where status is GEOM_VOID, GEOM_TISSUE or GEOM_BOUND and
    fibre values should be components of a unit vector for all non-void points.
    
    \author Ross McFarlane 
    \date 2009-09-23
    \param geomFile pointer to geometry file.
    \param geomFileName geometry file name.
    \param normaliseVectors indicates whether all vectors (fibre and boundary) should be normalised.
    \return True (1) on success, 0 otherwise.
**/
int populateMesh(FILE *geomFile, const char *geomFileName, int normaliseVectors, int checkVectors, FILE *outFile, const char *outFileName);

#endif /* end of include guard: _GEOMETRY_H_ */
