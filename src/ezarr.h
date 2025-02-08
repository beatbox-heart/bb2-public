/**
 * Copyright (C) (2010-2018) Vadim Biktashev, Irina Biktasheva et al. 
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

char filamentname[MAXPATH];
FILE *filament;
char images[MAXPATH];		/* generated image will converted to such files */
char imagescode[BUFLEN];	/* calculate a number to make part of the filename */
pp_fn imagescompiled;		/* k_code for the filename number calculation */
real plot_length[3];		/* determine shape of the box */
real plot_ini[3];		/* determine shape of the box */
GLfloat light_position0[4];
GLfloat lmodel_ambient[4];
CubeEdge        triangle_table[256][MAX_TRIANGLE_LIST_LENGTH];
                                /* triangle_table is the marching cubes lookup
                                 * table containing triples of edges which define
                                 * each triangle.  */
CubeEdge        edge_table[256][MAX_EDGE_LIST_LENGTH];
                                /* Same as triangle_table but without repetitions
                                 * of edges. */
CubeEdge        *triangle_list[2];
                                /* triangle_list[U_FIELD] and
                                 * triangle_list[V_FIELD] are used when filament
                                 * plotting to find common lines of the
                                 * triangles. */
GLReal          vertex[2][13][3];
                                /* Triangle vertices in real world space for each
                                 * edge in each volume. When just surface drawing
                                 * we only use the U or V part of the array. When
                                 * filament drawing we use both. */
GLReal          normals[13][3];
                                /* Surface normals at the triangle vertices. */
GLReal          pos[3];
                                /* Position of the 3rd vertex of the current
			         * marching cube in graphics coordinates. */
GLfloat mouse_down_mx[16];
unsigned int relative_edge_pos[13][3];

/* window title */
char title[1024];		/* (short) window title template */
char titlecode[1024];         	/* calculate a number to make part of the title */
pp_fn titlecompiled;		/* k_code for the title number calculation */
char shorttitle[MAXWINDOWTITLE]; /* short window title filled */
char longtitle[MAXWINDOWTITLE];  /* same + flags and angles */
