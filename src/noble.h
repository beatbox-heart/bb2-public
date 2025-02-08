/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
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

#include "extern.h"
 /*#include "p2c.h" */

#define neqn            60   /* neqn is max no. of differential equations */
typedef double Arneq[neqn];
EXTERN Arneq _F;	/* dy/dt for each variable in array */
EXTERN Arneq _Y;	/* array of variables to be integrated */
EXTERN Arneq alpha;	/* alpha rate coefficients */
EXTERN Arneq beta;	/* beta rate coefficients */
EXTERN FILE *Filin1, *Filout3;


/*#include "nobleini.h"*/
 /*#include "p2c.h" */
#include "nobleftn.h"
#include "noblestu.h"

