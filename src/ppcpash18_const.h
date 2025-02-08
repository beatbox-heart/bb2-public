/**
 * Copyright (C) (2010-2024) Vadim Biktashev, Irina Biktasheva et al. 
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

/* List of constants that could be used for tabulation of functions of V */

#define F  (96485.3415) /* Faraday constant (C/mol) */
#define R  (8314.472)   /* Gas constant (mJ/(mol*K)) */
#define T  (310.0)	/* temperature (K) */
#define RTF (R*T/F)	/* common subexpression (mV) */

#define Vc (8800.0e-15) /* cell volume (L) */
#define Vsr (583.73e-15)/* SR volume (L) */
#define Cm (9.87109e-11)/* cell membrane capacitance (F) */
#define CFV (Cm/(F*Vc)) /* common subexprresion (mM/mV) */

/* These must be constant to allow tabulation of ... */
#define Cao (1.8)	/* (mM) ... Xr1 */
#define Vh_hLate (87.61)	/* (mV) ... hL */
