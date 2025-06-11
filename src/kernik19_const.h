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

#define F  (96.4853415) /* Faraday constant (C/mmol) */
#define R  (8.314472)   /* Gas constant (J/(mol*K)=mJ/(mmol*K)) */
#define T  (310.0)	/* temperature (K) */
#define RTF ((R*T)/F)	/* common subexpression (mV) */
#define FRT (F/(R*T))
#define FFRT (F*FRT)

#define Vc_tenT (16404.0) /* && (um^3) */
#define VSR_tenT (1094.0) /* && (um^3) */
#define V_tot (3960.0) /* && (um^3) */
#define Cm (60.0) /* && cell membrane capacitance (pF) */
/* #define CFV (Cm/(F*Vc)) /\* common subexprresion (mM/mV) *\/ - not used */

#define NatoK_ratio (0.491) /* && */

/* Irel */
#define MaxSR (15.0) /* && */
#define MinSR (1.0) /* && */  
#define ec50SR (0.45) /* && (mM) */
/* #define ks (12.5) /\* && (1/ms) *\/ - this determines magnitude of j_rel so promoted to parameter */
#define kiCa (18.495) /* && (1/(mM ms)) = 54*0.3425 */
#define kim (5.571e-4)  /* && (1/ms) = 0.001*0.5571 */
#define koCa (643751.68) /* && (1/(mM^2 ms)) */
#define kom (0.21435) /* && (1/ms) = 1.5*0.1429 */

/* ICaL */
#define scale (1.2) /* && */
#define p_CaL (3.08027691378999990e-01) /* && (1/(mM^2 ms)) */
#define p_CaL_shannonCa (0.00054)  /* && */
#define p_CaL_shannonK (2.7e-07) /* && */
#define p_CaL_shannonNa (1.5e-08) /* && */

/* everything else */
/* typedef real GlobalData_t; */
#define tau_fCa (2.0)
#define Na_frac ((NatoK_ratio/(NatoK_ratio+1.0)))
