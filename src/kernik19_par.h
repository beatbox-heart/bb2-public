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

/* List of parameters values of which can be changed in the script */

/* Extracellular concentrations */
_(Nao		, 140.0)    	/* && (mM) */
_(Ko		, 5.4)      	/* && (mM) */
_(Cao		, 1.8)		/* && (mM) */

/* Other */
_(g_Na		, 9.72061340924100037)    /* && (1/ms) */
_(g_f		, 0.0435)	/* && (1/ms) */
_(g_CaT     ,     0.185) /* && */
_(g_to		, 1.17833333333299997e-01)	/* && (1/ms) */
_(g_Ks		, 0.0077)	/* && (1/ms) */
_(g_Kr		, 0.218025)	/* && (1/ms) */
_(g_K1		, 1.33785777797606004e-01)	/* && (1/ms) */
_(g_PCa		, 0.2625)	/* && (mV/ms) */ 
_(g_b_Na	, 4.35e-4)	/* && (1/ms) */
_(g_b_Ca	, 3.6704e-4)	/* && (1/ms) */

/* INaCa */
_(kNaCa		, 1100.0)	/* && (mV/ms) */
_(KmCa		, 1.38)  	/* && (mM) */
_(KmNai		, 87.5)		/* && (mM) */
_(Ksat		, 0.1)		/* && */
_(gamma		, 0.7)		/* && */
_(alpha     ,2.75)   /* && */


/* INaK */
_(Km_K		, 1.0)		/* && (mM) */
_(Km_Na		, 40.0)		/* && (mM) */
_(PNaK      ,2.476116) /* && (mV/ms) */

/* IpCa */
_(KPCa		, 0.0005)	/* && (mM) */

/* Ca2+ buffering */
_(Buf_C		, 0.06)		/* && (mM) */
_(Buf_SR	, 12.0)		/* && (mM) */
_(Kbuf_C	, 0.0006)	/* && (mM) */
_(Kbuf_SR	, 0.3)		/* && (mM) */

/* Ileak */
_(V_leak	, 1.6e-6)	/* &&  ~ magnitude of j_leak */


/* Irel */
_(ks		, 12.5)		/* && (1/ms) ~ magnitude of j_rel */

/* Iup */
_(Kup		, 1.755e-4)	/* && (mM) */
_(VmaxUp	, 1.105e-4)	/* && (mM/ms) = 0.000425 * 0.26 ~ magnitude of j_up */


