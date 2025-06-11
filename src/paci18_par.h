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

/* Parameters from optimizer */
_(j_up_max	,0.5113e-3)	/* (mM/ms) */
_(j_rel_max     ,62.5434e-3)	/* (1/ms) */
_(RyRa1		, 0.05354)  	/* (1) */
_(RyRa2		, 0.0488)   	/* (1) */
_(RyRahalf	, 0.02427)  	/* (1) */
_(RyRohalf	, 0.01042)  	/* (1) */
_(RyRchalf	, 0.00144)  	/* (1) */
_(kNaCa		, 3917.0463)	/* (mV/ms) */
_(PNaK		, 2.6351)   	/* (mV/ms) */
_(Kup		, 3.1928e-4)	/* (mM) */
_(j_leak_max	, 4.7279e-7)	/* (1/ms) */
_(alpha		, 2.5371)   	/* (1) */

/* Extracellular concentrations */
_(Nao		, 151.0)    	/* (mM) */
_(Ko		, 5.4)      	/* (mM) */
/*_(Cao		, 1.8)      (mM) - making this constant allows tabulation of Xr1 */

/* Intracellular concentrations */
_(Ki		, 150.0)     	/* (mM) */

/* Other */
_(g_Na		, 3.6712302)    /* (1/ms) */
_(g_NaL       	, 17.25e-3)	/* (1/ms) */
_(E_f         	, -17.0)	/* (mV) lit mV */
_(g_f		, 30.10312e-3)	/* (1/ms) lit 1/ms */
_(g_CaL		, 8.635702e-5)	/* (V/s)*(m^3/C) */
_(g_to		, 29.9038e-3)	/* (1/ms) */
_(g_Ks		, 2.041e-3)	/* (1/ms) */
_(g_Kr		, 29.8667e-3)	/* (1/ms) */
_(g_K1		, 28.1492e-3)	/* (1/ms) */
_(g_PCa		, 0.4125)	/* (mV/ms) */
_(g_b_Na	, 0.95e-3)	/* (1/ms) */
_(g_b_Ca	, 0.727272e-3)	/* (1/ms) */

/*% INaCa */
_(KmCa		, 1.38)  	/* (mM) */
_(KmNai		, 87.5)		/* (mM) */
_(Ksat		, 0.1)		/* (1) */
_(gamma		, 0.35)		/* (1) */

/*% INaK */
_(Km_K		, 1.0)		/* (mM) */
_(Km_Na		, 40.0)		/* (mM) */

/*% IpCa */
_(KPCa		, 0.0005)	/* (mM) */

/*% Ca2+ buffering */
_(Buf_C		, 0.25)		/* (mM) */
_(Buf_SR	, 10.0)		/* (mM) */
_(Kbuf_C	, 0.001)	/* (mM) */
_(Kbuf_SR	, 0.3)		/* (mM) */
