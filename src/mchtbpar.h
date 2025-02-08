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

/* Meta-parameters */
_(TB,	1,	"contribution of TB model",			    1)
_(CH,	1,	"contribution of CH model",			    1)

/* TB model constants */
_(Ktb,	1.08e5,	"TB phytoplankton saturation constant",		    ug N/m^3)
_(r,	0.3,	"TB phytoplankton maximal growth rate",		    1/d)
_(Rm,	0.7,	"TB zooplankton maximal grazing rate",		    1/d)
_(alpTB,5.7e3,	"TB zooplankton grazing satiation constant",	    ug N/m^3)    
_(gam,	0.05,	"TB zooplankton grazing efficiency",		    1)
_(mu,	0.012,	"TB zooplankton mortality and predation rate",	    1/d)

/* CH model constants */
_(k,	0.0154, "CH larvae grazing volume vs (weight/ug) coefficient",   m^3/d)
_(nu,	0.2234, "CH larvae grazing volume vs (weight/ug) exponent",	    1)
_(Mmax,	0.089,	"CH larvae initial mortality rate",		    1/d)
_(b,	0.005,	"CH larvae mortality decrease rate",		    1/d)
_(g,	0.0875,	"CH zooplankton instantaneous daily growth rate",   1/d)
_(n,	0.67,	"CH a larval metabolic constant",		    1)
_(j,	0.002,	"CH larvae metabolic change rate",		    1)
_(betmax,0.48,	"CH larvae maximal digestive coefficient (w acct of metabolic loss)",	    1)
/*_(betmin,0.135,	"CH larvae minimal (initial) digestive coefficient (w acct of metabolic loss)",1)*/
_(dbet,0.48-0.135,"CH larvae digestive coefficient increase (w acct of metabolic loss)",1)
_(Gmax,	0.12,	"CH larvae max daily growth rate",		    1)
_(metage, 100,	"CH larvae metamorphosis age",			    d)
_(metwgt, 3165,	"CH larvae metamorphosis weight",		    ug)
_(Wi,	33,	"CH init. larvavweight",			    ug)

#if (! VARKCH)
  _(Kch,	    2.6,	"CH larvae weight-metabolic cost coefficient",	    ug^{1-n})
#endif
_(Wc,	    3.3,"(new) small larvae weight",			    ug)
_(metagecoeff, 1,"(new) larvae metamorphosis coefficient",	    d^{-1})
_(metwgtcoeff, 1,"(new) larvae metamorphosis coefficient",	    d^{-1})
_(metagealw, 3,	"(new) larvae metamorphosis age allowance",	    d)
_(metwgtalw, 31,"(new) larvae metamorphosis weight allowance",	    ug)
_(starvcoeff, 0.001, "(new) starvation death coefficient",		    1)
_(starvexp, 1, "(new) starvation death exponent",		    1)
_(Glexp, 0.0 /*0.25*/,	"(new) larvae max daily growth rate decrease exponent",   1)

_(Hollingtype, 1,	"(new) Holling type of larva satiation kinetics (1/2)",   1)
