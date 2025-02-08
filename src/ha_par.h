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


_(Nab,		130.0*mM)	/* ??? sodium	*/
_(Kb,		5.4*mM)		/* ??? potassium	*/
_(Cab,		1.8*mM)		/* ??? calcium	*/
_(Mgi,		2.5*mM)		/* ? internal magnium */
_(ECaapp,	60.0*mV)	/* ? some voltage */
_(kCa,		WHICH(0.25,0.025)*mM)	/* ??? */
_(R,		8314*mJ/(mol*Kel)) /* universal gas constant */
_(T,		306.15*Kel)	/* temperature, 33C */
_(F,		96487*(mJ/mV/mol))/* Faraday's constant: Coulomb=J/V=mJ/mV != nF*mV !!! - but used correctly */
_(Cm,		0.05*nF)	/* total cell membrane capacitance */
_(Voli,		0.005884*nL)	/* ? cell internal volume */
_(Volc,		WHICH(0.138,0.136)*Voli)	/* ? cleft volume per cell */
_(Vold,		0.02*Voli)	/* ??? volume per cell */
_(Volrel,	4.41e-5*nL)
_(Volup,	3.969e-4*nL)
_(tNa,		14.3*s)
_(tK,		10.0*s)
_(tCa,		24.7*s)
_(tdi,		0.010*s)
_(_INaK,	70.8253*pA)
_(kNaKK,	1.0*mM)
_(kNaKNa,	11.0*mM)
_(_ICaP,	4.0*pA)
_(kCaP,		0.0002*mM)
_(kNaCa,	0.0374842*pA/(mM*mM*mM*mM))
_(gam,		0.45)
_(dNaCa,	0.0003*1/(mM*mM*mM*mM))
_(PhiNaen,	1.68*pA)	/* NB: sign inverted! - el-neutr Na+ influx for ionic homeostasis */
_(_Iup,		2800.0*pA)
_(kcyca,	0.0003*mM)
_(ksrca,	0.5*mM)
_(kxcs,		0.4)
_(ttr,		0.01*s)
_(alprel,	2e5*pA/mM)
_(kreli,	0.0003*mM)
_(kreld,	0.003*mM)
_(rrecov,	0.815*1/s)
_(PNa,		0.0016*nL/s)
_(_gCaL,	6.75*nS)
_(_gt,		7.5*nS)
_(_gsus,	2.75*nS)
_(_gKs,		1.0*nS)
_(_gKr,		0.5*nS)
_(_gK1,		3.0*nS)
_(_gBNa,	0.060599*nS)
_(_gBCa,	0.078681*nS)
_(Iext,		0.0*pA)
_(ht,		RNONE)
_(RTF,		WHICH(((R*T)/F),T*0.08554))
_(FRT,		(1./RTF))
