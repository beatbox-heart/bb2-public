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


/*_(ATP,		5.0)	// ATP concentration mM */
_(Cao,		2.0)	/* external Ca concentration */
_(C,		0.0002)	/* capacitance */
_(CTrop,	0.05)	/* c-troponin concentration */
_(DNaCa,	0.0)	/* denominator factor for sodium-calcium exchange */
_(gB,		0.0)	/* read by datain but not used ! - now used for fraction of persistently open h subunits in long-QT syndrome - vnb 951216 */
_(gbCa,		0.00025)	/* background calcium current conductance */
_(gbK,		0.0006)	/* background potassium conductance */
_(gK1,		1.0)	/* maximum conductance of inward rectifier */
_(gbNa,		0.0006)	/* background sodium conductance */
_(gNa,		2.5)	/* sodium channel conductance */
_(Gto,		0.005)	/* transient outward current conductance */
_(iKm,		1.0)	/* maximum delayed K current */
_(iNaCam,	1.e4)   /* maximum rate of Na-Ca exchange */
_(kATP,		0.1)	/* binding constant for ATP control of channels */
_(Ko,		4.0)	/* [K] in bulk extracellular space */
_(KCyCa,	0.0003)	/* see Hilgemann SR pump equations */
_(Km,		1.0)	/* binding constant for potassium activation of Na pump */
_(KmCa,		0.0005)	/* binding constant for calcium SR release */
_(KmCa2,	250.0)	/* rate constant for calcium release */
/*_(Kmf,		45.0)	// Km for K activation of i(f) */
_(KmK1,		10.0)	/* Km for K activation of iK1 */
_(KmNa,		40.0)	/* binding constant for sodium activation of sodium pump */
_(kNaCa,	0.0005)	/* scaling factor for INaCa */
_(KSRCa,	0.5)	/* see Hilgemann SR pump equations */
_(Kxcs,		0.4)	/* see Hilgemann SR pump equations */
_(MTrop,	0.02)	/* m-troponin concentration */
_(Nao,		140.0)	/* external sodium concentration */
_(nNaK,		1.5)	/* stoichiometry for sodium-potassium pump */
_(PCa,		0.25)	/* calcium channel permeability */
_(PCaK,		0.002)	/* permeability ratio K : Ca for Ca channel */
_(PNaK,		0.12)	/* K permeability of sodium channels */
_(Pump,		0.7)	/* maximum pump current */
_(ShiftK1,	20.0)	/* voltage shift for iK1 */
_(steepK1,	2.0)	/* steepness of iK1 rectification */
_(Temp,		37.0)	/* temperature Celcius */
_(V12,		0.01)	/* fraction occupied by uptake SR */
_(V13,		0.1)	/* fraction occupied by release SR */
_(VS,		50.0)	/* surface potential on calcium channel */
_(yNaCa,	0.5)	/* energy barrier position for Na-Ca exchange */
_(CellLength,	0.08)	/* length of cell, mm */
_(CellRadius,	0.015)	/* radius of cell, mm */
_(vol,		0.4)	/* fractional extracellular volume */

_(KKK,KCyCa*Kxcs/KSRCa)
_(KK_K,KCyCa*Kxcs+KCyCa)
_(Nao3,Nao*Nao*Nao)
_(RT_F,(Temp+273)*0.08554)
_(SRvol,M_PI*CellRadius*CellRadius*CellLength*V12)   /* volume of uptake compartment of SR */
_(VSFRT,VS/RT_F)
_(Vi,M_PI*CellRadius*CellRadius*CellLength*(1-vol-V12-V13)) /* intracellular volume */
_(ViF,Vi*Faraday)

_(krel,SRvol*V13/(Vi*V12))
_(kup,Vi/SRvol)
_(ktr,V12/V13)

_(Iext,0.0)		/* external current, nA/cell, negative as all trmbr currents */
_(IV,0.0)		/* external current, mV/ms, positive as produced by Laplacian */
