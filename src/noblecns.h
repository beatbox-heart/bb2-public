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


__AA(double,Scale)   /* general scaling parameter - USED FOR INPUT ONLY */
__AA(double,TimeScale) /* scaling factor for time, 1 = secs, 1000 = msec */
__AA(double,Kb)   /* [K] in bulk extracellular space */
__AA(double,V)   /* volume of cell or tissue */
__AA(double,VF)   /* volume x Faraday */
__AA(double,PCa3)   /* permeability constant for very slow channel */
__AA(double,PF)  /* rate constant for diffusion from cleft to bulk space */
__AA(double,Pump)   /* maximum pump current */
__AA(double,gNa)   /* sodium channel conductance */
__AA(double,gCa)   /* not used in present version */
__AA(double,gbK)   /* background potassium conductance */
__AA(double,gbNa)   /* background sodium conductance */
__AA(double,gB)   /* read by datain but not used ! - now used for fraction of persistently open h subunits in long-QT syndrome - vnb 951216 */
__AA(double,C)   /* capacitance */
__AA(double,gfNa)   /* sodium component of i(f) conductance */
__AA(double,gfK)   /* potassium component of i(f) conductance */
__AA(double,CAP)   /* reciprocal of C */
__AA(double,dVtest)   /* rate of voltage change for DT change in HEART */
__AA(double,KmCa)   /* binding constant for calcium SR release */
__AA(double,SRvol)   /* volume of uptake compartment of SR */
__AA(double,Relvol)   /* volume of release compartment of SR */
__AA(double,dx)   /* radial distance interval in cylindrical preps */
__AA(double,Km)    /* binding constant for potassium activation of Na pump */
__AA(double,gK1)   /* maximum conductance of inward rectifier */
__AA(double,Nai)   /* internal sodium concentration */
__AA(double,Vi)   /* intracellular volume */
__AA(double,ViF)   /* intracellular volume x Faraday */
__AA(double,KmNa)/* binding constant for sodium activation of sodium pump */
__AA(double,Cai)  /* intracellular calcium concentration   STILL USED ?? */
__AA(double,KmK1) /* Km for K activation of iK1 */
__AA(double,Kmf)   /* Km for K activation of i(f) */
__AA(double,Diff)   /* extracellular potassium diffusion constant */
__AA(double,Tort)   /* tortuosity factor for extracellular diffusion */
__AA(double,kNaCa)   /* scaling factor for INaCa */
__AA(double,Nao)   /* external sodium concentration */
__AA(double,gTO)   /* transient outward current conductance */
__AA(double,CaLim)   /* no longer used */
__AA(double,Ki)   /* intracellular potassium concentration */
__AA(double,ShiftTO)   /* voltage shift for i(to) */
__AA(double,ShiftK1)   /* voltage shift for iK1 */
__AA(double,PCa)   /* calcium channel permeability */
__AA(double,gbCa)   /* background calcium current conductance */
/* __AA(double,Em)   / * membrane potential */
__AA(double,gNaK)   /* non-specific cation channel conductance */
__AA(double,DNaCa)   /* denominator factor for sodium-calcium exchange */
__AA(double,nNaK)   /* stoichiometry for sodium-potassium pump */
__AA(double,Temp)   /* temperature */
__AA(double,RTonF)   /* RT/F */
__AA(double,Ca12m) /* maximum concentration of calcium in SR uptake store */
__AA(double,NaPE)   /* sodium concentration in perfusion electrode */
__AA(double,CaPE)   /* calcium concentration in perfusion electrode */
__AA(double,KPE)   /* potassium concentration in perfusion electrode */
__AA(double,TNaP)   /* time constant for sodium perfusion */
__AA(double,TCaP)   /* time constant for calcium perfusion */
__AA(double,TKP)   /* time constant for potassium perfusion */
__AA(double,T0sav)   /* saved value of t0 */
__AA(double,V12)   /* fraction occupied by uptake SR */
__AA(double,V12F)   /* x Faraday */
__AA(double,V13)   /* fraction occupied by release SR */
__AA(double,V13F)   /* x Faraday */
__AA(double,Tau12)   /* SR uptake constant */
__AA(double,Tau13)   /* SR release constant */
__AA(double,DiffCa)   /* diffusion constant for extracellular calcium */
__AA(double,CaB)   /* bulk [Ca]o, when cleft [Ca] is variable */
__AA(double,ACh)   /* acetyl choline concentration */
__AA(double,AChgK)   /* conductance for iACh -- new formulation */
__AA(double,AChRecKm) /* Km for receptor for ACh activation of K channel */
__AA(double,AChRecfKm)   /* Km for receptor for ACh shift of i(f) */
__AA(double,AChRecCKm)   /* Km for receptor for ACh inhibition of iCa */
__AA(double,iNaCam)   /* maximum rate of Na-Ca exchange */
__AA(double,Kminact) /* Km for Ca inactivation of iCa --- DF formulation */
__AA(double,KmTO)   /* Km for K activation of i(TO) */
__AA(double,MiB)   /* minimal fraction of ibNa in Ca activation mode */
__AA(double,MiK)   /* minimal fraction of iK   in Ca activation mode */
__AA(double,MiK1)   /* minimal fraction of iK1  in Ca activation mode */
__AA(double,MiTO)   /* minimal fraction of iTO  in Ca activation mode */
__AA(double,TauF2) /* limiting recovery time constant for iCa when [Ca]i->0 */
__AA(double,TauRel)/* SR release constant */
__AA(double,vol)   /* fractional extracellular volume */
__AA(double,yNaCa) /* energy barrier position for Na-Ca exchange */
__AA(double,Calmod)   /* calmodulin concentration */
__AA(double,CTrop)   /* c-troponin concentration */
__AA(double,MTrop)   /* m-troponin concentration */
__AA(double,CaVOL)   /* fractional cytosol space occupied by calcium */
__AA(double,BufVol)   /* fractional cytosol space occupied by buffers */
__AA(double,Mg)   /* intracellular magnesium concentration */
__AA(double,PrepLength)   /* length of cell or tissue */
__AA(double,VSurfCa)   /* surface potential on calcium channel */
__AA(double,PCaK)   /* permeability ratio K : Ca for Ca channel */
__AA(double,KcMRel)   /* mitochondrial calcium release constant */
__AA(double,KnMRel)   /* mitochondrial calcium uptake constant */
__AA(double,KcMup)   /* constant for binding to mitochondrial Ca pump */
__AA(double,NcMup)   /* Ca ions binding to mitochondrial pump */
__AA(double,NNMRel)   /* ions bound by uptake pump */
__AA(double,Vmit)   /* volume fraction occupied by mitochondria */
__AA(double,PCa2)   /* permeability of T type calcium channel */
__AA(double,gkATPm) /* maximum conductance for ATP dependent K channels */
__AA(double,kATP)   /* binding constant for ATP control of channels */
__AA(double,KCaCHoff)   /* Ca induced iCa inactivation constant */
__AA(double,KCyCa)   /* see Hilgemann SR pump equations */
__AA(double,Kxcs)   /*                 |               */
__AA(double,KSRCa)   /*_________________|_______________*/
__AA(double,SRleak)   /* calcium leak from SR */
__AA(double,KSLpump)   /* Sarcolemmal Ca pump scaling factor */
__AA(double,KMSLpump)   /* Ca binding constant for sarcolemmal Ca pump */
__AA(double,PNaK)   /* K permeability of sodium channels */
__AA(double,gK)   /* delayed rectifier potassium conductance */
__AA(double,KMNaCa) /* binding constant for Ca activation of Na Ca exchange */
__AA(double,SPvol)   /* fractional volume of subsarcolemmal space */
__AA(double,Vspace)   /* subsarcolemmal space volume */
__AA(double,iCafract) /* fraction of iCa flowing to subsarcolemmal space */
__AA(double,iNCfract) /* fraction of iNaCa flowing to subsarcolemmal space */
__AA(double,SRfract)  /* fraction of SR adjacent to subsarcolemmal space */
__AA(double,TauSPvol) /* time constant for diffusion from subsarcolemmal space */
__AA(double,FuraS)   /* calcium indicator in subsarcolemmal space */
__AA(double,VSpaceF)   /* space volume x Faraday */
__AA(double,Fura)   /* intracellular Ca indicator concentration */
__AA(double,KmCa2)   /* rate constant for calcium release */
__AA(double,CBden)   /* cross bridge density */
__AA(double,CBturn)   /* cross bridge turnover rate */
__AA(double,SarcoLength)   /* sarcomere length */
__AA(double,steepK1) /* steepness of iK1 rectification */
__AA(double,iKm)   /* maximum delayed K current */
__AA(double,iAChm)   /* maximum value of iACh  -- old formulation */
__AA(double,Stim)   /* amplitude of stimulus */
__AA(double,iPulse)   /* current pulse size */
__AA(double,AChgKmax) /* maximum conductance for iACh -- new formulation */

/* Kept here as constant parameters */
__AA(Arneq,Shift)   /* voltage shifts for y processes */
__AA(Arneq,Speed)   /* scaling factor for speed of y processes */
__AA(Ardsub,KC)   /* cleft values of [K] */
__AA(Ardsub,KCE)   /* extrapolated values of cleft [K] */
/* arrays for parameters changed at times T2 .. T9 */
__AA(Ar10,PumpCH)/* values of Na pump activity at change times,pump2 ... pump9 */
__AA(Ar10,KCH)   /* values of [K]o at change times,K2 .... K9 */
__AA(Ar10,NaCH)  /* values of [Na]o at change times,Na2 ... Na9 */
__AA(Ar10,CaCH)  /* values of [Ca]o at change times,Ca2 ... Ca9 */
__AA(Ar10,repCH)/* values of repetition intervals at change times,rep2 ... rep9 */
__AA(Ar10,dtCH) /* values of DT at change times,DT2 ... DT9 */
__AA(Ar10,CACT) /* parameters for Ca activation of ionic currents */

__AA(boolean,valid)   /* used to check validity on exit from various procedures */
__AA(boolean,spmimic) /* mimics subsarcolemmal space by speeding Ca release */
__AA(boolean,BufFast)  /* if true,buffer rates computed,otherwise steady state used */
__AA(boolean,ContractMode)   /* contraction computed if true */
__AA(boolean,SLPumpMode)   /* sarcolemmal calcium pump included if true */
__AA(boolean,Head)   /* used in TABLE to write column headings in non-graphic mode */
__AA(boolean,Fast)   /* integrate sodium activation equations if true */
__AA(boolean,Slow)    /* used in slow modes to set some parameters to steady state values */
__AA(boolean,GNaSS)   /* gNa set to steady state if true */
__AA(boolean,CaBuff)   /* `clamps' [Ca]i in DN equations */
__AA(boolean,ComputeCaIntegral)  /* determines whether this integral is computed */
__AA(boolean,CaSS)   /* if true,Ca set to steady state */
__AA(boolean,ComputeATP) /* determines whether [ATP] is computed */
__AA(boolean,hiAcc)   /* DTSS even further reduced if true */
__AA(boolean,ramp)   /* if true,uses T2 .. T9 for ramps,not steps */

__AA(short,out)   /* sets output parameters in TABLE */
__AA(short,Amode)   /* sets choice for Ca activation of currents */
__AA(short,Bmode)   /* sets choice of calcium buffer equations */
__AA(short,CaOmode)   /* sets choice of equations for [Ca]o */
__AA(short,Camode)   /* sets choice of equations for iCa */
__AA(short,CaNmode)   /* determines whether Na permeates Ca channel */
__AA(short,Iscal)   /* global current scaling parameter */
__AA(short,SRmode)   /* sets choice of equations for SR uptake and release */
__AA(short,Kmode)   /* sets choice of equations for iK */
__AA(short,Namode)   /* sets choice of equations for iNa */
__AA(short,NNaCa)   /* stoichiometry for sodium-calcium exchange */
__AA(short,Nmode)   /* sets choice of Na-Ca equations */
__AA(short,Nrel)   /* no of Ca ions bound to release site */
__AA(short,Prep)   /* sets model to be used  --- not fully up to date ! */
__AA(short,Rmode)   /* sets choice of repriming equations */
__AA(short,Space)   /* sets configuration of extracellular space */
__AA(short,TOmode)   /* sets choice of iTO equations */
__AA(short,Ymode)   /* sets choice of equations for i(f) */
__AA(short,Run)   /* increments following each run of calculation */
/*__AA(short,TIM[26])*/
__AA(short,Iflag)   /* indicates state of numerical computation */
__AA(short,NEQND)/* forces IFLAG = 0 when m equations first integrated */
__AA(short,Chan)   /* subscript for change variables at times T2 .. T9 */
__AA(short,Mmode)  /* sets choice of equations for mitochondrial calcium */
__AA(short,KPNaCa)  /* no of Ca ions bound to control site of Na-Ca exchange*/
__AA(short,Ncalmod)   /* no of Ca ions binding to calmodulin */
__AA(short,Nctrop)   /* no of Ca ions bidning to c troponin */
__AA(short,Nmtrop)   /* no of Ca ions binding to m troponin */

#define Infinite        100000L   /* Infinite is beyond likely Tend */
#define Faraday         965.0e2   /* Faraday constant */
#define delta3          0.0001   /* small value */
#define delta4          0.00001   /* smaller value */
#define delta6          0.0000001   /* very small value */

#define Version         4.0   /* Program version */

#define Update          "1st March 1993"


