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


/* SHOULD BE MADE LOCAL VARIABLES */

__AA(double,iCa3)   /* very slow calcium current */
/* __AA(double,iNa)   / * sodium channel current */
#define iNa (y[56])
/* __AA(double,iCa)   / * calcium channel net current */
#define iCa (y[59])
/* __AA(double,iK)   / * delayed rectifer K current */
#define iK (y[57])
__AA(double,ifNa)   /* sodium component of i(f) */
__AA(double,iTO)   /* transient outward current */
__AA(double,ibK)   /* background potassium current */
__AA(double,ibNa)   /* background sodium current */
__AA(double,ifK)   /* potassium component of i(f) */
__AA(double,iP)   /* sodium pump current */
/* __AA(double,iK1)   / * inward rectifier current */
#define iK1 (y[58])
__AA(double,itot)   /* total ionic current */
__AA(double,imK)   /* net membrane potassium current */
__AA(double,imNa)   /* net membrane sodium current */
__AA(double,iNaCa)   /* sodium-calcium exchange current */
__AA(double,ImCa)   /* net membrane calcium current */
__AA(double,iNaK)   /* non-specific cation channel current */
__AA(double,iSI)   /* net `second inward' current */
__AA(double,iCaNa)   /* sodium flow through calcium channel */
__AA(double,iCaK)   /* potassium flow through calcium channel */
__AA(double,iCaCa)   /* calcium flow through calcium channel */
__AA(double,ibCa)   /* background calcium current */
__AA(double,iUP)   /* calcium uptake to SR,expressed as a current */
__AA(double,iTran)   /* flow from uptake to release SR */
__AA(double,iRel)   /* calcium release,expressed as a current */
__AA(double,ikATP)   /* ATP dependent K channel current */
__AA(double,FCtrop)   /* occupied c troponin */
__AA(double,FMtrop)   /* occupied m troponin */
__AA(double,IrelS)   /* SR release within space */
__AA(double,IupS)   /* SR uptake within space */
__AA(double,ITranS)   /* SR transfer within space */
__AA(double,ImCaS)   /* net membrane calcium current in space */
__AA(double,INaCaS)   /* exchange current in space */
__AA(double,ICaCaS)   /* calcium component of iCa in space */
__AA(double,ICaNaS)   /* sodium component of iCa in space */
__AA(double,ICaKS)   /* potassium component of iCa in space */
__AA(double,iCas)   /* calcium current in space */
__AA(double,ENaCa)   /* equilibrium potential for INaCa */
__AA(double,ECa)   /* calcium equilibrium potential */
__AA(double,CaCHoff)   /* fraction of inactivation sites with Ca bound */
__AA(double,CaCHon)   /* fraction of calcium channels available */
__AA(double,Fcyt) /* see Hilgemann SR pump equations */
__AA(double,FSR)   /*                 |               */
/*__AA(double,FcytFree)  / *                 |               */
/*__AA(double,FSRFree)   / *_________________|_______________*/
__AA(double,ENa)   /* sodium equilibrium potential */
__AA(double,Cao)   /* external calcium concentration */
__AA(double,Emh)   /* sodium channel reversal potential */
__AA(double,EK)   /* potassium equilibrium potential */
__AA(double,CBdenA)   /* available cross bridge density */
__AA(double,CaFmit) /* rate of flow of Ca ions to and from mitochondria */
/*__AA(double,ECaS)   / * calcium equilibrium potential in space */
/*__AA(double,ENaCas)   / * NaCa exchange equilibrium potential in space */
/*__AA(double,IsoTen)   / * isometric tension */
__AA(double,K)   /* external potassium concentration */
__AA(double,iACh)   /* ACh activated K current */
__AA(double,iCa2)   /* T type calcium current */
#ifdef HEART
__AA(double,save1)   /* stored values between calls of DESOL */
__AA(double,save2)   /*                  |                   */
__AA(double,save3)   /*                  |                   */
__AA(double,save4)   /*                  |                   */
__AA(double,save5)   /*                  |                   */
__AA(double,save6)   /*                  |                   */
__AA(double,save7)   /*                  |                   */
__AA(double,save8)   /*                  |                   */
__AA(double,save9)   /*                  |                   */
__AA(double,save10)   /*                  |                   */
__AA(double,save11)   /*                  |                   */
__AA(double,save12)   /*__________________|___________________*/
__AA(double,Esav)   /* stored value of EC */
__AA(double,DCNP)   /* Sarcolemmal Calcium pump variable */
__AA(double,DCP)   /*_________________|_________________*/
#endif

/*__AA(Arneq,saven0)  / * save arrays used by DESOL in calls to ADAMS */
/*__AA(Arneq,saven1)  / *                       |                     */
/*__AA(Arneq,saven2)  / *                       |                     */
/*__AA(Arneq,saven3)  / *_______________________|_____________________*/
/*__AA(Arneq,savY)    / * stored values of y -- used dy DeSol */
/*__AA(Arneq,Ysave)   / * stored values of y -- used by main program */
/*__AA(Arneq,Yres)    / * restart values of y array */

/* arrays for values of parameters in extracellular cleft spaces */
#if 0
/* !!!  should move all this to local place */
__AA(Ardsub,EKC)   /* cleft values of EK */
/* __AA(Ardsub,IMKSAV)   / * saved values of imK */
__AA(Ardsub,IMKC)   /* cleft values of imK */
__AA(Ardsub,KC)   /* cleft values of [K] */
__AA(Ardsub,IKC)   /* cleft values of iK */
__AA(Ardsub,IK1C)   /* cleft values of iK1 */
__AA(Ardsub,IPC)   /* cleft values of ip */
__AA(Ardsub,KCE)   /* extrapolated values of cleft [K] */
__AA(Ardsub,KCSAV)   /* saved values of [K] cleft */
__AA(Ardsub,ITOC)   /* cleft values of iTO */
__AA(Ardsub,IFKC)   /* cleft values of K component of i(f) */
__AA(Ardsub,ICAKC)   /* cleft values of K component of iCa */
__AA(Ardsub,YSAV,IFNAC)   /* cleft values of Na component of i(f) */
__AA(Ardsub,IBKC)   /* cleft values of background K current */
__AA(Ardsub,IACHC)  /* cleft values of iACh */
#endif




/* Global INTEGER variables */

/*__AA(short,Counter)   / * counts number of integration steps */
/* __AA(short,Kflag)    / * used by desol for return from ADAMS and for interpolation */
/*__AA(short,Jflag)   / * used by ADAMS */
/*__AA(short,M)   / * generally available global integers */
/*__AA(short,J)   / *                    |                */
/* __AA(short,L)   / *____________________|________________*/
/*__AA(short,NREDO)/ * tests for discontinuity in ADAMS integration procedure */
/* __AA(short,Isave1)   / * parameters used by DESOL */
/* __AA(short,Isave2)   / *            |             */
/* __AA(short,Isave3)   / *            |             */
/* __AA(short,Isave4)   / *____________|_____________*/

/* Global BOOLEAN variables */
/* __AA(boolean,Done)   / * if true,event has occurred */
