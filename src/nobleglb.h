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


/* REALLY GLOBAL VARIABLES - now made contants		*/
/* Mostly reserved for compatibility 			*/
__AA(boolean,Finish)   /* if true,computation terminates */
__AA(short,NEQNN)   /* actual no of differential equations used in Adams */

__AA(double,T0)   /* time at beginning of integration step */
__AA(double,Tout)   /* time at end of integration step */
__AA(double,Tend)   /* time at end of run */
__AA(double,dt)   /* standard step length */
__AA(double,dts)   /* short step length */
__AA(double,dtl)   /* long step length */
__AA(double,TAB)   /* time for next tabulation */
__AA(double,TABT)   /* time interval between tabulations */
__AA(double,dtss)   /* very short integration step interval */
__AA(double,Tstart)   /* time at start of calculation */
/* __AA(double,dtSav)   / * saved value of DT */
__AA(double,on)  /* time for switching on current and voltage pulses */
__AA(double,off)  /* time for switching off current and voltage pulses */
__AA(double,Rep)   /* repetition interval between pulses */
__AA(double,TP)   /* time for activation of stimulus */
__AA(double,TPoff)   /* time for switch off of stimulus */
__AA(double,Restart)   /* time for creating a restart file */
__AA(double,StimSize)   /* stimulus amplitude */
__AA(double,PulseLength)   /* duration of pulse */
__AA(double,iPulseSize)   /* current pulse amplitude */
/* __AA(double,RestSav)   / * stored value of RESTART */
/*__AA(double,onSav)   / * stored value of ON */
/*__AA(double,offSav)   / * stored value of OFF */
/*__AA(double,TRsav)   / * stored value for time of restart */
__AA(double,TPS)   /* time of application of stimulus */
__AA(double,TPsav) /* saved value of TPS */
__AA(double,EC) /* clamp potential/ now used as initial potential only - VNB */
