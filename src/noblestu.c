/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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

#include <math.h>
#include "beatbox.h"

/* - how to do this on mips?
#if MAXREAL != MAXDOUBLE && defined(NOBLE)
  #error QHEART needs double precisioni !!
#endif
*/

#include "p2c.h"
#include "noble.h"
#include "nobleini.h"

#include "scr.h"


#include "noblest2.h"
#include "noblestu.h"

#define WindowInp       false
#define WinCounter      0
#define path            1
static void ReadKeyFile(void) {
  strcpy(licence, "00000");
  *site = '\0';
  *contact = '\0';
}

static double _Z2, _Z3, _Z4, _Z5, _Z6, _Z7;
static Char _ch;


void SetModels(void) {
  STANDARD();
  fprintf(Filout3, "%s\n\n", File0);
  /* look for starred input to set standard preparation parameter */
  _ch = getc(Filin1);
  if (_ch < '\001') {
    valid = false;
    return;
  }
  putc(_ch, Filout3);
  if (_ch != '*') {
    PURKINJE();
    return;
  }
  /* finds star then looks for colon */
  while (_ch != ':') {
    _ch = getc(Filin1);
    if (_ch < '\001') {
      valid = false;
      return;
    }
    putc(_ch, Filout3);
  }
  /* looks for capitals */
  while (_ch < 'A' || _ch > 'Z') {
    _ch = getc(Filin1);
    if (_ch < '\001') {
      valid = false;
      return;
    }
    putc(_ch, Filout3);
  }
  switch (_ch) {

  case 'P':
    PURKINJE();
    strcpy(model, "      Purkinje fibre");
    break;

  case 'S':
    SINUS();
    strcpy(model, "Sinus node");
    break;

  case 'E':
    SINUS();
    SINSIN();
    strcpy(model, "Single sinus cell");
    break;

  case 'H':
    HATRIUM();
    strcpy(model, "Rabbit atrium");
    break;

  case 'Y':
    HATRIUM();
    SINAT();
    strcpy(model, "Single atrial cell");
    break;

  case 'V':
    HATRIUM();
    SINAT();
    VENTRICLE();
    strcpy(model, "Guinea-pig ventricular cell");
    break;

  case 'G':
    HATRIUM();
    SINAT();
    GPCELL();
    spmimic = true;
    strcpy(model, "Simplified GP cell");
    break;

  case 'R':
    RATVCELL();
    strcpy(model, "Rat ventricular cell");
    break;

  case 'F':
    FROGSINUS();
    strcpy(model, "Frog sinus node cell");
    break;

  case 'A':
    ATRIUM();
    strcpy(model, "Frog atrial cell");
    break;

  default:
    PURKINJE();
    strcpy(model, "      Purkinje fibre");
    break;
  }
}  /* of procedure SetModels */



/*static*/ void ScaleParam(void) {
						#undef _L
  short _L;

  /* NB This procedure needs updating for more recent models */



  gNa = Scale * gNa;
  PCa = Scale * PCa;
  iKm = Scale * iKm;
  gK = Scale * gK;
  PCa3 = Scale * PCa3;
  PCa2 = Scale * PCa2;
  gK1 = Scale * gK1;
  gB = Scale * gB;
  gbNa = Scale * gbNa;
  gbK = Scale * gbK;
  gfK = Scale * gfK;
  gfNa = Scale * gfNa;
  gTO = Scale * gTO;
  C = Scale * C;
  gbCa = Scale * gbCa;
  Pump = Scale * Pump;
  StimSize = Scale * StimSize;

  kNaCa = Scale * kNaCa;
  iPulseSize = Scale * iPulseSize;
  for (_L = 1; _L <= 8; _L++)
    PumpCH[_L] = Scale * PumpCH[_L];

}  /* of procedure ScaleParam */



static void RemoveFText(void) {
  if (!*File1) return;
  if (Filout1 != NULL) {
    fclose(Filout1);
    erase(Filout1);
    Filout1 = NULL;
  }
/* p2c: stunit.pas, line 240:
 * Warning: Symbol 'ERASE' is not defined [221] */
  /*FileOpen(&Filout1, File1, true);*/  frewrite(Filout1,File1);
  if (out == 1)
    return;

  if (Filout2 != NULL) {
    fclose(Filout2);
    erase(Filout2);
    Filout2 = NULL;
  }
/* p2c: stunit.pas, line 244:
 * Warning: Symbol 'ERASE' is not defined [221] */
  /*FileOpen(&Filout2, File2, true);*/	frewrite(Filout2,File2);
}  /* of procedure RemoveFText */


/*static*/ void SetChangeVar(void) {
  short _L;


  for (_L = 0; _L <= 9; _L++) {
    KCH[_L] = Kb;
    PumpCH[_L] = Pump;
    NaCH[_L] = Nao;
    CaCH[_L] = _Y[23]/* =Cao */;
    repCH[_L] = Rep;
    dtCH[_L] = dt;

  }





  dts = 0.5 * dt;
  dtss = 0.1 * dts;
  dtl = dt;

}  /* of procedure SetChangeVar */


/*static*/ void SetUnchanVar(void) {
  short _L;


  /* reset unchanged variables to initial values to  ensure  that  if
     only initial value  is changed  from the  default value,  all  other
     values will be held at the new initial value */

  _Z2 = 0.0;
  _Z3 = 0.0;
  _Z4 = 0.0;
  _Z5 = 0.0;
  _Z6 = 0.0;
  _Z7 = 0.0;
  for (_L = 2; _L <= 9; _L++) {
    if (KCH[_L - 1] - KCH[_L - 2] != 0)
      _Z2 = 1.0;
    if (PumpCH[_L - 1] - PumpCH[_L - 2] != 0)
      _Z3 = 1.0;
    if (NaCH[_L - 1] - NaCH[_L - 2] != 0)
      _Z4 = 1.0;
    if (CaCH[_L - 1] - CaCH[_L - 2] != 0)
      _Z5 = 1.0;
    if (repCH[_L - 1] - repCH[_L - 2] != 0)
      _Z6 = 1.0;
    if (dtCH[_L - 1] - dtCH[_L - 2] != 0)
      _Z7 = 1.0;
  }


  if (_Z2 == 0) {
    for (_L = 1; _L <= 8; _L++)
      KCH[_L] = Kb;
  }
  if (_Z3 == 0) {
    for (_L = 1; _L <= 8; _L++)
      PumpCH[_L] = Pump;
  }
  if (_Z4 == 0) {
    for (_L = 1; _L <= 8; _L++)
      NaCH[_L] = Nao;
  }
  if (_Z5 == 0) {
    for (_L = 1; _L <= 8; _L++)
      CaCH[_L] = _Y[23]/* =Cao */;
  }
  if (_Z6 == 0) {
    for (_L = 1; _L <= 8; _L++)
      repCH[_L] = Rep;
  }
  if (_Z7 == 0) {
    for (_L = 1; _L <= 8; _L++)
      dtCH[_L] = dt;
  }

}  /* of procedure SetUnchanVar */






/*==========================================================================*/

void START(void) {

  /*     This procedure sets initial values of parameters.

         It is the driver procedure for a set of procedures that set and
         store parameters. These are:

         1. The standard set procedure, STANDARD

         2. The individual preparation procedures:
            PURKINJE,SINUS,SINSIN,HATRIUM,SINAT,VENTRICLE,RATVCELL,FROGSINUS

         3. The procedure for reading input data files: DATAIN

         4. The procedures for calculating and storing initial values:
            RECALL,LISTDATA,STORE,INIT

         5. The interactive menu procedure MENU

         6. The screen display procedure SDISPLAY

         7. The file input and output procedures FENTER & FOUT

         8. The procedure for listing parameters: LISTD

         9. The procedure (READKEYFILE) that reads HEART.KEY. This
            procedure is in the unit KEYUNIT.TPU


The procedure contains:

          SetModels, ScaleParam, RemoveFText,
          SetChangeVar, SetUnchanVar, ResetChangeVar,
          CheckSpeed.

===========================================================================*/
  valid = true;

  ReadKeyFile();
  if (Finish)
    return;

  /* no screen HEART display on reruns or command line runs */
  /*2.2*/

  /* reset initial values on rerun unless last run created a restart
     file or if mode = 0, close files */


  /* Open files for input and output */
  FENTER();

  /* Set standard parameters unless this is a rerun */
  SetModels();

  /* Set default values of variable changes and dt parameters */
  SetChangeVar();

  /* call DATAIN or MENU to read input data */
  DATAIN();





  if (out != 0)
    RemoveFText();

  SetUnchanVar();

  /* list input parameters in output file */
  LISTDATA();



  /* scale parameters if required */
  if (Scale != 1)
    ScaleParam();

  /* set initial values of various parameters */
  INI();

}  /* of procedure Start */


/* $I LISTD.PAS*/
/* ListD.pas */

/* lists parameters set in HEART */


/* ************************************************************************* */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*             +    +   +++++      +      ++++    +++++++                    */
/*             +    +   +         + +     +   +      +                       */
/*             ++++++   ++++     +   +    ++++       +                       */
/*             +    +   +       +++++++   +   +      +                       */
/*             +    +   +++++  ++     ++  +    +     +                       */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/* ************************************************************************* */


static void Rwrite(double r) {
  Char s[256];
  Char *TEMP1;


  fgets(s, 256, Filin3);
  TEMP1 = strchr(s, '\n');
  if (TEMP1 != NULL)
    *TEMP1 = 0;
  if (r != 0 && r != Infinite)
    fprintf(Filout3, "%s%12.7f\n", s, r);

}



static void Iwrite(short i) {
  Char s[256];
  Char *TEMP1;


  fgets(s, 256, Filin3);
  TEMP1 = strchr(s, '\n');
  if (TEMP1 != NULL)
    *TEMP1 = 0;
  fprintf(Filout3, "%s%6d\n", s, i);

}



static void Bwrite(boolean b) {
  Char s[256];
  Char *TEMP1;


  fgets(s, 256, Filin3);
  TEMP1 = strchr(s, '\n');
  if (TEMP1 != NULL)
    *TEMP1 = 0;
  fprintf(Filout3, "%s%s\n", s, b ? " TRUE" : "FALSE");

}



/*==========================================================================*/

void LISTDATA(void) {  /*of ListData*/

  /* This procedure lists the starting data in the third output file
     It uses text from LIST.TXT

     Calls: no other procedures

     Contains: Rwrite, Iwrite, Bwrite.

===========================================================================*/
  short l;





  strcpy(File5, "list.txt");
  /*FileOpen(&Filin3, File5, false);*/ freset(Filin3,File5);
  fprintf(Filout3, "%s\n", site);
  fprintf(Filout3, "Contact: %s\n", contact);
  fprintf(Filout3, "Program version: %3.1f\n\n\n", Version);
  Rwrite(ACh);
  Rwrite(AChgKmax);
  Rwrite(alpha[11]);
  Iwrite(Amode);   /*amp*/
  Rwrite(0.0);
  Rwrite(_Y[34]);
  Iwrite(Bmode);
  Bwrite(BufFast);
  Rwrite(C);
  Rwrite(Calmod);
  Rwrite(CTrop);
  Rwrite(_Y[3]);
  Rwrite(_Y[23]/* =Cao */);
  Iwrite(CaOmode);
  Rwrite(CaLim);
  Iwrite(Camode);
  Iwrite(CaNmode);
  if (Amode != 0) {   /*continue*/
    for (l = 0; l <= 3; l++)
      Rwrite(CACT[l]);
  } else {
    for (l = 1; l <= 4; l++) {
      fscanf(Filin3, "%*[^\n]");
      (void)getc(Filin3);
    }
  }
  Iwrite(0);
  Bwrite(ContractMode);
  Rwrite(_Y[12]);
  Rwrite(_Y[11]);
  Rwrite(Ca12m);
  for (l = 1; l <= 8; l++) {
    if (CaCH[l] != _Y[23]/* =Cao */)
      Rwrite(CaCH[l]);
    else
      Rwrite(0.0);
  }
  Rwrite(_Y[6]);
  Rwrite(_Y[15]);
  Rwrite(dx);
  Rwrite(DiffCa);   /*disp*/
  Iwrite(0);
  Rwrite(DNaCa);
  Rwrite(dVtest);
  Rwrite(dt);
  Rwrite(EC);
  for (l = 2; l <= 9; l++)
    Rwrite(0.0);
  Rwrite(_Y[7]);
  Rwrite(_Y[14]);
  if (PCa3 != 0)
    Rwrite(_Y[16]);
  else
    Rwrite(0.0);
  if (Fura != 0) {
    Rwrite(Fura);
    Rwrite(_Y[33]);
  } else {
    fscanf(Filin3, "%*[^\n]");
    (void)getc(Filin3);
    fscanf(Filin3, "%*[^\n]");
    (void)getc(Filin3);
  }
  Rwrite(gB);
  Rwrite(gbK);
  Rwrite(gbNa);
  Rwrite(gbCa);
  Rwrite(gfK);
  Rwrite(gfNa);
  Rwrite(gK1);
  Rwrite(gK);
  Rwrite(gTO);
  Rwrite(gNa);
  Rwrite(_Y[8]);
  Rwrite(iAChm);   /*idisp*/
  Iwrite(0);
  Rwrite(iKm);
  Rwrite(iPulseSize);
  Rwrite(iNaCam);
  Iwrite(Iscal);
  Rwrite(Kb);
  Rwrite(_Y[10]);
  Rwrite(_Y[2]/*=K which is now local variable VNB*/);
  if (SRmode == 2) {
    Rwrite(KCyCa);
    Rwrite(KSRCa);
    Rwrite(Kxcs);
  } else {
    for (l = 1; l <= 3; l++) {
      fscanf(Filin3, "%*[^\n]");
      (void)getc(Filin3);
    }
  }
  if (SLPumpMode) {
    Rwrite(KMSLpump);
    Rwrite(KSLpump);
  } else {
    for (l = 1; l <= 2; l++) {
      fscanf(Filin3, "%*[^\n]");
      (void)getc(Filin3);
    }
  }
  if (Camode == 8)
    Rwrite(KCaCHoff);
  else
    Rwrite(0.0);
  if (Bmode != 0) {
    for (l = 19; l <= 22; l++) {
      Rwrite(alpha[l]);
      Rwrite(beta[l]);
    }
  } else {
    for (l = 1; l <= 8; l++) {
      fscanf(Filin3, "%*[^\n]");
      (void)getc(Filin3);
    }
  }
  Rwrite(kNaCa);
  Rwrite(KMNaCa);
  Rwrite(KmCa);
  Rwrite(KmCa2);
  Rwrite(Kmf);
  Rwrite(KmK1);
  Rwrite(KmTO);
  Rwrite(Km);
  Iwrite(Kmode);
  for (l = 1; l <= 8; l++) {
    if (KCH[l] != Kb)
      Rwrite(KCH[l]);
    else
      Rwrite(0.0);
  }
  Rwrite(PrepLength);
  Rwrite(_Y[9]);   /*mode*/
  Iwrite(1);
  Rwrite(MTrop);
  Rwrite(MiB);
  Rwrite(MiK);
  Rwrite(MiK1);
  Rwrite(MiTO);
  Rwrite(_Y[1]);
  Iwrite(Namode);
  Rwrite(Nao);
  Iwrite(NNaCa);
  Iwrite(Nmode);
  Rwrite(nNaK);
  Iwrite(Nrel);
  for (l = 1; l <= 8; l++) {
    if (NaCH[l] != Nao)
      Rwrite(NaCH[l]);
    else
      Rwrite(0.0);
  }
  Rwrite(off);
  Rwrite(on);
  Iwrite(out);
  Rwrite(PF);   /*pmode*/
  Iwrite(0);
  Rwrite(PCa);
  Rwrite(PCa2);
  Rwrite(PCa3);
  Rwrite(PNaK);
  Rwrite(_Y[13]);
  Iwrite(Prep);   /*pulsesize*/
  Rwrite(0.0);
  Rwrite(Pump);
  for (l = 1; l <= 8; l++) {
    if (PumpCH[l] != Pump)
      Rwrite(PumpCH[l]);
    else
      Rwrite(0.0);
  }
  Rwrite(_Y[18]);
  Rwrite(_Y[17]);
  Rwrite(Rep);
  for (l = 1; l <= 8; l++) {
    if (repCH[l] != Rep)
      Rwrite(repCH[l]);
    else
      Rwrite(0.0);
  }
  Rwrite(Restart);
  Iwrite(Rmode);   /*re*/
  Rwrite(0.0);   /*rs*/
  Rwrite(0.0);
  Rwrite(Scale);
  for (l = 4; l <= 9; l++)
    Rwrite(Shift[l]);
  Rwrite(ShiftTO);
  Rwrite(ShiftK1);
  for (l = 13; l <= 18; l++)
    Rwrite(Shift[l]);
  for (l = 4; l <= 9; l++)
    Rwrite(Speed[l]);
  for (l = 13; l <= 18; l++)
    Rwrite(Speed[l]);
  Iwrite(Space);
  if (SPvol != 0) {
    Rwrite(SPvol);
    Rwrite(iCafract);
    Rwrite(iNCfract);
    Rwrite(SRfract);
    Rwrite(TauSPvol);
  } else {
    for (l = 1; l <= 5; l++) {
      fscanf(Filin3, "%*[^\n]");
      (void)getc(Filin3);
    }
  }
  Iwrite(SRmode);
  Rwrite(SRleak);
  Rwrite(steepK1);
  Rwrite(StimSize);   /*tswitch*/
  Rwrite(0.0);
  Rwrite(Tau12);
  Rwrite(Tau13);
  Rwrite(TauRel);
  Rwrite(TABT);
  Rwrite(Temp);
  Rwrite(Tend);
  Rwrite(TimeScale);
  Iwrite(TOmode);
  Rwrite(Tort);
  Rwrite(Tstart);
  for (l = 2; l <= 9; l++)
    Rwrite(0.0);
  Rwrite(TP);
  Rwrite(TPS);   /*tplot1*/
  Rwrite(0.0);   /*tplot2*/
  Rwrite(0.0);   /*tplot3*/
  Rwrite(0.0);
  Rwrite(vol);
  Rwrite(V12);
  Rwrite(V13);
  Rwrite(_Y[5]);
  Rwrite(_Y[4]);
  Rwrite(yNaCa);
  Iwrite(Ymode);

  putc('\n', Filout3);
  if (Filin3 != NULL)
    fclose(Filin3);
  Filin3 = NULL;

}  /* of procedure LISTDATA */



/* $I MENU.PAS*/

/* $I FENTER.PAS*/
/* FEnter.pas */

/* organises the entry of file names in HEART */


void FENTER(void) {
  Char CH;
  Run = 1;
  GetTime(&hour, &minute, &second, &sec100);
  fprintf(Filout3, "Run started at %u hours, %u minutes, %u seconds\n",
	  hour, minute, second);
  /*TIM[0] = ((long)hour * 60L + minute) * 60L + second;*/

  /* read and print heading information in input file */
#ifdef VERBATIM
  CH = '\001';
  while ((CH != '$') & (!P_eof(Filin1))) {
    while ((CH != '$') & (!P_eoln(Filin1))) {
      CH = getc(Filin1);
      if (CH != '$') {
	if (CH != '\n')
	  putc(CH, Filout3);
      }
    }
    if (CH != '$')
      CH = getc(Filin1);
    putc('\n', Filout3);
  }
#else
  while ( !P_eof(Filin1) && ('$'!=(CH=getc(Filin1))) ) putc(CH,Filout3);
  putc('\n',Filout3);
#endif
}  /*of procedure Fenter*/


/* $I FOUT.PAS*/
/* Fout.pas */

/* writes initial output to HEART files */


/* ************************************************************************* */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*             +    +   +++++      +      ++++    +++++++                    */
/*             +    +   +         + +     +   +      +                       */
/*             ++++++   ++++     +   +    ++++       +                       */
/*             +    +   +       +++++++   +   +      +                       */
/*             +    +   +++++  ++     ++  +    +     +                       */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/* ************************************************************************* */


/* $I INIT.PAS*/
/* init.pas */

/* carries out initial computations of parameters used in HEART */

/* NB: nothing to do with INITIAL.PAS ! This procedure was called
   INITIAL until Turbo decided to pre-empt the name ! */


/* ************************************************************************* */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*             +    +   +++++      +      ++++    +++++++                    */
/*             +    +   +         + +     +   +      +                       */
/*             ++++++   ++++     +   +    ++++       +                       */
/*             +    +   +       +++++++   +   +      +                       */
/*             +    +   +++++  ++     ++  +    +     +                       */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/* ************************************************************************* */



/*=========================================================================*/

void INI(void) {

  /* This procedure sets the initial values of various parameters after
     using Stand and Datain and writes some of them to the third output file

Calls: no other procedures
==========================================================================*/
  double Z1;
  short J, L;


  /* Initialise flags and booleans */

  Iflag = 0;   /* force integration to start from scratch */
  NEQND = 0;   /* initialise flag for iflag */

  Head = true;   /* tabulation headings to be written first time */
  Fast = false;   /* assume m not changing rapidly */




  Slow = false;
  if (GNaSS)   /* if GNaSS false then integrate m */
    Fast = false;
  else
    Fast = true;


  /* Initialize times and stimuli */


  Chan = 2;   /* change counter set to 2 */
  Tout = Tstart;   /* initial time */
  T0 = Tstart;
  TAB = Tstart - 0.00001 - TABT * TABT / dt;   /* initialise TAB */
  TPoff = Tend + 1;   /* ensure TPoff is initially off timescale */
  if (TPS > Tstart)   /* stim off at beginning */
    Stim = 0.0;
  if (TPS < Tstart)
    TPS = Tend + 1;
  if (Restart <= Tstart)
    Restart = Tend + 1;
  if (out == 0)
    Restart = Tstart + restrep;
  if (StimSize == 0)   /* no stimulus, so put TPS off scale */
    TPS = Tend + 1;
  PulseLength = off - on;



  /* Initialize voltage and cleft [K] */

  _Y[0] = EC;   /* initial voltage */
/*  Esav = _Y[0] + 1;   / * set ESAV <> y[1] to force rates calculation */
  /*_Y[2]=K*/   /* set external [K] values */
  for (L = 0; L < depth; L++) {
    KC[L] = _Y[2]/*=K*/;
    KCE[L] = _Y[2]/*=K*/;
/*    KCSAV[L] = K; */
  }
  KC[depth] = Kb;
  KCE[depth] = Kb;
/*  KCSAV[depth] = Kb; */
  KC[depth + 1] = Kb;
  KCE[depth + 1] = Kb;
/*  KCSAV[depth + 1] = Kb; */



  /* Initialise diffusion and space parameters */

  Diff = 1500 / (Tort * Tort);
  RTonF = (Temp + 273) * 0.08554;
  Z1 = 0.001 * depth * dx;
  if (Space == 2)
    V = 3.14 * 1.33333 * Z1 * Z1 * Z1;
  else
    V = 3.14 * PrepLength * Z1 * Z1;
  /* use cylindrical geometry except for space = 2 */


  /* For Hilgemann formulation allow for SR volume in calculating vols */

  if (SRmode == 2)
    Vi = V * (1 - vol - V12 - V13);
  else
    Vi = V * (1 - vol);
  V *= vol;
  VF = V * Faraday;

  ViF = Vi * Faraday;
  CAP = 1 / C;

  /* initialise CaCHoff when Hilgemann iCa equations used */
  /*
  if (Camode == 8)
    CaCHoff = _Y[3] / (KCaCHoff + _Y[3]);
	- removed, as CaCHoff made local variable VNB
  */


  /* initialise mitochondrial calcium system */

  if (Mmode == 1) {
    Vmit = 0.4 * Vi;
    Vi = 0.6 * Vi;
    ViF = Vi * Faraday;
    alpha[24] = Vi * alpha[24];
    beta[24] = Vi * beta[24];
  }


  /* initialisation of subsarcolemmal space parameters */

  if (SPvol != 0) {
    Vspace = SPvol * Vi;
    VSpaceF = Vspace * Faraday;
    Vi = (1 - SPvol) * Vi;
  }


  /* initialize reticulum parameters */

  if (SRmode == 2) {
    SRvol = V / vol * V12;
    Relvol = V / vol * V13;
  } else {
    SRvol = Vi * V12;
    Relvol = Vi * V13;
  }
  V12F = SRvol * Faraday;
  V13F = Relvol * Faraday;


  /* initialise parameters for DiFrancesco-Noble model */

  if (SRmode == 1) {
    alpha[11] = 2 * ViF / (Tau12 * Ca12m);
    beta[11] = 0.000001 * alpha[11] * Ca12m;
    alpha[12] = 2 * V13F / Tau13;
    beta[12] = 2 * V13F / TauRel;
  }
  putchar('\n');


  /* Initialize calcium buffers */

  switch (Bmode) {

  case 1:
    _Y[19] = alpha[19] * _Y[3] / (beta[19] + alpha[19] * _Y[3]);
    _Y[20] = alpha[20] * _Y[3] / (beta[20] + alpha[20] * _Y[3]);
    _Y[21] = alpha[21] * Mg - alpha[21] * Mg * _Y[3];
    _Y[21] /= alpha[21] * Mg + beta[21];
    _Y[22] = alpha[22] * _Y[3] + beta[22] * alpha[21] * Mg / beta[21] + beta[22];
    _Y[22] = alpha[22] * _Y[3] / _Y[22];

    break;


  case 2:
    _Y[30] = _Y[3];
    _Y[21] = 0.0;
    _Y[22] = 0.0;
    _Y[33] = 0.0;

    break;


  }/* of Bmode cases */


  /* Initialize external [Ca] */

  if (CaOmode == 1) {
    CaB = _Y[23]/* =Cao */;
  } /*else
    _Y[23] = 0.0;*/


  /* output relevant parameters to third output file */

  if (Run != 1)
    return;

  fprintf(Filout3, "  PARAMETERS COMPUTED AT START OF COMPUTATION: \n\n");
  fprintf(Filout3, "Total volume of preparation = %15.12f microlitres\n",
	  V + Vi);
  fprintf(Filout3, "Extracellular volume = %15.12f  microlitres \n", V);
  fprintf(Filout3, "Representation of extracellular space is ");
  switch (Space) {

  case 1:
    fprintf(Filout3, "cylindrical");
    break;

  case 2:
    fprintf(Filout3, "spherical");
    break;

  case 3:
    fprintf(Filout3, "uniform -- 3 COMPARTMENT MODEL");
    break;

  case 4:
    fprintf(Filout3, "not restricted -- SINGLE CELL");
    break;
  }
  putc('\n', Filout3);

  if (SRmode == 1) {
    fprintf(Filout3, "Ca release power = %d\n", Nrel);
    Z1 = PoN(KmCa, Nrel);
    fprintf(Filout3, "REAL KMCA = %15.12f\n", Z1);
    fprintf(Filout3, "Ca uptake rate = %12.5f\n", alpha[11]);
    fprintf(Filout3, "Ca repriming rate = %12.5f\n", alpha[12]);
    fprintf(Filout3, "Ca release rate = %12.5f\n", beta[12]);
  }

  if (Bmode > 0) {
    fprintf(Filout3, " Initial buffer states: ");
    for (J = 19; J <= 22; J++)
      fprintf(Filout3, "%8.4f", _Y[J]);
    putc('\n', Filout3);
  }

  putc('\n', Filout3);
}  /*of procedure INI  ==================================================*/




/* $I FIRST.PAS*/
/* First.pas */

/* computes and sets some initial parameters -- see also init.pas */


/* ************************************************************************* */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*             +    +   +++++      +      ++++    +++++++                    */
/*             +    +   +         + +     +   +      +                       */
/*             ++++++   ++++     +   +    ++++       +                       */
/*             +    +   +       +++++++   +   +      +                       */
/*             +    +   +++++  ++     ++  +    +     +                       */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/* ************************************************************************* */


/*==========================================================================*/

void FIRST(void) {

  /* Computes first values of concentrations, stimuli & pulses
Calls: NAISS, CAISS, IZERO

===========================================================================*/
  /*double Z1;*/


  /* check whether stim or ipulse starts at tstart */

  T0 = Tstart;
  NEQNN = neqn;

  if (StimSize != 0) {
    if (TPS == Tstart) {
      printf("\n             STIMULUS ON  >>>    \n\n");
      Stim = StimSize;
      TPoff = TP + Tstart;
      TPS = Tend + 1;
      dt = dtss;
    }
  }
  if (on == Tstart) {
    if (off - on < 5 * dt)
      dt = dtss;
    iPulse = iPulseSize;
    on = Tstart + Rep;
    printf("\n             CURRENT PULSE 0N  >>>   \n\n");
  }


  /* check whether [Ca]i is buffered to minimum value */

  if (CaLim >= _Y[3]) {
    _Y[3] = 0.99 * CaLim;
    CaBuff = true;
  }


}  /*of procedure First ==================================================*/







/* End. */
