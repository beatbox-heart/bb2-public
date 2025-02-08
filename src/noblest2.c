/**
 * Copyright (C) (2010-2025) Vadim Biktashev, Irina Biktasheva et al. 
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
#include "system.h"
#include "p2c.h"
#include "noble.h"
#include "nobleini.h"


#define OWN
#include "noblest2.h"


/* Controls the input of data to HEART from data files and interrupts */
static double _Z1/*, _Z2*/;
						#undef _L
						#undef _S
static short _J, _L;
static Char _CHI[60];
static Char _S[256];

static void convert(void) {

  /*This procedure converts up to three characters into decimal equivalent
and the result is stored in _L*/
  short k;

  k = 1;
  do {
    k++;
  } while (_CHI[k - 1] < '1' || _CHI[k - 1] > '9');
  _L = _CHI[k - 1] - 48;
  if (_CHI[k] == ' ' || _CHI[k] == ']')
    return;
  _L = _L * 10 + _CHI[k] - 48;
  if (_CHI[k + 1] != ' ' && _CHI[k + 1] != ']')
    _L = _L * 10 + _CHI[k + 1] - 48;
}  /*convert*/


static void SetBoolean(boolean *b) {
  if (_Z1 > 0)
    *b = true;
  else
    *b = false;
}

static double deltashift5(void) {
  return (7.2161 * PoN(ACh, 0.689928) /
	  (PoN(AChRecfKm, 0.689928) + PoN(ACh, 0.689928)));
}

static void GetDataString(Char *st) {
  _L = 1;
  *st = '\0';
  do {
    snprintf(st + strlen(st), 256, "%c", _CHI[_L - 1]);
    _L++;
  } while (_CHI[_L - 1] != ' ' && _CHI[_L - 1] != '[' && _CHI[_L - 1] >= 'A');
}  /* of procedure GetDataString */


static void chsearch(void) {
  Char _S[256];

  GetDataString(_S);

  /* find variable name and set to the input number */

  if (!Strcmp(_S, "ach")) {
    if (ACh != 0)
      Shift[4] += deltashift5();
    ACh = 0.000001 * _Z1;   /* Read as uM, convert to M */
    if (ACh != 0)
      Shift[4] -= deltashift5();
  }
  if (!Strcmp(_S, "achgkmax"))
    AChgKmax = _Z1;
  if (!Strcmp(_S, "achreckm"))
    AChRecKm = 0.000001 * _Z1;
  if (!Strcmp(_S, "achrecckm"))
    AChRecCKm = 0.000001 * _Z1;
  if (!Strcmp(_S, "achrecfkm")) {
    Shift[4] += deltashift5();
    AChRecfKm = 0.000001 * _Z1;
    Shift[4] -= deltashift5();
  }
  if (!Strcmp(_S, "amode"))
    Amode = _J;
  if (!Strcmp(_S, "atp"))
    _Y[34] = _Z1;
  if (!Strcmp(_S, "alpha")) {
    convert();
    alpha[_L - 1] = _Z1;
  }
  if (!Strcmp(_S, "beta")) {
    convert();
    beta[_L - 1] = _Z1;
  }
  if (!Strcmp(_S, "bmode"))
    Bmode = _J;
  if (!Strcmp(_S, "bufvol"))
    BufVol = _Z1;
  if (!Strcmp(_S, "buffast"))
    SetBoolean(&BufFast);

  if (!Strcmp(_S, "cbden"))
    CBden = _Z1;
  if (!Strcmp(_S, "cbturn"))
    CBturn = _Z1;
  if (!Strcmp(_S, "computecaintegral"))
    SetBoolean(&ComputeCaIntegral);
  if (!Strcmp(_S, "cape"))
    CaPE = _Z1;
  if (!Strcmp(_S, "capacitance"))
    C = _Z1;
  if (!Strcmp(_S, "cai"))
    _Y[3] = _Z1;
  if (!Strcmp(_S, "cao")) {
			/* Cao now made local variable VNB */
    _Y[23] = /*Cao =*/ _Z1;
  }
  if (!Strcmp(_S, "caomode"))
    CaOmode = _J;
  if (!Strcmp(_S, "calim"))
    CaLim = _Z1;
  if (!Strcmp(_S, "calmod"))
    Calmod = _Z1;
  if (!Strcmp(_S, "camode"))
    Camode = _J;
  if (!Strcmp(_S, "canmode"))
    CaNmode = _J;
  if (!Strcmp(_S, "cass"))
    CaSS = true;
  if (!Strcmp(_S, "cact")) {
    _L = _CHI[4];
    CACT[_L - 49] = _Z1;
  }
  if (!Strcmp(_S, "cavol"))
    CaVOL = _Z1;
  if (!Strcmp(_S, "computeatp"))
    SetBoolean(&ComputeATP);
  if (!Strcmp(_S, "contractmode"))
    SetBoolean(&ContractMode);
  if (!Strcmp(_S, "ctrop"))
    CTrop = _Z1;
  if (!Strcmp(_S, "ca")) {
    if (_CHI[3] == '3')
      _Y[12] = _Z1;
    if (_CHI[3] == '2') {
      if (_CHI[4] == 'M')
	Ca12m = _Z1;
      else
	_Y[11] = _Z1;
    }
    if (_CHI[2] > '1') {
      _L = _CHI[2];
      CaCH[_L - 49] = _Z1;
    }
    if (_CHI[2] == '[') {
      convert();
      _Y[_L - 1] = _Z1;
    }
  }
  if (!Strcmp(_S, "cycle"))
    Rep = _Z1;
  if (!Strcmp(_S, "d")) {
    if (_CHI[1] == ' ')
      _Y[6] = _Z1;
    else
      _Y[15] = _Z1;
  }
  if (!Strcmp(_S, "dx"))
    dx = _Z1;

  if (!Strcmp(_S, "diffca"))
    DiffCa = _Z1;
  if (!Strcmp(_S, "dnaca"))
    DNaCa = _Z1;
  if (!Strcmp(_S, "dvtest"))
    dVtest = _Z1;
  if (!Strcmp(_S, "dt")) {
    if (_CHI[2] > '1') {
      _L = _CHI[2];
      dtCH[_L - 49] = _Z1;
    } else {
      dt = _Z1;
      dtl = dt;
      dts = 0.5 * dt;
      dtss = 0.1 * dt;
    }
  }
  if (!Strcmp(_S, "dts"))
    dts = _Z1;
  if (!Strcmp(_S, "dtss"))
    dtss = _Z1;

  if (!Strcmp(_S, "e")) {
    if (_CHI[1] == 'c' || _CHI[1] == ' ')
      _Y[0] /*!!= EC!!*/ = _Z1;
  }

  if (!Strcmp(_S, "f")) {
    if (_CHI[1] == ' ')
      _Y[7] = _Z1;
    if (_CHI[1] == '2')
      _Y[14] = _Z1;
    if (_CHI[1] == '3')
      _Y[16] = _Z1;
  }
  if (!Strcmp(_S, "fura"))
    Fura = _Z1;
  if (!Strcmp(_S, "furas"))
    FuraS = _Z1;

  if (!Strcmp(_S, "gb"))
    gB = _Z1;
  if (!Strcmp(_S, "gbna"))
    gbNa = _Z1;
  if (!Strcmp(_S, "gbk"))
    gbK = _Z1;
  if (!Strcmp(_S, "gbca"))
    gbCa = _Z1;
  if (!Strcmp(_S, "gfk"))
    gfK = _Z1;
  if (!Strcmp(_S, "gfna"))
    gfNa = _Z1;
  if (!Strcmp(_S, "gk")) {
    if (_CHI[2] == ' ')
      gK = _Z1;
    else
      gK1 = _Z1;
  }
  if (!Strcmp(_S, "gkatpm"))
    gkATPm = _Z1;
  if (!Strcmp(_S, "gna"))
    gNa = _Z1;
  if (!Strcmp(_S, "gnak"))
    gNaK = _Z1;
  if (!Strcmp(_S, "gnass"))
    SetBoolean(&GNaSS);
  if (!Strcmp(_S, "gto"))
    gTO = _Z1;

  if (!Strcmp(_S, "h"))
    _Y[8] = _Z1;
  if (!Strcmp(_S, "hiacc"))
    SetBoolean(&hiAcc);

  if (!Strcmp(_S, "iachm"))
    iAChm = _Z1;
  if (!Strcmp(_S, "icafract"))
    iCafract = _Z1;
  if (!Strcmp(_S, "ikm"))
    iKm = _Z1;
  if (!Strcmp(_S, "ip") || !Strcmp(_S, "ipulse") || !Strcmp(_S, "ipulsesize"))
    iPulseSize = _Z1;
  if (!Strcmp(_S, "incfract"))
    iNCfract = _Z1;
  if (!Strcmp(_S, "inacam"))
    iNaCam = _Z1;

  if (!Strcmp(_S, "iscale"))
    Iscal = _J;

  if (!Strcmp(_S, "k")) {
    _L = _CHI[1];
    KCH[_L - 49] = _Z1;
  }
  if (!Strcmp(_S, "katp"))
    kATP = _Z1;
  if (!Strcmp(_S, "kb"))
    Kb = _Z1;
  if (!Strcmp(_S, "ki"))
    _Y[10] = _Z1;
  if (!Strcmp(_S, "kc")) {
    _Y[2]/*=K*/ = _Z1;
  }
  if (!Strcmp(_S, "kcalon"))
    alpha[19] = _Z1;
  if (!Strcmp(_S, "kcaloff"))
    beta[19] = _Z1;
  if (!Strcmp(_S, "kctropon"))
    alpha[20] = _Z1;
  if (!Strcmp(_S, "kctropoff"))
    beta[20] = _Z1;
  if (!Strcmp(_S, "kcachoff"))
    KCaCHoff = _Z1;
  if (!Strcmp(_S, "kcyca"))
    KCyCa = _Z1;
  if (!Strcmp(_S, "knaca"))
    kNaCa = _Z1;
  if (!Strcmp(_S, "kpe"))
    KPE = _Z1;
  if (!Strcmp(_S, "kpnaca"))
    KPNaCa = _J;
  if (!Strcmp(_S, "kmca")) {
    if (_CHI[4] == '2')
      KmCa2 = _Z1;
    else
      KmCa = _Z1;
  }
  if (!Strcmp(_S, "kminact"))
    Kminact = _Z1;
  if (!Strcmp(_S, "kmf"))
    Kmf = _Z1;
  if (!Strcmp(_S, "kmk"))
    KmK1 = _Z1;
  if (!Strcmp(_S, "kmnaca"))
    KMNaCa = _Z1;
  if (!Strcmp(_S, "kmslpump"))
    KMSLpump = _Z1;
  if (!Strcmp(_S, "kmto"))
    KmTO = _Z1;
  if (!Strcmp(_S, "kmtropmon"))
    alpha[21] = _Z1;
  if (!Strcmp(_S, "kmtropmoff"))
    beta[21] = _Z1;
  if (!Strcmp(_S, "kmtropcon"))
    alpha[22] = _Z1;
  if (!Strcmp(_S, "kmtropcoff"))
    beta[22] = _Z1;
  if (!Strcmp(_S, "km"))
    Km = _Z1;
  if (!Strcmp(_S, "kmode"))
    Kmode = _J;
  if (!Strcmp(_S, "kslpump"))
    KSLpump = _Z1;
  if (!Strcmp(_S, "ksrca"))
    KSRCa = _Z1;
  if (!Strcmp(_S, "kxcs"))
    Kxcs = _Z1;

  if (!Strcmp(_S, "_L"))
    PrepLength = _Z1;

  if (!Strcmp(_S, "m"))
    _Y[9] = _Z1;

  if (!Strcmp(_S, "mtrop"))
    MTrop = _Z1;
  if (!Strcmp(_S, "mib"))
    MiB = _Z1;
  if (!Strcmp(_S, "mik")) {
    if (_CHI[3] == '1')
      MiK1 = _Z1;
    else
      MiK = _Z1;
  }
  if (!Strcmp(_S, "mito"))
    MiTO = _Z1;

  if (!Strcmp(_S, "na")) {
    _L = _CHI[2];
    NaCH[_L - 49] = _Z1;
  }
  if (!Strcmp(_S, "nai"))
    _Y[1] = _Z1;
  if (!Strcmp(_S, "nao"))
    Nao = _Z1;
  if (!Strcmp(_S, "namode"))
    Namode = _J;
  if (!Strcmp(_S, "nmode"))
    Nmode = _J;
  if (!Strcmp(_S, "nape"))
    NaPE = _Z1;
  if (!Strcmp(_S, "nnaca"))
    NNaCa = _J;
  if (!Strcmp(_S, "nnak"))
    nNaK = _Z1;

  if (!Strcmp(_S, "nrel"))
    Nrel = _J;

  if (!Strcmp(_S, "out"))
    out = _J;
  if (!Strcmp(_S, "on"))
    on = _Z1;
  if (!Strcmp(_S, "off")) {
    off = _Z1;
    PulseLength = off - on;
  }

  if (!Strcmp(_S, "pf"))
    PF = _Z1;
  if (!Strcmp(_S, "pca")) {
    if (_CHI[3] == '3')
      PCa3 = _Z1;
    if (_CHI[3] == '2')
      PCa2 = _Z1;
    if (_CHI[3] == ' ')
      PCa = _Z1;
  }
  if (!Strcmp(_S, "pcak"))
    PCaK = _Z1;

  if (!Strcmp(_S, "pnak"))
    PNaK = _Z1;

  if (!Strcmp(_S, "prim"))
    _Y[13] = _Z1;
  if (!Strcmp(_S, "preplength"))
    PrepLength = _Z1;
  if (!Strcmp(_S, "prep"))
    Prep = _J;
  if (!Strcmp(_S, "pump")) {
    if (_CHI[4] == ' ')
      Pump = _Z1;
    else {
      _L = _CHI[4];
      PumpCH[_L - 49] = _Z1;
    }
  }
  if (!Strcmp(_S, "q"))
    _Y[18] = _Z1;
  if (!Strcmp(_S, "r"))
    _Y[17] = _Z1;
  if (!Strcmp(_S, "ramp"))
    SetBoolean(&ramp);
  if (!Strcmp(_S, "rep")) {
    if (_CHI[3] == ' ')
      Rep = _Z1;
    else {
      _L = _CHI[3];
      repCH[_L - 49] = _Z1;
    }
  }
  if (!Strcmp(_S, "restart"))
    Restart = _Z1;
  if (!Strcmp(_S, "restrep"))
    restrep = _Z1;
  if (!Strcmp(_S, "restfileno"))
    restfileno = _J;
  if (!Strcmp(_S, "rmode"))
    Rmode = _J;
  if (!Strcmp(_S, "sarcolength"))
    SarcoLength = _Z1;
  if (!Strcmp(_S, "scale"))
    Scale = _Z1;
  if (!Strcmp(_S, "shifty"))
    Shift[4] = _Z1;
  if (!Strcmp(_S, "shiftx"))
    Shift[5] = _Z1;
  if (!Strcmp(_S, "shiftd"))
    Shift[6] = _Z1;
  if (!Strcmp(_S, "shiftf"))
    Shift[7] = _Z1;
  if (!Strcmp(_S, "shifth"))
    Shift[8] = _Z1;
  if (!Strcmp(_S, "shiftm"))
    Shift[9] = _Z1;
  if (!Strcmp(_S, "shiftto"))
    ShiftTO = _Z1;
  if (!Strcmp(_S, "shiftk"))
    ShiftK1 = _Z1;
  if (!Strcmp(_S, "shift")) {
    convert();
    Shift[_L - 1] = _Z1;
  }

  if (!Strcmp(_S, "slpumpmode"))
    SetBoolean(&SLPumpMode);
  if (!Strcmp(_S, "speedy"))
    Speed[4] = _Z1;
  if (!Strcmp(_S, "speedx"))
    Speed[5] = _Z1;
  if (!Strcmp(_S, "speedd"))
    Speed[6] = _Z1;
  if (!Strcmp(_S, "speedf"))
    Speed[7] = _Z1;
  if (!Strcmp(_S, "speedh"))
    Speed[8] = _Z1;
  if (!Strcmp(_S, "speedm"))
    Speed[9] = _Z1;
  if (!Strcmp(_S, "speed")) {
    convert();
    Speed[_L - 1] = _Z1;
  }
  if (!Strcmp(_S, "spmimic"))
    SetBoolean(&spmimic);
  if (!Strcmp(_S, "space"))
    Space = _J;
  if (!Strcmp(_S, "spvol"))
    SPvol = _Z1;
  if (!Strcmp(_S, "srmode"))
    SRmode = _J;
  if (!Strcmp(_S, "srleak"))
    SRleak = _Z1;
  if (!Strcmp(_S, "srfract"))
    SRfract = _Z1;
  if (!Strcmp(_S, "stimsize"))
    StimSize = _Z1;
  if (!Strcmp(_S, "steepk"))
    steepK1 = _Z1;

  if (!Strcmp(_S, "tnap"))
    TNaP = _Z1;
  if (!Strcmp(_S, "tkp"))
    TKP = _Z1;
  if (!Strcmp(_S, "tcap"))
    TCaP = _Z1;
  if (!Strcmp(_S, "tauf"))
    alpha[14] = 1 / _Z1;
  if (!Strcmp(_S, "tauspvol"))
    TauSPvol = _Z1;
  if (!Strcmp(_S, "tauup"))
    Tau12 = _Z1;
  if (!Strcmp(_S, "taureprim"))
    Tau13 = _Z1;
  if (!Strcmp(_S, "taurelease"))
    TauRel = _Z1;
  if (!Strcmp(_S, "tabt")) {
    TABT = _Z1;
    if (dt > TABT) {
      dt = TABT;
      dts = 0.5 * dt;
      dtss = 0.1 * dt;
      dtl = dt;
    }
  }
  if (!Strcmp(_S, "temperature"))
    Temp = _Z1;
  if (!Strcmp(_S, "tend"))
    Tend = _Z1;
  if (!Strcmp(_S, "timescale"))
    TimeScale = _Z1;
  if (!Strcmp(_S, "tort"))
    Tort = _Z1;
  if (!Strcmp(_S, "tomode"))
    TOmode = _J;
  if (!Strcmp(_S, "tstart"))
    Tstart = _Z1;
  if (!Strcmp(_S, "tp"))
    TP = _Z1;
  if (!Strcmp(_S, "tps"))
    TPS = _Z1;



  if (!Strcmp(_S, "vsurfca"))
    VSurfCa = _Z1;
  if (!Strcmp(_S, "vol"))
    vol = _Z1;
  if (!Strcmp(_S, "v")) {
    if (_CHI[2] == '2')
      V12 = _Z1;
    if (_CHI[2] == '3')
      V13 = _Z1;
  }



  if (!Strcmp(_S, "x"))
    _Y[5] = _Z1;

  if (!Strcmp(_S, "y")) {
    if (_CHI[1] == ' ')
      _Y[4] = _Z1;
    else {
      convert();
      _Y[_L - 1] = _Z1;
    }
  }
  /*if (!Strcmp(_S, "ymax"))   // save to use with ymin 
    _Z2 = _Z1;*/

  if (!Strcmp(_S, "ynaca"))
    yNaCa = _Z1;



}  /*of Procedure CHSEARCH */


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void DATAIN(void) {  /*of procedure DATAIN*/

  /* this procedure reads changes in input data either from an input
     file or from data entered during an interrupt. The integer variable
     ii is set to:

     1 -- for data entered during interrupts (for HEART)
     2 -- for data entered during interrupts in graphics mode (HEART)
     3 -- for data read from an input file (for HEART)
     4 -- for data read from a graphics parameter file (for HPLOT)
     5 -- for data read from automatic interrupts (for HEART)


Contains: CharIn, CHsearch, Convert, SetBoolean, SetRange, Setynn,
             Deltashift5, GetDataString

Uses: FileOpen (Scr),

==========================================================================*/
  Char *TEMP1;


_L2:
  _L = 1;

  _CHI[0] = '\001';

  /* look for letters */

  while ((_CHI[0] > 'z' || _CHI[0] < 'A') && _CHI[0] != '\\' &&
	 _CHI[0] != '(' &&
	 _CHI[0] != '{' && _CHI[0] != '%' && _CHI[0] != '$') {
    *_CHI = getc(Filin1);
    if (isupper(_CHI[0]))
      _CHI[0] += 32;
    if (_CHI[0] != '\\' && _CHI[0] != '(' && _CHI[0] != '{' &&
	_CHI[0] != '%' && _CHI[0] != '$')
      putc(_CHI[0], Filout3);
  }
  if (_CHI[0] == '$')
    goto _L3;

  /* write comments to filout3 */

  if (_CHI[0] == '\\' || _CHI[0] == '(' || _CHI[0] == '{' || _CHI[0] == '%') {
    fgets(_S, 256, Filin1);
    TEMP1 = strchr(_S, '\n');
    if (TEMP1 != NULL)
      *TEMP1 = 0;
    goto _L2;
  }

  /* look for equals sign */

  while (_CHI[_L - 1] != '=') {
    _L++;
    _CHI[_L - 1] = getc(Filin1);
    if (isupper(_CHI[_L - 1]))
      _CHI[_L - 1] += 32;
    putc(_CHI[_L - 1], Filout3);
    if (_CHI[_L - 1] == '*')
      goto _L2;
  }

  /* read input parameter or string */


  fscanf(Filin1, "%lg", &_Z1);
  fprintf(Filout3, "%12.9f", _Z1);
  _J = (short)floor(_Z1 + 0.5);
  chsearch();



  if (_CHI[0] != '$')
    goto _L2;

  /* end of input of new data */


_L3:
  if (P_eof(Filin1))
    Iflag = 3;
  Iflag = 0;
}  /*of procedure datain ================================================*/


/* $I STAND.PAS*/
/* Stand.pas */

/* Sets standard values of parameters in HEART */


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


/*============================================================================*/

void STANDARD(void) {

  /*This procedure sets standard values of parameters usually common to all
preparations. These are the values used unless changed either by calling
a specific preparation procedure, by reading an input file, or setting a
parameter using the menus.

Calls: no other procedures

=============================================================================*/
  short /*J,*/ L;


  /* initialise all neqn arrays*/

  for (L = 0; L < neqn; L++) {
    _Y[L] = 0.0;
    _F[L] = 0.0;
    alpha[L] = 1.0;
    beta[L] = 1.0;
    Shift[L] = 0.0;
    Speed[L] = 1.0;
#ifdef HEART
    saven0[L] = 0.0;
    saven1[L] = 0.0;
    saven2[L] = 0.0;
    saven3[L] = 0.0;
    savY[L] = 0.0;
    Yres[L] = 0.0;
#endif
  }


  /* zero all save parameters*/
#ifdef HEART
  save1 = 0.0;
  save2 = 0.0;
  save3 = 0.0;
  save4 = 0.0;
  save5 = 0.0;
  save6 = 0.0;
  save7 = 0.0;
  save8 = 0.0;
  save9 = 0.0;
  save10 = 0.0;
  save11 = 0.0;
  save12 = 0.0;
  Isave1 = 0;
  Isave2 = 0;
  Isave3 = 0;
  Isave4 = 0;
  itot = 0.0;   /* zero total current */
#endif

  /* set modes to default values.
     note that many of these modes are reset by the preparation
     procedures. Some standard values are now only rarely used */


  Rmode = 0;   /* voltage dependent repriming */
  Camode = 1;   /* original Oxford equations for i Ca */
  Kmode = 0;   /* MNT model equations */
  Amode = 0;   /* no calcium activation */

  Nmode = 1;   /* standard DiFrancesco-Noble equations for i NaCa*/
  TOmode = 1;   /* inactivation of i to*/
  CaNmode = 1;   /* original Reuter-Scholz model of i Ca */
  SRmode = 1;   /* DiFrancesco-Noble SR equations */
  Ymode = 0;   /* DiFrancesco-Noble i f equations */
  Bmode = 0;   /* no calcium buffers */
  Namode = 1;   /* DiFrancesco-Noble i Na equations */
  Mmode = 0;   /* no mitochondria */
  CaOmode = 0;   /* don't compute [Ca]o */
  ContractMode = false;   /* no contraction */
  BufFast = true;   /* buffer rates computed in relevant models - VNB */
  SLPumpMode = false;   /* no sarcolemmal Ca pump */
  ComputeCaIntegral = false;   /* do not compute this integral */
  ComputeATP = false;   /* do not compute [ATP] */

  Prep = 1;   /* Purkinje fibre */
  Space = 3;   /* cleft space */
  out = 4;   /* output main ionic currents to second data file */
  Temp = 37.0;   /* temperature */
  CaLim = 0.000001;   /* minimal [Ca]i --- no longer used in most cases */
  Tort = 1.41;   /* external diffusion tortuosity factor */
  Kmf = 45.0;   /* Km for K activation of i f mM */
  gB = 0.0;   /* background conductance */
  dVtest = 200.0;   /* switch to full m equations at 200V/s */
  ACh = 0.0;   /* ACh concentration zero */




  CaSS = false;   /* no steady state assumptions for [Ca]i */

  GNaSS = false;   /* no steady state assumption for i Na */
  CaBuff = false;   /* fast (simplified) calcium buffer equations */


  ramp = false;   /* no ramp clamps */
  Run = 1;   /* initialise run counter */
  hiAcc = false;   /* do not use DTSS */


  /* initialise stimuli and pulses */

  Stim = 0.0;   /* stimulus amplitude zero */
  iPulse = 0.0;   /* I pulse zero */
  iPulseSize = 0.0;   /* I pulse zero */

  StimSize = 0.0;
  TP = 0.0;

  /* initialise times and related parameters.
     Note that parameters that are not to be used are set to
     'infinite' time. This is defined in initial.pas: it is a
     very large number, not actually infinite */

  Tstart = 0.0;   /* start at zero */
  T0 = 0.0;   /* first step from zero */

  Restart = Infinite;   /* restart time -> infinite */
  restfileno = 1;   /* first restart file for OUT = 0 */

  on = Infinite;   /* pulse on time -> infinite */
  off = Infinite;   /* pulse off time -> infinite */
  Rep = Infinite;   /* pulse repeat time -> infinite */
  TPS = 0.0;   /* stim time at zero */



  /* scaling parameters */

  TimeScale = 1.0;   /* tabulate time in secs */
  Iscal = 1;   /* tabulate current in nA */
  Scale = 1.0;   /* general scaling set to 1 */

  /* Calcium buffer parameters */

  BufVol = 0.5;   /* fractional space occupied by buffers = 0.5 */
  CaVOL = 1.0;   /* fractional space occupied by calcium = 1 */
  Ncalmod = 3;   /* 3 Ca sites on calmodulin */
  Nctrop = 1;   /* 1 Ca site on c troponin */
  Nmtrop = 2;   /* 2 Ca sites on m troponin */
  alpha[20] = 100000.0;   /* on rate for troponin 18^8 M^-1 sec^-1*/
  /* note: 10^5 actually used since concn in mM */
  beta[20] = 200.0;   /* off rate for troponin 200 sec ^-1 */
  _Y[20] = 0.0015;   /* fraction of troponin sites occupied */
  alpha[19] = 100000.0;
      /* on rate for calmodulin, same units as alpha[21] */
  beta[19] = 50.0;   /* off rate for calmodulin 50 sec ^-1 */
  _Y[19] = 0.0005;   /* fraction of calmodulin sites occupied */
  _Y[21] = 0.0;   /* Mg bound to troponin */
  _Y[22] = 0.0;   /* Ca bound to M troponin */
  Mg = 2.5;   /* [Mg]i = 2.5 mM */




  /* contraction parameters */

  alpha[28] = 12000.0;
  beta[28] = 100.0;
  _Y[28] = 0.0;   /* light chain conformation */
  alpha[29] = 60.0;
  beta[29] = 25.0;
  _Y[29] = 0.0;   /* cross bridge reaction */
  SarcoLength = 0.0;
  /* CBdenA = 0.0; - this is a strange object -
    const parameter OR intermediate variable depending on flags!!
    Let it be intermediate - see when it will not be enough VNB */
  /* IsoTen = 0.0; - this initial value not used actually VNB */

  /* original DiF-N SR parameters */

  Nrel = 2;   /* 2 Ca ions bind to release site */
  Ca12m = 5.0;   /* maximum [Ca] 5 mM */
  V12 = 0.05;   /* volume for uptake SR 5 % */
  V13 = 0.02;   /* volume for release SR 2 % */
  Tau12 = 0.005;   /* SR uptake constant */
  Tau13 = 0.2;   /* SR transfer constant */
  TauRel = 0.01;   /* SR release constant */



  /* Hilgemann SR parameters */

  alpha[11] = 3.0;   /* rate constant for uptake to SR */
  beta[11] = 0.23;   /* rate constant for release from SR */
  _Y[11] = 0.3;   /* uptake calcium 0.3 mM */
  alpha[12] = 50.0;   /* rate constant for transfer from UP to REL */
  _Y[12] = 0.3;   /* release calcium 0.3 mM */
  SRleak = 0.0;   /* no leak in original Hilgemann SR model */
  alpha[25] = 600.0;   /* rate constant for V dep formation of SR activator */
  beta[25] = 500.0;   /* rate constant for Ca-dep formation of SR activator */
  alpha[26] = 60.0;   /* constant for Ca-ind inactivation rate */
  beta[26] = 500.0;   /* constant for Ca-dep inactivation rate */
  _Y[26] = 0.0;   /* SR release activator fraction */
  beta[27] = 1.0;   /* rate for regeneration of precursor */
  _Y[27] = 0.0;   /* product fraction */
  /* precursor fraction y[26] = 1 - y[27] - y[28] = 1 */
  KmCa2 = 0.0;   /* rate constant for release set to zero to avoid overflow */

  /* i Na */

  PNaK = 0.12;   /* K permeability of Na channel 12 % */


  /* i K1 parameters */

  steepK1 = 2.0;   /* rectification steepness 2 */
  KmK1 = 10.0;   /* Km for K activation of i K1 */
  ShiftK1 = 0.0;   /* standard voltage dependence of rectifier */


  /* i TO parameters */

  ShiftTO = 0.0;   /* standard voltage dependence */
  KmTO = 10.0;   /* Km for K activation 10 mM */
  _Y[17] = 1.0;   /* inactivation r */
  _Y[18] = 0.0;   /* activation q */


  /* These variables are used when computing extracellular [Ca] */

  DiffCa = 0.0005;


  /* data for calcium buffering from Robertson, Johnson & Potter,
     Biophysical Journal, (1981), 34, 559-569 */

  Calmod = 0.00767;   /* [calmodulin] = 7.67 uM */
  CTrop = 0.049;   /* [ctroponin] = 49 uM */
  MTrop = 0.0429;   /* [mtroponin] = 42.9 uM */

  /* old buffer parameters: omit these while using versions 2 and 3 */

  /* alpha[20] := 100000.0;  beta[20] := 238;
     alpha[21] := 3900;      beta[21] := 19.6;
     alpha[22] := 10;        beta[22] := 33.3;
     alpha[23] := 10000;     beta[23] := 0.33;
    */

  /* Data for calcium probe __ see Jackson et al FEBS LETTERS 216, 35-39, 1*/

  Fura = 0.0;   /* [fura] = 0 mM */
  FuraS = 0.0;   /* SS space [fura] = 0 mM */
  alpha[33] = 100000.0;   /* on rate 100,000 */
  beta[33] = 84.0;   /* off rate 84 sec -1 */


  /* parameters for mitochondrial calcium */

  alpha[24] = 1.25;
  beta[24] = 0.15;
  KcMRel = 19.34;   /* Ca release constant */
  KnMRel = 11.36;   /* Ca uptake constant */
  KcMup = 0.0877;   /* constant for binding to M pump */
  NNMRel = 2.7;   /* ions bound by uptake pump */
  NcMup = 1.4;   /* Ca ions binding to M pump */
  _Y[24] = 0.005;   /* [Ca] in M = 5 uM */


  /* sodium pump */

  nNaK = 1.5;   /* sodium pump stoichiometry 3:2 Na:K */
  Km = 1.0;   /* Km for [K]o activation mM */
  KmNa = 40.0;   /* Km for [Na]i activation of sodium pump mM */


  /* sodium-calcium exchange */

  NNaCa = 3;   /* stoichiometry */
  DNaCa = 0.001;   /* denominator constant */
  yNaCa = 0.5;   /* position of peak energy barrier */
  KMNaCa = 0.0;   /* Km for calcium regulatory site */
  KPNaCa = 1;   /* no of Ca ions binding to regulatory site */
  iNaCam = 10000.0;   /* maximum i NaCa */


  /* i Ca,L channel */

  PCaK = 0.01;   /* K perm 1 % of calcium permeability */
  KCaCHoff = 0.001;   /* Ca induced i Ca inactivation constant */
  VSurfCa = 50.0;   /* Ca channel surface V = 50 mV */


  /* i Ca,T channel */

  PCa2 = 0.0;   /* T channel P set to zero */
  _Y[31] = 0.0;   /* T channel activation */
  _Y[32] = 1.0;   /* T channel inactivation */


  /* slowly inactivated Ca channel -- `KS Lee current' */

  PCa3 = 0.0;   /* no `slow' calcium channel permeability */
  _Y[15] = 0.0;   /* channel activation  d3 */
  _Y[16] = 1.0;   /* channel inactivation f3 */


  /* Ca inactivation of i Ca */

  Kminact = 0.001;   /* Km for [Ca]i inactivation of i Ca = 1 uM */
  alpha[14] = 10.0;   /* rate constant for recovery sec^-1 */



  /* default values for calcium-activation of currents */

  CACT[0] = 0.0;   /* no Ca activation of i K1 */
  CACT[1] = 0.0;   /* no Ca activation of i K */
  CACT[2] = 0.0;   /* no Ca activation of i b,Na */
  CACT[3] = 0.0;   /* no Ca activation of i to */
  MiB = 0.0;
  MiK = 0.0;
  MiTO = 0.0;
  MiK1 = 0.0;


  gK = 0.0;   /* Kmode 5 conductance zero */


  /* default values of subsarcolemmal (SS) space parameters */

  SPvol = 0.0;   /* no SS space */
  iCafract = 0.0;   /* no I Ca in space */
  iNCfract = 0.0;   /* no i NaCa in space */
  SRfract = 0.0;   /* no SR in space */
  spmimic = false;   /* no SS space mimic */



  iAChm = 0.0;   /* Ach channel zero */



  /* default values of ATP, KATP channel, CB density and CB turnover */

  CBden = 0.05;   /* cross bridge density 50 uM */
  CBturn = 0.02;   /* cross bridge turnover rate 20 msec */
  _Y[34] = 5.0;   /* [ATP] = 5 mM */
  gkATPm = 0.0;   /* Maximum conductance set to zero */
  /* Sensible values to set are in the range */
  /* 1-2 uS per cell --- Nichols & Lederer, 1990 */
  kATP = 0.1;   /* Binding constant for ATP 100 uM */

  /* Default parameters for Mark Boyett's formulation of iach */

  Shift[50] = 50.0;
  Shift[51] = 50.0;
  /* iACh = 0.0;    - this initial value is not used actually VNB */
  AChgKmax = 0.0;   /* set this to non-zero value to use new formulation */
  alpha[50] = 3.684211;
  alpha[51] = 73.09924;
  _Y[50] = 0.5;
  _Y[51] = 0.5;
  AChRecKm = 2.795e-7;   /* binding constant for iKACh */
  AChRecCKm = 2.0e-7;   /* binding constant for iCa */
  AChRecfKm = 1.256e-8;   /* binding constant for i(f) */


  /* default values of Mark Boyett's formulation of slow iCa inactivation */

  for (L = 52; L <= 55; L++)
    _Y[L] = 1.0;
  Speed[53] = 3.0;
  Speed[54] = 30.0;
  Speed[55] = 150.0;



  /* default values of graphics ranges */



}  /*of procedure STANDARD =================================================*/



/* $I PREPS.PAS*/
/* Preps.pas */

/* defines the various preparation models in HEART */


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



/* These procedures define the standard models so far developed.
   To determine the parameters used the program first sets parameters
   using procedure standard, then adds or alters parameters with
   the various preparation procedures (which are sometimes called
   in succession, e.g. where one model has been developed from another)
   and finally, some parameters are set by calculation in procedure init*/


/*==========================================================================*/

void PURKINJE(void) {

  /*  set standard parameters for Purkinje fibre

This is the original DiFrancesco-Noble (1985) multicellular model

calls: no other procedures

===========================================================================*/
  C = 0.075;   /* Capacitance 75 nF */
  gNa = 750.0;   /* sodium conductance uS */
  PCa = 15.0;   /* Ca permeability x F */
  iKm = 180.0;   /* maximum i K  nA */
  gK1 = 920.0;   /* maximum gK1  uS */
  KmK1 = 210.0;   /* K activation of i K1  mM */
  Camode = 5;   /* Ca-induced Ca inactivation set */
  Amode = 4;   /* Ca-dependent activation of i to */
  CACT[3] = 0.0005;   /* binding constant 0.5 uM */
  Shift[4] = 10.0;   /* i f activation curve shift 10 mV */
  gbNa = 0.18;   /* background sodium conductance 0.18 uS */
  gbK = 0.0;   /* background potassium conductance 0 */
  gfK = 3.0;   /* i f K conductance 3 uS */
  gfNa = 3.0;   /* i f Na conductance 3 uS */
  gTO = 1.0;   /* i to maximum conductance 1 uS */
  gbCa = 0.02;   /* background calcium conductance 0.02 uS */
  Pump = 125.0;   /* maximum pump current 125 nA */
  kNaCa = 0.02;   /* i NaCa scaling factor */
  DNaCa = 0.001;   /* i NaCa denominator factor */
  _Y[10] = 140.0;   /* [K]i  140 mM */
  _Y[1] = 8.0;   /* [Na]i 8 mM */
  _Y[3] = 0.00005;   /* [Ca]i 50 nM */
  Kb = 4.0;   /* [K]b  4 mM */
  _Y[2]/*=K*/ = 4.0;   /* [K]o  4 mM */
  Nao = 140.0;   /* [Na]o 140 mM */
  _Y[23] = /*Cao =*/ 2.0;   /* [Ca]o  2 mM */
  _Y[4] = 0.2;   /* i f activation */
  _Y[5] = 0.01;   /* i K activation */
  _Y[6] = 0.005;   /* i Ca activation */
  _Y[7] = 1.0;   /* i Ca inactivation */
  _Y[14] = 1.0;   /* i Ca Ca-induced inactivation */
  _Y[8] = 0.8;   /* i Na inactivation */
  _Y[9] = 0.01;   /* i Na activation */
  Shift[6] = 5.0;   /* voltage shift of i Ca activation */
  TOmode = 1;   /* i to equations include inactivation */
  dx = 5.0;   /* radial segment width um, total radius 10 x 5 = 50 um */
  PrepLength = 2.0;   /* fibre length 2 mm */
  vol = 0.1;   /* fractional extracellular space volume */
  PF = 0.7;   /* rate constant for diffusion to bulk space sec^-1 */
  Tend = 5.0;   /* terminate computation at 5 sec */
  dt = 0.005;   /* main step length 5 msec */
  TABT = 0.005;   /* main tabulation interval 5 msec */
  TP = 0.005;   /* stimulus duration 5 msec */
  TPS = 0.3;   /* applied at 300 msec */
  dVtest = 200.0;   /* switch to full i Na activation at 200 v/sec */
  EC = -87.0;   /* initial voltage */
  StimSize = -350.0;   /* stimulus 350 nA */
  TPoff = Tend + 1;   /* ensure TPoff never reached */
  _Y[11] = 2.0;   /* [Ca]sr uptake store 2 mM */
  _Y[12] = 1.0;   /* [Ca]sr release store 1 mM */
  _Y[13] = 1.0;   /* repriming */
  Tau12 = 0.025;   /* SR uptake constant */
  Tau13 = 2.0;   /* SR transfer constant */
  TauRel = 0.05;   /* SR release constant */
  KmCa = 0.001;   /* Ca binding constant to release site 1 uM */
  Prep = 1;   /* Purkinje fibre */
  CaNmode = 0;   /* no sodium permeability in i Ca,L channel */
  alpha[14] = 5.0;   /* rate constant for recovery from Ca inact 5 sec^-1 */


}  /*of procedure PURKINJE ==============================================*/





void SINUS(void) {


  /*========================================================================*/

  /*   set standard parameters for sinus node as at 20 November 1983

This is the original multicellular sinus model (Noble & Noble 1984)

calls: no other procedures

=========================================================================*/
  C = 0.006;   /* Capacitance  6 nF */
  gNa = 1.25;   /* Sodium conductance uS */
  PCa = 7.5;   /* Ca permeability x F */
  iKm = 20.0;   /* maximum iK  nA */
  gK1 = 0.75;   /* maximum gK1  uS */
  Camode = 5;   /* Ca induced inactivation set */
  CACT[3] = 0.0005;   /* binding constant 0.5 uM */
  gbNa = 0.07;   /* background sodium conductance uS */
  gbK = 0.0;   /* background potassium conductance */
  gfK = 6.0;   /* K conductance for i f  6 uS */
  gfNa = 6.0;   /* Na conductance for i f  6 uS */
  gTO = 0.0;   /* transient outward current zero */
  gbCa = 0.01;   /* background calcium conductance uS */
  Pump = 50.0;   /* maximum pump current nA */
  kNaCa = 0.002;   /* Na-Ca exchange scaling factor */
  DNaCa = 0.0001;   /* exchange denominator factor */
  _Y[10] = 140.0;   /* [K]i = 140 mM */
  _Y[1] = 7.5;   /* [Na]i = 7.5 mM */
  _Y[3] = 0.000058;   /* [Ca]i = 58 nM */
  Kb = 3.0;   /* bulk [K]o 3mM */
  _Y[2]/*=K*/ = 3.0;   /* [K]o = 3 mM */
  Nao = 140.0;   /* [Na]o = 140 mM */
  _Y[23] = /*Cao =*/ 2.0;   /* [Ca]o = 2 mM */
  _Y[4] = 0.007;   /* i f activation */
  _Y[5] = 0.54;   /* i K activation */
  _Y[6] = 0.0011;   /* i Ca activation */
  _Y[7] = 0.785;   /* i Ca inactivation */
  _Y[8] = 0.015;   /* i Na inactivation */
  _Y[9] = 0.076;   /* i Na activation */
  _Y[14] = _Y[7];   /* Ca-induced inactivation */
  Shift[6] = 5.0;   /* Voltage shift of iCa activation */
  TOmode = 1;   /* i TO equations include inactivation */
  dx = 8.0;   /* radial segment width um, 10 x 8 = 80 um */
  PrepLength = 0.08;   /* fibre length = 80 um */
  vol = 0.1;   /* fractional extracellular space 10 % */
  PF = 1.0;   /* rate constant for external diffusion sec^-1 */
  Tend = 1.25;   /* terminate computation at 1.25 sec */
  dt = 0.0025;   /* main step length 2.5 msec */
  TABT = 0.005;   /* main tabulation interval 5 msec */
  TP = 0.0;   /* pulse time zero, therefore no pulse */
  dVtest = 6000.0;   /* switch to full activation for i Na at 6000V/s */
  /* this is large enough never to be activated */
  EC = -60.0;   /* initial voltage -60 mV */
  StimSize = 0.0;   /* no stimuli necessary for spontaneous prep */
  _Y[11] = 1.98;   /* [Ca]sr uptake 1.98 mM */
  _Y[12] = 0.55;   /* [Ca]sr release 0.55 mM */
  _Y[13] = _Y[7];   /* repriming */
  Tau12 = 0.005;   /* SR uptake constant */
  Tau13 = 0.2;   /* SR transfer constant */
  TauRel = 0.01;   /* SR release constant */
  KmCa = 0.002;   /* Ca binding for release 2 uM */
  iNaCam = 10000.0;   /* maximum i NaCa 10,000 nA */
  Prep = 3;   /* preparation sinus */
  Kmode = 2;   /* iK equations: Oxford sinus iK */
  Speed[4] = 2.0;   /* i f activation speed x 2 */
  Kminact = 0.0005;   /* Km for Ca inactivation by Ca 0.5 uM */
  Shift[13] = -30.0;   /* shift repriming curve 30 mV negative */



}  /*of procedure SINUS =================================================*/




void SINSIN(void) {

  /*=========================================================================*/

  /*  set standard parameters for single sinus node as at
12 August 1989. Must be used after first calling SINUS

This is the first of the single sinus node cell models from
which others have been developed. See HEART.UPD for information
on more recent models.

calls: no other procedures

=============================================================================*/
  Space = 4;   /* single cell: no restricted space */
  EC = -69.1865;   /* initial voltage mV */
  _Y[3] = 0.000056;   /* [Ca]i 56 nM */
  _Y[4] = 0.0822;   /* i f activation */
  _Y[5] = 0.1231;   /* i K activation */
  _Y[6] = 0.0000;   /* i Ca activation */
  _Y[7] = 0.9997;   /* i Ca inactivation */
  _Y[8] = 0.1969;   /* i Na inactivation */
  _Y[9] = 0.0365;   /* i Na activation */
  _Y[11] = 2.3909;   /* [Ca]sr uptake 2.3909 mM */
  _Y[12] = 0.2207;   /* [Ca]sr release 0.2207 mM */
  _Y[13] = 0.2370;   /* repriming */
  _Y[14] = 0.5765;   /* Ca inactivation of iCa */
  _Y[15] = 1.0000;   /* i Ca3 activation -- but not used in this model!*/
  _Y[17] = 0.1028;   /* i to inactivation -- not used */
  _Y[18] = 1.0000;   /* i to activation -- not used */
  gNa = 0.0125;   /* sodium conductance uS */
  PCa = 0.12;   /* Ca permeability x F */
  iKm = 0.8;   /* maximum iK 800 pA */
  gK1 = 0.0075;   /* maximum gK1 uS */
  gbNa = 0.0007;   /* background sodium conductance 700 pS */
  gfK = 0.06;   /* K conductance of i f 60 nS */
  gfNa = 0.06;   /* Na conductance of i f 60 nS */
  gbCa = 0.0001;   /* background calcium conductance 100 pS */
  Pump = 0.45;   /* maximum pump current 450 pA */
  kNaCa = 0.00002;   /* i NaCa scaling factor */
  PrepLength = 0.028;   /* cell length 28 um */
  dx = 1.0;   /* cell radius 10 x dx = 10 um */
  C = 0.00006;   /* capacitance 60 pF */
  Ymode = 20;   /* i f kinetics use single cell data, squared */
  Kmode = 4;   /* i K kinetics use single cell data */
  TimeScale = 1000.0;   /* convert tabulation time scale to msec */
  Iscal = 1000;   /* convert tabulation current scale to pA */
  Tend = 1.0;   /* terminate computation at 1 sec */



}  /*of procedure SINSIN ==================================================*/




/*============================================================================*/

void HATRIUM(void) {

  /* sets Hilgemann atrium model as at 6 July 1986

Calls:

=============================================================================*/
  /* This part is a repeat of relevant parts of SINUS */

  C = 0.006;   /* capacitance 6 nF */
  CACT[3] = 0.0005;   /* Km for Ca activation of ito 500 nM */
  DNaCa = 0.0001;   /* denominator factor for sodium-calcium exchange */
  _Y[10] = 140.0;   /* [K]i = 140 mM */
  Nao = 140.0;   /* [Na]o = 140 mM */
  _Y[23] = /*Cao =*/ 2.0;   /* [Ca]o = 2 mM */
  _Y[6] = 0.0011;   /* iCa activation gate */
  Shift[6] = 5.0;   /* shift iCa activation 5 mV */
  _Y[7] = 0.785;   /* iCa inactivation gate */
  _Y[8] = 0.015;   /* iNa inactivation gate */
  _Y[9] = 0.076;   /* iNa activation gate */
  _Y[14] = _Y[7];   /* slow calcium inactivation */
  dx = 8.0;   /* radius 10 x DX = 80 um */
  PrepLength = 0.08;   /* length 80 um */
  PF = 1.0;   /* rate constant for diffusion to bulk space sec^-1 */
  Tend = 1.25;   /* end computation at 1.25 sec */
  TP = 0.0;   /* pulsetime at zero */
  StimSize = 0.0;   /* pulse set to zero */

  Prep = 3;
  Kminact = 0.0005;   /* Km for calcium inactivation of iCa 0.5 uM */

  CaOmode = 1;   /* calculate external calcium */
  Camode = 8;   /* use Hilgemann formulation for iCa */
  Bmode = 2;   /* use Hilgemann buffer equations */
  SLPumpMode = false;   /* no SL calcium pump */
  ContractMode = true;   /* cross-bridge formation calculated */
  SRmode = 2;   /* use Hilgemann SR equations */
  dt = 0.001;   /* integration steps 1 msec */
  TABT = 0.002;   /* tabulation steps 2 msec */
  out = 13;   /* output data includes contraction */
  gTO = 0.0;   /* no ito in this model */
  iKm = 0.0;   /* no iK in this model */
  Kb = 4.0;   /* [K]b = 4 mM */
  Pump = 14.0;   /* maximum pump current 14 nA */
  gNa = 50.0;   /* sodium conductance 50 uS */
  gfNa = 0.0;   /* no i f in this model */
  gbNa = 0.012;   /* background sodium conductance 12 nS */
  gbCa = 0.005;   /* background calcium conductance 5 nS */
  ShiftK1 = 20.0;   /* iK1 rectifier shifted 20 mV */
  gK1 = 1.7;   /* maximum gK1 1.7 uS */
  gbK = 0.17;   /* background K conductance 170 nS */
  _Y[0] = -88.0;   /* initial voltage -88 mV */
  EC = -88.0;   /* initial voltage -88 mV */
  Space = 4;   /* no external K gradients */
  _Y[2]/*=K*/ = 4.0;   /* [K]o = 4 mM */
  V12 = 0.01;   /* SR volume fraction 1 % */
  V13 = 0.1;   /* release volume fraction 10 % */
  vol = 0.4;   /* external volume fraction 40 % */
  kNaCa = 0.01;   /* sodium calcium exchange scaling factor */
  PCa = 5.0;   /* calcium permeability x F */
  Speed[6] = 3.0;   /* calcium activation speed x 3 */
  _Y[3] = 0.00001;   /* [Ca]i 10 nM */
  _Y[1] = 6.5;   /* [Na]i 6.5 mM */
  PCaK = 0.002;   /* K permeability of Ca channel 2 % */
  KCaCHoff = 0.001;   /* Ca binding to Ca regulatory site for iCa */
  CTrop = 0.15;   /* [ctopronin] = 150 uM */
  MTrop = 0.02;   /* [mtroponin] = 20 uM */
  BufFast = false;   /* use simplified buffer equations */
  alpha[20] = 100000.0;
      /* on rate for Ca binding to ctroponin mM^-1 sec^-1 */
  beta[20] = 200.0;   /* off rate for Ca binding to ctroponin sec^-1 */
  _Y[20] = 0.0015;   /* fraction of ctroponin with Ca bound */
  alpha[19] = 100000.0;
      /* on rate for Ca binding to calmodulin mM^-1 sec^-1 */
  beta[19] = 50.0;   /* off rate for Ca binding to calmodulin sec^-1 */
  _Y[19] = 0.0005;   /* fraction of calmodulin with Ca bound */
  KmCa = 0.0005;   /* Km for calcium-induced release 500 nM */
  _Y[21] = 0.0;   /* fraction of ctroponin with Mg bound */
  _Y[22] = 0.0;   /* fraction of mtroponin with Mg bound */

  alpha[11] = 3.0;   /* rate constant for uptake to SR */
  beta[11] = 0.23;   /* rate constant for leak from SR */
  _Y[11] = 0.3;   /* [Ca]sr,uptake = 0.3 mM */
  alpha[12] = 50.0;   /* rate constant for transfer from SR,up to SR,rel */
  _Y[12] = 0.3;   /* [Ca]sr,release = 0.3 mM */
  alpha[28] = 12000.0;   /* rate constant for light chain conformation */
  beta[28] = 100.0;   /* rate constant for decay of l c conformation */
  _Y[28] = 0.0;   /* fraction of light chain conformation formed */
  alpha[29] = 60.0;   /* rate constant for formation of cross bridges */
  beta[29] = 25.0;   /* rate constant for breaking of cross bridges */
  _Y[29] = 0.0;   /* fraction of cross-bridges formed */
  KmCa2 = 250.0;   /* rate constant for release by open SR channels */
  KCyCa = 0.0003;   /* rate coefficient KcyCa for SR Ca pump */
  KSRCa = 0.5;   /* rate coefficient KsrCa for SR Ca pump */
  Kxcs = 0.4;   /* rate coefficient Kxcs for SR Ca pump */
  KMSLpump = 0.0002;   /* Km for Ca activation of SL Ca pump = 200 nM */
  KSLpump = 22.0;
      /* maximum rate of SL pump, expressed as a Ca current nA */



}  /* OF PROCEDURE HATRIUM =================================================*/



void SINAT(void) {

  /*============================================================================*/

  /* set standard parameters for SINGLE ATRIAL CELL AS AT 12 AUGUST 1989.

Must be used after calling HATRIUM

calls: no other procedures

=============================================================================*/
  iPulseSize = -1.3;   /* current pulse -1.3 nA */
  on = 0.05;   /* applied at 50 msec */
  off = 0.052;   /* lasting 2 msec */
  Rep = 0.3;   /* repeated every 300 msec */
  dVtest = 200.0;
  Nai = 6.48;   /* [Na]i = 6.48 mM */
  Space = 4;   /* single cell: no restricted space */
  TOmode = 5;   /* Hilgemann-Noble ito */
  CaOmode = 0;   /* external calcium constant */
  gTO = 0.01;   /* g to 10 nS */
  _Y[18] = 0.0;   /* i to activation gate */
  _Y[17] = 1.0;   /* i to inactivation gate */

  gNa = 0.5;   /* g Na 500 nS */
  PCa = 0.05;   /* PCa x F */
  iKm = 0.0;   /* no iK */
  gbNa = 0.00012;   /* background sodium conductance 120 pS */
  gfK = 0.0;   /* no i f */
  gfNa = 0.0;
  gbCa = 0.00005;   /* background calcium conductance 50 pS */
  Pump = 0.14;   /* maximum pump current 140 pA */
  kNaCa = 0.0001;   /* Na-Ca exchange scaling factor */
  gK1 = 0.017;   /* maximum gK1 17 nS */
  gbK = 0.0017;   /* background K conductance 1.7 nS */
  EC = -91.6;   /* initial voltage -91.6 mV */
  PrepLength = 0.08;   /* cell length 80 um */
  dx = 1.0;   /* cell radius 10 x 1 = 10 um */
  C = 0.00004;   /* capacitance 40 pF */
  TimeScale = 1000.0;   /* convert time in output to msec */
  Iscal = 1000;   /* convert current in output to pA */
  Tend = 0.29;   /* end computation at 290 msec */

}  /*of procedure sinat ---------------------------------------------------*/





/*============================================================================*/

void VENTRICLE(void) {

  /* sets standard parameters for guinea-pig ventricle as at 20 April 1991

Should be used only after calling HATRIUM and SINAT

calls: no other procedures

=============================================================================*/
  iPulseSize = -6.0;   /* current pulse -6 nA */
  on = 0.1;   /* applied at 100 msec */
  off = 0.102;   /* lasting 2 msec */
  dVtest = 200.0;
  Rep = 1.0;   /* repeated every second */
  out = 4;   /* tabulate standard ionic currents */

  CaOmode = 0;   /* No external calcium space */
  Speed[7] = 0.5;   /* slow iCa inactivation 50% */
  BufFast = true;   /* use full Hilgemann buffer equations */

  CTrop = 0.05;   /* [ctroponin] 50 uM to increase [Ca]i transient */
  alpha[11] = 0.4;   /* reduce SR uptake */
  beta[11] = 0.03;   /* and leak to give shape of calcium transient */
  alpha[25] = 0.0;   /* no voltage dependent release */
  Camode = 5;   /* original iCa equations with Ca inactivation */
  TOmode = 5;   /* Hilgemann-Noble ito */
  gTO = 0.01;   /* transient outward conductance 10 nS */
  gNa = 2.5;   /* sodium conductance 2.5 uS */
  PCa = 0.25;   /* calcium channel permeability x F */
  iKm = 1.0;   /* maximum iK 1 nA */
  Kmode = 0;   /* original iK equations */
  gbNa = 0.0006;   /* background sodium conductance 600 pS */
  gfK = 0.0;   /* no i f */
  gfNa = 0.0;
  gbCa = 0.00025;   /* background calcium conductance 250 pS */
  Pump = 0.7;   /* maximum pump current 700 pA */
  kNaCa = 0.0002;   /* Na-Ca exchange scaling parameter */
  DNaCa = 0.0;   /* Na-Ca exchange denominator factor */
  gK1 = 1.0;   /* maximum gK1 1 uS */
  gbK = 0.0006;   /* background K conductance 600 pS */
  PrepLength = 0.12;   /* cell length 120 um */
  dx = 1.5;   /* cell radius 10 x 1.5 = 15 um */
  C = 0.0002;   /* capacitance 200 pF */
  TimeScale = 1000.0;   /* convert time to msec */
  dt = 0.001;   /* integration steps 1 msec */
  TABT = 0.001;   /* tabulation interval 1 msec */
  Iscal = 1000;   /* current in pA */
  steepK1 = 2.1;   /* steepness of iK1 rectifier 2.1 in place of 2 */
  ShiftK1 = 11.0;   /* shift rectifier voltage dependence 11 mV */
  KmK1 = 14.0;   /* Km for [K]o activation of iK1 14 mM */
  Tend = 0.4;   /* end calculation at 400 msec */


  SPvol = 0.1;   /* subsarcolemmal space 10% of internal space */
  iCafract = 0.5;   /* 50% of iCa flowing into it */
  iNCfract = 0.5;   /* 50% of iNaCa flowing to and from it */
  SRfract = 0.2;   /* 20% of SR sensing [Ca]i in space */
  TauSPvol = 0.01;   /* time constant for exchange with cytosol 10 msec */

  _Y[0] = -93.0094;   /* initial voltage  mV */
  _Y[1] = 5.0027;   /* initial [Na]i mM */
  _Y[2] = 4.0;   /* initial [K]o mM */
  _Y[3] = 0.00000995;   /* initial [Ca]i mM */
  _Y[5] = 0.0010;   /* iK activation */
  _Y[6] = 0.0;   /* iCa activation */
  _Y[7] = 1.0;   /* iCa inactivation */
  _Y[8] = 0.9946;   /* iNa inactivation */
  _Y[9] = 0.0016;   /* iNa activation */
  _Y[10] = 139.9982;   /* [K]i mM */
  _Y[11] = 0.3024;   /* [Ca]sr,uptake mM */
  _Y[12] = 0.2967;   /* [Ca]sr,release mM */
  _Y[13] = 0.7850;   /* repriming ????? */
  _Y[14] = 0.8974;   /* slow calcium inactivation -- Ca dependent */
  _Y[17] = 0.9999;   /* ito inactivation */
  _Y[18] = 0.0;   /* ito activation */
  _Y[19] = 0.0005;   /* Ca bound to calmodulin */
  _Y[20] = 0.0003;   /* Ca bound to troponin */
  _Y[21] = 0.0;   /* Mg bound to troponin ???? */
  _Y[22] = 0.0;   /* Ca bound to m-troponin ???? */
  _Y[23] = 2.00091257;   /* [Ca]o ??????? */
  _Y[25] = 0.9539;   /* SR precursor */
  _Y[26] = 0.0042;   /* SR activator */
  _Y[27] = 0.0420;   /* SR product */
  _Y[28] = 0.0001;   /* light chain conformation */
  _Y[29] = 0.0006;   /* cross bridge formation */
  _Y[30] = 0.0007;   /* total cytosol calcium */
  _Y[34] = 4.997660;   /* ATP ????? */
  _Y[36] = 0.000003;   /* [Ca]space */
  _Y[37] = 0.299014;   /* [Ca]sr,up,space */
  _Y[38] = 0.299825;   /* [Ca]sr,rel,space */
  _Y[39] = 0.000129;   /* buffer 1 */
  _Y[40] = 0.000074;   /* buffer 2 */
  _Y[41] = 5.000000;   /* [Na]space */
  _Y[42] = 0.785000;   /* Ca inactivation of iCa in space */
  _Y[43] = 0.0;   /* calcium probe */
  _Y[44] = 0.000206;   /* total cytosol calcium in space */


}  /*of procedure ventricle =============================================*/






void RATVCELL(void) {

  /*=========================================================================*/

  /*  sets standard parameters for rat ventricular cell at 4 July 1985

calls: no other procedures

==========================================================================*/
  C = 0.0002;   /* capacitance 200 pF */
  gNa = 4.0;   /* gNa 4 uS */
  PCa = 0.4;   /* PCa x F */
  iKm = 1.0;   /* maximum iK 1 nA */
  gK1 = 0.7;   /* maximum gK1 0.7 uS */
  gbNa = 0.005;   /* background sodium conductance 5 nS */
  gbK = 0.03;   /* background K conductance 30 nS */
  gfK = 0.0;   /* no i f */
  gfNa = 0.0;
  gTO = 0.01;   /* g to 10 nS */
  gbCa = 0.0025;   /* background calcium conductance 2.5 nS */
  Pump = 4.0;   /* maximum pump current 4 nA */
  kNaCa = 0.0005;   /* iNaCa scaling factor */
  _Y[10] = 140.0;   /* [K]i 140 mM */
  _Y[1] = 10.0;   /* [Na]i 10 mM */
  _Y[3] = 0.00005;   /* [Ca]i 50 nM */
  Kb = 3.8;   /* [K]o 3.8 mM */
  _Y[2]/*=K*/ = 3.8;
  Nao = 135.0;   /* [Na]o 135 mM */
  _Y[23] = /*Cao =*/ 2.5;   /* [Ca]o 2.5 mM */
  _Y[5] = 0.01;   /* iK activation */
  _Y[6] = 0.005;   /* iCa activation */
  _Y[7] = 1.0;   /* iCa inactivation */
  _Y[14] = _Y[7];   /* slow Ca inactivation */
  _Y[8] = 0.8;   /* i Na inactivation */
  _Y[9] = 0.01;   /* i Na activation */
  Shift[6] = 5.0;   /* shift iCa activation 5 mV */
  TOmode = 1;   /* include ito inactivation */
  dx = 1.0;   /* radius of cell 10 x 1 = 10 um */
  PrepLength = 0.11;   /* cell length 110 um */
  vol = 0.1;   /* external volume 10 % ????? */
  PF = 0.7;   /* ?????? */
  Tend = 0.25;   /* end calculation at 250 msec */
  dt = 0.002;   /* integration steps 2 msec */
  TABT = 0.002;   /* tabulation 2 msec */
  dVtest = 200.0;   /* test for m equations 200 V/sec */
  EC = -80.0;   /* initial voltage mV */
  TPoff = Tend + 1;
  Speed[7] = 3.5;   /* speed calcium inactivation */
  _Y[11] = 1.6;   /* [Ca]sr,uptake mM */
  _Y[12] = 1.4;   /* [Ca]sr,release mM */
  _Y[13] = _Y[7];   /* SR repriming */
  Tau12 = 0.01;   /* uptake time constant */
  Tau13 = 0.1;   /* transfer time constant */
  TauRel = 0.05;   /* release time constant */
  KmCa = 0.005;   /* Km for Ca dependente Ca release 5 uM */
  Prep = 4;   /* Rat cell */
  Space = 4;   /* single cell space */
  out = 4;   /* standard currents in output file */

  TimeScale = 1000.0;   /* time in msec */
  Iscal = 1000;   /* current in pA */
  iPulseSize = -5.0;   /* current pulse -5 nA */
  on = 0.05;   /* applied at 50 msec */
  off = 0.052;   /* lasting 2 msec */



}  /* of procedure RATVCELL =============================================*/







/*=========================================================================*/

void GPCELL(void) {

  /* sets standard parameters for guinea pig cell - April 1991

This is the model described by Noble et al (1991)
Annals of the New York Academy of Sciences.
Proceedings of 2nd International Meeting on Sodium-Calcium exchange

The major deficiency of this preliminary model is that it has
a notch at the beginning of the plateau. This is seen in the
ventricle of some species, but not usually the guinea-pig.

In this version, you can eliminate this problem if you wish
by setting SPMIMIC = 1. This speeds the calcium release.


must be used after calling Hatrium and Sinat

==========================================================================*/
  strcpy(model, "Simplified GP cell model");

  iPulseSize = -6.0;   /* 6 nA stimulus */
  on = 0.1;   /* at 100 msec */
  off = 0.102;   /* lasting 2 msec */
  Rep = 1.0;   /* repeat after 1 sec */
  dVtest = 200.0;
  _Y[1] = 5.0;   /* [Na]i = 5 mM */

  Speed[7] = 0.5;   /* slow iCa inactivation */
  BufFast = true;   /* use full buffer equations */

  CTrop = 0.05;   /* ctroponin 50 uM */
  alpha[11] = 0.4;   /* rate constant for uptake to SR */
  beta[11] = 0.03;   /* rate constant for leak from SR */
  alpha[25] = 0.0;   /* no voltage dependence of Ca release */
  Camode = 1;   /* 'Oxford' iCa equations */
  TOmode = 5;   /* Hilgemann-Noble i TO */
  gTO = 0.005;   /* g to 5 nS */
  gNa = 2.5;   /* g Na 2.5 uS */
  PCa = 0.25;   /* P Ca */
  iKm = 1.0;   /* maximum iK 1 nA */
  Kmode = 0;   /* MNT K equations */
  gbNa = 0.0006;   /* background sodium conductance 600 pS */
  gfK = 0.0;   /* no i f */
  gfNa = 0.0;   /* no i f */
  gbCa = 0.00025;   /* background calcium conductance 250 pS */
  Pump = 0.7;   /* maximum pump current 700 pA */
  kNaCa = 0.0005;   /* sodium-calcium exchange scaling factor */
  DNaCa = 0.0;   /* sodium-calcium exchange denominator factor */
  gK1 = 1.0;   /* maximum gK1 1 uS */
  gbK = 0.0006;   /* background K conductance 600 pS */
  PrepLength = 0.08;   /* cell length 80 um */
  dx = 1.5;   /* cell radius = 10 x DX = 15 um */
  C = 0.0002;   /* capacitance 200 pF */
  TimeScale = 1000.0;   /* convert time scale to msec in output file */
  Iscal = 1000;   /* convert current scale to pA in output file */
  Tend = 0.4;   /* end computation at 400 msec */

  EC = -93.0;   /* initial voltage -93 mV */
  Kb = 4.0;
  _Y[2]/*=K*/ = 4.0;   /* [K]o 4 mM */
  _Y[3] = 0.00000763;   /* [Ca]i mM */
  _Y[5] = 0.001;   /* initial i K activation */
  _Y[6] = 0.0;   /* initial i Ca activation */
  _Y[7] = 1.0;   /* initial i Ca inactivation */
  _Y[8] = 0.9951;   /* initial i Na inactivation */
  _Y[9] = 0.0015;   /* initial i Na activation */
  _Y[10] = 140.0;   /* [K]i */
  _Y[11] = 0.3013;   /* SR store calcium mM */
  _Y[12] = 0.2989;   /* SR release calcium mM */
  _Y[13] = 0.7850;   /* repriming variable */
  _Y[17] = 1.0;   /* i TO inactivation gate */
  _Y[18] = 0.0;   /* i TO activation gate */
  _Y[19] = 0.0003;   /* calmodulin sites occupied */
  _Y[20] = 0.0002;   /* troponin sites occupied */
  _Y[21] = 0.0;   /* Mg bound to troponin */
  _Y[22] = 0.0;   /* Ca bound to m troponin */
  _Y[23] = 2.0;   /* extracellular calcium mM */
  _Y[25] = 0.9687;   /* precursor fraction */
  _Y[26] = 0.0023;   /* activator fraction */
  _Y[27] = 0.0291;   /* product fraction */
  _Y[28] = 0.0;   /* light chain conformation */
  _Y[29] = 0.0003;   /* cross bridge formation */
  _Y[30] = 0.0005;   /* total cytosol calcium */
  _Y[34] = 5.0;   /* ATP concentration mM */

}  /*of procedure GPcell */




/*============================================================================*/

void FROGSINUS(void) {

  /* sets standard parameters for single frog sinus venosus based on
                Calgary equations as at 19 July 1984

calls: no other procedures

===========================================================================*/
  printf("N.B. This parameter set is available only on the Oxford and\n");
  printf("Calgary versions of the program\n");
  printf("The model is not yet finished\n\n\n");

}  /* of procedure FROGSINUS ===============================================*/





void ATRIUM(void) {

  /*           Sets standard parameters for frog atrium as at 18 July 1984
Calls: no other procedures

===========================================================================*/
  printf("Frog atrium - 18-july-1984\n\n");
  printf("This parameter set is available only in Calgary and Oxford versions\n\n");

}  /* of procedure ATRIUM ==================================================*/
