/*
This is the CellML export to C of the OHara Rudy code. The aim of this exercise
is to put the cell model into Beatbox format with minimal effort and modification.
Along the way, the way in which it is done should also be documented so that
it gives rise to a MSc project in Liverpool or Cesena.
*/
//==============================================================================
// CellML file:   C:\Users\sanjay\Desktop\Ohara_Rudy_2011.cellml
// CellML model:  MyModel1
// Date and time: 13/12/2011 at 12:58:02
//------------------------------------------------------------------------------
// Conversion from CellML 1.0 to C was done using COR (0.9.31.1409)
//    Copyright 2002-2011 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================

#include "Ohara_Rudy_2011.h"

//------------------------------------------------------------------------------

#include <math.h>
#include <string.h>

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

double Y[_NB_OF_STATE_VARIABLES_];
double dY[_NB_OF_STATE_VARIABLES_];
// 0: CaMKt (millimolar) (in CaMK)
// 1: d (dimensionless) (in ICaL)
// 2: fcaf (dimensionless) (in ICaL)
// 3: fcafp (dimensionless) (in ICaL)
// 4: fcas (dimensionless) (in ICaL)
// 5: ff (dimensionless) (in ICaL)
// 6: ffp (dimensionless) (in ICaL)
// 7: fs (dimensionless) (in ICaL)
// 8: jca (dimensionless) (in ICaL)
// 9: nca (dimensionless) (in ICaL)
// 10: xk1 (dimensionless) (in IK1)
// 11: xrf (dimensionless) (in IKr)
// 12: xrs (dimensionless) (in IKr)
// 13: xs1 (dimensionless) (in IKs)
// 14: xs2 (dimensionless) (in IKs)
// 15: hL (dimensionless) (in INaL)
// 16: hLp (dimensionless) (in INaL)
// 17: mL (dimensionless) (in INaL)
// 18: hf (dimensionless) (in INa)
// 19: hs (dimensionless) (in INa)
// 20: hsp (dimensionless) (in INa)
// 21: j (dimensionless) (in INa)
// 22: jp (dimensionless) (in INa)
// 23: m (dimensionless) (in INa)
// 24: a (dimensionless) (in Ito)
// 25: ap (dimensionless) (in Ito)
// 26: iF (dimensionless) (in Ito)
// 27: iFp (dimensionless) (in Ito)
// 28: iS (dimensionless) (in Ito)
// 29: iSp (dimensionless) (in Ito)
// 30: cai (millimolar) (in intracellular_ions)
// 31: cajsr (millimolar) (in intracellular_ions)
// 32: cansr (millimolar) (in intracellular_ions)
// 33: cass (millimolar) (in intracellular_ions)
// 34: ki (millimolar) (in intracellular_ions)
// 35: kss (millimolar) (in intracellular_ions)
// 36: nai (millimolar) (in intracellular_ions)
// 37: nass (millimolar) (in intracellular_ions)
// 38: v (millivolt) (in membrane)
// 39: Jrelnp (dimensionless) (in ryr)
// 40: Jrelp (dimensionless) (in ryr)

char YNames[_NB_OF_STATE_VARIABLES_][7];
char YUnits[_NB_OF_STATE_VARIABLES_][14];
char YComponents[_NB_OF_STATE_VARIABLES_][19];

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

double CaMKo;   // dimensionless (in CaMK)
double KmCaM;   // millimolar (in CaMK)
double KmCaMK;   // millimolar (in CaMK)
double aCaMK;   // per_millimolar_per_millisecond (in CaMK)
double bCaMK;   // per_millisecond (in CaMK)
double Kmn;   // millimolar (in ICaL)
double PCa_b;   // dimensionless (in ICaL)
double k2n;   // per_millisecond (in ICaL)
double PCab;   // milliS_per_microF (in ICab)
double GK1_b;   // milliS_per_microF (in IK1)
double GKb_b;   // milliS_per_microF (in IKb)
double GKr_b;   // milliS_per_microF (in IKr)
double GKs_b;   // milliS_per_microF (in IKs)
double Gncx_b;   // milliS_per_microF (in INaCa_i)
double KmCaAct;   // millimolar (in INaCa_i)
double kasymm;   // dimensionless (in INaCa_i)
double kcaoff;   // per_millisecond (in INaCa_i)
double kcaon;   // per_millisecond (in INaCa_i)
double kna1;   // per_millisecond (in INaCa_i)
double kna2;   // per_millisecond (in INaCa_i)
double kna3;   // per_millisecond (in INaCa_i)
double qca;   // dimensionless (in INaCa_i)
double qna;   // dimensionless (in INaCa_i)
double wca;   // dimensionless (in INaCa_i)
double wna;   // dimensionless (in INaCa_i)
double wnaca;   // dimensionless (in INaCa_i)
double H;   // millimolar (in INaK)
double Khp;   // millimolar (in INaK)
double Kki;   // per_millisecond (in INaK)
double Kko;   // per_millisecond (in INaK)
double Kmgatp;   // millimolar (in INaK)
double Knai0;   // millimolar (in INaK)
double Knao0;   // millimolar (in INaK)
double Knap;   // millimolar (in INaK)
double Kxkur;   // millimolar (in INaK)
double MgADP;   // millimolar (in INaK)
double MgATP;   // millimolar (in INaK)
double Pnak_b;   // milliS_per_microF (in INaK)
double delta;   // millivolt (in INaK)
double eP;   // dimensionless (in INaK)
double k1m;   // per_millisecond (in INaK)
double k1p;   // per_millisecond (in INaK)
double k2m;   // per_millisecond (in INaK)
double k2p;   // per_millisecond (in INaK)
double k3m;   // per_millisecond (in INaK)
double k3p;   // per_millisecond (in INaK)
double k4m;   // per_millisecond (in INaK)
double k4p;   // per_millisecond (in INaK)
double GNaL_b;   // milliS_per_microF (in INaL)
double thL;   // millisecond (in INaL)
double PNab;   // milliS_per_microF (in INab)
double Ahf;   // dimensionless (in INa)
double GNa;   // milliS_per_microF (in INa)
double hssV1;   // millivolt (in INa)
double hssV2;   // millivolt (in INa)
double mssV1;   // millivolt (in INa)
double mssV2;   // millivolt (in INa)
double mtD1;   // dimensionless (in INa)
double mtD2;   // dimensionless (in INa)
double mtV1;   // millivolt (in INa)
double mtV2;   // millivolt (in INa)
double mtV3;   // millivolt (in INa)
double mtV4;   // millivolt (in INa)
double GpCa;   // milliS_per_microF (in IpCa)
double KmCap;   // millimolar (in IpCa)
double Gto_b;   // milliS_per_microF (in Ito)
double L;   // centimeter (in cell_geometry)
double rad;   // centimeter (in cell_geometry)
double celltype;   // dimensionless (in environment)
double cao;   // millimolar (in extracellular)
double ko;   // millimolar (in extracellular)
double nao;   // millimolar (in extracellular)
double BSLmax;   // millimolar (in intracellular_ions)
double BSRmax;   // millimolar (in intracellular_ions)
double KmBSL;   // millimolar (in intracellular_ions)
double KmBSR;   // millimolar (in intracellular_ions)
double cm;   // microF_per_centimeter_squared (in intracellular_ions)
double cmdnmax_b;   // millimolar (in intracellular_ions)
double csqnmax;   // millimolar (in intracellular_ions)
double kmcmdn;   // millimolar (in intracellular_ions)
double kmcsqn;   // millimolar (in intracellular_ions)
double kmtrpn;   // millimolar (in intracellular_ions)
double trpnmax;   // millimolar (in intracellular_ions)
double amp;   // microA_per_microF (in membrane)
double duration;   // millisecond (in membrane)
double F;   // coulomb_per_mole (in physical_constants)
double R;   // joule_per_kilomole_kelvin (in physical_constants)
double T;   // kelvin (in physical_constants)
double zca;   // dimensionless (in physical_constants)
double zk;   // dimensionless (in physical_constants)
double zna;   // dimensionless (in physical_constants)
double PKNa;   // dimensionless (in reversal_potentials)
double bt;   // millisecond (in ryr)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

double CaMKa;   // millimolar (in CaMK)
double CaMKb;   // millimolar (in CaMK)
double Afcaf;   // dimensionless (in ICaL)
double Afcas;   // dimensionless (in ICaL)
double Aff;   // dimensionless (in ICaL)
double Afs;   // dimensionless (in ICaL)
double ICaK;   // microA_per_microF (in ICaL)
double ICaL;   // microA_per_microF (in ICaL)
double ICaNa;   // microA_per_microF (in ICaL)
double PCa;   // dimensionless (in ICaL)
double PCaK;   // dimensionless (in ICaL)
double PCaKp;   // dimensionless (in ICaL)
double PCaNa;   // dimensionless (in ICaL)
double PCaNap;   // dimensionless (in ICaL)
double PCap;   // dimensionless (in ICaL)
double PhiCaK;   // dimensionless (in ICaL)
double PhiCaL;   // dimensionless (in ICaL)
double PhiCaNa;   // dimensionless (in ICaL)
double anca;   // dimensionless (in ICaL)
double dss;   // dimensionless (in ICaL)
double f;   // dimensionless (in ICaL)
double fICaLp;   // dimensionless (in ICaL)
double fca;   // dimensionless (in ICaL)
double fcap;   // dimensionless (in ICaL)
double fcass;   // dimensionless (in ICaL)
double fp;   // dimensionless (in ICaL)
double fss;   // dimensionless (in ICaL)
double km2n;   // per_millisecond (in ICaL)
double td;   // millisecond (in ICaL)
double tfcaf;   // millisecond (in ICaL)
double tfcafp;   // millisecond (in ICaL)
double tfcas;   // millisecond (in ICaL)
double tff;   // millisecond (in ICaL)
double tffp;   // millisecond (in ICaL)
double tfs;   // millisecond (in ICaL)
double tjca;   // millisecond (in ICaL)
double ICab;   // microA_per_microF (in ICab)
double GK1;   // milliS_per_microF (in IK1)
double IK1;   // microA_per_microF (in IK1)
double rk1;   // millisecond (in IK1)
double txk1;   // millisecond (in IK1)
double xk1ss;   // dimensionless (in IK1)
double GKb;   // milliS_per_microF (in IKb)
double IKb;   // microA_per_microF (in IKb)
double xkb;   // dimensionless (in IKb)
double Axrf;   // dimensionless (in IKr)
double Axrs;   // dimensionless (in IKr)
double GKr;   // milliS_per_microF (in IKr)
double IKr;   // microA_per_microF (in IKr)
double rkr;   // dimensionless (in IKr)
double txrf;   // millisecond (in IKr)
double txrs;   // millisecond (in IKr)
double xr;   // dimensionless (in IKr)
double xrss;   // dimensionless (in IKr)
double GKs;   // milliS_per_microF (in IKs)
double IKs;   // microA_per_microF (in IKs)
double KsCa;   // dimensionless (in IKs)
double txs1;   // millisecond (in IKs)
double txs2;   // millisecond (in IKs)
double xs1ss;   // dimensionless (in IKs)
double xs2ss;   // dimensionless (in IKs)
double E1_i;   // dimensionless (in INaCa_i)
double E1_ss;   // dimensionless (in INaCa_i)
double E2_i;   // dimensionless (in INaCa_i)
double E2_ss;   // dimensionless (in INaCa_i)
double E3_i;   // dimensionless (in INaCa_i)
double E3_ss;   // dimensionless (in INaCa_i)
double E4_i;   // dimensionless (in INaCa_i)
double E4_ss;   // dimensionless (in INaCa_i)
double Gncx;   // milliS_per_microF (in INaCa_i)
double INaCa_i;   // microA_per_microF (in INaCa_i)
double INaCa_ss;   // microA_per_microF (in INaCa_i)
double JncxCa_i;   // millimolar_per_millisecond (in INaCa_i)
double JncxCa_ss;   // millimolar_per_millisecond (in INaCa_i)
double JncxNa_i;   // millimolar_per_millisecond (in INaCa_i)
double JncxNa_ss;   // millimolar_per_millisecond (in INaCa_i)
double allo_i;   // dimensionless (in INaCa_i)
double allo_ss;   // dimensionless (in INaCa_i)
double h10_i;   // dimensionless (in INaCa_i)
double h10_ss;   // dimensionless (in INaCa_i)
double h11_i;   // dimensionless (in INaCa_i)
double h11_ss;   // dimensionless (in INaCa_i)
double h12_i;   // dimensionless (in INaCa_i)
double h12_ss;   // dimensionless (in INaCa_i)
double h1_i;   // dimensionless (in INaCa_i)
double h1_ss;   // dimensionless (in INaCa_i)
double h2_i;   // dimensionless (in INaCa_i)
double h2_ss;   // dimensionless (in INaCa_i)
double h3_i;   // dimensionless (in INaCa_i)
double h3_ss;   // dimensionless (in INaCa_i)
double h4_i;   // dimensionless (in INaCa_i)
double h4_ss;   // dimensionless (in INaCa_i)
double h5_i;   // dimensionless (in INaCa_i)
double h5_ss;   // dimensionless (in INaCa_i)
double h6_i;   // dimensionless (in INaCa_i)
double h6_ss;   // dimensionless (in INaCa_i)
double h7_i;   // dimensionless (in INaCa_i)
double h7_ss;   // dimensionless (in INaCa_i)
double h8_i;   // dimensionless (in INaCa_i)
double h8_ss;   // dimensionless (in INaCa_i)
double h9_i;   // dimensionless (in INaCa_i)
double h9_ss;   // dimensionless (in INaCa_i)
double hca;   // dimensionless (in INaCa_i)
double hna;   // dimensionless (in INaCa_i)
double k1_i;   // dimensionless (in INaCa_i)
double k1_ss;   // dimensionless (in INaCa_i)
double k2_i;   // dimensionless (in INaCa_i)
double k2_ss;   // dimensionless (in INaCa_i)
double k3_i;   // dimensionless (in INaCa_i)
double k3_ss;   // dimensionless (in INaCa_i)
double k3p_i;   // dimensionless (in INaCa_i)
double k3p_ss;   // dimensionless (in INaCa_i)
double k3pp_i;   // dimensionless (in INaCa_i)
double k3pp_ss;   // dimensionless (in INaCa_i)
double k4_i;   // dimensionless (in INaCa_i)
double k4_ss;   // dimensionless (in INaCa_i)
double k4p_i;   // dimensionless (in INaCa_i)
double k4p_ss;   // dimensionless (in INaCa_i)
double k4pp_i;   // dimensionless (in INaCa_i)
double k4pp_ss;   // dimensionless (in INaCa_i)
double k5_i;   // dimensionless (in INaCa_i)
double k5_ss;   // dimensionless (in INaCa_i)
double k6_i;   // dimensionless (in INaCa_i)
double k6_ss;   // dimensionless (in INaCa_i)
double k7_i;   // dimensionless (in INaCa_i)
double k7_ss;   // dimensionless (in INaCa_i)
double k8_i;   // dimensionless (in INaCa_i)
double k8_ss;   // dimensionless (in INaCa_i)
double x1_i;   // dimensionless (in INaCa_i)
double x1_ss;   // dimensionless (in INaCa_i)
double x2_i;   // dimensionless (in INaCa_i)
double x2_ss;   // dimensionless (in INaCa_i)
double x3_i;   // dimensionless (in INaCa_i)
double x3_ss;   // dimensionless (in INaCa_i)
double x4_i;   // dimensionless (in INaCa_i)
double x4_ss;   // dimensionless (in INaCa_i)
double E1;   // dimensionless (in INaK)
double E2;   // dimensionless (in INaK)
double E3;   // dimensionless (in INaK)
double E4;   // dimensionless (in INaK)
double INaK;   // microA_per_microF (in INaK)
double JnakK;   // millimolar_per_millisecond (in INaK)
double JnakNa;   // millimolar_per_millisecond (in INaK)
double Knai;   // millimolar (in INaK)
double Knao;   // millimolar (in INaK)
double P;   // dimensionless (in INaK)
double Pnak;   // milliS_per_microF (in INaK)
double a1;   // dimensionless (in INaK)
double a2;   // dimensionless (in INaK)
double a3;   // dimensionless (in INaK)
double a4;   // dimensionless (in INaK)
double b1;   // dimensionless (in INaK)
double b2;   // dimensionless (in INaK)
double b3;   // dimensionless (in INaK)
double b4;   // dimensionless (in INaK)
double x1;   // dimensionless (in INaK)
double x2;   // dimensionless (in INaK)
double x3;   // dimensionless (in INaK)
double x4;   // dimensionless (in INaK)
double GNaL;   // milliS_per_microF (in INaL)
double INaL;   // microA_per_microF (in INaL)
double fINaLp;   // dimensionless (in INaL)
double hLss;   // dimensionless (in INaL)
double hLssp;   // dimensionless (in INaL)
double mLss;   // dimensionless (in INaL)
double thLp;   // millisecond (in INaL)
double tmL;   // millisecond (in INaL)
double INab;   // microA_per_microF (in INab)
double Ahs;   // dimensionless (in INa)
double INa;   // microA_per_microF (in INa)
double fINap;   // dimensionless (in INa)
double h;   // dimensionless (in INa)
double hp;   // dimensionless (in INa)
double hss;   // dimensionless (in INa)
double hssp;   // dimensionless (in INa)
double jss;   // dimensionless (in INa)
double mss;   // dimensionless (in INa)
double thf;   // millisecond (in INa)
double ths;   // millisecond (in INa)
double thsp;   // millisecond (in INa)
double tj;   // millisecond (in INa)
double tjp;   // millisecond (in INa)
double tm;   // millisecond (in INa)
double IpCa;   // microA_per_microF (in IpCa)
double AiF;   // dimensionless (in Ito)
double AiS;   // dimensionless (in Ito)
double Gto;   // milliS_per_microF (in Ito)
double Ito;   // microA_per_microF (in Ito)
double ass;   // dimensionless (in Ito)
double assp;   // dimensionless (in Ito)
double delta_epi;   // dimensionless (in Ito)
double dti_develop;   // dimensionless (in Ito)
double dti_recover;   // dimensionless (in Ito)
double fItop;   // dimensionless (in Ito)
double i;   // dimensionless (in Ito)
double ip;   // dimensionless (in Ito)
double iss;   // dimensionless (in Ito)
double ta;   // millisecond (in Ito)
double tiF;   // millisecond (in Ito)
double tiF_b;   // millisecond (in Ito)
double tiFp;   // millisecond (in Ito)
double tiS;   // millisecond (in Ito)
double tiS_b;   // millisecond (in Ito)
double tiSp;   // millisecond (in Ito)
double Jleak;   // millimolar_per_millisecond (in SERCA)
double Jup;   // millimolar_per_millisecond (in SERCA)
double Jupnp;   // millimolar_per_millisecond (in SERCA)
double Jupp;   // millimolar_per_millisecond (in SERCA)
double fJupp;   // dimensionless (in SERCA)
double upScale;   // dimensionless (in SERCA)
double Acap;   // centimeter_squared (in cell_geometry)
double Ageo;   // centimeter_squared (in cell_geometry)
double vcell;   // microliter (in cell_geometry)
double vjsr;   // microliter (in cell_geometry)
double vmyo;   // microliter (in cell_geometry)
double vnsr;   // microliter (in cell_geometry)
double vss;   // microliter (in cell_geometry)
double Jdiff;   // millimolar_per_millisecond (in diff)
double JdiffK;   // millimolar_per_millisecond (in diff)
double JdiffNa;   // millimolar_per_millisecond (in diff)
double Bcai;   // dimensionless (in intracellular_ions)
double Bcajsr;   // dimensionless (in intracellular_ions)
double Bcass;   // dimensionless (in intracellular_ions)
double cmdnmax;   // millimolar (in intracellular_ions)
double Istim;   // microA_per_microF (in membrane)
double vffrt;   // coulomb_per_mole (in membrane)
double vfrt;   // dimensionless (in membrane)
double EK;   // millivolt (in reversal_potentials)
double EKs;   // millivolt (in reversal_potentials)
double ENa;   // millivolt (in reversal_potentials)
double Jrel;   // millimolar_per_millisecond (in ryr)
double Jrel_inf;   // dimensionless (in ryr)
double Jrel_inf_temp;   // dimensionless (in ryr)
double Jrel_infp;   // dimensionless (in ryr)
double Jrel_temp;   // dimensionless (in ryr)
double a_rel;   // millisecond (in ryr)
double a_relp;   // millisecond (in ryr)
double btp;   // millisecond (in ryr)
double fJrelp;   // dimensionless (in ryr)
double tau_rel;   // millisecond (in ryr)
double tau_rel_temp;   // millisecond (in ryr)
double tau_relp;   // millisecond (in ryr)
double tau_relp_temp;   // millisecond (in ryr)
double Jtr;   // millimolar_per_millisecond (in trans_flux)

//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void init()
{
   //---------------------------------------------------------------------------
   // State variables
   //---------------------------------------------------------------------------

   Y[0] = 0.0;   // CaMKt (millimolar) (in CaMK)
   Y[1] = 0.0;   // d (dimensionless) (in ICaL)
   Y[2] = 1.0;   // fcaf (dimensionless) (in ICaL)
   Y[3] = 1.0;   // fcafp (dimensionless) (in ICaL)
   Y[4] = 1.0;   // fcas (dimensionless) (in ICaL)
   Y[5] = 1.0;   // ff (dimensionless) (in ICaL)
   Y[6] = 1.0;   // ffp (dimensionless) (in ICaL)
   Y[7] = 1.0;   // fs (dimensionless) (in ICaL)
   Y[8] = 1.0;   // jca (dimensionless) (in ICaL)
   Y[9] = 0.0;   // nca (dimensionless) (in ICaL)
   Y[10] = 1.0;   // xk1 (dimensionless) (in IK1)
   Y[11] = 0.0;   // xrf (dimensionless) (in IKr)
   Y[12] = 0.0;   // xrs (dimensionless) (in IKr)
   Y[13] = 0.0;   // xs1 (dimensionless) (in IKs)
   Y[14] = 0.0;   // xs2 (dimensionless) (in IKs)
   Y[15] = 1.0;   // hL (dimensionless) (in INaL)
   Y[16] = 1.0;   // hLp (dimensionless) (in INaL)
   Y[17] = 0.0;   // mL (dimensionless) (in INaL)
   Y[18] = 1.0;   // hf (dimensionless) (in INa)
   Y[19] = 1.0;   // hs (dimensionless) (in INa)
   Y[20] = 1.0;   // hsp (dimensionless) (in INa)
   Y[21] = 1.0;   // j (dimensionless) (in INa)
   Y[22] = 1.0;   // jp (dimensionless) (in INa)
   Y[23] = 0.0;   // m (dimensionless) (in INa)
   Y[24] = 0.0;   // a (dimensionless) (in Ito)
   Y[25] = 0.0;   // ap (dimensionless) (in Ito)
   Y[26] = 1.0;   // iF (dimensionless) (in Ito)
   Y[27] = 1.0;   // iFp (dimensionless) (in Ito)
   Y[28] = 1.0;   // iS (dimensionless) (in Ito)
   Y[29] = 1.0;   // iSp (dimensionless) (in Ito)
   Y[30] = 1.0e-4;   // cai (millimolar) (in intracellular_ions)
   Y[31] = 1.2;   // cajsr (millimolar) (in intracellular_ions)
   Y[32] = 1.2;   // cansr (millimolar) (in intracellular_ions)
   Y[33] = 1.0e-4;   // cass (millimolar) (in intracellular_ions)
   Y[34] = 145.0;   // ki (millimolar) (in intracellular_ions)
   Y[35] = 145.0;   // kss (millimolar) (in intracellular_ions)
   Y[36] = 7.0;   // nai (millimolar) (in intracellular_ions)
   Y[37] = 7.0;   // nass (millimolar) (in intracellular_ions)
   Y[38] = -87.0;   // v (millivolt) (in membrane)
   Y[39] = 0.0;   // Jrelnp (dimensionless) (in ryr)
   Y[40] = 0.0;   // Jrelp (dimensionless) (in ryr)

   strcpy(YNames[0], "CaMKt");
   strcpy(YNames[1], "d");
   strcpy(YNames[2], "fcaf");
   strcpy(YNames[3], "fcafp");
   strcpy(YNames[4], "fcas");
   strcpy(YNames[5], "ff");
   strcpy(YNames[6], "ffp");
   strcpy(YNames[7], "fs");
   strcpy(YNames[8], "jca");
   strcpy(YNames[9], "nca");
   strcpy(YNames[10], "xk1");
   strcpy(YNames[11], "xrf");
   strcpy(YNames[12], "xrs");
   strcpy(YNames[13], "xs1");
   strcpy(YNames[14], "xs2");
   strcpy(YNames[15], "hL");
   strcpy(YNames[16], "hLp");
   strcpy(YNames[17], "mL");
   strcpy(YNames[18], "hf");
   strcpy(YNames[19], "hs");
   strcpy(YNames[20], "hsp");
   strcpy(YNames[21], "j");
   strcpy(YNames[22], "jp");
   strcpy(YNames[23], "m");
   strcpy(YNames[24], "a");
   strcpy(YNames[25], "ap");
   strcpy(YNames[26], "iF");
   strcpy(YNames[27], "iFp");
   strcpy(YNames[28], "iS");
   strcpy(YNames[29], "iSp");
   strcpy(YNames[30], "cai");
   strcpy(YNames[31], "cajsr");
   strcpy(YNames[32], "cansr");
   strcpy(YNames[33], "cass");
   strcpy(YNames[34], "ki");
   strcpy(YNames[35], "kss");
   strcpy(YNames[36], "nai");
   strcpy(YNames[37], "nass");
   strcpy(YNames[38], "v");
   strcpy(YNames[39], "Jrelnp");
   strcpy(YNames[40], "Jrelp");

   strcpy(YUnits[0], "millimolar");
   strcpy(YUnits[1], "dimensionless");
   strcpy(YUnits[2], "dimensionless");
   strcpy(YUnits[3], "dimensionless");
   strcpy(YUnits[4], "dimensionless");
   strcpy(YUnits[5], "dimensionless");
   strcpy(YUnits[6], "dimensionless");
   strcpy(YUnits[7], "dimensionless");
   strcpy(YUnits[8], "dimensionless");
   strcpy(YUnits[9], "dimensionless");
   strcpy(YUnits[10], "dimensionless");
   strcpy(YUnits[11], "dimensionless");
   strcpy(YUnits[12], "dimensionless");
   strcpy(YUnits[13], "dimensionless");
   strcpy(YUnits[14], "dimensionless");
   strcpy(YUnits[15], "dimensionless");
   strcpy(YUnits[16], "dimensionless");
   strcpy(YUnits[17], "dimensionless");
   strcpy(YUnits[18], "dimensionless");
   strcpy(YUnits[19], "dimensionless");
   strcpy(YUnits[20], "dimensionless");
   strcpy(YUnits[21], "dimensionless");
   strcpy(YUnits[22], "dimensionless");
   strcpy(YUnits[23], "dimensionless");
   strcpy(YUnits[24], "dimensionless");
   strcpy(YUnits[25], "dimensionless");
   strcpy(YUnits[26], "dimensionless");
   strcpy(YUnits[27], "dimensionless");
   strcpy(YUnits[28], "dimensionless");
   strcpy(YUnits[29], "dimensionless");
   strcpy(YUnits[30], "millimolar");
   strcpy(YUnits[31], "millimolar");
   strcpy(YUnits[32], "millimolar");
   strcpy(YUnits[33], "millimolar");
   strcpy(YUnits[34], "millimolar");
   strcpy(YUnits[35], "millimolar");
   strcpy(YUnits[36], "millimolar");
   strcpy(YUnits[37], "millimolar");
   strcpy(YUnits[38], "millivolt");
   strcpy(YUnits[39], "dimensionless");
   strcpy(YUnits[40], "dimensionless");

   strcpy(YComponents[0], "CaMK");
   strcpy(YComponents[1], "ICaL");
   strcpy(YComponents[2], "ICaL");
   strcpy(YComponents[3], "ICaL");
   strcpy(YComponents[4], "ICaL");
   strcpy(YComponents[5], "ICaL");
   strcpy(YComponents[6], "ICaL");
   strcpy(YComponents[7], "ICaL");
   strcpy(YComponents[8], "ICaL");
   strcpy(YComponents[9], "ICaL");
   strcpy(YComponents[10], "IK1");
   strcpy(YComponents[11], "IKr");
   strcpy(YComponents[12], "IKr");
   strcpy(YComponents[13], "IKs");
   strcpy(YComponents[14], "IKs");
   strcpy(YComponents[15], "INaL");
   strcpy(YComponents[16], "INaL");
   strcpy(YComponents[17], "INaL");
   strcpy(YComponents[18], "INa");
   strcpy(YComponents[19], "INa");
   strcpy(YComponents[20], "INa");
   strcpy(YComponents[21], "INa");
   strcpy(YComponents[22], "INa");
   strcpy(YComponents[23], "INa");
   strcpy(YComponents[24], "Ito");
   strcpy(YComponents[25], "Ito");
   strcpy(YComponents[26], "Ito");
   strcpy(YComponents[27], "Ito");
   strcpy(YComponents[28], "Ito");
   strcpy(YComponents[29], "Ito");
   strcpy(YComponents[30], "intracellular_ions");
   strcpy(YComponents[31], "intracellular_ions");
   strcpy(YComponents[32], "intracellular_ions");
   strcpy(YComponents[33], "intracellular_ions");
   strcpy(YComponents[34], "intracellular_ions");
   strcpy(YComponents[35], "intracellular_ions");
   strcpy(YComponents[36], "intracellular_ions");
   strcpy(YComponents[37], "intracellular_ions");
   strcpy(YComponents[38], "membrane");
   strcpy(YComponents[39], "ryr");
   strcpy(YComponents[40], "ryr");

   //---------------------------------------------------------------------------
   // Constants
   //---------------------------------------------------------------------------

   CaMKo = 0.05;   // dimensionless (in CaMK)
   KmCaM = 0.0015;   // millimolar (in CaMK)
   KmCaMK = 0.15;   // millimolar (in CaMK)
   aCaMK = 0.05;   // per_millimolar_per_millisecond (in CaMK)
   bCaMK = 0.00068;   // per_millisecond (in CaMK)
   Kmn = 0.002;   // millimolar (in ICaL)
   PCa_b = 0.0001;   // dimensionless (in ICaL)
   k2n = 1000.0;   // per_millisecond (in ICaL)
   PCab = 2.5e-8;   // milliS_per_microF (in ICab)
   GK1_b = 0.1908;   // milliS_per_microF (in IK1)
   GKb_b = 0.003;   // milliS_per_microF (in IKb)
   GKr_b = 0.046;   // milliS_per_microF (in IKr)
   GKs_b = 0.0034;   // milliS_per_microF (in IKs)
   Gncx_b = 0.0008;   // milliS_per_microF (in INaCa_i)
   KmCaAct = 150.0e-6;   // millimolar (in INaCa_i)
   kasymm = 12.5;   // dimensionless (in INaCa_i)
   kcaoff = 5.0e3;   // per_millisecond (in INaCa_i)
   kcaon = 1.5e6;   // per_millisecond (in INaCa_i)
   kna1 = 15.0;   // per_millisecond (in INaCa_i)
   kna2 = 5.0;   // per_millisecond (in INaCa_i)
   kna3 = 88.12;   // per_millisecond (in INaCa_i)
   qca = 0.167;   // dimensionless (in INaCa_i)
   qna = 0.5224;   // dimensionless (in INaCa_i)
   wca = 6.0e4;   // dimensionless (in INaCa_i)
   wna = 6.0e4;   // dimensionless (in INaCa_i)
   wnaca = 5.0e3;   // dimensionless (in INaCa_i)
   H = 1.0e-7;   // millimolar (in INaK)
   Khp = 1.698e-7;   // millimolar (in INaK)
   Kki = 0.5;   // per_millisecond (in INaK)
   Kko = 0.3582;   // per_millisecond (in INaK)
   Kmgatp = 1.698e-7;   // millimolar (in INaK)
   Knai0 = 9.073;   // millimolar (in INaK)
   Knao0 = 27.78;   // millimolar (in INaK)
   Knap = 224.0;   // millimolar (in INaK)
   Kxkur = 292.0;   // millimolar (in INaK)
   MgADP = 0.05;   // millimolar (in INaK)
   MgATP = 9.8;   // millimolar (in INaK)
   Pnak_b = 30.0;   // milliS_per_microF (in INaK)
   delta = -0.155;   // millivolt (in INaK)
   eP = 4.2;   // dimensionless (in INaK)
   k1m = 182.4;   // per_millisecond (in INaK)
   k1p = 949.5;   // per_millisecond (in INaK)
   k2m = 39.4;   // per_millisecond (in INaK)
   k2p = 687.2;   // per_millisecond (in INaK)
   k3m = 79300.0;   // per_millisecond (in INaK)
   k3p = 1899.0;   // per_millisecond (in INaK)
   k4m = 40.0;   // per_millisecond (in INaK)
   k4p = 639.0;   // per_millisecond (in INaK)
   GNaL_b = 0.0075;   // milliS_per_microF (in INaL)
   thL = 200.0;   // millisecond (in INaL)
   PNab = 3.75e-10;   // milliS_per_microF (in INab)
   Ahf = 0.99;   // dimensionless (in INa)
   GNa = 75.0;   // milliS_per_microF (in INa)
   hssV1 = 82.9;   // millivolt (in INa)
   hssV2 = 6.086;   // millivolt (in INa)
   mssV1 = 39.57;   // millivolt (in INa)
   mssV2 = 9.871;   // millivolt (in INa)
   mtD1 = 6.765;   // dimensionless (in INa)
   mtD2 = 8.552;   // dimensionless (in INa)
   mtV1 = 11.64;   // millivolt (in INa)
   mtV2 = 34.77;   // millivolt (in INa)
   mtV3 = 77.42;   // millivolt (in INa)
   mtV4 = 5.955;   // millivolt (in INa)
   GpCa = 0.0005;   // milliS_per_microF (in IpCa)
   KmCap = 0.0005;   // millimolar (in IpCa)
   Gto_b = 0.02;   // milliS_per_microF (in Ito)
   L = 0.01;   // centimeter (in cell_geometry)
   rad = 0.0011;   // centimeter (in cell_geometry)
   celltype = 0.0;   // dimensionless (in environment)
   cao = 1.8;   // millimolar (in extracellular)
   ko = 5.4;   // millimolar (in extracellular)
   nao = 140.0;   // millimolar (in extracellular)
   BSLmax = 1.124;   // millimolar (in intracellular_ions)
   BSRmax = 0.047;   // millimolar (in intracellular_ions)
   KmBSL = 0.0087;   // millimolar (in intracellular_ions)
   KmBSR = 0.00087;   // millimolar (in intracellular_ions)
   cm = 1.0;   // microF_per_centimeter_squared (in intracellular_ions)
   cmdnmax_b = 0.05;   // millimolar (in intracellular_ions)
   csqnmax = 10.0;   // millimolar (in intracellular_ions)
   kmcmdn = 0.00238;   // millimolar (in intracellular_ions)
   kmcsqn = 0.8;   // millimolar (in intracellular_ions)
   kmtrpn = 0.0005;   // millimolar (in intracellular_ions)
   trpnmax = 0.07;   // millimolar (in intracellular_ions)
   amp = -80.0;   // microA_per_microF (in membrane)
   duration = 0.5;   // millisecond (in membrane)
   F = 96485.0;   // coulomb_per_mole (in physical_constants)
   R = 8314.0;   // joule_per_kilomole_kelvin (in physical_constants)
   T = 310.0;   // kelvin (in physical_constants)
   zca = 2.0;   // dimensionless (in physical_constants)
   zk = 1.0;   // dimensionless (in physical_constants)
   zna = 1.0;   // dimensionless (in physical_constants)
   PKNa = 0.01833;   // dimensionless (in reversal_potentials)
   bt = 4.75;   // millisecond (in ryr)

   //---------------------------------------------------------------------------
   // Computed variables
   //---------------------------------------------------------------------------

   Aff = 0.6;
   Afs = 1.0-Aff;
   tjca = 75.0;

   if (celltype == 1.0)
      PCa = PCa_b*1.2;
   else if (celltype == 2.0)
      PCa = PCa_b*2.5;
   else
      PCa = PCa_b;

   PCap = 1.1*PCa;
   PCaNa = 0.00125*PCa;
   PCaK = 3.574e-4*PCa;
   PCaNap = 0.00125*PCap;
   PCaKp = 3.574e-4*PCap;

   if (celltype == 1.0)
      GK1 = GK1_b*1.2;
   else if (celltype == 2.0)
      GK1 = GK1_b*1.3;
   else
      GK1 = GK1_b;

   if (celltype == 1.0)
      GKb = GKb_b*0.6;
   else
      GKb = GKb_b;

   if (celltype == 1.0)
      GKr = GKr_b*1.3;
   else if (celltype == 2.0)
      GKr = GKr_b*0.8;
   else
      GKr = GKr_b;

   if (celltype == 1.0)
      GKs = GKs_b*1.4;
   else
      GKs = GKs_b;

   Ahs = 1.0-Ahf;
   h10_i = kasymm+1.0+nao/kna1*(1.0+nao/kna2);
   h11_i = nao*nao/(h10_i*kna1*kna2);
   h12_i = 1.0/h10_i;
   k1_i = h12_i*cao*kcaon;
   k2_i = kcaoff;
   k5_i = kcaoff;

   if (celltype == 1.0)
      Gncx = Gncx_b*1.1;
   else if (celltype == 2.0)
      Gncx = Gncx_b*1.4;
   else
      Gncx = Gncx_b;

   h10_ss = kasymm+1.0+nao/kna1*(1.0+nao/kna2);
   h11_ss = nao*nao/(h10_ss*kna1*kna2);
   h12_ss = 1.0/h10_ss;
   k1_ss = h12_ss*cao*kcaon;
   k2_ss = kcaoff;
   k5_ss = kcaoff;
   b1 = k1m*MgADP;
   a2 = k2p;
   a4 = k4p*MgATP/Kmgatp/(1.0+MgATP/Kmgatp);

   if (celltype == 1.0)
      Pnak = Pnak_b*0.9;
   else if (celltype == 2.0)
      Pnak = Pnak_b*0.7;
   else
      Pnak = Pnak_b;

   thLp = 3.0*thL;

   if (celltype == 1.0)
      GNaL = GNaL_b*0.6;
   else
      GNaL = GNaL_b;

   if (celltype == 1.0)
      Gto = Gto_b*4.0;
   else if (celltype == 2.0)
      Gto = Gto_b*4.0;
   else
      Gto = Gto_b;

   if (celltype == 1.0)
      upScale = 1.3;
   else
      upScale = 1.0;

   vcell = 1000.0*3.14*rad*rad*L;
   Ageo = 2.0*3.14*rad*rad+2.0*3.14*rad*L;
   Acap = 2.0*Ageo;
   vmyo = 0.68*vcell;
   vnsr = 0.0552*vcell;
   vjsr = 0.0048*vcell;
   vss = 0.02*vcell;

   if (celltype == 1.0)
      cmdnmax = cmdnmax_b*1.3;
   else
      cmdnmax = cmdnmax_b;

   a_rel = 0.5*bt;
   btp = 1.25*bt;
   a_relp = 0.5*btp;
}

//------------------------------------------------------------------------------
// Computation
//------------------------------------------------------------------------------

void compute(double time)
{
   // time: time (millisecond)

   CaMKb = CaMKo*(1.0-Y[0])/(1.0+KmCaM/Y[33]);
   CaMKa = CaMKb+Y[0];
   dY[0] = aCaMK*CaMKb*(CaMKb+Y[0])-bCaMK*Y[0];
   dss = 1.0/(1.0+exp(-(Y[38]+3.94)/4.23));
   td = 0.6+1.0/(exp(-0.05*(Y[38]+6.0))+exp(0.09*(Y[38]+14.0)));
   dY[1] = (dss-Y[1])/td;
   fss = 1.0/(1.0+exp((Y[38]+19.58)/3.696));
   tff = 7.0+1.0/(0.0045*exp(-(Y[38]+20.0)/10.0)+0.0045*exp((Y[38]+20.0)/10.0));
   tfs = 1000.0+1.0/(0.000035*exp(-(Y[38]+5.0)/4.0)+0.000035*exp((Y[38]+5.0)/6.0));
   dY[5] = (fss-Y[5])/tff;
   dY[7] = (fss-Y[7])/tfs;
   f = Aff*Y[5]+Afs*Y[7];
   fcass = fss;
   tfcaf = 7.0+1.0/(0.04*exp(-(Y[38]-4.0)/7.0)+0.04*exp((Y[38]-4.0)/7.0));
   tfcas = 100.0+1.0/(0.00012*exp(-Y[38]/3.0)+0.00012*exp(Y[38]/7.0));
   Afcaf = 0.3+0.6/(1.0+exp((Y[38]-10.0)/10.0));
   Afcas = 1.0-Afcaf;
   dY[2] = (fcass-Y[2])/tfcaf;
   dY[4] = (fcass-Y[4])/tfcas;
   fca = Afcaf*Y[2]+Afcas*Y[4];
   dY[8] = (fcass-Y[8])/tjca;
   tffp = 2.5*tff;
   dY[6] = (fss-Y[6])/tffp;
   fp = Aff*Y[6]+Afs*Y[7];
   tfcafp = 2.5*tfcaf;
   dY[3] = (fcass-Y[3])/tfcafp;
   fcap = Afcaf*Y[3]+Afcas*Y[4];
   km2n = Y[8]*1.0;
   anca = 1.0/(k2n/km2n+pow(1.0+Kmn/Y[33], 4.0));
   dY[9] = anca*k2n-Y[9]*km2n;
   vffrt = Y[38]*F*F/(R*T);
   vfrt = Y[38]*F/(R*T);
   PhiCaL = 4.0*vffrt*(Y[33]*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
   PhiCaNa = 1.0*vffrt*(0.75*Y[37]*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
   PhiCaK = 1.0*vffrt*(0.75*Y[35]*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
   fICaLp = 1.0/(1.0+KmCaMK/CaMKa);
   ICaL = (1.0-fICaLp)*PCa*PhiCaL*Y[1]*(f*(1.0-Y[9])+Y[8]*fca*Y[9])+fICaLp*PCap*PhiCaL*Y[1]*(fp*(1.0-Y[9])+Y[8]*fcap*Y[9]);
   ICaNa = (1.0-fICaLp)*PCaNa*PhiCaNa*Y[1]*(f*(1.0-Y[9])+Y[8]*fca*Y[9])+fICaLp*PCaNap*PhiCaNa*Y[1]*(fp*(1.0-Y[9])+Y[8]*fcap*Y[9]);
   ICaK = (1.0-fICaLp)*PCaK*PhiCaK*Y[1]*(f*(1.0-Y[9])+Y[8]*fca*Y[9])+fICaLp*PCaKp*PhiCaK*Y[1]*(fp*(1.0-Y[9])+Y[8]*fcap*Y[9]);
   ICab = PCab*4.0*vffrt*(Y[30]*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
   xk1ss = 1.0/(1.0+exp(-(Y[38]+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
   txk1 = 122.2/(exp(-(Y[38]+127.2)/20.36)+exp((Y[38]+236.8)/69.33));
   dY[10] = (xk1ss-Y[10])/txk1;
   rk1 = 1.0/(1.0+exp((Y[38]+105.8-2.6*ko)/9.493));
   EK = R*T/F*log(ko/Y[34]);
   IK1 = GK1*sqrt(ko)*rk1*Y[10]*(Y[38]-EK);
   xkb = 1.0/(1.0+exp(-(Y[38]-14.48)/18.34));
   IKb = GKb*xkb*(Y[38]-EK);
   xrss = 1.0/(1.0+exp(-(Y[38]+8.337)/6.789));
   txrf = 12.98+1.0/(0.3652*exp((Y[38]-31.66)/3.869)+4.123e-5*exp(-(Y[38]-47.78)/20.38));
   txrs = 1.865+1.0/(0.06629*exp((Y[38]-34.7)/7.355)+1.128e-5*exp(-(Y[38]-29.74)/25.94));
   Axrf = 1.0/(1.0+exp((Y[38]+54.81)/38.21));
   Axrs = 1.0-Axrf;
   dY[11] = (xrss-Y[11])/txrf;
   dY[12] = (xrss-Y[12])/txrs;
   xr = Axrf*Y[11]+Axrs*Y[12];
   rkr = 1.0/(1.0+exp((Y[38]+55.0)/75.0))*1.0/(1.0+exp((Y[38]-10.0)/30.0));
   IKr = GKr*sqrt(ko/5.4)*xr*rkr*(Y[38]-EK);
   xs1ss = 1.0/(1.0+exp(-(Y[38]+11.6)/8.932));
   txs1 = 817.3+1.0/(2.326e-4*exp((Y[38]+48.28)/17.8)+0.001292*exp(-(Y[38]+210.0)/230.0));
   dY[13] = (xs1ss-Y[13])/txs1;
   xs2ss = xs1ss;
   txs2 = 1.0/(0.01*exp((Y[38]-50.0)/20.0)+0.0193*exp(-(Y[38]+66.54)/31.0));
   dY[14] = (xs2ss-Y[14])/txs2;
   KsCa = 1.0+0.6/(1.0+pow(3.8e-5/Y[30], 1.4));
   EKs = R*T/F*log((ko+PKNa*nao)/(Y[34]+PKNa*Y[36]));
   IKs = GKs*KsCa*Y[13]*Y[14]*(Y[38]-EKs);
   mss = 1.0/(1.0+exp(-(Y[38]+mssV1)/mssV2));
   tm = 1.0/(mtD1*exp((Y[38]+mtV1)/mtV2)+mtD2*exp(-(Y[38]+mtV3)/mtV4));
   dY[23] = (mss-Y[23])/tm;
   hss = 1.0/(1.0+exp((Y[38]+hssV1)/hssV2));
   thf = 1.0/(1.432e-5*exp(-(Y[38]+1.196)/6.285)+6.149*exp((Y[38]+0.5096)/20.27));
   ths = 1.0/(0.009794*exp(-(Y[38]+17.95)/28.05)+0.3343*exp((Y[38]+5.73)/56.66));
   dY[18] = (hss-Y[18])/thf;
   dY[19] = (hss-Y[19])/ths;
   h = Ahf*Y[18]+Ahs*Y[19];
   jss = hss;
   tj = 2.038+1.0/(0.02136*exp(-(Y[38]+100.6)/8.281)+0.3052*exp((Y[38]+0.9941)/38.45));
   dY[21] = (jss-Y[21])/tj;
   hssp = 1.0/(1.0+exp((Y[38]+89.1)/6.086));
   thsp = 3.0*ths;
   dY[20] = (hssp-Y[20])/thsp;
   hp = Ahf*Y[18]+Ahs*Y[20];
   tjp = 1.46*tj;
   dY[22] = (jss-Y[22])/tjp;
   fINap = 1.0/(1.0+KmCaMK/CaMKa);
   ENa = R*T/F*log(nao/Y[36]);
   INa = GNa*(Y[38]-ENa)*pow(Y[23], 3.0)*((1.0-fINap)*h*Y[21]+fINap*hp*Y[22]);
   hca = exp(qca*Y[38]*F/(R*T));
   hna = exp(qna*Y[38]*F/(R*T));
   h1_i = 1.0+Y[36]/kna3*(1.0+hna);
   h2_i = Y[36]*hna/(kna3*h1_i);
   h3_i = 1.0/h1_i;
   h4_i = 1.0+Y[36]/kna1*(1.0+Y[36]/kna2);
   h5_i = Y[36]*Y[36]/(h4_i*kna1*kna2);
   h6_i = 1.0/h4_i;
   h7_i = 1.0+nao/kna3*(1.0+1.0/hna);
   h8_i = nao/(kna3*hna*h7_i);
   h9_i = 1.0/h7_i;
   k3p_i = h9_i*wca;
   k3pp_i = h8_i*wnaca;
   k3_i = k3p_i+k3pp_i;
   k4p_i = h3_i*wca/hca;
   k4pp_i = h2_i*wnaca;
   k4_i = k4p_i+k4pp_i;
   k6_i = h6_i*Y[30]*kcaon;
   k7_i = h5_i*h2_i*wna;
   k8_i = h8_i*h11_i*wna;
   x1_i = k2_i*k4_i*(k7_i+k6_i)+k5_i*k7_i*(k2_i+k3_i);
   x2_i = k1_i*k7_i*(k4_i+k5_i)+k4_i*k6_i*(k1_i+k8_i);
   x3_i = k1_i*k3_i*(k7_i+k6_i)+k8_i*k6_i*(k2_i+k3_i);
   x4_i = k2_i*k8_i*(k4_i+k5_i)+k3_i*k5_i*(k1_i+k8_i);
   E1_i = x1_i/(x1_i+x2_i+x3_i+x4_i);
   E2_i = x2_i/(x1_i+x2_i+x3_i+x4_i);
   E3_i = x3_i/(x1_i+x2_i+x3_i+x4_i);
   E4_i = x4_i/(x1_i+x2_i+x3_i+x4_i);
   allo_i = 1.0/(1.0+pow(KmCaAct/Y[30], 2.0));
   JncxNa_i = 3.0*(E4_i*k7_i-E1_i*k8_i)+E3_i*k4pp_i-E2_i*k3pp_i;
   JncxCa_i = E2_i*k2_i-E1_i*k1_i;
   INaCa_i = 0.8*Gncx*allo_i*(zna*JncxNa_i+zca*JncxCa_i);
   h1_ss = 1.0+Y[37]/kna3*(1.0+hna);
   h2_ss = Y[37]*hna/(kna3*h1_ss);
   h3_ss = 1.0/h1_ss;
   h4_ss = 1.0+Y[37]/kna1*(1.0+Y[37]/kna2);
   h5_ss = Y[37]*Y[37]/(h4_ss*kna1*kna2);
   h6_ss = 1.0/h4_ss;
   h7_ss = 1.0+nao/kna3*(1.0+1.0/hna);
   h8_ss = nao/(kna3*hna*h7_ss);
   h9_ss = 1.0/h7_ss;
   k3p_ss = h9_ss*wca;
   k3pp_ss = h8_ss*wnaca;
   k3_ss = k3p_ss+k3pp_ss;
   k4p_ss = h3_ss*wca/hca;
   k4pp_ss = h2_ss*wnaca;
   k4_ss = k4p_ss+k4pp_ss;
   k6_ss = h6_ss*Y[33]*kcaon;
   k7_ss = h5_ss*h2_ss*wna;
   k8_ss = h8_ss*h11_ss*wna;
   x1_ss = k2_ss*k4_ss*(k7_ss+k6_ss)+k5_ss*k7_ss*(k2_ss+k3_ss);
   x2_ss = k1_ss*k7_ss*(k4_ss+k5_ss)+k4_ss*k6_ss*(k1_ss+k8_ss);
   x3_ss = k1_ss*k3_ss*(k7_ss+k6_ss)+k8_ss*k6_ss*(k2_ss+k3_ss);
   x4_ss = k2_ss*k8_ss*(k4_ss+k5_ss)+k3_ss*k5_ss*(k1_ss+k8_ss);
   E1_ss = x1_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   E2_ss = x2_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   E3_ss = x3_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   E4_ss = x4_ss/(x1_ss+x2_ss+x3_ss+x4_ss);
   allo_ss = 1.0/(1.0+pow(KmCaAct/Y[33], 2.0));
   JncxNa_ss = 3.0*(E4_ss*k7_ss-E1_ss*k8_ss)+E3_ss*k4pp_ss-E2_ss*k3pp_ss;
   JncxCa_ss = E2_ss*k2_ss-E1_ss*k1_ss;
   INaCa_ss = 0.2*Gncx*allo_ss*(zna*JncxNa_ss+zca*JncxCa_ss);
   Knai = Knai0*exp(delta*Y[38]*F/(3.0*R*T));
   Knao = Knao0*exp((1.0-delta)*Y[38]*F/(3.0*R*T));
   P = eP/(1.0+H/Khp+Y[36]/Knap+Y[34]/Kxkur);
   a1 = k1p*pow(Y[36]/Knai, 3.0)/(pow(1.0+Y[36]/Knai, 3.0)+pow(1.0+Y[34]/Kki, 2.0)-1.0);
   b2 = k2m*pow(nao/Knao, 3.0)/(pow(1.0+nao/Knao, 3.0)+pow(1.0+ko/Kko, 2.0)-1.0);
   a3 = k3p*pow(ko/Kko, 2.0)/(pow(1.0+nao/Knao, 3.0)+pow(1.0+ko/Kko, 2.0)-1.0);
   b3 = k3m*P*H/(1.0+MgATP/Kmgatp);
   b4 = k4m*pow(Y[34]/Kki, 2.0)/(pow(1.0+Y[36]/Knai, 3.0)+pow(1.0+Y[34]/Kki, 2.0)-1.0);
   x1 = a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
   x2 = b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
   x3 = a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
   x4 = b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
   E1 = x1/(x1+x2+x3+x4);
   E2 = x2/(x1+x2+x3+x4);
   E3 = x3/(x1+x2+x3+x4);
   E4 = x4/(x1+x2+x3+x4);
   JnakNa = 3.0*(E1*a3-E2*b3);
   JnakK = 2.0*(E4*b1-E3*a1);
   INaK = Pnak*(zna*JnakNa+zk*JnakK);
   mLss = 1.0/(1.0+exp(-(Y[38]+42.85)/5.264));
   tmL = tm;
   dY[17] = (mLss-Y[17])/tmL;
   hLss = 1.0/(1.0+exp((Y[38]+87.61)/7.488));
   dY[15] = (hLss-Y[15])/thL;
   hLssp = 1.0/(1.0+exp((Y[38]+93.81)/7.488));
   dY[16] = (hLssp-Y[16])/thLp;
   fINaLp = 1.0/(1.0+KmCaMK/CaMKa);
   INaL = GNaL*(Y[38]-ENa)*Y[17]*((1.0-fINaLp)*Y[15]+fINaLp*Y[16]);
   INab = PNab*vffrt*(Y[36]*exp(vfrt)-nao)/(exp(vfrt)-1.0);
   IpCa = GpCa*Y[30]/(KmCap+Y[30]);
   ass = 1.0/(1.0+exp(-(Y[38]-14.34)/14.82));
   ta = 1.0515/(1.0/(1.2089*(1.0+exp(-(Y[38]-18.4099)/29.3814)))+3.5/(1.0+exp((Y[38]+100.0)/29.3814)));
   dY[24] = (ass-Y[24])/ta;
   iss = 1.0/(1.0+exp((Y[38]+43.94)/5.711));

   if (celltype == 1.0)
      delta_epi = 1.0-0.95/(1.0+exp((Y[38]+70.0)/5.0));
   else
      delta_epi = 1.0;

   tiF_b = 4.562+1.0/(0.3933*exp(-(Y[38]+100.0)/100.0)+0.08004*exp((Y[38]+50.0)/16.59));
   tiS_b = 23.62+1.0/(0.001416*exp(-(Y[38]+96.52)/59.05)+1.78e-8*exp((Y[38]+114.1)/8.079));
   tiF = tiF_b*delta_epi;
   tiS = tiS_b*delta_epi;
   AiF = 1.0/(1.0+exp((Y[38]-213.6)/151.2));
   AiS = 1.0-AiF;
   dY[26] = (iss-Y[26])/tiF;
   dY[28] = (iss-Y[28])/tiS;
   i = AiF*Y[26]+AiS*Y[28];
   assp = 1.0/(1.0+exp(-(Y[38]-24.34)/14.82));
   dY[25] = (assp-Y[25])/ta;
   dti_develop = 1.354+1.0e-4/(exp((Y[38]-167.4)/15.89)+exp(-(Y[38]-12.23)/0.2154));
   dti_recover = 1.0-0.5/(1.0+exp((Y[38]+70.0)/20.0));
   tiFp = dti_develop*dti_recover*tiF;
   tiSp = dti_develop*dti_recover*tiS;
   dY[27] = (iss-Y[27])/tiFp;
   dY[29] = (iss-Y[29])/tiSp;
   ip = AiF*Y[27]+AiS*Y[29];
   fItop = 1.0/(1.0+KmCaMK/CaMKa);
   Ito = Gto*(Y[38]-EK)*((1.0-fItop)*Y[24]*i+fItop*Y[25]*ip);
   Jupnp = upScale*0.004375*Y[30]/(Y[30]+0.00092);
   Jupp = upScale*2.75*0.004375*Y[30]/(Y[30]+0.00092-0.00017);
   fJupp = 1.0/(1.0+KmCaMK/CaMKa);
   Jleak = 0.0039375*Y[32]/15.0;
   Jup = (1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
   JdiffNa = (Y[37]-Y[36])/2.0;
   JdiffK = (Y[35]-Y[34])/2.0;
   Jdiff = (Y[33]-Y[30])/0.2;
   dY[36] = -(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap*cm/(F*vmyo)+JdiffNa*vss/vmyo;
   dY[37] = -(ICaNa+3.0*INaCa_ss)*cm*Acap/(F*vss)-JdiffNa;

   if (time <= duration)
      Istim = amp;
   else
      Istim = 0.0;

   dY[34] = -(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*cm*Acap/(F*vmyo)+JdiffK*vss/vmyo;
   dY[35] = -ICaK*cm*Acap/(F*vss)-JdiffK;
   Bcai = 1.0/(1.0+cmdnmax*kmcmdn/pow(kmcmdn+Y[30], 2.0)+trpnmax*kmtrpn/pow(kmtrpn+Y[30], 2.0));
   dY[30] = Bcai*(-(IpCa+ICab-2.0*INaCa_i)*cm*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);
   Bcass = 1.0/(1.0+BSRmax*KmBSR/pow(KmBSR+Y[33], 2.0)+BSLmax*KmBSL/pow(KmBSL+Y[33], 2.0));
   fJrelp = 1.0/(1.0+KmCaMK/CaMKa);
   Jrel = (1.0-fJrelp)*Y[39]+fJrelp*Y[40];
   dY[33] = Bcass*(-(ICaL-2.0*INaCa_ss)*cm*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);
   Jtr = (Y[32]-Y[31])/100.0;
   dY[32] = Jup-Jtr*vjsr/vnsr;
   Bcajsr = 1.0/(1.0+csqnmax*kmcsqn/pow(kmcsqn+Y[31], 2.0));
   dY[31] = Bcajsr*(Jtr-Jrel);
   dY[38] = -(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim);
   Jrel_inf_temp = a_rel*-ICaL/(1.0+1.0*pow(1.5/Y[31], 8.0));

   if (celltype == 2.0)
      Jrel_inf = Jrel_inf_temp*1.7;
   else
      Jrel_inf = Jrel_inf_temp;

   tau_rel_temp = bt/(1.0+0.0123/Y[31]);

   if (tau_rel_temp < 0.001)
      tau_rel = 0.001;
   else
      tau_rel = tau_rel_temp;

   dY[39] = (Jrel_inf-Y[39])/tau_rel;
   Jrel_temp = a_relp*-ICaL/(1.0+pow(1.5/Y[31], 8.0));

   if (celltype == 2.0)
      Jrel_infp = Jrel_temp*1.7;
   else
      Jrel_infp = Jrel_temp;

   tau_relp_temp = btp/(1.0+0.0123/Y[31]);

   if (tau_relp_temp < 0.001)
      tau_relp = 0.001;
   else
      tau_relp = tau_relp_temp;

   dY[40] = (Jrel_infp-Y[40])/tau_relp;
}

//==============================================================================
// End of file
//==============================================================================
