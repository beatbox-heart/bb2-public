//==============================================================================
// CellML file:   C:\Users\sanjay\Desktop\Ohara_Rudy_2011.cellml
// CellML model:  MyModel1
// Date and time: 13/12/2011 at 12:58:02
//------------------------------------------------------------------------------
// Conversion from CellML 1.0 to C (header) was done using COR (0.9.31.1409)
//    Copyright 2002-2011 Dr Alan Garny
//    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
//------------------------------------------------------------------------------
// http://www.cellml.org/
//==============================================================================

#ifndef __OHARA_RUDY_2011_H__
#define __OHARA_RUDY_2011_H__

//------------------------------------------------------------------------------
// State variables
//------------------------------------------------------------------------------

#define _NB_OF_STATE_VARIABLES_ 41

extern double Y[_NB_OF_STATE_VARIABLES_];
extern double dY[_NB_OF_STATE_VARIABLES_];
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

extern char YNames[_NB_OF_STATE_VARIABLES_][7];
extern char YUnits[_NB_OF_STATE_VARIABLES_][14];
extern char YComponents[_NB_OF_STATE_VARIABLES_][19];

//------------------------------------------------------------------------------
// Constants
//------------------------------------------------------------------------------

extern double CaMKo;   // dimensionless (in CaMK)
extern double KmCaM;   // millimolar (in CaMK)
extern double KmCaMK;   // millimolar (in CaMK)
extern double aCaMK;   // per_millimolar_per_millisecond (in CaMK)
extern double bCaMK;   // per_millisecond (in CaMK)
extern double Kmn;   // millimolar (in ICaL)
extern double PCa_b;   // dimensionless (in ICaL)
extern double k2n;   // per_millisecond (in ICaL)
extern double PCab;   // milliS_per_microF (in ICab)
extern double GK1_b;   // milliS_per_microF (in IK1)
extern double GKb_b;   // milliS_per_microF (in IKb)
extern double GKr_b;   // milliS_per_microF (in IKr)
extern double GKs_b;   // milliS_per_microF (in IKs)
extern double Gncx_b;   // milliS_per_microF (in INaCa_i)
extern double KmCaAct;   // millimolar (in INaCa_i)
extern double kasymm;   // dimensionless (in INaCa_i)
extern double kcaoff;   // per_millisecond (in INaCa_i)
extern double kcaon;   // per_millisecond (in INaCa_i)
extern double kna1;   // per_millisecond (in INaCa_i)
extern double kna2;   // per_millisecond (in INaCa_i)
extern double kna3;   // per_millisecond (in INaCa_i)
extern double qca;   // dimensionless (in INaCa_i)
extern double qna;   // dimensionless (in INaCa_i)
extern double wca;   // dimensionless (in INaCa_i)
extern double wna;   // dimensionless (in INaCa_i)
extern double wnaca;   // dimensionless (in INaCa_i)
extern double H;   // millimolar (in INaK)
extern double Khp;   // millimolar (in INaK)
extern double Kki;   // per_millisecond (in INaK)
extern double Kko;   // per_millisecond (in INaK)
extern double Kmgatp;   // millimolar (in INaK)
extern double Knai0;   // millimolar (in INaK)
extern double Knao0;   // millimolar (in INaK)
extern double Knap;   // millimolar (in INaK)
extern double Kxkur;   // millimolar (in INaK)
extern double MgADP;   // millimolar (in INaK)
extern double MgATP;   // millimolar (in INaK)
extern double Pnak_b;   // milliS_per_microF (in INaK)
extern double delta;   // millivolt (in INaK)
extern double eP;   // dimensionless (in INaK)
extern double k1m;   // per_millisecond (in INaK)
extern double k1p;   // per_millisecond (in INaK)
extern double k2m;   // per_millisecond (in INaK)
extern double k2p;   // per_millisecond (in INaK)
extern double k3m;   // per_millisecond (in INaK)
extern double k3p;   // per_millisecond (in INaK)
extern double k4m;   // per_millisecond (in INaK)
extern double k4p;   // per_millisecond (in INaK)
extern double GNaL_b;   // milliS_per_microF (in INaL)
extern double thL;   // millisecond (in INaL)
extern double PNab;   // milliS_per_microF (in INab)
extern double Ahf;   // dimensionless (in INa)
extern double GNa;   // milliS_per_microF (in INa)
extern double hssV1;   // millivolt (in INa)
extern double hssV2;   // millivolt (in INa)
extern double mssV1;   // millivolt (in INa)
extern double mssV2;   // millivolt (in INa)
extern double mtD1;   // dimensionless (in INa)
extern double mtD2;   // dimensionless (in INa)
extern double mtV1;   // millivolt (in INa)
extern double mtV2;   // millivolt (in INa)
extern double mtV3;   // millivolt (in INa)
extern double mtV4;   // millivolt (in INa)
extern double GpCa;   // milliS_per_microF (in IpCa)
extern double KmCap;   // millimolar (in IpCa)
extern double Gto_b;   // milliS_per_microF (in Ito)
extern double L;   // centimeter (in cell_geometry)
extern double rad;   // centimeter (in cell_geometry)
extern double celltype;   // dimensionless (in environment)
extern double cao;   // millimolar (in extracellular)
extern double ko;   // millimolar (in extracellular)
extern double nao;   // millimolar (in extracellular)
extern double BSLmax;   // millimolar (in intracellular_ions)
extern double BSRmax;   // millimolar (in intracellular_ions)
extern double KmBSL;   // millimolar (in intracellular_ions)
extern double KmBSR;   // millimolar (in intracellular_ions)
extern double cm;   // microF_per_centimeter_squared (in intracellular_ions)
extern double cmdnmax_b;   // millimolar (in intracellular_ions)
extern double csqnmax;   // millimolar (in intracellular_ions)
extern double kmcmdn;   // millimolar (in intracellular_ions)
extern double kmcsqn;   // millimolar (in intracellular_ions)
extern double kmtrpn;   // millimolar (in intracellular_ions)
extern double trpnmax;   // millimolar (in intracellular_ions)
extern double amp;   // microA_per_microF (in membrane)
extern double duration;   // millisecond (in membrane)
extern double F;   // coulomb_per_mole (in physical_constants)
extern double R;   // joule_per_kilomole_kelvin (in physical_constants)
extern double T;   // kelvin (in physical_constants)
extern double zca;   // dimensionless (in physical_constants)
extern double zk;   // dimensionless (in physical_constants)
extern double zna;   // dimensionless (in physical_constants)
extern double PKNa;   // dimensionless (in reversal_potentials)
extern double bt;   // millisecond (in ryr)

//------------------------------------------------------------------------------
// Computed variables
//------------------------------------------------------------------------------

extern double CaMKa;   // millimolar (in CaMK)
extern double CaMKb;   // millimolar (in CaMK)
extern double Afcaf;   // dimensionless (in ICaL)
extern double Afcas;   // dimensionless (in ICaL)
extern double Aff;   // dimensionless (in ICaL)
extern double Afs;   // dimensionless (in ICaL)
extern double ICaK;   // microA_per_microF (in ICaL)
extern double ICaL;   // microA_per_microF (in ICaL)
extern double ICaNa;   // microA_per_microF (in ICaL)
extern double PCa;   // dimensionless (in ICaL)
extern double PCaK;   // dimensionless (in ICaL)
extern double PCaKp;   // dimensionless (in ICaL)
extern double PCaNa;   // dimensionless (in ICaL)
extern double PCaNap;   // dimensionless (in ICaL)
extern double PCap;   // dimensionless (in ICaL)
extern double PhiCaK;   // dimensionless (in ICaL)
extern double PhiCaL;   // dimensionless (in ICaL)
extern double PhiCaNa;   // dimensionless (in ICaL)
extern double anca;   // dimensionless (in ICaL)
extern double dss;   // dimensionless (in ICaL)
extern double f;   // dimensionless (in ICaL)
extern double fICaLp;   // dimensionless (in ICaL)
extern double fca;   // dimensionless (in ICaL)
extern double fcap;   // dimensionless (in ICaL)
extern double fcass;   // dimensionless (in ICaL)
extern double fp;   // dimensionless (in ICaL)
extern double fss;   // dimensionless (in ICaL)
extern double km2n;   // per_millisecond (in ICaL)
extern double td;   // millisecond (in ICaL)
extern double tfcaf;   // millisecond (in ICaL)
extern double tfcafp;   // millisecond (in ICaL)
extern double tfcas;   // millisecond (in ICaL)
extern double tff;   // millisecond (in ICaL)
extern double tffp;   // millisecond (in ICaL)
extern double tfs;   // millisecond (in ICaL)
extern double tjca;   // millisecond (in ICaL)
extern double ICab;   // microA_per_microF (in ICab)
extern double GK1;   // milliS_per_microF (in IK1)
extern double IK1;   // microA_per_microF (in IK1)
extern double rk1;   // millisecond (in IK1)
extern double txk1;   // millisecond (in IK1)
extern double xk1ss;   // dimensionless (in IK1)
extern double GKb;   // milliS_per_microF (in IKb)
extern double IKb;   // microA_per_microF (in IKb)
extern double xkb;   // dimensionless (in IKb)
extern double Axrf;   // dimensionless (in IKr)
extern double Axrs;   // dimensionless (in IKr)
extern double GKr;   // milliS_per_microF (in IKr)
extern double IKr;   // microA_per_microF (in IKr)
extern double rkr;   // dimensionless (in IKr)
extern double txrf;   // millisecond (in IKr)
extern double txrs;   // millisecond (in IKr)
extern double xr;   // dimensionless (in IKr)
extern double xrss;   // dimensionless (in IKr)
extern double GKs;   // milliS_per_microF (in IKs)
extern double IKs;   // microA_per_microF (in IKs)
extern double KsCa;   // dimensionless (in IKs)
extern double txs1;   // millisecond (in IKs)
extern double txs2;   // millisecond (in IKs)
extern double xs1ss;   // dimensionless (in IKs)
extern double xs2ss;   // dimensionless (in IKs)
extern double E1_i;   // dimensionless (in INaCa_i)
extern double E1_ss;   // dimensionless (in INaCa_i)
extern double E2_i;   // dimensionless (in INaCa_i)
extern double E2_ss;   // dimensionless (in INaCa_i)
extern double E3_i;   // dimensionless (in INaCa_i)
extern double E3_ss;   // dimensionless (in INaCa_i)
extern double E4_i;   // dimensionless (in INaCa_i)
extern double E4_ss;   // dimensionless (in INaCa_i)
extern double Gncx;   // milliS_per_microF (in INaCa_i)
extern double INaCa_i;   // microA_per_microF (in INaCa_i)
extern double INaCa_ss;   // microA_per_microF (in INaCa_i)
extern double JncxCa_i;   // millimolar_per_millisecond (in INaCa_i)
extern double JncxCa_ss;   // millimolar_per_millisecond (in INaCa_i)
extern double JncxNa_i;   // millimolar_per_millisecond (in INaCa_i)
extern double JncxNa_ss;   // millimolar_per_millisecond (in INaCa_i)
extern double allo_i;   // dimensionless (in INaCa_i)
extern double allo_ss;   // dimensionless (in INaCa_i)
extern double h10_i;   // dimensionless (in INaCa_i)
extern double h10_ss;   // dimensionless (in INaCa_i)
extern double h11_i;   // dimensionless (in INaCa_i)
extern double h11_ss;   // dimensionless (in INaCa_i)
extern double h12_i;   // dimensionless (in INaCa_i)
extern double h12_ss;   // dimensionless (in INaCa_i)
extern double h1_i;   // dimensionless (in INaCa_i)
extern double h1_ss;   // dimensionless (in INaCa_i)
extern double h2_i;   // dimensionless (in INaCa_i)
extern double h2_ss;   // dimensionless (in INaCa_i)
extern double h3_i;   // dimensionless (in INaCa_i)
extern double h3_ss;   // dimensionless (in INaCa_i)
extern double h4_i;   // dimensionless (in INaCa_i)
extern double h4_ss;   // dimensionless (in INaCa_i)
extern double h5_i;   // dimensionless (in INaCa_i)
extern double h5_ss;   // dimensionless (in INaCa_i)
extern double h6_i;   // dimensionless (in INaCa_i)
extern double h6_ss;   // dimensionless (in INaCa_i)
extern double h7_i;   // dimensionless (in INaCa_i)
extern double h7_ss;   // dimensionless (in INaCa_i)
extern double h8_i;   // dimensionless (in INaCa_i)
extern double h8_ss;   // dimensionless (in INaCa_i)
extern double h9_i;   // dimensionless (in INaCa_i)
extern double h9_ss;   // dimensionless (in INaCa_i)
extern double hca;   // dimensionless (in INaCa_i)
extern double hna;   // dimensionless (in INaCa_i)
extern double k1_i;   // dimensionless (in INaCa_i)
extern double k1_ss;   // dimensionless (in INaCa_i)
extern double k2_i;   // dimensionless (in INaCa_i)
extern double k2_ss;   // dimensionless (in INaCa_i)
extern double k3_i;   // dimensionless (in INaCa_i)
extern double k3_ss;   // dimensionless (in INaCa_i)
extern double k3p_i;   // dimensionless (in INaCa_i)
extern double k3p_ss;   // dimensionless (in INaCa_i)
extern double k3pp_i;   // dimensionless (in INaCa_i)
extern double k3pp_ss;   // dimensionless (in INaCa_i)
extern double k4_i;   // dimensionless (in INaCa_i)
extern double k4_ss;   // dimensionless (in INaCa_i)
extern double k4p_i;   // dimensionless (in INaCa_i)
extern double k4p_ss;   // dimensionless (in INaCa_i)
extern double k4pp_i;   // dimensionless (in INaCa_i)
extern double k4pp_ss;   // dimensionless (in INaCa_i)
extern double k5_i;   // dimensionless (in INaCa_i)
extern double k5_ss;   // dimensionless (in INaCa_i)
extern double k6_i;   // dimensionless (in INaCa_i)
extern double k6_ss;   // dimensionless (in INaCa_i)
extern double k7_i;   // dimensionless (in INaCa_i)
extern double k7_ss;   // dimensionless (in INaCa_i)
extern double k8_i;   // dimensionless (in INaCa_i)
extern double k8_ss;   // dimensionless (in INaCa_i)
extern double x1_i;   // dimensionless (in INaCa_i)
extern double x1_ss;   // dimensionless (in INaCa_i)
extern double x2_i;   // dimensionless (in INaCa_i)
extern double x2_ss;   // dimensionless (in INaCa_i)
extern double x3_i;   // dimensionless (in INaCa_i)
extern double x3_ss;   // dimensionless (in INaCa_i)
extern double x4_i;   // dimensionless (in INaCa_i)
extern double x4_ss;   // dimensionless (in INaCa_i)
extern double E1;   // dimensionless (in INaK)
extern double E2;   // dimensionless (in INaK)
extern double E3;   // dimensionless (in INaK)
extern double E4;   // dimensionless (in INaK)
extern double INaK;   // microA_per_microF (in INaK)
extern double JnakK;   // millimolar_per_millisecond (in INaK)
extern double JnakNa;   // millimolar_per_millisecond (in INaK)
extern double Knai;   // millimolar (in INaK)
extern double Knao;   // millimolar (in INaK)
extern double P;   // dimensionless (in INaK)
extern double Pnak;   // milliS_per_microF (in INaK)
extern double a1;   // dimensionless (in INaK)
extern double a2;   // dimensionless (in INaK)
extern double a3;   // dimensionless (in INaK)
extern double a4;   // dimensionless (in INaK)
extern double b1;   // dimensionless (in INaK)
extern double b2;   // dimensionless (in INaK)
extern double b3;   // dimensionless (in INaK)
extern double b4;   // dimensionless (in INaK)
extern double x1;   // dimensionless (in INaK)
extern double x2;   // dimensionless (in INaK)
extern double x3;   // dimensionless (in INaK)
extern double x4;   // dimensionless (in INaK)
extern double GNaL;   // milliS_per_microF (in INaL)
extern double INaL;   // microA_per_microF (in INaL)
extern double fINaLp;   // dimensionless (in INaL)
extern double hLss;   // dimensionless (in INaL)
extern double hLssp;   // dimensionless (in INaL)
extern double mLss;   // dimensionless (in INaL)
extern double thLp;   // millisecond (in INaL)
extern double tmL;   // millisecond (in INaL)
extern double INab;   // microA_per_microF (in INab)
extern double Ahs;   // dimensionless (in INa)
extern double INa;   // microA_per_microF (in INa)
extern double fINap;   // dimensionless (in INa)
extern double h;   // dimensionless (in INa)
extern double hp;   // dimensionless (in INa)
extern double hss;   // dimensionless (in INa)
extern double hssp;   // dimensionless (in INa)
extern double jss;   // dimensionless (in INa)
extern double mss;   // dimensionless (in INa)
extern double thf;   // millisecond (in INa)
extern double ths;   // millisecond (in INa)
extern double thsp;   // millisecond (in INa)
extern double tj;   // millisecond (in INa)
extern double tjp;   // millisecond (in INa)
extern double tm;   // millisecond (in INa)
extern double IpCa;   // microA_per_microF (in IpCa)
extern double AiF;   // dimensionless (in Ito)
extern double AiS;   // dimensionless (in Ito)
extern double Gto;   // milliS_per_microF (in Ito)
extern double Ito;   // microA_per_microF (in Ito)
extern double ass;   // dimensionless (in Ito)
extern double assp;   // dimensionless (in Ito)
extern double delta_epi;   // dimensionless (in Ito)
extern double dti_develop;   // dimensionless (in Ito)
extern double dti_recover;   // dimensionless (in Ito)
extern double fItop;   // dimensionless (in Ito)
extern double i;   // dimensionless (in Ito)
extern double ip;   // dimensionless (in Ito)
extern double iss;   // dimensionless (in Ito)
extern double ta;   // millisecond (in Ito)
extern double tiF;   // millisecond (in Ito)
extern double tiF_b;   // millisecond (in Ito)
extern double tiFp;   // millisecond (in Ito)
extern double tiS;   // millisecond (in Ito)
extern double tiS_b;   // millisecond (in Ito)
extern double tiSp;   // millisecond (in Ito)
extern double Jleak;   // millimolar_per_millisecond (in SERCA)
extern double Jup;   // millimolar_per_millisecond (in SERCA)
extern double Jupnp;   // millimolar_per_millisecond (in SERCA)
extern double Jupp;   // millimolar_per_millisecond (in SERCA)
extern double fJupp;   // dimensionless (in SERCA)
extern double upScale;   // dimensionless (in SERCA)
extern double Acap;   // centimeter_squared (in cell_geometry)
extern double Ageo;   // centimeter_squared (in cell_geometry)
extern double vcell;   // microliter (in cell_geometry)
extern double vjsr;   // microliter (in cell_geometry)
extern double vmyo;   // microliter (in cell_geometry)
extern double vnsr;   // microliter (in cell_geometry)
extern double vss;   // microliter (in cell_geometry)
extern double Jdiff;   // millimolar_per_millisecond (in diff)
extern double JdiffK;   // millimolar_per_millisecond (in diff)
extern double JdiffNa;   // millimolar_per_millisecond (in diff)
extern double Bcai;   // dimensionless (in intracellular_ions)
extern double Bcajsr;   // dimensionless (in intracellular_ions)
extern double Bcass;   // dimensionless (in intracellular_ions)
extern double cmdnmax;   // millimolar (in intracellular_ions)
extern double Istim;   // microA_per_microF (in membrane)
extern double vffrt;   // coulomb_per_mole (in membrane)
extern double vfrt;   // dimensionless (in membrane)
extern double EK;   // millivolt (in reversal_potentials)
extern double EKs;   // millivolt (in reversal_potentials)
extern double ENa;   // millivolt (in reversal_potentials)
extern double Jrel;   // millimolar_per_millisecond (in ryr)
extern double Jrel_inf;   // dimensionless (in ryr)
extern double Jrel_inf_temp;   // dimensionless (in ryr)
extern double Jrel_infp;   // dimensionless (in ryr)
extern double Jrel_temp;   // dimensionless (in ryr)
extern double a_rel;   // millisecond (in ryr)
extern double a_relp;   // millisecond (in ryr)
extern double btp;   // millisecond (in ryr)
extern double fJrelp;   // dimensionless (in ryr)
extern double tau_rel;   // millisecond (in ryr)
extern double tau_rel_temp;   // millisecond (in ryr)
extern double tau_relp;   // millisecond (in ryr)
extern double tau_relp_temp;   // millisecond (in ryr)
extern double Jtr;   // millimolar_per_millisecond (in trans_flux)

//------------------------------------------------------------------------------
// Procedures
//------------------------------------------------------------------------------

extern void init();
extern void compute(double time);

//------------------------------------------------------------------------------

#endif

//==============================================================================
// End of file
//==============================================================================
