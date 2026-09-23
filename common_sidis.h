/* /home/prokudin/cpp/collins08 my headers for fitaut */
/* Alexei Prokudin 31/10/2008        */

#ifndef __COMMON_SIDIS_H__
#define __COMMON_SIDIS_H__


#include "TFile.h"
#include "TTree.h"

using namespace std;

namespace sidis_namespace {

  extern double x_used, z_used, pt_used, y_used, xf_used, W_used, Q2_used, pt2_used;


const int COMPASS_SIDIS    = 1;
const int HERMES           = 2;
const int HERMES_LP        = 3;
const int HERMES_UNP       = 4;
const int COMPASS_SIDIS_UNP= 5;
const int EIC_LP           = 6;
const int MEIC1             = 7; // EIC JLAB Design min energy 400 GeV^2
const int MEIC2             = 8; // EIC JLAB Design max energy 1200 GeV^2
const int JLAB12            = 9;
const int HERMES09          = 10; // HERMES cuts 2009 for Sivers...
const int E143              = 11; // SlAC E143 experiment g2 measurement...
const int JLAB6             = 12;
const int E155              = 13; // SlAC E155 experiment g2 measurement...
const int SANE              = 14; // JLAB SANE experiment g2 measurement...
const int BELLE             = 15; // BELLE e+e- experiment collins measurement...
const int RHIC_STAR               = 16; //  RHIC PP
const int RHIC_BRAHMS             = 17; //  RHIC PP
const int BELLE12             = 18; // BELLE e+e- experiment collins measurement sin(phi1+phi2) method...
const int ZEUS               = 19; // ZEUS unpolarised cross sections...
const int EMC                = 20; // EMC unpolarised cross sections...
const int RHIC_ANDY             = 21; //  RHIC PP
const int RHIC_PHENIX             = 22; //  RHIC PP
const int BABAR              = 23; //  BABAR
extern int EXPERIMENT;

extern bool cuba_verbose; // iether printing or not verbose results of Cuba



extern TFile* vegas_root;

extern TTree* tree; // we will produce a tree with variables...

extern double xa_tree;
extern double xb_tree;
extern double xf_tree;
extern double zc_tree;
extern double jt_tree;
extern double pt_tree;
extern double y_tree;




//const bool VALENCE_SIDIS = true; // NOTE used only in Gamberg:2013kla parametrizations!

const bool VALENCE_SIDIS = false; // NOTE used in most other parametrizations

// For SIVERS x, z, pt dependence...
enum DEPENDENCE
   {
     x_dependence = 1,
     z_dependence,
     y_dependence,
     pt_dependence,
     m_dependence,
     xf_dependence,
     kinematical_dependence, // to integrate all
     W_dependence,
     Q2_dependence,
     pt2_dependence, // dependence on pt^2
     jt_dependence // dependence on jt of the pion in jet
   };



// Integration limits
struct LIMITS {
  //Experimental cuts

  double xF_min;
  double xF_max;

  bool xF_cut_used;

  double x_min;
  double x_max;

  bool x_cut_used;

  double y_min;
  double y_max;

  bool y_cut_used;

  double z_min;
  double z_max;

  bool z_cut_used;

  double pt_min;
  double pt_max;

  bool pt_cut_used;

  double kt_min;
  double kt_max;

  bool kt_cut_used;

  double Q2_min;
  double Q2_max;

  bool Q2_cut_used;

  double W2_min;
  double W2_max;

  bool W2_cut_used;

  double Eh_min;
  double Eh_max;

  bool Eh_cut_used;

};





// Switch to sivers, collins effects......
enum EFFECT {
  sivers_effect = 1,
  cahn_effect,
  boermulders_effect,
  collins_effect,
  pretzelosity_effect,
  g1t_effect_weighted_ww, // A_LT with WW approximation PT weighted
  g1t_effect_ww, // A_LT with WW approximation
  g1t_effect, // A_LT
  lp_full,    // full calculation l P --> PI X
  lpjet_full, // full calculation l P --> JET X
  pt_average, // average PT calculation for boer mulders etc...
  g1_effect, // A_LL
  g2_structure,
  h1l_effect_weighted_ww, // A_UL^sin(2 Phi) with WW approximation PT weighted
  h1l_effect_ww, // A_UL^sin(2 Phi) with WW approximation
  h1l_effect, // A_UL^sin(2 Phi)
  no_effect,
  sivers_effect_pp, // twist-3 formalizm in pp
  collins_effect_pp, // twist-3 formalizm in pp
  collins_type_effect, // for Andreas Metz --> H\hat contributions
  collins_type_effect1, // for Andreas Metz --> HUF  contributions
  sum_effect_pp, // twist-3 formalizm in pp
  sivers_effect_nlo, // NLO TMD evolution with Mert Aybat and Ted Rogers
  normalization // this one to calculate normalization for DIS cross section unpolarised
};






}

#endif //#ifndef  __COMMON_SIDIS_H__
