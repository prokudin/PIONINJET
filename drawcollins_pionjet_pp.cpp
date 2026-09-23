// ./drawcollins_pionjet_sidis.cpp
// Draw results of pion in jet
// Alexei Prokudin 13/4/2016
//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cfortran.h>
#include <stdio.h>
#include <sidis.h>
#include <hadron.h>
#include <common_sidis.h>
#include <setexperiment_sidis.h>
#include <draw_result_sidis.h>
#include <experimental_data.h>
//#include <draw_sidis.h>
#include <pionjet_pp.h>

// --- This is for root use in C++
#include "TApplication.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TPaveLabel.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraph.h"
#include "TFile.h"
#include "TObject.h"
#include "TPaveStats.h"
#include "TText.h"
//#include <gsl/gsl_integration.h>
#include "reader_sidis.h"

// Cuba library
//#include "cuba.h"
#include <cuba_sidis.h>


using namespace sidis_namespace;

using namespace std;


bool Coridor            = false;

bool sidis_namespace::cuba_verbose = false; // iether printing or not verbose results of Cuba

bool NLO_JET = false; // ANSELMINO 13
//bool NLO_JET = true; // KANG 15

double sidis_namespace::x_used, sidis_namespace::xf_used, sidis_namespace::z_used, sidis_namespace::pt_used, sidis_namespace::pt2_used, sidis_namespace::y_used, sidis_namespace::W_used, sidis_namespace::Q2_used;

TFile* sidis_namespace::vegas_root        = new TFile("vegas_root_sidis.root","recreate");

double sidis_namespace::xa_tree = 0.;
double sidis_namespace::xb_tree  = 0.;
double sidis_namespace::xf_tree  = 0.;
double sidis_namespace::zc_tree  = 0.;
double sidis_namespace::pt_tree  = 0.;
double sidis_namespace::jt_tree = 0.;
double sidis_namespace::y_tree = 0.;
TTree* sidis_namespace::tree = new TTree("ntuple","Ntuple");
//

int sidis_namespace::EXPERIMENT  = COMPASS_SIDIS;



int NDIM_JET = 3;
int NCOMP_JET = 1;


// Plot the results pion in jet using VEGAS
void draw_result_pionjet(int npoints, const int ncomputations, const int ndim, int step_var, bool Coridor, TGraph* graph[], double zc, double pt, double jt,  double y_min, double y_max,  double jt_min, double jt_max){
  sidis& SIDIS = *sidis::Instance();
  double x, x_min, x_max;


//  double zc_min_draw =  0.1;
//  double zc_max_draw =  0.7; // main
  double zc_min_draw =  0.35;
  double zc_max_draw =  0.85; // Star 2020 predictions


  double jt_min_draw =  0.1;
  double jt_max_draw =  1.;

  double pt_min_draw =  1.;
  double pt_max_draw =  30.;


  double Draw_min = -0.2;
  double Draw_max = 0.;

  double asymmetry;


  double xx[npoints];
  double yy[npoints];

  double min[npoints]; // For Coridor
  double max[npoints];
  double corrx[2*npoints + 1]; // For Coridor
  double corry[2*npoints + 1];


  NCOMP_JET = ncomputations + 1; // Important!!! Check the subroutine  Integrand_SIDIS ncomputations defines how many effects we want to sum + 1 = also unpolarised
  NDIM_JET  = ndim;

  int       neval, fail;
  double integral[NCOMP_JET], error[NCOMP_JET], prob[NCOMP_JET];
  double numerator[NCOMP_JET];
  double denominator = 0;


  for(int n = 0; n <= npoints; n++){
      yy[n] = 0.;
  }



  switch(EXPERIMENT){
   case RHIC_STAR:
    Draw_min = -0.08;
    Draw_max = 0.08;
    break;
  default:
    cout << "Define Experiment! DRAW RESULT PION JET" << endl;
    break;
  }


  switch(SIDIS.dependence){
   case z_dependence:
    x_min = zc_min_draw;
    x_max = zc_max_draw;
    break;
   case pt_dependence:
    x_min = pt_min_draw;
    x_max = pt_max_draw;
    break;
   case jt_dependence: // it is jt here
    x_min = jt_min_draw;
    x_max = jt_max_draw;
    break;
   default:
    printf("Unknown dependence in draw results approximate_pionjet_pp.cpp !\n");
    break;
  }

  double x_step  = (x_max - x_min)/(npoints - 1);

  if( step_var == step_log) { // log scale
    x_step  = (log( x_max ) - log( x_min ))/float(npoints - 1);
  }

  for(int i = 0; i < npoints; ++i){

    if(i == 0)
      x = x_min;
    else {
      x += x_step;
      if( step_var == step_log)  // log scale
	x = exp( log( x_min ) + float(i) * x_step );
    }

    xx[i] = x;

    cout << "x = " << x << endl;

  switch(SIDIS.dependence){
   case z_dependence:
    z_used = x;
    if (!NLO_JET) asymmetry = calc_pionjet( z_used, pt, jt, y_min, y_max, jt_min, jt_max);
    if (NLO_JET) asymmetry = calc_pionjet_NLO( z_used, pt, jt, y_min, y_max, jt_min, jt_max);
    break;
   case pt_dependence:
    pt_used = x;
    if (!NLO_JET) asymmetry = calc_pionjet( zc, pt_used, jt, y_min, y_max, jt_min, jt_max);
    if (NLO_JET) asymmetry = calc_pionjet_NLO( zc, pt_used, jt, y_min, y_max, jt_min, jt_max);
    break;
   case jt_dependence: // it is jt here
    W_used = x;
    if (!NLO_JET) asymmetry = calc_pionjet( zc, pt, W_used, y_min, y_max, jt_min, jt_max);
    if (NLO_JET) asymmetry = calc_pionjet_NLO( zc, pt, W_used, y_min, y_max, jt_min, jt_max);
    break;
   default:
    printf("Unknown dependence in draw results approximate_pionjet_pp.cpp !\n");
    break;
  }


    cout << "asymmetry = " << asymmetry << endl;


    yy[i] =  asymmetry; // asymmetry


    // WORK ON IT !!!!!!
    if(Coridor) { // We draw also a coridor with errors WORK on IT

      //COLLINS
      if( SIDIS.calculation == collins_effect_pp ){
	for(int k = 0; k < set_number; k++){

	  SIDIS.Params.Collins.UseSet(k); // we do loop over all sets...
	  SIDIS.Params.Transversity.UseSet(k); // we do loop over all sets...


//     SIDIS.effect = no_effect;
//     double denominator = calc_pionjet( zc, pt, jt, y_min, y_max, jt_min, jt_max);
//
//
//
//     //Unpolarized
//     SIDIS.effect = sivers_effect;
//     double numerator = calc_pionjet( zc, pt, jt, y_min, y_max, jt_min, jt_max);
//
//
//     double asymmetry = numerator/denominator; // asymmetry
    if (NLO_JET) reset_hoppet();

    switch(SIDIS.dependence){
   case z_dependence:
    z_used = x;
    if (!NLO_JET) asymmetry = calc_pionjet( z_used, pt, jt, y_min, y_max, jt_min, jt_max);
    if (NLO_JET) asymmetry = calc_pionjet_NLO( z_used, pt, jt, y_min, y_max, jt_min, jt_max);
    break;
   case pt_dependence:
    pt_used = x;
    if (!NLO_JET) asymmetry = calc_pionjet( zc, pt_used, jt, y_min, y_max, jt_min, jt_max);
    if (NLO_JET) asymmetry = calc_pionjet_NLO( zc, pt_used, jt, y_min, y_max, jt_min, jt_max);
    break;
   case jt_dependence: // it is jt here
    W_used = x;
    if (!NLO_JET) asymmetry = calc_pionjet( zc, pt, W_used, y_min, y_max, jt_min, jt_max);
    if (NLO_JET) asymmetry = calc_pionjet_NLO( zc, pt, W_used, y_min, y_max, jt_min, jt_max);
    break;
   default:
    printf("Unknown dependence in draw results approximate_pionjet_pp.cpp !\n");
    break;
  }


    double sechenie =  asymmetry; // asymmetry

	  if( k == 0 ){
	    min[i] = sechenie;
	    max[i] = sechenie;
	  } else {
	    if( (sechenie ) > max[i] ){
	      max[i] = sechenie;
	    }
	    if(  sechenie < min[i] ){
	      min[i] = sechenie;
	    }
	  }
	}
	SIDIS.Params.Collins.UseMainSet(); // main set is restored after the loop ...
	SIDIS.Params.Transversity.UseMainSet(); // main set is restored after the loop ...
	    if (NLO_JET) reset_hoppet();

      }
    }


  }



  //Draw Coridor if needed
  if(Coridor){
    for (int i = 0; i < npoints; ++i) {
      corry[i] = min[i];
      corrx[i] = xx[i];
    }

    for (int i = 0; i < npoints; ++i) {
      corry[npoints + i] = max[npoints - i - 1];
      corrx[npoints + i] = xx[npoints - i - 1];
    }

    corry[2*npoints] = min[0];
    corrx[2*npoints] = xx[0];

    TGraph* grcor;
    grcor = new TGraph(2*npoints + 1,corrx,corry);
    grcor->SetFillColor(29);
    grcor->SetLineColor(1);
    grcor->SetLineWidth(0.1);

    graph[ncomputations + 2] = grcor; // this will become a coridor...
  }



    graph[0] = new TGraph(npoints, xx, yy);
    graph[0]->SetLineColor(2);
    graph[0]->SetLineWidth(3);
    graph[0]->SetTitle("");
    graph[0]->SetMinimum(Draw_min);
    graph[0]->SetMaximum(Draw_max);
}






//========================================================= main
int main(int argc, char **argv)
{
  sidis& SIDIS = *sidis::Instance();
  char* file_data   = " ";
  ofstream out;
  int read;
  //int npoints = 30; // npoints to draw...
  int npoints = 7; // npoints to draw...
  //int npoints = 15; // npoints to draw...

//   double zc = 0.37;
//   //double zc = 0.23;
//   //double zc = 0.13;
//   //double pt = 12.9;
//   //double pt = 31.;
//   double pt = 15.; // for JIM
//   double jt = 1.;
//   double y_min = 0.;
//   double y_max = 1.;
//   //double y_min = -1.;
//   //double y_max = 0.;
//   double jt_min = 0.0;
//   double jt_max = 3.;

  // PREDICTIONS FOR RENEE, 200 GeV paper
  double zc = 0.2;
  double pt = 12.46; // for RENEE
  double jt = 0.33;
  double y_min = 0.;
  double y_max = 1.;
  double jt_min = 0.05;
  double jt_max = 4.5;

  // DATA RHIC_STAR
  EXPERIMENTAL_DATA RHIC_STAR_Z_PIP200("./data/pp_jetpip_star200.dat");
  RHIC_STAR_Z_PIP200.X(1); // x
  RHIC_STAR_Z_PIP200.Y(3); // y
  RHIC_STAR_Z_PIP200.dY(4,5); // dy

  EXPERIMENTAL_DATA RHIC_STAR_Z_PIP500("./data/pp_jetpip_star500.dat");
  RHIC_STAR_Z_PIP500.X(1); // x
  RHIC_STAR_Z_PIP500.Y(3); // y
  RHIC_STAR_Z_PIP500.dY(4,5); // dy


  EXPERIMENTAL_DATA RHIC_STAR_JT_013_PIP500("./data/pp_jetpip_jt_z013_star500.dat");
  RHIC_STAR_JT_013_PIP500.X(1); // x
  RHIC_STAR_JT_013_PIP500.Y(3); // y
  RHIC_STAR_JT_013_PIP500.dY(4,5); // dy


  EXPERIMENTAL_DATA RHIC_STAR_Z_PIM200("./data/pp_jetpim_star200.dat");
  RHIC_STAR_Z_PIM200.X(1); // x
  RHIC_STAR_Z_PIM200.Y(3); // y
  RHIC_STAR_Z_PIM200.dY(4,5); // dy

  EXPERIMENTAL_DATA RHIC_STAR_Z_PIM500("./data/pp_jetpim_star500.dat");
  RHIC_STAR_Z_PIM500.X(1); // x
  RHIC_STAR_Z_PIM500.Y(3); // y
  RHIC_STAR_Z_PIM500.dY(4,5); // dy

  EXPERIMENTAL_DATA RHIC_STAR_JT_013_PIM500("./data/pp_jetpim_jt_z013_star500.dat");
  RHIC_STAR_JT_013_PIM500.X(1); // x
  RHIC_STAR_JT_013_PIM500.Y(3); // y
  RHIC_STAR_JT_013_PIM500.dY(4,5); // dy



  EXPERIMENTAL_DATA RHIC_STAR_JT_023_PIP500("./data/pp_jetpip_jt_z023_star500.dat");
  RHIC_STAR_JT_023_PIP500.X(1); // x
  RHIC_STAR_JT_023_PIP500.Y(3); // y
  RHIC_STAR_JT_023_PIP500.dY(4,5); // dy

  EXPERIMENTAL_DATA RHIC_STAR_JT_023_PIM500("./data/pp_jetpim_jt_z023_star500.dat");
  RHIC_STAR_JT_023_PIM500.X(1); // x
  RHIC_STAR_JT_023_PIM500.Y(3); // y
  RHIC_STAR_JT_023_PIM500.dY(4,5); // dy


  EXPERIMENTAL_DATA RHIC_STAR_JT_037_PIP500("./data/pp_jetpip_jt_z037_star500.dat");
  RHIC_STAR_JT_037_PIP500.X(1); // x
  RHIC_STAR_JT_037_PIP500.Y(3); // y
  RHIC_STAR_JT_037_PIP500.dY(4,5); // dy

  EXPERIMENTAL_DATA RHIC_STAR_JT_037_PIM500("./data/pp_jetpim_jt_z037_star500.dat");
  RHIC_STAR_JT_037_PIM500.X(1); // x
  RHIC_STAR_JT_037_PIM500.Y(3); // y
  RHIC_STAR_JT_037_PIM500.dY(4,5); // dy


  // ---  Initialize random generator
  srand (time(NULL));

  set_hoppet(); // DON'T FORGET TO SET HOPPET

  set_hoppet(); // DON'T FORGET TO SET HOPPET


  // Branches for the tree:
  tree->Branch("xa",&xa_tree,"xa_tree/D");
  tree->Branch("xb" ,&xb_tree, "xb_tree/D");
  tree->Branch("xf" ,&xf_tree, "xf_tree/D");
  tree->Branch("zc" ,&zc_tree, "zc_tree/D");
  tree->Branch("jt",&jt_tree,"jt_tree/D");
  tree->Branch("pt",&pt_tree,"pt_tree/D");
  tree->Branch("y",&y_tree,"y_tree/D");
  // Choose experiment name!!!


  // --- Settings for EXPERIMENT
  if( SetExperiment_SIDIS(argv[1]) != 0) // EXPERIMENT
    {
      cerr << "Error in SetExperiment_SIDIS (" << argv[1] << ")" << endl;
      return 1;
    }


  if( SetBeam_SIDIS(argv[2]) != 0)  // Beam
    {
      cerr << "Error in SetBeam_SIDIS (" << argv[2] << ")" << endl;
      return 1;
    }

  if( SetTarget_SIDIS(argv[3]) != 0)  // TARGET
    {
      cerr << "Error in SetTarget_SIDIS (" << argv[3] << ")" << endl;
      return 1;
    }

  if( SetHadron_SIDIS(argv[4]) != 0) // PRODUCED HADRON
    {
      cerr << "Error in SetHadron_SIDIS (" << argv[4] << ")" << endl;
      return 1;
    }

    if( SetDependence_SIDIS(argv[5]) != 0 ) // DEPENDENCE
      {
	cerr << "Error in SetDependence_SIDIS (" << argv[5] << ")" << endl;
	return 1;
      }

  if(strcmp(argv[6],"-CORIDOR")== 0)
    {
      Coridor = true;
    } else {
      Coridor = false;
    }



  SIDIS.calculation = collins_effect_pp; // we want to calculate Collins effect in PP collision!!!

  int ncomputations = 1;

  int ndim = 3; // if  z, pt
  if (SIDIS.dependence == jt_dependence) ndim = 2; //jt dependence

  if (NLO_JET) ndim = 5; // if  z, pt at NLO -> bt and zhat integration
  if (SIDIS.dependence == jt_dependence && NLO_JET) ndim = 4; //jt dependence


//   NCOMP_JET = 2;
//   NDIM_JET = 3;
//   z_used = zc;
//
//   double asymmetry = calc_pionjet1( );


//  double asymmetry = calc_pionjet( zc, pt, jt, y_min, y_max, jt_min, jt_max);

//  hold(true);
// --- Root
  TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetTitleColor(0);
  plain->SetStatColor(0);
  plain->SetPalette(1);


  gROOT->SetStyle("Plain");
  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132,"xyz");

  TApplication* theApp = new TApplication("App", &argc, argv);
  TCanvas* canvas1 = new TCanvas("result","result");
  canvas1->SetFillColor(0);
  canvas1->SetFrameFillStyle(0);
  canvas1->SetFillStyle(0);
  canvas1->Draw();



  TLatex *t = new TLatex();
  t->SetTextAngle(90);
  t->SetTextSize(0.08);


  t->DrawLatex(0.08,0.6,"A_{UT}^{sin(#phi_{S} - #phi_{H})}");

  t->SetTextAngle(0);
  t->SetTextSize(0.05);
  switch( SIDIS.dependence ){
  case z_dependence:
    t->DrawLatex(0.85,0.05,"z");
    break;
  case pt_dependence:
    t->DrawLatex(0.85,0.05,"P_{T}^{jet} (GeV)");
    break;
  case jt_dependence:
    t->DrawLatex(0.85,0.05,"j_{T} (GeV)");
    break;
  default:
    break;
  }

  TPad* tpad = new TPad("nemo","This is nemo",0.1,0.1  ,1.,1.,0,0,0);
  tpad->Range(0,0,10,10);
  tpad->SetFillColor(0);
  tpad->SetFrameFillStyle(0);
  tpad->SetFillStyle(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameLineWidth(0);
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetLabelOffset(0.04,"x");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetNdivisions(8,"xy");

  Float_t small = 0;

  gPad->SetFillColor(0);


  tpad->Draw();
  tpad->cd();

  int step = step_lin;
  if(SIDIS.dependence == jt_dependence) step = step_log;


//   TGraph* graph[ncomputations + 1];
  TGraph* graph[ncomputations + 2]; // the last graph is for Coridor drawing!!!!

  draw_result_pionjet(npoints, ncomputations, ndim, step, Coridor, graph, zc, pt, jt, y_min, y_max, jt_min, jt_max);



  graph[0]->SetLineStyle(1);
  graph[0]->Draw("ACP"); // ? STRANGE LOOK IN BOER MULDERS>>>>>

  //CORIDOR
  if(Coridor){
    graph[ncomputations + 2]->Draw("F");
    graph[ncomputations + 2]->Draw("L");
  }

  graph[0]->SetLineStyle(1);
  graph[0]->Draw("CP"); // WORK ON IT !!!!!!!!!!!!!!!!!!!!!!!!!!


  double minimum_x = graph[0]->GetXaxis()->GetXmin();
  double maximum_x = graph[0]->GetXaxis()->GetXmax();

  TLine *line = new TLine(minimum_x, 0., maximum_x, 0.);
  line->Draw();
  //Description
  SIDIS.description()->Draw();

    if(  SIDIS.dependence == jt_dependence ){
    tpad->SetLogx();
    tpad->Update();
  }


  int red = 2;

  if( EXPERIMENT == RHIC_STAR ) {
    if( SIDIS.Target.GetHadronType() == PROTON && SIDIS.Produced.GetHadronType() == PION) {
      switch( SIDIS.dependence ){
      case pt_dependence:
      case x_dependence:
      case xf_dependence:
      case y_dependence:
 	break;
      case z_dependence:
  	  if( SIDIS.Produced.GetCharge() == +1. && SIDIS.GetEnergy() == 200.*200.) RHIC_STAR_Z_PIP200.plot()->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == +1. && SIDIS.GetEnergy() == 500.*500.) RHIC_STAR_Z_PIP500.plot()->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == -1. && SIDIS.GetEnergy() == 200.*200.) RHIC_STAR_Z_PIM200.plot(4,24)->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == -1. && SIDIS.GetEnergy() == 500.*500.) RHIC_STAR_Z_PIM500.plot(4,24)->Draw("PZS"); // Draw the data
	break;
      case jt_dependence:
  	  if( SIDIS.Produced.GetCharge() == +1. && SIDIS.GetEnergy() == 500.*500. && zc == 0.13) RHIC_STAR_JT_013_PIP500.plot()->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == -1. && SIDIS.GetEnergy() == 500.*500. && zc == 0.13) RHIC_STAR_JT_013_PIM500.plot(4,24)->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == +1. && SIDIS.GetEnergy() == 500.*500. && zc == 0.23) RHIC_STAR_JT_023_PIP500.plot()->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == -1. && SIDIS.GetEnergy() == 500.*500. && zc == 0.23) RHIC_STAR_JT_023_PIM500.plot(4,24)->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == +1. && SIDIS.GetEnergy() == 500.*500. && zc == 0.37) RHIC_STAR_JT_037_PIP500.plot()->Draw("PZS"); // Draw the data
  	  if( SIDIS.Produced.GetCharge() == -1. && SIDIS.GetEnergy() == 500.*500. && zc == 0.37) RHIC_STAR_JT_037_PIM500.plot(4,24)->Draw("PZS"); // Draw the data
	break;
       default:
	break;
      }
    }
}


  tree->Write();
  vegas_root->Close();


  canvas1->Update();


  // Now the kinematics...
  TFile* file_reader = new TFile("./vegas_root_sidis.root");
  //reader_sidis( file_reader  );

  //Now we draw functions...
  double Q2 = 2.41;
   //SIDIS.DrawTransversityDistribution(SIDIS.Target, Q2, Coridor);
  // Draw Collins Fragmentation as well...
  //Q2 = 2.41;
  //SIDIS.DrawCollinsDistribution(SIDIS.Produced, Q2, Coridor);

  theApp->Run(true);
  canvas1->Close();




  return 0;
}
