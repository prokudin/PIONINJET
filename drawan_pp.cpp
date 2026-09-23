// Alexei Prokudin 21/12/2016
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
#include <approximate_collins_pp.h>
//#include <draw_sidis.h>

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
#include "reader_sidis.h"

// Cuba library
//#include "cuba.h"
#include <cuba_sidis.h>


using namespace sidis_namespace;

using namespace std;


bool Coridor            = false;

bool sidis_namespace::cuba_verbose = true; // iether printing or not verbose results of Cuba

double sidis_namespace::x_used, sidis_namespace::xf_used, sidis_namespace::z_used, sidis_namespace::pt_used, sidis_namespace::pt2_used, sidis_namespace::y_used, sidis_namespace::W_used, sidis_namespace::Q2_used;

TFile* sidis_namespace::vegas_root        = new TFile("vegas_root.root","recreate");

double sidis_namespace::xa_tree = 0.;
double sidis_namespace::xb_tree  = 0.;
double sidis_namespace::xf_tree  = 0.;
double sidis_namespace::zc_tree  = 0.;
double sidis_namespace::pt_tree  = 0.;
double sidis_namespace::jt_tree = 0.;
double sidis_namespace::y_tree = 0.;
TTree* sidis_namespace::tree = new TTree("ntuple","Ntuple");


int sidis_namespace::EXPERIMENT  = COMPASS_SIDIS;

int NDIM_JET = 3;
int NCOMP_JET = 1;

 
// Plot the results AN using gauss Holoborodko...
void draw_result_an(int npoints, const int ncomputations, const int ndim, int step_var, bool Coridor, TGraph* graph[]){
  sidis& SIDIS = *sidis::Instance();
  double x, x_min, x_max;
  double xF;

 
  double xff_min = -0.7;
  double xff_max =  0.7;


  double pt_min =  1.;
  double pt_max =  50.;


  double Draw_min = -0.2;
  double Draw_max = 0.;

  double asymmetry;


  double xx[npoints];
  double yy[npoints];

  double min[npoints]; // For Coridor
  double max[npoints];
  double corrx[2*npoints + 1]; // For Coridor
  double corry[2*npoints + 1];


  double numerator   = 0;
  double denominator = 0;

 
  for(int n = 0; n <= npoints; n++){
      yy[n] = 0.;
  }
  


  switch(EXPERIMENT){
  case RHIC_BRAHMS:
//     xff_min = 0.17;
//     xff_max = 0.35; 
    xff_min = 0.25;
    xff_max = 0.35; 
    pt_min  = 1.;
    pt_max  = 10.;
    Draw_min = -0.1;
    Draw_max = 0.1;
    break;
  case RHIC_STAR:
      xff_min = 0.16;
      xff_max = 0.55; 
//    xff_min = 0.22;
//    xff_max = 0.75; 
//    xff_min = 0.22;
//    xff_max = 0.65; 
//      xff_min = -0.6;
//      xff_max = -0.27; 
    pt_min  = 1.;
    pt_max  = 15.;
    Draw_min = -0.1;
    Draw_max = 0.15;
    break;
  case RHIC_ANDY:
    xff_min = 0.14;
    xff_max = 0.6; 
    pt_min  = 1.;
    pt_max  = 2.e+10;
    Draw_min = -0.1;
    Draw_max = 0.1;
    break;
   default:
    cout << "Define Experiment! DRAW RESULT SIDIS cpp" << endl;
    break;
  }


  switch(SIDIS.dependence){
   case xf_dependence:
    x_min = xff_min;
    x_max = xff_max;
    break;
   case pt_dependence:
    x_min = pt_min;
    x_max = pt_max;
    break;
     break;
  default:
    printf("Unknown dependence in draw results drawan_pp.cpp !\n");
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
   case xf_dependence:
    xF = x;

    SIDIS.effect = collins_effect; // numerator
    
    numerator = an_numerator_denominator_gauss( xF );

    SIDIS.effect = no_effect; // denominator
    
    denominator = an_numerator_denominator_gauss( xF );
    
    asymmetry = numerator/denominator;

    break;
   case pt_dependence:
    pt_used = x;

 
    SIDIS.effect = collins_effect; // numerator
    
    numerator = an_numerator_denominator_gauss_pt( pt_used );

    SIDIS.effect = no_effect; // denominator
    
    denominator = an_numerator_denominator_gauss_pt( pt_used );
    
    asymmetry = numerator/denominator;


    break;
   default:
    printf("Unknown dependence in drawan_pp.cpp !\n");
    break;
  }

    
    cout << "asymmetry = " << asymmetry << endl;

 
    yy[i] =  asymmetry; // asymmetry

    // ATTENTION ATTENTION ATTENTION
    // WORK ON IT -- Not effective!!! It should reset hoppet and then go other all points at once!
    // No resetting hoppet for every point !!!!!!
    if(Coridor) { // We draw also a coridor with errors WORK on IT

      //COLLINS
      if( SIDIS.calculation == collins_effect_pp ){
	for(int k = 0; k < set_number; k++){

	  SIDIS.Params.Collins.UseSet(k); // we do loop over all sets...
	  SIDIS.Params.Transversity.UseSet(k); // we do loop over all sets...
	  
      reset_hoppet();
      cout << k << endl;

  switch(SIDIS.dependence){
   case xf_dependence:
    xF = x;

    SIDIS.effect = collins_effect; // numerator
    
    numerator = an_numerator_denominator_gauss( xF );

    SIDIS.effect = no_effect; // denominator
    
    denominator = an_numerator_denominator_gauss( xF );
    
    asymmetry = numerator/denominator;

    break;
   case pt_dependence:
    pt_used = x;

 
    SIDIS.effect = collins_effect; // numerator
    
    numerator = an_numerator_denominator_gauss_pt( pt_used );

    SIDIS.effect = no_effect; // denominator
    
    denominator = an_numerator_denominator_gauss_pt( pt_used );
    
    asymmetry = numerator/denominator;


    break;
   default:
    printf("Unknown dependence in drawan_pp.cpp !\n");
    break;
  }

// WORK HERE!!!!   

    double sechenie =  asymmetry; // asymmetry
    
    cout << sechenie << endl;
	  
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
    reset_hoppet(); // DO NOT FORGET TO RESET HOPPET!
	
	
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
    grcor->SetFillColorAlpha(29,0.55);
    grcor->SetLineColor(1);
    grcor->SetLineWidth(0.1);

    graph[ncomputations+2] = grcor; // this will become a coridor...
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
  int npoints = 15; // npoints to draw...
  //int npoints = 2; // npoints to draw...
  

 
  set_hoppet(); // DON'T FORGET TO SET HOPPET

  // DATA RHIC_STAR 
//  EXPERIMENTAL_DATA RHIC_STAR_37("./data/pp_an_star2004.dat");
  EXPERIMENTAL_DATA RHIC_STAR_37("./data/pp_an_star3.7.dat");
  RHIC_STAR_37.X(1); // x
  RHIC_STAR_37.Y(2); // y
  RHIC_STAR_37.dY(3,4); // dy

  EXPERIMENTAL_DATA RHIC_STAR_33("./data/pp_an_star3.3.dat");
  RHIC_STAR_33.X(1); // x
  RHIC_STAR_33.Y(2); // y
  RHIC_STAR_33.dY(3,4); // dy
 
  EXPERIMENTAL_DATA RHIC_STAR_368("./data/pp_an_star3.68.dat");
  RHIC_STAR_368.X(1); // x
  RHIC_STAR_368.Y(2); // y
  RHIC_STAR_368.dY(3,4); // dy

  // DATA RHIC_BRAHMS 2.3grad
  EXPERIMENTAL_DATA RHIC_BRAHMS_PIM_23("./data/pp_an_brahms_pim_2.3.dat");
  RHIC_BRAHMS_PIM_23.X(1); // x
  RHIC_BRAHMS_PIM_23.Y(2); // y
  RHIC_BRAHMS_PIM_23.dY(3,4); // dy

  EXPERIMENTAL_DATA RHIC_BRAHMS_PIP_23("./data/pp_an_brahms_pip_2.3.dat");
  RHIC_BRAHMS_PIP_23.X(1); // x
  RHIC_BRAHMS_PIP_23.Y(2); // y
  RHIC_BRAHMS_PIP_23.dY(3,4); // dy
 

  // DATA RHIC_BRAHMS 4grad
  EXPERIMENTAL_DATA RHIC_BRAHMS_PIM_4("./data/pp_an_brahms_pim_4.dat");
  RHIC_BRAHMS_PIM_4.X(1); // x
  RHIC_BRAHMS_PIM_4.Y(2); // y
  RHIC_BRAHMS_PIM_4.dY(3,4); // dy

  EXPERIMENTAL_DATA RHIC_BRAHMS_PIP_4("./data/pp_an_brahms_pip_4.dat");
  RHIC_BRAHMS_PIP_4.X(1); // x
  RHIC_BRAHMS_PIP_4.Y(2); // y
  RHIC_BRAHMS_PIP_4.dY(3,4); // dy
  
  // data pt 500 gev
  EXPERIMENTAL_DATA RHIC_STAR_PI0_02("./data/pp_an_star500_pt_16_24_70mr.dat");
  RHIC_STAR_PI0_02.X(1); // x
  RHIC_STAR_PI0_02.Y(2); // y
  RHIC_STAR_PI0_02.dY(4); // dy
  RHIC_STAR_PI0_02.dX(3); // dy
  
  EXPERIMENTAL_DATA RHIC_STAR_PI0_028("./data/pp_an_star500_pt_24_32_70mr.dat");
  RHIC_STAR_PI0_028.X(1); // x
  RHIC_STAR_PI0_028.Y(2); // y
  RHIC_STAR_PI0_028.dY(4); // dy
  RHIC_STAR_PI0_028.dX(3); // dy


  EXPERIMENTAL_DATA RHIC_STAR_PI0_036("./data/pp_an_star500_pt_32_40_70mr.dat");
  RHIC_STAR_PI0_036.X(1); // x
  RHIC_STAR_PI0_036.Y(2); // y
  RHIC_STAR_PI0_036.dY(4); // dy
  RHIC_STAR_PI0_036.dX(3); // dy

  // ---  Initialize random generator
  srand (time(NULL));
 

  // Branches for the tree:
  tree->Branch("xa",&xa_tree,"xa_tree/D");
  tree->Branch("xb" ,&xb_tree, "xb_tree/D");
  tree->Branch("xf" ,&xf_tree, "xf_tree/D");
  tree->Branch("zc" ,&zc_tree, "zc_tree/D");
  tree->Branch("jt",&jt_tree,"jt_tree/D");
  tree->Branch("pt",&pt_tree,"pt_tree/D");
  tree->Branch("y",&pt_tree,"y_tree/D");

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

  if(strcmp(argv[7],"-4")== 0)
    {
      SIDIS.Average.y = 4.;
    }  
  if(strcmp(argv[7],"-2.3")== 0)
    {
      SIDIS.Average.y = 2.3;
    }  
  if(strcmp(argv[7],"-3.7")== 0)
    {
      SIDIS.Average.y = 3.7;
    }  
  if(strcmp(argv[7],"--3.7")== 0)
    {
      SIDIS.Average.y = -3.7;
    }  
  if(strcmp(argv[7],"-3.3")== 0)
    {
      SIDIS.Average.y = 3.3;
    }  
  if(strcmp(argv[7],"-3.68")== 0)
    {
      SIDIS.Average.y = 3.68;
    }  
    
    
  if( SIDIS.dependence == pt_dependence ){
    
    if(strcmp(argv[7],"-0.2")== 0)
    {
      SIDIS.Average.x = 0.2; // it is xF
    }  
    if(strcmp(argv[7],"-0.28")== 0)
    {
      SIDIS.Average.x = 0.28; // it is xF
    }  
      if(strcmp(argv[7],"-0.36")== 0)
    {
      SIDIS.Average.x = 0.36; // xF
    }  
      if(strcmp(argv[7],"-0.5")== 0)
    {
      SIDIS.Average.x = 0.5; // xF
    }  

    
    }
 
 
 
   double xmin = 0.0001;
   double xmax = 0.9999;
   double Q2min = 1.0000;
   double Q2max = 90.;

   int nd = 100;
   double stepx = (xmax - xmin)/nd;
   double stepQ2 = (Q2max - Q2min)/nd;
   double x = xmin, Q2 = Q2min;
// 
//   
// // write g1 in a file:
//    FILE *out1 = fopen( "g1.dat", "w" );
// 
//    for(int i = 0; i < nd; i++){
//    	x = xmin + stepx*i;
//    	for(int j = 0; j < nd; j++){
//     Q2 = Q2min + stepQ2*j;
//    
//     PARTONCONTENT test = SIDIS.helicity(x, Q2);
//     SIDIS.Unpolarised(SIDIS.Target, x, Q2);
//     
//     //cout <<  test.up/SIDIS.Target.partcontent.up << endl;
//  
//     fprintf(out1,"%8.6e  %8.6e  ",x, Q2);
//     fprintf(out1,"%8.6e  %8.6e ",test.up, test.down);
//     fprintf(out1,"%8.6e  %8.6e ",test.strange, test.anti_up);
//     fprintf(out1,"%8.6e  %8.6e\n",test.anti_down, test.anti_strange);
//     }
//    }
//    fclose(out1);
//    
//    // write g1 in a file:
//    out1 = fopen( "SofferBound.dat", "w" );
// 
//    for(int i = 0; i < nd; i++){
//    	x = xmin + stepx*i;
//    	for(int j = 0; j < nd; j++){
//     Q2 = Q2min + stepQ2*j;
//    
//     PARTONCONTENT test = SIDIS.SofferBound(x, Q2);
//  
//     fprintf(out1,"%8.6e  %8.6e  ",x, Q2);
//     fprintf(out1,"%8.6e  %8.6e ",test.up, -test.down);
//     fprintf(out1,"%8.6e  %8.6e ",test.strange, test.anti_up);
//     fprintf(out1,"%8.6e  %8.6e\n",-test.anti_down, test.anti_strange);
//     }
//    }
//    fclose(out1);
// 
// 
//    // write g1 in a file:
//    out1 = fopen( "fragmentationpiminus.dat", "w" );
// 
//    for(int i = 0; i < nd; i++){
//    	x = xmin + stepx*i;
//    	for(int j = 0; j < nd; j++){
//     Q2 = Q2min + stepQ2*j;
//    
//     SIDIS.Fragmentation(SIDIS.Produced,x, Q2);
//  
//     fprintf(out1,"%8.6e  %8.6e  ",x, Q2);
//     fprintf(out1,"%8.6e  %8.6e ",SIDIS.Produced.fragmentation.up, SIDIS.Produced.fragmentation.down);
//     fprintf(out1,"%8.6e  %8.6e ",SIDIS.Produced.fragmentation.strange, SIDIS.Produced.fragmentation.anti_up);
//     fprintf(out1,"%8.6e  %8.6e\n",SIDIS.Produced.fragmentation.anti_down, SIDIS.Produced.fragmentation.anti_strange);
//     }
//    }
//    fclose(out1);


   // // write g1 in a file:
   FILE* out1 = fopen( "pion_mrss.dat", "w" );

   for(int i = 0; i < nd; i++){
   	x = xmin + stepx*i;
   	for(int j = 0; j < nd; j++){
    Q2 = Q2min + stepQ2*j;
    PARTONCONTENT test = SIDIS.unpolarised_pion(x, Q2);
    
 
    fprintf(out1,"%8.6e  %8.6e  ",x, Q2);
    fprintf(out1,"%8.6e  %8.6e ",test.up, test.down);
    fprintf(out1,"%8.6e  %8.6e ",test.strange, test.anti_up);
    fprintf(out1,"%8.6e  %8.6e\n",test.anti_down, test.anti_strange);
   
     }
   }
   fclose(out1);
   


   hold(true);

 

//  SIDIS.calculation = sivers_effect_pp; // we want to calculate Sivers effect in PP collision!!!
  SIDIS.calculation = collins_effect_pp; // we want to calculate Collins effect in PP collision!!!
//  SIDIS.calculation = sum_effect_pp; // we want to calculate Sivers+Collins effect in PP collision!!!
  int ncomputations = 1; // this depend on how many effects we want to sum
                         // for BM effect for example: BM + CAHN -> NCOMPUTATIONS = 2
 
  int ndim = 2; // if  xa and xb

 
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
  //gStyle->Reset();
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
    
    
  t->DrawLatex(0.08,0.6,"A_{N}");
    
  t->SetTextAngle(0);
  t->SetTextSize(0.05);
  switch( SIDIS.dependence ){
  case xf_dependence:
    t->DrawLatex(0.85,0.05,"x_{F}");
    break;
  case pt_dependence:
    t->DrawLatex(0.85,0.05,"P_{hT} (GeV)");
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


//   TGraph* graph[ncomputations + 1];
  TGraph* graph[ncomputations + 2]; // the last graph is for Coridor drawing!!!!

  set_hoppet(); // DON'T FORGET TO SET HOPPET

  //draw_result_sidis(npoints, ncomputations, ndim, step, Coridor, Integrand_PP, graph);
  
  draw_result_an( npoints, ncomputations, ndim, step, Coridor, graph);



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
  
  int red = 2;

  if( EXPERIMENT == RHIC_STAR ) {

    if( SIDIS.Beam.GetHadronType() == PROTON && SIDIS.Produced.GetHadronType() == PION) {
      switch( SIDIS.dependence ){
      case pt_dependence:
        	if( SIDIS.Produced.GetCharge() == 0. && SIDIS.Average.x == 0.2 ) RHIC_STAR_PI0_02.plot()->Draw("PZS"); // Draw the data
        	if( SIDIS.Produced.GetCharge() == 0. && SIDIS.Average.x == 0.28 ) RHIC_STAR_PI0_028.plot()->Draw("PZS"); // Draw the data
        	if( SIDIS.Produced.GetCharge() == 0. && SIDIS.Average.x == 0.36 ) RHIC_STAR_PI0_036.plot()->Draw("PZS"); // Draw the data
      case x_dependence:
      case z_dependence:
      case y_dependence:
 	break;
      case xf_dependence:
  	if( SIDIS.Produced.GetCharge() == 0. && SIDIS.Average.y == 3.7 ) RHIC_STAR_37.plot()->Draw("PZS"); // Draw the data
  	if( SIDIS.Produced.GetCharge() == 0. && SIDIS.Average.y == 3.3 ) RHIC_STAR_33.plot()->Draw("PZS"); // Draw the data
  	if( SIDIS.Produced.GetCharge() == 0. && SIDIS.Average.y == 3.68 ) RHIC_STAR_368.plot()->Draw("PZS"); // Draw the data
	break;
       default:
	break;
      }
    }
}

 if( EXPERIMENT == RHIC_BRAHMS ) {
    if( SIDIS.Beam.GetHadronType() == PROTON && SIDIS.Produced.GetHadronType() == PION) {
      switch( SIDIS.dependence ){
      case pt_dependence:
      case x_dependence:
      case z_dependence:
      case y_dependence:
 	break;
      case xf_dependence:
  	if( SIDIS.Produced.GetCharge() > 0. && SIDIS.Average.y == 2.3 ) RHIC_BRAHMS_PIP_23.plot()->Draw("PZS"); // Draw the data
  	if( SIDIS.Produced.GetCharge() < 0. && SIDIS.Average.y == 2.3 ) RHIC_BRAHMS_PIM_23.plot()->Draw("PZS"); // Draw the data
  	if( SIDIS.Produced.GetCharge() > 0. && SIDIS.Average.y == 4. ) RHIC_BRAHMS_PIP_4.plot()->Draw("PZS"); // Draw the data
  	if( SIDIS.Produced.GetCharge() < 0. && SIDIS.Average.y == 4. ) RHIC_BRAHMS_PIM_4.plot()->Draw("PZS"); // Draw the data
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
  TFile* file_reader = new TFile("./vegas_root.root");
  reader_sidis( file_reader  );

  //Now we draw functions...
  //double Q2 = 2.41;
  //SIDIS.DrawSiversDistribution(SIDIS.Target, Q2);
  //SIDIS.DrawSiversDistribution(SIDIS.Target, Q2, Coridor);
  //SIDIS.DrawF1TDistribution(SIDIS.Target, Q2, Coridor);
  //SIDIS.DrawF1TDistributionBound(SIDIS.Target, Q2, Coridor);
  //SIDIS.DrawSiversDistributionBound(SIDIS.Target, Q2, Coridor);


  //SIDIS.DrawSivers3DDistribution(SIDIS.Target, Q2);

  //SIDIS.DrawSivers3DDistribution(SIDIS.Target, Q2, Coridor);
  //SIDIS.DrawF1T3DDistribution(SIDIS.Target, Q2);
  //SIDIS.DrawF1T3DDistributionSum(SIDIS.Target, Q2);

  //SIDIS.DrawFragmentationFunctions(SIDIS.Produced, Q2);

  //SIDIS.Params.Sivers.Print();

  theApp->Run(true);
  canvas1->Close();




  return 0; 
}
