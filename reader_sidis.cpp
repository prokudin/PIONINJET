#include <TROOT.h>
#include <TStyle.h>
#include <TObject.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>

#include <TMath.h>
#include <TLorentzVector.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TF2.h>

#include <TProfile.h>

#include <TLatex.h>
#include <TLegend.h> 
#include <TLegendEntry.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TMinuit.h> 
#include <TApplication.h>

#include <TRandom3.h>


#include <TColor.h>


#include <TApplication.h>
#include <TMinuit.h> 

#include <common_sidis.h>
#include "reader_sidis.h"

using namespace sidis_namespace;

int reader_sidis( TFile *file ) 
{



  TStyle *plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetTitleColor(1);
  plain->SetStatColor(0);
  plain->SetPalette(1);
  gStyle->SetOptStat(1); // Statistics?! No -> 0
   

  gROOT->SetStyle("Plain");
  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132,"xyz");

  gStyle->SetTitleSize(0.09,"xy");
  gStyle->SetLabelSize(0.06,"xy");
  gStyle->SetLabelOffset(0.02,"x");
  gStyle->SetLabelOffset(0.01,"y"); 
  gStyle->SetNdivisions(8,"xy");


//   // COLOR PALETTES!!!!
//   UInt_t Number = 3;
//   Double_t Red[3]    = { 1.00, 0.00, 0.00};
//   Double_t Green[3]  = { 0.00, 1.00, 0.00};
//   Double_t Blue[3]   = { 1.00, 0.00, 1.00};
//   Double_t Length[3] = { 0.00, 0.00, 1.00 };
//   Int_t nb=50;
//   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
//   //  f2->SetContour(nb);


    const Int_t NRGBs = 5;
//     const Int_t NCont = 255;
//    const Int_t NCont = 25;
   const Int_t NCont = 25;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
//     Double_t stops[NRGBs] = { 0.00, 0.10, 0.16, 0.84, 1.00 };
//     Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
//     Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
//     Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);




  TH1F* hist1[100];
//   TH2F* hist2[100];
  TH2* hist2[100];
  TTree* tree(NULL);

  hist1[1] = new TH1F("xa","x_{a}",500,0,1);
  hist1[2] = new TH1F("xb","x_{b}",500,0,1);
  hist2[1] = new TH2F("yvsxa","y vs x_{a}",500,-1.,1.,500,0.,1.); // 
  hist2[2] = new TH2F("yvsxb","y vs x_{b}",1500,-1.,1.,1500,0.,1.); 
  hist2[3] = new TH2F("xavsxb","x_{a} vs x_{b}",1500,0.,1.,1500,0.,1.);
  hist2[4] = new TH2F("xfvsxa","x_{F} vs x_{a}",1500,-1.,1.,1500,0.,1.);
  hist2[5] = new TH2F("xavszc","x_{a} vs z_{c}",1500,0.,1.,1500,0.,1.);
  hist2[6] = new TH2F("xbvszc","x_{b} vs z_{c}",1500,0.,1.,1500,0.,1.);

double  xa_tree = 0.;
double  xb_tree  = 0.;
double  xf_tree  = 0.;
double  zc_tree  = 0.;
double  pt_tree  = 0.;
double  jt_tree = 0.;
double  y_tree = 0.;


  tree = (TTree*)file->Get("ntuple");
  if(!tree)
    {
      printf("READER.CPP: Not able to open ROOT file!!!!\n");
      abort();
    } 
  tree->SetBranchAddress("xa",&xa_tree);
  tree->SetBranchAddress("xb",&xb_tree);
  tree->SetBranchAddress("xf",&xf_tree);
  tree->SetBranchAddress("zc",&zc_tree);
  tree->SetBranchAddress("jt",&jt_tree);
  tree->SetBranchAddress("pt",&pt_tree);
  tree->SetBranchAddress("y",&y_tree);

  int N=(int)tree->GetEntries();// number of entries en the tree (ntuple)

  printf("N = %d\n",N);

  for(int i=0;i<N;i++)
    {
      tree->GetEntry(i);// get entry number "i"

      hist1[1]->Fill(xa_tree);
      hist1[2]->Fill(xb_tree);
      hist2[1]->Fill(y_tree,xa_tree); //FOR PP
      hist2[2]->Fill(y_tree,xb_tree);
      hist2[3]->Fill(xa_tree,xb_tree);
      hist2[4]->Fill(xf_tree,xa_tree);
      hist2[5]->Fill(xa_tree,zc_tree);
      hist2[6]->Fill(xb_tree,zc_tree);

    }



  TCanvas* canvas1 = new TCanvas("canvas1","histograms",500,500);
  canvas1->Divide(2,4);
  canvas1->Draw();

  canvas1->cd(1);
  hist1[1]->GetXaxis()->SetTitle("x_{a}");
  canvas1->GetPad(1)->SetLogy();
  hist1[1]->GetXaxis()->SetLabelOffset(0.04);
  hist1[1]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(1)->SetLeftMargin(0.2);
  canvas1->GetPad(1)->SetBottomMargin(0.2);
  canvas1->Update();
  hist1[1]->Draw();
 
  canvas1->cd(2);
  hist1[2]->GetXaxis()->SetTitle("x_{b}");
  hist1[2]->GetXaxis()->SetLabelOffset(0.04);
  hist1[2]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(2)->SetLeftMargin(0.2);
  canvas1->GetPad(2)->SetBottomMargin(0.2);
  hist1[2]->Draw();
 
 
  canvas1->cd(3); //PP!!!
  hist2[1]->GetXaxis()->SetTitle("y"); //PP!!!
  hist2[1]->GetYaxis()->SetTitle("x_{a}");
  hist2[1]->GetXaxis()->SetLabelOffset(0.04);
  hist2[1]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(3)->SetLeftMargin(0.2);
  canvas1->GetPad(3)->SetBottomMargin(0.2);
  canvas1->Update();
  //hist2[1]->SetContour(nb);
  hist2[1]->Draw("colz");

  canvas1->cd(4);
  hist2[2]->GetXaxis()->SetTitle("y");
  hist2[2]->GetYaxis()->SetTitle("x_{b}");
  hist2[2]->GetXaxis()->SetLabelOffset(0.04);
  hist2[2]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(4)->SetLeftMargin(0.2);
  canvas1->GetPad(4)->SetBottomMargin(0.2);
  hist2[2]->Draw("colz");

  canvas1->cd(5);
  hist2[3]->GetXaxis()->SetTitle("x_{a}");
  hist2[3]->GetYaxis()->SetTitle("x_{b}");
  hist2[3]->GetXaxis()->SetLabelOffset(0.04);
  hist2[3]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(5)->SetLeftMargin(0.2);
  canvas1->GetPad(5)->SetBottomMargin(0.2);
  hist2[3]->Draw("colz");

  canvas1->cd(6);
  hist2[4]->GetXaxis()->SetTitle("x_{F}");
  hist2[4]->GetYaxis()->SetTitle("x_{a}");
  hist2[4]->GetXaxis()->SetLabelOffset(0.04);
  hist2[4]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(6)->SetLeftMargin(0.2);
  canvas1->GetPad(6)->SetBottomMargin(0.2);
  hist2[4]->Draw("colz");


  canvas1->cd(7);
  hist2[5]->GetXaxis()->SetTitle("x_{a}");
  hist2[5]->GetYaxis()->SetTitle("z_{c}");
  hist2[5]->GetXaxis()->SetLabelOffset(0.04);
  hist2[5]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(7)->SetLeftMargin(0.2);
  canvas1->GetPad(7)->SetBottomMargin(0.2);
  hist2[5]->Draw("colz");

  canvas1->cd(8);
  hist2[6]->GetXaxis()->SetTitle("x_{b}");
  hist2[6]->GetYaxis()->SetTitle("z_{c}");
  hist2[6]->GetXaxis()->SetLabelOffset(0.04);
  hist2[6]->GetYaxis()->SetLabelOffset(0.04);
  canvas1->GetPad(8)->SetLeftMargin(0.2);
  canvas1->GetPad(8)->SetBottomMargin(0.2);
  hist2[6]->Draw("colz");

 
  canvas1->Update();
  //hist2[1]->SetContour(nb);


  return 0;
}
