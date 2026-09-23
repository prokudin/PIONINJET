#ifndef _DRAW_RESULT_SIDIS__
#define _DRAW_RESULT_SIDIS__

#include <sidis.h>
#include <cuba.h>
#include <experimental_data.h>

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

const int step_log = 0;
const int step_lin = 1;
// linear or logariphmic step....

void draw_result_sidis(int npoints, const int ncomputations, const int ndim, int step_var, bool Coridor, integrand_t integrand, TGraph* graphs[]);

#endif // #ifndef _DRAW_RESULT_SIDIS__
