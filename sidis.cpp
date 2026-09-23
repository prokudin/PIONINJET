// Functions for SIDIS.
// Alexei Prokudin 31/01/2009

//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sidis.h>
#include <memory.h>
#include <sstream>

#include "TLatex.h"
#include "TROOT.h"
#include "TColor.h"
#include "TColor.h"
#include "TAttImage.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLine.h"
#include "TAxis.h"
#include "TStyle.h"
// this is fo my palette...
#include "TAttImage.h"
#include "TASPaletteEditor.h"
#include "THStack.h"




using namespace std;



extern  int transversity_order;
extern  int polarised_order;

int polarised_order = standart_lo; // 'STANDARD' SCENARIO, LEADING ORDER

int transversity_order = transversity_lo; // LO

char name[160] = "cteq6l1"; // test on UBUNTU
//char name[160] = "CT10"; // CT10
//char name[160] = "CT10nlo"; // CT10nlo


const int kretzer_ff = 1;
const int kkp_ff     = 2;
const int bfgw_ff    = 3;
const int fdss_ff    = 4; // Daniel de Florian
const int hkns_ff    = 6; // Japan team...
const int akk_ff     = 5; // AKK08...
//int fragmentation_functions = hkns_ff; // Japan team...
//int fragmentation_functions = kretzer_ff;
 int fragmentation_functions = fdss_ff; // Daniel de Florian
//int fragmentation_functions = akk_ff; // AKK08

const int  LO_FF = 0;
const int NLO_FF = 1;
int fragmetation_functions_order = LO_FF; // Leading order fragmentation...
//int fragmetation_functions_order = NLO_FF; // Next To Leading order fragmentation...



void hold(bool keepgoing){
   while(keepgoing){
     cout << "Execute (x) to proceed." << endl;
     string command;
     getline(cin,command);
     if(command == "x") keepgoing = false;
   }
}


sidis* sidis::instance = NULL;
//========================================================= Instance...
sidis* sidis::Instance()
{
  if(!instance)
    instance = new sidis;
  return instance;
}

//description
string sidis::process(void){

  string process        = "test";
  string starget        = "";
  string starget_charge = "";
  string starget_pol    = "";
  string sbeam          = "";
  string sbeam_charge   = "";
  string sbeam_pol      = "";
  string shadron        = "";
  string shadron_charge = "";
  string shadron_pol    = "";




  switch( Beam.GetHadronType() ){
  case PROTON:
    sbeam =  "P";
    break;
  case ANTIPROTON:
    sbeam =  "#bar{P}";
    break;
  case PION:
    sbeam =  "#pi";
    if( Beam.GetCharge() == negative )
      sbeam_charge =  "^{-}";
    if( Beam.GetCharge() == positive )
      sbeam_charge =  "^{+}";
    if( Beam.GetCharge() == neutral )
      starget_charge =  "^{0}";
    break;
  case DEUTRON:
    sbeam =  "D";
    break;
  case NEUTRON:
    sbeam =  "N";
    break;
  default:
    break;
  }

  switch( Beam.GetPolarization() ){
  case longitudinal:
    sbeam_pol =    "^{#rightarrow}";
    break;
  case transverse:
    sbeam_pol =     "^{#uparrow}";
    break;
  case unp:
    break;
  default:
    break;
  }


  switch( Target.GetHadronType() ){
  case PROTON:
    starget =  "P";
    break;
  case ANTIPROTON:
    starget =  "#bar{P}";
    break;
  case PION:
    starget =  "#pi";
    if( Target.GetCharge() == negative )
      starget_charge =  "^{-}";
    if( Target.GetCharge() == positive )
      starget_charge =  "^{+}";
    if( Target.GetCharge() == neutral )
      starget_charge =  "^{0}";
    break;
  case DEUTRON:
    starget =  "D";
    break;
  case NEUTRON:
    starget =  "N";
    break;
  default:
    break;
  }


  switch( Target.GetPolarization() ){
  case longitudinal:
    starget_pol =    "^{#rightarrow}";
    break;
  case transverse:
    starget_pol =     "^{#uparrow}";
    break;
  case unp:
    break;
  default:
    break;
  }


  switch( Produced.GetHadronType() ){
  case PROTON:
    shadron =  "P";
    break;
  case ANTIPROTON:
    shadron =  "#bar{P}";
    break;
  case HADRONS:
    shadron =  "h";
    if( Produced.GetCharge() == negative )
      shadron_charge =  "^{-}";
    if( Produced.GetCharge() == positive )
      shadron_charge =  "^{+}";
    if( Produced.GetCharge() == neutral )
      shadron_charge =  "^{0}";
    break;
  case PION:
    shadron =  "#pi";
    if( Produced.GetCharge() == negative )
      shadron_charge =  "^{-}";
    if( Produced.GetCharge() == positive )
      shadron_charge =  "^{+}";
    if( Produced.GetCharge() == neutral )
      shadron_charge =  "^{0}";
    break;
  case KAON:
    shadron =  "K";
    if( Produced.GetCharge() == negative )
      shadron_charge =  "^{-}";
    if( Produced.GetCharge() == positive )
      shadron_charge =  "^{+}";
    if( Produced.GetCharge() == neutral )
      shadron_charge =  "^{0}";
    break;
  case JET:
    shadron =  "Jet";
    break;
  case PHOTON:
    shadron =  "#gamma";
    break;
  case DEUTRON:
    shadron =  "D";
    break;
  case NEUTRON:
    shadron =  "N";
    break;
  default:
    break;
  }


  switch( Produced.GetPolarization() ){
  case longitudinal:
    shadron_pol =    "^{#rightarrow}";
    break;
  case transverse:
    shadron_pol =     "^{#uparrow}";
    break;
  case unp:
    break;
  default:
    break;
  }

  process = sbeam + sbeam_charge + sbeam_pol + starget + starget_charge + starget_pol;

  process += "#rightarrow " + shadron + shadron_charge + shadron_pol + " X";

  return process;
}

//description
string sidis::print_energy(void){
  stringstream en;

  en << sqrt(GetEnergy());

  string process = " #sqrt{s} = " + en.str() + " (GeV) ";

  return process;
}

//hadron description
string sidis::HadronDescription( HADRON h ){
  string shadron        = "";
  string shadron_charge = "";

  switch( h.GetHadronType() ){
  case PROTON:
    shadron =  "P";
    break;
  case ANTIPROTON:
    shadron =  "#bar{P}";
    break;
  case PION:
    shadron =  "#pi";
    if( h.GetCharge() == negative )
      shadron_charge =  "^{-}";
    if( h.GetCharge() == positive )
      shadron_charge =  "^{+}";
    if( h.GetCharge() == neutral )
      shadron_charge =  "^{0}";
    break;
  case KAON:
    shadron =  "K";
    if( h.GetCharge() == negative )
      shadron_charge =  "^{-}";
    if( h.GetCharge() == positive )
      shadron_charge =  "^{+}";
    if( h.GetCharge() == neutral )
      shadron_charge =  "^{0}";
    break;
  case DEUTRON:
    shadron =  "D";
    break;
  case NEUTRON:
    shadron =  "N";
    break;
  default:
    break;
  }

  return shadron = shadron + shadron_charge;
}



//description
string sidis::limit_x(void){
  stringstream min1, max1;

  min1 << Limits.x_min;
  max1 << Limits.x_max;

  string process = min1.str() + " < x < " + max1.str();

  return process;
}

//description
string sidis::limit_Q2(void){
  stringstream min1, max1;

  min1 << Limits.Q2_min;
  max1 << Limits.Q2_max;

  string process = min1.str() + " < Q^{2} < " + max1.str()+ " (GeV^{2}) ";

  return process;
}

//description
string sidis::limit_W2(void){
  stringstream min1, max1;

  min1 << Limits.W2_min;
  max1 << Limits.W2_max;

  string process = min1.str() + " < W^{2} < " + max1.str()+ " (GeV^{2}) ";

  return process;
}

//description
string sidis::limit_Eh(void){
  stringstream min1, max1;

  min1 << Limits.Eh_min;
  max1 << Limits.Eh_max;

  string process = min1.str() + " < E_{h} < " + max1.str()+ " (GeV) ";

  return process;
}

//description
string sidis::limit_z(void){
  stringstream min1, max1;

  min1 << Limits.z_min;
  max1 << Limits.z_max;

  string process = min1.str() + " < z < " + max1.str();

  return process;
}

//description
string sidis::limit_xF(void){
  stringstream min1, max1;

  min1 << Limits.xF_min;
  max1 << Limits.xF_max;

  string process = min1.str() + " < x_{F} < " + max1.str();

  return process;
}

//description
string sidis::limit_Y(void){
  stringstream min, max;

  min << Limits.y_min;
  max << Limits.y_max;

  string process = min.str() + " < y < " + max.str();

  return process;
}

//description
string sidis::limit_PT(void){
  stringstream min, max;

  min << Limits.pt_min;
  max << Limits.pt_max;

  string process = min.str() + " < P_{T} < " + max.str() + " (GeV) ";

  return process;
}


//description
string sidis::limit_kt(void){
  stringstream min, max;

  min << Limits.kt_min;
  max << Limits.kt_max;

  string process = min.str() + " < k_{#perp} < " + max.str() + " (GeV) ";

  return process;
}

// description as Pave WORK ON IT !!!!!!!!!!!!
TPaveText* sidis::description(){
  //TPaveText* ptstats = new TPaveText(0.78,0.755,0.98,0.995,"brNDC");
  TPaveText* ptstats = new TPaveText(0.2,0.65,0.4,0.85,"brNDC");

  ptstats->SetName("description");
  ptstats->SetBorderSize(0);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);

//   TText *text = ptstats->AddText(name);
  TText *text = ptstats->AddText(experiment_name.c_str());
  text = ptstats->AddText(process().c_str()); // process
  text =  ptstats->AddText(print_energy().c_str()); //energy
  if( Limits.y_cut_used )  text = ptstats->AddText(limit_Y().c_str());
  if( Limits.x_cut_used )  text = ptstats->AddText(limit_x().c_str());
  if( Limits.z_cut_used )  text = ptstats->AddText(limit_z().c_str());
  if( Limits.Q2_cut_used ) text = ptstats->AddText(limit_Q2().c_str());
  if( Limits.W2_cut_used ) text = ptstats->AddText(limit_W2().c_str());
  if( Limits.Eh_cut_used ) text = ptstats->AddText(limit_Eh().c_str());
  if( Limits.xF_cut_used ) text = ptstats->AddText(limit_xF().c_str());
  if( Limits.pt_cut_used ) text = ptstats->AddText(limit_PT().c_str());
  if( Limits.kt_cut_used ) text = ptstats->AddText(limit_kt().c_str());

  return ptstats;
}

//========================================================= DY intitialization...
sidis::sidis()
{
  const int SUBSET = 0;




 if( strcmp( name, "GRV98lo.LHgrid") == 0 ){
  LHAPDF::initPDFSet("MSTW2008nlo68cl.LHgrid", SUBSET); // init mstw instead!
 }
 else {
  LHAPDF::initPDFSet(name, SUBSET); //
 };


  double x = 0.2;
  double Q2 = pow2(20.);
  double UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU;
  double Q = sqrt(Q2);

//  vector<double> dens = pdf.xfx(x, Q); // NOTE LINUX
  vector<double> dens = LHAPDF::xfx(x, Q); // NOTE MAC ONLY

  UPV =  dens[8] / x; // LHAPDF return already u = uvalence + u_bar
  DNV =  dens[7] / x; // LHAPDF return already d = dvalence + d_bar
  USEA = dens[4] / x; //
  DSEA = dens[5] / x; //
  STR =  dens[9] / x; //
  CHM =  dens[10]/ x; //
  BOT =  dens[11]/ x; //
  TOP =  dens[12]/ x; //
  GLU =  dens[6] / x; //



  printf("=====TEST STRUCTM\n  %f %f %f %f %f %f %f %f %f\n",
	 UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU);
  cout << "=====End TEST STRUCTM" << endl;



  double U, D, UB, DB, ST, GL, G1P, G1N;
  int iset = polarised_order;
  PARPOL(iset, x, Q2, U, D, UB, DB, ST, GL, G1P, G1N);
  printf("=====TEST PARPOL\n  %f %f %f %f %f %f %f %f\n",
	 U, D, UB, DB, ST, GL, G1P, G1N);
  cout << "=====End TEST PARPOL" << endl;


  int ISET_DSSV = 0; //BEST fit for DSSV
  DSSVINI(ISET_DSSV);
  double DUV,DDV,DUBAR,DDBAR,DSTR,DGLU;
  DSSVFIT(x,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU);
  printf("=====TEST DSSVFIT\n  %f %f %f %f %f %f \n",
	 DUV,DDV,DUBAR,DDBAR,DSTR,DGLU);
  cout << "=====End TEST DSSVFIT" << endl;




//hold(true);
  // Cuts are not used by default...
  Limits.xF_cut_used  = false;
  Limits.x_cut_used  = false;
  Limits.y_cut_used   = false;
  Limits.z_cut_used   = false;
  Limits.Q2_cut_used   = false;
  Limits.W2_cut_used   = false;
  Limits.pt_cut_used  = false;
  Limits.kt_cut_used  = false;

  // x, and xF, kt have default values:
  Limits.x_min =  0.;
  Limits.x_max =  1.;
  Limits.xF_min = -1.;
  Limits.xF_max =  1.;
  Limits.kt_min =  -3.*sqrt(Params.GetKt2Average());
  Limits.kt_max =  3.*sqrt(Params.GetKt2Average());





}


 //Settings for hoppet:
 void set_hoppet( void ){ // this function will set hoppet correctly

 int    dummy = 0;


//      ! start the dglap evolution/convolution package
 double   ymax  = 12.;        //      ! max value of ln 1/x
 double   dy    = 0.1;        //     ! the internal grid spacing (smaller->higher accuarcy)
                                // ! 0.1 should provide at least 10^{-3} accuracy
//   double   Qmin  = sqrt(0.3);  //   ! smallest Q value in tabulation
//   double   Qmin  = 3.2;  //   ! smallest Q value in tabulation
//   double   Qmin  = C1/bmax;  //   ! smallest Q value in tabulation
 double   Qmin  = 0.5;  //   ! smallest Q value in tabulation
 double   Qmax  = 1.e+5;       // ! largest Q value in tabulation
 double   dlnlnQ = dy/4.;     // ! tabulation spacing in dlnlnQ (dy/4 recommended)
 int   nloop  = 1;            // ! the number of loops to initialise (max=3!)
 int   order  = -6;           // ! numerical interpolation order (-6 is a good choice)

//   ! change the ipdf and scheme to select different pdfs
//   int  ipdf = 1;               // ! 1=f1, 2=f1Tp1, 3=g1, 4=g1T1, 5=h1,
                                // ! 6=h1p1, 7=h1Lp1, 8=h1Tp1, 9=H1p1
 int   scheme = 5;            // ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar,
                                // ! 4=frag, 5=TransMsbar, 6=SiversMsbar
//     int   scheme = 1;            // ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar,
                                // ! 4=frag, 5=TransMsbar, 6=SiversMsbar

// NOTE NO DUMMY VARIABLE NEEDED IN HOPPET-1-15!!!! CHECK THIS
//      ! call this once at the beginning of your program
//   hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,
//                              scheme,dummy);
// NOTE: new hoppet
   hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,
                              scheme);


   cout <<  "Splitting functions initialised!" << endl;

//      ! tell hoppet to use a variable flavour number scheme with
//      ! the following c,b,t quark masses
//      hoppetSetVFN(1.4,4.5,175);

//      ! set the initial scale and the coupling there (in general the PDF
//      ! and coupling may be specified at different scales -- this is not
//      ! done here) and
//      double asQ0 =  0.706862, Q0=sqrt(0.3);
//      double asQ0 =  0.13939, Q0=91.187;
//      double asQ0 = 0.127, Q0=91.1876;

      double Q0 = sqrt(2.4), asQ0 = AlphaS(Q0);
      //cout << "alphas0 " << asQ0 << endl;
      //cout << "alphaMz" << AlphaS(91.1876) << endl;
      //hold(true);

//      ! the ratio xmu = mu_F/mu_R to be used in the evolution.
      double xmu  = 1.0;

//      ! carry out the evolution to create a tabulation, corresponding
//      ! to the initial condition pdf_init(...) given below
//      hoppetEvolve(asQ0,Q0,nloop,xmu,pdf_init,Qmin,dummy);
//NOTE new hoppet for TRANSVERSITY and COLLINS:
      hoppetEvolve(asQ0,Q0,nloop,xmu,pdf_init,Q0); // FENG 2014 COLLINS



//      ! alternatively if you need to repeat the same evolution on very
//      ! many pdf sets, used a cached evolution (set up once, use many times
//      ! and gain a factor 3-4 in speed.
//      !call hoppetPreEvolve(asQ0,Q0,nloop,xmu,Q0)
//      !call hoppetCachedEvolve(evolvePDF)
      cout<< "Evolution done!" << endl;
}

 //Settings for hoppet:->reset evolution when parameters are changed
 void reset_hoppet( void ){ // this function will set hoppet correctly

      int   nloop  = 1;            // ! the number of loops to initialise (max=3!)

      double Q0 = sqrt(2.4), asQ0 = AlphaS(Q0);
//      ! the ratio xmu = mu_F/mu_R to be used in the evolution.
      double xmu  = 1.0;

//NOTE new hoppet for TRANSVERSITY and COLLINS:
      hoppetEvolve(asQ0,Q0,nloop,xmu,pdf_init,Q0); // FENG 2014 COLLINS
      cout<< "Evolution reset done!" << endl;
}


//----------------------------------------------------------------------
// the initial condition for COLLINS AND TRANSVERSITY!
void  pdf_init(const double & x, const double & Q,  double * pdf) {


  // MAIN WITH FENG 2014
  // we will implement the following: 1/2 (f1(x,Q0) + g1(x,Q0))
  sidis& SIDIS = *sidis::Instance();
  double U,D,USEA,DSEA,STR,CHM,BOT,TOP,GLU;



  vector<double> dens = LHAPDF::xfx(x, Q); //

  U    =  dens[8]  ; // LHAPDF return already x * u = uvalence + u_bar
  D    =  dens[7]  ; // LHAPDF return already x * d = dvalence + d_bar
  USEA = dens[4]  ; //
  DSEA = dens[5]  ; //
  STR =  dens[9]  ; //
  CHM =  dens[10] ; //
  BOT =  dens[11] ; //
  TOP =  dens[12] ; //
  GLU =  dens[6]  ; //


  double x1;
  if (x <= 1.01e-5 ) // DSSV does not work below 1.e-5 --> freeze!
  {x1 = 1.01e-5;}
  else
  {x1 = x;};

  double DUV,DDV,DUBAR,DDBAR,DSTR,DGLU;
  DSSVFIT(x1,Q*Q,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU);
  /* C           ALWAYS X*DISTRIBUTION IS RETURNED !!!                            C */
  /* C       OUTPUT:                                                              C */
/* C           DUV   :  X * U VALENCE DISTRIBUTION                              C */
/* C           DDV   :  X * D VALENCE DISTRIBUTION                              C */
/* C           DUBAR :  X * UBAR DISTRIBUTION                                   C */
/* C           DDBAR :  X * DBAR DISTRIBUTION                                   C */
/* C           DSTR  :  X * STRANGE DISTRIBUTION                                C */
/* C           DGLU  :  X * GLUON DISTRIBUTION                                  C */


// MAIN CHOICE FOR 2015
  pdf[-3+6] = 0.; //strange
  pdf[ 3+6] = 0.;//anti_strange
  pdf[ 2+6] = 0.5 * SIDIS.transversity(x).up * (U + DUV + DUBAR) ;  // x*UP quark times n_up x^a (1-x)^b etc
  pdf[-2+6] = 0.5 * SIDIS.transversity(x).down * (D + DDV + DDBAR); //x*DOWN quark
// MAIN CHOICE FOR 2015

// FRAGMENTATION IS HERE:


   // ffset   1,2,3 means K, KKP, BFGW
   // fforder 0,1 is LO, NLO(MSbar)
   // ihadron 1,2,3,4,5 is pi,K,h,p,n
   // icp     1,2,3 chooses between particle, anti-particle or sum of both
   // ipi     1,2,3 is a flag for BFGW; inactive for K, KKP
   // icharge 0,1,2,3 is 0,+,-,+&-
   // int ffset  = 1; // Kretzer
  int ffset  = 4; // DSS
  int fforder= 1; // NLO
  int ihadron= 1; // PION
  int icp    = 1; //fragmentation.icp;
  int ipi    = 1;
  int icharge= 1; // positive

  double dff[11];

  dff[0] = 0.;
  dff[1] = 0.;
  dff[2] = 0.;
  dff[3] = 0.;
  dff[4] = 0.;
  dff[5] = 0.;
  dff[6] = 0.;
  dff[7] = 0.;
  dff[8] = 0.;
  dff[9] = 0.;
  dff[10]= 0.;

// * Subroutine "dlib" returns D(z,Q2) [NOT z*D(z,Q2)] as an array        *
  DLIB(x,Q*Q,dff,ffset,fforder,ihadron,icharge,icp,ipi);



// this part is for Collins FF:
  pdf[ 6+6] = 0.;
  pdf[-6+6] = 0.;
  pdf[ 5+6] = 0.;
  pdf[-5+6] = 0.;
  pdf[ 4+6] = 0.;
  pdf[-4+6] = 0.;


  pdf[ 1+6] = SIDIS.Params.Collins.parameters.n_up * pow(x,SIDIS.Params.Collins.parameters.a_up) *
  pow(1.-x, SIDIS.Params.Collins.parameters.b_up)  * x * dff[1+5]; // Collins H(3) up == FAV(z) NB collins = 2 N x^a (1-x)^b  //U->pi+

  pdf[-1+6] = SIDIS.Params.Collins.parameters.n_down * pow(x, SIDIS.Params.Collins.parameters.a_down) *
  pow(1.-x, SIDIS.Params.Collins.parameters.b_down) * x * dff[-1+5]; // Collins H(3) down == UNFAV(z) UBAR --> pi+
  pdf[3+6] = SIDIS.Params.Collins.parameters.n_down * pow(x, SIDIS.Params.Collins.parameters.a_down) *
  pow(1.-x, SIDIS.Params.Collins.parameters.b_down) * x * dff[3+5]; // Collins H(3) strange == UNFAV1(z) strange --> pi+
//NOTE HOPPET EVOLVES X * DISTRIBUTION(X) !!!!!
  pdf[-3+6] = 0.;
//ATTENTION NOTE NOTE:
//ATTENTION!!!!!!
// HOPPEN DOES NOT EVOLVE pdf[ 6+6], I DO NOT KNOW WHY!!!!!!
//ATTENTION
}






//========================================================= Gaussian
double sidis::gauss( double kt, double width ){
  double result = 0.;

  if(width != 0.){
  result = 1./( PI*width ) * exp( -kt*kt / width );
  } else
  if(width == 0.){
  result = 1.; // in order to implement transversity  kt dependence!
  };

  return result;
}



//========================================================= Unpolarized partcontent
void sidis::Unpolarised(HADRON& h, double x, double Q2)
  // Returns the parton content at x,Q2
{

  switch( h.GetHadronType() ){
  case PROTON:
    h.partcontent = unpolarised(x, Q2);
    break;
  case ANTIPROTON:
    h.partcontent = unpolarised(x, Q2);
    Antiproton(h.partcontent);
    break;
  case PION:
    h.partcontent = unpolarised_pion(x, Q2);
    if( h.GetCharge() == negative )
      PionNegative(h.partcontent);
    break;
  case DEUTRON:
    h.partcontent = unpolarised(x, Q2);
    Deutron(h.partcontent);
    break;
  case NEUTRON:
    h.partcontent = unpolarised(x, Q2);
    Neutron(h.partcontent);
    break;
  default:
    break;
  }

};



//========================================================= Unpolarized partcontent
PARTONCONTENT sidis::unpolarised(double x, double Q2)
  // Returns the parton content at x,Q2
{
  PARTONCONTENT partcontent;

  double Q = sqrt(Q2);

  vector<double> dens;


  dens = LHAPDF::xfx(x, Q); //NOTE MAC ONLY

//   // 0, 1, 2, 3, 4, 5 = tbar, barb, barchm, bars, ubar, dbar;
//   // 6 = g;
//   // 7, 8, 9, 10, 11, 12 = d, u, s, chm, b, t

  //Assign parton densities PROTON:
  if(!VALENCE_SIDIS){ // WE USE STANDARD PARAMETRIZATION!!!!
  partcontent.up     =  dens[8] / x; // LHAPDF return already u = uvalence + u_bar
  partcontent.down   =  dens[7] / x; // LHAPDF return already d = dvalence + d_bar
  partcontent.anti_up  =  dens[4] / x; //
  partcontent.anti_down=  dens[5] / x; //
  partcontent.strange    =  dens[9] / x; //
  partcontent.anti_strange =  dens[3] / x; //
  partcontent.charm    =  dens[10]/ x; //
  partcontent.anti_charm =  dens[2] / x; //
  partcontent.bottom    =  dens[11]/ x; //
  partcontent.anti_bottom =  dens[1] / x; //
  partcontent.top    =  dens[12]/ x; //
  partcontent.anti_top =  dens[0] / x; //
  partcontent.glu    =  dens[6] / x; //
  }
  if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
  partcontent.up     =  (dens[8]-dens[4]) / x; // LHAPDF return already u = uvalence + u_bar
  partcontent.down   =  (dens[7]-dens[5]) / x; // LHAPDF return already d = dvalence + d_bar
  partcontent.anti_up  =  dens[4] / x; //
  partcontent.anti_down=  dens[5] / x; //
  partcontent.strange    =  dens[9] / x; //
  partcontent.anti_strange =  dens[3] / x; //
  partcontent.charm    =  dens[10]/ x; //
  partcontent.anti_charm =  dens[2] / x; //
  partcontent.bottom    =  dens[11]/ x; //
  partcontent.anti_bottom =  dens[1] / x; //
  partcontent.top    =  dens[12]/ x; //
  partcontent.anti_top =  dens[0] / x; //
  partcontent.glu    =  dens[6] / x; //
  }

  return partcontent;

}


//========================================================= Unpolarized partcontent PION
PARTONCONTENT sidis::unpolarised_pion(double x, double Q2)
  // Returns the parton content at x,Q2
{
  PARTONCONTENT partcontent;


  double ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL;
  double Q = sqrt(Q2);
   //
  int MODE = 2;
  //GRVPILO(x,Q,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZTB,ZGL); // RETURNS DISTRIBUTIONS!!!
  PION(x,Q,MODE,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL); // RETURS x * Distributions(x)
  //grvpilo(x,Q,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL);
  //grvpilo(x,Q,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL);
   //Assign parton densities PION:
  if(!VALENCE_SIDIS){ // WE USE STANDARD PARAMETRIZATION!!!!
  partcontent.up     = ( ZUV + ZUDB ) / x; //
  partcontent.down   = ZUDB / x; //
  partcontent.anti_up  = ZUDB / x; //
  partcontent.anti_down= ( ZUV + ZUDB ) / x; //
  partcontent.strange    = ZSB / x; //
  partcontent.anti_strange = ZSB / x; //
  partcontent.charm    = ZCB / x; //
  partcontent.anti_charm = ZCB / x; //
  partcontent.bottom    = ZBB / x; //
  partcontent.anti_bottom = ZBB / x; //
  partcontent.top    = 0 / x; //
  partcontent.anti_top = 0 / x; //
  partcontent.glu    = ZGL / x; //
  }
  if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
  partcontent.up     = ZUV / x; //
  partcontent.down   = ZUDB / x; //
  partcontent.anti_up  = ZUDB / x; //
  partcontent.anti_down= ZUV / x; //
  partcontent.strange    = ZSB / x; //
  partcontent.anti_strange = ZSB / x; //
  partcontent.charm    = ZCB / x; //
  partcontent.anti_charm = ZCB / x; //
  partcontent.bottom    = ZBB / x; //
  partcontent.anti_bottom = ZBB / x; //
  partcontent.top    = 0 / x; //
  partcontent.anti_top = 0 / x; //
  partcontent.glu    = ZGL / x; //
  }

  return partcontent;
}


 // definitions of DEUTRON
void sidis::Deutron( struct PARTONCONTENT& partcontent ){

    partcontent.up      += partcontent.down;
    partcontent.down     = partcontent.up;
    partcontent.anti_up   += partcontent.anti_down;
    partcontent.anti_down  = partcontent.anti_up;
    partcontent.strange     *= 2.;
    partcontent.anti_strange  *= 2.;
    partcontent.charm     *= 2.;
    partcontent.anti_charm  *= 2.;
    partcontent.glu     *= 2.;
    partcontent.bottom     *= 2.;
    partcontent.anti_bottom  *= 2.;
    partcontent.top     *= 2.;
    partcontent.anti_top  *= 2.;
    // IF WE USE VALENCE PARAMETRIZATION --> FORMULAS ARE THE SAME!

};

// definitions of NEUTRON
void sidis::Neutron( struct PARTONCONTENT& partcontent ){
      double neutron_pdf       = partcontent.up;
      partcontent.up    = partcontent.down;
      partcontent.down  = neutron_pdf;
      neutron_pdf       = partcontent.anti_up;
      partcontent.anti_up = partcontent.anti_down;
      partcontent.anti_down  = neutron_pdf;
    // IF WE USE VALENCE PARAMETRIZATION --> FORMULAS ARE THE SAME!
 };

// definitions of ANTIPROTON
void sidis::Antiproton( struct PARTONCONTENT& partcontent ){
      double antiproton_pdf       = partcontent.up;
      partcontent.up    = partcontent.anti_up;
      partcontent.anti_up = antiproton_pdf;

      antiproton_pdf    = partcontent.down;
      partcontent.down  = partcontent.anti_down;
      partcontent.anti_down= antiproton_pdf;

      antiproton_pdf    = partcontent.strange;
      partcontent.strange   = partcontent.anti_strange;
      partcontent.anti_strange  = antiproton_pdf;
};


// definitions of NEGATIVE PION
void sidis::PionNegative( struct PARTONCONTENT& partcontent ){
  double tmp =  partcontent.up;
  double tmp1 = partcontent.anti_down;

  partcontent.up   = partcontent.anti_up;
  partcontent.anti_down  = partcontent.down;

  partcontent.anti_up = tmp;
  partcontent.down   = tmp1;

};


const int grv_g1 = 1; // GLUK REYA
const int dssv_g1 = 2; // DE FLORIAN ETC
//int helicity_model = grv_g1;
int helicity_model = dssv_g1;


//========================================================= Helicity partcontent
PARTONCONTENT sidis::helicity(double x, double Q2)
  // Returns the helicity content at x,Q2
{
  PARTONCONTENT partcontent;

//  double Q = sqrt(Q2);

  vector<double> dens;



  for (int i=0; i<13; i++) dens.push_back(0.);
// *********************************************************************
// *                                                                   *
// *    POLARIZED RADIATIVELY GENERATED LO AND NLO PARTON DENSITIES    *
// *                                                                   *
// *         M. GLUCK, E. REYA, M. STRATMANN AND W. VOGELSANG,         *
// *              Phys.Rev.D63:094005,2001,    hep-ph/0011215          *
// *                                                                   *
// *               PROBLEMS/QUESTIONS TO wvogelsang@bnl.gov            *
// *            OR TO marco.stratmann@physik.uni-regensburg.de         *
// *                                                                   *
// *   INPUT:   ISET = number of the parton set :                      *
// *            ISET = 1  'STANDARD' SCENARIO, NEXT-TO-LEADING ORDER   *
// *                      (MS-bar)                                     *
// *                      (DATA FILE 'std2000_nlo.grid' UNIT=11, TO BE *
// *                       DEFINED BY THE USER )                       *
// *            ISET = 2  'VALENCE' SCENARIO,  NEXT-TO-LEADING ORDER   *
// *                      (MS-bar)                                     *
// *                      (DATA FILE 'val2000_nlo.grid' UNIT=22, TO BE *
// *                       DEFINED BY THE USER )                       *
// *            ISET = 3  'STANDARD' SCENARIO, LEADING ORDER           *
// *                      (DATA FILE 'std2000_lo.grid' UNIT=33, TO BE  *
// *                       DEFINED BY THE USER )                       *
// *            ISET = 4  'VALENCE' SCENARIO,  LEADING ORDER           *
// *                      (DATA FILE 'val2000_lo.grid' UNIT=44, TO BE  *
// *                       DEFINED BY THE USER )                       *
// *                                                                   *
// *            X  = Bjorken-x       (between  1.E-4  and  1)          *
// *            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        *
// *                                                                   *
// *   OUTPUT:  U = x * DELTA u                                        *
// *            D = x * DELTA d                                        *
// *            UB = x * DELTA ubar                                    *
// *            DB = x * DELTA dbar                                    *
// *            ST = x * DELTA STRANGE                                 *
// *            GL = x * DELTA GLUON                                   *
// *            G1P = g_1^proton                                       *
// *            G1N = g_1^neutron                                      *
// *                                                                   *
// *          (  For the parton distributions always x times           *
// *                   the distribution is returned .                  *
// *                 This is NOT the case for g1(p,n)  )               *
// *                                                                   *
// *            The sets are the result of a combined fit to           *
// *            data for the spin asymmetries A_1 (p,n,d)              *
// *                                                                   *
// *            Note: No charm is included                             *
  if( helicity_model == grv_g1 ){
  int ISET = polarised_order; //LO
  double U, D, UB, DB, ST, GL, G1P, G1N;

  PARPOL(ISET, x, Q2, U, D, UB, DB, ST, GL, G1P, G1N);


  dens[8] = U;
  dens[7] = D;
  dens[4] = UB;
  dens[5] = DB;
  dens[9] = ST;
  dens[3] = ST;
  dens[10]= 0.;
  dens[2] = 0.;
  dens[11]= 0.;
  dens[1] = 0.;
  dens[12]= 0.;
  dens[0] = 0.;
  dens[6] = GL;




  } else if ( helicity_model == dssv_g1 ){
  double DUV,DDV,DUBAR,DDBAR,DSTR,DGLU;

// C       OUTPUT:                                                              C
// C           DUV   :  X * U VALENCE DISTRIBUTION                              C
// C           DDV   :  X * D VALENCE DISTRIBUTION                              C
// C           DUBAR :  X * UBAR DISTRIBUTION                                   C
// C           DDBAR :  X * DBAR DISTRIBUTION                                   C
// C           DSTR  :  X * STRANGE DISTRIBUTION                                C
// C           DGLU  :  X * GLUON DISTRIBUTION                                  C
// C                                                                            C
// C       IMPORTANT:                                                           C
// C           ALWAYS X*DISTRIBUTION IS RETURNED !!!                            C
// C           ALL PDFs ARE IN THE MSbar SCHEME                                 C

  DSSVFIT(x,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU);
  dens[8] = DUV + DUBAR;
  dens[7] = DDV + DDBAR;
  dens[4] = DUBAR;
  dens[5] = DDBAR;
  dens[9] = DSTR;
  dens[3] = DSTR;
  dens[10]= 0.;
  dens[2] = 0.;
  dens[11]= 0.;
  dens[1] = 0.;
  dens[12]= 0.;
  dens[0] = 0.;
  dens[6] = DGLU;


  }


//   // 0, 1, 2, 3, 4, 5 = tbar, barb, barchm, bars, ubar, dbar;
//   // 6 = g;
//   // 7, 8, 9, 10, 11, 12 = d, u, s, chm, b, t

  //Assign parton densities PROTON:
  partcontent.up     =  dens[8] / x; // LHAPDF return already u = uvalence + u_bar
  partcontent.down   =  dens[7] / x; // LHAPDF return already d = dvalence + d_bar
  partcontent.anti_up  =  dens[4] / x; //
  partcontent.anti_down=  dens[5] / x; //
  partcontent.strange    =  dens[9] / x; //
  partcontent.anti_strange =  dens[3] / x; //
  partcontent.charm    =  dens[10]/ x; //
  partcontent.anti_charm =  dens[2] / x; //
  partcontent.bottom    =  dens[11]/ x; //
  partcontent.anti_bottom =  dens[1] / x; //
  partcontent.top    =  dens[12]/ x; //
  partcontent.anti_top =  dens[0] / x; //
  partcontent.glu    =  dens[6] / x; //

  return partcontent;

}

//========================================================= Helicity Distribution
void sidis::HelicityDistribution( HADRON& h, double x, double Q2)
  // Returns the parton content at x,Q2
{

  h.partcontent = helicity(x, Q2);

  switch( h.GetHadronType() ){
  case PROTON:
    h.partcontent = helicity( x, Q2);
    break;
  case ANTIPROTON:
    h.partcontent = helicity( x, Q2);
    Antiproton(h.partcontent);
    break;
  case PION:
    break;
  case DEUTRON:
    h.partcontent = helicity( x, Q2);
    Deutron(h.partcontent);
    break;
  case NEUTRON:
    h.partcontent = helicity( x, Q2);
    Neutron(h.partcontent);
    break;
  default:
    break;
  }

};



//========================================================= fragmentation functions
void sidis::fragmentation(HADRON& h, int charge, double z, double Q2)
  // Returns fragmentation content at z,Q2
{
  double dff[11];

  dff[0] = 0.;
  dff[1] = 0.;
  dff[2] = 0.;
  dff[3] = 0.;
  dff[4] = 0.;
  dff[5] = 0.;
  dff[6] = 0.;
  dff[7] = 0.;
  dff[8] = 0.;
  dff[9] = 0.;
  dff[10]= 0.;


  if( fragmentation_functions ==  kretzer_ff || fragmentation_functions ==  kkp_ff
      || fragmentation_functions ==  bfgw_ff || fragmentation_functions == fdss_ff
      || fragmentation_functions ==  akk_ff || fragmentation_functions ==  hkns_ff){


// * Library of FF parametrizations updated to be posted on               *
// * http://www.pv.infn.it/~radici/FFdatabase/                            *
// * site maintained by Marco Radici (radici@pv.infn.it)                  *
// * and Rainer Jakob (rainer@theorie.physik.uni-wuppertal.de)            *
// * within the EU network:                                               *
// * "Hadronic Physics with High Energy Electromagnetic Probes"           *
// *                                                                      *
// * Also posted on that site are the individual FF sets included in this *
// * library as distributed by the corresponding authors.                 *
// ************************************************************************
// * K, KKP and BFGW fragmentation functions for                          *
// * iparton -> ihadron                                                   *
// *                                                                      *
// * Subroutine "dlib" returns D(z,Q2) [NOT z*D(z,Q2)] as an array        *
// * "dff(iparton)" of dimension dff(-5:5)                                *
// *                                                                      *
// * iparton = 5,4,3,2,1,0,-1,...,-5 means b,c,s,d,u,g,ubar,...,bbar      *
// *                                                                      *
// * [Q2] = GeV^2                                                         *
// *                                                                      *
// * ffset = 1,2,3 means K, KKP, BFGW                                     *
// ************************************************************************
// *        ALEXEI PROKUDIN prokudin@to.infn.it                           *
// *        making the same interface FF fDSS as for DLIB                 *
// ************************************************************************
// * ffset = 4 means                                                      *
// *        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                     *
// *  D.de Florian, R.Sassot, M.Stratmann   hep-ph/0703242)               *
// *    Phys.Rev.D.75:114010,2007                                         *
// *                                                                      *
// ************************************************************************
// *        ALEXEI PROKUDIN prokudin@to.infn.it                           *
// *        making the same interface FF AKK as for DLIB                  *
// *                    26/03/2008                                        *
// * ffset = 5 means                                                      *
// c----------------------------------------------------------------------
// *
// c     AKK ROUTINES 2008
// *
// c----------------------------------------------------------------------*
// *
// *      SUBROUTINE AKK(IH,Z,Q,DH)*
// *
// * AKK Update: Improvements from New Theoretical Input and Experimental Data.
// * S. Albino, B.A. Kniehl, G. Kramer . Mar 2008.
// * e-Print: arXiv:0803.2768 [hep-ph]
// * Details are to be found in the corresponding references:             *
// * S.Kretzer, Phys.Rev.D62, 054001 (2000)                               *
// * B.A.Kniehl, G.Kramer, B.Potter, Nucl.Phys.B582, 514 (2000)           *
// * L.Bourhis, M.Fontannaz, J.P.Guillet, M.Werlen,Eur.Phys.J.C19,89(2001)*
// *                                                                      *
// * fforder=0,1 is LO, NLO(MSbar)                                        *
// * no LO set for BFGW                                                   *
// *                                                                      *
// * ihadron=1,2,3,4,5 is pi,K,h,p,n                                      *
// * no pi,K sets for BFGW                                                *
// * no p,n sets for K, BFGW                                              *
// * no n set for fDSS                                                    * AP
// *                                                                      *
// * icharge=0,1,2,3 is 0,+,-,+&-                                         *
// * note: 3 = +&- = charge sum (NOT average) in this library             *
// * for the neutral (icharge=0) particles K^0, n:                        *
// * icp = 1,2,3 chooses between particle, anti-particle or sum of both   *
// * icp is inactive for pi,h,p and charged Kaons                         *
// *                                                                      *
// * Following KKP, FFs into {pi^0;(anti-)K^0} and into {pi^+/-;K^+/-}    *
// * are related by isospin.                                              *
// *                                                                      *
// * ipi = 1,2,3 is a flag for BFGW; inactive for K, KKP:                 *
// * ipi = 1: best fit     (formula (8))  of BFGW                         *
// * ipi = 2: large Ng set (formula (9))  of BFGW                         *
// * ipi = 3: low Ng set   (formula (10)) of BFGW                         *
// ************************************************************************
// * If options (ihadron, icharge etc.) are chosen which do not exist for *
// * a given parametrization (ffset) or if a flag is chosen outside its   *
// * range as defined above then the code returns a corresponding warning *
// * and stops.                                                           *
// ************************************************************************
// ************************************************************************
// **************************************************************
// *        ALEXEI PROKUDIN prokudin@to.infn.it    2009         *
// *        making the same interface FF fHKNS as for DLIB      *
// **************************************************************
// ************************************************************************
// *        ALEXEI PROKUDIN prokudin@to.infn.it                           *
// *        making the same interface FF AKK as for DLIB                  *
// *                    26/03/2008                                        *
// * ffset = 6 means                                                      *
// C *********************************************************************
// C  hknsff07.f  Version 1.0                                2006/FEB/27
// C
// C  [Package for the HKNS fragmentaion functions]
// C  Reference:
// C    Determination of fragmentation functions and their uncertainties
// C    M. Hirai, S. Kumano, T.-H. Nagai, and K. Sudoh
// C    hep-ph/0702250 (KEK-TH-1108), submitted for publication.
// C *********************************************************************
// C ---------------------------------------------------------------------
// C  SUBROUTINE HKNSFF(Q2,X,ISET,ICHARGE,FF,GRAD):
// C
// C   Subroutine HKNSFF returns the values of fragmentaion functions
// C   and their gradient terms at specified Q^2 and x point
// C   by interpolating the grid data.
// C   [ Log(Q^2): LINEAR INTERPOLATION, x: CUBIC SPLINE INTERPOLATION ]
// C
// C   INPUT:
// C     Q2, X ... Q^2 and x values at which the functions are calculated.
// C               Available range: 10^-2 <= X <= 1.0,
// C                           1.0 GeV^2 <= Q^2 <= 10^8 GeV^2.
// C     ISET=1: Pion LO Fragmentation functions and their gradient terms
// C          2: Pion NLO
// C          3: Kaon LO
// C          4: Kaon NLO
// C          5: Proton LO
// C          6: Proton NLO
// C
// C     ICHARGE=1: pi^+, K^+, or proton
// C     ICHARGE=2: pi^-, K^-, or neutron
// C     ICHARGE=3: pi^0=[pi^+ + pi^-]/2, [K^0+K^0b]/2, or [p+pb]/2
// C        If you want to obtain the fragmentaion functions for each
// C        K0, K0b, pb, or nb, you may use the relations in Appendix of
// C        hep-ph/0702250.
// C
// C   OUTPUT: Arrays FF(-5:5) & GRAD(I,J)
// C
// C     FF(I) --> HKNS fragmentation functions (FFs).
// C      I = -5 ... b-bar quark (D_b-bar = D_b)
// C          -4 ... c-bar quark (D_c-bar = D_c)
// C          -3 ... s-bar quark
// C          -2 ... d-bar quark
// C          -1 ... u-bar quark
// C           0 ... gluon D_g(x)
// C           1 ... u quark
// C           2 ... d quark
// C           3 ... s quark
// C           4 ... c quark
// C           5 ... b quark
// C
// C     GRAD(I,J) --> Gradient terms of HKSN FFs
// C      I is the same index as the one in FF(I).
// C      J indicates the parameter index for a gradient term dFF(I)/da_J
// C      (a_J = parameter).
// C
// C Pion,J= 1..2: g (2ndM, alpha)           2ndM = second moment
// C         3..5: u (2ndM, alpha, beta)
// C         6..8: d (2ndM, alpha, beta)
// C        9..11: c (2ndM, alpha, beta)
// C       12..14: b (2ndM, alpha, beta)
// C   For example, the above J=14 indicates d D_b^{pi^+}/d beta_b^{pi^+}.
// C
// C Kaon,J= 1..2: g (2ndM, beta)
// C         3..5: u (2ndM, alpha, beta)
// C         6..8: d (2ndM, alpha, beta)
// C        9..11: sb(2ndM, alpha, beta)
// C       12..14: c (2ndM, alpha, beta)
// C       15..17: b (2ndM, alpha, beta)
// C
// CProton,J=1..2: g (2ndM, beta)
// C         3..5: u (2nsM, alpha, beta)
// C         6..8: qb(2ndM, alpha, beta)
// C        9..11: c (2ndM, alpha, beta)
// C       11..13: b (2ndM, alpha, beta)
// C
// C   NOTE: The returned values are not multiplied by x.
// C
// C      *  Error matrix can be used by declaring a common block:
// C         COMMON/ERRM/EM(17,17). This matrix is defined as
// C         the inverse matrix of Hessian multiplied by Delta chi^2:
// C         EM(i,j)=Delta chi^2*H_ij^-1.
// C         The values of Delta chi^2 are as follows:
// C         15.9359730(pion), 19.1977555(kaon), 14.8470228(proton).
// C *********************************************************************



    int positive    = 1;
    int negative    = 2;
    int neutral     = 0;


   // ffset   1,2,3 means K, KKP, BFGW
   // fforder 0,1 is LO, NLO(MSbar)
   // ihadron 1,2,3,4,5 is pi,K,h,p,n
   // icp     1,2,3 chooses between particle, anti-particle or sum of both
   // ipi     1,2,3 is a flag for BFGW; inactive for K, KKP
   // icharge 0,1,2,3 is 0,+,-,+&-
   // int ffset  = 1; // Kretzer
  int ffset  = fragmentation_functions; // Kretzer&ko
  int fforder= fragmetation_functions_order; // LO or NLO
  int ihadron= 0;
  int icp    = 1; //fragmentation.icp;
  int ipi    = 1;
  int icharge= charge;

  switch( h.GetHadronType() ){
  case PION:
    ihadron = 1;
    break;
  case KAON:
    ihadron = 2;
    break;
  case HADRONS:
    ihadron = 3;
    break;
  case PROTON:
    ihadron = 4;
    break;
  case NEUTRON:
    ihadron = 5;
    break;
  default:
    break;
  }


  DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi);

  if( h.GetHadronType() == KAON && icharge == 0 ){ // Our assumptions
    // for K0_S = 1/sqrt{2} {d {bar s} + {bar d} s) !!!!

//    double up_in_kplus   = dff[1+5];
//    double down_in_kplus = dff[2+5];
//    double sbar_in_kplus = dff[-3+5];


    double SB_Kp, DB_Kp, UB_Kp, GL_Kp, U_Kp, D_Kp, S_Kp,
      SB_Km, DB_Km, UB_Km, GL_Km, U_Km, D_Km, S_Km;


    icharge = positive;
    DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi); // K^+ is produced!!!


    SB_Kp  = dff[5-3];
    DB_Kp  = dff[5-2];
    UB_Kp  = dff[5-1];
    GL_Kp  = dff[5+0];
    U_Kp   = dff[5+1];
    D_Kp   = dff[5+2];
    S_Kp   = dff[5+3];

    icharge = negative;
    DLIB(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi); // K^- is produced!!!

    SB_Km  = dff[5-3];
    DB_Km  = dff[5-2];
    UB_Km  = dff[5-1];
    GL_Km  = dff[5+0];
    U_Km   = dff[5+1];
    D_Km   = dff[5+2];
    S_Km   = dff[5+3];




//????????????????????????????????

    dff[5-5] = 0.  ;
    dff[5-4] = 0.  ;
    dff[5-3] = 0.5 * ( SB_Kp + SB_Km ); // D_{s-bar}^{K0} =   1/2[ D_{s-bar}^{K+} + D_{s-bar}^{K-}] ]
    dff[5-2] = 0.5 * ( UB_Kp + UB_Km ); // D_{d-bar}^{K0} =  1/2[ D_{u-bar}^{K+} + D_{u-bar}^{K-}] ] ???
    dff[5-1] = 0.5 * ( DB_Kp + D_Km );  // D_{u-bar}^{K0} =  1/2[ D_{d-bar}^{K+} + D_{d}^{K-}] ]     ???
    dff[5+0] = 0.5 * ( GL_Kp + GL_Km );
    dff[5+1] = 0.5 * ( D_Kp  + DB_Km  );// D_{u}^{K0} =   1/2[ D_{d}^{K+} + D_{d-bar}^{K-}] ]
    dff[5+2] = 0.5 * ( U_Kp  + U_Km  ); // D_{d}^{K0} =   1/2[ D_{u}^{K+} + D_{u-bar}^{K-}] ]
    dff[5+3] = 0.5 * ( S_Kp  + S_Km  ); // D_{s}^{K0} =   1/2[ D_{s}^{K+} + D_{s}^{K-}] ]
    dff[5+4] = 0.  ;
    dff[5+5] = 0.  ;


//     fragmentation.dff[1+5]  = down_in_kplus;  // up --> K^0   = down --> K^+
//     fragmentation.dff[2+5]  = up_in_kplus;    // down --> K^0 = up --> K^+
//     fragmentation.dff[-2+5] = up_in_kplus;    // dsea --> K^0 = up --> K^+
//     fragmentation.dff[3+5]  = sbar_in_kplus;  // s --> K^0    = sbar --> K^+
//     fragmentation.dff[-3+5] = sbar_in_kplus;  // sbar --> K^0 = sbar --> K^+

    icharge = neutral;  // STATUS QUO IS RESTORED.
   }

  }

// * iparton = 5,4,3,2,1,0,-1,  -2,  -3,  -4,  -5 means
//             b,c,s,d,u,g,ubar,dbar,sbar,cbar,bbar      *

    h.fragmentation.up       = dff[1+5];
    h.fragmentation.down     = dff[2+5];
    h.fragmentation.anti_up    = dff[-1+5];
    h.fragmentation.anti_down  = dff[-2+5];
    h.fragmentation.strange      = dff[3+5];
    h.fragmentation.anti_strange   = dff[-3+5];
    h.fragmentation.charm      = dff[4+5];
    h.fragmentation.anti_charm   = dff[-4+5];
    h.fragmentation.glu      = dff[0+5];
    h.fragmentation.bottom      = dff[5+5];
    h.fragmentation.anti_bottom   = dff[-5+5];
    h.fragmentation.top      = 0.;
    h.fragmentation.anti_top   = 0.;
}





//========================================================= Unpolarized fragmentation
void sidis::Fragmentation(HADRON& h, double z, double Q2)
  // Returns the parton content at x,Q2
{
   // icharge 0,1,2,3 is 0,+,-,+&-
  double charge = 1;

  if( h.GetCharge() == neutral || h.GetCharge() == 0. ){ charge = 0;} // the same as 3
  else if( h.GetCharge() == positive || h.GetCharge() > 0.  ){ charge = 1;}
  else if( h.GetCharge() == negative || h.GetCharge() < 0. ){ charge = 2;}

  fragmentation(h, charge, z, Q2);


};




//========================================================= fragmentation functions ala transversity
void sidis::fragmentationCollins(HADRON& h, int charge, double z, double Q2)
  // Returns fragmentation content at z,Q2
{
     fragmentation(h,  charge,   z,  Q2);
}


//========================================================= fragmentation Collins NLO with Feng 2014
void sidis::FragmentationCollinsNLO(HADRON& h, double z, double Q2)
  // Returns fragmentation content at z,Q2
{
    // Note ethat it works only for pions now!

     double pdf_a[13];
  hoppetEval(z, sqrt(Q2), pdf_a); // evaluate transversity&collins for pi+ u--> pi+, down --> pi+
//HOPPET GIVES Z * D(Z) !

    h.fragmentation.up       =  pdf_a[1+6]/z;
    h.fragmentation.down     =  pdf_a[-1+6]/z;
    h.fragmentation.anti_up    = pdf_a[-1+6]/z;
    h.fragmentation.anti_down  = pdf_a[1+6]/z;
    h.fragmentation.strange      = pdf_a[3+6]/z;   // strange functions are evolved separately
    h.fragmentation.anti_strange   = pdf_a[3+6]/z; // strange functions are evolved separately
    h.fragmentation.charm      = 0.;
    h.fragmentation.anti_charm   = 0.;
    h.fragmentation.glu      = 0.;
    h.fragmentation.bottom      = 0.;
    h.fragmentation.anti_bottom   = 0.;
    h.fragmentation.top      = 0.;
    h.fragmentation.anti_top   = 0.;

}



//========================================================= Collins fragmentation with Feng 2014....
// Work on it !!!!!!!!!!!!!!!!!!!!!
void sidis::CollinsDistributionNLO( HADRON& h, double z, double Q2)
  // Returns the parton content at x,Q2
{


    FragmentationCollinsNLO(h, z, Q2); // u --> pi+, down --> pi+


    PARTONCONTENT coll = h.fragmentation;
    h.fragmentation.up = 1.;
    h.fragmentation.down = 1.;
    h.fragmentation.anti_up = 1.;
    h.fragmentation.anti_down = 1.;
    h.fragmentation.strange = 1.;
    h.fragmentation.anti_strange = 1.;
    h.fragmentation.charm = 0.;
    h.fragmentation.anti_charm = 0.;
    h.fragmentation.bottom= 0.;
    h.fragmentation.anti_bottom = 0.;
    h.fragmentation.glu = 0.;



    collins_distribution(h, coll); // This function does correct "reshuffling" -- check the logic!

};

//========================================================= Collins fragmentation ala transversity
void sidis::FragmentationCollins(HADRON& h, double z, double Q2)
  // Returns the parton content at x,Q2
{
   // icharge 0,1,2,3 is 0,+,-,+&-
  double charge = 1;

  if( h.GetCharge() == neutral || h.GetCharge() == 0. ){ charge = 0;} // the same as 3
  else if( h.GetCharge() == positive || h.GetCharge() > 0.  ){ charge = 1;}
  else if( h.GetCharge() == negative || h.GetCharge() < 0. ){ charge = 2;}

  fragmentationCollins(h, charge, z, Q2);
};




//========================================================= Sivers effect...x dependence
double sidis::collins_z_dependence(double z, double a, double b, double n)
{

  double collins_z = 0.;
  if( b == 0. && a == 0. )
    collins_z = n;
  else if( b == 0. )
    collins_z = n * pow(z,a);
  else if( a == 0. )
    collins_z = n * pow((1.-z),b);
  else if( a < 0. || b < 0.)
    collins_z = n * pow(z,a) * pow((1.-z),b);
  else collins_z = n * pow(z,a) * pow((1.-z),b) * pow((a+b),(a+b)) /(pow(a,a) * pow(b,b));

//    collins_z = n * pow(z,a) * pow((1.-z),b); // NOTE TEST KANG ONLY PP FITTING!!

  return collins_z;
}



//========================================================= Collins fragmentation functions
PARTONCONTENT sidis::collins( double z )
  // Returns fragmentation content at z
{

  // Defaults are for pi^+(u, bard)
  // u --> pi^+ favoured
  // d --> pi^+ unfavoured
  PARTONCONTENT partcontent;


// IF PIONS OR HADRONS
if( Produced.GetHadronType() == PION || Produced.GetHadronType() == HADRONS){

  partcontent.up  = 2.*
    collins_z_dependence(z, Params.Collins.parameters.a_up,
			Params.Collins.parameters.b_up,
			Params.Collins.parameters.n_up); // Favoured fragmentation

// cout << Params.Collins.parameters.n_up << " n_up" << endl;

  partcontent.down  = 2.*
    collins_z_dependence(z, Params.Collins.parameters.a_down,
			Params.Collins.parameters.b_down,
			Params.Collins.parameters.n_down); // UnFavoured fragmentation

  partcontent.anti_up = 2.*
    collins_z_dependence(z, Params.Collins.parameters.a_anti_up,
			Params.Collins.parameters.b_anti_up,
			Params.Collins.parameters.n_anti_up); // UnFavoured fragmentation

  partcontent.anti_down = 2.*
    collins_z_dependence(z, Params.Collins.parameters.a_anti_down,
			Params.Collins.parameters.b_anti_down,
			Params.Collins.parameters.n_anti_down); //  AntiFavoured fragmentation = Fav1



  partcontent.strange = 2.*
    collins_z_dependence(z, Params.Collins.parameters.a_strange,
			Params.Collins.parameters.b_strange,
			Params.Collins.parameters.n_strange); // StrangeUnFavoured fragmentation


  partcontent.anti_strange = 2.*
    collins_z_dependence(z, Params.Collins.parameters.a_anti_strange,
			Params.Collins.parameters.b_anti_strange,
			Params.Collins.parameters.n_anti_strange); //    StrangeUnFavoured fragmentation
}


// IF KAON
if( Produced.GetHadronType() == KAON ){

  partcontent.up  = 2.*
    collins_z_dependence(z, Params.CollinsKaon.parameters.a_up,
			Params.CollinsKaon.parameters.b_up,
			Params.CollinsKaon.parameters.n_up); // Favoured fragmentation


  partcontent.down  = 2.*
    collins_z_dependence(z, Params.CollinsKaon.parameters.a_down,
			Params.CollinsKaon.parameters.b_down,
			Params.CollinsKaon.parameters.n_down); // UnFavoured fragmentation

  partcontent.anti_up = 2.*
    collins_z_dependence(z, Params.CollinsKaon.parameters.a_anti_up,
			Params.CollinsKaon.parameters.b_anti_up,
			Params.CollinsKaon.parameters.n_anti_up); // UnFavoured fragmentation

  partcontent.anti_down = 2.*
    collins_z_dependence(z, Params.CollinsKaon.parameters.a_anti_down,
			Params.CollinsKaon.parameters.b_anti_down,
			Params.CollinsKaon.parameters.n_anti_down); //  AntiFavoured fragmentation = Fav1



  partcontent.strange = 2.*
    collins_z_dependence(z, Params.CollinsKaon.parameters.a_strange,
			Params.CollinsKaon.parameters.b_strange,
			Params.CollinsKaon.parameters.n_strange); // StrangeFavoured fragmentation


  partcontent.anti_strange = 2.*
    collins_z_dependence(z, Params.CollinsKaon.parameters.a_anti_strange,
			Params.CollinsKaon.parameters.b_anti_strange,
			Params.CollinsKaon.parameters.n_anti_strange); //    StrangeUnFavoured fragmentation
}

  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;

  return partcontent;

}




// USEFUL FUNCTION TO RE_SHUFFLE THINGS FOR COLLINS DISTRIBUTION...
void sidis::collins_distribution( HADRON& h, PARTONCONTENT a){
  // reshuffle fragmentation according to hadron type
  // a is Collins function always defined for pi^+ up = favoured, down = unfavoured etc
//  double temp = 0;


  switch( h.GetHadronType() ){
  case HADRONS:
  case PION:
    //fragmentation(h, 1, z, Q2); // It works badly DLIB reads grid file each time I change the sign
    // a solution "not elegant" is to u <-> d two times: so that Fav and Unfav Collins are realised...



//    if( h.GetCharge() == neutral || h.GetCharge() == 0. ){// NEUTRAL pi0 = 1/2 (pi+ + pi-) !!!
//      // pi = 1/2 (pi+ + pi-)
//      // For Neutral we need 1/2 (Hfav + Hunf)
//      h.fragmentation.up       -= 0.5*h.fragmentation.strange; // 1/2 (FAV + UNFAV)  - 1/2 UNF(strange = UNF)
//      h.fragmentation.down      = 0.5*h.fragmentation.strange; //  1/2 UNF(strange = UNF)
//      h.fragmentation.anti_up   = h.fragmentation.down; //  1/2 UNF(strange = UNF)
//      h.fragmentation.anti_down = h.fragmentation.up; // 1/2 FAV
//    } // CHECK IT

// NOTE Neutral pion we get fragmentation quark->pi0 = 1/2 (quark->pi+ + quark->pi-)
// NOTE we need the following combination collins_pi0 = 1/2 (collins_pi+ + collins_pi-)
// NOTE  how to do it??????????????????


    if( h.GetCharge() == negative || h.GetCharge() < 0. ){ // pi- (u-bar d)
      // a here is collins(z)
      double fav           = a.up;
      double antifav       = a.anti_down;
      double unfav         = a.down;
      double unfav_strange = a.strange; // MAIN TORINO

      a.anti_up = antifav; // a.up --> favoured Collins MAIN TORINO
      a.down    = fav; // a.up --> favoured Collins

      // all the rest unfavoured:
      a.up           = unfav; // a.down --> unfavoured Collins
      a.anti_down    = unfav; // a.down --> unfavoured Collins
      a.strange      = unfav_strange; // a.down --> unfavoured Collins
      a.anti_strange = unfav_strange; // a.down --> unfavoured Collins


// NOTE test stefano only fav and unfav
//       a.anti_up = fav; // a.up --> favoured Collins
//       a.down    = fav; // a.up --> favoured Collins
//
//       // all the rest unfavoured:
//       a.up           = unfav; // a.down --> unfavoured Collins
//       a.anti_down    = unfav; // a.down --> unfavoured Collins
//       a.strange      = unfav; // a.down --> unfavoured Collins
//       a.anti_strange = unfav; // a.down --> unfavoured Collins

    }

    if( h.GetCharge() == positive || h.GetCharge() > 0. ){ // pi+ (u d-bar)
      // a here is collins(z)
      double fav           = a.up;
      double antifav       = a.anti_down;
      double unfav         = a.down;
      double unfav_strange = a.strange; // MAIN TORINO

      a.up           = fav; // a.up --> favoured Collins
      a.anti_down    = antifav; // a.up --> favoured Collins MAIN TORINO

      // all the rest unfavoured:
      a.anti_up      = unfav; // a.down --> unfavoured Collins
      a.down         = unfav; // a.down --> unfavoured Collins
      a.anti_strange = unfav_strange; // a.down --> unfavoured Collins
      a.strange      = unfav_strange; // a.down --> unfavoured Collins

// NOTE test stefano only fav and unfav
//       a.up           = fav; // a.up --> favoured Collins
//       a.anti_down    = fav; // a.up --> favoured Collins
//
//       // all the rest unfavoured:
//       a.anti_up      = unfav; // a.down --> unfavoured Collins
//       a.down         = unfav; // a.down --> unfavoured Collins
//       a.strange      = unfav; // a.down --> unfavoured Collins
//       a.anti_strange = unfav; // a.down --> unfavoured Collins


    }



    if( h.GetCharge() == neutral || h.GetCharge() == 0. ){ // pi0 = 1/2 (pi^+ + pi^-)
      // a here is collins(z)
      double fav           = 0.5 * (a.up + a.down);
//      double unfav         = 0.5 * (a.up + a.down);

      double unfav_strange = a.strange;




      a.up           = fav; //
      a.anti_down    = fav; //
      a.anti_up      = fav; //
      a.down         = fav; //

      a.anti_strange = unfav_strange; //  //NOTE TEST IT !!!!!!!!!
      a.strange      = unfav_strange; //  //NOTE TEST IT !!!!!!!!!
     }

    h.fragmentation = product( a, h.fragmentation); // IMPORTANT STEP!!!!




//    if( h.GetCharge() == neutral || h.GetCharge() == 0. ){// NEUTRAL pi0 = 1/2 (pi+ + pi-) !!!
//      // NOW we have up = 1/2 hfav, down = 1/2 Hunf etc, let's sum them:
//      h.fragmentation.up       += h.fragmentation.down;
//      h.fragmentation.down      = h.fragmentation.up; //
//      h.fragmentation.anti_up   = h.fragmentation.up; //
//      h.fragmentation.anti_down = h.fragmentation.up; //

//      h.fragmentation.strange      = h.fragmentation.strange;
//      h.fragmentation.anti_strange = h.fragmentation.strange;

//       // OR ONLY FAVOURED!!!
//       h.fragmentation.up    += h.fragmentation.up; //
//       h.fragmentation.down   = h.fragmentation.up; //
//       h.fragmentation.anti_up= h.fragmentation.up; //
//       h.fragmentation.anti_down= h.fragmentation.up; //

//    } // CHECK IT


    break;
   case KAON:
//   we suppose that
// H1^bar s/D^bar s| K+ = H1^u/D^u| K+ = H1^u/D^u| pi+
//
    //fragmentation(h, 1, z, Q2); // It works badly DLIB reads grid file each time I change the sign
    // a solution "not elegant" is to u <-> d two times: so that Fav and Unfav Collins are realised...
    // K^+ u, sbar --> favoured, all the rest --> unfavoured

    if( h.GetCharge() == negative || h.GetCharge() < 0. ){ // K- (u-bar s)
      double fav   = a.up;
      double unfav = a.down;
      double fav_strange = a.strange;
      double unfav_strange = a.anti_strange;
      //double fav_strange = fav;
      //double unfav_strange = unfav;

      a.anti_up = fav; // a.up --> favoured Collins
      a.strange = fav_strange; // a.up --> favoured Collins

      // all the rest unfavoured:
      a.up           = unfav; // a.down --> unfavoured Collins
      a.down         = unfav; // a.down --> unfavoured Collins
      a.anti_down    = unfav; // a.down --> unfavoured Collins
      a.anti_strange = unfav_strange; // a.down --> unfavoured Collins
    }

    if( h.GetCharge() == positive || h.GetCharge() > 0. ){ // K+ (u s-bar)
      double fav   = a.up;
      double unfav = a.down;
      double fav_strange = a.strange;
      double unfav_strange = a.anti_strange;
      //double fav_strange = fav;
      //double unfav_strange = unfav;

      a.up           = fav; // a.up --> favoured Collins
      a.anti_strange = fav_strange; // a.up --> favoured Collins

      // all the rest unfavoured:
      a.anti_up      = unfav; // a.down --> unfavoured Collins
      a.down         = unfav; // a.down --> unfavoured Collins
      a.anti_down    = unfav; // a.down --> unfavoured Collins
      a.strange      = unfav_strange; // a.down --> unfavoured Collins
    }


// for K0_S = 1/sqrt{2} {d {bar s} + {bar d} s) !!!!
// IMPLEMENT K0 NOT DONE YET !!!!!!!!!!!!!!!!!!!!!!!


    h.fragmentation = product( a, h.fragmentation); // IMPORTANT STEP!!!!

    break;
  default:
    h.fragmentation.up      = 0.;
    h.fragmentation.down    = 0.;
    h.fragmentation.anti_up   = 0.;
    h.fragmentation.anti_down = 0.;
    h.fragmentation.strange     = 0.;
    h.fragmentation.anti_strange  = 0.;
    h.fragmentation.charm     = 0.;
    h.fragmentation.anti_charm  = 0.;
    h.fragmentation.bottom     = 0.;
    h.fragmentation.anti_bottom  = 0.;
    h.fragmentation.top     = 0.;
    h.fragmentation.anti_top  = 0.;
    break;
  }


}

//========================================================= Collins fragmentation....
// Work on it !!!!!!!!!!!!!!!!!!!!!
void sidis::CollinsDistribution( HADRON& h, double z, double Q2)
  // Returns the parton content at x,Q2
{

    //Fragmentation(h, z, Q2);
    FragmentationCollins(h, z, Q2);
    collins_distribution(h, collins(z));
    //if( h.GetHadronType() == PION || h.GetHadronType() == HADRONS) collins_distribution(h, collins(z));
    //if( h.GetHadronType() == KAON) collins_distribution(h, collinsKaon(z));

};



//========================================================= Collins fragmentation functions kt integrated dependence FIRST MOMENT d2pt pt^2/2zM_h^2 H1^\perp(z,pt^2) ...
PARTONCONTENT sidis::collinsFirstMoment( double z )
  // Returns fragmentation content at z
{

  // Defaults are for pi^+(u, bard)
  // u --> pi^+ favoured
  // d --> pi^+ unfavoured
  PARTONCONTENT partcontent;

  partcontent = collins( z ); // Collins z dependence...

  double ptq2_average = Params.GetPtq2Average();


// IF PIONS OR HADRONS
if( Produced.GetHadronType() == PION || Produced.GetHadronType() == HADRONS){


  double collins_width_up = pow(Params.Collins.parameters.m2_up,1.5) *ptq2_average/pow2(ptq2_average + Params.Collins.parameters.m2_up);

  double integrand_up = sqrt( e/2. )/2. * collins_width_up / (z * mpion); // see  SIDIS_Collins.nb
// /2. because collins(z) gives 2 N z^a(1-z)^b

  double collins_width_down = pow(Params.Collins.parameters.m2_down,1.5)*ptq2_average/pow2(ptq2_average + Params.Collins.parameters.m2_down);

  double integrand_down = sqrt( e/2. )/2. * collins_width_down / (z * mpion);

  double collins_width_anti_up = pow(Params.Collins.parameters.m2_anti_up,1.5)*ptq2_average/pow2(ptq2_average + Params.Collins.parameters.m2_anti_up);

  double integrand_anti_up = sqrt( e/2. )/2. *   collins_width_anti_up / (z * mpion );

  double collins_width_anti_down = pow(Params.Collins.parameters.m2_anti_down,1.5)*ptq2_average/pow2(ptq2_average + Params.Collins.parameters.m2_anti_down);

  double integrand_anti_down = sqrt( e/2. )/2. *  collins_width_anti_down  / (z * mpion);

  double collins_width_strange = pow(Params.Collins.parameters.m2_strange,1.5)*ptq2_average/pow2(ptq2_average + Params.Collins.parameters.m2_strange);

  double integrand_strange = sqrt( e/2. )/2. *  collins_width_strange  / (z * mpion);

  double collins_width_anti_strange = pow(Params.Collins.parameters.m2_anti_strange,1.5)*ptq2_average/pow2(ptq2_average + Params.Collins.parameters.m2_anti_strange);

  double integrand_anti_strange = sqrt( e/2. )/2. *  collins_width_anti_strange  / (z * mpion);


  partcontent.up  *=  integrand_up; //

  partcontent.down  *= integrand_down; //

  partcontent.anti_up  *=  integrand_anti_up; //

  partcontent.anti_down  *= integrand_anti_down; //

  partcontent.strange *= integrand_strange; //

  partcontent.anti_strange *= integrand_anti_strange; //
  }

// IF KAON
if( Produced.GetHadronType() == KAON){


  double collins_width_up = pow(Params.CollinsKaon.parameters.m2_up,1.5) *ptq2_average/pow2(ptq2_average + Params.CollinsKaon.parameters.m2_up);

  double integrand_up = sqrt( e/2. )/2. * collins_width_up / (z * mpion); // see  SIDIS_Collins.nb
// /2. because collins(z) gives 2 N z^a(1-z)^b

  double collins_width_down = pow(Params.CollinsKaon.parameters.m2_down,1.5)*ptq2_average/pow2(ptq2_average + Params.CollinsKaon.parameters.m2_down);

  double integrand_down = sqrt( e/2. )/2. * collins_width_down / (z * mpion);

  double collins_width_anti_up = pow(Params.CollinsKaon.parameters.m2_anti_up,1.5)*ptq2_average/pow2(ptq2_average + Params.CollinsKaon.parameters.m2_anti_up);

  double integrand_anti_up = sqrt( e/2. )/2. *   collins_width_anti_up / (z * mpion );

  double collins_width_anti_down = pow(Params.CollinsKaon.parameters.m2_anti_down,1.5)*ptq2_average/pow2(ptq2_average + Params.CollinsKaon.parameters.m2_anti_down);

  double integrand_anti_down = sqrt( e/2. )/2. *  collins_width_anti_down  / (z * mpion);

  double collins_width_strange = pow(Params.CollinsKaon.parameters.m2_strange,1.5)*ptq2_average/pow2(ptq2_average + Params.CollinsKaon.parameters.m2_strange);

  double integrand_strange = sqrt( e/2. )/2. *  collins_width_strange  / (z * mpion);

  double collins_width_anti_strange = pow(Params.CollinsKaon.parameters.m2_anti_strange,1.5)*ptq2_average/pow2(ptq2_average + Params.CollinsKaon.parameters.m2_anti_strange);

  double integrand_anti_strange = sqrt( e/2. )/2. *  collins_width_anti_strange  / (z * mpion);


  partcontent.up  *=  integrand_up; //

  partcontent.down  *= integrand_down; //

  partcontent.anti_up  *=  integrand_anti_up; //

  partcontent.anti_down  *= integrand_anti_down; //

  partcontent.strange *= integrand_strange; //

  partcontent.anti_strange *= integrand_anti_strange; //
  }


  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;

  return partcontent;

}


//========================================================= Collins FF First moment
void sidis::CollinsDistributionFirstMoment( HADRON& h, double z, double Q2)
// Returns the parton content at x,Q2
{

    //Fragmentation(h, z, Q2);
    FragmentationCollins(h, z, Q2);
    //h.fragmentation.up = 1;
  collins_distribution( h, collinsFirstMoment( z ));

};




 // transversity....
//========================================================= transversity...x dependence
double sidis::transversity_x_dependence(double x, double a, double b, double n, double gamma, double eta)
{

  double trans_x = 0.;
  if( n == 0. )
    trans_x = 0.;
  else if( b == 0. && a == 0. )
    trans_x = n;
  else if( b == 0. )
    trans_x = n * pow(x,a);
  else if( a == 0. )
    trans_x = n * pow((1.-x),b);
  else
    trans_x = n * pow(x,a) * pow((1.-x),b) * pow((a+b),(a+b)) /(pow(a,a) * pow(b,b)); // MAIN TORINO GROUP



  return trans_x;
}

//=========================================================Soffer Bound partcontent
PARTONCONTENT sidis::SofferBound(double x, double Q2)
  // Returns the parton content at x,Q2
{
  PARTONCONTENT partcontent;

//   int iset = 2; //  LEADING ORDER
  int iset = transversity_order; //  LEADING ORDER
  double  u = 0., d = 0., ub = 0., db = 0., st = 0.;

  PARPOLT(iset, x, Q2, u, d, ub, db, st); // We invoke x*h1(x,Q)

  //Assign parton densities PROTON:
  if(!VALENCE_SIDIS){
  partcontent.up           =  u / x; //
  partcontent.down         =  (-1.) * d / x; //
  partcontent.anti_up      =  ub / x; //
  partcontent.anti_down    =  (-1.) * db / x; //
  partcontent.strange      =  st / x; //
  partcontent.anti_strange =  st / x; //
  partcontent.charm        =  0.; //
  partcontent.anti_charm   =  0.; //
  partcontent.bottom       =  0.; //
  partcontent.anti_bottom  =  0.; //
  partcontent.top          =  0.; //
  partcontent.anti_top     =  0.; //
  partcontent.glu          =  0.; //
  } else if(VALENCE_SIDIS){
  partcontent.up           =  (u - ub) / x; //
  partcontent.down         =  (-1.) * (d - db) / x; //
  partcontent.anti_up      =  ub / x; //
  partcontent.anti_down    =  (-1.) * db / x; //
  partcontent.strange      =  st / x; //
  partcontent.anti_strange =  st / x; //
  partcontent.charm        =  0.; //
  partcontent.anti_charm   =  0.; //
  partcontent.bottom       =  0.; //
  partcontent.anti_bottom  =  0.; //
  partcontent.top          =  0.; //
  partcontent.anti_top     =  0.; //
  partcontent.glu          =  0.; //
  }

  return partcontent;

};


//========================================================= transversity...
PARTONCONTENT sidis::transversity(double x ) {
  // Returns the parton content for Sivers effect...
  PARTONCONTENT partcontent;

  partcontent.up  =
  transversity_x_dependence(x, Params.Transversity.parameters.a_up,
      		Params.Transversity.parameters.b_up,
      		Params.Transversity.parameters.n_up,
      		Params.Transversity.parameters.gamma_up,
      		Params.Transversity.parameters.eta_up);


  partcontent.down  =
  transversity_x_dependence(x, Params.Transversity.parameters.a_down,
      		Params.Transversity.parameters.b_down,
      		Params.Transversity.parameters.n_down,
      		Params.Transversity.parameters.gamma_down,
      		Params.Transversity.parameters.eta_down); //

  partcontent.anti_up =
  transversity_x_dependence(x, Params.Transversity.parameters.a_anti_up,
      		Params.Transversity.parameters.b_anti_up,
      		Params.Transversity.parameters.n_anti_up,
      		Params.Transversity.parameters.gamma_anti_up,
      		Params.Transversity.parameters.eta_anti_up); //

  partcontent.anti_down =
  transversity_x_dependence(x, Params.Transversity.parameters.a_anti_down,
      		Params.Transversity.parameters.b_anti_down,
      		Params.Transversity.parameters.n_anti_down,
      		Params.Transversity.parameters.gamma_anti_down,
      		Params.Transversity.parameters.eta_anti_down); //


  partcontent.strange =
  transversity_x_dependence(x, Params.Transversity.parameters.a_strange,
      		Params.Transversity.parameters.b_strange,
      		Params.Transversity.parameters.n_strange,
      		Params.Transversity.parameters.gamma_strange,
      		Params.Transversity.parameters.eta_strange); //


  partcontent.anti_strange =
  transversity_x_dependence(x, Params.Transversity.parameters.a_anti_strange,
      		Params.Transversity.parameters.b_anti_strange,
      		Params.Transversity.parameters.n_anti_strange,
      		Params.Transversity.parameters.gamma_anti_strange,
      		Params.Transversity.parameters.eta_anti_strange); //

  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;

  return partcontent;

};


//========================================================= Transversity partcontent
void sidis::TransversityDistribution( HADRON& h, double x, double Q2)
  // Returns the parton content at x,Q2
{

  switch( h.GetHadronType() ){
  case PROTON:
    h.partcontent = SofferBound( x, Q2);
    h.partcontent = product(transversity(x), h.partcontent);
    break;
  case ANTIPROTON:
    h.partcontent = SofferBound( x, Q2);
    h.partcontent = product(transversity(x), h.partcontent);
    Antiproton(h.partcontent);
    break;
  case PION:
    break;
  case DEUTRON:
    h.partcontent = SofferBound( x, Q2);
    h.partcontent = product(transversity(x), h.partcontent);
    Deutron(h.partcontent);
    break;
  case NEUTRON:
    h.partcontent = SofferBound( x, Q2);
    h.partcontent = product(transversity(x), h.partcontent);
    Neutron(h.partcontent);
    break;
  default:
    break;
  }

};


//========================================================= Transversity partcontent using hoppen at NLO!!! FENG 2014
void sidis::TransversityDistributionNLO( HADRON& h, double x, double Q2)
  // Returns the parton content at x,Q2
{


  double pdf_a[13];
  hoppetEval(x, sqrt(Q2), pdf_a); // evaluate transversity

  h.partcontent.up     = pdf_a[6+2]/x; // UP, see pdf_init
  //h.partcontent.anti_up  = 0.; // MAIN 2015 no sea quarks
  h.partcontent.anti_up  = pdf_a[6+4]/x; // MAIN 2016 - sea quarks


  h.partcontent.down     = pdf_a[6-2]/x; // DOWN, see pdf_init
  //h.partcontent.anti_down  = 0.; // MAIN 2015 no sea quarks
  h.partcontent.anti_down  = pdf_a[6-4]/x; // MAIN 2016 - sea quarks

  h.partcontent.strange     = 0.;
  h.partcontent.anti_strange  = 0.;

  h.partcontent.charm     = 0.;
  h.partcontent.anti_charm  = 0.;
  h.partcontent.bottom     = 0.;
  h.partcontent.anti_bottom  = 0.;
  h.partcontent.top     = 0.;
  h.partcontent.anti_top  = 0.;
  h.partcontent.glu   = 0.;


  switch( h.GetHadronType() ){
  case PROTON:
    break;
  case ANTIPROTON:
    Antiproton(h.partcontent);
    break;
  case PION:
    break;
  case DEUTRON:
    Deutron(h.partcontent);
    break;
  case NEUTRON:
     Neutron(h.partcontent);
    break;
  default:
    break;
  }

};






//========================================================= Transversity partcontent First moment
void sidis::TransversityDistributionFirstMoment( HADRON& h, double x, double Q2)
// Returns the parton content at x,Q2
{
  TransversityDistribution( h, x, Q2);
};



//========================================================= Sivers partcontent First moment
void sidis::F1TDistributionFirstMoment( HADRON& h,  double x, double Q2)
// Returns the parton content at x,Q2
{
  SiversDistributionFirstMoment( h,  x, Q2);
//   Relation OUr first moment = -1. f1tperp(1)
  h.partcontent = product(-1., h.partcontent);

};



//========================================================= Sivers partcontent First moment
void sidis::SiversDistributionFirstMoment( HADRON& h,  double x, double Q2)
// Returns the parton content at x,Q2
{
  h.partcontent = unpolarised(x, Q2);



  h.partcontent = product( sivers(x), h.partcontent);



//  double integrand = 1;
//   double mpr = 0.93827203 ; // the proton mass
  double kt2_average = Params.GetKt2Average();



  h.partcontent.up  *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_up, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_up +kt2_average, 2); // Eq (22.) hep-ph/0507181

  h.partcontent.down  *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_down, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_down +kt2_average, 2); //

  h.partcontent.anti_up *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_up, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_up +kt2_average, 2); //

  h.partcontent.anti_down *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_anti_down, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_anti_down +kt2_average, 2); //


  h.partcontent.strange *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_strange, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_strange +kt2_average, 2); //


  h.partcontent.anti_strange *= sqrt( 2.*exp(1.) )/( 4.*mpr ) * pow(Params.Sivers.parameters.m2_anti_strange, 3./2.)*kt2_average/
    pow(Params.Sivers.parameters.m2_anti_strange +kt2_average, 2); //


  h.partcontent.charm     = 0.;
  h.partcontent.anti_charm  = 0.;
  h.partcontent.bottom     = 0.;
  h.partcontent.anti_bottom  = 0.;
  h.partcontent.top     = 0.;
  h.partcontent.anti_top  = 0.;
  h.partcontent.glu  = 0.;




  switch( h.GetHadronType() ){
  case PROTON:
    break;
  case ANTIPROTON:
    Antiproton(h.partcontent);
    break;
  case PION:
    break;
  case DEUTRON:
    Deutron(h.partcontent);
    break;
  case NEUTRON:
    Neutron(h.partcontent);
    break;
  default:
    break;
  }

};



//========================================================= Sivers effect...
PARTONCONTENT  sidis::sivers(double x ) {
  // Returns the parton content for Sivers effect...
  // // // //TEST 2 BROKEN SEA MAIN PAPER 07...

  // Delta N_Sivers < 2. * f1(x)
  PARTONCONTENT partcontent;

  partcontent.up  = 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_up,
			Params.Sivers.parameters.b_up,
			Params.Sivers.parameters.n_up,
			Params.Sivers.parameters.gamma_up,
			Params.Sivers.parameters.eta_up);

  //  partcontent.up  *= sqrt( Params.Sivers.parameters.m2_up ); // test Stefano

  partcontent.down  = 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_down,
			Params.Sivers.parameters.b_down,
			Params.Sivers.parameters.n_down,
			Params.Sivers.parameters.gamma_down,
			Params.Sivers.parameters.eta_down); //

  //  partcontent.down  *= sqrt( Params.Sivers.parameters.m2_down ); // test Stefano


  partcontent.anti_up = 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_anti_up,
			Params.Sivers.parameters.b_anti_up,
			Params.Sivers.parameters.n_anti_up ,
			Params.Sivers.parameters.gamma_anti_up,
			Params.Sivers.parameters.eta_anti_up); //

  //  partcontent.anti_up  *= sqrt( Params.Sivers.parameters.m2_anti_up ); // test Stefano


  partcontent.anti_down = 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_anti_down,
			Params.Sivers.parameters.b_anti_down,
			Params.Sivers.parameters.n_anti_down,
			Params.Sivers.parameters.gamma_anti_down,
			Params.Sivers.parameters.eta_anti_down); //

  //  partcontent.anti_down  *= sqrt( Params.Sivers.parameters.m2_anti_down ); // test Stefano

  partcontent.strange = 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_strange,
			Params.Sivers.parameters.b_strange,
			Params.Sivers.parameters.n_strange,
			Params.Sivers.parameters.gamma_strange,
			Params.Sivers.parameters.eta_strange); //


  //  partcontent.strange  *= sqrt( Params.Sivers.parameters.m2_strange ); // test Stefano

  partcontent.anti_strange = 2.*
    sivers_x_dependence(x, Params.Sivers.parameters.a_anti_strange,
			Params.Sivers.parameters.b_anti_strange,
			Params.Sivers.parameters.n_anti_strange,
			Params.Sivers.parameters.gamma_anti_strange,
			Params.Sivers.parameters.eta_anti_strange); //

  //  partcontent.anti_strange  *= sqrt( Params.Sivers.parameters.m2_anti_strange ); // test Stefano

  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;

  return partcontent;

};


//========================================================= Sivers effect...x dependence
double sidis::sivers_x_dependence(double x, double a, double b, double n, double gamma, double eta)
{

  double sivers_x = 0.;
   if( n == 0. )
     sivers_x = 0.;
   else if( b == 0. && a == 0. )
     sivers_x = n;
   else if( b == 0. )
     sivers_x = n * pow(x,a);
   else if( a == 0. )
     sivers_x = n * pow((1.-x),b);
    else if( a < 0. || b < 0. )
      sivers_x = n * pow(x,a) * pow((1.-x),b);
    else
      sivers_x = n * pow(x,a) * pow((1.-x),b) * pow((a+b),(a+b)) /(pow(a,a) * pow(b,b)); // NOTE MAIN TORINO GROUP


//    sivers_x = n * pow(x,a) * pow((1.-x),b); //NOTE TEST MAIN JET AN NOTE NOTE NOTE test KANG assumption SEATTLE ...



  //sivers with nodes!

  //sivers_x += gamma * pow( (1-x), 7.); // SIVERS WITH NODES test SEATTLE...

  //  sivers_x *= tanh( (x - gamma) * 10. ); // test SEATTLE we fix the turnover point as gamma ... SPEAK WITH DANIEL -NODE ON SIVERS FUNCTION !!!!
//     sivers_x *= 1. + gamma * sqrt( x ) + eta * x;
//   sivers_x *= 1. + 0.*gamma * sqrt( x ) * pow2(1. -x ) + eta * x * pow(1. -x, gamma );
//    sivers_x *= 1. + gamma * pow( x, eta ); // Boer 4 fit
//   sivers_x *= 1. + gamma * pow( 1. - x, eta ); // Boer 5 fit
//   sivers_x *= -(1. + gamma * ( 1. - x )/ x); // Boer 6 fit alla Bacchetta
//   sivers_x *= 1. + gamma * exp( eta * log( x ) ); // Boer 7 fit without antiquarks...

//   if( x >= 0.01 )
//     sivers_x = n/x; // Boer 8 ... :)
//   else
//     sivers_x = gamma/x;

//   double x_1 = 0.05, x_2 = 0.3;
//   double x_1 = 0.05, x_2 = 1.;
//   double x_1 = eta, x_2 = 0.95; // Boer 10 ... :)

//   if( x >= x_1 && x <= x_2)
//     sivers_x = n/x; // Boer 9 ... :)
//   else if( x < x_1 )
//     sivers_x = gamma/x;
//   else
//     sivers_x = 0.;

//   sivers_x = n; // Boer 9 ... :)




  return sivers_x;
}
