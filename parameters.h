#ifndef  __PARAMETERS_H_
#define  __PARAMETERS_H_

#include <vector>
#include <cmath>

#include <common_sidis.h>


using namespace std;
//using namespace dy_namespace;
using namespace sidis_namespace;


// Parameters class
class VALUES {
 private:
 public:
  VALUES() {};
  ~VALUES() {};

  // Down quark
  double a_down;
  double b_down;
  double n_down;
  double gamma_down;
  double eta_down;
  double m2_down; // kt dependence


  // Up quark
  double a_up;
  double b_up;
  double n_up;
  double gamma_up;
  double eta_up;
  double m2_up; // kt dependence

  // anti Down quark
  double a_anti_down;
  double b_anti_down;
  double n_anti_down;
  double gamma_anti_down;
  double eta_anti_down;
  double m2_anti_down; // kt dependence


  // anti Up quark
  double a_anti_up;
  double b_anti_up;
  double n_anti_up;
  double gamma_anti_up;
  double eta_anti_up;
  double m2_anti_up; // kt dependence

  // Strange quark
  double a_strange;
  double b_strange;
  double n_strange;
  double gamma_strange;
  double eta_strange;
  double m2_strange;


  // anti Strange quark
  double a_anti_strange;
  double b_anti_strange;
  double n_anti_strange;
  double gamma_anti_strange;
  double eta_anti_strange;
  double m2_anti_strange;
};




const int set_number = 1500; // we will generate 1500 sets
//const int set_number = 200; // we will generate 200 sets MAIN
//const int set_number = 290; // we will generate 200 sets MAIN TEST SCAN ANDY!!!!
//const int set_number = 500; // we will generate 500 sets MAIN
//const int set_number = 10000; // number of sets for sivers coridor TEST FIT5
//const int set_number = 200; // number of sets for sivers coridor COLLINS
//const int set_number = 100; // number of sets for sivers coridor

//const int set_number = 86; // we will generate 86 sets MY PENG 2015
//const int set_number = 135; // we will generate 135 setsPENG 2015

//const int set_number = 290; // we will generate 290 sets Sivers KANG, GAMBERG 2013

//const int set_number = 10000; // we will generate 10000 sets test 2016

//const int set_number = 20000; // we will generate 20000 sets test1 2016

//const int set_number = 1848; // we will generate STEFANO 2015 TRANSVERSITY ANSELMINO 13


//const int set_number = 100; // we will generate TEST

class SETS{ // organize all sets here params, errors, sets of params
 private:
  VALUES main_parameters;
  const char* FILE_READ; // file to read parameters
  const char* SETS_FILE_READ; // file to read set parameters
  const char* FILE_WRITE; // file to write parameters
  const char* SETS_FILE_WRITE; // file to read set parameters
  const char* FILE_READ_OLD; // file to read parameters OLD FITS < 09
  const char* SETS_FILE_READ_OLD; // file to read set parameters OLD FITS < 09

 public:
  SETS() {}; // constructor...
  ~SETS() {}; // deconstructor...

  VALUES parameters;
  VALUES errors;

  VALUES set[set_number];


//  const char * deals with warning: deprecated conversion from string constant to ?char*?
//  http://en.wikibooks.org/wiki/GCC_Debugging/g%2B%2B/Warnings/deprecated_conversion_from_string_constant


  void SetReadFile(const char* file);
  void SetSetsReadFile(const char* file);
  void SetReadFileOld(const char* file);
  void SetSetsReadFileOld(const char* file);
  void SetWriteFile(const char* file);
  void SetSetsWriteFile(const char* file);
  void ReadParameters();
  void ReadSets();
  void ReadParametersOld();
  void ReadSetsOld();
  void WriteParameters();
  void WriteSets();
  void UseMainSet(void);
  void UseSet(int k);
  void Print();

  vector<double> MinMax( double (*f)(double x), double x ); // evaluate min and max of a function f(x) on the parameters sets


};


// HERE WE DEFINE ALL PARAMETERS...
class PARAMETERS{
 private:
  // All params...
  double kt2_average;
  double ptq2_average;


 public:
  PARAMETERS(); // constructor...
  //Get params...
  double GetKt2Average( void ) const { return kt2_average;};
  double GetPtq2Average( void ) const { return ptq2_average;};

/*   double GetKt2Average( double x ) const { return 0.25; }; // STANDARD */
/*   double GetPtq2Average( double z ) const { return 0.2; }; // STANDARD */


  //double GetKt2Average( double x ) const { return 0.25; }; // COMPASS...
/*   double GetKt2Average( double x ) const { return 0.18; }; // HERMES... */
/*   double GetPtq2Average( double z ) const { return 0.42*pow(1.-z,0.54)*pow(z,0.37);}; // HERMES... */


  // MY TESTS:
/*   double GetKt2Average( double x ) const { return 0.1*pow(1.-x,0.5)*pow(x,-0.3); }; // MY TESTS */
/*   double GetPtq2Average( double z ) const { return 0.42*pow(1.-z,0.54)*pow(z,0.37);}; // MY TES */

/*   double GetKt2Average( double x ) const { return 0.08*pow(1.-x,-0.2)*pow(x,-0.3); }; // MY TESTS */
/*   double GetKt2Average( double x ) const { return 0.25 - 0.14*pow(1.-x,0.1)*pow(x,0.3); }; // MY TESTS Goog!!! */
/*   double GetKt2Average( double x ) const { return 0.3 - 1.43*pow(1.-x,2.4)*pow(x,0.96); }; // MY TESTS */
  double GetKt2Average( double x ) const { return 0.3 - 0.25*pow(1.-x,0.1)*pow(x,0.3); }; // MY TESTS BEST GUESS !!!
  double GetPtq2Average( double z ) const { return 0.42*pow(1.-z,0.54)*pow(z,0.37);}; // MY TEST BEST GUESS !!!

  void SetKt2Average( double x ) { if( x >= 0.) kt2_average = x;};
  void SetPtq2Average( double x ) { if( x >= 0.) ptq2_average = x;};


  double GetKt2Average( int EXPERIMENT, double x ) {
    if( sidis_namespace::EXPERIMENT == HERMES || sidis_namespace::EXPERIMENT == HERMES_UNP)
      //return  GetKt2Average( x ); //
      return 0.18;
    else if( sidis_namespace::EXPERIMENT == COMPASS_SIDIS || sidis_namespace::EXPERIMENT == COMPASS_SIDIS_UNP)
      //return GetKt2Average( x ); // 0.25;
      return 0.25;
    else
      return GetKt2Average();
  };  // EXPERIMENT DEPENDENT KT...

  double GetPtq2Average( int EXPERIMENT, double z ) {
    if( EXPERIMENT == sidis_namespace::HERMES || EXPERIMENT == sidis_namespace::HERMES_UNP)
      //return GetPtq2Average( z ); // return 0.42*pow(1.-z,0.54)*pow(z,0.37);
    return 0.2;
    else if( EXPERIMENT == sidis_namespace::COMPASS_SIDIS || EXPERIMENT == sidis_namespace::COMPASS_SIDIS_UNP)
      //return GetPtq2Average( z ); // return 0.2;
    return 0.2;
    else
      return GetPtq2Average();
  };  // EXPERIMENT DEPENDENT KT...



  SETS Sivers;
  SETS Collins;
  SETS CollinsKaon;
  SETS Transversity;
  SETS BoerMulders;
  SETS Pretzelosity;
  SETS G1T;
  SETS G2;
  SETS H1L;
  SETS G1;
// UNpolarised PDF and Unpolarised FF
  SETS UnpolarisedPDF;
  SETS UnpolarisedFF;
  SETS UnpolarisedGK;
};

#endif // #ifndef  __PARAMETERS_H_
