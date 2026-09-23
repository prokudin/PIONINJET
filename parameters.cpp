#ifndef __PARAMETERS_H_
#include "parameters.h"
#endif

#include <iostream>
#include <fstream>

using namespace std;

PARAMETERS::PARAMETERS() {
//   // init all parameters
//    kt2_average = 0.25; // MAIN
//    ptq2_average = 0.2; // MAIN
//   // init all parameters
//  kt2_average = 0.025; // TEST
//  ptq2_average = 0.02; // TEST


  kt2_average = 0.57; // MAIN HERMES 2016
  ptq2_average = 0.12; // MAIN HERMES 2016

//   kt2_average = 0.25; // MAIN HERMES 2016
//   ptq2_average = 0.1502; // MAIN HERMES 2016

//   kt2_average = 0.4; // MAIN HERMES 2016
//   ptq2_average = 0.1327; // MAIN HERMES 2016
//

//  kt2_average = 0.15; // MAIN HERMES 2016
//  ptq2_average = 0.167; //MAIN HERMES 2016

//   // init all parameters
//   kt2_average = 0.2602; // results of the fit for BM
//   ptq2_average = 0.1886;

   Sivers.SetReadFile("./parameters/sivers_parameters.dat"); // USE VALENCE_SIDIS=True if sivers_parametersMAIN_ZHONGBO_2013.dat is used and False if others
   Sivers.ReadParameters();
//
   Sivers.SetSetsReadFile("./parameters/sivers_sets.dat");
   Sivers.ReadSets();
//
//   Pretzelosity.SetReadFile("./parameters/pretzelosity_parameters.dat");
//   Pretzelosity.ReadParameters();
//
//   Pretzelosity.SetSetsReadFile("./parameters/pretzelosity_sets.dat");
//   Pretzelosity.ReadSets();
//
//   BoerMulders.SetReadFile("./parameters/boermulders_parameters.dat");
//   BoerMulders.ReadParameters();
//
//   BoerMulders.SetSetsReadFile("./parameters/boermulders_sets.dat");
//   BoerMulders.ReadSets();
//
//
//   CollinsKaon.SetReadFile("./parameters/collins_kaons_parameters.dat");
//   CollinsKaon.ReadParameters();
//
//   CollinsKaon.SetSetsReadFile("./parameters/collins_kaons_sets.dat");
//   CollinsKaon.ReadSets();


  Collins.SetReadFile("./parameters/collins_parameters.dat");
  Collins.ReadParameters();

  Collins.SetSetsReadFile("./parameters/collins_sets.dat");
  Collins.ReadSets();

  Transversity.SetReadFile("./parameters/transversity_parameters.dat");
  Transversity.ReadParameters();

  Transversity.SetSetsReadFile("./parameters/transversity_sets.dat");
  Transversity.ReadSets();


//   G1T.SetReadFile("./parameters/g1t_parameters.dat");
//   G1T.ReadParameters();
//
//   G2.SetReadFile("./parameters/g2_parameters.dat");
//   G2.ReadParameters();
//
//   G1.SetReadFile("./parameters/g1_parameters.dat");
//   G1.ReadParameters();
//
//   H1L.SetReadFile("./parameters/h1L_parameters.dat");
//   H1L.ReadParameters();


  UnpolarisedPDF.SetReadFile("./parameters/unpolarisedpdf_parameters.dat");
  UnpolarisedPDF.ReadParameters();

  UnpolarisedFF.SetReadFile("./parameters/unpolarisedff_parameters.dat");
  UnpolarisedFF.ReadParameters();

  UnpolarisedGK.SetReadFile("./parameters/unpolarisedgk_parameters.dat");
  UnpolarisedGK.ReadParameters();


};






//========================================================= parameters
void SETS::ReadParameters(){

   ifstream in;
   in.open(FILE_READ); // fit

   if(!in) {
     cerr << "PARAMETERS: Unable to open file with parameters! " << FILE_READ << endl;
   } else {
     cout << "PARAMETERS: Open " << FILE_READ << endl;
   }

   in >>  parameters.a_up >>  errors.a_up;
   in >>  parameters.b_up >>  errors.b_up;
   in >>  parameters.n_up >>  errors.n_up;
   in >>  parameters.gamma_up >>  errors.gamma_up;
   in >>  parameters.eta_up >>  errors.eta_up;
   in >>  parameters.m2_up>>  errors.m2_up;

   in >>  parameters.a_down     >>  errors.a_down;
   in >>  parameters.b_down     >>  errors.b_down;
   in >>  parameters.n_down     >>  errors.n_down;
   in >>  parameters.gamma_down >>  errors.gamma_down;
   in >>  parameters.eta_down >>  errors.eta_down;
   in >>  parameters.m2_down    >>  errors.m2_down;

   in >>  parameters.a_anti_up >>  errors.a_anti_up;
   in >>  parameters.b_anti_up >>  errors.b_anti_up;
   in >>  parameters.n_anti_up >>  errors.n_anti_up;
   in >>  parameters.gamma_anti_up >>  errors.gamma_anti_up;
   in >>  parameters.eta_anti_up >>  errors.eta_anti_up;
   in >>  parameters.m2_anti_up>>  errors.m2_anti_up;

   in >>  parameters.a_anti_down     >>  errors.a_anti_down;
   in >>  parameters.b_anti_down     >>  errors.b_anti_down;
   in >>  parameters.n_anti_down     >>  errors.n_anti_down;
   in >>  parameters.gamma_anti_down >>  errors.gamma_anti_down;
   in >>  parameters.eta_anti_down >>  errors.eta_anti_down;
   in >>  parameters.m2_anti_down    >>  errors.m2_anti_down;


   in >>  parameters.a_strange >>  errors.a_strange;
   in >>  parameters.b_strange >>  errors.b_strange;
   in >>  parameters.n_strange >>  errors.n_strange;
   in >>  parameters.gamma_strange >>  errors.gamma_strange;
   in >>  parameters.eta_strange >>  errors.eta_strange;
   in >>  parameters.m2_strange>>  errors.m2_strange;

   in >>  parameters.a_anti_strange >>  errors.a_anti_strange;
   in >>  parameters.b_anti_strange >>  errors.b_anti_strange;
   in >>  parameters.n_anti_strange >>  errors.n_anti_strange;
   in >>  parameters.gamma_anti_strange >>  errors.gamma_anti_strange;
   in >>  parameters.eta_anti_strange >>  errors.eta_anti_strange;
   in >>  parameters.m2_anti_strange>>  errors.m2_anti_strange;


   in.close();

   main_parameters = parameters; // remember the main set of parameters
};


//========================================================= parameters OLD results Sivers
void SETS::ReadParametersOld(){

   ifstream in;
   in.open(FILE_READ_OLD); // fit

   if(!in) {
     cerr << "PARAMETERS: Unable to open file with parameters! " << FILE_READ_OLD << endl;
   } else {
     cout << "PARAMETERS: Open " << FILE_READ_OLD << endl;
   }


   double a_up, a_up_error,
     b_up, b_up_error,
     n_up, n_up_error,
     a_down, a_down_error,
     b_down, b_down_error,
     n_down, n_down_error,
     a_up_sea, a_up_sea_error,
     b_up_sea, b_up_sea_error,
     n_up_sea, n_up_sea_error,
     a_down_sea, a_down_sea_error,
     b_down_sea, b_down_sea_error,
     n_down_sea, n_down_sea_error,
     m2, m2_error, kt2_averageO, kt2_average_error,
     ptq2_averageO,  ptq2_average_error;

   in >> a_up >> a_up_error;
   in >> b_up >> b_up_error;
   in >> n_up >> n_up_error;

   in >> a_down >> a_down_error;
   in >> b_down >> b_down_error;
   in >> n_down >> n_down_error;

   in >> a_up_sea >> a_up_sea_error;
   in >> b_up_sea >> b_up_sea_error;
   in >> n_up_sea >> n_up_sea_error;

   in >> a_down_sea >> a_down_sea_error;
   in >> b_down_sea >> b_down_sea_error;
   in >> n_down_sea >> n_down_sea_error;

   in >> m2 >> m2_error;

   in >> kt2_averageO >> kt2_average_error;
   in >> ptq2_averageO>> ptq2_average_error;

   //// Parametrization FIT08 or FIT06 Sivers aka Stefano......

   parameters.a_up = a_up;  errors.a_up = a_up_error;
   parameters.b_up = b_up;  errors.b_up = b_up_error;
   parameters.n_up = n_up;  errors.n_up = n_up_error;
   parameters.m2_up= m2;    errors.m2_up= m2_error;

   parameters.a_down  = a_down;   errors.a_down  = a_down_error;
   parameters.b_down  = b_up;     errors.b_down  = b_up_error;
   parameters.n_down  = n_down;   errors.n_down  = n_down_error;
   parameters.m2_down = m2;       errors.m2_down = m2_error;

   parameters.a_anti_up = a_up_sea;  errors.a_anti_up = a_up_sea_error;
   parameters.b_anti_up = b_up;      errors.b_anti_up = b_up_error;
   parameters.n_anti_up = n_up_sea;  errors.n_anti_up = n_up_sea_error;
   parameters.m2_anti_up= m2;        errors.m2_anti_up= m2_error;

   parameters.a_anti_down  = a_up_sea;    errors.a_anti_down = a_up_sea_error;
   parameters.b_anti_down  = b_up;        errors.b_anti_down = b_up_error;
   parameters.n_anti_down  =  n_down_sea; errors.n_anti_down = n_down_sea_error;
   parameters.m2_anti_down = m2;          errors.m2_anti_down= m2_error;


   parameters.a_strange  = a_up_sea;    errors.a_strange   = a_up_sea_error;
   parameters.b_strange  = b_up;        errors.b_strange   = b_up_error;
   parameters.n_strange  = a_down_sea;  errors.n_strange = a_down_sea_error;
   parameters.m2_strange = m2;          errors.m2_strange= m2_error;

   parameters.a_anti_strange  = a_up_sea;   errors.a_anti_strange = a_up_sea_error;
   parameters.b_anti_strange  = b_up;       errors.b_anti_strange = b_up_error;
   parameters.n_anti_strange  = b_down_sea; errors.n_anti_strange = b_down_sea_error;
   parameters.m2_anti_strange = m2;         errors.m2_anti_strange= m2_error;


   in.close();

   main_parameters = parameters; // remember the main set of parameters

};



//========================================================= parameters
void SETS::UseSet(int k){
   parameters =  set[k];
}

//========================================================= parameters
void SETS::UseMainSet(void){
   parameters = main_parameters;
}




//========================================================= Constructor
void SETS::SetReadFile(const char* file){
  FILE_READ = file;
};

//========================================================= Constructor
void SETS::SetSetsReadFile(const char* file){
  SETS_FILE_READ = file;
};

//========================================================= Constructor
void SETS::SetReadFileOld(const char* file){
  FILE_READ_OLD = file;
};

//========================================================= Constructor
void SETS::SetSetsReadFileOld(const char* file){
  SETS_FILE_READ_OLD = file;
};

//========================================================= Constructor
void SETS::SetWriteFile(const char* file){
  FILE_WRITE = file;
};

//========================================================= Constructor
void SETS::SetSetsWriteFile(const char* file){
  SETS_FILE_WRITE = file;
};

//========================================================= parameters
void SETS::ReadSets(){

   ifstream in;
   in.open(SETS_FILE_READ); // fit

   if(!in) {
     cerr << "PARAMETERS: SETS: Unable to open file with parameters! " << SETS_FILE_READ << endl;
   } else {
     cout << "PARAMETERS: SETS: Open " << SETS_FILE_READ << endl;
   }

   for(int i = 0; i < set_number; i++){
   in >>  set[i].a_up;
   in >>  set[i].b_up;
   in >>  set[i].n_up;
   in >>  set[i].gamma_up;
   in >>  set[i].eta_up;
   in >>  set[i].m2_up;

   in >>  set[i].a_down;
   in >>  set[i].b_down;
   in >>  set[i].n_down;
   in >>  set[i].gamma_down;
   in >>  set[i].eta_down;
   in >>  set[i].m2_down;

   in >>  set[i].a_anti_up;
   in >>  set[i].b_anti_up;
   in >>  set[i].n_anti_up;
   in >>  set[i].gamma_anti_up;
   in >>  set[i].eta_anti_up;
   in >>  set[i].m2_anti_up;

   in >>  set[i].a_anti_down;
   in >>  set[i].b_anti_down;
   in >>  set[i].n_anti_down;
   in >>  set[i].gamma_anti_down;
   in >>  set[i].eta_anti_down;
   in >>  set[i].m2_anti_down;


   in >>  set[i].a_strange;
   in >>  set[i].b_strange;
   in >>  set[i].n_strange;
   in >>  set[i].gamma_strange;
   in >>  set[i].eta_strange;
   in >>  set[i].m2_strange;

   in >>  set[i].a_anti_strange;
   in >>  set[i].b_anti_strange;
   in >>  set[i].n_anti_strange;
   in >>  set[i].gamma_anti_strange;
   in >>  set[i].eta_anti_strange;
   in >>  set[i].m2_anti_strange;
   }

   in.close();
};

//========================================================= parameters < 09
void SETS::ReadSetsOld(){

   ifstream in;
   in.open(SETS_FILE_READ_OLD); // fit

   if(!in) {
     cerr << "PARAMETERS: SETS: Unable to open file with parameters! " << SETS_FILE_READ_OLD << endl;
   } else {
     cout << "PARAMETERS: SETS: Open " << SETS_FILE_READ_OLD << endl;
   }

   double a_up,
     b_up,
     n_up,
     a_down,
     b_down,
     n_down,
     a_up_sea,
     b_up_sea,
     n_up_sea,
     a_down_sea,
     b_down_sea,
     n_down_sea,
     m2;



   for(int i = 0; i < set_number; i++){

   in >> a_up;
   in >> b_up;
   in >> n_up;

   in >> a_down;
   in >> b_down;
   in >> n_down;

   in >> a_up_sea;
   in >> b_up_sea;
   in >> n_up_sea;

   in >> a_down_sea;
   in >> b_down_sea;
   in >> n_down_sea;

   in >> m2;

   //// Parametrization FIT08 or FIT06 Sivers aka Stefano......

   set[i].a_up = a_up;
   set[i].b_up = b_up;
   set[i].n_up = n_up;
   set[i].m2_up= m2;

   set[i].a_down  = a_down;
   set[i].b_down  = b_up;
   set[i].n_down  = n_down;
   set[i].m2_down = m2;

   set[i].a_anti_up = a_up_sea;
   set[i].b_anti_up = b_up;
   set[i].n_anti_up = n_up_sea;
   set[i].m2_anti_up= m2;

   set[i].a_anti_down  = a_up_sea;
   set[i].b_anti_down  = b_up;
   set[i].n_anti_down  =  n_down_sea;
   set[i].m2_anti_down = m2;


   set[i].a_strange  = a_up_sea;
   set[i].b_strange  = b_up;
   set[i].n_strange  = a_down_sea;
   set[i].m2_strange = m2;

   set[i].a_anti_strange  = a_up_sea;
   set[i].b_anti_strange  = b_up;
   set[i].n_anti_strange  = b_down_sea;
   set[i].m2_anti_strange = m2;
   }

   in.close();
};



//========================================================= parameters
void SETS::WriteParameters(){

//    ofstream out;
//    out.open(FILE_WRITE); // fit

//    if(!out) {
//      cerr << "PARAMETERS: Unable to open file to write parameters! " << FILE_WRITE << endl;
//    } else {
//      cout << "PARAMETERS: Open " << FILE_WRITE << endl;
//    }

   FILE *out = fopen( FILE_WRITE, "w" );

   if( out == NULL ){
     cerr << "PARAMETERS: Unable to open file to write parameters! " << FILE_WRITE << endl;
   } else {
     cout << "PARAMETERS: Open " << FILE_WRITE << endl;
   }


   fprintf(out,"%8.6e  %8.6e\n",parameters.a_up, errors.a_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.b_up, errors.b_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.n_up, errors.n_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.gamma_up, errors.gamma_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.eta_up, errors.eta_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.m2_up, errors.m2_up);

   fprintf(out,"%8.6e  %8.6e\n",parameters.a_down, errors.a_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.b_down, errors.b_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.n_down, errors.n_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.gamma_down, errors.gamma_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.eta_down, errors.eta_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.m2_down, errors.m2_down);

   fprintf(out,"%8.6e  %8.6e\n",parameters.a_anti_up, errors.a_anti_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.b_anti_up, errors.b_anti_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.n_anti_up, errors.n_anti_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.gamma_anti_up, errors.gamma_anti_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.eta_anti_up, errors.eta_anti_up);
   fprintf(out,"%8.6e  %8.6e\n",parameters.m2_anti_up, errors.m2_anti_up);

   fprintf(out,"%8.6e  %8.6e\n",parameters.a_anti_down, errors.a_anti_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.b_anti_down, errors.b_anti_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.n_anti_down, errors.n_anti_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.gamma_anti_down, errors.gamma_anti_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.eta_anti_down, errors.eta_anti_down);
   fprintf(out,"%8.6e  %8.6e\n",parameters.m2_anti_down, errors.m2_anti_down);

   fprintf(out,"%8.6e  %8.6e\n",parameters.a_strange, errors.a_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.b_strange, errors.b_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.n_strange, errors.n_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.gamma_strange, errors.gamma_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.eta_strange, errors.eta_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.m2_strange, errors.m2_strange);


   fprintf(out,"%8.6e  %8.6e\n",parameters.a_anti_strange, errors.a_anti_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.b_anti_strange, errors.b_anti_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.n_anti_strange, errors.n_anti_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.gamma_anti_strange, errors.gamma_anti_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.eta_anti_strange, errors.eta_anti_strange);
   fprintf(out,"%8.6e  %8.6e\n",parameters.m2_anti_strange, errors.m2_anti_strange);


//    out.close();

};

//========================================================= parameters
void SETS::WriteSets(){

//    ofstream out;
//    out.open(SETS_FILE_WRITE); // fit

//    if(!out) {
//      cerr << "PARAMETERS: SETS: Unable to open file to write parameters! " << SETS_FILE_WRITE << endl;
//    } else {
//      cout << "PARAMETERS: SETS: Open " << SETS_FILE_WRITE << endl;
//    }



   FILE *out = fopen( SETS_FILE_WRITE, "w" );

   if( out == NULL ){
     cerr << "PARAMETERS: SETS: Unable to open file to write parameters! " << SETS_FILE_WRITE << endl;
   } else {
     cout << "PARAMETERS: SETS: Open " << SETS_FILE_WRITE << endl;
   }




   for(int i = 0; i < set_number; i++){

   fprintf(out,"%8.6e  ",set[i].a_up);
   fprintf(out,"%8.6e  ",set[i].b_up);
   fprintf(out,"%8.6e  ",set[i].n_up);
   fprintf(out,"%8.6e  ",set[i].gamma_up);
   fprintf(out,"%8.6e  ",set[i].eta_up);
   fprintf(out,"%8.6e  ",set[i].m2_up);

   fprintf(out,"%8.6e  ",set[i].a_down);
   fprintf(out,"%8.6e  ",set[i].b_down);
   fprintf(out,"%8.6e  ",set[i].n_down);
   fprintf(out,"%8.6e  ",set[i].gamma_down);
   fprintf(out,"%8.6e  ",set[i].eta_down);
   fprintf(out,"%8.6e  ",set[i].m2_down);

   fprintf(out,"%8.6e  ",set[i].a_anti_up);
   fprintf(out,"%8.6e  ",set[i].b_anti_up);
   fprintf(out,"%8.6e  ",set[i].n_anti_up);
   fprintf(out,"%8.6e  ",set[i].gamma_anti_up);
   fprintf(out,"%8.6e  ",set[i].eta_anti_up);
   fprintf(out,"%8.6e  ",set[i].m2_anti_up);

   fprintf(out,"%8.6e  ",set[i].a_anti_down);
   fprintf(out,"%8.6e  ",set[i].b_anti_down);
   fprintf(out,"%8.6e  ",set[i].n_anti_down);
   fprintf(out,"%8.6e  ",set[i].gamma_anti_down);
   fprintf(out,"%8.6e  ",set[i].eta_anti_down);
   fprintf(out,"%8.6e  ",set[i].m2_anti_down);

   fprintf(out,"%8.6e  ",set[i].a_strange );
   fprintf(out,"%8.6e  ",set[i].b_strange );
   fprintf(out,"%8.6e  ",set[i].n_strange );
   fprintf(out,"%8.6e  ",set[i].gamma_strange);
   fprintf(out,"%8.6e  ",set[i].eta_strange);
   fprintf(out,"%8.6e  ",set[i].m2_strange);


   fprintf(out,"%8.6e  ",set[i].a_anti_strange);
   fprintf(out,"%8.6e  ",set[i].b_anti_strange);
   fprintf(out,"%8.6e  ",set[i].n_anti_strange);
   fprintf(out,"%8.6e  ",set[i].gamma_anti_strange);
   fprintf(out,"%8.6e  ",set[i].eta_anti_strange);
   fprintf(out,"%8.6e  \n",set[i].m2_anti_strange);

   }

//    out.close();

};


// Evaluate min and max of a function f(x) on the parameters sets
std::vector<double> SETS::MinMax( double (*f)(double ), double x ){
// std::vector<double> SETS::MinMax( void (*f)(const int ncomputations, const int ndim, double x, integrand_t integrand, double asymmetry[]), double x ){
  vector<double> r;
  double min, max;

//   for(int k = 0; k < set_number; k++){
  for(int k = 0; k < 2; k++){

    SETS::UseSet(k); // we do loop over all sets...

    double sechenie = (*f)( x );

    cout << k << endl;

    if( k == 0 ){
      min = sechenie;
      max = sechenie;
    } else {
      if( (sechenie = (*f)( x ) ) > max ){
	max = sechenie;
      }
      if(  sechenie < min ){
	min = sechenie;
      }
    }
  }

  SETS::UseMainSet(); // main set is restored after the loop ...

  r.push_back(min);
  r.push_back(max);

  return r;

};

void SETS::Print(){ //Just prints values of parameters
  cout << "U quark " << endl;
  cout <<  parameters.a_up << " " << errors.a_up << endl;
  cout <<  parameters.b_up << " " <<  errors.b_up << endl;
  cout <<  parameters.n_up << " " <<  errors.n_up << endl;
  cout <<  parameters.gamma_up << " " <<  errors.gamma_up << endl;
  cout <<  parameters.eta_up << " " <<  errors.eta_up << endl;
  cout <<  parameters.m2_up<< " " <<  errors.m2_up<< endl;

  cout << "D quark " << endl;
  cout <<  parameters.a_down     << " " <<  errors.a_down << endl;
  cout <<  parameters.b_down     << " " <<  errors.b_down << endl;
  cout <<  parameters.n_down     << " " <<  errors.n_down << endl;
  cout <<  parameters.gamma_down << " " <<  errors.gamma_down << endl;
  cout <<  parameters.eta_down << " " <<  errors.eta_down << endl;
  cout <<  parameters.m2_down    << " " <<  errors.m2_down << endl;

  cout << "UBAR quark " << endl;
  cout <<  parameters.a_anti_up << " " <<  errors.a_anti_up << endl;
  cout <<  parameters.b_anti_up << " " <<  errors.b_anti_up << endl;
  cout <<  parameters.n_anti_up << " " <<  errors.n_anti_up << endl;
  cout <<  parameters.gamma_anti_up << " " <<  errors.gamma_anti_up << endl;
  cout <<  parameters.eta_anti_up << " " <<  errors.eta_anti_up << endl;
  cout <<  parameters.m2_anti_up<< " " <<  errors.m2_anti_up << endl;

  cout << "D bar quark " << endl;
  cout <<  parameters.a_anti_down     << " " <<  errors.a_anti_down << endl;
  cout <<  parameters.b_anti_down     << " " <<  errors.b_anti_down << endl;
  cout <<  parameters.n_anti_down     << " " <<  errors.n_anti_down << endl;
  cout <<  parameters.gamma_anti_down << " " <<  errors.gamma_anti_down << endl;
  cout <<  parameters.eta_anti_down << " " <<  errors.eta_anti_down << endl;
  cout <<  parameters.m2_anti_down    << " " <<  errors.m2_anti_down << endl;


  cout << "S quark " << endl;
  cout <<  parameters.a_strange << " " <<  errors.a_strange << endl;
  cout <<  parameters.b_strange << " " <<  errors.b_strange << endl;
  cout <<  parameters.n_strange << " " <<  errors.n_strange << endl;
  cout <<  parameters.gamma_strange << " " <<  errors.gamma_strange << endl;
  cout <<  parameters.eta_strange << " " <<  errors.eta_strange << endl;
  cout <<  parameters.m2_strange<< " " <<  errors.m2_strange << endl;

  cout << "S BAR quark " << endl;
  cout <<  parameters.a_anti_strange << " " <<  errors.a_anti_strange << endl;
  cout <<  parameters.b_anti_strange << " " <<  errors.b_anti_strange << endl;
  cout <<  parameters.n_anti_strange << " " <<  errors.n_anti_strange << endl;
  cout <<  parameters.gamma_anti_strange << " " <<  errors.gamma_anti_strange << endl;
  cout <<  parameters.eta_anti_strange << " " <<  errors.eta_anti_strange << endl;
  cout <<  parameters.m2_anti_strange<< " " <<  errors.m2_anti_strange << endl;

};
