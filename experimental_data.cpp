#include <iostream>
#include "experimental_data.h"
#include <fstream>

//--- This is for root use in C++
#include "TApplication.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLatex.h"
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
#include "TROOT.h"
#include "TFile.h" 
#include "TObject.h"


using namespace std;


EXPERIMENTAL_DATA::EXPERIMENTAL_DATA() : max_number(30), max_length(1000), ntotal(0){
  // allocate memory for data
  data = new double * [max_number];

  for(int i = 0; i < max_number; i++)
    data[i] = new double[ max_length ];
} //constructor


// EXPERIMENTAL_DATA::EXPERIMENTAL_DATA(int ntotal, int nraws) : nraws(nraws), ntotal(ntotal){
// 
// 
//   data = new double * [ntotal];
// 
//   for(int i = 0; i < ntotal; i++)
//     data[i] = new double[ nraws ];
// 
// //     if( i <= ntotal && j <= nraws){ 
// //       return data[i][j];
// //     }
// //   // allocate memory for data
// //   data = new double * [max_number];
// // 
// //   for(int i = 0; i < max_number; i++)
// //     data[i] = new double[ max_length ];
// } //constructor

//constructor read file, skip lines, obtain ndata etc...
EXPERIMENTAL_DATA::EXPERIMENTAL_DATA(char const * file) : max_number(30), max_length(1000), ntotal(0){
  filename = file;
  process_file( file );
}; 

//constructor read file, x -> data[i1], y-> data[i2], dy -> i3
EXPERIMENTAL_DATA::EXPERIMENTAL_DATA(char const * file, int i1, int i2, int i3) : max_number(30), max_length(1000){
  filename = file;
  process_file( file );
  cout << "New constructor" << endl;

  X(i1);  // assign x (remember that in c++ we count from 0 not 1 => i1 -1)
  Y(i2);  // assign y
  dY(i3); // assign dy

  // dX is zero
  dx = new double[ntotal];
  for (int i = 0; i < ntotal; i++)
    dx[i] = 0.;

  //  cout << x[1] << endl;
}; 

EXPERIMENTAL_DATA::~EXPERIMENTAL_DATA() {
/*  delete [] data;
  delete [] x;
  delete [] y;
  delete [] dx;
  delete [] dy;*/
} //deconstructor


void EXPERIMENTAL_DATA::write_file(char const * file){
  ofstream out;
  out.open(file, ios::out);
  cout << "Open file " << file << endl;

  if(!out){
    cerr << "File " << file << " cannot be opened.";
  }; 
  
  for (int i = 0; i < ntotal; i++){
    for (int j = 0; j < nraws; j++){
      out << GetData( i, j) << " " ;
    }
    out << x[i]  << " " << y[i]  << " " << dy[i]  << "\n";
  }
  out.close();
};


void EXPERIMENTAL_DATA::process_file(char const * file){
  ifstream in;
  char text[256], comment[] = "#";
  in.open(file, ios::in);
  cout << "Open file " << file << endl;

  if(!in){
    cerr << "File " << file << " cannot be opened.";
  }; 
  
  long int position = 0;

  // read comments
  do {
    in.getline(text,256);
    //cout << text << endl;

    position += strlen(text) + 1;

  } while( strncmp (text,comment,1) == 0);

  position -= strlen(text) + 1;

  // set file position pointer to the beggining of the line
  in.seekg(position, ios::beg);


  // read data from the file, determine how many raws and columns we have...
  int counter_raws = 0;
  int counter_lines = 0;
  double temp_data;

  do {
    in >> temp_data;
    //cout << temp_data << " ";
    while(in.peek() == ' ') in.seekg(1, ios::cur);
    if( counter_lines == 0 ) // we count data in a first line
      counter_raws++;
  } while (in.peek() != '\n' );

  //cout << endl;

  //cout << "Number of raws: " << counter_raws << endl;

  in.seekg(position, ios::beg);

  // now I read n = counter_raws numbers from file
  while (1) {
    double temp_data[counter_raws];
    for( int i = 0; i < counter_raws; i ++){
      in >> temp_data[i];
      vdata.push_back(temp_data[i]);
      //cout << temp_data[i] << " ";
    }
    //cout << endl;

    if (in.eof()) break;

    counter_lines++;

    
    if (in.fail()) { 
      cout << "Error in reading" << file << endl; 
    }
  }

  //cout << endl << "Number of lines: " << counter_lines << endl;

  // now  create table with the data:
  data = new double * [counter_lines];

  for(int i = 0; i < counter_lines; i++)
    data[i] = new double[ counter_raws ];

  // fill the data
  for(int i = 0; i < counter_lines; i++){
    for(int j = 0; j < counter_raws; j++){

      data[i][j] = vdata[ i*counter_raws + j ];
      //cout << data[i][j] << " " ;
    }
    //cout << endl;
  }

  ntotal = counter_lines; // number of experimental points
  nraws  = counter_raws; // number of raws in data file
  //system("ls");

  in.close();
}; // read file, skip lines, obtain ndata etc...

// assign x array
void EXPERIMENTAL_DATA::X(int i) {
  i = i - 1; 

//   x = new double[sizeof( data[i] )/sizeof( double )];
  x = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    x[j] = data[j][i];
    //cout << "x[" << j << "] = " << x[j] << endl;
  }
 
}

// assign x array
void EXPERIMENTAL_DATA::X(int i, int j) {
  i = i - 1; 
  j = j - 1;
 
  x = new double[ntotal];

  for(int k = 0 ; k < ntotal; k++){
    x[k] = ( data[k][i] + data[k][j] ) / 2.;
    //cout << "x[" << k << "] = " << x[k] << endl;
  }
}

// assign y array
void EXPERIMENTAL_DATA::Y(int i) {
  i = i - 1; 

  y = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    y[j] = data[j][i];
    //cout << "y[" << j << "] = " << y[j] << endl;
  }

}

// assign dx array
void EXPERIMENTAL_DATA::dX(int i) {
  i = i - 1; 

  //  dx = new double[sizeof( data[i] )/sizeof( double )];

  dx = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    dx[j] = data[j][i];
    //cout << "dx[" << j << "] = " << dx[j] << endl;
  }

}

// assign dy array
void EXPERIMENTAL_DATA::dY(int i) {
  i = i - 1;

  dy = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    dy[j] = data[j][i];
    //cout << "dy[" << j << "] = " << dy[j] << endl;
  }
 
}

// assign dy array
void EXPERIMENTAL_DATA::dY(int i1, int i2) {
  i1 = i1 - 1;
  i2 = i2 - 1;

  dy = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    dy[j] = sqrt( data[j][i1] * data[j][i1]  + data[j][i2] * data[j][i2] );
    //cout << "dx[" << j << "] = " << dy[j] << endl;
  }
}

// assign dy array
void EXPERIMENTAL_DATA::dY(int i1, int i2, int i3) {
  i1 = i1 - 1;
  i2 = i2 - 1;
  i3 = i3 - 1;

  dy = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    dy[j] = sqrt( data[j][i1] * data[j][i1]  + data[j][i2] * data[j][i2]  + data[j][i3] * data[j][i3]);
    //cout << "dx[" << j << "] = " << dy[j] << endl;
  }
}

// assign dy array COMPASS 0.4 sigma ecc
void EXPERIMENTAL_DATA::dY(int i1, int i2, double coeff) {
  i1 = i1 - 1;
  i2 = i2 - 1;

  dy = new double[ntotal];

  for(int j = 0 ; j < ntotal; j++){
    dy[j] = sqrt( data[j][i1] * data[j][i1]  + coeff * coeff * data[j][i2] * data[j][i2] );
    //cout << "dx[" << j << "] = " << dy[j] << endl;
  }
}

TGraphErrors* EXPERIMENTAL_DATA::plot(){

  TGraphErrors* graph = plot(1,8);

  return graph;
}

TGraphErrors* EXPERIMENTAL_DATA::plot(int color, int marker){
  TGraphErrors* graph = new TGraphErrors(1); // in case data does not exist
  if( ntotal != 0 ){ // data exist
      double x_error[ntotal];
      for(int j = 0 ; j < ntotal; j++){
	x_error[j] = 0.;
	//cout << x[j] << " " << y[j] << " " << dy[j] <<endl;
      }

      graph = new TGraphErrors(ntotal, x, y, x_error, dy);
            
      graph->SetMarkerColor( color );
      graph->SetMarkerStyle( marker); 
      graph->SetMarkerSize(1.);
      graph->SetLineWidth(2.);
      graph->SetLineColor(color);
      graph->SetTitle("");
  }
  return graph;
}
