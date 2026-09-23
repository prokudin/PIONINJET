#ifndef __EXPERIMENTAL_DATA_
#define __EXPERIMENTAL_DATA_

#include <vector>
//#include <stdlib>
#include <iostream>
#include <string>

#include "TGraphErrors.h"

using namespace std;

// class experimental data
class EXPERIMENTAL_DATA {
 private:
  const int max_number; // lines in a file
  const int max_length; // number of exp data
  int ntotal; // number of data
  int nraws; // number of raws in data files
  double **data; // arrays with data 
  char const * filename; // filename to read the data from
  double *x, *y, *dx, *dy; // pointers to x, y, dx, dy  arrays
  vector<double> vdata;
  void process_file(char const * file);
  void write_file(char const * file);
 public:
  EXPERIMENTAL_DATA(); //constructor
/*   EXPERIMENTAL_DATA(char * file, double X( int i), double Y( int i1), double dY( int i2, int i3 )); //constructor read file, skip lines, obtain ndata etc... */
  EXPERIMENTAL_DATA(char const * file); //constructor read file, skip lines, obtain ndata etc...
  EXPERIMENTAL_DATA(char const * file, int i1, int i2, int i3); //constructor read file, x -> data[i1], y-> data[i2], dy -> i3

//  EXPERIMENTAL_DATA(int ntotal, int nraws); //constructor  of empty data with ntotal&nraws...

    ~EXPERIMENTAL_DATA(); //deconstructor
    // copy constructor...

  void X( int i );
  void X( int i , int j);
  void Y( int i );
  void dX(int i ); 
  void dY( int i1); // dy  
  void dY( int i1, int i2); // dy 
  void dY( int i1, int i2, int i3); // dy 
  void dY( int i1, int i2, double coeff); // dy for 2009 COMPASS data, they use 0.4 sigma for systematics.. 


  void read(char const * file){
    process_file(file);
  };   


  void write(char const * file){
    write_file(file);
  };    

  void print(){
     for (int i = 0; i < Ntotal(); i++){
     for (int j = 0; j < Nraws(); j++){
      cout << GetData( i, j) << " ";
    }
   cout << endl;
  }
  };    

  double GetData( int i, int j) {
    if( i <= ntotal && j <= nraws){ 
      return data[i][j];
    } else {
      if( i > ntotal) cerr << "experimental_data.h: GetData:" << i << " is outside of ntotal = " << ntotal << endl;
      if( j > nraws) cerr << "experimental_data.h: GetData:" << j << " is outside of nraws in data file = " << nraws << endl;
      return 0;
    }
  };

  double SetData( int i, int j, double value) {
    if( i <= ntotal && j <= nraws){ 
      return data[i][j] = value;
    } else {
      if( i > ntotal) cerr << "experimental_data.h: SetData:" << i << " is outside of ntotal = " << ntotal << endl;
      if( j > nraws) cerr << "experimental_data.h: SetData:" << j << " is outside of nraws in data file = " << nraws << endl;
      return 0;
    }
  };

  char const * GetFilename() {
    return filename;
  };


  int SetX( int i, double value ) { 
    if( i <= ntotal  && i >= 0){ 
      x[i] = value;
      return 0;
    } else {
      cerr << "experimental_data.h: SetX:" << i << " is outside of ntotal = " << ntotal << endl;
      return 1;
    }
  };

  double GetX( int i ) { 
    if( i <= ntotal  && i >= 0){ 
      return x[i];
    } else {
      cerr << "experimental_data.h: GetX:" << i << " is outside of ntotal = " << ntotal << endl;
      return 0;
    }
  };

  int SetY( int i, double value ) {
    if( i <= ntotal && i >= 0){ 
      y[i] = value;
      return 0;
    } else {
      cerr << "experimental_data.h: SetY:" << i << " is outside of ntotal = " << ntotal << endl;
      return 1;
    }
  }; 

  double GetY( int i) { 
    if( i <= ntotal  && i >= 0){ 
      return y[i];
    } else {
      cerr << "experimental_data.h: GetY:" << i << " is outside of ntotal = " << ntotal << endl;
      return 0;
    }
  };


  int SetdY( int i, double value ) { 
    if( i <= ntotal && i >= 0){ 
      dy[i] = value;
      return 0;
    } else {
      cerr << "experimental_data.h: GetdY:" << i << " is outside of ntotal = " << ntotal << endl;
      return 1;
    }
  };

  double GetdY( int i ) { 
    if( i <= ntotal && i >= 0){ 
      return dy[i];
    } else {
      cerr << "experimental_data.h: GetdY:" << i << " is outside of ntotal = " << ntotal << endl;
      return 0;
    }
  };

  TGraphErrors* plot(); 
  TGraphErrors* plot(int color, int marker);

  int Ntotal() {return ntotal;}; // how many entries...
  int Nraws() {return nraws;}; // how many raws...

};

#endif // #ifndef __EXPERIMENTAL_DATA_
