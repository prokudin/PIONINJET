#ifndef __PARTON_H
#define __PARTON_H

#include "particle.h"
#include <iostream>

using namespace std;


enum PARTON_TYPE {
  UP  = 1,
  DOWN,
  STRANGE,
  TOP,
  BOTTOM,
  GLUON,
  ANTI_UP,
  ANTI_DOWN,
  ANTI_STRANGE,
  ANTI_TOP,
  QUARK, // generik quark...
};



class PARTON: public PARTICLE {
 private:
  enum PARTON_TYPE parton_type;

 public:
  PARTON(); // Constructor...
  virtual ~PARTON() {}; // Virtual deconstructor, Child classes might change it in future


  void   SetParton( enum PARTON_TYPE parton );

  virtual void SetCharge(double _Charge); 
  virtual void  SetSpin(double _Spin); 
  virtual void  SetMass(double _Mass); 
  virtual void  SetName(char* _Name); 


  enum PARTON_TYPE GetPartonType() const {return parton_type;};

  virtual void print( std::ostream& kuda ) { // Prints out the description...
    kuda << " *************************************" << endl;
    kuda << " *           PARTON                  *" << endl;
    kuda << " *           " << PARTICLE::GetName() << "                 *"<< endl;
    kuda << " *************************************" << endl;
    PARTICLE::print( kuda );
  }; 

};



#endif // #ifndef __HADRON_H
