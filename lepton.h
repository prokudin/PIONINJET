#ifndef __LEPTON_H
#define __LEPTON_H

#include "particle.h"
#include <iostream>

using namespace std;


enum LEPTON_TYPE {
  ELECTRON  = 1,
  MUON
};



class LEPTON: public PARTICLE {
 private:
  enum LEPTON_TYPE lepton_type;

 public:
  LEPTON(); // Constructor...
  virtual ~LEPTON() {}; // Virtual deconstructor, Child classes might change it in future


  void   SetLepton( enum LEPTON_TYPE lepton );

  virtual void SetCharge(double _Charge); 
  virtual void  SetSpin(double _Spin); 
  virtual void  SetMass(double _Mass); 
  virtual void  SetName(char* _Name); 


  enum LEPTON_TYPE GetLeptonType() const {return lepton_type;};

  virtual void print( std::ostream& kuda ) { // Prints out the description...
    kuda << " *************************************" << endl;
    kuda << " *           LEPTON                  *" << endl;
    kuda << " *           " << PARTICLE::GetName() << "                 *"<< endl;
    kuda << " *************************************" << endl;
    PARTICLE::print( kuda );
  }; 

};



#endif // #ifndef __HADRON_H
