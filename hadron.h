#ifndef __HADRON_H
#define __HADRON_H

//#include "particle.h"
#include <particle.h>
#include <iostream>

using namespace std;


enum HADRON_TYPE {
  PROTON  = 1,
  DEUTRON,
  NEUTRON,
  ANTIPROTON,
  PION,
  KAON,
  HADRONS, // NB: used for fragmentation only
  COPPER,   // USED FOR DY EXPERIMENTS --. 29 protons, 35 neutrons....
  JET, // NOTE Let us use jet as a hadron
  PHOTON // NOTE Let us use photon as a hadron
};


// Partonic content DIS
struct PARTONCONTENT {
  double up,down,anti_up,anti_down,strange,anti_strange,charm,anti_charm,bottom,anti_bottom,top,anti_top,glu;
};


inline PARTONCONTENT division(PARTONCONTENT a, PARTONCONTENT b){
  PARTONCONTENT c;

  c.up          = a.up / b.up;
  c.down        = a.down / b.down;
  c.anti_up     = a.anti_up / b.anti_up;
  c.anti_down   = a.anti_down / b.anti_down;
  c.strange     = a.strange / b.strange;
  c.anti_strange= a.anti_strange / b.anti_strange;
  if( c.charm != 0. ){c.charm       = a.charm / b.charm; } else { c.charm = 0.; } ;
  if( c.anti_charm != 0. ){c.anti_charm  = a.anti_charm / b.anti_charm;} else { c.anti_charm = 0.; }
  if( b.bottom != 0. ){c.bottom      = a.bottom / b.bottom;} else { c.bottom     = 0.;};
  if( b.anti_bottom != 0. ){c.anti_bottom = a.anti_bottom / b.anti_bottom;} else { c.anti_bottom     = 0.;};
  if( b.top != 0. ){c.top         = a.top / b.top;} else { c.top     = 0.;};
  if( b.anti_top != 0. ){c.anti_top    = a.anti_top / b.anti_top;} else { c.anti_top     = 0.;};
  c.glu         = a.glu / b.glu;

  return c;
}


inline PARTONCONTENT product(PARTONCONTENT a, PARTONCONTENT b){
  PARTONCONTENT c;

  c.up          = a.up * b.up;
  c.down        = a.down * b.down;
  c.anti_up     = a.anti_up * b.anti_up;
  c.anti_down   = a.anti_down * b.anti_down;
  c.strange     = a.strange * b.strange;
  c.anti_strange= a.anti_strange * b.anti_strange;
  c.charm       = a.charm * b.charm;
  c.anti_charm  = a.anti_charm * b.anti_charm;
  c.bottom      = a.bottom * b.bottom;
  c.anti_bottom = a.anti_bottom * b.anti_bottom;
  c.top         = a.top * b.top;
  c.anti_top    = a.anti_top * b.anti_top;
  c.glu         = a.glu * b.glu;

  return c;
}



inline PARTONCONTENT product(double a, PARTONCONTENT b){
  PARTONCONTENT c;

  c.up          = a * b.up;
  c.down        = a * b.down;
  c.anti_up     = a * b.anti_up;
  c.anti_down   = a * b.anti_down;
  c.strange     = a * b.strange;
  c.anti_strange= a * b.anti_strange;
  c.charm       = a * b.charm;
  c.anti_charm  = a * b.anti_charm;
  c.bottom      = a * b.bottom;
  c.anti_bottom = a * b.anti_bottom;
  c.top         = a * b.top;
  c.anti_top    = a * b.anti_top;
  c.glu         = a * b.glu;

  return c;
}

inline PARTONCONTENT sum(double a, PARTONCONTENT b){
  PARTONCONTENT c;

  c.up          = a + b.up;
  c.down        = a + b.down;
  c.anti_up     = a + b.anti_up;
  c.anti_down   = a + b.anti_down;
  c.strange     = a + b.strange;
  c.anti_strange= a + b.anti_strange;
  c.charm       = a + b.charm;
  c.anti_charm  = a + b.anti_charm;
  c.bottom      = a + b.bottom;
  c.anti_bottom = a + b.anti_bottom;
  c.top         = a + b.top;
  c.anti_top    = a + b.anti_top;
  c.glu         = a + b.glu;

  return c;
}

inline PARTONCONTENT sum(PARTONCONTENT a, PARTONCONTENT b){
  PARTONCONTENT c;

  c.up          = a.up + b.up;
  c.down        = a.down + b.down;
  c.anti_up     = a.anti_up + b.anti_up;
  c.anti_down   = a.anti_down + b.anti_down;
  c.strange     = a.strange + b.strange;
  c.anti_strange= a.anti_strange + b.anti_strange;
  c.charm       = a.charm + b.charm;
  c.anti_charm  = a.anti_charm + b.anti_charm;
  c.bottom      = a.bottom + b.bottom;
  c.anti_bottom = a.anti_bottom + b.anti_bottom;
  c.top         = a.top + b.top;
  c.anti_top    = a.anti_top + b.anti_top;
  c.glu         = a.glu + b.glu;

  return c;
}

inline PARTONCONTENT difference(PARTONCONTENT a, PARTONCONTENT b){
  PARTONCONTENT c;

  c.up          = a.up - b.up;
  c.down        = a.down - b.down;
  c.anti_up     = a.anti_up - b.anti_up;
  c.anti_down   = a.anti_down - b.anti_down;
  c.strange     = a.strange - b.strange;
  c.anti_strange= a.anti_strange - b.anti_strange;
  c.charm       = a.charm - b.charm;
  c.anti_charm  = a.anti_charm - b.anti_charm;
  c.bottom      = a.bottom - b.bottom;
  c.anti_bottom = a.anti_bottom - b.anti_bottom;
  c.top         = a.top - b.top;
  c.anti_top    = a.anti_top - b.anti_top;
  c.glu         = a.glu - b.glu;

  return c;
}



class HADRON: public PARTICLE {
 private:
  enum HADRON_TYPE hadron_type;

 public:
  HADRON(); // Constructor...
  virtual ~HADRON() {}; // Virtual deconstructor, Child classes might change it in future

  struct PARTONCONTENT partcontent;
  struct PARTONCONTENT fragmentation;

  void   SetHadron( enum HADRON_TYPE hadron );

  virtual void SetCharge(double _Charge);
  virtual void  SetSpin(double _Spin);
  virtual void  SetMass(double _Mass);
  virtual void  SetName(char* _Name);


  enum HADRON_TYPE GetHadronType() const {return hadron_type;};

  virtual void print( std::ostream& kuda ) { // Prints out the description...
    kuda << " *************************************" << endl;
    kuda << " *           HADRON                  *" << endl;
    kuda << " *           " << PARTICLE::GetName() << "                 *"<< endl;
    kuda << " *************************************" << endl;
    PARTICLE::print( kuda );
  };

};



#endif // #ifndef __HADRON_H
