#ifndef __PARTICLE_H
#define __PARTICLE_H
//class particle

// --- This is for root use in C++
#include "TVector3.h"
#include "TLorentzVector.h"
#include "constants.h"


// Switch to longitudinal, transverse polarization......
enum POLARIZATION {
  unp = 1,
  longitudinal,
  transverse,
  NA
};


class PARTICLE {
 private:
  char     Name[32];
  double   Charge;
  double   Spin;
  double   Mass;
  TLorentzVector Momentum4;
  TVector3 Momentum3;
  TVector3 PolarizationVector;
  enum POLARIZATION Polarization;

 public:
  PARTICLE(); 
  virtual ~PARTICLE() {}; 
  PARTICLE(char* Name, double Charge, double Spin, double Mass);

  PARTICLE& operator=(const PARTICLE& P) {
    strcpy(Name, P.Name);
    Charge = P.Charge;
    Spin = P.Spin;
    Mass = P.Mass;
    Momentum3 = P.Momentum3;
    Momentum4 = P.Momentum4;
    PolarizationVector = P.PolarizationVector;
    Polarization = P.Polarization;
    // by convention, always return *this
    return *this;
  }



  virtual void SetCharge(double _Charge) { Charge = _Charge;}; 
  virtual void  SetSpin(double _Spin) { Spin = _Spin;}; 
  virtual void  SetMass(double _Mass) { Mass = _Mass;}; 

//  const char * deals with warning: deprecated conversion from string constant to ?char*?
//  http://en.wikibooks.org/wiki/GCC_Debugging/g%2B%2B/Warnings/deprecated_conversion_from_string_constant


  virtual void  SetName(const char* _Name) { strcpy(Name, _Name);}; 



  double GetCharge() const {return Charge;}; // these are to read protected values...
  double GetSpin() const {return Spin;}; 
  double GetMass() const {return Mass;}; 

  enum POLARIZATION GetPolarization() const {return Polarization;};
  void SetPolarization( enum POLARIZATION p ) { 
    if ( GetSpin() == 0 ) Polarization = NA; 
    else Polarization = p;
  };

  virtual char*  GetName() {return Name;}
  virtual TVector3 GetMomentum3() const {return Momentum3;};
  virtual TLorentzVector GetMomentum4() const {return Momentum4;};
  virtual TVector3 GetPolarizationVector() const {return PolarizationVector;};

  //  virtual double GetEnergy() const { return mod3(Mass, Momentum);}

  virtual double GetEnergy() const { return sqrt(Mass*Mass + Momentum3*Momentum3);}

  void SetMomentum3( const TVector3&  p) { Momentum3 = p; };

  void SetMomentum4( const TLorentzVector&  p) { 
    Momentum4 = p; 
    Momentum3.SetX( p.Px() );
    Momentum3.SetY( p.Py() );
    Momentum3.SetZ( p.Pz() );
  };

  void SetPolarizationVector( const TVector3& s) { PolarizationVector = s; };

   virtual void print( std::ostream& kuda ); 
//  virtual void print();
}; 

#endif // #ifndef __PARTICLE_H
