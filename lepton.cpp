#include "lepton.h"
#include "constants.h"

using namespace std;

LEPTON::LEPTON() : PARTICLE() {};

 
void LEPTON::SetLepton( enum LEPTON_TYPE lepton_type ){ 


  switch( lepton_type ){
  case ELECTRON:
    PARTICLE::SetName("ELECTRON");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(0.);
    this->lepton_type = ELECTRON;
    break;
  case MUON:
    PARTICLE::SetName("MUON");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(mmuon);
    this->lepton_type = ELECTRON;
    break;
  default:
    PARTICLE::SetName("ELECTRON");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(0.);
    this->lepton_type = ELECTRON;
    break;
  }
}



void LEPTON::SetCharge(double _Charge){
  switch( lepton_type ){
  case ELECTRON:
    PARTICLE::SetCharge(+1.);
    break;
  default:
    PARTICLE::SetCharge(_Charge);
    break;
  }
}; 

void  LEPTON::SetSpin(double _Spin){
  switch( lepton_type ){
  case ELECTRON:
    PARTICLE::SetSpin(0.5);
    break;
  case MUON:
    PARTICLE::SetSpin(0.5);
    break;
  default:
    PARTICLE::SetSpin(_Spin);
    break;
  }
}; 

void  LEPTON::SetMass(double _Mass){
  switch( lepton_type ){
  case ELECTRON:
    PARTICLE::SetMass(0.);
    break;
  case MUON:
    PARTICLE::SetMass(mmuon);
    break;
  default:
    PARTICLE::SetMass(_Mass);
    break;
  }
};
 
void  LEPTON::SetName(char* _Name){
  switch( lepton_type ){
  case ELECTRON:
    PARTICLE::SetName("ELECTRON");
    break;
  case MUON:
    PARTICLE::SetName("MUON");
    break;
  default:
    PARTICLE::SetName(_Name);
    break;
  }
}; 


