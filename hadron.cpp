#include "hadron.h"

using namespace std;

HADRON::HADRON() : PARTICLE() {};

 void HADRON::SetHadron( enum HADRON_TYPE hadron_type ){


  switch( hadron_type ){
  case PROTON:
    PARTICLE::SetName("PROTON");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(mpr);
    this->hadron_type = PROTON;
    break;
  case ANTIPROTON:
    PARTICLE::SetName("ANTIPROTON");
    PARTICLE::SetCharge(-1.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(mpr);
    this->hadron_type = ANTIPROTON;
    break;
  case PION:
    PARTICLE::SetName("PION");
    PARTICLE::SetCharge(0.);
    PARTICLE::SetSpin(0.);
    PARTICLE::SetMass(mpion);
    this->hadron_type = PION;
    break;
  case KAON:
    PARTICLE::SetName("KAON");
    PARTICLE::SetCharge(0.);
    PARTICLE::SetSpin(0.);
    PARTICLE::SetMass(mkaon);
    this->hadron_type = KAON;
    break;
  case HADRONS:
    PARTICLE::SetName("HADRONS");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(0.);
    PARTICLE::SetMass(mpion);
    this->hadron_type = HADRONS;
    break;
  case COPPER: // FOR DY EXPERIMENTS
    PARTICLE::SetName("COPPER");
    PARTICLE::SetCharge(+29.);
    PARTICLE::SetSpin(0.5); // effective target...
    PARTICLE::SetMass(63.546);
    this->hadron_type = COPPER;
    break;
  case DEUTRON:
    PARTICLE::SetName("DEUTRON");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(1.);
    PARTICLE::SetMass(mdeutron);
    this->hadron_type = DEUTRON;
    break;
  case NEUTRON:
    PARTICLE::SetName("NEUTRON");
    PARTICLE::SetCharge(0.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(mneutron);
    this->hadron_type = NEUTRON;
    break;
  case JET:
    PARTICLE::SetName("JET");
    PARTICLE::SetCharge(0.);
    PARTICLE::SetSpin(0);
    PARTICLE::SetMass(0.);
    this->hadron_type = JET;
    break;
 case PHOTON:
    PARTICLE::SetName("PHOTON");
    PARTICLE::SetCharge(0.);
    PARTICLE::SetSpin(0);
    PARTICLE::SetMass(0.);
    this->hadron_type = PHOTON;
    break;
  default:
    PARTICLE::SetName("PROTON");
    PARTICLE::SetCharge(+1.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(mpr);
    this->hadron_type = PROTON;
    break;
  }
}



void HADRON::SetCharge(double _Charge){
  switch( hadron_type ){
  case PROTON:
    PARTICLE::SetCharge(+1.);
    break;
  case ANTIPROTON:
    PARTICLE::SetCharge(-1.);
    break;
  case PION:
    PARTICLE::SetCharge(_Charge);
    break;
  case HADRONS:
    PARTICLE::SetCharge(_Charge);
    break;
  case KAON:
    PARTICLE::SetCharge(_Charge);
    break;
  case COPPER:
    PARTICLE::SetCharge(+29.);
    break;
  case DEUTRON:
    PARTICLE::SetCharge(+1.);
    break;
  case NEUTRON:
    PARTICLE::SetCharge(0.);
    break;
  default:
    PARTICLE::SetCharge(_Charge);
    break;
  }
};

void  HADRON::SetSpin(double _Spin){
  switch( hadron_type ){
  case PROTON:
    PARTICLE::SetSpin(0.5);
    break;
  case ANTIPROTON:
    PARTICLE::SetSpin(0.5);
    break;
  case PION:
    PARTICLE::SetSpin(0.);
    break;
  case COPPER:
    PARTICLE::SetSpin(0.5);
    break;
  case HADRONS:
    PARTICLE::SetSpin(0.);
    break;
  case DEUTRON:
    PARTICLE::SetSpin(1.);
    break;
  case NEUTRON:
    PARTICLE::SetSpin(0.5);
    break;
  default:
    PARTICLE::SetSpin(_Spin);
    break;
  }
};

void  HADRON::SetMass(double _Mass){
  switch( hadron_type ){
  case PROTON:
    PARTICLE::SetMass(mpr);
    break;
  case ANTIPROTON:
    PARTICLE::SetMass(mpr);
    break;
  case PION:
    PARTICLE::SetMass(mpion);
    break;
  case HADRONS:
    PARTICLE::SetMass(mpion);
    break;
  case DEUTRON:
    PARTICLE::SetMass(mdeutron);
    break;
  case COPPER:
    PARTICLE::SetMass(63.546);
    break;
  case NEUTRON:
    PARTICLE::SetMass(mneutron);
    break;
  default:
    PARTICLE::SetMass(_Mass);
    break;
  }
};

void  HADRON::SetName(char* _Name){
  switch( hadron_type ){
  case PROTON:
    PARTICLE::SetName((char *)"PROTON");
    break;
  case ANTIPROTON:
    PARTICLE::SetName((char *)"ANTIPROTON");
    break;
  case PION:
    PARTICLE::SetName((char *)"PION");
    break;
  case HADRONS:
    PARTICLE::SetName((char *)"HADRONS");
    break;
  case DEUTRON:
    PARTICLE::SetName((char *)"DEUTRON");
    break;
  case COPPER:
    PARTICLE::SetName((char *)"COPPER");
    break;
  case NEUTRON:
    PARTICLE::SetName((char *)"NEUTRON");
    break;
  default:
    PARTICLE::SetName(_Name);
    break;
  }
};
