#include "parton.h"
#include "constants.h"

using namespace std;

PARTON::PARTON() : PARTICLE() {};

 
void PARTON::SetParton( enum PARTON_TYPE parton_type ){ 

  switch( parton_type ){
  case UP:
    PARTICLE::SetName((char *)"UP");
    PARTICLE::SetCharge(+2./3.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(0.);
    this->parton_type = UP;
    break; // write others after...
  case QUARK:
    PARTICLE::SetName((char *)"QUARK");
    PARTICLE::SetCharge(+2./3.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(0.);
    this->parton_type = QUARK;
    break; // write others after...
  default:
    PARTICLE::SetName((char *)"QUARK");
    PARTICLE::SetCharge(+2./3.);
    PARTICLE::SetSpin(0.5);
    PARTICLE::SetMass(0.);
    this->parton_type = QUARK;
    break;
  }
}



void PARTON::SetCharge(double _Charge){
  switch( parton_type ){
  case UP:
    PARTICLE::SetCharge(+2./3.);
    break;
  case QUARK:
  default:
    PARTICLE::SetCharge(_Charge);
    break;
  }
}; 

void  PARTON::SetSpin(double _Spin){
  switch( parton_type ){
  case UP:
    PARTICLE::SetSpin(0.5);
    break;
  case QUARK:
    PARTICLE::SetSpin(0.5);
    break;
  default:
    PARTICLE::SetSpin(_Spin);
    break;
  }
}; 

void  PARTON::SetMass(double _Mass){
  switch( parton_type ){
  case UP:
    PARTICLE::SetMass(0.);
    break;
  case QUARK:
  default:
    PARTICLE::SetMass(_Mass);
    break;
  }
};
 
void  PARTON::SetName(char* _Name){
  switch( parton_type ){
  case UP:
    PARTICLE::SetName((char *)"UP");
    break;
  case QUARK:
    PARTICLE::SetName((char *)"QUARK");
    break;
  default:
    PARTICLE::SetName(_Name);
    break;
  }
}; 


