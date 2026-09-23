//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <test.h>
#include <memory.h>
#include <sstream>


using namespace std;


//========================================================= main
int main(int argc, char **argv)
{
// FRAGMENTATION IS HERE:
    
   // ffset   1,2,3 means K, KKP, BFGW
   // fforder 0,1 is LO, NLO(MSbar)
   // ihadron 1,2,3,4,5 is pi,K,h,p,n
   // icp     1,2,3 chooses between particle, anti-particle or sum of both 
   // ipi     1,2,3 is a flag for BFGW; inactive for K, KKP  
   // icharge 0,1,2,3 is 0,+,-,+&-   
   // int ffset  = 1; // Kretzer
  int ffset  = 4; // DSS
  int fforder= 1; // NLO
  int ihadron= 1; // PION
  int icp    = 1; //fragmentation.icp;
  int ipi    = 1;  
  int icharge= 1; // positive

  double dff[11];

  double z = 0.5, Q = 10.;

  dff[0] = 0.; 
  dff[1] = 0.; 
  dff[2] = 0.;
  dff[3] = 0.;
  dff[4] = 0.; 
  dff[5] = 0.;
  dff[6] = 0.; 
  dff[7] = 0.;
  dff[8] = 0.; 
  dff[9] = 0.; 
  dff[10]= 0.;

// * Subroutine "dlib" returns D(z,Q2) [NOT z*D(z,Q2)] as an array        *  
  DLIB(z,Q*Q,dff,ffset,fforder,ihadron,icharge,icp,ipi);

  printf("=====TEST DSS \n  %f %f %f %f %f %f \n",
	 dff[0],dff[1],dff[2],dff[3],dff[4],dff[5]);
  cout << "=====End TEST DSS " << endl;
}