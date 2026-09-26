// /home/prokudin/cpp/dy/ 
// Functions for DY.
// Alexei Prokudin 26/04/2011

//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cfortran.h>

#include <sidis.h>
//#include <gsl/gsl_sf_erf.h>
//#include <gauss_legendre.h>
#include <pionjet_pp.h>
//#include <evolution_sivers_sidis.h>
 

using namespace std;

double Pqq( double x) {
return CF*(1.+pow2(x))/(1.-x);
};

double Pgq( double x) {
return CF*(1.+pow2(1.-x))/x;
};

double deltaPqq( double x) { // splitting for transversity
return CF*(2.*x)/(1.-x);
};



//========================================================= Unpolarised Width...x dependence
double sidis::unpolarisedwidth_x_dependence_Feng(double x, double a, double b, double n, double gamma, double eta)
{

  double dep_x = 0.;
 

// PENG AND FENG 2014 g1+ g3 (x0/x)^lambda
   dep_x =   n  + a * pow( b / x, gamma)  ;  

  return dep_x;
}


//g2 part common to all functions
double NonPertEvolFeng(double bt, double mu, double Q0) {
 sidis& SIDIS = *sidis::Instance(); 
  
 double result;

 
 double gk_bt  = SIDIS.Params.UnpolarisedGK.parameters.a_up;
 double gk_bt2 = SIDIS.Params.UnpolarisedGK.parameters.n_up;

 double bstar = bt/sqrt(1.0+ pow(bt,2)/pow(bmax,2));

 

// FENG AND PENG 2014 
 if( (result = -gk_bt/2. * log( bt/bstar ) * log(mu/Q0) )<= -20.) 
 result = -20.;

 return exp(result);
}


//========================================================= Unpolarised PDF or FF
PARTONCONTENT  sidis::unpolarised_btpart_nonpert_Feng(double bt, double x) {
   
   
  PARTONCONTENT partcontent;
// 
   partcontent.up  =  exp( -bt * bt  *
    unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedPDF.parameters.a_up, 
			Params.UnpolarisedPDF.parameters.b_up, 
			Params.UnpolarisedPDF.parameters.n_up,
			Params.UnpolarisedPDF.parameters.gamma_up, 
			Params.UnpolarisedPDF.parameters.eta_up)    );//  *

 
  partcontent.down  =  exp( -bt * bt  *
    unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedPDF.parameters.a_down, 
			Params.UnpolarisedPDF.parameters.b_down, 
			Params.UnpolarisedPDF.parameters.n_down,
			Params.UnpolarisedPDF.parameters.gamma_down, 
			Params.UnpolarisedPDF.parameters.eta_down)    );// *

 

  partcontent.anti_up =  exp( -bt * bt  *
    unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedPDF.parameters.a_anti_up, 
			Params.UnpolarisedPDF.parameters.b_anti_up, 
			Params.UnpolarisedPDF.parameters.n_anti_up ,
			Params.UnpolarisedPDF.parameters.gamma_anti_up, 
			Params.UnpolarisedPDF.parameters.eta_anti_up)    ); // *

 

  partcontent.anti_down =   exp( -bt * bt  *
    unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedPDF.parameters.a_anti_down, 
			Params.UnpolarisedPDF.parameters.b_anti_down, 
			Params.UnpolarisedPDF.parameters.n_anti_down,
			Params.UnpolarisedPDF.parameters.gamma_anti_down, 
			Params.UnpolarisedPDF.parameters.eta_anti_down)    ); // *

 
  partcontent.strange =   exp( -bt * bt  *
    unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedPDF.parameters.a_strange, 
			Params.UnpolarisedPDF.parameters.b_strange, 
			Params.UnpolarisedPDF.parameters.n_strange,
			Params.UnpolarisedPDF.parameters.gamma_strange, 
			Params.UnpolarisedPDF.parameters.eta_strange)    );//  *



  partcontent.anti_strange =   exp( -bt  *bt  *
    unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedPDF.parameters.a_anti_strange, 
			Params.UnpolarisedPDF.parameters.b_anti_strange, 
			Params.UnpolarisedPDF.parameters.n_anti_strange,
			Params.UnpolarisedPDF.parameters.gamma_anti_strange, 
			Params.UnpolarisedPDF.parameters.eta_anti_strange)    );//  *

 
  partcontent.charm     = 0.;
  partcontent.anti_charm  = 0.;
  partcontent.bottom     = 0.;
  partcontent.anti_bottom  = 0.;
  partcontent.top     = 0.;
  partcontent.anti_top  = 0.;

  return partcontent;

};




//========================================================= Unpolarised FF
PARTONCONTENT  sidis::unpolarisedff_btpart_nonpert_Feng( HADRON& h, double bt, double x) {
   
   
   
  PARTONCONTENT fragmentation;


// HERE WE NEED TO DECIDE WHICH ONE IS FAVOURED, WHICH ONE IS UNFAVOURED
   // icharge 0,1,2,3 is 0,+,-,+&-  
  double charge = 1; 

  if( h.GetCharge() == neutral || h.GetCharge() == 0. ){ charge = 0;} // the same as 3
  else if( h.GetCharge() == positive || h.GetCharge() > 0.  ){ charge = 1;}
  else if( h.GetCharge() == negative || h.GetCharge() < 0. ){ charge = 2;}


  switch( h.GetHadronType() ){
  case HADRONS:
  case PION:
  	if( h.GetCharge() == positive || h.GetCharge() >= 0.  ){ 
   	//pi+ (u bar d) only fav--> u and unvaf-->d, more complicated --> program later!
      	double fav           = exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_up, 
					Params.UnpolarisedFF.parameters.b_up, 
					Params.UnpolarisedFF.parameters.n_up,
					Params.UnpolarisedFF.parameters.gamma_up, 
					Params.UnpolarisedFF.parameters.eta_up)    );//  *
       	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)    );//  *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  unfav; // 
	fragmentation.anti_up        =  unfav; //
	fragmentation.anti_down    =   fav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   unfav; //   
  	}
 	if( h.GetCharge() == negative || h.GetCharge() <= 0.  ){ 
   	//pi- (bar u d) only fav--> d and unvaf-->u, more complicated --> program later!
      	double fav           = exp( -bt  * bt / (x*x)*
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_up, 
					Params.UnpolarisedFF.parameters.b_up, 
					Params.UnpolarisedFF.parameters.n_up,
					Params.UnpolarisedFF.parameters.gamma_up, 
					Params.UnpolarisedFF.parameters.eta_up)    );//  *
       	double unfav         = exp( -bt  *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)   );//  *
 
	fragmentation.up               =  unfav;
	fragmentation.down           =  fav; // 
	fragmentation.anti_up        =  fav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   unfav; //   
  	}
 	if( h.GetCharge() == neutral || h.GetCharge() == 0.  ){ 
   	//// pi0 = 1/2 (pi^+ + pi^-), more complicated --> program later!
      	double fav           = 0.5 *( exp( -bt *bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_up, 
					Params.UnpolarisedFF.parameters.b_up, 
					Params.UnpolarisedFF.parameters.n_up,
					Params.UnpolarisedFF.parameters.gamma_up, 
					Params.UnpolarisedFF.parameters.eta_up)    ) +
					exp( -bt  * bt / (x*x) * 
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)    )  );
      	double unfav         = exp( -bt  *bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)    );//   *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  fav; // 
	fragmentation.anti_up        =  fav; //
	fragmentation.anti_down    =   fav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   unfav; //   
  	}

  break;
  case KAON:
 	if( h.GetCharge() == positive || h.GetCharge() >= 0.  ){ 
   	//K+ (u bar s) only fav--> u and unvaf-->d, more complicated --> program later!
   	      	double fav           = exp( -bt  * bt / (x*x)*
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_up, 
					Params.UnpolarisedFF.parameters.b_up, 
					Params.UnpolarisedFF.parameters.n_up,
					Params.UnpolarisedFF.parameters.gamma_up, 
					Params.UnpolarisedFF.parameters.eta_up)    );//  *

        	double unfav         = exp( -bt  * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)     );//  *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  unfav; // 
	fragmentation.anti_up        =  unfav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   fav; //   
  	}
 	if( h.GetCharge() == negative || h.GetCharge() <= 0.  ){ 
   	//K- (bar s d) only fav--> d and unvaf-->u, more complicated --> program later!
      	double fav           = exp( -bt  *bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_up, 
					Params.UnpolarisedFF.parameters.b_up, 
					Params.UnpolarisedFF.parameters.n_up,
					Params.UnpolarisedFF.parameters.gamma_up, 
					Params.UnpolarisedFF.parameters.eta_up)   );//  *
       	double unfav         = exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)   );//  *
 
	fragmentation.up               =  unfav;
	fragmentation.down           =  fav; // 
	fragmentation.anti_up        =  unfav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   fav; // 
	fragmentation.anti_strange =   unfav; //   
  	}
 	if( h.GetCharge() == neutral || h.GetCharge() == 0.  ){ 
   	//// pi0 = 1/2 (K^+ + K^-), more complicated --> program later!
      	double fav           = 0.5 *( exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_up, 
					Params.UnpolarisedFF.parameters.b_up, 
					Params.UnpolarisedFF.parameters.n_up,
					Params.UnpolarisedFF.parameters.gamma_up, 
					Params.UnpolarisedFF.parameters.eta_up)    ) +
					exp( -bt *bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)    ) );
      	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, Params.UnpolarisedFF.parameters.a_down, 
					Params.UnpolarisedFF.parameters.b_down, 
					Params.UnpolarisedFF.parameters.n_down,
					Params.UnpolarisedFF.parameters.gamma_down, 
					Params.UnpolarisedFF.parameters.eta_down)    );//  *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  unfav; // 
	fragmentation.anti_up        =  fav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   fav; // 
	fragmentation.anti_strange =   fav; //   
  	}




  break;
  default:
    cout << " No case for hadron type " << h.GetHadronType() << " is programmed in sidis::unpolarisedff_btpart_nonpert( HADRON& h, double bt, double x) evolution_sivers_sidis.cpp "<< endl;
  break;
  }

    
  fragmentation.charm     = 0.;
  fragmentation.anti_charm  = 0.;
  fragmentation.bottom     = 0.;
  fragmentation.anti_bottom  = 0.;
  fragmentation.top     = 0.;
  fragmentation.anti_top  = 0.;

  return fragmentation;

};




//========================================================= Collins FF
PARTONCONTENT  sidis::collinsff_btpart_nonpert_Feng( HADRON& h, double bt, double x) {
   
   
   
  PARTONCONTENT fragmentation;


// HERE WE NEED TO DECIDE WHICH ONE IS FAVOURED, WHICH ONE IS UNFAVOURED
   // icharge 0,1,2,3 is 0,+,-,+&-  
  double charge = 1; 

  if( h.GetCharge() == neutral || h.GetCharge() == 0. ){ charge = 0;} // the same as 3
  else if( h.GetCharge() == positive || h.GetCharge() > 0.  ){ charge = 1;}
  else if( h.GetCharge() == negative || h.GetCharge() < 0. ){ charge = 2;}


  switch( h.GetHadronType() ){
  case HADRONS:
  case PION:
  	if( h.GetCharge() == positive || h.GetCharge() >= 0.  ){ 
   	//pi+ (u bar d) only fav--> u and unvaf-->d, more complicated --> program later!
      	double fav           = exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_up - Params.Collins.parameters.m2_up,
					0., 
					0.)    );//  *
       	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    );//  *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  unfav; // 
	fragmentation.anti_up        =  unfav; //
	fragmentation.anti_down    =   fav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   unfav; //   
  	}
 	if( h.GetCharge() == negative || h.GetCharge() <= 0.  ){ 
   	//pi- (bar u d) only fav--> d and unvaf-->u, more complicated --> program later!
      	double fav           = exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_up - Params.Collins.parameters.m2_up,
					0., 
					0.)    );//  *
       	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    );//  *
 
	fragmentation.up               =  unfav;
	fragmentation.down           =  fav; // 
	fragmentation.anti_up        =  fav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   unfav; //   
  	}
 	if( h.GetCharge() == neutral || h.GetCharge() == 0.  ){ 
   	//// pi0 = 1/2 (pi^+ + pi^-), more complicated --> program later!
      	double fav           = 0.5 *( exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_up - Params.Collins.parameters.m2_up,
					0., 
					0.)    ) +
					exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    )  );
      	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    );//   *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  fav; // 
	fragmentation.anti_up        =  fav; //
	fragmentation.anti_down    =   fav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   unfav; //   
  	}

  break;
  case KAON:
 	if( h.GetCharge() == positive || h.GetCharge() >= 0.  ){ 
   	//K+ (u bar s) only fav--> u and unvaf-->d, more complicated --> program later!
   	      	double fav           = exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_up - Params.Collins.parameters.m2_up,
					0., 
					0.)    );//  *

        	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    );//  *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  unfav; // 
	fragmentation.anti_up        =  unfav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   unfav; // 
	fragmentation.anti_strange =   fav; //   
  	}
 	if( h.GetCharge() == negative || h.GetCharge() <= 0.  ){ 
   	//K- (bar s d) only fav--> d and unvaf-->u, more complicated --> program later!
      	double fav           = exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_up - Params.Collins.parameters.m2_up,
					0., 
					0.)    );//  *
       	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    );//  *
 
	fragmentation.up               =  unfav;
	fragmentation.down           =  fav; // 
	fragmentation.anti_up        =  unfav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   fav; // 
	fragmentation.anti_strange =   unfav; //   
  	}
 	if( h.GetCharge() == neutral || h.GetCharge() == 0.  ){ 
   	//// pi0 = 1/2 (K^+ + K^-), more complicated --> program later!
      	double fav           = 0.5 *( exp( -bt * bt / (x*x) *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_up - Params.Collins.parameters.m2_up,
					0., 
					0.)    ) +
					exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    ) );
      	double unfav         = exp( -bt *bt / (x*x)  *
    					unpolarisedwidth_x_dependence_Feng(x, 0., 
					0., 
					Params.UnpolarisedFF.parameters.n_down - Params.Collins.parameters.m2_down,
					0., 
					0.)    );//  *
 
	fragmentation.up               =  fav;
	fragmentation.down           =  unfav; // 
	fragmentation.anti_up        =  fav; //
	fragmentation.anti_down    =   unfav; //  
	fragmentation.strange        =   fav; // 
	fragmentation.anti_strange =   fav; //   
  	}




  break;
  default:
    cout << " No case for hadron type " << h.GetHadronType() << " is programmed in sidis::unpolarisedff_btpart_nonpert( HADRON& h, double bt, double x) evolution_sivers_sidis.cpp "<< endl;
  break;
  }

    
  fragmentation.charm     = 0.;
  fragmentation.anti_charm  = 0.;
  fragmentation.bottom     = 0.;
  fragmentation.anti_bottom  = 0.;
  fragmentation.top     = 0.;
  fragmentation.anti_top  = 0.;

  return fragmentation;

};
 


// c--------------------------------------------------------------
// c     alphas(q) -- strong coupling constant
// c--------------------------------------------------------------      
// checked with zhongbo's code 4/10/16
double alphasQCD_JET( double q){
      int nf;
      double lambda, alphasQCD;
      double mb=4.5;

      if(q <= mb){
         nf=4;
         lambda = 0.326;}
      else{
         nf=5;
         lambda=0.226;
      };

      double b0=11.-2./3. * nf;
      double b1=51.-19./3. * nf;
      double q2=q*q;

      double lambda2=lambda*lambda;
      double tt=log(q2/lambda2);

      alphasQCD=4. * PI/(b0*tt)*(1.-2.*b1/(b0*b0)*log(tt)/tt);

      return alphasQCD;
};


// C-------------------------------------------------------------------
// checked with zhongbo's code 4/10/16
int Upp_JET(double s,double t,double u, double* WQ1,double* WQ2,double* WQ3,
		double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8){
       
      double Nc=3.;

// c     qq' ->qq'
      *WQ1=(pow2(Nc)-1.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t);
// c     qq  ->qq
      *WQ2=(pow2(Nc)-1.)/(2.*pow2(Nc))*( (s*s+u*u)/(t*t)+(s*s+t*t)/(u*u) )
          -(pow2(Nc)-1.)/pow(Nc,3)*(s*s)/(t*u);
// c     qqb ->q'qb'
      *WQ3=(pow2(Nc)-1.)/(2.*pow2(Nc))*(t*t+u*u)/(s*s); 
// c     qqb ->qqb
      *WQ4=(pow2(Nc)-1.)/(2.*pow2(Nc))*( (s*s+u*u)/(t*t)+(t*t+u*u)/(s*s) )
          -(pow2(Nc)-1.)/pow(Nc,3)*(u*u)/(s*t);
// c     qqb ->gg
      *WQ5=pow2(pow2(Nc)-1.)/(2.*pow(Nc,3))*(u/t+t/u)
          -(pow2(Nc)-1.)/Nc*(t*t+u*u)/(s*s); // 
// c     gg  ->qqb
      *WQ6=1./(2.*Nc)*(t/u+u/t)-Nc/(pow2(Nc)-1.)*(t*t+u*u)/(s*s);
// c     qg  ->qg
      *WQ7=(pow2(Nc)-1.)/(2.*pow2(Nc))*(-s/u-u/s)+(s*s+u*u)/(t*t);
// c     gg  ->gg
      *WQ8=4.*pow2(Nc)/(pow2(Nc)-1.)*(3.-t*u/(s*s)-s*u/(t*t)-s*t/(u*u));


      return 0;
};





// C-------------------------------------------------------------------
// C     unpolarized partonic cross section
// C-------------------------------------------------------------------
// checked with zhongbo's code 4/10/16
      double UUSIG_JET(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
){
       
      //real*8 HU1,HUB1,HD1,HDB1,HS1,HSB1
      double SIG1,SIG2,SIG3,SIG4;
      double WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8;
      double WT1,WT2,WT3,WT4,WT5,WT6,WT7,WT8;
 

      Upp_JET(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4,&WQ5,&WQ6,&WQ7,&WQ8);
      Upp_JET(s,u,t,&WT1,&WT2,&WT3,&WT4,&WT5,&WT6,&WT7,&WT8);

      SIG1=WQ2*(U1*U2*fu+UB1*UB2*fub+D1*D2*fd+DB1*DB2*fdb
              +S1*S2*fs+SB1*SB2*fsb)
         +WQ3*(U1*UB2*(fd+fs)+D1*DB2*(fu+fs)+S1*SB2*(fu+fd)
              +UB1*U2*(fdb+fsb)+DB1*D2*(fub+fsb)
              +SB1*S2*(fub+fdb))
         +WT3*(UB1*U2*(fd+fs)+DB1*D2*(fu+fs)+SB1*S2*(fu+fd)
              +U1*UB2*(fdb+fsb)+D1*DB2*(fub+fsb)
              +S1*SB2*(fub+fdb));

      SIG2=WQ4*(U1*UB2*fu+D1*DB2*fd+S1*SB2*fs
              +UB1*U2*fub+DB1*D2*fdb+SB1*S2*fsb)
         +WT4*(UB1*U2*fu+DB1*D2*fd+SB1*S2*fs
              +U1*UB2*fub+D1*DB2*fdb+S1*SB2*fsb)
         +WQ1*(U1*(D2+DB2+S2+SB2)*fu
              +D1*(U2+UB2+S2+SB2)*fd
              +UB1*(D2+DB2+S2+SB2)*fub
              +DB1*(U2+UB2+S2+SB2)*fdb
              +S1*(U2+UB2+D2+DB2)*fs
              +SB1*(U2+UB2+D2+DB2)*fsb);

      SIG3=WT1*((D1+DB1+S1+SB1)*U2*fu
              +(U1+UB1+S1+SB1)*D2*fd
              +(D1+DB1+S1+SB1)*UB2*fub
              +(U1+UB1+S1+SB1)*DB2*fdb
              +(U1+UB1+D1+DB1)*S2*fs
              +(U1+UB1+D1+DB1)*SB2*fsb)
         +WQ5*(U1*UB2+D1*DB2+S1*SB2)*fg
         +WT5*(UB1*U2+DB1*D2+SB1*S2)*fg;

      SIG4=WQ6*GL1*GL2*(fu+fd+fs)
         +WT6*GL1*GL2*(fub+fdb+fsb)
         +WQ7*(GL2*(U1*fu+D1*fd+S1*fs+UB1*fub+DB1*fdb+SB1*fsb)
              +GL1*fg*(U2+UB2+D2+DB2+S2+SB2))
         +WT7*(GL1*(U2*fu+D2*fd+S2*fs+UB2*fub+DB2*fdb+SB2*fsb)
              +GL2*fg*(U1+UB1+D1+DB1+S1+SB1))
         +WQ8*GL1*GL2*fg;

      double SIG = SIG1+SIG2+SIG3+SIG4;

      return SIG;
};    




//C-------------------------------------------------------------------
// checked with zhongbo's code 4/10/16
      int Collpp(double s,double t,double u,double* WQ1,double* WQ2,double* WQ4,
	  double* WT4,double* WQ7){
      double Nc;
        
      Nc  = 3.;

//c     q+q'->q+q'
      *WQ1= (pow2(Nc)-1.)/(pow2(Nc))*(s*u)/(-t*t);
//c     q+q->q+q
      *WQ2 = (pow2(Nc)-1.)/(pow2(Nc))*( (s*u)/(-t*t)-1./Nc*s/(-t)) ;
//c     q+qb->q+qb, qb+q->qb+q
      *WQ4 = (pow2(Nc)-1.)/(pow2(Nc))*( (s*u)/(-t*t)+1./Nc*u/t); 
//c     qb+q->q+qb, q+qb->qb+q
      *WT4 =  -(pow2(Nc)-1.)/(pow(Nc,3)) ;
//c     q+g->q+g
      *WQ7 = (pow2(Nc)-1.)/(pow2(Nc))+2.*(s*u)/(-t*t);
//-----all other channels are equal to zero
//c     q'q->qq', qqb->q'qb', qqb->gg, gg->qqb, gg->gg

 
      return 0;
};



 

//c-------------------------------------------------------------------
//c     Collins partonic cross section
//c-------------------------------------------------------------------
// checked with zhongbo's code 4/10/16
double CollSIG(double s,double t,double u,
		double hU1,double hD1,double hUB1,double hDB1,double hS1,double hSB1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double Hu,double Hub,double Hd,double Hdb,double Hs,double Hsb
){
       
      
      double SIG1,SIG2,SIG3,SIG4;
      double WQ1,WQ2,WQ3,WQ4,WT4,WQ7;
 
 
     Collpp(s,t,u,&WQ1,&WQ2,&WQ4,&WT4,&WQ7);
 
 
 
     SIG1=WQ2*(hU1*U2*Hu+hUB1*UB2*Hub+hD1*D2*Hd+hDB1*DB2*Hdb
              +hS1*S2*Hs+hSB1*SB2*Hsb);
     
     SIG2=WQ1*(hU1*(D2+DB2+S2+SB2)*Hu
              +hD1*(U2+UB2+S2+SB2)*Hd
              +hUB1*(D2+DB2+S2+SB2)*Hub
              +hDB1*(U2+UB2+S2+SB2)*Hdb
              +hS1*(U2+UB2+D2+DB2)*Hs
              +hSB1*(U2+UB2+D2+DB2)*Hsb);
     
     SIG3=WQ4*(hU1*UB2*Hu+hD1*DB2*Hd+hUB1*U2*Hub
              +hDB1*D2*Hdb+hS1*SB2*Hs+hSB1*S2*Hsb)
         +WT4*(hUB1*U2*Hu+hDB1*D2*Hd+hU1*UB2*Hub
              +hD1*DB2*Hdb+hSB1*S2*Hs+hS1*SB2*Hsb);
	 
     SIG4=WQ7*GL2*(hU1*Hu+hD1*Hd+hUB1*Hub+hDB1*Hdb
              +hS1*Hs+hSB1*Hsb);
     
     
      
      double SIG=SIG1+SIG2+SIG3+SIG4;
 
      return SIG;
       
}







//========================================================= Approximate PION IN JET in PP for jet production...
// checked with zhongbo's code 4/10/16
int sidis::approximate_pionjet_pp_calculation( double xb, double y, double zc, double pt, double jt, double Q2 )
{

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);


  double Q  = sqrt( Q2 );

  double FU =0.,FD =0.,FUB =0.,FDB =0.,FS =0.,FSB =0.,FGL =0.; 
  double U1 =0.,D1 =0.,UB1 =0.,DB1 =0.,S1 =0.,SB1 =0.,GL1 =0.; 
  double U2 =0.,D2 =0.,UB2 =0.,DB2 =0.,S2 =0.,SB2 =0.,GL2 =0.; 

  

 
  double cross_section;
 
  double SS = s;
  double UU = - pt * sqrts * exp(y);
  double TT = - pt * sqrts * exp(-y);

  double xa = -xb*UU/(xb*SS +TT);


  if (xa >= 0.99 || xb >= 0.99) {sidis::cross_section = 0.; return 0;}; // otherwise we will hit x=1 point and sometimes (1-x)^-N --> infinity

  double  sh = xa * xb * SS;
  double  uh = xb * UU;
  double  th = xa * TT;

   
  
  double coeff, coeff_coll, coeff_unp;

 
  switch( effect ){
  case collins_effect: // COLLINS

    switch( Beam.GetPolarization() ){
    case transverse:
      TransversityDistribution( Beam , xa, Q2 ); //  Transversity & LO

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
  
  
      Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 
 

      CollinsDistribution( Produced, zc, Q2); // Collins Fragmentation ??? What function? pt dependent or not? ALEXEI PROKUDIN 4/10/16
      Produced.fragmentation = product(0.5,Produced.fragmentation); // I have factor "2" in collins(z) -> N(z) = 2 z^a (1-z)^b (a+b)^(a+b)/a^a/b^b Zhongbo uses collins(z)/2


      FU  = Produced.fragmentation.up;
      FD  = Produced.fragmentation.down;     
      FUB = Produced.fragmentation.anti_up;   
      FDB = Produced.fragmentation.anti_down;  
      FS  = Produced.fragmentation.strange;   
      FSB = Produced.fragmentation.anti_strange;  

       


       
      break;
    default:
	cerr << " No unp or long here!! int sidis::approximate_pionjet_pp_calculation(  )   " << endl;
    break;

    };
    
    if(!VALENCE_SIDIS){ 
      cross_section = CollSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1, 
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB  );
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
       cross_section = CollSIG( sh,  th,  uh,
		  U1 + UB1,  D1 + DB1, UB1, DB1,  S1,  SB1,  
		  U2 + UB2,  D2 + DB2, UB2, DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB  );
    }


    switch( dependence ){
    case pt_dependence:
    default: // x, z dependence
       coeff_coll = 2.*pt/(xb*SS+TT)/(xb*xa) ;
       coeff_coll *= 2.*jt*sqrt(2.*exp(1.))*jt/(sqrt(Params.Collins.parameters.m2_up))*
        exp(-pow2(jt)/Params.Collins.parameters.m2_up)*1./(PI*Params.GetPtq2Average())*
        exp(-pow2(jt)/Params.GetPtq2Average()); // jt is from d^2 jt
     break;
    }
            
    coeff = coeff_coll;

    break;

  case no_effect: // UNPOLARISED


    Unpolarised( Beam , xa, Q2 );  // 1st hadron

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
      GL1 = Beam.partcontent.glu; 
 
    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 
 

      Fragmentation( Produced, zc, Q2);


      FU  = Produced.fragmentation.up;
      FD  = Produced.fragmentation.down;     
      FUB = Produced.fragmentation.anti_up;   
      FDB = Produced.fragmentation.anti_down;  
      FS  = Produced.fragmentation.strange;   
      FSB = Produced.fragmentation.anti_strange;  
      FGL = Produced.fragmentation.glu; 


    if(!VALENCE_SIDIS){  
      cross_section = UUSIG_JET( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);  
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
      cross_section = UUSIG_JET( sh,  th,  uh,
		  U1 + UB1,  D1 + DB1, UB1, DB1,  S1,  SB1,  GL1,
		  U2 + UB2,  D2 + DB2, UB2, DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    }

    switch( dependence ){
    case pt_dependence:
    default: // x, z dependence
      coeff_unp = 2.*pt/(xb*SS+TT)/(xb*xa);
      coeff_unp *= 2.*jt/(PI*Params.GetPtq2Average())*
        exp(-pow2(jt)/Params.GetPtq2Average()) ; // jt is from d^2 jt
     break;
    }

    coeff = coeff_unp;

    break;
  default: // 
    
    cout << "Error in sidis::approximate_pionjet_pp_calculation(   ) no effect is defined" << endl;
    return 1;
    
    break;
  };

  
  sidis::cross_section = cross_section;  

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) ) / s ; 

  return 0;
}









//Integrand for PP --> PION in JET X twist-3 approach AP 4/13/2016
static int Integrand_PP_PIONJET(const int *ndim, const double xxx[],
                     const int *ncomp, double ff[], const double userdata[])
{
  sidis& SIDIS = *sidis::Instance();

  double xa = 1., xb = 1., y = 1., xF = 1.; // y here stands for rapidity!, W will stand for jT!
 
  double SS = 0, TT = 0, UU = 0, sh = 0, uh = 0, th = 0;

  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s); 

  double jacobian = 1.;

  double Q2 = 0.;

  

  double zc = userdata[0];
  double pt = userdata[1];
  double jt = userdata[2];
  double y_min = userdata[3];
  double y_max = userdata[4];
  double jt_min = userdata[5];
  double jt_max = userdata[6];

  double xb_min = 0.;
  double xb_max = 1.;
  

  // cuts in rapidity will be used for integration in rapidity!
 
  for( int i = 0; i < *ncomp; i++) ff[i] = 1.e-30; // Important!!!! otherwise the noise will count to the integral
 
  // xxx is n dim cube!!!
  switch(SIDIS.dependence){
  case pt_dependence: // pt_jet dependence -> integrate over xb, jt
    pt = pt_used;
    

    y = y_min + (y_max - y_min) * xxx[0];
    xb = xb_min + (xb_max - xb_min) * xxx[1];
    jt = jt_min + (jt_max - jt_min) * xxx[2];
     

    Q2 = pt*pt;
 
    

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    xa = -xb*UU/(xb*SS +TT);

    sh = xa * xb * SS;
    uh = xb * UU;
    th = xa * TT;


    jacobian = (y_max - y_min)*(xb_max - xb_min)*(jt_max - jt_min);


    break;
  case jt_dependence: // jt dependence
    jt = W_used;


    y = y_min + (y_max - y_min) * xxx[0];
    xb = xb_min + (xb_max - xb_min) * xxx[1];
    
     

    Q2 = pt*pt;
 
    

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    xa = -xb*UU/(xb*SS +TT);

    sh = xa * xb * SS;
    uh = xb * UU;
    th = xa * TT;


    jacobian = (y_max - y_min)*(xb_max - xb_min);
 
 
    break;
  case z_dependence: // zc dependence
    zc = z_used;

    
    y = y_min + (y_max - y_min) * xxx[0];
    xb = xb_min + (xb_max - xb_min) * xxx[1];
    jt = jt_min + (jt_max - jt_min) * xxx[2];
     

    Q2 = pt*pt;
 
    

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    xa = -xb*UU/(xb*SS +TT);

    sh = xa * xb * SS;
    uh = xb * UU;
    th = xa * TT;


    jacobian = (y_max - y_min)*(xb_max - xb_min)*(jt_max - jt_min);
 
 
    break;
  default:
    printf("Unknown dependence in Integrand_PP_PIONJET !\n");
    break;
  }


  if(
     (-UU/(SS+TT) < xa && xa  < 1.) &&
     (-TT/(SS+UU) < xb && xb  < 1.) &&
     SIDIS.Q2_cut(Q2) ){

    if ( cuba_verbose ) {
      switch(SIDIS.dependence){
      case kinematical_dependence:
      case pt_dependence:
      case z_dependence:
      case x_dependence:
      case y_dependence:
      case jt_dependence:

 	xa_tree = xa;
 	xb_tree  = xb;
 	xf_tree  = xa - xb; // ????
 	zc_tree  = zc;
 	pt_tree = pt;
 	jt_tree = jt;
 	y_tree = y;
 	tree->Fill();
	break;
      default:
	printf("Unknown dependence in Integrand_Sivers_sidis cuba_sivers_sidis.cpp !\n");
	break;
      }
    }
  


    if(SIDIS.calculation ==  collins_effect_pp){ //Collins pion in jet
 
       //Collins
      SIDIS.effect = collins_effect; // NEED TO BE MORE FLEXIBLE
      if( SIDIS.approximate_pionjet_pp_calculation( xb, y, zc, pt, jt, Q2 )  || 0){
	cerr << "approximate_pionjet_pp.cpp: ERROR wrong result!" << endl;     
      }
      
      ff[1] =  SIDIS.GetCrossSection() * jacobian;
       
      
      //Unpolarized
      SIDIS.effect = no_effect;
      if( SIDIS.approximate_pionjet_pp_calculation( xb, y, zc, pt, jt, Q2 )  || 0){
	cerr << "approximate_pionjet_pp.cpp: ERROR wrong result!" << endl;     
      }
      ff[0] =  SIDIS.GetCrossSection() * jacobian;
      
      

    }
 

  }
return 0;
}



 

#define USERDATA NULL
double EPSREL_JET = 1.e-2; // MAIN
double EPSABS_JET = 1.e-9;
 
 
//double EPSREL_JET = 1.e-1; // MAIN COURSE RESULT?
//double EPSABS_JET = 1.e-7;
 

int VERBOSE_JET = 0;
 
#define LAST 1
#define SEED 0
#define NVEC 1

 
 
#define NBATCH 100000
#define GRIDNO 0
#define STATEFILE NULL
void *SPIN = NULL;

#define NNEW 1000
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

 int MINEVAL_JET = 100; // MAIN
 int MAXEVAL_JET = 1500000; // MAIN
 int NSTART_JET  = 1000; // MAIN
 int NINCREASE_JET = 200; // MAIN


// Routine that calls the Vegas integration, calculates asymmetry as function of zc or pt
double calc_pionjet(double zc, double pt, double jt, double y_min, double y_max,  double jt_min, double jt_max){ 


int       neval, fail;
double integral[NCOMP_JET], error[NCOMP_JET], prob[NCOMP_JET];


double mydata[7];
// [0] is x , [1] is Q2, [2] is kt
mydata[0]=zc;
mydata[1]=pt;
mydata[2]=jt;
mydata[3]=y_min;
mydata[4]=y_max;
mydata[5]=jt_min;
mydata[6]=jt_max;

cubacores(0,1000); // uses only ONE CORE FOR CALCULATIONS
 
Vegas(NDIM_JET, NCOMP_JET, (integrand_t) Integrand_PP_PIONJET, mydata, NVEC,
      EPSREL_JET, EPSABS_JET, VERBOSE_JET, SEED, 
      MINEVAL_JET, MAXEVAL_JET, NSTART_JET, NINCREASE_JET, NBATCH,
      GRIDNO, STATEFILE, SPIN, 
      &neval, &fail, integral, error, prob); // 2016, NVEC should follow mydata 


//cubawait(&SPIN);

 return integral[1]/integral[0];//  THIS IS ASYMMETRY
//
// return integral[0];//  unpolarised

  

}







//========================================================= Approximate PION IN JET in PP for jet production at NLO accuracy...
// This function is in bt, must be integrated in order to have answer in jt
int sidis::approximate_pionjet_pp_calculation_NLO( double xb, double y, double zc, double zhat, double bt, double pt, double jt, double Q2 )
{

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);


  double Q  = sqrt( Q2 );

  double FU =0.,FD =0.,FUB =0.,FDB =0.,FS =0.,FSB =0.,FGL =0.; 
  double U1 =0.,D1 =0.,UB1 =0.,DB1 =0.,S1 =0.,SB1 =0.,GL1 =0.; 
  double U2 =0.,D2 =0.,UB2 =0.,DB2 =0.,S2 =0.,SB2 =0.,GL2 =0.; 

  

 
  double cross_section;
 
  double SS = s;
  double UU = - pt * sqrts * exp(y);
  double TT = - pt * sqrts * exp(-y);

  double xa = -xb*UU/(xb*SS +TT);


  if (xa >= 0.99 || xb >= 0.99) {sidis::cross_section = 0.; return 0;}; // otherwise we will hit x=1 point and sometimes (1-x)^-N --> infinity

  double  sh = xa * xb * SS;
  double  uh = xb * UU;
  double  th = xa * TT;

     
  double mu = Q ;
  double xi = mu * mu;

  double bstar = bt/sqrt(1.0+ pow(bt,2)/pow(bmax,2));
  double mub2 = pow(C1/bstar,2);
  double mub  = C1/bstar;

  
  
  double as = AlphaS(mub);
 
 
  PARTONCONTENT ff_z;
  PARTONCONTENT ff_zhat;

 
  
  double z1 = zc/zhat;



  double evolution;   

  double Q2ff = mub2; 
  
  double coeff, coeff_coll, coeff_unp;
  
  double coefficient;  

  double relation_to_trento;

 
  switch( effect ){
  case collins_effect: // COLLINS

    switch( Beam.GetPolarization() ){
    case transverse:
      TransversityDistributionNLO( Beam , xa, Q2 ); //  Transversity & NLO

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
  
  
      Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 
 
// Collins fragmentation at NLO is here:
      CollinsDistributionNLO( Produced, zc, Q2ff); // Fragmentation EVOLUTION mub   
   
      ff_z = Produced.fragmentation;
   
      CollinsDistributionNLO( Produced, zhat, Q2ff); // Collins Fragmentation EVOLUTION mub  in z hat  
   
      ff_zhat = Produced.fragmentation;
    
//CONVOLUTIONS ARE HERE   
      Produced.fragmentation = product( (1.+ as/PI*(-2.*CF))/(1.-zc), ff_z); // the part with delta function
      //I need to divide by 1-z in order to have correct result after integration
      ff_zhat = product( as/PI/zhat *( deltaPqq(z1)*log(z1)),ff_zhat); //part with deltaPqq
      Produced.fragmentation = sum(ff_zhat,  Produced.fragmentation); //done
//END CONVOLUTIONS  
   
  
      evolution = PertEvolFeng(mu,bt)*NonPertEvolFeng(bt,mu,Q0);  


    
      Produced.fragmentation = product( collinsff_btpart_nonpert_Feng(Produced, bt, zc), Produced.fragmentation);
 
    
 
      Produced.fragmentation = product( evolution , Produced.fragmentation); // FULL EVOLUTION


      coefficient = bt*bt/(2.*PI*zc*zc)/zc * TMath::BesselJ1( bt * jt/zc );  

      relation_to_trento = - zc;   
      Produced.fragmentation = product( relation_to_trento*coefficient, Produced.fragmentation); // FULL evolution



      FU  = Produced.fragmentation.up;
      FD  = Produced.fragmentation.down;     
      FUB = Produced.fragmentation.anti_up;   
      FDB = Produced.fragmentation.anti_down;  
      FS  = Produced.fragmentation.strange;   
      FSB = Produced.fragmentation.anti_strange;  
       
      break;
    default:
	cerr << " No unp or long here!! int sidis::approximate_pionjet_pp_calculation_NLO(  )   " << endl;
    break;

    };
    
    if(!VALENCE_SIDIS){ 
      cross_section = CollSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1, 
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB  );
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
       cross_section = CollSIG( sh,  th,  uh,
		  U1 + UB1,  D1 + DB1, UB1, DB1,  S1,  SB1,  
		  U2 + UB2,  D2 + DB2, UB2, DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB  );
    }


    switch( dependence ){
    case pt_dependence:
    default: // x, z dependence
       coeff_coll = 2.*pt/(xb*SS+TT)/(xb*xa) ;
       coeff_coll *= 2.*jt;
      break;
    }
            
    coeff = coeff_coll;

    break;

  case no_effect: // UNPOLARISED


    Unpolarised( Beam , xa, Q2 );  // 1st hadron

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
      GL1 = Beam.partcontent.glu; 
 
    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 
 
// Fragmentation at NLO
      Fragmentation( Produced, zc, Q2ff); // Fragmentation EVOLUTION mub   
   
      ff_z = Produced.fragmentation;
      
      Fragmentation( Produced, zhat, Q2ff); // Fragmentation EVOLUTION mub  in z hat  
   
      ff_zhat = Produced.fragmentation;    
//CONVOLUTIONS ARE HERE
      Produced.fragmentation = product( (1.+ as/PI*(-2.*CF))/(1.-zc), ff_z); // the part with delta function
      //I need to divide by 1-z in order to have correct result after integration
      Produced.fragmentation = sum(as/PI * ff_zhat.glu/zhat *( CF/2. * z1     + Pgq(z1)*log(z1)),  Produced.fragmentation); // with gluon splitting    
      ff_zhat = product( as/PI/zhat *( CF/2. * (1.-z1) + Pqq(z1)*log(z1)),ff_zhat); //part with Pqq
      Produced.fragmentation = sum(ff_zhat,  Produced.fragmentation); 
//END CONVOLUTIONS  
  

// EVOLUTION PART
 
     evolution = PertEvolFeng(mu,bt)*NonPertEvolFeng(bt,mu,Q0);  

     
  
     Produced.fragmentation = product( unpolarisedff_btpart_nonpert_Feng(Produced, bt, zc), Produced.fragmentation);
 
      
 
     Produced.fragmentation = product( evolution , Produced.fragmentation); // FULL EVOLUTION
 
 
     coefficient = bt/(2.*PI*zc*zc) * TMath::BesselJ0( bt * jt/zc );  

   
     Produced.fragmentation = product( coefficient, Produced.fragmentation); // FULL evolution


      FU  = Produced.fragmentation.up;
      FD  = Produced.fragmentation.down;     
      FUB = Produced.fragmentation.anti_up;   
      FDB = Produced.fragmentation.anti_down;  
      FS  = Produced.fragmentation.strange;   
      FSB = Produced.fragmentation.anti_strange;  
      FGL = Produced.fragmentation.glu; 


    if(!VALENCE_SIDIS){  
      cross_section = UUSIG_JET( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);  
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
      cross_section = UUSIG_JET( sh,  th,  uh,
		  U1 + UB1,  D1 + DB1, UB1, DB1,  S1,  SB1,  GL1,
		  U2 + UB2,  D2 + DB2, UB2, DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    }

    switch( dependence ){
    case pt_dependence:
    default: // x, z dependence
      coeff_unp = 2.*pt/(xb*SS+TT)/(xb*xa);
      coeff_unp *= 2.*jt;
    break;
    }

    coeff = coeff_unp;

    break;
  default: // 
    
    cout << "Error in sidis::approximate_pionjet_pp_calculation_NLO(   ) no effect is defined" << endl;
    return 1;
    
    break;
  };

  

  sidis::cross_section = cross_section;  

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) ) / s ; 

  return 0;
}


static int Integrand_PP_PIONJET_NLO(const int *ndim, const double xxx[],
                     const int *ncomp, double ff[], const double userdata[])
{
  sidis& SIDIS = *sidis::Instance();

  double xa = 1., xb = 1., y = 1., xF = 1.; // y here stands for rapidity!, W will stand for jT!
 
  double SS = 0, TT = 0, UU = 0, sh = 0, uh = 0, th = 0;
  
  double zhat, bt;

  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s); 

  double jacobian = 1.;

  double Q2 = 0.;

  

  double zc = userdata[0];
  double pt = userdata[1];
  double jt = userdata[2];
  double y_min = userdata[3];
  double y_max = userdata[4];
  double jt_min = userdata[5];
  double jt_max = userdata[6];

  double xb_min = 0.;
  double xb_max = 1.;
 
  // minimum and maximum for integration in bt:
  //double bt_integration_min = 0.00001;
  double bt_integration_min = 0.001;
  double bt_integration_max = 10.;
 

  // cuts in rapidity will be used for integration in rapidity!
 
  for( int i = 0; i < *ncomp; i++) ff[i] = 1.e-30; // Important!!!! otherwise the noise will count to the integral
 
  // xxx is n dim cube!!!
  switch(SIDIS.dependence){
  case pt_dependence: // pt_jet dependence -> integrate over xb, jt
    pt = pt_used;
    

    y = y_min + (y_max - y_min) * xxx[0];
    xb = xb_min + (xb_max - xb_min) * xxx[1];
    jt = jt_min + (jt_max - jt_min) * xxx[2];
    zhat = zc + (1 - zc) * xxx[3];
    bt =  bt_integration_min + (bt_integration_max - bt_integration_min) * xxx[4];

    Q2 = pt*pt;
 
    

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    xa = -xb*UU/(xb*SS +TT);

    sh = xa * xb * SS;
    uh = xb * UU;
    th = xa * TT;


    jacobian = (y_max - y_min)*(xb_max - xb_min)*(jt_max - jt_min)*(1-zc)*(bt_integration_max - bt_integration_min);


    break;
  case jt_dependence: // jt dependence
    jt = W_used;


    y = y_min + (y_max - y_min) * xxx[0];
    xb = xb_min + (xb_max - xb_min) * xxx[1];
    zhat = zc + (1 - zc) * xxx[2];
    bt =  bt_integration_min + (bt_integration_max - bt_integration_min) * xxx[3];
    
     

    Q2 = pt*pt;
 
    

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    xa = -xb*UU/(xb*SS +TT);

    sh = xa * xb * SS;
    uh = xb * UU;
    th = xa * TT;


    jacobian = (y_max - y_min)*(xb_max - xb_min)*(1-zc)*(bt_integration_max - bt_integration_min);
 
 
    break;
  case z_dependence: // zc dependence
    zc = z_used;

    
    y = y_min + (y_max - y_min) * xxx[0];
    xb = xb_min + (xb_max - xb_min) * xxx[1];
    jt = jt_min + (jt_max - jt_min) * xxx[2];
    zhat = zc + (1 - zc) * xxx[3];
    bt =  bt_integration_min + (bt_integration_max - bt_integration_min) * xxx[4];

    Q2 = pt*pt;
 
    

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    xa = -xb*UU/(xb*SS +TT);

    sh = xa * xb * SS;
    uh = xb * UU;
    th = xa * TT;


    jacobian = (y_max - y_min)*(xb_max - xb_min)*(jt_max - jt_min)*(1-zc)*(bt_integration_max - bt_integration_min);
 
 
    break;
  default:
    printf("Unknown dependence in Integrand_PP_PIONJET !\n");
    break;
  }


  if(
     (-UU/(SS+TT) < xa && xa  < 1.) &&
     (-TT/(SS+UU) < xb && xb  < 1.) &&
     SIDIS.Q2_cut(Q2) //&&
     //SIDIS.xF_cut(xa-xb) // Shall we use cut in xF???
     ){

  
    if ( cuba_verbose ) {
      switch(SIDIS.dependence){
      case kinematical_dependence:
      case pt_dependence:
      case z_dependence:
      case x_dependence:
      case y_dependence:
      case jt_dependence:

 	xa_tree = xa;
 	xb_tree  = xb;
 	xf_tree  = xa - xb; // ????
 	zc_tree  = zc;
 	pt_tree = pt;
 	jt_tree = jt;
 	y_tree = y;
 	tree->Fill();
	break;
      default:
	printf("Unknown dependence in Integrand_Sivers_sidis cuba_sivers_sidis.cpp !\n");
	break;
      }
    }


    if(SIDIS.calculation ==  collins_effect_pp){ //Collins pion in jet
 
       //Collins
      SIDIS.effect = collins_effect; // NEED TO BE MORE FLEXIBLE
      if( SIDIS.approximate_pionjet_pp_calculation_NLO( xb, y, zc, zhat, bt, pt, jt, Q2 )  || 0){
	cerr << "approximate_pionjet_pp.cpp: ERROR wrong result!" << endl;     
      }
      
      ff[1] =  SIDIS.GetCrossSection() * jacobian;
       
      
      //Unpolarized
      SIDIS.effect = no_effect;
      if( SIDIS.approximate_pionjet_pp_calculation_NLO( xb, y, zc,zhat, bt, pt, jt, Q2 )  || 0){
	cerr << "approximate_pionjet_pp.cpp: ERROR wrong result!" << endl;     
      }
      ff[0] =  SIDIS.GetCrossSection() * jacobian;
      
      

    }
 

  }
return 0;
}


// Routine that calls the Vegas integration, calculates asymmetry as function of zc or pt
double calc_pionjet_NLO(double zc, double pt, double jt, double y_min, double y_max,  double jt_min, double jt_max){ 


int       neval, fail;
double integral[NCOMP_JET], error[NCOMP_JET], prob[NCOMP_JET];


double mydata[7];
// [0] is x , [1] is Q2, [2] is kt
mydata[0]=zc;
mydata[1]=pt;
mydata[2]=jt;
mydata[3]=y_min;
mydata[4]=y_max;
mydata[5]=jt_min;
mydata[6]=jt_max;

cubacores(0,1000); // uses only ONE CORE FOR CALCULATIONS
 
Vegas(NDIM_JET, NCOMP_JET, (integrand_t) Integrand_PP_PIONJET_NLO, mydata, NVEC,
      EPSREL_JET, EPSABS_JET, VERBOSE_JET, SEED, 
      MINEVAL_JET, MAXEVAL_JET, NSTART_JET, NINCREASE_JET, NBATCH,
      GRIDNO, STATEFILE, SPIN, 
      &neval, &fail, integral, error, prob); // 2016, NVEC should follow mydata 


//cubawait(&SPIN);
 return integral[1]/integral[0];// asymmetry

  
// return integral[0];// Unpolarized

}
 


