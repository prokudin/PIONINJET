#include <draw_result_sidis.h>
#include <cuba_sidis.h>



// Plot the results..
void draw_result_sidis(int npoints, const int ncomputations, const int ndim, int step_var, bool Coridor, integrand_t integrand, TGraph* graph[]){
  sidis& SIDIS = *sidis::Instance();

  double x = 0. , x_min =0., x_max= 0.;

  double xb_min = 0.05;
  double xb_max = 0.8;

  double xff_min = -0.7;
  double xff_max =  0.7;

  double z_min = 0.2;
  double z_max = 0.8;

  double y_min = 0.2;
  double y_max = 0.8;

  double pt_min = 0.1;
  double pt_max = 1.;

  double W_min = 6.;
  double W_max = 15.;

  double Q2_min = 2.4;
  double Q2_max = 50.;

  double xx[npoints];
  double yy[ncomputations + 1][npoints];


  double min[npoints]; // For Coridor
  double max[npoints];
  double corrx[2*npoints + 1]; // For Coridor
  double corry[2*npoints + 1];


  for(int i = 0; i <= ncomputations + 1; i++){
    for(int n = 0; n <= npoints; n++){
      yy[i][n] = 0.;
    }
  }

  double Draw_min = -0.15;
  double Draw_max = 0.15;

  double asymmetry[ncomputations];

  switch(EXPERIMENT){
  case RHIC_BRAHMS:
    xff_min = 0.25;
    xff_max = 0.35; 
    pt_min  = 1.;
    pt_max  = 2.e+10;
    Draw_min = -0.1;
    Draw_max = 0.1;
    break;
  case RHIC_STAR:
    xff_min = 0.22;
    xff_max = 0.65; 
//      xff_min = -0.6;
//      xff_max = -0.27; 
    pt_min  = 2.;
    pt_max  = 10.;
    Draw_min = -0.1;
    Draw_max = 0.15;
    y_min = -1; 
    y_max = 2;
    break;
  case RHIC_ANDY:
    xff_min = 0.14;
    xff_max = 0.6; 
    pt_min  = 1.;
    pt_max  = 2.e+10;
    Draw_min = -0.1;
    Draw_max = 0.1;
    y_min = -1; 
    y_max = 2;
    break;
   default:
    cout << "Define Experiment! DRAW RESULT SIDIS cpp" << endl;
    break;
  }

  if( SIDIS.calculation == no_effect) { // UNPOLARISED
    Draw_min = 1.e-3;
    Draw_max = 1.e+2;
  }


 

 


  switch(SIDIS.dependence){
   case xf_dependence:
    x_min = xff_min;
    x_max = xff_max;
    break;
   case pt_dependence:
    x_min = pt_min;
    x_max = pt_max;
    break;
     break;
  default:
    printf("Unknown dependence in draw results draw_results_sidis.cpp !\n");
    break;
  }

  double x_step  = (x_max - x_min)/(npoints - 1);

  if( step_var == step_log) { // log scale
    x_step  = (log( x_max ) - log( x_min ))/float(npoints - 1);
  }


// MAIN computation
  for(int i = 0; i < npoints; ++i){
    
    if(i == 0)
      x = x_min;
    else {
      x += x_step; 
      if( step_var == step_log)  // log scale
	x = exp( log( x_min ) + float(i) * x_step );
    }
    
    xx[i] = x;
    
//        SIDIS.Params.Transversity.UseSet(1); // we do not loop over all sets...
// 	    SIDIS.Params.Collins.UseSet(1); // we do not loop over all sets...
// 	    cout << "Set used in CORIDOR = "<< 27 << endl; 

//  reset_hoppet(); // DON'T FORGET TO SET HOPPET
    
    Asymmetry_SIDIS(ncomputations, ndim, x, integrand, asymmetry);
 
 

    for(int k = 0; k < ncomputations; k++){
      yy[k][i] = asymmetry[k];
      yy[ncomputations][i] += asymmetry[k]; // The sum of all effects...
    }

  }

 


 
 
 
 
 
// CORRIDOR computation COLLINS_PP
    // WORK ON IT !!!!!!
   if(Coridor) { // We draw also a coridor with errors WORK on IT

      //COLLINS:
   if( SIDIS.calculation == collins_effect_pp ){



  for(int k = 0; k < set_number; k++){
  
	    SIDIS.Params.Transversity.UseSet(k); // we do not loop over all sets...
	    SIDIS.Params.Collins.UseSet(k); // we do not loop over all sets...
  
  reset_hoppet(); // DON'T FORGET TO RESET HOPPET
  
  for(int i = 0; i < npoints; ++i){
    
    if(i == 0)
      x = x_min;
    else {
      x += x_step; 
      if( step_var == step_log)  // log scale
	x = exp( log( x_min ) + float(i) * x_step );
    }
    
    xx[i] = x;
	  
	    Asymmetry_SIDIS(ncomputations, ndim, x, integrand, asymmetry);
	    cout << "Set used in CORIDOR = "<< k << endl; 

	    double sechenie = asymmetry[0];
	  
	    //cout << "HERE" << k << endl;
	    //cout << "ASY" << sechenie << endl;
	  
	   
	  if( k == 0 ){
	    min[i] = sechenie;
	    max[i] = sechenie;
	  } else {
	    if( (sechenie ) > max[i] ){
	      max[i] = sechenie;
	    }
	    if(  sechenie < min[i] ){
	      min[i] = sechenie;
	    }
	  }

    }
    
    }

	SIDIS.Params.Transversity.UseMainSet(); // main set is restored after the loop ...
	SIDIS.Params.Collins.UseMainSet(); // main set is restored after the loop ...
    }
    
    }      
// END   CORRIDOR computation COLLINS_PP 
 



 


  //Draw Coridor if needed
  if(Coridor){
    for (int i = 0; i < npoints; ++i) {
      corry[i] = min[i];
      corrx[i] = xx[i];
    }

    for (int i = 0; i < npoints; ++i) {
      corry[npoints + i] = max[npoints - i - 1];
      corrx[npoints + i] = xx[npoints - i - 1];
    }

    corry[2*npoints] = min[0];
    corrx[2*npoints] = xx[0];

    TGraph* grcor;
    grcor = new TGraph(2*npoints + 1,corrx,corry);  
    grcor->SetFillColor(16);
    grcor->SetFillColor(16);
    //   grcor->SetFillColor(19);
    //   grcor->SetFillColor(19);
    //   grcor->SetFillColor(29);
    //   grcor->SetFillColor(29);
    grcor->SetLineColor(1);
    grcor->SetLineWidth(0.1);
//     grcor->Draw("F");
//     grcor->Draw("L");

    graph[ncomputations + 2] = grcor; // this will become a coridor...
  }


  for(int k = 0; k <= ncomputations + 1; k++){
    graph[k] = new TGraph(npoints, xx, yy[k]);
    graph[k]->SetLineColor(2);
    graph[k]->SetLineWidth(3);
    graph[k]->SetTitle("");
    graph[k]->SetMinimum(Draw_min);
    graph[k]->SetMaximum(Draw_max);
  }

 
}

