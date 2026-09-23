#include <setexperiment_sidis.h>

int SetDependence_SIDIS(const char * argv){
  sidis& SIDIS = *sidis::Instance();


  if(strcmp(argv,"-PT") == 0)
    {
      SIDIS.dependence = pt_dependence;
    }
  else if(strcmp(argv,"-XF") == 0)
    {
      SIDIS.dependence = xf_dependence;
    }
  else if(strcmp(argv,"-Z") == 0)
    {
      SIDIS.dependence = z_dependence;
    }
  else if(strcmp(argv,"-JT") == 0)
    {
      SIDIS.dependence = jt_dependence;
    }
  else
    {
      cout << " USAGE of SetDependence_SIDIS: " << endl;
      cout << " available dependence " << endl;
      cout << "-PT"<< endl;
      cout << "-XF"<< endl;
      cout << "-Z"<< endl;
      cout << "-JT"<< endl;
      return 1;
  }

  return 0;
};

int SetHadron_SIDIS(const char * argv){
  sidis& SIDIS = *sidis::Instance();

  if(strcmp(argv,"-H-") == 0)
    {
      SIDIS.Produced.SetHadron(HADRONS);
      SIDIS.Produced.SetCharge(negative);
    }
  else if(strcmp(argv,"-H+") == 0)
    {
      SIDIS.Produced.SetHadron(HADRONS);
      SIDIS.Produced.SetCharge(positive);
    }
  else if(strcmp(argv,"-H0") == 0)
    {
      SIDIS.Produced.SetHadron(HADRONS);
      SIDIS.Produced.SetCharge(neutral);
    }
  else if(strcmp(argv,"-PI-") == 0)
    {
      SIDIS.Produced.SetHadron(PION);
      SIDIS.Produced.SetCharge(negative);
    }
  else if(strcmp(argv,"-PI+") == 0)
    {
      SIDIS.Produced.SetHadron(PION);
      SIDIS.Produced.SetCharge(positive);
    }
  else if(strcmp(argv,"-PI0") == 0)
    {
      SIDIS.Produced.SetHadron(PION);
      SIDIS.Produced.SetCharge(neutral);
    }
  else if(strcmp(argv,"-K-") == 0)
    {
      SIDIS.Produced.SetHadron(KAON);
      SIDIS.Produced.SetCharge(negative);
    }
  else if(strcmp(argv,"-K+") == 0)
    {
      SIDIS.Produced.SetHadron(KAON);
      SIDIS.Produced.SetCharge(positive);
    }
  else if(strcmp(argv,"-K0") == 0)
    {
      SIDIS.Produced.SetHadron(KAON);
      SIDIS.Produced.SetCharge(neutral);
    }
  else if(strcmp(argv,"-PIONS") == 0)
    {
      SIDIS.Produced.SetHadron(PION);
      SIDIS.Produced.SetCharge(positive); // I'll change the charge in the program
    }
  else if(strcmp(argv,"-KAONS") == 0)
    {
      SIDIS.Produced.SetHadron(KAON);
      SIDIS.Produced.SetCharge(positive); // I'll change the charge in the program
    }
  else if(strcmp(argv,"-JET") == 0)
    {
      SIDIS.Produced.SetHadron(JET);
      SIDIS.Produced.SetCharge(neutral);
    }
  else if(strcmp(argv,"-PHOTON") == 0)
    {
      SIDIS.Produced.SetHadron(PHOTON);
      SIDIS.Produced.SetCharge(neutral);
    }
  else {
      cout << " USAGE of SetHadron_SIDIS: " << endl;
      cout << " available hadrons " << endl;
      cout << "-H-"<< endl;
      cout << "-H+"<< endl;
      cout << "-H0"<< endl;
      cout << "-PI-"<< endl;
      cout << "-PI+"<< endl;
      cout << "-PI0"<< endl;
      cout << "-K-"<< endl;
      cout << "-K+"<< endl;
      cout << "-K0"<< endl;
      cout << "-PIONS"<< endl;
      cout << "-KAONS"<< endl;
      cout << "-JET"<< endl;
      cout << "-PHOTON"<< endl;
      return 1;

  }

  return 0;
};

int SetTarget_SIDIS(const char* argv){
  sidis& SIDIS = *sidis::Instance();

  if(strcmp(argv,"-PROTON") == 0)
    {
      SIDIS.Target.SetHadron(PROTON);
    }
  else if(strcmp(argv,"-DEUTRON") == 0)
    {
      SIDIS.Target.SetHadron(DEUTRON);
    }
  else if(strcmp(argv,"-NEUTRON") == 0)
    {
      SIDIS.Target.SetHadron(NEUTRON);
    }
  else
    {
      cout << " USAGE of SetTarget_SIDIS: " << endl;
      cout << " available targets " << endl;
      cout << "-PROTON" << endl;
      cout << "-DEUTRON" << endl;
      cout << "-NEUTRON" << endl;


      return 1;
    }


  return 0;
};

int SetBeam_SIDIS(const char* argv){
  sidis& SIDIS = *sidis::Instance();

  if(strcmp(argv,"-PROTON") == 0)
    {
      SIDIS.Beam.SetHadron(PROTON);
    }
  else if(strcmp(argv,"-DEUTRON") == 0)
    {
      SIDIS.Beam.SetHadron(DEUTRON);
    }
  else if(strcmp(argv,"-NEUTRON") == 0)
    {
      SIDIS.Beam.SetHadron(NEUTRON);
    }
  else
    {
      cout << " USAGE of SetBeam_SIDIS: " << endl;
      cout << " available targets " << endl;
      cout << "-PROTON" << endl;
      cout << "-DEUTRON" << endl;
      cout << "-NEUTRON" << endl;


      return 1;
    }


  return 0;
};

// void SetExperiment_SIDIS(int Experiment){ // We use it for vegas
int SetExperiment_SIDIS(const char * argv){ // We use it for vegas
  sidis& SIDIS = *sidis::Instance();

    double EE ;
    double EP;
    double energy;

  if(strcmp(argv,"-RHIC_STAR")== 0)
    {
      EXPERIMENT = RHIC_STAR;
    }
  else if(strcmp(argv,"-RHIC_BRAHMS")== 0)
    {
      EXPERIMENT = RHIC_BRAHMS;
    }
  else if(strcmp(argv,"-RHIC_ANDY")== 0)
    {
      EXPERIMENT = RHIC_ANDY;
    }
  else if(strcmp(argv,"-RHIC_PHENIX")== 0)
    {
      EXPERIMENT = RHIC_PHENIX;
    }
  else
    {
      cout << " USAGE of SetExperiment_SIDIS: " << endl;
      cout << " available experiments " << endl;
      cout << "-RHIC_STAR" << endl;
      cout << "-RHIC_BRAHMS" << endl;
      cout << "-RHIC_ANDY" << endl;
      cout << "-RHIC_PHENIX" << endl;

      return 1;
    }


  int Experiment = EXPERIMENT;

  switch (Experiment){
  case RHIC_ANDY: // RHIC PP ANDY


    EXPERIMENT = RHIC_ANDY;

    SIDIS.Beam.SetHadron(PROTON); //
    SIDIS.Beam.SetPolarization( transverse );

    SIDIS.Target.SetHadron(PROTON); //
    SIDIS.Target.SetPolarization( unp  );

    SIDIS.Produced.SetHadron(PION);
    SIDIS.Produced.SetCharge(positive);
    SIDIS.Produced.SetPolarization( unp );

    SIDIS.SetEnergy( 510.*510. );
//    SIDIS.SetEnergy( 200.*200. );


    SIDIS.Average.y = 3.25; // This is rapidity!

    SIDIS.Limits.Q2_min  = 1.;
    SIDIS.Limits.Q2_max  = 1.e+10;


    SIDIS.SetName("ANDY");

    break;
  case RHIC_PHENIX: // RHIC PP PHENIX


    EXPERIMENT = RHIC_PHENIX;

    SIDIS.Beam.SetHadron(PROTON); //
    SIDIS.Beam.SetPolarization( transverse );

    SIDIS.Target.SetHadron(PROTON); //
    SIDIS.Target.SetPolarization( unp  );

    SIDIS.Produced.SetHadron(PION);
    SIDIS.Produced.SetCharge(positive);
    SIDIS.Produced.SetPolarization( unp );
    //SIDIS.Produced.print( cout );


    SIDIS.SetEnergy( 200.*200. );


    SIDIS.Average.y = 3.5; // This is rapidity!

    SIDIS.Limits.Q2_min  = 1.;
    SIDIS.Limits.Q2_max  = 1.e+10;


    SIDIS.SetName("PHENIX");

    break;
  case RHIC_STAR: // RHIC PP STAR


    EXPERIMENT = RHIC_STAR;

    SIDIS.Beam.SetHadron(PROTON); //
    SIDIS.Beam.SetPolarization( transverse );

    SIDIS.Target.SetHadron(PROTON); //
    SIDIS.Target.SetPolarization( unp   );
    //SIDIS.Target.print( cout );

    SIDIS.Produced.SetHadron(PION);
    SIDIS.Produced.SetCharge(positive);
    SIDIS.Produced.SetPolarization( unp );
    //SIDIS.Produced.print( cout );

    SIDIS.SetEnergy( 200.*200. );
    //SIDIS.SetEnergy( 500.*500. );


    //SIDIS.Average.y = 3.7; // This is rapidity!
    SIDIS.Average.y = 3.; // This is rapidity for JIM!

    //SIDIS.Limits.Q2_min  = 1.; // MAIN
    SIDIS.Limits.Q2_min  = 2.; // NOTE STAR 2020
    SIDIS.Limits.Q2_max  = 1.e+10;

    SIDIS.Limits.xF_min  = 0.;
    SIDIS.Limits.xF_max  = 1.;

    SIDIS.SetName("STAR");

    break;
  case RHIC_BRAHMS: // RHIC PP BRAHMS


    EXPERIMENT = RHIC_BRAHMS;

    SIDIS.Beam.SetHadron(PROTON); //
    SIDIS.Beam.SetPolarization( transverse );

    SIDIS.Target.SetHadron(PROTON); //
    SIDIS.Target.SetPolarization( unp  );
    //SIDIS.Target.print( cout );

    SIDIS.Produced.SetHadron(PION);
    SIDIS.Produced.SetCharge(positive);
    SIDIS.Produced.SetPolarization( unp );
    //SIDIS.Produced.print( cout );

    SIDIS.SetEnergy( 200.*200. );


    //SIDIS.Average.y = 3.7; // This is rapidity!

    SIDIS.Limits.Q2_min  = 1.;
    SIDIS.Limits.Q2_max  = 1.e+10;


    SIDIS.SetName("BRAHMS");

    break;
   default:
    break;
  }


      return 0;

};
