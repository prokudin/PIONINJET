// Validation adapter, not used by the production quadrature.
// The full original C++ NLO integrand is extracted by prepare.py.
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;
constexpr double PI=3.14159265358979323846, CF=4./3., bmax=1.5;
constexpr double C1=1.1229189671337703, Q0=1.5491933384829668;
inline double pow2(double x) { return x*x; }
double AlphaS(double), PertEvolFeng(double,double);
double UUSIG_JET(double,double,double,double,double,double,double,double,double,double,
                double,double,double,double,double,double,double,
                double,double,double,double,double,double,double);
double CollSIG(double,double,double,double,double,double,double,double,double,
               double,double,double,double,double,double,double,
               double,double,double,double,double,double);
extern "C" void pj_pdf(double,double,double*,double*);
extern "C" void pj_ff(double,double,int,double*);
extern "C" void pj_coll(double,double,double*);
#include "reference_struct.inc"
enum {negative=-1,neutral=0,positive=1,transverse=2,PION=3,HADRONS=4,KAON=5,
      collins_effect=6,no_effect=7,pt_dependence=8};
constexpr bool VALENCE_SIDIS=false;
namespace TMath { // The only ROOT symbols used by this routine.
    double BesselJ0(double x) { return ::j0(x); }
    double BesselJ1(double x) { return ::j1(x); }
}
struct HADRON {
    PARTONCONTENT partcontent{},fragmentation{};
    int charge=1;
    int GetCharge() const { return charge; }
    int GetHadronType() const { return PION; }
    int GetPolarization() const { return transverse; }
};
struct Parameters {
    double a_up=0,b_up=0,n_up=0,gamma_up=0,eta_up=0,m2_up=0;
    double a_down=0,b_down=0,n_down=0,gamma_down=0,eta_down=0,m2_down=0;
};
struct Fit { Parameters parameters; };
static PARTONCONTENT content(const double* f) {
    return {f[0],f[1],f[2],f[3],f[4],f[5],0,0,0,0,0,0,f[6]};
}
static void conjugate(double* f,int charge) {
    if(charge<0) { swap(f[0],f[2]); swap(f[1],f[3]); swap(f[4],f[5]); }
}
class sidis {
public:
    struct { Fit Collins,UnpolarisedFF,UnpolarisedGK; } Params;
    HADRON Beam,Target,Produced;
    int effect=no_effect,dependence=pt_dependence,ff_order=1;
    double cross_section=0,energy_squared=250000;
    static sidis* Instance() { static sidis instance; return &instance; }
    double GetEnergy() const { return energy_squared; }
    void Unpolarised(HADRON& h,double x,double q2) {
        double f[7],t[7]; pj_pdf(x,sqrt(q2),f,t); h.partcontent=content(f);
    }
    void TransversityDistributionNLO(HADRON& h,double x,double q2) {
        double f[7],t[7]; pj_pdf(x,sqrt(q2),f,t); h.partcontent=content(t);
    }
    void Fragmentation(HADRON& h,double z,double q2) {
        double f[7]; pj_ff(z,sqrt(q2),ff_order,f); conjugate(f,h.charge);
        h.fragmentation=content(f);
    }
    void CollinsDistributionNLO(HADRON& h,double z,double q2) {
        double f[7]; pj_coll(z,sqrt(q2),f); conjugate(f,h.charge);
        h.fragmentation=content(f);
    }
    double unpolarisedwidth_x_dependence_Feng(double,double,double,double,double,double);
    PARTONCONTENT unpolarisedff_btpart_nonpert_Feng(HADRON&,double,double);
    PARTONCONTENT collinsff_btpart_nonpert_Feng(HADRON&,double,double);
    int approximate_pionjet_pp_calculation_NLO(double,double,double,double,double,double,double,double);
};
#include "reference_integrand.inc"

extern "C" void pj_reference_config(double gh,double gc,double g2,int ff_order) {
    auto& s=*sidis::Instance();
    s.Params.UnpolarisedFF.parameters.n_up=gh;
    s.Params.UnpolarisedFF.parameters.n_down=gh;
    s.Params.Collins.parameters.m2_up=gc;
    s.Params.Collins.parameters.m2_down=gc;
    s.Params.UnpolarisedGK.parameters.a_up=g2;
    s.ff_order=ff_order;
}
// rows: energy,pt,y,xb,z,zhat,b,jt,charge; output: FUU,FUT in driver convention.
extern "C" void pj_reference_batch(int n,const double* rows,double* output) {
    auto& s=*sidis::Instance();
    for(int i=0;i<n;++i) {
        const double* r=rows+9*i;
        s.energy_squared=r[0]*r[0]; s.Produced.charge=int(r[8]);
        for(int k=0;k<2;++k) {
            s.effect=k?collins_effect:no_effect;
            s.approximate_pionjet_pp_calculation_NLO(r[3],r[2],r[4],r[5],r[6],r[1],r[7],r[1]*r[1]);
            output[2*i+k]=s.cross_section;
        }
    }
}
