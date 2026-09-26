// Minimal C ABI for Python. Numerical hard factors/Sudakov are extracted verbatim.
#include <LHAPDF/LHAPDF.h>
#include <algorithm>
#include <memory>
#include "kernels.inc"

extern "C" {
void fdss_(int*,int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
void dssvini_(int*);
void dssvfit_(double*,double*,double*,double*,double*,double*,double*,double*);
void hoppetstartextended_(double*,double*,double*,double*,double*,int*,int*,int*);
void hoppetevolve_(double*,double*,int*,double*,void (*)(const double&,const double&,double*),double*);
void hoppeteval_(const double*,const double*,double*);
extern struct {int initialized;} fragini_deflorian_;
}

static std::unique_ptr<LHAPDF::PDF> proton;
static double trans[12], coll[12];
static int current_ff_order=-1;

// Flavor order throughout the Python API: u,d,ubar,dbar,s,sbar,g.
extern "C" void pj_ff(double z, double q, int order, double* f) {
    if(current_ff_order!=order) {
        fragini_deflorian_.initialized=0;
        current_ff_order=order;
    }
    int ih=1, charge=1;
    double c,b,q2=q*q;
    fdss_(&ih,&charge,&order,&z,&q2,f,f+2,f+1,f+3,f+4,f+5,&c,&b,f+6);
    for(int k=0;k<7;++k) f[k]/=z;
}

double shape(double x,const double* p) {
    const double a=p[0],b=p[1],n=p[2];
    double norm=1.;
    if(a!=0 && b!=0) norm=std::pow(a+b,a+b)/(std::pow(a,a)*std::pow(b,b));
    return n*std::pow(x,a)*std::pow(1-x,b)*norm;
}

void initial(const double& x,const double& q,double* f) {
    std::fill(f,f+13,0.); // Original callback omitted the gluon slot.
    double xx=std::max(x,1.01e-5),q2=q*q,uv,dv,ub,db,s,g,ff[7];
    dssvfit_(&xx,&q2,&uv,&dv,&ub,&db,&s,&g);
    f[8]=.5*shape(x,trans)*(proton->xfxQ(2,x,q)+uv+ub);
    f[4]=.5*shape(x,trans+6)*(proton->xfxQ(1,x,q)+dv+db);
    pj_ff(x,q,1,ff);
    f[7]=coll[2]*std::pow(x,coll[0])*std::pow(1-x,coll[1])*x*ff[0];
    f[5]=coll[8]*std::pow(x,coll[6])*std::pow(1-x,coll[7])*x*ff[2];
    f[9]=coll[8]*std::pow(x,coll[6])*std::pow(1-x,coll[7])*x*ff[4];
}

extern "C" int pj_init(const double* t,const double* c,double dy,const char* pdfset) {
    try {
        proton.reset(LHAPDF::mkPDF(pdfset,0));
        std::copy(t,t+12,trans); std::copy(c,c+12,coll);
        int set=0; dssvini_(&set);
        double ymax=12,qmin=.5,qmax=1.e5,dln=dy/4.;
        int loops=1,order=-6,scheme=5;
        hoppetstartextended_(&ymax,&dy,&qmin,&qmax,&dln,&loops,&order,&scheme);
        double q0=std::sqrt(2.4),as=AlphaS(q0),ratio=1.;
        hoppetevolve_(&as,&q0,&loops,&ratio,initial,&q0);
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<std::endl; return 1; }
}

extern "C" void pj_pdf(double x,double q,double* f,double* h) {
    const int ids[7]={2,1,-2,-1,3,-3,21};
    for(int k=0;k<7;++k) f[k]=proton->xfxQ(ids[k],x,q)/x;
    double v[13]; hoppeteval_(&x,&q,v);
    std::fill(h,h+7,0.);
    h[0]=v[8]/x; h[1]=v[4]/x; h[2]=v[10]/x; h[3]=v[2]/x;
}

extern "C" void pj_coll(double z,double q,double* f) {
    double v[13]; hoppeteval_(&z,&q,v);
    f[0]=v[7]/z; f[1]=v[5]/z; f[2]=v[5]/z; f[3]=v[7]/z;
    f[4]=v[9]/z; f[5]=v[9]/z; f[6]=0.;
}

double unpolarized(double s,double t,double u,const double* a,const double* b,const double* f) {
    return UUSIG_JET(s,t,u,a[0],a[1],a[2],a[3],a[4],a[5],a[6],
        b[0],b[1],b[2],b[3],b[4],b[5],b[6],f[0],f[2],f[1],f[3],f[4],f[5],f[6]);
}
double polarized(double s,double t,double u,const double* a,const double* b,const double* f) {
    return CollSIG(s,t,u,a[0],a[1],a[2],a[3],a[4],a[5],
        b[0],b[1],b[2],b[3],b[4],b[5],b[6],f[0],f[2],f[1],f[3],f[4],f[5]);
}

extern "C" void pj_hard(double energy,double pt,double y,double xb,double* uu,double* ut) {
    std::fill(uu,uu+7,0.); std::fill(ut,ut+7,0.);
    double ss=energy*energy,tt=-pt*energy*std::exp(-y),uu0=-pt*energy*std::exp(y);
    double xa=-xb*uu0/(xb*ss+tt);
    if(xa<=0 || xa>=.99 || xb<=0 || xb>=.99) return;
    double a[7],b[7],h[7],hb[7];
    pj_pdf(xa,pt,a,h); pj_pdf(xb,pt,b,hb);
    double pref=2*pt/(xb*ss+tt)/(xb*xa)/ss;
    // Common alpha_s(pt)^2 cancels in A_UT; Python reports reduced F_UU.
    for(int k=0;k<7;++k) {
        double f[7]={}; f[k]=1.;
        uu[k]=pref*unpolarized(xa*xb*ss,xa*tt,xb*uu0,a,b,f);
        ut[k]=pref*polarized(xa*xb*ss,xa*tt,xb*uu0,h,b,f);
    }
}
extern "C" double pj_alpha(double q) {return AlphaS(q);}
extern "C" double pj_alpha_hard(double q) {return alphasQCD_JET(q);}
extern "C" double pj_sudakov(double q,double b) {return PertEvolFeng(q,b);}
extern "C" double pj_j0(double x) {return ::j0(x);}
extern "C" double pj_j1(double x) {return ::j1(x);}
extern "C" void pj_initial(double x,double q,double* f) {initial(x,q,f);}
extern "C" void pj_evolved(double x,double q,double* f) {hoppeteval_(&x,&q,f);}
extern "C" void pj_direct(double energy,double pt,double y,double xb,
                           const double* d,const double* c,double* out) {
    double ss=energy*energy,tt=-pt*energy*std::exp(-y),uu=-pt*energy*std::exp(y);
    double xa=-xb*uu/(xb*ss+tt),a[7],b[7],h[7],hb[7];
    pj_pdf(xa,pt,a,h); pj_pdf(xb,pt,b,hb);
    double pref=2*pt/(xb*ss+tt)/(xb*xa)/ss;
    out[0]=pref*unpolarized(xa*xb*ss,xa*tt,xb*uu,a,b,d);
    out[1]=pref*polarized(xa*xb*ss,xa*tt,xb*uu,h,b,c);
}
