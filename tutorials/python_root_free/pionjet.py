"""ROOT-free central-fit pp -> (jet pi) X calculation for arXiv:1707.00913.

Python handles integration, configuration, validation and output. The tiny native
backend retains the supplied hard factors, Sudakov expression, DSS/DSSV and
custom transversity HOPPET. No ROOT, Cuba, Pythia, fitting framework or cfortran.
All scales/momenta are in GeV; b is in GeV^-1.
"""
from pathlib import Path
import argparse
import ctypes as ct
import hashlib
import json
import os
import time
import platform
from contextlib import contextmanager
import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.special import j0, j1
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
CF = 4 / 3
C1 = 2 * np.exp(-np.euler_gamma)
Q0 = np.sqrt(2.4)
BMAX = 1.5


def gauss(n, lo, hi):
    x, w = leggauss(n)
    return (hi+lo)/2 + (hi-lo)*x/2, (hi-lo)*w/2


def parameters(path):
    v = np.loadtxt(path)
    if v.shape != (36, 2) or not np.isfinite(v).all():
        raise ValueError(f'Invalid fit parameter file: {path}')
    return np.ascontiguousarray(v[:, 0])


class Backend:
    _initialized = False

    def __init__(self, build=HERE/'build', dy=.1, pdfset='CT10nlo', ff_order=1):
        if Backend._initialized:
            raise RuntimeError('HOPPET has global state: use a separate process for another Backend.')
        if ff_order not in (0,1) or not 0 < dy <= .2:
            raise ValueError('Require ff_order=0 or 1 and 0<dy<=0.2')
        self.build = Path(build).resolve()
        self.ff_order=ff_order
        for name in ('PINLO.GRID', 'PILO.GRID', 'DSSV-GRID.NLO', 'libpionjet.so'):
            if not (self.build/name).is_file():
                raise FileNotFoundError(self.build/name)
        self.t = parameters(self.build/'parameters/transversity_parametersKANG2015.dat')
        self.c = parameters(self.build/'parameters/collins_parametersKANG2015.dat')
        self.gh = parameters(self.build/'parameters/unpolarisedff_parameters.dat')[2]
        self.g2 = parameters(self.build/'parameters/unpolarisedgk_parameters.dat')[0]
        widths = parameters(self.build/'parameters/unpolarisedff_parameters.dat').reshape(6,6)
        if np.any(widths[:,[0,1,3,4]] != 0) or np.any(widths[:,2] != self.gh):
            raise ValueError('This tutorial implements a constant common unpolarized width.')
        if self.c[5] != self.c[11] or not 0 < self.c[5] < self.gh:
            raise ValueError('Require equal Collins widths with 0<gc<gh.')
        self.lib = ct.CDLL(str(self.build/'libpionjet.so'))
        arr = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        d = ct.c_double
        signatures = {
            'pj_init': ([arr, arr, d, ct.c_char_p], ct.c_int),
            'pj_ff': ([d,d,ct.c_int,arr], None),
            'pj_coll': ([d,d,arr], None),
            'pj_pdf': ([d,d,arr,arr], None),
            'pj_hard': ([d,d,d,d,arr,arr], None),
            'pj_direct': ([d,d,d,d,arr,arr,arr], None),
            'pj_initial': ([d,d,arr], None), 'pj_evolved': ([d,d,arr], None),
            'pj_sudakov': ([d,d], d), 'pj_alpha': ([d], d),
            'pj_alpha_hard': ([d],d),
            'pj_j0': ([d], d), 'pj_j1': ([d], d),
        }
        for name, (args, ret) in signatures.items():
            fun = getattr(self.lib, name)
            fun.argtypes, fun.restype = args, ret
        # The original Fortran routines open tables relative to the current cwd.
        # This backend is process-global and must not be used concurrently in threads.
        with self.native_directory():
            if self.lib.pj_init(self.t, self.c, dy, pdfset.encode()):
                raise RuntimeError('Native physics initialization failed')
        Backend._initialized = True

    @contextmanager
    def native_directory(self):
        previous = Path.cwd()
        try:
            os.chdir(self.build)
            yield
        finally:
            os.chdir(previous)

    def ff(self,z,q):
        f=np.empty(7)
        with self.native_directory():
            self.lib.pj_ff(z,q,self.ff_order,f)
        return f

    def coll(self,z,q):
        f=np.empty(7)
        self.lib.pj_coll(z,q,f)
        return f

    def hard(self,energy,pt,n=32):
        """Integrate y and xb once per jet kinematics; keep outgoing flavors.

        The original integrand is linear in outgoing FFs, independent of xa/xb/y.
        Thus this factorization is algebraically exact at fixed jet pt and z.
        """
        result=np.zeros((2,7))
        if not (np.isfinite(energy) and np.isfinite(pt) and energy>2*pt and .75<=pt<=1e5):
            raise ValueError('Require energy>2*pt and .75<=pt<=1e5 (GeV).')
        ys,wy=gauss(n,0,1)
        for y,w_y in zip(ys,wy):
            ss=energy**2; tt=-energy*pt*np.exp(-y); uu=-energy*pt*np.exp(y)
            # Solve xa < .99 exactly (the original numerical endpoint cut).
            if .99*ss+uu <= 0: continue
            xbmin=-.99*tt/(.99*ss+uu)
            if xbmin>=.99: continue
            logs,wx=gauss(3*n,np.log(xbmin),np.log(.99))
            for logx,w_x in zip(logs,wx):
                xb=np.exp(logx)
                pair=np.empty((2,7))
                self.lib.pj_hard(energy,pt,y,xb,pair[0],pair[1])
                result+=pair*(w_y*w_x*xb)
        return result

    def tmd(self,z,pt,jt,n=32,bmin=.001,bmax=10,matching='repository'):
        """NLL TMD evolution and matching in the original conventions.

        Outgoing gluon TMD is excluded as specified in the paper. Gluon mixing
        in quark matching coefficients is retained (Eq. 10).
        """
        jt=np.atleast_1d(jt)
        if matching not in ('repository','unit-delta','tree'):
            raise ValueError('Unknown matching option')
        if not (.01<=z<.95 and C1/1e5<=bmin<bmax and n>=8
                and np.isfinite(jt).all() and np.all(jt>=0)):
            raise ValueError('TMD arguments outside the supported numerical domain')
        zs,wz=gauss(2*n,z,1)
        zz=z/zs
        # log(z/z')/(1-z/z') is regular at the upper endpoint, which Gauss avoids.
        logz=np.log(zz)
        qq=CF/2*(1-zz)+CF*(1+zz**2)/(1-zz)*logz
        gq=CF/2*zz+CF*(1+(1-zz)**2)/zz*logz
        hq=CF*2*zz/(1-zz)*logz
        boundaries=sorted(set([bmin,bmax]+[b for b in (.01,.05,.1,.2,.5,1,2,4,7) if bmin<b<bmax]))
        bs,ws=zip(*(gauss(n,a,b) for a,b in zip(boundaries[:-1],boundaries[1:])))
        bs,ws=np.concatenate(bs),np.concatenate(ws)
        out=np.zeros((2,len(jt),7))
        for b,wb in zip(bs,ws):
            bstar=b/np.sqrt(1+(b/BMAX)**2)
            mu=C1/bstar; a=self.lib.pj_alpha(mu)/np.pi
            if matching=='tree': a=0.
            fz,cz=self.ff(z,mu),self.coll(z,mu)
            fhat=np.array([self.ff(zh,mu) for zh in zs])
            chat=np.array([self.coll(zh,mu) for zh in zs])
            # Equations (8),(10) in 1707.00913 have a unit delta coefficient.
            # The uploaded source instead has 1-2*CF*alpha_s/pi in both.
            delta=1-2*CF*a if matching=='repository' else 1.
            f=delta*fz + a*np.sum((wz*qq/zs)[:,None]*fhat,axis=0)
            f+=a*np.dot(wz*gq/zs,fhat[:,6])
            c=delta*cz + a*np.sum((wz*hq/zs)[:,None]*chat,axis=0)
            # f_g is not a gluon TMD matching calculation; exclude it explicitly.
            f[6]=0.; c[6]=0.
            evol=self.lib.pj_sudakov(pt,b)*np.exp(max(-20.,-self.g2/2*np.log(b/bstar)*np.log(pt/Q0)))
            f*=evol*np.exp(-self.gh*b*b/(z*z))
            # KANG2015 has equal g_c for favored/unfavored. Assert rather than
            # silently generalizing this simplified scalar-width implementation.
            c*=evol*np.exp(-(self.gh-self.c[5])*b*b/(z*z))
            out[0]+=wb*b/(2*np.pi*z*z)*j0(b*jt[:,None]/z)*f
            # Original relation_to_trento=-z, times b^2/(2*pi*z^3).
            out[1]+=-wb*b*b/(2*np.pi*z*z)*j1(b*jt[:,None]/z)*c
        if not np.isfinite(out).all():
            raise FloatingPointError('Nonfinite TMD integral')
        return out

    def validate(self):
        checks={}
        max_rel=0
        input_errors=[]
        for x in (.01,.05,.13,.2,.37,.7):
            initial,evolved=np.empty(13),np.empty(13)
            with self.native_directory():
                self.lib.pj_initial(x,Q0,initial)
            self.lib.pj_evolved(x,Q0,evolved)
            active=np.abs(initial)>1.e-9
            max_rel=max(max_rel,float(np.max(np.abs(evolved[active]/initial[active]-1))))
            input_errors.append({'x':x,'initial':initial.tolist(),'evolved':evolved.tolist()})
        checks['hoppet_input_scale_max_relative_error']=max_rel
        checks['hoppet_input_scale_samples']=input_errors
        # Interpolation of the tabulated input functions need not be exact at Q0.
        # Require sub-percent closure and record the full residuals for scrutiny.
        assert max_rel<.01, checks
        # Independently integrate the NLL Sudakov in log(mu). The extracted
        # antiderivative uses d(mu^2)/mu^2 and then an explicit factor 1/2.
        residual=[]
        for b in (.03,.2,1.,3.):
            mu=C1*np.sqrt(1+(b/BMAX)**2)/b
            def kernel(logmu):
                a=self.lib.pj_alpha(np.exp(logmu))/np.pi
                A=CF*a+CF*((67/36-np.pi**2/12)*3-5/18*5)*a*a
                B=-1.5*CF*a
                return A*np.log(31**2/np.exp(2*logmu))+B
            analytic=self.lib.pj_sudakov(31,b)
            numerical=np.exp(-quad(kernel,np.log(mu),np.log(31),epsabs=1e-11)[0])
            residual.append(abs(analytic/numerical-1))
        checks['sudakov_max_relative_error']=max(residual)
        assert max(residual)<1e-9, checks
        rng=np.random.default_rng(170700913)
        residual=[]
        for _ in range(12):
            y=float(rng.uniform(0,1)); xb=float(rng.uniform(.2,.8))
            weights=np.empty((2,7)); f=rng.uniform(.1,1,7); c=rng.uniform(-1,1,7)
            self.lib.pj_hard(500,31,y,xb,weights[0],weights[1])
            direct=np.empty(2)
            self.lib.pj_direct(500,31,y,xb,f,c,direct)
            contracted=np.array([weights[0]@f,weights[1]@c])
            residual.append(float(np.max(np.abs(contracted/direct-1))))
        checks['hard_factorization_max_relative_error']=max(residual)
        assert max(residual)<1e-12, checks
        checks['bessel_max_absolute_error']=max(
            abs(getattr(self.lib,'pj_'+name)(x)-fun(x))
            for name,fun in [('j0',j0),('j1',j1)] for x in (.001,.1,1.,3.,10.,80.))
        return checks


def asymmetry(hard,tmd,charge):
    # Charge conjugation; s=sbar for both pion charges.
    order=np.array([0,1,2,3,4,5,6] if charge==1 else [2,3,0,1,5,4,6])
    uu=tmd[0][:,order]@hard[0]
    ut=tmd[1][:,order]@hard[1]
    if np.any(uu<=0): raise FloatingPointError('Nonpositive denominator')
    return ut/uu,uu,ut


def run(args):
    start=time.monotonic()
    output=args.output.resolve(); output.mkdir(parents=True,exist_ok=True)
    engine=Backend(args.build,args.dy,args.pdfset,args.ff_order)
    validation=engine.validate()
    print('Validation:',{k:v for k,v in validation.items() if not isinstance(v,list)},flush=True)
    jt=np.geomspace(.1,1.5,65)
    hard=engine.hard(args.energy,args.pt,args.order)
    hcheck=engine.hard(args.energy,args.pt,args.order*2)
    records=[]; comparisons=[]
    for z in args.z:
        tmd=engine.tmd(z,args.pt,jt,args.order,bmin=args.bmin,bmax=args.bmax,matching=args.matching)
        refined=engine.tmd(z,args.pt,jt,args.order*2,bmin=args.bmin,bmax=args.bmax,matching=args.matching)
        for charge in (1,-1):
            a,uu,ut=asymmetry(hcheck,refined,charge)
            coarse=asymmetry(hard,tmd,charge)[0]
            comparisons.append({'z':z,'charge':charge,'max_absolute_A_change':float(np.max(abs(a-coarse))),
                                'max_relative_A_change':float(np.max(abs(a/coarse-1)))})
            for i,j in enumerate(jt):
                records.append((z,charge,j,a[i],coarse[i],uu[i],ut[i]))
            print(f'z={z} charge={charge:+}: A at j=0.1,~0.5,1.5: {a[[0,38,-1]]}',flush=True)
    np.savetxt(output/'figure4_curves.csv',records,delimiter=',',
               header='z,charge,jt_GeV,A_UT,A_UT_coarse,FUU_reduced,FUT_reduced',comments='')
    validation['quadrature_refinement']=comparisons
    metadata=vars(args).copy()
    # Published metadata contains logical folder names, not workstation paths.
    metadata={k:v.name if isinstance(v,Path) else v for k,v in metadata.items()}
    metadata.update({'fit':'KANG2015 (central)','PDF':args.pdfset+' member 0',
        'FF':('DSS NLO; repository root PINLO.GRID, first 816 records' if args.ff_order else 'DSS LO; repository root PILO.GRID'),
        'Collins_initial_FF':'DSS NLO, independently of unpolarized FF choice',
        'dy':args.dy,'outgoing_gluon_TMD':False,'gluon_to_quark_matching':True,
        'rapidity':[0,1],'x_cut':.99,'b_range':[args.bmin,args.bmax],
        'reduced_cross_section':'Omits common alpha_s(pt)^2 and 2*jt; do not compare absolute normalization to Fig. 5.',
        'elapsed_seconds':time.monotonic()-start,
        'upstream_revision':(engine.build/'upstream_revision.txt').read_text().strip(),
        'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'wrapper_hashes':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                          for name in ['pionjet.py','bridge.cpp','prepare.py','build.sh','validation/reference.cpp']},
        'validation':validation,'source_hashes':json.loads((engine.build/'source_hashes.json').read_text()),
        'input_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
             [engine.build/'parameters'/n for n in ['collins_parametersKANG2015.dat',
              'transversity_parametersKANG2015.dat','unpolarisedff_parameters.dat','unpolarisedgk_parameters.dat']]}})
    (output/'run.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Finished in {metadata["elapsed_seconds"]:.1f} s: {output}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,default=HERE/'build')
    p.add_argument('--output',type=Path,default=HERE/'results')
    p.add_argument('--order',type=int,default=24)
    p.add_argument('--dy',type=float,default=.1)
    p.add_argument('--pdfset',default='CT10nlo')
    p.add_argument('--matching',choices=['unit-delta','repository','tree'],default='repository',
                   help='repository retains its finite delta term; unit-delta follows Eqs. 8,10; tree drops one-loop matching')
    p.add_argument('--ff-order',type=int,choices=[0,1],default=1)
    p.add_argument('--energy',type=float,default=500.)
    p.add_argument('--pt',type=float,default=31.)
    p.add_argument('--z',type=float,nargs='+',default=[.13,.37])
    p.add_argument('--bmin',type=float,default=.001)
    p.add_argument('--bmax',type=float,default=10.)
    args=p.parse_args()
    if not (args.order>=8 and 0<args.dy<=.2 and args.energy>2*args.pt>0
            and .75<=args.pt<=1e5 and all(.01<=z<.95 for z in args.z)
            and C1/1e5<=args.bmin<args.bmax):
        p.error('Require order>=8, 0<dy<=0.2, energy>2*pt, .75<=pt<=1e5, .01<=z<.95, C1/1e5<=bmin<bmax')
    run(args)
