"""Compare the Python quadrature with the original four-dimensional integrand.

This is deliberately separate from the production algorithm. Scrambled Sobol
replicates integrate y, xb, zhat and b together, using the original C++ matching,
nonperturbative factors, Bessel factors and hard contraction. The two paths share
the PDF/FF providers and HOPPET library; this is not an independent fit or theory.
"""
from pathlib import Path
import argparse
import ctypes as ct
import json
import sys
import time
import numpy as np
from scipy.stats import qmc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pionjet import Backend, asymmetry


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--power',type=int,default=16,help='2**power points per replicate')
    p.add_argument('--replicates',type=int,default=8)
    p.add_argument('--output',type=Path,default=Path('results/audit.json'))
    args=p.parse_args()
    if not (10<=args.power<=22 and 4<=args.replicates<=32):
        p.error('Require 10<=power<=22 and 4<=replicates<=32')
    engine=Backend()
    d=ct.c_double
    arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    engine.lib.pj_reference_config.argtypes=[d,d,d,ct.c_int]
    engine.lib.pj_reference_config.restype=None
    engine.lib.pj_reference_batch.argtypes=[ct.c_int,arr,arr]
    engine.lib.pj_reference_batch.restype=None
    engine.lib.pj_reference_config(engine.gh,engine.c[5],engine.g2,engine.ff_order)
    # Ensure fDSS has opened the requested grid before the batch callback.
    engine.ff(.37,31.)
    start=time.monotonic()
    hard=engine.hard(500.,31.,64)
    cases=[(.13,.1,1),(.37,.5,1),(.37,.5,-1),(.37,1.5,1)]
    results=[]
    for case,(z,jt,charge) in enumerate(cases):
        target=asymmetry(hard,engine.tmd(z,31.,[jt],64),charge)[0][0]
        integrals=[]
        for rep in range(args.replicates):
            unit=qmc.Sobol(4,scramble=True,seed=170700913+rep).random_base2(args.power)
            y=unit[:,0]
            tt=-31*500*np.exp(-y); uu=-31*500*np.exp(y)
            xbmin=-.99*tt/(.99*500**2+uu)
            logrange=np.log(.99/xbmin)
            xb=xbmin*np.exp(unit[:,1]*logrange)
            b=.001*np.exp(unit[:,3]*np.log(10000))
            rows=np.empty((len(unit),9))
            rows[:]=[500,31,0,0,z,0,0,jt,charge]
            rows[:,2]=y; rows[:,3]=xb
            rows[:,5]=z+(1-z)*unit[:,2]; rows[:,6]=b
            out=np.empty((len(unit),2))
            with engine.native_directory():
                engine.lib.pj_reference_batch(len(unit),rows,out)
            jacobian=xb*logrange*(1-z)*b*np.log(10000)
            integrals.append(np.mean(out*jacobian[:,None],axis=0))
        integrals=np.array(integrals)
        ratios=integrals[:,1]/integrals[:,0]
        estimate=float(np.mean(ratios))
        error=float(np.std(ratios,ddof=1)/np.sqrt(args.replicates))
        result={'z':z,'jt_GeV':jt,'charge':charge,'python_gauss':float(target),
                'original_cpp_sobol':estimate,'replicate_standard_error':error,
                'difference_in_standard_errors':float((estimate-target)/error),
                'relative_difference':float(estimate/target-1),
                'replicate_integrals_FUU_FUT':integrals.tolist()}
        # Stochastic check: 5 standard errors plus the separately measured
        # quadrature tolerance. This does not certify the historical paper.
        result['passed']=bool(abs(estimate-target) <= 5*error+1e-6)
        results.append(result)
        print(json.dumps({k:v for k,v in result.items() if not isinstance(v,list)}),flush=True)
    report={'method':'Original C++ NLO integrand; independent scrambled Sobol 4D integration',
            'shared_inputs':'LHAPDF, DSS, DSSV, HOPPET and initial conditions',
            'reference_adaptations':['ROOT Bessel calls mapped to C math',
                                     'uninitialized NP fragmentation slots zeroed',
                                     'sidis state adapter; equal widths and charged pions'],
            'points_per_replicate':2**args.power,'replicates':args.replicates,
            'seed_base':170700913,'elapsed_seconds':time.monotonic()-start,'cases':results}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    if not all(r['passed'] for r in results):
        raise RuntimeError('Original integrand comparison failed; inspect the saved audit.')


if __name__=='__main__': main()
