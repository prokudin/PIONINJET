"""Check Figure 5 in both the driver and Cartesian differential conventions.

The driver's integrands multiply by (2 pt)(2 jt). Eq. (1) uses d^2pt d^2jt.
Store both, without fitting a normalization to the published curve.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pionjet import Backend

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reference',type=Path,required=True)
p.add_argument('--evolved',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
args=p.parse_args()
ref=np.genfromtxt(args.reference.resolve(),delimiter=',',names=True)
evolved=np.genfromtxt(args.evolved.resolve()/'figure4_curves.csv',delimiter=',',names=True)
output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
engine=Backend(pdfset='cteq6l1',ff_order=0)
hard=engine.hard(500.,31.,64)[0]
ff=engine.ff(.37,31.)
jt=np.linspace(.1,1.5,200)
alpha2=engine.lib.pj_alpha_hard(31.)**2
gaussian=np.exp(-jt**2/.12)/(np.pi*.12)
no_gluon=2*jt*alpha2*(hard[:6]@ff[:6])*gaussian
with_gluon=2*jt*alpha2*(hard@ff)*gaussian
er=evolved[(evolved['z']==.37)&(evolved['charge']==1)]
ev=2*er['jt_GeV']*alpha2*er['FUU_reduced']
np.savetxt(output/'figure5_gaussian.csv',np.column_stack([jt,no_gluon,with_gluon,
           no_gluon/(4*31*jt),with_gluon/(4*31*jt)]),delimiter=',',comments='',
           header='jt,driver_no_gluon,driver_with_gluon,cartesian_no_gluon,cartesian_with_gluon')
np.savetxt(output/'figure5_evolved.csv',np.column_stack([er['jt_GeV'],ev,ev/(4*31*er['jt_GeV'])]),
           delimiter=',',comments='',header='jt,driver_no_gluon,cartesian_no_gluon')
stats=[]
for evolution,label,xx,yy in [(0,'Gaussian, quark jets',jt,no_gluon),
                              (0,'Gaussian, quark and gluon jets',jt,with_gluon),
                              (1,'TMD evolution, quark jets',er['jt_GeV'],ev)]:
    rr=ref[ref['evolution']==evolution]
    if evolution:
        interp=np.exp(np.interp(rr['jt_GeV'],xx,np.log(yy)))
    else:
        amplitude=yy[0]/(xx[0]*np.exp(-xx[0]**2/.12))
        interp=amplitude*rr['jt_GeV']*np.exp(-rr['jt_GeV']**2/.12)
    ratio=interp/rr['cross_section_paper']
    stats.append({'evolution':evolution,'calculation':label,
                  'driver_to_paper_ratio_min':float(min(ratio)),
                  'driver_to_paper_ratio_max':float(max(ratio)),
                  'max_relative_difference':float(max(abs(ratio-1)))})
print(json.dumps(stats,indent=2))
(output/'figure5_checks.json').write_text(json.dumps({'comparison':stats,
    'alpha_hard_squared':alpha2,'LO_pdf':'cteq6l1','LO_ff':'DSS LO, PILO.GRID',
    'gaussian_width_GeV2':.12,
    'convention':'Driver includes 4*pt*jt relative to Cartesian FUU of paper Eq. 1.',
    'outgoing_gluon_fraction_LO':float(hard[6]*ff[6]/(hard@ff))},indent=2)+'\n')
fig,ax=plt.subplots(figsize=(7,5),layout='constrained')
ax.semilogy(jt,no_gluon,color='#245eb1',label='Gaussian, quark jets only')
ax.semilogy(jt,with_gluon,color='#245eb1',ls=':',label='Gaussian, including gluon jets')
ax.semilogy(er['jt_GeV'],ev,color='#be342e',label='TMD evolution, quark jets only')
for e,c in [(0,'#245eb1'),(1,'#be342e')]:
    rr=ref[ref['evolution']==e]
    ax.plot(rr['jt_GeV'],rr['cross_section_paper'],'o',color=c,mfc='white',ms=4)
ax.plot([],[],'o',color='0.3',mfc='white',label='Published curves')
ax.set(xlabel=r'$j_\perp$ [GeV]',ylabel='Cross section in repository driver convention',
       title=r'Figure 5 check: $\pi^+$, $z_h=0.37$, $p_T=31$ GeV',
       xlim=(.08,1.55),ylim=(1.e-13,1.e-5))
ax.grid(alpha=.15);ax.legend(frameon=False,loc='lower left',fontsize=9)
fig.savefig(output/'figure5_comparison.png',dpi=180)
fig.savefig(output/'figure5_comparison.svg')
