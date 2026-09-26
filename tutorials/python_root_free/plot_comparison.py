"""Plot fresh calculations against vector-extracted Figure 4 central curves."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    r=np.genfromtxt(args.results/'figure4_curves.csv',delimiter=',',names=True)
    ref=np.genfromtxt(args.reference,delimiter=',',names=True)
    meta=json.loads((args.results/'run.json').read_text())
    args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'figure.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex='col',
                          gridspec_kw={'height_ratios':[3,1]},layout='constrained')
    comparisons=[]
    for col,z in enumerate((.13,.37)):
        ax,delta=axes[:,col]
        for charge,color in [(1,'#be342e'),(-1,'#245eb1')]:
            rr=r[(r['z']==z)&(r['charge']==charge)]
            pp=ref[(ref['z']==z)&(ref['charge']==charge)]
            # Last reference x differs from 1.5 by vector rounding only.
            fresh=np.interp(np.log(pp['jt_GeV']),np.log(rr['jt_GeV']),rr['A_UT'])
            residual=fresh-pp['A_UT_paper']
            ax.plot(rr['jt_GeV'],rr['A_UT'],color=color,lw=2)
            ax.plot(pp['jt_GeV'],pp['A_UT_paper'],'o',mfc='white',mec=color,ms=4)
            delta.plot(pp['jt_GeV'],1e3*residual,'o-',color=color,ms=3,lw=1)
            comparisons.append({'z':z,'charge':charge,
                'max_absolute_difference':float(np.max(abs(residual))),
                'rms_absolute_difference':float(np.sqrt(np.mean(residual**2))),
                'max_relative_difference':float(np.max(abs(residual/pp['A_UT_paper']))),
                'max_quadrature_change':float(np.max(abs(rr['A_UT']-rr['A_UT_coarse'])))})
        ax.set_title(rf'$z_h={z}$',loc='left',fontweight='bold')
        ax.axhline(0,color='0.5',lw=.6)
        delta.axhline(0,color='0.5',lw=.6)
        for a in (ax,delta):
            a.set_xscale('log'); a.grid(alpha=.15)
            a.set_xlim(.095,1.6)
        delta.set_xticks([.1,.2,.5,1,1.5],['0.1','0.2','0.5','1','1.5'])
        delta.set_xlabel(r'$j_\perp$ [GeV]')
    axes[0,0].set_ylabel(r'$A_{UT}^{\sin(\phi_S-\phi_H)}$')
    axes[1,0].set_ylabel(r'$10^3\,(A_{\rm run}-A_{\rm paper})$')
    fig.legend(handles=[
        Line2D([],[],color='#be342e',lw=2,label=r'Calculated $\pi^+$'),
        Line2D([],[],color='#245eb1',lw=2,label=r'Calculated $\pi^-$'),
        Line2D([],[],color='0.3',marker='o',mfc='white',ls='',label='Published central curves')],loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Pion in jet: ROOT-free calculation versus Figure 4\n'
                 r'$\sqrt{s}=500$ GeV, $p_T=31$ GeV, $0<\eta<1$; TMD evolution',fontsize=14)
    fig.savefig(args.output/'figure4_comparison.png')
    fig.savefig(args.output/'figure4_comparison.svg')
    (args.output/'comparison.json').write_text(json.dumps({'run':args.results.name,
        'reference':'1707.00913v2, Fig. 4; PDF vector paths, 15 points per curve',
        'configuration':{'PDF':meta['PDF'],'FF':meta['FF'],'matching':meta.get('matching','repository')},
        'curves':comparisons},indent=2)+'\n')
    print(json.dumps(comparisons,indent=2))


if __name__=='__main__': main()
