"""Run the tutorial from any working directory; generated files go to results/."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit',action='store_true',help='Also integrate the original C++ routine independently')
    p.add_argument('--convergence',action='store_true',help='Also vary the HOPPET spacing and b limits')
    args=p.parse_args()
    if not (HERE/'build/libpionjet.so').exists():
        p.error('Build first: python3 prepare.py, then bash build.sh (see README.md)')

    def run(*parts):
        subprocess.run([sys.executable,*parts],cwd=HERE,check=True)

    run('pionjet.py','--order','32','--output','results/central')
    run('plot_comparison.py','--results','results/central',
        '--reference','reference/paper_figure4_vector.csv','--output','results/comparison')
    run('figure5.py','--evolved','results/central',
        '--reference','reference/paper_figure5_vector.csv','--output','results/figure5')
    if args.convergence:
        run('pionjet.py','--order','32','--dy','.05','--output','results/fine_hoppet')
        run('pionjet.py','--order','32','--bmin','.0001','--bmax','15','--output','results/b_range')
        base=np.genfromtxt(HERE/'results/central/figure4_curves.csv',delimiter=',',names=True)
        changes={}
        for name in ('fine_hoppet','b_range'):
            other=np.genfromtxt(HERE/f'results/{name}/figure4_curves.csv',delimiter=',',names=True)
            changes[name]={'max_absolute_A_change':float(np.max(abs(other['A_UT']-base['A_UT']))),
                           'max_relative_A_change':float(np.max(abs(other['A_UT']/base['A_UT']-1)))}
        (HERE/'results/convergence.json').write_text(json.dumps(changes,indent=2)+'\n')
        print(json.dumps(changes,indent=2),flush=True)
    if args.audit:
        run('validation/audit.py','--output','results/audit.json')


if __name__=='__main__': main()
