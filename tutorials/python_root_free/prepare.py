"""Extract the unchanged numerical kernels; leave the upstream checkout untouched."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import re

UPSTREAM = '3945b91086dd153b4d940531e2f09f8f70c0496f'
PARAMETERS = ('transversity_parametersKANG2015.dat', 'collins_parametersKANG2015.dat',
              'unpolarisedff_parameters.dat', 'unpolarisedgk_parameters.dat')


def function(text, signature):
    # Ignore comments and literals while locating braces. Keep the returned
    # source byte-for-byte, including comments, for auditability.
    start = text.index(signature)
    masked = re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
                    lambda m: ' '*len(m.group()), text, flags=re.S)
    body = masked.index('{', start)
    depth = 1
    end = body + 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return text[start:end] + '\n'


def prepare(repo, build):
    repo, build = repo.resolve(), build.resolve()
    build.mkdir(parents=True, exist_ok=True)
    physics = (repo / 'pionjet_pp.cpp').read_text()
    evolution = (repo / 'evolution_common.cpp').read_text()
    kernels = '#include <cmath>\nusing namespace std;\n'
    kernels += 'constexpr double PI=3.14159265358979323846, CF=4./3., nf=5., LamQCD=.225;\n'
    kernels += 'constexpr double C1=1.1229189671337703, bmax=1.5, fourpi=4*PI;\n'
    kernels += 'inline double pow2(double x){return x*x;}\n'
    for sig in ('double AlphaS(', 'double Sudakovkernel(', 'double PertEvolFeng('):
        kernels += function(evolution, sig)
    for sig in ('double alphasQCD_JET(', 'int Upp_JET(', 'double UUSIG_JET(', 'int Collpp(', 'double CollSIG('):
        kernels += function(physics, sig)
    (build / 'kernels.inc').write_text(kernels)
    hoppet = build / 'hoppet'
    hoppet.mkdir(exist_ok=True)
    for src in (repo / 'hoppet-1.1.5/src').iterdir():
        if src.suffix.lower() in ('.f', '.f90') or src.name == 'Makefile':
            shutil.copy2(src, hoppet / src.name)
    for src in ('FORTRAN/frag/fDSS.f', 'pdf_polarized/DSSV-2008.f'):
        shutil.copy2(repo / src, build / Path(src).name)
    # Use precisely the tables found by the original driver from the repo root.
    # PINLO.GRID contains three identical 816-record blocks; fDSS reads the first.
    for filename in ('PINLO.GRID', 'PILO.GRID', 'DSSV-GRID.NLO'):
        shutil.copy2(repo / filename, build / filename)
    (build/'parameters').mkdir(exist_ok=True)
    for name in PARAMETERS:
        shutil.copy2(repo/'parameters'/name, build/'parameters'/name)
    # Validation only: compile the complete original NLO integrand against a
    # small sidis adapter. The adapter supplies the same PDF/FF providers.
    hadron = (repo/'hadron.h').read_text()
    reference = function(hadron, 'struct PARTONCONTENT {').rstrip()+';\n'
    for sig in ('inline PARTONCONTENT product(PARTONCONTENT',
                'inline PARTONCONTENT product(double',
                'inline PARTONCONTENT sum(double',
                'inline PARTONCONTENT sum(PARTONCONTENT'):
        reference += function(hadron, sig)
    (build/'reference_struct.inc').write_text(reference)
    reference = function((repo/'twist3.cpp').read_text(), 'double alphasQCD(')
    for sig in ('double Pqq(', 'double Pgq(', 'double deltaPqq(',
                'double sidis::unpolarisedwidth_x_dependence_Feng(',
                'double NonPertEvolFeng(',
                'PARTONCONTENT  sidis::unpolarisedff_btpart_nonpert_Feng(',
                'PARTONCONTENT  sidis::collinsff_btpart_nonpert_Feng(',
                'int sidis::approximate_pionjet_pp_calculation_NLO('):
        reference += function(physics, sig)
    # The source left the NP gluon slot indeterminate. Zero it explicitly, as
    # in the production wrapper; do not compare against undefined behavior.
    assert reference.count('PARTONCONTENT fragmentation;') == 2
    reference = reference.replace('PARTONCONTENT fragmentation;', 'PARTONCONTENT fragmentation{};')
    (build/'reference_integrand.inc').write_text(reference)
    inputs = ['pionjet_pp.cpp', 'evolution_common.cpp', 'PINLO.GRID', 'PILO.GRID',
              'DSSV-GRID.NLO', 'FORTRAN/frag/fDSS.f', 'pdf_polarized/DSSV-2008.f',
              'sidis.cpp', 'hadron.h', 'twist3.cpp', 'hoppet-1.1.5/src/Makefile']
    manifest = {s: hashlib.sha256((repo/s).read_bytes()).hexdigest() for s in inputs}
    for s in (repo/'hoppet-1.1.5/src').iterdir():
        if s.suffix.lower() in ('.f','.f90'):
            manifest[str(s.relative_to(repo))]=hashlib.sha256(s.read_bytes()).hexdigest()
    (build / 'source_hashes.json').write_text(json.dumps(manifest, indent=2)+'\n')
    # Check the numerical input files against the audited upstream commit,
    # even when HEAD is a later tutorial commit. Never silently label a
    # modified source tree as the audited revision.
    for name, digest in manifest.items():
        original = subprocess.check_output(['git','-C',str(repo),'show',f'{UPSTREAM}:{name}'])
        if hashlib.sha256(original).hexdigest() != digest:
            raise ValueError(f'Input differs from audited upstream {UPSTREAM}: {name}')
    (build/'upstream_revision.txt').write_text(UPSTREAM+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--build', type=Path, default=Path(__file__).resolve().parent / 'build')
    args = parser.parse_args()
    prepare(args.repo, args.build)
