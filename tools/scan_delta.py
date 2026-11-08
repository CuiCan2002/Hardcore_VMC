#!/usr/bin/env python3
"""Prepare/submit fixed-delta runs; collect final energies and long-range correlations."""
import argparse
import csv
import datetime
import math
import os
from pathlib import Path
import subprocess


def read_input(path):
    values = {}
    for line in path.read_text().splitlines():
        parts = line.split('#', 1)[0].replace('=', ' ').split()
        if not parts:
            continue
        if len(parts) != 2 or parts[0] in values:
            raise ValueError('Invalid/duplicate input: ' + line)
        values[parts[0]] = parts[1]
    if 'lambda_eff' in values or 'stability_margin' in values:
        raise ValueError('Migrate old input to delta first; lambda_eff is now derived.')
    return values


def collect(manifest):
    rows = list(csv.DictReader(manifest.open()))
    report = manifest.parent / 'scan_summary.csv'
    fields = ['delta_requested', 'job_id', 'status', 'run_dir', 'delta', 'chi1', 'Delta',
              'lambda_eff', 'min_Ek', 'E_per_site', 'blocked_error_per_site',
              'sr_iterations_completed', 'sr_stop_reason', 'G1_Lx2_Ly2', 'G2xx_Lx2_Ly2']
    with report.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            result = {k: row.get(k, '') for k in fields}
            if row.get('run_dir'):
                run = Path(row['run_dir'])
                status = run / 'logs/status.txt'
                result['status'] = status.read_text().strip() if status.exists() else 'UNKNOWN'
                summary = run / 'results/summary.txt'
                if summary.exists():
                    values = {}
                    for line in summary.read_text().splitlines():
                        if ' = ' in line:
                            k, v = line.split(' = ', 1)
                            values[k] = v
                    for k in fields:
                        if k in values:
                            result[k] = values[k]
                    corr = run / 'results/correlations.csv'
                    if corr.exists():
                        for r in csv.DictReader(corr.open()):
                            if int(r['dx']) == int(values['Lx'])//2 and int(r['dy']) == int(values['Ly'])//2:
                                result['G1_Lx2_Ly2'] = r['G1']
                                result['G2xx_Lx2_Ly2'] = r['G2_xx']
            writer.writerow(result)
    print('Report: ' + str(report))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=Path('examples/delta_scan.inp'))
    parser.add_argument('--deltas', type=float, nargs='+', default=[1e-2, 1e-4, 1e-6, 1e-8, 1e-10])
    parser.add_argument('--samples', type=int)
    parser.add_argument('--sr-iterations', type=int)
    parser.add_argument('--submit', action='store_true', help='submit sequential Slurm jobs')
    parser.add_argument('--collect', type=Path, help='collect an existing manifest.csv, without submitting')
    args = parser.parse_args()
    if args.collect:
        collect(args.collect.resolve())
        return
    project = Path(__file__).resolve().parents[1]
    base = args.base.resolve()
    values = read_input(base)
    if args.samples is not None:
        values['samples'] = str(args.samples)
    if args.sr_iterations is not None:
        values['sr_iterations'] = str(args.sr_iterations)
    if len(set(args.deltas)) != len(args.deltas) or any(not math.isfinite(x) or x <= 0 for x in args.deltas):
        parser.error('deltas must be finite, positive and distinct')
    stamp = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y%m%d_%H%M%S_%f')
    folder = project / 'scan_inputs' / stamp
    folder.mkdir(parents=True)
    (folder/'base.original.inp').write_bytes(base.read_bytes())
    manifest = folder / 'manifest.csv'
    previous = ''
    with manifest.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['delta_requested', 'input', 'job_id', 'run_dir', 'status'])
        writer.writeheader()
        for delta in args.deltas:
            v = dict(values)
            v.update(delta=format(delta, '.17g'), optimize_delta='0', correlations='1')
            inp = folder / ('delta_{}.inp'.format(format(delta, '.12g')))
            inp.write_text('# Independent fixed-delta optimization; same initial chi1 and Delta.\n' +
                           ''.join('{} = {}\n'.format(k, value) for k, value in v.items()))
            row = dict(delta_requested=v['delta'], input=str(inp), job_id='', run_dir='', status='PREPARED')
            if args.submit:
                env = dict(os.environ)
                if previous:
                    env['VMC_DEPENDENCY'] = 'afterany:' + previous
                proc = subprocess.run(['bash', str(project/'run_optimize.sh'), str(inp)],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      universal_newlines=True, env=env)
                print(proc.stdout, end='', flush=True)
                if proc.returncode:
                    row['status'] = 'SUBMISSION_FAILED'
                    writer.writerow(row); out.flush()
                    raise RuntimeError(proc.stderr or proc.stdout)
                for line in proc.stdout.splitlines():
                    if line.startswith('Submitted job: '):
                        row['job_id'] = line.split(': ', 1)[1]
                        previous = row['job_id'].split(';')[0]
                    if line.startswith('Run directory: '):
                        row['run_dir'] = line.split(': ', 1)[1]
                row['status'] = 'SUBMITTED'
            writer.writerow(row);out.flush()
    print('Manifest: ' + str(manifest), flush=True)


if __name__ == '__main__':
    main()
