#!/usr/bin/env python3
"""Analyze a final fixed-parameter diagnostics.csv. Warmup is always excluded.

NumPy enables FFT autocorrelation. Without NumPy, blocking and time-segment
reports still work. Neither method certifies equilibration of a trapped chain.
tau convention = 1 + 2 sum rho; ESS = N/tau.
"""
import argparse
import array
import csv
import json
import math
from pathlib import Path
try:
    import numpy as np
except ImportError:
    np = None


def acf_window(x, c):
    n=len(x);centered=x-x.mean()
    spectrum=np.fft.rfft(centered,n=1<<((2*n-1).bit_length()))
    cov=np.fft.irfft(spectrum*spectrum.conj())[:n]/np.arange(n,0,-1)
    rho=cov/cov[0];taus=2*np.cumsum(rho)-1
    candidates=np.flatnonzero((np.arange(n)>=c*taus)&(taus>0)&(np.arange(n)>=5)&(np.arange(n)<=n//20))
    if len(candidates)==0:
        return None,None,rho
    window=int(candidates[0]);return float(taus[window]),window,rho


def analyze(path,out,discard=0):
    columns={};sweeps=array.array('q');warmup=0;seen=0
    with path.open() as stream:
        reader=csv.DictReader(stream)
        names=[x for x in reader.fieldnames if x not in ['phase','sample','sweep']]
        columns={k:array.array('d') for k in names}
        for row in reader:
            if row['phase']=='warmup':warmup+=1;continue
            if row['phase']!='measurement':raise ValueError('Unknown phase')
            seen+=1
            if seen<=discard:continue
            sweeps.append(int(row['sweep']))
            for key in names:
                v=float(row[key])
                if not math.isfinite(v):raise ValueError('Nonfinite '+key)
                columns[key].append(v)
    n=len(sweeps)
    if n<100:raise ValueError('Need at least 100 retained measurement records')
    interval=sweeps[1]-sweeps[0]
    if interval<=0 or any(sweeps[i]-sweeps[i-1]!=interval for i in range(1,n)):
        raise ValueError('Measurement sweep spacing must be uniform')
    out.mkdir(parents=True,exist_ok=False)
    summaries=[];blocks=[];segments=[]
    for name,values in columns.items():
        mean=sum(values)/n;var=sum((v-mean)**2 for v in values)/n
        row=dict(observable=name,N=n,mean=mean,tau=None,tau_sweeps=None,ESS=None,
                 mcse=None,window=None,tau_c10=None,status='')
        for b in [10,50,100,200,500,1000,5000,10000,50000]:
            count=n//b
            if count<20:continue
            means=[sum(values[i*b:(i+1)*b])/b for i in range(count)]
            m=sum(means)/count;variance=sum((x-m)**2 for x in means)/(count-1)
            blocks.append(dict(observable=name,block_records=b,block_sweeps=b*interval,
                               blocks=count,used_records=count*b,mean=m,mcse=math.sqrt(variance/count),
                               tau_from_blocking=b*variance/var if var>0 else None))
        for i in range(10):
            a=i*n//10;b=(i+1)*n//10
            segments.append(dict(observable=name,segment=i+1,start_sweep=sweeps[a],
                                 end_sweep=sweeps[b-1],N=b-a,mean=sum(values[a:b])/(b-a)))
        if var==0:
            row['status']='constant: conserved OR unvisited sectors; no tau estimate'
        elif np is None:
            row['status']='NumPy unavailable: inspect blocking.csv; FFT not computed'
        else:
            x=np.asarray(values);tau,window,rho=acf_window(x,5)
            tau10,_,_=acf_window(x,10)
            row['tau_c10']=tau10
            if tau is None:
                row['status']='unresolved autocorrelation window; increase chain length'
            else:
                row.update(tau=tau,tau_sweeps=tau*interval,ESS=n/tau,
                           mcse=math.sqrt(var*tau/n),window=window,status='estimated; not proof of equilibration')
                if n<50*tau or tau10 is None or abs(tau10-tau)>0.2*tau:
                    row['status']='window/length sensitive; inspect long tails and independent chains'
            upto=min(n,2001)
            np.savetxt(str(out/(name+'_acf.dat')),np.column_stack((np.arange(upto),np.arange(upto)*interval,rho[:upto])),
                       header='lag_records lag_sweeps normalized_autocorrelation')
        summaries.append(row)
    for filename,rows in [('summary.csv',summaries),('blocking.csv',blocks),('segments.csv',segments)]:
        if rows:
            with (out/filename).open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (out/'metadata.json').write_text(json.dumps(dict(input=str(path),warmup_rows_excluded=warmup,
        measurement_rows_discarded=discard,N=n,interval_sweeps=interval,numpy_available=np is not None,
        tau_convention='1 + 2 sum rho; ESS=N/tau; positive self-consistent window W >= c*tau',
        warning='Single-chain diagnostic. Different observables can have different slow modes. Compare independent chains, burn-in lengths and block-size plateaus. No automatic equilibrium verdict.'),indent=2))
    print('Reports: '+str(out))
    if np is None:print('NumPy not installed: blocking/time-segment reports only. Install NumPy or run on a NumPy-equipped machine for FFT ACF/ESS.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trace',type=Path)
    parser.add_argument('--discard',type=int,default=0,help='additional measurement records to discard after warmup')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    if args.discard<0:parser.error('discard must be nonnegative')
    analyze(args.trace.resolve(),(args.output_dir or args.trace.parent/'diagnostics_analysis').resolve(),args.discard)
