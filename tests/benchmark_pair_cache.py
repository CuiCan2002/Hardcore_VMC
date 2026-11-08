#!/usr/bin/env python3
"""Same-seed end-to-end benchmark; run inside a 32-task Slurm allocation."""
import csv, json, pathlib, subprocess, sys, time

project=pathlib.Path(sys.argv[1]); root=pathlib.Path(sys.argv[2]); old=pathlib.Path(sys.argv[3])
root.mkdir(parents=True,exist_ok=False)
base=dict(Lx=12,Ly=12,Nb=18,JB=0,VB=2.5,K=1,
          chi1='0.00383304127420137124461',Delta='0.36754805197563652518',delta='0.0781742545329976973556',
          warmup_sweeps=100,sweeps_between=2,block_size=1,seed=983001,sr_iterations=0,
          Ifjastrow=0,diagnostics=1,diagnostic_stride=1)
times={}
def run(name,binary,**changes):
    c=dict(base);c.update(changes);c['output']=str(root/name)
    inp=root/(name+'.inp');inp.write_text(''.join('{} = {}\n'.format(k,v) for k,v in c.items()))
    start=time.monotonic()
    with (root/(name+'.log')).open('w') as f:
        subprocess.run(['mpirun','-np','32',str(binary),str(inp)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=600)
    times[name]=time.monotonic()-start
    print(name,times[name],flush=True)
    return root/name

a=run('old_full',old,samples=320,correlations=1)
b=run('new_full',project/'hardcore_vmc_mpi',samples=320,correlations=1,pair_measure_stride=1)
c=run('new_stride10',project/'hardcore_vmc_mpi',samples=320,correlations=1,pair_measure_stride=10)
for rank in range(32):
    rel=pathlib.Path('chains')/('chain_%04d'%rank)/'samples.csv'
    assert (a/rel).read_bytes()==(b/rel).read_bytes()==(c/rel).read_bytes(),str(rel)
print('PASS: old/new/stride10 have identical 32-chain energy and derivative trajectories',flush=True)
for name in ['one_body_correlation.dat','nn_pair_correlation.dat']:
    def values(p):return [float(x.split()[-1]) for x in p.read_text().splitlines() if x and not x.startswith('#')]
    x,y=values(a/name),values(b/name)
    err=max(abs(u-v) for u,v in zip(x,y))
    assert len(x)==len(y) and err<1e-10,(name,err)
    print('PASS:',name,'max abs error',err,flush=True)
assert (b/'one_body_correlation.dat').read_bytes()==(c/'one_body_correlation.dat').read_bytes()
for name,binary in [('old_sampling',old),('new_sampling',project/'hardcore_vmc_mpi')]:
    run(name,binary,samples=16000,correlations=0,diagnostics=0)
for rank in range(32):
    rel=pathlib.Path('chains')/('chain_%04d'%rank)/'samples.csv'
    assert (root/'old_sampling'/rel).read_bytes()==(root/'new_sampling'/rel).read_bytes()
print('PASS: persistent cache preserves 32-chain sampling trajectories',flush=True)
(root/'timings.json').write_text(json.dumps(times,indent=2))
print('ALL PERFORMANCE REGRESSION TESTS PASSED',flush=True)
