#!/usr/bin/env python3
"""Integration tests inside a Slurm allocation; Python 3.6+, no NumPy needed."""
import csv, math, pathlib, subprocess, sys

project=pathlib.Path(sys.argv[1]).resolve()
root=pathlib.Path(sys.argv[2]).resolve();root.mkdir(exist_ok=False)
base=dict(Lx=4,Ly=4,Nb=4,JB=0.06,VB=2.5,K=1,chi1=0.05,Delta=0.2,delta=0.66,
          samples=4600,block_size=100,warmup_sweeps=100,sweeps_between=2,seed=12345,
          sr_iterations=0,sr_step=0.02,sr_shift=0.001,correlations=1,diagnostics=1)
def run(name,np=None,**changes):
    c=dict(base);c.update(changes);out=root/name;c['output']=str(out)
    inp=root/(name+'.inp');inp.write_text(''.join('{} = {}\n'.format(k,v) for k,v in c.items()))
    cmd=[str(project/'hardcore_vmc'),str(inp)] if np is None else ['mpirun','-np',str(np),str(project/'hardcore_vmc_mpi'),str(inp)]
    with (root/(name+'.log')).open('w') as f: subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=180)
    return out
def rows(path):
    with path.open() as f:return list(csv.DictReader(f))
def near(a,b,label,tol=2e-10):
    assert abs(a-b)<=tol*max(1,abs(a),abs(b)),(label,a,b)
def matrix(path):
    return [float(x.split()[-1]) for x in path.read_text().splitlines() if x and not x.startswith('#')]

a=run('serial',sr_iterations=2);b=run('mpi1',1,sr_iterations=2)
for name in ['samples.csv','diagnostics.csv','optimization.csv','one_body_correlation.dat','nn_pair_correlation.dat']:
    assert (a/name).read_bytes()==(b/name).read_bytes(),name
print('PASS: MPI rank=1 reproduces serial trajectories, SR, diagnostics and matrices',flush=True)

for ranks in [2,4]:
    out=run('mpi'+str(ranks),ranks)
    summaries=rows(out/'chain_summary.csv');raw=[];replays=[];counts=[]
    for row in summaries:
        rank=int(row['rank']);n=int(row['samples']);counts.append(n)
        chain=out/'chains'/('chain_%04d'%rank)
        replay=run('replay{}_{}'.format(ranks,rank),samples=n,seed=row['seed']);replays.append(replay)
        for name in ['samples.csv','diagnostics.csv']:
            assert (chain/name).read_bytes()==(replay/name).read_bytes(),name
        raw.extend(rows(chain/'samples.csv'))
    assert sum(counts)==base['samples']
    n=len(raw);e=[float(x['E_total']) for x in raw];em=sum(e)/n
    o=[[float(x['O_'+key]) for x in raw] for key in ['chi1','Delta','delta']]
    om=[sum(x)/n for x in o];pooled={k:float(v) for k,v in rows(out/'pooled_statistics.csv')[0].items()}
    near(pooled['E'],em,'pooled E');near(pooled['variance'],sum((v-em)**2 for v in e)/n,'variance')
    for p in range(3):
        near(pooled['O'+str(p)],om[p],'mean O')
        near(pooled['F'+str(p)],sum((x-om[p])*(y-em) for x,y in zip(o[p],e))/n,'force')
        for q in range(3):near(pooled['S{}{}'.format(p,q)],sum((x-om[p])*(y-om[q]) for x,y in zip(o[p],o[q]))/n,'metric')
    near(pooled['error'],math.sqrt(sum((int(r['samples'])/n*float(r['blocked_error']))**2 for r in summaries)),'error')
    for name in ['one_body_correlation.dat','nn_pair_correlation.dat']:
        mats=[matrix(r/name) for r in replays];m=matrix(out/name)
        for j,x in enumerate(m):near(x,sum(counts[r]*mats[r][j] for r in range(ranks))/n,name)
    print('PASS: {} chains, unequal counts {}, independent replay and pooled S/F/correlations'.format(ranks,counts),flush=True)

out=run('mpi_sr',4,sr_iterations=2)
for it,row in enumerate(rows(out/'optimization.csv')):
    fixed=run('sr_replay'+str(it),4,seed=base['seed']+it,chi1=row['chi1'],Delta=row['Delta'],delta=row['delta'])
    stats=rows(fixed/'pooled_statistics.csv')[0]
    near(float(row['E_total']),float(stats['E']),'SR pooled E',1e-8)
    for p,key in enumerate(['chi1','Delta','delta']):near(float(row['F_'+key]),float(stats['F'+str(p)]),'SR pooled force',1e-8)
print('PASS: synchronized SR iterations reproduce independent pooled evaluations',flush=True)

bad=root/'bad.inp';bad.write_text((root/'mpi4.inp').read_text().replace('samples = 4600','samples = 400'))
with (root/'expected_failure.log').open('w') as f:
    ret=subprocess.run(['mpirun','-np','4',str(project/'hardcore_vmc_mpi'),str(bad)],stdout=f,stderr=subprocess.STDOUT,timeout=30)
assert ret.returncode!=0
print('PASS: invalid per-chain sample budget fails collectively without hanging',flush=True)
print('ALL MPI INTEGRATION TESTS PASSED',flush=True)
