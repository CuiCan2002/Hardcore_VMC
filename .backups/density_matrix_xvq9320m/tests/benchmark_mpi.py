#!/usr/bin/env python3
"""Small fixed-total-work benchmark; run in an allocated 32-core Slurm job."""
import csv,pathlib,subprocess,sys,time
p=pathlib.Path(sys.argv[1]);root=pathlib.Path(sys.argv[2]);root.mkdir(exist_ok=False)
rows=[]
for ranks in [1,4,16,32]:
    inp=root/('p%d.inp'%ranks)
    inp.write_text('''Lx = 8
Ly = 8
Nb = 8
JB = 0
VB = 2.5
K = 1
chi1 = 0.0053733039987898254
Delta = 0.41591949129532056
delta = 0.12297667789763767
samples = 64000
block_size = 100
warmup_sweeps = 500
sweeps_between = 2
sr_iterations = 2
sr_step = 0.02
sr_shift = 0.001
diagnostics = 1
correlations = 1
output = '''+str(root/('p%d'%ranks))+'\n')
    start=time.monotonic()
    with (root/('p%d.log'%ranks)).open('w') as f:
        subprocess.run(['mpirun','-np',str(ranks),str(p/'hardcore_vmc_mpi'),str(inp)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=600)
    seconds=time.monotonic()-start
    rows.append(dict(ranks=ranks,seconds=seconds,speedup=rows[0]['seconds']/seconds if rows else 1))
    print(rows[-1],flush=True)
with (root/'timings.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['ranks','seconds','speedup']);w.writeheader();w.writerows(rows)
