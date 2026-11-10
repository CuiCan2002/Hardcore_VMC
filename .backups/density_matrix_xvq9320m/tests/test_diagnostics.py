"""Run with Python + NumPy; checks analysis against a known AR(1) process."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import numpy as np

script=Path(__file__).resolve().parents[1]/'tools/analyze_diagnostics.py'
spec=importlib.util.spec_from_file_location('diagnostics',script)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
rng=np.random.RandomState(31415)
x=np.zeros(120000);noise=rng.normal(size=len(x))
for i in range(1,len(x)):x[i]=0.8*x[i-1]+noise[i]
x=x[20000:]
tau,window,_=module.acf_window(x,5)
assert abs(tau-9)<1.8,(tau,window)  # (1+phi)/(1-phi)
with tempfile.TemporaryDirectory() as folder:
    p=Path(folder);trace=p/'diagnostics.csv'
    with trace.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['phase','sample','sweep','X','constant'])
        writer.writerow(['warmup',0,0,1e9,99])
        for i,v in enumerate(x,1):writer.writerow(['measurement',i,500+4*i,v,4])
    module.analyze(trace,p/'report',discard=100)
    rows={r['observable']:r for r in csv.DictReader((p/'report/summary.csv').open())}
    assert abs(float(rows['X']['mean'])-np.mean(x[100:]))<1e-12
    assert rows['constant']['tau']=='' and 'constant:' in rows['constant']['status']
    meta=json.loads((p/'report/metadata.json').read_text())
    assert meta['interval_sweeps']==4 and meta['warmup_rows_excluded']==1 and meta['N']==99900
print('PASS: AR(1) autocorrelation, warmup exclusion, extra discard, time units and constant-series handling.')
