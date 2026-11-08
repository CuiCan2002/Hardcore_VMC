
import argparse
import bisect
import itertools
import math
import random
from pathlib import Path

import numpy as np

PARAMS = {
    74747: (.003736187223468507,3.346737716136566,9.45182214597509,
            .529151807277713,.3134820432651763,.1925557513348328),
    74740: (.002831103997648657,2.387473057346181,8.902565115769077,
            .5228164248999664,.2967780116443726,.1829371569008239),
}

def model(L, n, params):
    chi, gap, delta, v1, v2, v3 = params
    x,y=np.indices((L,L))
    cx,cy=np.cos(2*np.pi*x/L),np.cos(2*np.pi*y/L)
    gamma=cx+cy
    B=-3*gap*gamma
    bound=3*chi*gamma+8*n/(L*L)*cx*cy+abs(B)
    d=bound.max()-bound+delta
    f=np.fft.ifft2(-B/(abs(B)+d+np.sqrt(d*(d+2*abs(B))))).real
    f=(f+f[(-x)%L,(-y)%L])/2
    assert f.ravel()[1:].min()>0, 'Pair sampler requires positive f'
    f[0,0]=0
    i=np.arange(L*L)
    dx=(i[:,None]%L-i[None,:]%L)%L
    dy=(i[:,None]//L-i[None,:]//L)%L
    F=f[dx,dy]
    dx,dy=np.minimum(dx,L-dx),np.minimum(dy,L-dy)
    V=v1*(dx+dy==1)+v2*((dx==1)&(dy==1))+v3*((dx+dy==2)&((dx==0)|(dy==0)))
    bonds=[(i,j) for i in range(L*L)
           for j in (i//L*L+(i+1)%L,(i+L)%(L*L))]
    return F.tolist(),V.tolist(),bonds

def pair_bin(L,a,b):
    i,j=a//2,b//2
    return a%2,b%2,(j//L-i//L)%L,(j%L-i%L)%L

class Chain:
    def __init__(self,L,n,data,bias,seed):
        self.L,self.n,self.F,self.V,self.bonds=L,n,*data
        self.bias=bias
        self.rng=random.Random(seed)
        a=self.rng.randrange(len(self.bonds))
        self.ends=[a,a]
        self.bg=set(self.rng.sample(list(set(range(L*L))-set(self.bonds[a])),n-2))
        sites=list(sorted(self.bg))+list(self.bonds[a])
        self.sites=[sites[:],sites[:]]
        self.match=[[i^1 for i in range(n)] for _ in range(2)]
        self.tables=[None,None]
        self.attempts=[0]*3
        self.changes=[0]*3
        self.caps=0

    def log_ratio(self,r,changes):
        sites,m=self.sites[r],self.match[r]
        value=0.
        for i,a in changes.items():
            old=sites[i]
            j=m[i]
            if j not in changes or i<j:
                value+=math.log(self.F[a][changes.get(j,sites[j])]/self.F[old][sites[j]])
            for j,b in enumerate(sites):
                if j not in changes or i<j:
                    value+=self.V[old][b]-self.V[a][changes.get(j,b)]
        return value

    def accept(self,value):
        return value>=0 or self.rng.random()<math.exp(value)

    def apply(self,r,changes):
        for i,a in changes.items():
            self.sites[r][i]=a
        self.tables[r]=None

    def move_end(self):
        self.attempts[0]+=1
        r=self.rng.randrange(2)
        old,other=self.ends[r],self.ends[1-r]
        # Mixture of a global proposal and an explicit proposal to close.
        # Its reverse probability is evaluated with the other endpoint fixed.
        close=.2
        count=len(self.bonds)
        new=other if self.rng.random()<close else self.rng.randrange(count)
        if new==old:return
        a,b=self.bonds[new]
        if a in self.bg or b in self.bg:return
        changes={self.n-2:a,self.n-1:b}
        q_old=(1-close)/count+close*(old==other)
        q_new=(1-close)/count+close*(new==other)
        log_ratio=self.log_ratio(r,changes)+math.log(q_old/q_new)
        log_ratio+=math.log(self.bias)*((new==other)-(old==other))
        if self.accept(log_ratio):
            self.ends[r]=new
            self.apply(r,changes)
            self.changes[0]+=1

    def move_background(self):
        self.attempts[1]+=1
        p=self.rng.randrange(self.n-2)
        a=self.rng.randrange(self.L*self.L)
        if a in self.bg or any(a in self.bonds[b] for b in self.ends):return
        changes={p:a}
        if self.accept(sum(self.log_ratio(r,changes) for r in range(2))):
            old=self.sites[0][p]
            self.bg.remove(old)
            self.bg.add(a)
            for r in range(2):self.apply(r,changes)
            self.changes[1]+=1

    def worm(self):
        self.attempts[2]+=1
        r=self.rng.randrange(2)
        sites=self.sites[r]
        if self.tables[r] is None:
            rows=[list(itertools.accumulate(self.F[a][b] for b in sites)) for a in sites]
            self.tables[r]=(rows,max(row[-1] for row in rows))
        rows,normalizer=self.tables[r]
        m=self.match[r][:]
        head=self.rng.randrange(self.n)
        tail=m[head]
        for _ in range(10000):
            target=bisect.bisect_right(rows[tail],self.rng.random()*normalizer)
            if target==self.n or target==tail:continue
            if target==head:
                m[head]=tail
                m[tail]=head
                self.changes[2]+=m!=self.match[r]
                self.match[r]=m
                return
            displaced=m[target]
            m[tail]=target
            m[target]=tail
            tail=displaced
        self.caps+=1

    def step(self):
        move=self.rng.randrange(6)
        if move<2:self.move_end()
        elif move<4:self.move_background()
        else:self.worm()

def exact(L,n,data):
    F,V,bonds=data
    def haf(c):
        if not c:return 1.
        return sum(F[c[0]][p]*haf(c[1:i]+c[i+1:]) for i,p in enumerate(c) if i)
    norm=0.
    psi={}
    for C in itertools.combinations(range(L*L),n):
        J=sum(V[i][j] for i,j in itertools.combinations(C,2))
        psi[C]=math.exp(-J)*haf(C)
        norm+=psi[C]**2
    h=np.zeros((2,2,L,L))
    trace=0.
    for R in itertools.combinations(range(L*L),n-2):
        candidates=[(a,psi[tuple(sorted(R+bond))]) for a,bond in enumerate(bonds)
                    if not set(R).intersection(bond)]
        for a,pa in candidates:
            trace+=pa**2
            for b,pb in candidates:h[pair_bin(L,a,b)]+=pa*pb
    return trace/norm, h/h.sum()

def simulate(L,n,params,bias,steps,burn,seed,chains):
    data=model(L,n,params)
    hist=np.zeros((chains,20,2,2,L,L),dtype=np.int64)
    assert steps%20==0
    for k in range(chains):
        c=Chain(L,n,data,bias,seed+k)
        for _ in range(burn):c.step()
        c.attempts=[0]*3
        c.changes=[0]*3
        c.caps=0
        for i in range(steps):
            c.step()
            a,b=c.ends
            hist[(k,i//(steps//20))+pair_bin(L,a,b)]+=1
        diag=sum(hist[k,:,o,o,0,0].sum() for o in range(2))
        print('chain',k,'diag',diag,'accept/change endpoint/background/worm',
              [v/a for v,a in zip(c.changes,c.attempts)],'capped worms',c.caps,flush=True)
    return data,hist

def report(L,hist,bias,trace):
    blocks=hist.reshape(-1,2,2,L,L).astype(float)
    for o in range(2):blocks[:,o,o,0,0]/=bias
    def fraction(h):
        diag=h[0,0,0,0]+h[1,1,0,0]
        return h.sum()/(2*L*L*diag)
    total=blocks.sum(0)
    value=fraction(total)
    print('pooled Lambda_s/trace',value,flush=True)
    for size in (1,2,5,10,20):
        grouped=blocks.reshape(-1,size,2,2,L,L).sum(1)
        jk=np.array([fraction(total-b) for b in grouped])
        se=math.sqrt((len(jk)-1)*np.mean((jk-jk.mean())**2))
        print('jackknife blocks',len(jk),'SE',se,flush=True)
    for i,h in enumerate(blocks.reshape(len(hist),20,2,2,L,L).sum(1)):
        print('chain fraction',i,fraction(h))
    print('unweighted diagonal probability',1/(2*L*L*value))
    if trace is not None:
        print('Lambda_s',trace*value)
    return total/total.sum()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--job',type=int,choices=PARAMS,default=74747)
    p.add_argument('--size',type=int,default=12)
    p.add_argument('--particles',type=int,default=18)
    p.add_argument('--steps',type=int,default=400000,help='Recorded proposals per chain')
    p.add_argument('--burn',type=int,default=100000)
    p.add_argument('--chains',type=int,default=4)
    p.add_argument('--bias',type=float,default=20,help='Extra weight for a=b')
    p.add_argument('--seed',type=int,default=72700)
    p.add_argument('--trace',type=float,help='Independent measurement of Tr G')
    p.add_argument('--exact',action='store_true',help='Enumerate small systems only')
    p.add_argument('--output',type=Path,help='Optional new NPZ file with raw block histograms')
    a=p.parse_args()
    if not (a.size>=3 and 4<=a.particles<a.size**2 and a.particles%2==0
            and a.bias>0 and a.chains>=2 and a.steps>=20 and a.steps%20==0 and a.burn>=0):
        p.error('invalid size, particles, chain count, bias or steps')
    if a.exact and (a.size>4 or a.particles>4):p.error('exact mode is restricted to L<=4, Nb=4')
    data,hist=simulate(a.size,a.particles,PARAMS[a.job],a.bias,a.steps,a.burn,a.seed,a.chains)
    prob=report(a.size,hist,a.bias,a.trace)
    if a.exact:
        trace,reference=exact(a.size,a.particles,data)
        diag=reference[0,0,0,0]+reference[1,1,0,0]
        print('exact trace',trace,'exact Lambda_s/trace',1/(2*a.size**2*diag),
              'histogram total variation',abs(prob-reference).sum()/2)
    if a.output:
        with a.output.open('xb') as f:
            np.savez_compressed(f,hist=hist,bias=a.bias,prob=prob,job=a.job,particles=a.particles)

if __name__=='__main__':
    main()
