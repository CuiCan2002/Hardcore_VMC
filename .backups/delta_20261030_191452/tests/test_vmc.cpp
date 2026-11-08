#include "vmc.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace vmc;
void check(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
void near(Real a,Real b,Real tol,const char* text){check(std::abs(a-b)<=tol*(1+std::abs(b)),text);}
Real brute(const Kernel& k,const std::vector<int>& c) {
    if(c.empty())return 1;
    Real s=0;for(std::size_t b=1;b<c.size();++b){std::vector<int> rest;for(std::size_t j=1;j<c.size();++j)if(j!=b)rest.push_back(c[j]);s+=k.at(c[0],c[b])*brute(k,rest);}return s;
}
// Independently build H using diagonal bonds and their two spectator corners,
// instead of production code's center/delta enumeration.
std::vector<std::vector<Real>> dense_h(const Config& c,const Lattice& l,const std::vector<std::vector<int>>& basis) {
    std::map<std::vector<int>,int> lookup;for(int i=0;i<int(basis.size());++i)lookup[basis[i]]=i;
    std::vector<std::vector<Real>> H(basis.size(),std::vector<Real>(basis.size()));
    for(int col=0;col<int(basis.size());++col) {
        std::vector<int> n(l.n);for(int i:basis[col])n[i]=1;
        auto hop=[&](int from,int to,Real coeff){if(!n[from]||n[to])return;auto next=basis[col];*std::find(next.begin(),next.end(),from)=to;std::sort(next.begin(),next.end());H[lookup.at(next)][col]+=coeff;};
        for(int i=0;i<l.n;++i) {
            for(auto d:std::vector<std::array<int,2>>{{1,0},{0,1}}) {
                int j=l.shift(i,d[0],d[1]);H[col][col]+=c.VB*n[i]*n[j];hop(i,j,-c.JB);hop(j,i,-c.JB);
            }
            for(int dy:{-1,1}) {
                int j=l.shift(i,1,dy),s1=l.shift(i,1,0),s2=l.shift(i,0,dy);
                Real coeff=-c.K*(n[s1]+n[s2]);hop(i,j,coeff);hop(j,i,coeff);
            }
        }
    }return H;
}
int main()try {
    Config c;c.Lx=c.Ly=4;c.Nb=2;c.validate();Lattice l(c);Kernel k(c,l,c.alpha);
    check(l.bonds.size()==32&&l.corners.size()==64,"lattice multiplicity");
    for(int i=0;i<l.n;++i){near(k.at(i,i),0,1e-14L,"kernel diagonal");for(int j=0;j<l.n;++j){near(k.at(i,j),k.at(j,i),1e-14L,"kernel symmetry");near(k.at(i,j),k.at(l.shift(i,1,1),l.shift(j,1,1)),1e-13L,"translation");}}
    Config big;Lattice lb(big);Kernel kb(big,lb,big.alpha);
    std::vector<int> occupied{0,1,8,9,17,22,38,49};
    Hafnian h(kb,occupied);Real amp=h.amplitude();near(amp,brute(kb,occupied),1e-18L,"8-particle Hafnian pairing enumeration");
    for(int p=0;p<8;++p) {auto moved=occupied;moved[p]=63;Hafnian direct(kb,moved);near(h.moved(p,63)/amp,direct.amplitude()/amp,1e-11L,"minor-based ratio");}
    auto deriv=h.log_derivative();
    for(int p=0;p<3;++p) {
        Vec plus=big.alpha,minus=big.alpha;Real eps=1e-6L;plus[p]+=eps;minus[p]-=eps;
        Kernel kp(big,lb,plus),km(big,lb,minus);Hafnian hp(kp,occupied),hm(km,occupied);
        near(deriv[p],(hp.amplitude()-hm.amplitude())/(2*eps*amp),1e-6L,"analytic logarithmic derivative");
    }
    Kernel rank_one=kb;for(int i=0;i<lb.n;++i)for(int j=0;j<lb.n;++j)rank_one.f[i*lb.n+j]=(i==j?0:1);
    Hafnian hr(rank_one,occupied);near(hr.amplitude(),105,1e-14L,"rank-one (2M-1)!! identity");
    std::vector<std::vector<int>> basis;for(int i=0;i<l.n;++i)for(int j=i+1;j<l.n;++j)basis.push_back({i,j});
    auto H=dense_h(c,l,basis);std::vector<Real> psi;Real norm=0,numer=0;
    for(auto state:basis){Hafnian hs(k,state);psi.push_back(hs.amplitude());norm+=psi.back()*psi.back();}
    for(int i=0;i<int(basis.size());++i) {
        Real action=0;for(int j=0;j<int(basis.size());++j){near(H[i][j],H[j][i],1e-14L,"Hamiltonian Hermiticity");action+=H[i][j]*psi[j];}
        numer+=psi[i]*action;
        if(psi[i]!=0){Hafnian hs(k,basis[i]);near(local_energy(c,l,basis[i],hs).total(),action/psi[i],1e-11L,"local energy vs explicit Hamiltonian");}
    }
    Real exact=numer/norm;
    c.samples=40000;c.warmup_sweeps=1000;c.block_size=200;c.sweeps_between=2;
    auto s=sample(c,l,k,987654);
    std::cout<<"4x4 Nb=2 exact variational energy = "<<double(exact)<<"; MC = "<<double(s.energy)<<" +/- "<<double(s.error)<<'\n';
    check(std::abs(s.energy-exact)<6*s.error+0.005L,"MC energy vs exact variational sum");
    // Verify the full force against a finite difference of exact variational energy.
    Vec mean{},oe{};for(int i=0;i<int(basis.size());++i){Hafnian hs(k,basis[i]);if(psi[i]==0)continue;auto o=hs.log_derivative();Real w=psi[i]*psi[i]/norm,E=local_energy(c,l,basis[i],hs).total();for(int p=0;p<3;++p){mean[p]+=w*o[p];oe[p]+=w*o[p]*E;}}
    auto exact_energy=[&](const Vec& alpha){Kernel kt(c,l,alpha);std::vector<Real> v;Real z=0,en=0;for(auto state:basis){Hafnian hh(kt,state);v.push_back(hh.amplitude());z+=v.back()*v.back();}for(int i=0;i<int(v.size());++i)for(int j=0;j<int(v.size());++j)en+=v[i]*H[i][j]*v[j];return en/z;};
    for(int p=0;p<3;++p){Vec ap=c.alpha,am=c.alpha;ap[p]+=1e-6L;am[p]-=1e-6L;near(2*(oe[p]-mean[p]*exact),(exact_energy(ap)-exact_energy(am))/2e-6L,1e-6L,"energy gradient factor and sign");}
    Statistics fake;fake.metric={Vec{2,0.1L,0},Vec{0.1L,1,0.2L},Vec{0,0.2L,3}};fake.force={1,2,3};auto d=sr_direction(fake,0.05L,0.001L);
    for(int p=0;p<3;++p){Real sum=0.001L*d[p];for(int q=0;q<3;++q)sum+=fake.metric[p][q]*d[q];near(sum,-0.05L*fake.force[p],1e-12L,"SR linear solve");}
    bool rejected=false;try{Vec bad{0.05L,0.05L,-1};Kernel invalid(c,l,bad);}catch(const std::exception&){rejected=true;}check(rejected,"unstable reference rejected");
    rejected=false;try{Config odd=c;odd.Nb=3;odd.validate();}catch(const std::exception&){rejected=true;}check(rejected,"odd sector rejected");
    // Full matrices vs independent Fock-space operator action on 3x3, Nb=4.
    // This includes x/y cross channels and overlapping, equal, disjoint bonds.
    Config cc;cc.Lx=cc.Ly=3;cc.Nb=4;cc.alpha[2]=8;
    Lattice lc(cc);Kernel kc(cc,lc,cc.alpha);CorrelationAccumulator measured(lc);
    std::map<unsigned,Real> amplitudes;Real z=0;
    for(unsigned mask=0;mask<(1U<<lc.n);++mask) {
        std::vector<int> state;for(int i=0;i<lc.n;++i)if(mask&(1U<<i))state.push_back(i);
        if(state.size()!=4)continue;
        Real a=brute(kc,state);amplitudes[mask]=a;z+=a*a;
        if(a!=0){Hafnian hh(kc,state);measured.add(lc,kc,state,hh,a*a);}
    }
    auto expectation=[&](const std::vector<int>& remove,const std::vector<int>& create) {
        Real sum=0;
        for(auto entry:amplitudes) {
            unsigned next=entry.first;bool valid=true;
            for(int i:remove){if(!(next&(1U<<i))){valid=false;break;}next^=1U<<i;}
            if(!valid)continue;
            for(int i:create){if(next&(1U<<i)){valid=false;break;}next|=1U<<i;}
            if(valid)sum+=entry.second*amplitudes.at(next);
        }
        return sum/z;
    };
    Real trace=0;
    for(int i=0;i<lc.n;++i)for(int j=0;j<lc.n;++j) {
        near(measured.one[i*lc.n+j]/measured.weight,expectation({j},{i}),1e-12L,"full one-body matrix");
        if(i==j)trace+=measured.one[i*lc.n+j]/measured.weight;
    }
    near(trace,cc.Nb,1e-12L,"one-body trace Nb");
    for(std::size_t a=0;a<lc.bonds.size();++a)for(std::size_t b=0;b<lc.bonds.size();++b) {
        auto target=lc.bonds[a],source=lc.bonds[b];
        near(measured.pair[a*lc.bonds.size()+b]/measured.weight,
             expectation({source[1],source[0]},{target[1],target[0]}),1e-12L,"full NN pair matrix");
    }
    std::cout<<"Full one-body and NN pair matrices match exact Fock-space sums.\n";
    std::cout<<"All physics and numerical tests passed.\n";return 0;
}catch(const std::exception& e){std::cerr<<"TEST FAILED: "<<e.what()<<'\n';return 1;}
