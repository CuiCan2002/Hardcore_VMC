#include "vmc.hpp"
#include <algorithm>
#include <cmath>
#if __has_include(<filesystem>)
#include <filesystem>
namespace vmc_fs = std::filesystem;
#else
// GCC 7: C++17 filesystem is provided by the experimental library.
#include <experimental/filesystem>
namespace vmc_fs = std::experimental::filesystem;
#endif
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#ifdef VMC_MPI
#include <mpi.h>
#endif

namespace vmc {
namespace {
int rank_id() {
#ifdef VMC_MPI
    int r=0;MPI_Comm_rank(MPI_COMM_WORLD,&r);return r;
#else
    return 0;
#endif
}
int world_size() {
#ifdef VMC_MPI
    int n=1;MPI_Comm_size(MPI_COMM_WORLD,&n);return n;
#else
    return 1;
#endif
}
void synchronize() {
#ifdef VMC_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
#ifdef VMC_MPI
std::uint64_t chain_seed(std::uint64_t seed,int rank) {
    if(rank==0)return seed; // preserve the original serial trajectory
    auto z=seed+0x9e3779b97f4a7c15ULL*std::uint64_t(rank);
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
    z=(z^(z>>27))*0x94d049bb133111ebULL;
    return z^(z>>31);
}
#endif
void require(bool ok,const std::string& message) {if(!ok) throw std::runtime_error(message);}
template<class T> void parse(const std::string& s,T& target) {
    std::istringstream in(s); T value{}; std::string extra;
    require(bool(in>>value)&&!(in>>extra),"Invalid value: "+s); target=value;
}
std::ofstream output(const std::string& path) {
    std::ofstream f(path);require(bool(f),"Cannot write "+path); f<<std::setprecision(17);return f;
}
bool usable(Real x) {return x!=0 && std::isfinite(x);}
int first_bit(std::uint32_t mask) {int i=0;while(!(mask&1U)){mask>>=1;++i;}return i;}
std::vector<int> occupancy(const Lattice& l,const std::vector<int>& c) {
    std::vector<int> n(l.n);for(int i:c)n[i]=1;return n;
}
}
Config Config::read(const std::string& path) {
    std::ifstream f(path);require(bool(f),"Cannot read "+path);
    Config c;std::string line;std::set<std::string> seen;int lineno=0;
    while(std::getline(f,line)) {
        ++lineno;line=line.substr(0,line.find('#'));
        std::replace(line.begin(),line.end(),'=',' ');
        std::istringstream in(line);std::string key,value,extra;
        if(!(in>>key))continue;
        require(bool(in>>value)&&!(in>>extra),"Expected key = value at line "+std::to_string(lineno));
        require(seen.insert(key).second,"Duplicate key: "+key);
        if(key=="Lx")parse(value,c.Lx);
        else if(key=="Ly")parse(value,c.Ly);
        else if(key=="Nb")parse(value,c.Nb);
        else if(key=="JB")parse(value,c.JB);
        else if(key=="VB")parse(value,c.VB);
        else if(key=="K")parse(value,c.K);
        else if(key=="chi1")parse(value,c.alpha[0]);
        else if(key=="Delta")parse(value,c.alpha[1]);
        else if(key=="delta")parse(value,c.alpha[2]);
        else if(key=="lambda_eff")throw std::runtime_error("lambda_eff is derived now; input delta = lambda_eff - lambda_c instead");
        else if(key=="stability_margin")throw std::runtime_error("stability_margin was replaced by strictly positive delta");
        else if(key=="local_probability")parse(value,c.local_probability);
        else if(key=="warmup_sweeps")parse(value,c.warmup_sweeps);
        else if(key=="Ifjastrow") {int b;parse(value,b);require(b==0||b==1,"Ifjastrow must be 0 or 1");c.Ifjastrow=b;}
        else if(key=="v1")parse(value,c.alpha[3]);
        else if(key=="v2")parse(value,c.alpha[4]);
        else if(key=="v3")parse(value,c.alpha[5]);
        else if(key=="samples")parse(value,c.samples);
        else if(key=="samples_sr")parse(value,c.samples_sr);
        else if(key=="sweeps_between")parse(value,c.sweeps_between);
        else if(key=="block_size")parse(value,c.block_size);
        else if(key=="sr_iterations")parse(value,c.sr_iterations);
        else if(key=="sr_step")parse(value,c.sr_step);
        else if(key=="sr_shift")parse(value,c.sr_shift);
        else if(key=="max_log_delta_step")parse(value,c.max_log_delta_step);
        else if(key=="optimize_delta") {int b;parse(value,b);require(b==0||b==1,"optimize_delta must be 0 or 1");c.optimize_delta=b!=0;}
        else if(key=="seed") {require(!value.empty()&&value[0]!='-',"seed must be nonnegative");parse(value,c.seed);}
        else if(key=="correlations") {int b;parse(value,b);require(b==0||b==1,"correlations must be 0 or 1");c.correlations=b!=0;}
        else if(key=="diagnostics") {int b;parse(value,b);require(b==0||b==1,"diagnostics must be 0 or 1");c.diagnostics=b!=0;}
        else if(key=="diagnostic_stride")parse(value,c.diagnostic_stride);
        else if(key=="pair_measure_stride")parse(value,c.pair_measure_stride);
        else if(key=="output")c.output=value;
        else throw std::runtime_error("Unknown input key: "+key);
    }
    c.validate();return c;
}
void Config::validate() const {
    require(Lx>=3&&Ly>=3&&Lx<=128&&Ly<=128&&Lx*Ly<=4096,"Require Lx,Ly >= 3 and Ns <= 4096 (nondegenerate PBC)");
    require(Nb>=2&&Nb<=20&&Nb%2==0&&Nb<Lx*Ly,"Require even 2 <= Nb <= 20 and Nb < Ns (exact Hafnian limit)");
    for(Real x:alpha)require(std::isfinite(x),"Nonfinite variational parameter");
    for(Real x:{JB,VB,K,alpha[0],alpha[1],alpha[2],local_probability,sr_step,sr_shift,max_log_delta_step})
        require(std::isfinite(x),"All real parameters must be finite");
    require(alpha[1]!=0&&(VB-4*K)!=0,"Bk identically zero: pair-projected state vanishes");
    require(alpha[2]>0,"delta must be strictly positive");
    require(max_log_delta_step>0,"max_log_delta_step must be positive");
    require(local_probability>=0&&local_probability<1,"Require 0 <= local_probability < 1 to retain global moves");
    require(warmup_sweeps>=0&&samples>=2&&sweeps_between>=1,"Invalid sampling controls");
    require(block_size>=1&&samples%block_size==0&&samples/block_size>=2,"samples must contain at least two complete blocks");
    require(samples_sr>=0,"samples_sr must be nonnegative (0 inherits samples)");
    if(sr_iterations>0)require(sr_samples()>=2&&sr_samples()%block_size==0&&sr_samples()/block_size>=2,"samples_sr must contain at least two complete blocks");
    require(sr_iterations>=0&&sr_step>0&&sr_shift>0,"Invalid SR controls");
    require(diagnostic_stride>=1,"diagnostic_stride must be positive");
    require(pair_measure_stride>=1,"pair_measure_stride must be positive");
    require(!output.empty(),"Empty output path");
}
Lattice::Lattice(const Config& c):lx(c.Lx),ly(c.Ly),n(lx*ly) {
    for(int i=0;i<n;++i) {
        bonds.push_back({i,shift(i,1,0)});bonds.push_back({i,shift(i,0,1)});
        for(int dx:{-1,1})for(int dy:{-1,1})corners.push_back({i,shift(i,dx,0),shift(i,0,dy)});
    }
}
int Lattice::site(int x,int y)const{return ((x%lx+lx)%lx)+lx*((y%ly+ly)%ly);}
int Lattice::shift(int i,int dx,int dy)const{return site(i%lx+dx,i/lx+dy);}
Kernel::Kernel(const Config& c,const Lattice& l,const Vec& a):n(l.n),f(n*n),
    min_margin(std::numeric_limits<Real>::infinity()),min_gap(min_margin) {
    jastrow=c.Ifjastrow;v=a;for(auto x:a)require(std::isfinite(x),"Nonfinite variational parameter");
    {shell.resize(n*n,-1);for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){int dx=std::abs(i%l.lx-j%l.lx),dy=std::abs(i/l.lx-j/l.lx);dx=std::min(dx,l.lx-dx);dy=std::min(dy,l.ly-dy);int t=-1;if(dx+dy==1)t=0;else if(dx==1&&dy==1)t=1;else if((dx==2&&dy==0)||(dx==0&&dy==2))t=2;shell[i*n+j]=shell[j*n+i]=t;}}
    for(auto& d:df)d.resize(n*n);
    const Real pi=std::acos(Real(-1)),g=c.VB-4*c.K,nu=Real(c.Nb)/n;
    std::vector<Real> fr(n);std::array<std::vector<Real>,3> dr;
    for(auto& d:dr)d.resize(n);
    require(g*a[1]!=0,"Zero pair kernel");
    require(std::isfinite(a[0])&&std::isfinite(a[1])&&std::isfinite(a[2])&&a[2]>0,"Require finite chi1, Delta and delta > 0");
    // Use exactly the same discrete grid for the boundary and the Fourier kernel.
    std::vector<Real> gammas(n),bs(n),boundaries(n);
    lambda_c=-std::numeric_limits<Real>::infinity();
    for(int ky=0;ky<l.ly;++ky)for(int kx=0;kx<l.lx;++kx) {
        int q=kx+l.lx*ky;Real x=2*pi*kx/l.lx,y=2*pi*ky/l.ly;
        Real gamma=std::cos(x)+std::cos(y),eta=std::cos(x)*std::cos(y);
        gammas[q]=gamma;bs[q]=2*g*a[1]*gamma;
        boundaries[q]=2*(c.JB-g*a[0])*gamma+8*c.K*nu*eta+std::abs(bs[q]);
        require(std::isfinite(boundaries[q]),"Nonfinite stability boundary");
        lambda_c=std::max(lambda_c,boundaries[q]);
    }
    lambda_eff=lambda_c+a[2];
    require(std::isfinite(lambda_eff),"Nonfinite derived lambda_eff");
    // Average active gradients at a numerical tie: a subgradient of max().
    // Symmetry-related maximizers have identical gradients. At a genuine branch
    // crossing the function is not differentiable; this convention is explicit.
    Real tie=16*std::numeric_limits<Real>::epsilon()*std::max(Real(1),std::abs(lambda_c));
    for(int q=0;q<n;++q)if(lambda_c-boundaries[q]<=tie) {
        ++boundary_modes;
        boundary_derivative[0]+=-2*g*gammas[q];
        boundary_derivative[1]+=(bs[q]>0?1:(bs[q]<0?-1:0))*2*g*gammas[q];
    }
    for(int p=0;p<2;++p)boundary_derivative[p]/=boundary_modes;
    for(int ky=0;ky<l.ly;++ky)for(int kx=0;kx<l.lx;++kx) {
        int q=kx+l.lx*ky;Real x=2*pi*kx/l.lx,y=2*pi*ky/l.ly;
        Real gamma=gammas[q],B=bs[q];
        // Do not subtract nearly equal A and |B|: retain small positive delta.
        Real margin=(lambda_c-boundaries[q])+a[2];
        Real A=std::abs(B)+margin,E=std::sqrt(margin)*std::sqrt(margin+2*std::abs(B));
        require(std::isfinite(E)&&E>0,"Gaussian gap not representable at this delta");
        Real den=A+E,fk=-B/den;
        min_margin=std::min(min_margin,margin);min_gap=std::min(min_gap,E);
        Vec da{boundary_derivative[0]+2*g*gamma,boundary_derivative[1],1};
        Vec db{0,2*g*gamma,0},dk{};
        for(int p=0;p<3;++p) {
            // Stable version of (A*dA-B*dB)/E near the soft mode.
            Real de=(margin*da[p]+std::abs(B)*(da[p]-(B>0?1:(B<0?-1:0))*db[p]))/E;
            dk[p]=-db[p]/den+(B/den)*(da[p]+de)/den;
            require(std::isfinite(dk[p]),"Nonfinite pairing derivative");
        }
        for(int r=0;r<n;++r) {
            Real phase=std::cos(x*(r%l.lx)+y*(r/l.lx))/n;
            fr[r]+=phase*fk;for(int p=0;p<3;++p)dr[p][r]+=phase*dk[p];
        }
    }
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j) {
        int r=l.site(i%l.lx-j%l.lx,i/l.lx-j/l.lx);
        f[i*n+j]=f[j*n+i]=fr[r];
        for(int p=0;p<3;++p)df[p][i*n+j]=df[p][j*n+i]=dr[p][r];
    } // diagonal and its derivatives remain zero after hard-core projection
}
std::array<Real,3> Kernel::counts(const std::vector<int>& c) const {
    std::array<Real,3> b{};
    for(std::size_t i=0;i<c.size();++i)for(std::size_t j=i+1;j<c.size();++j){int t=shell[c[i]*n+c[j]];if(t>=0)++b[t];}
    return b;
}
Real Kernel::change(const std::vector<int>& a,const std::vector<int>& b) const {
    if(!jastrow)return 0;
    auto x=counts(a),y=counts(b);Real d=0;
    for(int i=0;i<3;++i)d+=v[i+3]*(y[i]-x[i]);
    return d;
}
Hafnian::Hafnian(const Kernel& k,const std::vector<int>& c):kernel(k),occupied(c),
    memo(std::size_t(1)<<c.size(),std::numeric_limits<Real>::quiet_NaN()){memo[0]=1;}
Real Hafnian::solve(std::uint32_t mask) {
    Real& value=memo[mask];if(!std::isnan(value))return value;
    int a=first_bit(mask);auto rest=mask^(1U<<a);value=0;
    for(auto bits=rest;bits;bits&=bits-1) {int b=first_bit(bits);value+=kernel.at(occupied[a],occupied[b])*solve(rest^(1U<<b));}
    return value;
}
Real Hafnian::amplitude(){return solve((1U<<occupied.size())-1);}
Real Hafnian::moved(int particle,int destination) {
    auto rest=((1U<<occupied.size())-1)^(1U<<particle);Real value=0;
    for(auto bits=rest;bits;bits&=bits-1){int b=first_bit(bits);value+=kernel.at(destination,occupied[b])*solve(rest^(1U<<b));}
    return value;
}
Real Hafnian::moved_pair(int first,int second,int i,int j) {
    require(first!=second,"Pair move requires distinct source particles");
    const auto rest=((1U<<occupied.size())-1)^(1U<<first)^(1U<<second);
    Real value=kernel.at(i,j)*solve(rest);
    for(int a=0;a<int(occupied.size());++a) {
        if(!(rest&(1U<<a)))continue;
        for(int b=a+1;b<int(occupied.size());++b) {
            if(!(rest&(1U<<b)))continue;
            const Real pairing=kernel.at(i,occupied[a])*kernel.at(j,occupied[b])
                              +kernel.at(i,occupied[b])*kernel.at(j,occupied[a]);
            value+=pairing*solve(rest^(1U<<a)^(1U<<b));
        }
    }
    return value;
}
void Hafnian::accept_move(int particle,int destination) {
    occupied[particle]=destination;
    // Subsets excluding the changed particle retain their exact cached value.
    // Keep particle indices fixed so bit masks retain the same interpretation.
    const std::size_t bit=std::size_t(1)<<particle;
    for(std::size_t start=bit;start<memo.size();start+=2*bit)
        std::fill(memo.begin()+start,memo.begin()+start+bit,
                  std::numeric_limits<Real>::quiet_NaN());
    // Re-evaluate the full amplitude lazily using its original summation order.
}
Real Hafnian::move_change(int particle,int destination) const {
    if(!kernel.jastrow)return 0;
    Real d=0;int old=occupied[particle];
    for(int b:occupied)if(b!=old){int x=kernel.shell[old*kernel.n+b],y=kernel.shell[destination*kernel.n+b];if(x>=0)d-=kernel.v[x+3];if(y>=0)d+=kernel.v[y+3];}
    return d;
}
Vec Hafnian::log_derivative() {
    Vec result{};Real psi=amplitude();require(usable(psi),"Zero/nonfinite amplitude in logarithmic derivative");
    auto full=(1U<<occupied.size())-1;
    for(int a=0;a<int(occupied.size());++a)for(int b=a+1;b<int(occupied.size());++b) {
        Real minor=solve(full^(1U<<a)^(1U<<b))/psi;
        for(int p=0;p<3;++p)result[p]+=kernel.df[p][occupied[a]*kernel.n+occupied[b]]*minor;
    }
    if(kernel.jastrow){auto b=kernel.counts(occupied);for(int p=0;p<3;++p)result[p+3]=-b[p];}
    return result;
}
Energy local_energy(const Config& cfg,const Lattice& l,const std::vector<int>& c,Hafnian& h) {
    auto n=occupancy(l,c);std::vector<int> index(l.n,-1);
    for(int p=0;p<int(c.size());++p)index[c[p]]=p;
    Real psi=h.amplitude();require(usable(psi),"Zero/nonfinite sampled amplitude");Energy e;
    auto hop=[&](int a,int b)->Real {
        if(n[a]==n[b])return 0;
        if(n[b])std::swap(a,b);
        return h.moved(index[a],b)/psi*std::exp(-h.move_change(index[a],b));
    };
    for(auto b:l.bonds){e.interaction+=cfg.VB*n[b[0]]*n[b[1]];e.hopping-=cfg.JB*hop(b[0],b[1]);}
    // Eq. (1) directly: four oriented (delta_x,delta_y) per center, plus h.c.
    // hop() selects the single allowed direction; no additional factor of two.
    for(auto t:l.corners)if(n[t[0]])e.correlated-=cfg.K*hop(t[1],t[2]);
    require(std::isfinite(e.total()),"Nonfinite local energy");return e;
}
Diagnostics measure_diagnostics(const Lattice& l,const Kernel& k,const std::vector<int>& c,Hafnian& h) {
    Diagnostics d;auto n=occupancy(l,c);std::vector<int> index(l.n,-1);
    for(int p=0;p<int(c.size());++p){int i=c[p];index[i]=p;if((i%l.lx+i/l.lx)%2==0)++d.na;else ++d.nb;}
    const int translated=l.shift(0,1,0),far=l.site(l.lx/2,l.ly/2);
    d.density={Real(n[0]),Real(n[translated]),Real(n[far])};
    Real psi=h.amplitude();require(usable(psi),"Invalid diagnostic amplitude");
    auto single=[&](int i,int j)->Real {
        if(i==j)return n[i];
        if(n[i]==n[j])return 0;
        if(n[j])std::swap(i,j);
        return h.moved(index[i],j)/psi*std::exp(-h.move_change(index[i],j))/2; // (b_i^dagger b_j + h.c.)/2
    };
    auto directed_pair=[&](int i,int j,int a,int b)->Real {
        if(!n[a]||!n[b]||(n[i]&&i!=a&&i!=b)||(n[j]&&j!=a&&j!=b))return 0;
        auto moved=c;
        moved.erase(std::find(moved.begin(),moved.end(),a));
        moved.erase(std::find(moved.begin(),moved.end(),b));
        moved.push_back(i);moved.push_back(j);
        return h.moved_pair(index[a],index[b],i,j)/psi*std::exp(-k.change(c,moved));
    };
    for(int i=0;i<l.n;++i) {
        int j=l.shift(i,l.lx/2,l.ly/2),ip=l.shift(i,1,0),jp=l.shift(j,1,0);
        Real one=single(i,j),pair=(directed_pair(i,ip,j,jp)+directed_pair(j,jp,i,ip))/2;
        d.one[3]+=one/l.n;d.pair[3]+=pair/l.n;
        d.one_other+=single(i,l.shift(i,l.lx/2,l.ly/2+1))/l.n;
        if(i==0){d.one[0]=one;d.pair[0]=pair;}
        if(i==translated){d.one[1]=one;d.pair[1]=pair;}
    }
    d.one[2]=d.one[0]-d.one[1];d.pair[2]=d.pair[0]-d.pair[1];return d;
}
CorrelationAccumulator::CorrelationAccumulator(const Lattice& l):
    one(std::size_t(l.n)*l.n),pair(l.bonds.size()*l.bonds.size()),density(l.n) {}
void CorrelationAccumulator::add(const Lattice& l,const Kernel& k,const std::vector<int>& c,Hafnian& h,Real w,bool measure_pair) {
    auto n=occupancy(l,c);Real psi=h.amplitude();
    require(usable(psi),"Invalid amplitude in correlation measurement");
    weight+=w;
    for(int p=0;p<int(c.size());++p) {
        int j=c[p];one[j*l.n+j]+=w;
        for(int i=0;i<l.n;++i)if(!n[i]) {
            Real value=w*h.moved(p,i)/psi*std::exp(-h.move_change(p,i))/2;
            // For a real state, averaging O and O^dagger is unbiased.
            one[i*l.n+j]+=value;one[j*l.n+i]+=value;
        }
    }
    const std::size_t nb=l.bonds.size();
    if(measure_pair) {
        pair_weight+=w;
        std::vector<int> index(l.n,-1);
        for(int p=0;p<int(c.size());++p)index[c[p]]=p;
        for(std::size_t source=0;source<nb;++source) {
            int a=l.bonds[source][0],b=l.bonds[source][1];
            if(!n[a]||!n[b])continue;
            auto remaining=c;
            remaining.erase(std::find(remaining.begin(),remaining.end(),a));
            remaining.erase(std::find(remaining.begin(),remaining.end(),b));
            // Annihilate first, then create: includes coincident and overlapping bonds.
            for(std::size_t target=0;target<nb;++target) {
                int i=l.bonds[target][0],j=l.bonds[target][1];
                if((n[i]&&i!=a&&i!=b)||(n[j]&&j!=a&&j!=b))continue;
                auto moved=remaining;moved.push_back(i);moved.push_back(j);
                Real value=w*h.moved_pair(index[a],index[b],i,j)/psi*std::exp(-k.change(c,moved))/2;
                pair[target*nb+source]+=value;pair[source*nb+target]+=value;
            }
        }
    }
    for(int i:c)for(int j:c) {
        int r=l.site(j%l.lx-i%l.lx,j/l.lx-i/l.lx);density[r]+=w/l.n;
    }
}
void CorrelationAccumulator::write(const Lattice& l,const std::string& path) const {
    require(weight>0,"Empty correlation measurement");
    require(pair_weight>0,"No pair samples: reduce pair_measure_stride");
    auto parent=vmc_fs::path(path).parent_path();
    auto g=output((parent/"one_body_correlation.dat").string());
    g<<"# i j <b_i^dagger b_j>; site=x+Lx*y, zero-based; real Hermitian estimator\n";
    for(int i=0;i<l.n;++i)for(int j=0;j<l.n;++j)
        g<<i<<' '<<j<<' '<<one[i*l.n+j]/weight<<'\n';
    auto p=output((parent/"nn_pair_correlation.dat").string());
    p<<"# i j k l <b_i^dagger b_j^dagger b_k b_l>; zero-based; each undirected NN bond listed once (+x,+y)\n";
    std::size_t nb=l.bonds.size();
    for(std::size_t a=0;a<nb;++a)for(std::size_t b=0;b<nb;++b)
            p<<l.bonds[a][0]<<' '<<l.bonds[a][1]<<' '<<l.bonds[b][0]<<' '<<l.bonds[b][1]<<' '<<pair[a*nb+b]/pair_weight<<'\n';
    // Retain the original translation-averaged outputs for compatibility.
    auto f=output(path);f<<"dx,dy,G1,G2_xx,nn\n";
    for(int r=0;r<l.n;++r) {
        Real g1=0,g2=0;
        for(int i=0;i<l.n;++i) {
            int j=l.shift(i,r%l.lx,r/l.lx);
            g1+=one[i*l.n+j];g2+=pair[std::size_t(2*i)*nb+2*j];
        }
        f<<r%l.lx<<','<<r/l.lx<<','<<g1/(weight*l.n)<<','<<g2/(pair_weight*l.n)<<','<<density[r]/weight<<'\n';
    }
    auto sq=output(path+".structure.csv");sq<<"kx_index,ky_index,S_connected\n";
    Real pi=std::acos(Real(-1)),nu=density[0]/weight;
    for(int ky=0;ky<l.ly;++ky)for(int kx=0;kx<l.lx;++kx) {
        Real value=0;
        for(int r=0;r<l.n;++r)value+=std::cos(2*pi*(Real(kx)*(r%l.lx)/l.lx+Real(ky)*(r/l.lx)/l.ly))*(density[r]/weight-nu*nu);
        sq<<kx<<','<<ky<<','<<value<<'\n';
    }
    auto counts=output((parent/"correlation_counts.txt").string());
    counts<<"one_body_samples = "<<weight<<"\npair_samples = "<<pair_weight<<'\n';
}
Sampler::Sampler(const Config& c,const Lattice& l,const Kernel& k,std::uint64_t seed):cfg(c),lat(l),kernel(k),rng(seed) {
    std::vector<int> sites(l.n);std::iota(sites.begin(),sites.end(),0);
    for(int tries=0;tries<10000;++tries) {
        std::shuffle(sites.begin(),sites.end(),rng);occupied.assign(sites.begin(),sites.begin()+c.Nb);
        cached=std::make_unique<Hafnian>(k,occupied);
        if(usable(cached->amplitude()))return;
    }throw std::runtime_error("No nonzero initial configuration after 10000 attempts");
}
void Sampler::sweep() {
    static const int dx[]={1,-1,0,0},dy[]={0,0,1,-1};
    for(int attempt=0;attempt<cfg.Nb;++attempt) {
        ++proposed;int p=std::uniform_int_distribution<int>(0,cfg.Nb-1)(rng),j;
        if(uniform(rng)<cfg.local_probability) {
            int d=std::uniform_int_distribution<int>(0,3)(rng);j=lat.shift(occupied[p],dx[d],dy[d]);
            if(std::find(occupied.begin(),occupied.end(),j)!=occupied.end())continue;
        } else {
            do {j=std::uniform_int_distribution<int>(0,lat.n-1)(rng);}
            while(std::find(occupied.begin(),occupied.end(),j)!=occupied.end());
        }
        auto& h=*cached;
        Real old=h.amplitude(),next=h.moved(p,j);
        require(usable(old)&&std::isfinite(next),"Nonfinite Metropolis amplitude");
        if(next==0)continue;
        // Log acceptance avoids squaring a potentially huge ratio.
        Real log_accept=2*(std::log(std::abs(next))-std::log(std::abs(old)))-2*h.move_change(p,j);
        if(log_accept>=0||std::log(uniform(rng))<log_accept){
            h.accept_move(p,j);
            occupied[p]=j;
            ++accepted;
        }
    }
}
Statistics sample_local(const Config& c,const Lattice& l,const Kernel& k,std::uint64_t seed,
                  const std::string& trace_path,const std::string& corr_path,
                  CorrelationAccumulator* external=nullptr) {
    Sampler chain(c,l,k,seed);
    std::ofstream diagnostic;
    if(c.diagnostics&&!trace_path.empty()) {
        auto parent=vmc_fs::path(trace_path).parent_path();
        diagnostic=output((parent/"diagnostics.csv").string());
        diagnostic<<"phase,sample,sweep,E_total,n0,nx,nfar,density_diff_x,density_diff_far,NA,NB,NA_minus_NB,NA_minus_NB_sq,G1_far_0,G1_far_x,G1_far_diff,G1_far_avg,G1_other_avg,G2xx_far_0,G2xx_far_x,G2xx_far_diff,G2xx_far_avg,B1,B2,B3,J\n";
        auto info=output((parent/"diagnostics_info.txt").string());
        info<<"Lx = "<<l.lx<<"\nLy = "<<l.ly<<"\nNb = "<<c.Nb
            <<"\nfar_dx = "<<l.lx/2<<"\nfar_dy = "<<l.ly/2
            <<"\nother_dx = "<<l.lx/2<<"\nother_dy = "<<(l.ly/2+1)%l.ly
            <<"\nsite_0 = 0\nsite_x = "<<l.shift(0,1,0)<<"\nsite_far = "<<l.site(l.lx/2,l.ly/2)
            <<"\nmeasurement_interval_sweeps = "<<std::int64_t(c.sweeps_between)*c.diagnostic_stride
            <<"\nmeasurement_seed = "<<seed
            <<"\nCoordinates: zero-based site=x+Lx*y. A: (x+y)%2=0.\n"
            <<"NA/NB are a true bipartition for even Lx and Ly; odd periodic sizes only have checkerboard labels.\n"
            <<"one/pair use real Hermitian estimators, identical to full correlation outputs.\n"
            <<"G2xx: pair P_i=b_i*b_(i+x); displacement far joins pair origins.\n"
            <<"_x translates BOTH endpoints/bonds by +x; _diff is _0 minus _x; _avg averages all origins.\n"
            <<"phase=warmup includes initial configuration at sweep=0; exclude warmup from stationary statistics.\n"
            <<"These are instantaneous estimators, not running averages. Only final fixed-parameter chain is recorded.\n";
    }
    auto record_diagnostic=[&](const char* phase,int sample,std::int64_t sweep,Hafnian& h,Real energy) {
        if(!diagnostic.is_open())return;
        auto d=measure_diagnostics(l,k,chain.occupied,h);int imbalance=d.na-d.nb;
        diagnostic<<phase<<','<<sample<<','<<sweep<<','<<energy
            <<','<<d.density[0]<<','<<d.density[1]<<','<<d.density[2]
            <<','<<d.density[0]-d.density[1]<<','<<d.density[0]-d.density[2]
            <<','<<d.na<<','<<d.nb<<','<<imbalance<<','<<imbalance*imbalance;
        for(auto v:d.one)diagnostic<<','<<v;
        diagnostic<<','<<d.one_other;
        for(auto v:d.pair)diagnostic<<','<<v;
        auto bc=k.counts(chain.occupied);Real jc=0;for(int p=0;p<3;++p){diagnostic<<','<<bc[p];jc+=bc[p]*k.v[p+3];}diagnostic<<','<<jc;
        diagnostic<<'\n';
    };
    if(diagnostic.is_open()){auto& h=chain.current_hafnian();record_diagnostic("warmup",0,0,h,local_energy(c,l,chain.occupied,h).total());}
    const std::int64_t interval=std::int64_t(c.sweeps_between)*c.diagnostic_stride;
    for(int i=0;i<c.warmup_sweeps;++i) {
        chain.sweep();
        if(diagnostic.is_open()&&(i+1)%interval==0){auto& h=chain.current_hafnian();record_diagnostic("warmup",0,i+1,h,local_energy(c,l,chain.occupied,h).total());}
    }
    if(diagnostic.is_open())diagnostic.flush();
    chain.proposed=chain.accepted=0;
    std::ofstream trace;if(!trace_path.empty()){trace=output(trace_path);trace<<"sample,E_total,E_J,E_V,E_K,O_chi1,O_Delta,O_delta,O_v1,O_v2,O_v3\n";}
    Statistics s;Real m2=0,block_sum=0;std::vector<Real> blocks;
    // Online centered covariances reduce cancellation in nearly redundant directions.
    Vec cov_oe{};Mat cov_oo{};
    // Allocate dense correlation matrices only for the final requested measurement.
    std::unique_ptr<CorrelationAccumulator> owned;
    if(!corr_path.empty()&&!external)owned=std::make_unique<CorrelationAccumulator>(l);
    auto corr=external?external:owned.get();
    if(corr)require(c.samples>=c.pair_measure_stride,"Each chain needs at least one pair sample; reduce pair_measure_stride");
    for(int t=1;t<=c.samples;++t) {
        for(int i=0;i<c.sweeps_between;++i)chain.sweep();
        auto& h=chain.current_hafnian();auto e=local_energy(c,l,chain.occupied,h);auto o=h.log_derivative();Real E=e.total();
        Real de=E-s.energy;Vec d{};for(int p=0;p<6;++p)d[p]=o[p]-s.mean_o[p];
        s.energy+=de/t;for(int p=0;p<6;++p)s.mean_o[p]+=d[p]/t;
        m2+=de*(E-s.energy);
        for(int p=0;p<6;++p){cov_oe[p]+=d[p]*(E-s.energy);for(int q=0;q<6;++q)cov_oo[p][q]+=d[p]*(o[q]-s.mean_o[q]);}
        s.components.interaction+=e.interaction/c.samples;s.components.hopping+=e.hopping/c.samples;s.components.correlated+=e.correlated/c.samples;
        block_sum+=E;if(t%c.block_size==0){blocks.push_back(block_sum/c.block_size);block_sum=0;}
        if(trace)trace<<t<<','<<E<<','<<e.hopping<<','<<e.interaction<<','<<e.correlated<<','<<o[0]<<','<<o[1]<<','<<o[2]<<','<<o[3]<<','<<o[4]<<','<<o[5]<<'\n';
        if(corr)corr->add(l,k,chain.occupied,h,1,t%c.pair_measure_stride==0);
        if(t%c.diagnostic_stride==0)record_diagnostic("measurement",t,std::int64_t(c.warmup_sweeps)+std::int64_t(t)*c.sweeps_between,h,E);
        if(diagnostic.is_open()&&t%1000==0)diagnostic.flush();
    }
    s.variance=m2/c.samples;s.acceptance=Real(chain.accepted)/chain.proposed;
    Real bvar=0;for(Real b:blocks)bvar+=(b-s.energy)*(b-s.energy);
    s.error=std::sqrt(bvar/(blocks.size()*(blocks.size()-1)));
    for(int p=0;p<6;++p){s.force[p]=cov_oe[p]/c.samples;for(int q=0;q<6;++q)s.metric[p][q]=cov_oo[p][q]/c.samples;}
    if(corr&&!corr_path.empty())corr->write(l,corr_path);
    return s;
}
Statistics merge_statistics(const std::vector<Statistics>& parts,const std::vector<int>& counts) {
    require(!parts.empty()&&parts.size()==counts.size(),"Invalid chain statistics");
    Real total=0;for(int n:counts){require(n>0,"Empty chain");total+=n;}
    Statistics out;
    for(std::size_t r=0;r<parts.size();++r) {
        const auto& s=parts[r];Real w=counts[r]/total;
        out.energy+=w*s.energy;
        for(int p=0;p<6;++p)out.mean_o[p]+=w*s.mean_o[p];
        out.acceptance+=w*s.acceptance;
        out.components.hopping+=w*s.components.hopping;
        out.components.interaction+=w*s.components.interaction;
        out.components.correlated+=w*s.components.correlated;
        out.error+=w*w*s.error*s.error;
    }
    out.error=std::sqrt(out.error);
    for(std::size_t r=0;r<parts.size();++r) {
        const auto& s=parts[r];Real w=counts[r]/total,de=s.energy-out.energy;
        out.variance+=w*(s.variance+de*de);
        for(int p=0;p<6;++p) {
            Real dp=s.mean_o[p]-out.mean_o[p];
            out.force[p]+=w*(s.force[p]+dp*de);
            for(int q=0;q<6;++q)out.metric[p][q]+=w*(s.metric[p][q]+dp*(s.mean_o[q]-out.mean_o[q]));
        }
    }
    return out;
}
Statistics sample(const Config& c,const Lattice& l,const Kernel& k,std::uint64_t seed,
                  const std::string& trace_path,const std::string& corr_path) {
    const int size=world_size();
    if(size==1)return sample_local(c,l,k,seed,trace_path,corr_path);
#ifdef VMC_MPI
    const int rank=rank_id();
    // Distribute whole blocks so no block straddles chains or drops observations.
    const int blocks=c.samples/c.block_size;
    require(blocks/size>=2,"MPI requires at least two full blocks per chain; increase samples or reduce ranks/block_size");
    std::vector<int> counts(size);
    for(int r=0;r<size;++r)counts[r]=(blocks/size+(r<blocks%size))*c.block_size;
    Config local=c;local.samples=counts[rank];
    std::string local_trace;
    const auto local_seed=chain_seed(seed,rank);
    if(!trace_path.empty()) {
        std::ostringstream name;name<<"chain_"<<std::setw(4)<<std::setfill('0')<<rank;
        auto dir=vmc_fs::path(trace_path).parent_path()/"chains"/name.str();
        vmc_fs::create_directories(dir);local_trace=(dir/"samples.csv").string();
        auto meta=output((dir/"chain_info.txt").string());
        meta<<"rank = "<<rank<<"\nmpi_ranks = "<<size<<"\nseed = "<<local_seed
            <<"\nsamples = "<<local.samples<<"\ntotal_samples = "<<c.samples
            <<"\npair_measure_stride = "<<c.pair_measure_stride
            <<"\npair_samples = "<<(c.correlations?local.samples/c.pair_measure_stride:0)<<"\n";
    }
    std::unique_ptr<CorrelationAccumulator> corr;
    if(!corr_path.empty())corr=std::make_unique<CorrelationAccumulator>(l);
    auto s=sample_local(local,l,k,local_seed,local_trace,"",corr.get());
    // Explicit packing avoids struct padding and preserves long-double precision.
    auto pack=[](const Statistics& v) {
        std::array<Real,55> a{};int j=0;
        for(Real x:{v.energy,v.error,v.variance,v.acceptance,v.components.hopping,v.components.interaction,v.components.correlated})a[j++]=x;
        for(Real x:v.mean_o)a[j++]=x;
        for(Real x:v.force)a[j++]=x;
        for(const auto& row:v.metric)for(Real x:row)a[j++]=x;
        return a;
    };
    auto packed=pack(s);std::vector<Real> gathered(rank==0?size*55:0);
    MPI_Gather(packed.data(),55,MPI_LONG_DOUBLE,gathered.data(),55,MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
    Statistics merged;
    if(rank==0) {
        std::vector<Statistics> parts(size);
        for(int r=0;r<size;++r) {
            auto& v=parts[r];int j=r*55;
            v.energy=gathered[j++];v.error=gathered[j++];v.variance=gathered[j++];v.acceptance=gathered[j++];
            v.components.hopping=gathered[j++];v.components.interaction=gathered[j++];v.components.correlated=gathered[j++];
            for(auto& x:v.mean_o)x=gathered[j++];
            for(auto& x:v.force)x=gathered[j++];
            for(auto& row:v.metric)for(auto& x:row)x=gathered[j++];
        }
        merged=merge_statistics(parts,counts);
        if(!trace_path.empty()) {
            auto stats=output((vmc_fs::path(trace_path).parent_path()/"pooled_statistics.csv").string());
            stats<<"E,error,variance,acceptance,E_J,E_V,E_K";for(int p=0;p<6;++p)stats<<",O"<<p;for(int p=0;p<6;++p)stats<<",F"<<p;for(int p=0;p<6;++p)for(int q=0;q<6;++q)stats<<",S"<<p<<q;stats<<'\n';
            auto values=pack(merged);
            for(std::size_t j=0;j<values.size();++j)stats<<(j?",":"")<<values[j];
            stats<<'\n';
            auto info=output((vmc_fs::path(trace_path).parent_path()/"chain_summary.csv").string());
            info<<"rank,seed,samples,E_total,blocked_error,acceptance\n";
            for(int r=0;r<size;++r)info<<r<<','<<chain_seed(seed,r)<<','<<counts[r]<<','<<parts[r].energy<<','<<parts[r].error<<','<<parts[r].acceptance<<'\n';
        }
    }
    if(corr) {
        for(auto* array:{&corr->one,&corr->pair,&corr->density}) {
            require(array->size()<=std::size_t(std::numeric_limits<int>::max()),"Correlation MPI count overflow");
            if(rank==0)MPI_Reduce(MPI_IN_PLACE,array->data(),int(array->size()),MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            else MPI_Reduce(array->data(),nullptr,int(array->size()),MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        }
        Real weight=0;MPI_Reduce(&corr->weight,&weight,1,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        Real pair_weight=0;
        MPI_Reduce(&corr->pair_weight,&pair_weight,1,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        if(rank==0){corr->weight=weight;corr->pair_weight=pair_weight;corr->write(l,corr_path);}
    }
    return merged; // only root consumes the pooled statistics
#else
    return {}; // serial world_size() is always one
#endif
}
Vec sr_direction(const Statistics& s,Real step,Real shift) {
    Mat a=s.metric;Vec b{};
    for(int i=0;i<6;++i){a[i][i]+=shift;b[i]=-step*s.force[i];}
    for(int i=0;i<6;++i) {
        int pivot=i;for(int j=i+1;j<6;++j)if(std::abs(a[j][i])>std::abs(a[pivot][i]))pivot=j;
        std::swap(a[i],a[pivot]);std::swap(b[i],b[pivot]);require(usable(a[i][i]),"Singular SR matrix");
        for(int j=i+1;j<6;++j){Real factor=a[j][i]/a[i][i];for(int q=i;q<6;++q)a[j][q]-=factor*a[i][q];b[j]-=factor*b[i];}
    }
    Vec d{};for(int i=5;i>=0;--i){Real v=b[i];for(int j=i+1;j<6;++j)v-=a[i][j]*d[j];d[i]=v/a[i][i];require(std::isfinite(d[i]),"Nonfinite SR update");}return d;
}
Vec sr_parameters(const Config& c,const Vec& alpha,const Statistics& s,Real& scale) {
    Statistics transformed=s;
    Vec jac{1,1,alpha[2],1,1,1}; // O_logdelta = delta * O_delta
    for(int p=0;p<6;++p) {
        transformed.force[p]=jac[p]*s.force[p];
        for(int q=0;q<6;++q)transformed.metric[p][q]=jac[p]*jac[q]*s.metric[p][q];
    }
    if(!c.optimize_delta) {
        // Solve the restricted two-parameter problem, not the full solve with d3 discarded.
        for(int p=0;p<6;++p)transformed.metric[p][2]=transformed.metric[2][p]=0;
        transformed.force[2]=0;
    }
    if(!c.Ifjastrow)for(int p=3;p<6;++p){transformed.force[p]=0;for(int q=0;q<6;++q)transformed.metric[p][q]=transformed.metric[q][p]=0;}
    Vec direction=sr_direction(transformed,c.sr_step,c.sr_shift);
    scale=1;
    if(std::abs(direction[2])>c.max_log_delta_step)scale=c.max_log_delta_step/std::abs(direction[2]);
    Lattice l(c);
    for(int retry=0;retry<30;++retry,scale*=0.5L) {
        Vec next{alpha[0]+scale*direction[0],alpha[1]+scale*direction[1],alpha[2]};
        for(int p=3;p<6;++p)next[p]=c.Ifjastrow?alpha[p]+scale*direction[p]:0;
        if(c.optimize_delta)next[2]=alpha[2]*std::exp(scale*direction[2]);
        try {Kernel trial(c,l,next);return next;}catch(const std::runtime_error&){}
    }
    throw std::runtime_error("Nonfinite/zero pair state in SR update after 30 backtracks");
}
namespace {
void save_parameters(const std::string& path,const Config& c,const Vec& a) {
    auto f=output(path);
    f<<std::setprecision(std::numeric_limits<Real>::max_digits10);
    f<<"# Physical coordinates: chi1, Delta, delta. SR internally uses log(delta).\n"
     <<"Lx = "<<c.Lx<<"\nLy = "<<c.Ly<<"\nNb = "<<c.Nb<<"\nJB = "<<c.JB<<"\nVB = "<<c.VB<<"\nK = "<<c.K
     <<"\nchi1 = "<<a[0]<<"\nDelta = "<<a[1]<<"\ndelta = "<<a[2]
     <<"\nIfjastrow = "<<c.Ifjastrow<<"\nv1 = "<<a[3]<<"\nv2 = "<<a[4]<<"\nv3 = "<<a[5]
     <<"\noptimize_delta = "<<c.optimize_delta<<"\nmax_log_delta_step = "<<c.max_log_delta_step
     <<"\nwarmup_sweeps = "<<c.warmup_sweeps<<"\nsamples_sr = "<<c.sr_samples()<<"\nsamples = "<<c.samples<<"\nsweeps_between = "<<c.sweeps_between
     <<"\nblock_size = "<<c.block_size<<"\nseed = "<<c.seed<<"\nlocal_probability = "<<c.local_probability
     <<"\nsr_iterations = 0\nsr_step = "<<c.sr_step<<"\nsr_shift = "<<c.sr_shift
     <<"\npair_measure_stride = "<<c.pair_measure_stride
     <<"\ndiagnostics = "<<c.diagnostics<<"\ndiagnostic_stride = "<<c.diagnostic_stride
     <<"\ncorrelations = "<<c.correlations<<"\noutput = "<<c.output<<"_continued\n";
}
}
void run(const Config& c) {
    c.validate();Lattice l(c);Vec alpha=c.alpha;if(!c.Ifjastrow)for(int p=3;p<6;++p)alpha[p]=0;Kernel initial(c,l,alpha);
    const bool root=rank_id()==0;
    if(world_size()>1)require(c.samples/c.block_size/world_size()>=2,"MPI requires at least two blocks per chain");
    Config sr_config=c;sr_config.samples=c.sr_samples();
    if(c.sr_iterations>0&&world_size()>1)require(sr_config.samples/c.block_size/world_size()>=2,"samples_sr requires at least two blocks per MPI chain");
    if(c.correlations)require((c.samples/c.block_size/world_size())*c.block_size>=c.pair_measure_stride,
        "Each MPI chain needs at least one complete pair measurement interval");
    std::ofstream log;
    if(root) {
    require(!vmc_fs::exists(c.output),"Output path already exists; choose a new output in .inp: "+c.output);
    vmc_fs::create_directories(c.output);
    save_parameters(c.output+"/checkpoint.inp",c,alpha);
    log=output(c.output+"/optimization.csv");
    log<<std::setprecision(std::numeric_limits<Real>::max_digits10);
    log<<"iteration,chi1,Delta,delta,lambda_c,lambda_eff,E_total,blocked_error,acceptance,min_margin,min_gap,step_scale,delta_log_change,boundary_modes,v1,v2,v3";
    for(const char* name:{"chi1","Delta","delta","v1","v2","v3"})log<<",F_"<<name;
    for(int p=0;p<6;++p)for(int q=p;q<6;++q)log<<",S"<<p<<q;
    log<<'\n';
    }
    synchronize();
    if(root)std::cout<<"Samples: SR="<<c.sr_samples()<<", final="<<c.samples<<" (total across chains)"<<std::endl;
    std::string stop_reason="requested_iterations_completed";int completed=0;
    for(int it=0;it<c.sr_iterations;++it) {
        Kernel k(c,l,alpha);auto s=sample(sr_config,l,k,c.seed+std::uint64_t(it));Real scale=1;Vec next{};int update_ok=1;
        if(root) {
        try {next=sr_parameters(c,alpha,s,scale);}
        catch(const std::runtime_error& e) {
            stop_reason="update_failed_using_last_valid_state";
            std::cerr<<"SR stopped: "<<e.what()<<". Measuring the checkpoint state.\n";update_ok=0;
        }
        }
#ifdef VMC_MPI
        MPI_Bcast(&update_ok,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(next.data(),6,MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
#endif
        if(!update_ok){stop_reason="update_failed_using_last_valid_state";break;}
        if(root) {
        log<<it<<','<<alpha[0]<<','<<alpha[1]<<','<<alpha[2]<<','<<k.lambda_c<<','<<k.lambda_eff
           <<','<<s.energy<<','<<s.error<<','<<s.acceptance<<','<<k.min_margin<<','<<k.min_gap<<','<<scale
           <<','<<std::log(next[2]/alpha[2])<<','<<k.boundary_modes;
        for(int p=3;p<6;++p)log<<','<<alpha[p];
        for(auto v:s.force)log<<','<<v;
        for(int p=0;p<6;++p)for(int q=p;q<6;++q)log<<','<<s.metric[p][q];
        log<<'\n';log.flush();
        }
        alpha=next;++completed;
        if(root) {
        save_parameters(c.output+"/checkpoint.inp",c,alpha);
        std::cout<<"SR "<<it<<" E/Ns="<<s.energy/l.n<<" +/- "<<s.error/l.n<<" delta_next="<<alpha[2]<<" scale="<<scale<<std::endl;
        }
    }
    Kernel final_kernel(c,l,alpha);
    if(root)std::cout<<"Final measurement: "<<c.samples<<" total samples on "<<world_size()<<" chain(s), delta="<<alpha[2]<<std::endl;
    auto s=sample(c,l,final_kernel,c.seed+std::uint64_t(completed),c.output+"/samples.csv",c.correlations?c.output+"/correlations.csv":"");
    if(!root)return;
    auto summary=output(c.output+"/summary.txt");
    summary<<"Ifjastrow = "<<c.Ifjastrow<<"\nv1 = "<<alpha[3]<<"\nv2 = "<<alpha[4]<<"\nv3 = "<<alpha[5]<<"\n";
    summary<<"mpi_ranks = "<<world_size()<<"\nsamples_scope = total_across_chains\n";
    summary<<"Lx = "<<c.Lx<<"\nLy = "<<c.Ly<<"\nNb = "<<c.Nb<<"\nnu = "<<Real(c.Nb)/l.n
        <<"\nJB = "<<c.JB<<"\nVB = "<<c.VB<<"\nK = "<<c.K<<"\nchi1 = "<<alpha[0]<<"\nDelta = "<<alpha[1]<<"\ndelta = "<<alpha[2]
        <<"\nlambda_c = "<<final_kernel.lambda_c<<"\nlambda_eff = "<<final_kernel.lambda_eff
        <<"\noptimize_delta = "<<c.optimize_delta<<"\nmax_log_delta_step = "<<c.max_log_delta_step
        <<"\nE_total = "<<s.energy<<"\nE_per_site = "<<s.energy/l.n<<"\nblocked_error_total = "<<s.error<<"\nblocked_error_per_site = "<<s.error/l.n
        <<"\nlocal_energy_variance = "<<s.variance<<"\nE_J = "<<s.components.hopping<<"\nE_V = "<<s.components.interaction<<"\nE_K = "<<s.components.correlated
        <<"\nacceptance = "<<s.acceptance<<"\nmin_A_minus_absB = "<<final_kernel.min_margin<<"\nmin_Ek = "<<final_kernel.min_gap
        <<"\nboundary_modes = "<<final_kernel.boundary_modes
        <<"\nseed = "<<c.seed<<"\nmeasurement_seed = "<<c.seed+std::uint64_t(completed)<<"\nsamples_sr = "<<c.sr_samples()<<"\nsamples = "<<c.samples<<"\nblock_size = "<<c.block_size<<"\nsweeps_between = "<<c.sweeps_between
        <<"\nwarmup_sweeps = "<<c.warmup_sweeps<<"\nlocal_probability = "<<c.local_probability<<"\nsr_iterations = "<<c.sr_iterations
        <<"\nsr_iterations_completed = "<<completed<<"\nsr_stop_reason = "<<stop_reason
        <<"\nsr_step = "<<c.sr_step<<"\nsr_shift = "<<c.sr_shift<<"\ncorrelations = "<<c.correlations
        <<"\npair_measure_stride = "<<c.pair_measure_stride
     <<"\ndiagnostics = "<<c.diagnostics<<"\ndiagnostic_stride = "<<c.diagnostic_stride
        <<"\ninitial_chi1 = "<<c.alpha[0]<<"\ninitial_Delta = "<<c.alpha[1]<<"\ninitial_delta = "<<c.alpha[2]
        <<"\ninitial_lambda_eff = "<<initial.lambda_eff
        <<"\nError is a block-mean estimate; check block-size saturation and independent seeds.\n";
    save_parameters(c.output+"/final_parameters.inp",c,alpha);
    std::cout<<std::setprecision(12)<<"Final E/Ns = "<<s.energy/l.n<<" +/- "<<s.error/l.n<<"; acceptance = "<<s.acceptance<<"\nOutput: "<<c.output<<std::endl;
}
} // namespace vmc
