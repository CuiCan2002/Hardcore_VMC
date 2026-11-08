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

namespace vmc {
namespace {
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
        else if(key=="samples")parse(value,c.samples);
        else if(key=="sweeps_between")parse(value,c.sweeps_between);
        else if(key=="block_size")parse(value,c.block_size);
        else if(key=="sr_iterations")parse(value,c.sr_iterations);
        else if(key=="sr_step")parse(value,c.sr_step);
        else if(key=="sr_shift")parse(value,c.sr_shift);
        else if(key=="max_log_delta_step")parse(value,c.max_log_delta_step);
        else if(key=="optimize_delta") {int b;parse(value,b);require(b==0||b==1,"optimize_delta must be 0 or 1");c.optimize_delta=b!=0;}
        else if(key=="seed") {require(!value.empty()&&value[0]!='-',"seed must be nonnegative");parse(value,c.seed);}
        else if(key=="correlations") {int b;parse(value,b);require(b==0||b==1,"correlations must be 0 or 1");c.correlations=b!=0;}
        else if(key=="output")c.output=value;
        else throw std::runtime_error("Unknown input key: "+key);
    }
    c.validate();return c;
}
void Config::validate() const {
    require(Lx>=3&&Ly>=3&&Lx<=128&&Ly<=128&&Lx*Ly<=4096,"Require Lx,Ly >= 3 and Ns <= 4096 (nondegenerate PBC)");
    require(Nb>=2&&Nb<=20&&Nb%2==0&&Nb<Lx*Ly,"Require even 2 <= Nb <= 20 and Nb < Ns (exact Hafnian limit)");
    for(Real x:{JB,VB,K,alpha[0],alpha[1],alpha[2],local_probability,sr_step,sr_shift,max_log_delta_step})
        require(std::isfinite(x),"All real parameters must be finite");
    require(alpha[1]!=0&&(VB-4*K)!=0,"Bk identically zero: pair-projected state vanishes");
    require(alpha[2]>0,"delta must be strictly positive");
    require(max_log_delta_step>0,"max_log_delta_step must be positive");
    require(local_probability>=0&&local_probability<1,"Require 0 <= local_probability < 1 to retain global moves");
    require(warmup_sweeps>=0&&samples>=2&&sweeps_between>=1,"Invalid sampling controls");
    require(block_size>=1&&samples%block_size==0&&samples/block_size>=2,"samples must contain at least two complete blocks");
    require(sr_iterations>=0&&sr_step>0&&sr_shift>0,"Invalid SR controls");
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
Vec Hafnian::log_derivative() {
    Vec result{};Real psi=amplitude();require(usable(psi),"Zero/nonfinite amplitude in logarithmic derivative");
    auto full=(1U<<occupied.size())-1;
    for(int a=0;a<int(occupied.size());++a)for(int b=a+1;b<int(occupied.size());++b) {
        Real minor=solve(full^(1U<<a)^(1U<<b))/psi;
        for(int p=0;p<3;++p)result[p]+=kernel.df[p][occupied[a]*kernel.n+occupied[b]]*minor;
    }return result;
}
Energy local_energy(const Config& cfg,const Lattice& l,const std::vector<int>& c,Hafnian& h) {
    auto n=occupancy(l,c);std::vector<int> index(l.n,-1);
    for(int p=0;p<int(c.size());++p)index[c[p]]=p;
    Real psi=h.amplitude();require(usable(psi),"Zero/nonfinite sampled amplitude");Energy e;
    auto hop=[&](int a,int b)->Real {
        if(n[a]==n[b])return 0;
        if(n[b])std::swap(a,b);
        return h.moved(index[a],b)/psi;
    };
    for(auto b:l.bonds){e.interaction+=cfg.VB*n[b[0]]*n[b[1]];e.hopping-=cfg.JB*hop(b[0],b[1]);}
    // Eq. (1) directly: four oriented (delta_x,delta_y) per center, plus h.c.
    // hop() selects the single allowed direction; no additional factor of two.
    for(auto t:l.corners)if(n[t[0]])e.correlated-=cfg.K*hop(t[1],t[2]);
    require(std::isfinite(e.total()),"Nonfinite local energy");return e;
}
CorrelationAccumulator::CorrelationAccumulator(const Lattice& l):
    one(std::size_t(l.n)*l.n),pair(l.bonds.size()*l.bonds.size()),density(l.n) {}
void CorrelationAccumulator::add(const Lattice& l,const Kernel& k,const std::vector<int>& c,Hafnian& h,Real w) {
    auto n=occupancy(l,c);Real psi=h.amplitude();
    require(usable(psi),"Invalid amplitude in correlation measurement");
    weight+=w;
    for(int p=0;p<int(c.size());++p) {
        int j=c[p];one[j*l.n+j]+=w;
        for(int i=0;i<l.n;++i)if(!n[i]) {
            Real value=w*h.moved(p,i)/psi/2;
            // For a real state, averaging O and O^dagger is unbiased.
            one[i*l.n+j]+=value;one[j*l.n+i]+=value;
        }
    }
    const std::size_t nb=l.bonds.size();
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
            Hafnian next(k,moved);Real value=w*next.amplitude()/psi/2;
            pair[target*nb+source]+=value;pair[source*nb+target]+=value;
        }
    }
    for(int i:c)for(int j:c) {
        int r=l.site(j%l.lx-i%l.lx,j/l.lx-i/l.lx);density[r]+=w/l.n;
    }
}
void CorrelationAccumulator::write(const Lattice& l,const std::string& path) const {
    require(weight>0,"Empty correlation measurement");
    auto parent=vmc_fs::path(path).parent_path();
    auto g=output((parent/"one_body_correlation.dat").string());
    g<<"# i j <b_i^dagger b_j>; site=x+Lx*y, zero-based; real Hermitian estimator\n";
    for(int i=0;i<l.n;++i)for(int j=0;j<l.n;++j)
        g<<i<<' '<<j<<' '<<one[i*l.n+j]/weight<<'\n';
    auto p=output((parent/"nn_pair_correlation.dat").string());
    p<<"# i j k l <b_i^dagger b_j^dagger b_k b_l>; zero-based; each undirected NN bond listed once (+x,+y)\n";
    std::size_t nb=l.bonds.size();
    for(std::size_t a=0;a<nb;++a)for(std::size_t b=0;b<nb;++b)
        p<<l.bonds[a][0]<<' '<<l.bonds[a][1]<<' '<<l.bonds[b][0]<<' '<<l.bonds[b][1]<<' '<<pair[a*nb+b]/weight<<'\n';
    // Retain the original translation-averaged outputs for compatibility.
    auto f=output(path);f<<"dx,dy,G1,G2_xx,nn\n";
    for(int r=0;r<l.n;++r) {
        Real g1=0,g2=0;
        for(int i=0;i<l.n;++i) {
            int j=l.shift(i,r%l.lx,r/l.lx);
            g1+=one[i*l.n+j];g2+=pair[std::size_t(2*i)*nb+2*j];
        }
        f<<r%l.lx<<','<<r/l.lx<<','<<g1/(weight*l.n)<<','<<g2/(weight*l.n)<<','<<density[r]/weight<<'\n';
    }
    auto sq=output(path+".structure.csv");sq<<"kx_index,ky_index,S_connected\n";
    Real pi=std::acos(Real(-1)),nu=density[0]/weight;
    for(int ky=0;ky<l.ly;++ky)for(int kx=0;kx<l.lx;++kx) {
        Real value=0;
        for(int r=0;r<l.n;++r)value+=std::cos(2*pi*(Real(kx)*(r%l.lx)/l.lx+Real(ky)*(r/l.lx)/l.ly))*(density[r]/weight-nu*nu);
        sq<<kx<<','<<ky<<','<<value<<'\n';
    }
}
Sampler::Sampler(const Config& c,const Lattice& l,const Kernel& k,std::uint64_t seed):cfg(c),lat(l),kernel(k),rng(seed) {
    std::vector<int> sites(l.n);std::iota(sites.begin(),sites.end(),0);
    for(int tries=0;tries<10000;++tries) {
        std::shuffle(sites.begin(),sites.end(),rng);occupied.assign(sites.begin(),sites.begin()+c.Nb);
        Hafnian h(k,occupied);if(usable(h.amplitude()))return;
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
        Hafnian h(kernel,occupied);Real old=h.amplitude(),next=h.moved(p,j);
        require(usable(old)&&std::isfinite(next),"Nonfinite Metropolis amplitude");
        if(next==0)continue;
        // Log acceptance avoids squaring a potentially huge ratio.
        Real log_accept=2*(std::log(std::abs(next))-std::log(std::abs(old)));
        if(log_accept>=0||std::log(uniform(rng))<log_accept){occupied[p]=j;++accepted;}
    }
}
Statistics sample(const Config& c,const Lattice& l,const Kernel& k,std::uint64_t seed,
                  const std::string& trace_path,const std::string& corr_path) {
    Sampler chain(c,l,k,seed);for(int i=0;i<c.warmup_sweeps;++i)chain.sweep();
    chain.proposed=chain.accepted=0;
    std::ofstream trace;if(!trace_path.empty()){trace=output(trace_path);trace<<"sample,E_total,E_J,E_V,E_K,O_chi1,O_Delta,O_delta\n";}
    Statistics s;Real m2=0,block_sum=0;std::vector<Real> blocks;
    // Online centered covariances reduce cancellation in nearly redundant directions.
    Vec cov_oe{};Mat cov_oo{};
    // Allocate dense correlation matrices only for the final requested measurement.
    std::unique_ptr<CorrelationAccumulator> corr;
    if(!corr_path.empty())corr=std::make_unique<CorrelationAccumulator>(l);
    for(int t=1;t<=c.samples;++t) {
        for(int i=0;i<c.sweeps_between;++i)chain.sweep();
        Hafnian h(k,chain.occupied);auto e=local_energy(c,l,chain.occupied,h);auto o=h.log_derivative();Real E=e.total();
        Real de=E-s.energy;Vec d{};for(int p=0;p<3;++p)d[p]=o[p]-s.mean_o[p];
        s.energy+=de/t;for(int p=0;p<3;++p)s.mean_o[p]+=d[p]/t;
        m2+=de*(E-s.energy);
        for(int p=0;p<3;++p){cov_oe[p]+=d[p]*(E-s.energy);for(int q=0;q<3;++q)cov_oo[p][q]+=d[p]*(o[q]-s.mean_o[q]);}
        s.components.interaction+=e.interaction/c.samples;s.components.hopping+=e.hopping/c.samples;s.components.correlated+=e.correlated/c.samples;
        block_sum+=E;if(t%c.block_size==0){blocks.push_back(block_sum/c.block_size);block_sum=0;}
        if(trace)trace<<t<<','<<E<<','<<e.hopping<<','<<e.interaction<<','<<e.correlated<<','<<o[0]<<','<<o[1]<<','<<o[2]<<'\n';
        if(corr)corr->add(l,k,chain.occupied,h);
    }
    s.variance=m2/c.samples;s.acceptance=Real(chain.accepted)/chain.proposed;
    Real bvar=0;for(Real b:blocks)bvar+=(b-s.energy)*(b-s.energy);
    s.error=std::sqrt(bvar/(blocks.size()*(blocks.size()-1)));
    for(int p=0;p<3;++p){s.force[p]=cov_oe[p]/c.samples;for(int q=0;q<3;++q)s.metric[p][q]=cov_oo[p][q]/c.samples;}
    if(corr)corr->write(l,corr_path);
    return s;
}
Vec sr_direction(const Statistics& s,Real step,Real shift) {
    Mat a=s.metric;Vec b{};
    for(int i=0;i<3;++i){a[i][i]+=shift;b[i]=-step*s.force[i];}
    for(int i=0;i<3;++i) {
        int pivot=i;for(int j=i+1;j<3;++j)if(std::abs(a[j][i])>std::abs(a[pivot][i]))pivot=j;
        std::swap(a[i],a[pivot]);std::swap(b[i],b[pivot]);require(usable(a[i][i]),"Singular SR matrix");
        for(int j=i+1;j<3;++j){Real factor=a[j][i]/a[i][i];for(int q=i;q<3;++q)a[j][q]-=factor*a[i][q];b[j]-=factor*b[i];}
    }
    Vec d{};for(int i=2;i>=0;--i){Real v=b[i];for(int j=i+1;j<3;++j)v-=a[i][j]*d[j];d[i]=v/a[i][i];require(std::isfinite(d[i]),"Nonfinite SR update");}return d;
}
Vec sr_parameters(const Config& c,const Vec& alpha,const Statistics& s,Real& scale) {
    Statistics transformed=s;
    Vec jac{1,1,alpha[2]}; // O_logdelta = delta * O_delta
    for(int p=0;p<3;++p) {
        transformed.force[p]=jac[p]*s.force[p];
        for(int q=0;q<3;++q)transformed.metric[p][q]=jac[p]*jac[q]*s.metric[p][q];
    }
    if(!c.optimize_delta) {
        // Solve the restricted two-parameter problem, not the full solve with d3 discarded.
        for(int p=0;p<3;++p)transformed.metric[p][2]=transformed.metric[2][p]=0;
        transformed.force[2]=0;
    }
    Vec direction=sr_direction(transformed,c.sr_step,c.sr_shift);
    scale=1;
    if(std::abs(direction[2])>c.max_log_delta_step)scale=c.max_log_delta_step/std::abs(direction[2]);
    Lattice l(c);
    for(int retry=0;retry<30;++retry,scale*=0.5L) {
        Vec next{alpha[0]+scale*direction[0],alpha[1]+scale*direction[1],alpha[2]};
        if(c.optimize_delta)next[2]=alpha[2]*std::exp(scale*direction[2]);
        try {Kernel trial(c,l,next);return next;}catch(const std::runtime_error&){}
    }
    throw std::runtime_error("Nonfinite/zero pair state in SR update after 30 backtracks");
}
namespace {
void save_parameters(const std::string& path,const Config& c,const Vec& a) {
    auto f=output(path);
    f<<"# Physical coordinates: chi1, Delta, delta. SR internally uses log(delta).\n"
     <<"Lx = "<<c.Lx<<"\nLy = "<<c.Ly<<"\nNb = "<<c.Nb<<"\nJB = "<<c.JB<<"\nVB = "<<c.VB<<"\nK = "<<c.K
     <<"\nchi1 = "<<a[0]<<"\nDelta = "<<a[1]<<"\ndelta = "<<a[2]
     <<"\noptimize_delta = "<<c.optimize_delta<<"\nmax_log_delta_step = "<<c.max_log_delta_step
     <<"\nwarmup_sweeps = "<<c.warmup_sweeps<<"\nsamples = "<<c.samples<<"\nsweeps_between = "<<c.sweeps_between
     <<"\nblock_size = "<<c.block_size<<"\nseed = "<<c.seed<<"\nlocal_probability = "<<c.local_probability
     <<"\nsr_iterations = 0\nsr_step = "<<c.sr_step<<"\nsr_shift = "<<c.sr_shift
     <<"\ncorrelations = "<<c.correlations<<"\noutput = "<<c.output<<"_continued\n";
}
}
void run(const Config& c) {
    c.validate();Lattice l(c);Vec alpha=c.alpha;Kernel initial(c,l,alpha);
    require(!vmc_fs::exists(c.output),"Output path already exists; choose a new output in .inp: "+c.output);
    vmc_fs::create_directories(c.output);
    save_parameters(c.output+"/checkpoint.inp",c,alpha);
    auto log=output(c.output+"/optimization.csv");
    log<<"iteration,chi1,Delta,delta,lambda_c,lambda_eff,E_total,blocked_error,acceptance,min_margin,min_gap,step_scale,delta_log_change,boundary_modes,F_chi1,F_Delta,F_delta,S00,S01,S02,S11,S12,S22\n";
    std::string stop_reason="requested_iterations_completed";int completed=0;
    for(int it=0;it<c.sr_iterations;++it) {
        Kernel k(c,l,alpha);auto s=sample(c,l,k,c.seed+std::uint64_t(it));Real scale=1;Vec next;
        try {next=sr_parameters(c,alpha,s,scale);}
        catch(const std::runtime_error& e) {
            stop_reason="update_failed_using_last_valid_state";
            std::cerr<<"SR stopped: "<<e.what()<<". Measuring the checkpoint state.\n";break;
        }
        log<<it<<','<<alpha[0]<<','<<alpha[1]<<','<<alpha[2]<<','<<k.lambda_c<<','<<k.lambda_eff
           <<','<<s.energy<<','<<s.error<<','<<s.acceptance<<','<<k.min_margin<<','<<k.min_gap<<','<<scale
           <<','<<std::log(next[2]/alpha[2])<<','<<k.boundary_modes;
        for(auto v:s.force)log<<','<<v;
        for(int p=0;p<3;++p)for(int q=p;q<3;++q)log<<','<<s.metric[p][q];
        log<<'\n';log.flush();alpha=next;++completed;
        save_parameters(c.output+"/checkpoint.inp",c,alpha);
        std::cout<<"SR "<<it<<" E/Ns="<<s.energy/l.n<<" +/- "<<s.error/l.n<<" delta_next="<<alpha[2]<<" scale="<<scale<<std::endl;
    }
    Kernel final_kernel(c,l,alpha);
    std::cout<<"Final measurement: "<<c.samples<<" samples, delta="<<alpha[2]<<std::endl;
    auto s=sample(c,l,final_kernel,c.seed+std::uint64_t(completed),c.output+"/samples.csv",c.correlations?c.output+"/correlations.csv":"");
    auto summary=output(c.output+"/summary.txt");
    summary<<"Lx = "<<c.Lx<<"\nLy = "<<c.Ly<<"\nNb = "<<c.Nb<<"\nnu = "<<Real(c.Nb)/l.n
        <<"\nJB = "<<c.JB<<"\nVB = "<<c.VB<<"\nK = "<<c.K<<"\nchi1 = "<<alpha[0]<<"\nDelta = "<<alpha[1]<<"\ndelta = "<<alpha[2]
        <<"\nlambda_c = "<<final_kernel.lambda_c<<"\nlambda_eff = "<<final_kernel.lambda_eff
        <<"\noptimize_delta = "<<c.optimize_delta<<"\nmax_log_delta_step = "<<c.max_log_delta_step
        <<"\nE_total = "<<s.energy<<"\nE_per_site = "<<s.energy/l.n<<"\nblocked_error_total = "<<s.error<<"\nblocked_error_per_site = "<<s.error/l.n
        <<"\nlocal_energy_variance = "<<s.variance<<"\nE_J = "<<s.components.hopping<<"\nE_V = "<<s.components.interaction<<"\nE_K = "<<s.components.correlated
        <<"\nacceptance = "<<s.acceptance<<"\nmin_A_minus_absB = "<<final_kernel.min_margin<<"\nmin_Ek = "<<final_kernel.min_gap
        <<"\nboundary_modes = "<<final_kernel.boundary_modes
        <<"\nseed = "<<c.seed<<"\nmeasurement_seed = "<<c.seed+std::uint64_t(completed)<<"\nsamples = "<<c.samples<<"\nblock_size = "<<c.block_size<<"\nsweeps_between = "<<c.sweeps_between
        <<"\nwarmup_sweeps = "<<c.warmup_sweeps<<"\nlocal_probability = "<<c.local_probability<<"\nsr_iterations = "<<c.sr_iterations
        <<"\nsr_iterations_completed = "<<completed<<"\nsr_stop_reason = "<<stop_reason
        <<"\nsr_step = "<<c.sr_step<<"\nsr_shift = "<<c.sr_shift<<"\ncorrelations = "<<c.correlations
        <<"\ninitial_chi1 = "<<c.alpha[0]<<"\ninitial_Delta = "<<c.alpha[1]<<"\ninitial_delta = "<<c.alpha[2]
        <<"\ninitial_lambda_eff = "<<initial.lambda_eff
        <<"\nError is a block-mean estimate; check block-size saturation and independent seeds.\n";
    save_parameters(c.output+"/final_parameters.inp",c,alpha);
    std::cout<<std::setprecision(12)<<"Final E/Ns = "<<s.energy/l.n<<" +/- "<<s.error/l.n<<"; acceptance = "<<s.acceptance<<"\nOutput: "<<c.output<<std::endl;
}
} // namespace vmc
