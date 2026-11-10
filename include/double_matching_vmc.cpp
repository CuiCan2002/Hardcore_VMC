
#include "vmc.hpp"
#include <algorithm>
#include <cmath>
#include <experimental/filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#ifdef DOUBLE_MATCHING_MPI
#include <mpi.h>
#endif

namespace fs = std::experimental::filesystem;
using namespace vmc;

namespace {
int rank_id=0;
int ranks=1;

int chain_samples(int total, int block_size) {
    const int blocks=total/block_size;
    return (blocks/ranks+(rank_id<blocks%ranks))*block_size;
}

void reduce_sum(Real* values, int count) {
#ifdef DOUBLE_MATCHING_MPI
    if (ranks==1) return;
    if (rank_id==0) MPI_Reduce(MPI_IN_PLACE,values,count,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    else MPI_Reduce(values,nullptr,count,MPI_LONG_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
#else
    (void)values;
    (void)count;
#endif
}

void reduce_correlations(CorrelationAccumulator& corr) {
    reduce_sum(corr.one.data(),int(corr.one.size()));
    reduce_sum(corr.density.data(),int(corr.density.size()));
    reduce_sum(corr.density_matrix.data(),int(corr.density_matrix.size()));
    reduce_sum(&corr.weight,1);
}

void check(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

Config read_config(const std::string& path) {
    Config c;
    c.output = "results_double_matching";
    c.diagnostics = false;
    std::ifstream input(path);
    check(bool(input), "Cannot open " + path);
    std::string line;
    std::vector<std::string> keys;
    while (std::getline(input, line)) {
        line.resize(std::min(line.size(), line.find('#')));
        std::replace(line.begin(), line.end(), '=', ' ');
        std::istringstream in(line);
        std::string key, extra;
        if (!(in >> key)) continue;
        check(std::find(keys.begin(), keys.end(), key) == keys.end(), "Duplicate key: " + key);
        keys.push_back(key);
        auto read = [&](const char* name, auto& value) {
            if (key != name) return false;
            check(bool(in >> value) && !(in >> extra), "Invalid value for " + key);
            return true;
        };
        if (read("Lx", c.Lx) || read("Ly", c.Ly) || read("Nb", c.Nb) ||
            read("JB", c.JB) || read("VB", c.VB) || read("K", c.K) ||
            read("chi1", c.alpha[0]) || read("Delta", c.alpha[1]) || read("delta", c.alpha[2]) ||
            read("v1", c.alpha[3]) || read("v2", c.alpha[4]) || read("v3", c.alpha[5]) ||
            read("Ifjastrow", c.Ifjastrow) || read("optimize_delta", c.optimize_delta) ||
            read("local_probability", c.local_probability) || read("warmup_sweeps", c.warmup_sweeps) ||
            read("samples", c.samples) || read("samples_sr", c.samples_sr) || read("samples_sr2", c.samples_sr2) ||
            read("sr_iterations", c.sr_iterations) || read("sr_iterations2", c.sr_iterations2) ||
            read("sr_step", c.sr_step) || read("sr_shift", c.sr_shift) ||
            read("max_log_delta_step", c.max_log_delta_step) || read("sweeps_between", c.sweeps_between) ||
            read("block_size", c.block_size) || read("seed", c.seed) || read("correlations", c.correlations) ||
            read("pair_measure_stride", c.pair_measure_stride) || read("diagnostics", c.diagnostics) ||
            read("diagnostic_stride", c.diagnostic_stride) || read("output", c.output)) continue;
        throw std::runtime_error("Unknown input key: " + key);
    }
    check(c.Lx >= 3 && c.Ly >= 3 && c.Lx <= 128 && c.Ly <= 128 && c.Lx*c.Ly <= 4096,
          "Invalid lattice size");
    check(c.Nb >= 2 && c.Nb % 2 == 0 && c.Nb < c.Lx*c.Ly, "Nb must be even and less than Ns");
    // Reuse input checks that do not depend on the exact Hafnian limit.
    Config validation = c;
    validation.Nb = std::min(c.Nb, 20);
    validation.validate();
    return c;
}

bool positive_kernel(const Kernel& k) {
    for (int i = 0; i < k.n; ++i)
        for (int j = i+1; j < k.n; ++j)
            if (!(k.at(i,j) > 0) || !std::isfinite(k.at(i,j))) return false;
    return true;
}

std::ofstream output(const fs::path& path) {
    std::ofstream out(path.string());
    check(bool(out), "Cannot write " + path.string());
    out << std::setprecision(std::numeric_limits<Real>::max_digits10);
    return out;
}

Real potential(const Kernel& k, int a, int b) {
    const int shell=k.shell[a*k.n+b];
    return k.jastrow && shell>=0 ? k.v[3+shell] : 0;
}

Real jastrow_change(const Kernel& k, const std::vector<int>& sites,
                    int p, int a, int q=-1, int b=-1) {
    Real d=0;
    for (int i=0;i<int(sites.size());++i) {
        if (i==p || i==q) continue;
        d+=potential(k,a,sites[i])-potential(k,sites[p],sites[i]);
        if (q>=0) d+=potential(k,b,sites[i])-potential(k,sites[q],sites[i]);
    }
    if (q>=0) d+=potential(k,a,b)-potential(k,sites[p],sites[q]);
    return d;
}

Real edge_ratio(const Kernel& k, const std::vector<int>& sites, const std::vector<int>& m,
                int p, int a, int q=-1, int b=-1) {
    if (m[p]==q) return k.at(a,b)/k.at(sites[p],sites[q]);
    Real ratio=k.at(a,sites[m[p]])/k.at(sites[p],sites[m[p]]);
    if (q>=0) ratio*=k.at(b,sites[m[q]])/k.at(sites[q],sites[m[q]]);
    return ratio;
}

class MatchingWorm {
    std::vector<std::vector<Real>> cumulative;
    Real total_weight=0;
    bool dirty=true;
public:
    explicit MatchingWorm(int n):cumulative(n,std::vector<Real>(n)) {}
    void invalidate() { dirty=true; }

    // Return 1 for a changed matching, 0 for an unchanged return, -1 for a capped excursion.
    int update(const Kernel& k, const std::vector<int>& sites, std::vector<int>& m,
               int head, std::mt19937_64& rng, std::uniform_real_distribution<Real>& unit) {
        const int n=int(sites.size());
        if (dirty) {
            total_weight=0;
            for (int i=0;i<n;++i) {
                Real sum=0;
                for (int j=0;j<n;++j) cumulative[i][j]=sum+=k.at(sites[i],sites[j]);
                total_weight=std::max(total_weight,sum);
            }
            dirty=false;
        }
        auto next=m;
        int tail=next[head];
        // A common normalization leaves unused mass as self-loops.
        // Rejecting paths above the length cap also rejects their reverses.
        for (int step=0;step<100000;++step) {
            const auto& row=cumulative[tail];
            const int target=int(std::upper_bound(row.begin(),row.end(),unit(rng)*total_weight)-row.begin());
            if (target==n || target==tail) continue;
            if (target==head) {
                next[head]=tail;
                next[tail]=head;
                const bool changed=next!=m;
                m.swap(next);
                return changed;
            }
            const int displaced=next[target];
            next[tail]=target;
            next[target]=tail;
            tail=displaced;
        }
        return -1;
    }
};

class Chain {
    const Config& cfg;
    const Lattice& lat;
    const Kernel& kernel;
    std::mt19937_64 rng;
    std::uniform_real_distribution<Real> unit{0,1};
    std::array<std::vector<int>,2> match;
    MatchingWorm matching_worm;

    int draw(int n) { return std::uniform_int_distribution<int>(0,n-1)(rng); }

    void move() {
        ++position_trials;
        const int p = draw(cfg.Nb), old = sites[p];
        const int axis = draw(4);
        const int dest = unit(rng) < cfg.local_probability
            ? lat.shift(old,(axis==0)-(axis==1),(axis==2)-(axis==3)) : draw(lat.n);
        if (index[dest] >= 0) return;
        const Real log_ratio = std::log(edge_ratio(kernel,sites,match[0],p,dest))+std::log(edge_ratio(kernel,sites,match[1],p,dest))
                               -2*jastrow_change(kernel,sites,p,dest);
        if (std::log(unit(rng)) >= log_ratio) return;
        index[old] = -1;
        index[dest] = p;
        sites[p] = dest;
        matching_worm.invalidate();
        ++position_accepts;
    }

    void worm() {
        ++worm_trials;
        auto& m=match[draw(2)];
        const int head=draw(cfg.Nb);
        const int result=matching_worm.update(kernel,sites,m,head,rng,unit);
        worm_changes+=result>0;
        worm_capped+=result<0;
    }

public:
    std::vector<int> sites, index;
    std::uint64_t position_trials=0, position_accepts=0, worm_trials=0, worm_changes=0, worm_capped=0;

    Chain(const Config& c, const Lattice& l, const Kernel& k, std::uint64_t seed)
        : cfg(c), lat(l), kernel(k), rng(seed), matching_worm(c.Nb),
          sites(l.n), index(l.n,-1) {
        std::iota(sites.begin(),sites.end(),0);
        std::shuffle(sites.begin(),sites.end(),rng);
        sites.resize(c.Nb);
        for (int i = 0; i < c.Nb; ++i) index[sites[i]] = i;
        for (auto& m : match) {
            m.resize(c.Nb);
            for (int i = 0; i < c.Nb; ++i) m[i] = i^1;
        }
    }

    void sweep() {
        for (int i = 0; i < cfg.Nb; ++i) {
            if (draw(2)) move();
            else worm();
        }
    }

    Real ratio(int replica, int p, int a) const {
        return edge_ratio(kernel,sites,match[replica],p,a)*std::exp(-jastrow_change(kernel,sites,p,a));
    }

    std::array<Energy,2> energies() const {
        std::array<Energy,2> e{};
        auto hop = [&](int a, int b, Real coupling, bool correlated) {
            if (coupling == 0 || (index[a]>=0)==(index[b]>=0)) return;
            if (index[b]>=0) std::swap(a,b);
            for (int r = 0; r < 2; ++r) {
                auto& component = correlated ? e[r].correlated : e[r].hopping;
                component -= coupling*ratio(r,index[a],b);
            }
        };
        for (auto bond : lat.bonds) {
            const int a = bond[0], b = bond[1];
            if (index[a]>=0 && index[b]>=0)
                for (auto& value : e) value.interaction += cfg.VB;
            hop(a,b,cfg.JB,false);
        }
        for (auto t : lat.corners) if (index[t[0]]>=0) hop(t[1],t[2],cfg.K,true);
        return e;
    }

    Vec derivative(int replica) const {
        Vec o{};
        for (int i = 0; i < cfg.Nb; ++i) {
            int j = match[replica][i];
            if (i > j) continue;
            const int edge = sites[i]*lat.n+sites[j];
            for (int p = 0; p < 3; ++p) o[p] += kernel.df[p][edge]/kernel.f[edge];
        }
        if (kernel.jastrow) {
            auto count = kernel.counts(sites);
            for (int p = 0; p < 3; ++p) o[p+3] = -count[p];
        }
        return o;
    }

    void measure(CorrelationAccumulator& corr, Real& pair_trace) const {
        ++corr.weight;
        for (int i : sites) {
            ++corr.one[i*lat.n+i];
            for (int j = 0; j < lat.n; ++j) if (index[j]<0) {
                Real w = (ratio(0,index[i],j)+ratio(1,index[i],j))/4;
                corr.one[i*lat.n+j] += w;
                corr.one[j*lat.n+i] += w;
            }
        }
        corr.add_density(lat,sites);
        for (const auto& bond:lat.bonds)
            pair_trace+=index[bond[0]]>=0 && index[bond[1]]>=0;
    }
};

class PairChain {
    const Lattice& lat;
    const Kernel& kernel;
    int background_size;
    static constexpr Real close_proposal=0.2L;
    std::mt19937_64 rng;
    std::uniform_real_distribution<Real> unit{0,1};
    std::array<std::vector<int>,2> sites, matching;
    std::array<MatchingWorm,2> worms;
    std::vector<bool> in_background;

    int draw(int n) { return std::uniform_int_distribution<int>(0,n-1)(rng); }

    Real log_ratio(int side, int p, int a, int q=-1, int b=-1) const {
        return std::log(edge_ratio(kernel,sites[side],matching[side],p,a,q,b))
               -jastrow_change(kernel,sites[side],p,a,q,b);
    }

    void move_endpoint() {
        ++trials[0];
        const int side=draw(2), other=end[1-side], old=end[side];
        const int candidate=unit(rng)<close_proposal ? other : draw(int(lat.bonds.size()));
        if (old==candidate) return;
        const auto bond=lat.bonds[candidate];
        if (in_background[bond[0]] || in_background[bond[1]]) return;
        const Real uniform=(1-close_proposal)/lat.bonds.size();
        const Real forward=uniform+close_proposal*(candidate==other);
        const Real reverse=uniform+close_proposal*(old==other);
        const Real bias=std::log(diagonal_bias)*((candidate==other)-(old==other));
        const Real log_weight=log_ratio(side,background_size,bond[0],background_size+1,bond[1])
                              +std::log(reverse/forward)+bias;
        if (std::log(unit(rng))>=log_weight) return;
        end[side]=candidate;
        sites[side][background_size]=bond[0];
        sites[side][background_size+1]=bond[1];
        worms[side].invalidate();
        ++changes[0];
    }

    void move_background() {
        if (background_size==0) return;
        ++trials[1];
        const int p=draw(background_size), dest=draw(lat.n);
        if (in_background[dest]) return;
        for (int bond:end)
            if (dest==lat.bonds[bond][0] || dest==lat.bonds[bond][1]) return;
        if (std::log(unit(rng))>=log_ratio(0,p,dest)+log_ratio(1,p,dest)) return;
        in_background[sites[0][p]]=false;
        in_background[dest]=true;
        for (int side=0;side<2;++side) {
            sites[side][p]=dest;
            worms[side].invalidate();
        }
        ++changes[1];
    }

public:
    static constexpr Real diagonal_bias=20;
    std::array<int,2> end{};
    std::array<Real,3> trials{}, changes{};
    Real capped=0;

    PairChain(const Lattice& l, const Kernel& k, int nb, std::uint64_t seed)
        :lat(l),kernel(k),background_size(nb-2),rng(seed),worms{MatchingWorm(nb),MatchingWorm(nb)},
         in_background(l.n,false) {
        end.fill(draw(int(l.bonds.size())));
        const auto bond=l.bonds[end[0]];
        auto& first=sites[0];
        for (int site=0;site<l.n;++site)
            if (site!=bond[0] && site!=bond[1]) first.push_back(site);
        std::shuffle(first.begin(),first.end(),rng);
        first.resize(background_size);
        for (int site:first) in_background[site]=true;
        first.insert(first.end(),bond.begin(),bond.end());
        sites[1]=first;
        for (auto& m:matching) {
            m.resize(nb);
            for (int i=0;i<nb;++i) m[i]=i^1;
        }
    }

    void sweep() {
        for (int i=0;i<background_size+2;++i) {
            switch (draw(3)) {
            case 0: move_endpoint(); break;
            case 1: move_background(); break;
            default:
                ++trials[2];
                const int side=draw(2), head=draw(background_size+2);
                const int result=worms[side].update(kernel,sites[side],matching[side],head,rng,unit);
                changes[2]+=result>0;
                capped+=result<0;
            }
        }
    }
};

int pair_bin(const Lattice& l, int a, int b) {
    // Lattice::bonds stores +x and +y in that order for each site.
    const int i=a/2, j=b/2;
    return (2*(a%2)+b%2)*l.n+l.site(j%l.lx-i%l.lx,j/l.lx-i/l.lx);
}

void sample_pairs(const Config& c, const Lattice& l, const Vec& a, std::uint64_t seed,
                  Real trace_sum, CorrelationAccumulator& corr, const fs::path& dir) {
    const Kernel k(c,l,a);
    PairChain chain(l,k,c.Nb,seed+std::uint64_t(rank_id)*UINT64_C(0x9e3779b97f4a7c15));
    for (int i=0;i<c.warmup_sweeps;++i) chain.sweep();
    chain.trials.fill(0);
    chain.changes.fill(0);
    chain.capped=0;
    const int n=chain_samples(c.samples,c.block_size)/c.pair_measure_stride;
    constexpr int block_count=20;
    std::array<Real,2*block_count> blocks{};
    std::vector<Real> hist(4*l.n);
    for (int i=0;i<n;++i) {
        for (int j=0;j<c.sweeps_between;++j) chain.sweep();
        ++hist[pair_bin(l,chain.end[0],chain.end[1])];
        const int block=int(std::int64_t(i)*block_count/n);
        ++blocks[2*block];
        blocks[2*block+1]+=chain.end[0]==chain.end[1];
    }
    std::vector<Real> all_blocks;
#ifdef DOUBLE_MATCHING_MPI
    if (ranks>1) {
        if (rank_id==0) all_blocks.resize(ranks*blocks.size());
        MPI_Gather(blocks.data(),int(blocks.size()),MPI_LONG_DOUBLE,all_blocks.data(),
                   int(blocks.size()),MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
    } else
#endif
        all_blocks.assign(blocks.begin(),blocks.end());
    reduce_sum(hist.data(),int(hist.size()));
    reduce_sum(chain.trials.data(),int(chain.trials.size()));
    reduce_sum(chain.changes.data(),int(chain.changes.size()));
    reduce_sum(&chain.capped,1);
    reduce_sum(&trace_sum,1);
    if (rank_id!=0) return;

    auto block_log=output(dir/"pair_sector_blocks.csv");
    block_log<<"rank,block,samples,biased_diagonal_visits\n";
    for (std::size_t i=0;i<all_blocks.size();i+=2)
        if (all_blocks[i]>0)
            block_log<<i/blocks.size()<<','<<(i%blocks.size())/2<<','<<all_blocks[i]<<','<<all_blocks[i+1]<<'\n';

    const Real records=std::accumulate(hist.begin(),hist.end(),Real(0));
    const Real diagonal=hist[0]+hist[3*l.n];
    check(diagonal>0,"No diagonal pair-sector visits; increase pair sampling");
    check(trace_sum>0,"No occupied nearest-neighbor bonds in physical samples");
    hist[0]/=PairChain::diagonal_bias;
    hist[3*l.n]/=PairChain::diagonal_bias;
    const Real corrected_diagonal=diagonal/PairChain::diagonal_bias;
    const Real trace=trace_sum/corr.weight;
    const Real mode_fraction=std::accumulate(hist.begin(),hist.end(),Real(0))
                             /(l.bonds.size()*corrected_diagonal);

    // Each translation bin sums Ns equivalent entries. Average opposite ordered bins
    // to enforce Hermiticity and retain the existing matrix-file normalization.
    corr.pair_weight=records;
    const auto nb=l.bonds.size();
    const Real factor=records*trace/(2*l.n*corrected_diagonal);
    for (std::size_t i=0;i<nb;++i) for (std::size_t j=0;j<nb;++j)
        corr.pair[i*nb+j]=factor*(hist[pair_bin(l,int(i),int(j))]+hist[pair_bin(l,int(j),int(i))]);
    auto info=output(dir/"pair_sector_summary.txt");
    info<<"pair_estimator = direct_pair_sector\ntranslation_averaged = 1\n"
        <<"pair_sector_samples = "<<records<<"\ndiagonal_bias = "<<PairChain::diagonal_bias
        <<"\nbiased_diagonal_visits = "<<diagonal<<"\npair_trace = "<<trace
        <<"\nuniform_pair_occupation = "<<trace*mode_fraction
        <<"\nuniform_pair_fraction = "<<mode_fraction<<"\n";
    const char* names[]={"endpoint_acceptance","background_acceptance","worm_change_rate"};
    for (int i=0;i<3;++i) info<<names[i]<<" = "<<chain.changes[i]/std::max(Real(1),chain.trials[i])<<'\n';
    info<<"capped_worms = "<<chain.capped<<'\n';
}

Statistics sample_matching(const Config& c, const Lattice& l, const Vec& a, int total, std::uint64_t seed,
                           bool optimize, CorrelationAccumulator* corr, std::ostream* trace, Real* pair_trace=nullptr) {
    Kernel k(c,l,a);
    check(positive_kernel(k),"Double matching requires strictly positive off-diagonal f");
    const int n=chain_samples(total,c.block_size);
    Chain chain(c,l,k,seed+std::uint64_t(rank_id)*UINT64_C(0x9e3779b97f4a7c15));
    for (int i = 0; i < c.warmup_sweeps; ++i) chain.sweep();
    chain.position_trials = chain.position_accepts = chain.worm_trials = chain.worm_changes = chain.worm_capped = 0;
    Statistics s;
    Real block=0, blocks_squared=0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < c.sweeps_between; ++j) chain.sweep();
        const auto e = chain.energies();
        const Real e0=e[0].total(), e1=e[1].total(), value=(e0+e1)/2;
        check(std::isfinite(value),"Nonfinite energy estimator");
        s.energy += value;
        s.components.interaction += (e[0].interaction+e[1].interaction)/2;
        s.components.hopping += (e[0].hopping+e[1].hopping)/2;
        s.components.correlated += (e[0].correlated+e[1].correlated)/2;
        block += value;
        if ((i+1)%c.block_size==0) {
            const Real mean=block/c.block_size;
            blocks_squared+=mean*mean;
            block=0;
        }
        if (optimize) {
            const Vec o0=chain.derivative(0), o1=chain.derivative(1);
            // The replicas are conditionally independent at fixed C.
            for (int p=0;p<6;++p) {
                s.mean_o[p]+=(o0[p]+o1[p])/2;
                s.force[p]+=(o0[p]*e1+o1[p]*e0)/2;
                for (int q=0;q<6;++q) s.metric[p][q]+=(o0[p]*o1[q]+o1[p]*o0[q])/2;
            }
        }
        if (corr) chain.measure(*corr,*pair_trace);
        if (trace && (i+1)%c.diagnostic_stride==0) {
            int imbalance=0;
            for (int site : chain.sites) imbalance += (site%l.lx+site/l.lx)%2 ? -1 : 1;
            *trace << i+1 << ',' << value/l.n << ',' << imbalance << '\n';
        }
    }
    Real position_trials=chain.position_trials, position_accepts=chain.position_accepts;
    Real worm_trials=chain.worm_trials, worm_changes=chain.worm_changes, worm_capped=chain.worm_capped;
#ifdef DOUBLE_MATCHING_MPI
    if (ranks>1) {
        // Reduce uncentered moments before centering across all chains.
        auto fields={&s.energy,&s.components.interaction,&s.components.hopping,&s.components.correlated,
                     &blocks_squared,&position_trials,&position_accepts,&worm_trials,&worm_changes,&worm_capped};
        std::vector<Real> moments;
        for (auto ptr:fields) moments.push_back(*ptr);
        if (optimize) {
            moments.insert(moments.end(),s.mean_o.begin(),s.mean_o.end());
            moments.insert(moments.end(),s.force.begin(),s.force.end());
            for (const auto& row:s.metric) moments.insert(moments.end(),row.begin(),row.end());
        }
        MPI_Allreduce(MPI_IN_PLACE,moments.data(),int(moments.size()),MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        auto value=moments.begin();
        for (auto ptr:fields) *ptr=*value++;
        if (optimize) {
            for (auto& x:s.mean_o) x=*value++;
            for (auto& x:s.force) x=*value++;
            for (auto& row:s.metric) for (auto& x:row) x=*value++;
        }
    }
#endif
    s.energy/=total;
    s.components.interaction/=total;
    s.components.hopping/=total;
    s.components.correlated/=total;
    const int blocks=total/c.block_size;
    s.error=std::sqrt(std::max(Real(0),blocks_squared/blocks-s.energy*s.energy)/(blocks-1));
    for (auto& o:s.mean_o) o/=total;
    for (int p=0;p<6;++p) {
        s.force[p]=s.force[p]/total-s.mean_o[p]*s.energy;
        for (int q=0;q<6;++q) s.metric[p][q]=s.metric[p][q]/total-s.mean_o[p]*s.mean_o[q];
    }
    s.acceptance=position_accepts/std::max(Real(1),position_trials);
    if (rank_id==0)
        std::cout << "position_acceptance=" << s.acceptance << " worm_change_rate="
                  << worm_changes/std::max(Real(1),worm_trials)
                  << " capped_worms=" << worm_capped << std::endl;
    return s;
}

bool positive_sr_metric(const Config& c, const Vec& a, const Statistics& s) {
    Mat lower{};
    Vec jac{1,1,a[2],1,1,1};
    auto active = [&](int i) { return (i != 2 || c.optimize_delta) && (i < 3 || c.Ifjastrow); };
    for (int i=0;i<6;++i) for (int j=0;j<=i;++j) {
        Real v = active(i) && active(j) ? jac[i]*jac[j]*s.metric[i][j] : 0;
        if (i==j) v+=c.sr_shift;
        for (int k=0;k<j;++k) v-=lower[i][k]*lower[j][k];
        if (i==j) {
            if (!(v>0) || !std::isfinite(v)) return false;
            lower[i][i]=std::sqrt(v);
        } else lower[i][j]=v/lower[j][j];
    }
    return true;
}

void save_state(const fs::path& path, const Config& c, const Vec& a) {
    auto out=output(path);
    out<<"Lx = "<<c.Lx<<"\nLy = "<<c.Ly<<"\nNb = "<<c.Nb<<"\nJB = "<<c.JB
       <<"\nVB = "<<c.VB<<"\nK = "<<c.K<<"\nIfjastrow = "<<c.Ifjastrow
       <<"\nchi1 = "<<a[0]<<"\nDelta = "<<a[1]<<"\ndelta = "<<a[2]
       <<"\nv1 = "<<a[3]<<"\nv2 = "<<a[4]<<"\nv3 = "<<a[5]
       <<"\nsr_iterations = 0\nsr_iterations2 = 0\nsamples = "<<c.samples
       <<"\nwarmup_sweeps = "<<c.warmup_sweeps<<"\nsweeps_between = "<<c.sweeps_between
       <<"\nblock_size = "<<c.block_size<<"\nseed = "<<c.seed
       <<"\npair_measure_stride = "<<c.pair_measure_stride<<"\ncorrelations = "<<c.correlations<<'\n';
}
}

int main(int argc, char** argv) {
#ifdef DOUBLE_MATCHING_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
#endif
    try {
        check(argc==2,"Usage: hardcore_vmc_double_matching input.inp");
        const Config c=read_config(argv[1]);
        auto validate_budget=[&](int n, const char* name) {
            check(n/c.block_size/ranks>=2,std::string(name)+" needs at least two complete blocks per chain");
        };
        validate_budget(c.samples,"samples");
        if(c.sr_iterations>0) validate_budget(c.sr_samples(),"samples_sr");
        if(c.sr_iterations2>0) validate_budget(c.sr2_samples(),"samples_sr2");
        if(c.correlations)
            check(c.samples/c.block_size/ranks*c.block_size>=c.pair_measure_stride,
                  "Each chain needs at least one complete pair measurement interval");
        const Lattice l(c);
        Vec a=c.alpha;
        if (!c.Ifjastrow) for(int i=3;i<6;++i) a[i]=0;
        check(positive_kernel(Kernel(c,l,a)),"Initial kernel has nonpositive off-diagonal weights");
        const fs::path dir(c.output);
        std::ofstream log;
        if(rank_id==0) {
            check(!fs::exists(dir),"Output already exists: "+c.output);
            fs::create_directories(dir);
            if(c.diagnostics && ranks>1) fs::create_directory(dir/"chains");
            save_state(dir/"checkpoint.inp",c,a);
            log=output(dir/"optimization.csv");
            log<<"iteration,samples,chi1,Delta,delta,v1,v2,v3,E_per_site,error,F_chi1,F_Delta,F_log_delta,F_v1,F_v2,F_v3,sr_shift,step_scale\n";
            std::cout<<"Double matching: "<<ranks<<" chain(s), sample budgets are totals"<<std::endl;
        }
#ifdef DOUBLE_MATCHING_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        const auto iterations=std::int64_t(c.sr_iterations)+c.sr_iterations2;
        for (std::int64_t it=0; it<iterations; ++it) {
            const int n=it<c.sr_iterations ? c.sr_samples() : c.sr2_samples();
            const auto s=sample_matching(c,l,a,n,c.seed+it,true,nullptr,nullptr);
            if(rank_id==0) {
                Config regularized=c;
                for(int retry=0;retry<30 && !positive_sr_metric(regularized,a,s);++retry)
                    regularized.sr_shift*=2;
                check(positive_sr_metric(regularized,a,s),"SR metric could not be regularized");
                Vec next{};
                Real scale=0;
                bool valid=false;
                for (int retry=0;retry<30;++retry,regularized.sr_step/=2) {
                    next=sr_parameters(regularized,a,s,scale);
                    if (positive_kernel(Kernel(c,l,next))) { valid=true; break; }
                }
                check(valid,"SR update left the positive-kernel domain");
                log<<it<<','<<n;
                for (Real v:a) log<<','<<v;
                log<<','<<s.energy/l.n<<','<<s.error/l.n;
                for(int i=0;i<6;++i) log<<','<<s.force[i]*(i==2?a[2]:1);
                log<<','<<regularized.sr_shift<<','<<scale*regularized.sr_step/c.sr_step<<std::endl;
                a=next;
                save_state(dir/"checkpoint.inp",c,a);
            }
#ifdef DOUBLE_MATCHING_MPI
            MPI_Bcast(a.data(),int(a.size()),MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);
#endif
        }
        std::unique_ptr<CorrelationAccumulator> corr;
        if(c.correlations) corr=std::make_unique<CorrelationAccumulator>(l);
        std::ofstream trace;
        if(c.diagnostics) {
            const auto path=ranks==1 ? dir/"samples.csv" :
                dir/"chains"/("samples_rank"+std::to_string(rank_id)+".csv");
            trace=output(path);
            trace<<"sample,E_per_site,NA_minus_NB\n";
        }
        Real pair_trace=0;
        const auto s=sample_matching(c,l,a,c.samples,c.seed+iterations,false,corr.get(),c.diagnostics?&trace:nullptr,&pair_trace);
        if(corr) {
            reduce_correlations(*corr);
            if(rank_id==0) std::cout<<"Measuring direct pair sector..."<<std::endl;
            sample_pairs(c,l,a,c.seed+iterations+UINT64_C(0xd1b54a32d192ed03),pair_trace,*corr,dir);
        }
        if(rank_id==0) {
            if(corr) corr->write(l,(dir/"correlations.csv").string());
            save_state(dir/"final_parameters.inp",c,a);
            auto summary=output(dir/"summary.txt");
            summary<<"E_per_site = "<<s.energy/l.n<<"\nblocked_error_per_site = "<<s.error/l.n
                   <<"\nE_J = "<<s.components.hopping<<"\nE_V = "<<s.components.interaction
                   <<"\nE_K = "<<s.components.correlated<<"\nposition_acceptance = "<<s.acceptance
                   <<"\nsamples = "<<c.samples<<"\nsamples_scope = total_across_chains"
                   <<"\nmpi_ranks = "<<ranks<<"\nsr_iterations_completed = "<<iterations<<'\n';
            std::cout<<"Final E/Ns = "<<s.energy/l.n<<" +/- "<<s.error/l.n<<'\n';
        }
    } catch(const std::exception& e) {
        std::cerr<<"double matching (rank "<<rank_id<<"): "<<e.what()<<std::endl;
#ifdef DOUBLE_MATCHING_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
#endif
        return 1;
    }
#ifdef DOUBLE_MATCHING_MPI
    MPI_Finalize();
#endif
}


