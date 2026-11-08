#pragma once
#include <array>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace vmc {
using Real = long double;
using Vec = std::array<Real,3>;
using Mat = std::array<Vec,3>;
struct Config {
    int Lx=8, Ly=8, Nb=8;
    Real JB=0.06L, VB=2.5L, K=1;
    Vec alpha{0.05L,0.05L,0.66L}; // chi1, Delta, delta = lambda_eff - lambda_c
    Real local_probability=0.5L;
    int warmup_sweeps=500, samples=2000, sweeps_between=2, block_size=100;
    int sr_iterations=0;
    Real sr_step=0.05L, sr_shift=1e-3L;
    bool optimize_delta=true;
    Real max_log_delta_step=1; // cap |change in log(delta)|, rescaling the whole SR step
    std::uint64_t seed=12345;
    bool correlations=true;
    bool diagnostics=true;
    int diagnostic_stride=1; // record every N final measurement samples
    std::string output="results";
    static Config read(const std::string& path);
    void validate() const;
};
struct Lattice {
    int lx,ly,n;
    std::vector<std::array<int,2>> bonds;
    std::vector<std::array<int,3>> corners; // center, x neighbor, y neighbor
    explicit Lattice(const Config& c);
    int site(int x,int y) const;
    int shift(int i,int dx,int dy) const;
};
struct Kernel {
    int n;
    std::vector<Real> f;
    std::array<std::vector<Real>,3> df;
    Real min_margin, min_gap;
    Real lambda_c, lambda_eff;
    Vec boundary_derivative{};
    int boundary_modes=0;
    Kernel(const Config&,const Lattice&,const Vec&);
    Real at(int i,int j) const { return f[i*n+j]; }
};
// Exact subset recurrence. Each instance caches sub-Hafnians of one configuration.
class Hafnian {
    const Kernel& kernel;
    std::vector<int> occupied;
    std::vector<Real> memo;
    Real solve(std::uint32_t mask);
public:
    Hafnian(const Kernel&,const std::vector<int>&);
    Real amplitude();
    Real moved(int particle,int destination);
    Vec log_derivative();
};
struct Energy { Real interaction=0, hopping=0, correlated=0;
    Real total() const {return interaction+hopping+correlated;}
};
Energy local_energy(const Config&,const Lattice&,const std::vector<int>&,Hafnian&);
struct Diagnostics {
    int na=0,nb=0;
    std::array<Real,3> density{}; // sites 0, +x, and r=(Lx/2,Ly/2)
    std::array<Real,4> one{},pair{}; // origin, translated +x, difference, spatial mean
    Real one_other=0; // spatial mean at (Lx/2,Ly/2+1)
};
Diagnostics measure_diagnostics(const Lattice&,const Kernel&,const std::vector<int>&,Hafnian&);
// Real-state Hermitian estimators. Pair indices run over the undirected NN bonds.
struct CorrelationAccumulator {
    std::vector<Real> one, pair, density;
    Real weight=0;
    explicit CorrelationAccumulator(const Lattice&);
    void add(const Lattice&,const Kernel&,const std::vector<int>&,Hafnian&,Real sample_weight=1);
    void write(const Lattice&,const std::string& legacy_csv) const;
};
class Sampler {
    const Config& cfg;
    const Lattice& lat;
    const Kernel& kernel;
    std::mt19937_64 rng;
    std::uniform_real_distribution<Real> uniform{0,1};
public:
    std::vector<int> occupied;
    std::uint64_t proposed=0,accepted=0;
    Sampler(const Config&,const Lattice&,const Kernel&,std::uint64_t seed);
    void sweep();
};
struct Statistics {
    Real energy=0,error=0,variance=0,acceptance=0;
    Energy components;
    Vec mean_o{},force{};
    Mat metric{};
};
// Pool centered statistics, including between-chain mean offsets.
Statistics merge_statistics(const std::vector<Statistics>&,const std::vector<int>&);
Statistics sample(const Config&,const Lattice&,const Kernel&,std::uint64_t,
                  const std::string& trace_path="", const std::string& correlation_path="");
Vec sr_direction(const Statistics&,Real step,Real shift);
Vec sr_parameters(const Config&,const Vec&,const Statistics&,Real& scale);
void run(const Config&);
} // namespace vmc
