#!/bin/bash
#SBATCH --job-name=pair-cache-test
#SBATCH --partition=norm2
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/DATA/cuic/Hardcore_boson_VMC/analysis/pair_cache_%j.log
set -euo pipefail
cd /DATA/cuic/Hardcore_boson_VMC
export OMP_NUM_THREADS=1
export I_MPI_PIN=1
./vmc_tests
python3 tests/test_mpi.py "$PWD" "$PWD/analysis/pair_cache_mpi_${SLURM_JOB_ID}"
python3 tests/test_jastrow_mpi.py "$PWD" "$PWD/analysis/pair_cache_jastrow_${SLURM_JOB_ID}"
python3 tests/benchmark_pair_cache.py "$PWD" "$PWD/analysis/pair_cache_benchmark_${SLURM_JOB_ID}" "$PWD/.backups/pair_cache_20261105/hardcore_vmc_mpi"
