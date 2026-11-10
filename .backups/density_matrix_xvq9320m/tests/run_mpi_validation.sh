#!/bin/bash
#SBATCH --job-name=vmc_mpi_check
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=norm2
#SBATCH --time=00:20:00
#SBATCH --mem=16G
set -eo pipefail
source /opt/intel/oneapi/mpi/latest/env/vars.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /DATA/cuic/Hardcore_boson_VMC
make all mpi test LDFLAGS=-lstdc++fs
python3 tests/test_mpi.py "$PWD" "$PWD/validation_mpi_20261101/check_${SLURM_JOB_ID}"
python3 tests/benchmark_mpi.py "$PWD" "$PWD/validation_mpi_20261101/bench_${SLURM_JOB_ID}"
