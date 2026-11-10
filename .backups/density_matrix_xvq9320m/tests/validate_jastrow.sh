#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=norm2
#SBATCH --time=00:15:00
set -e
source /opt/intel/oneapi/mpi/latest/env/vars.sh
cd /DATA/cuic/Hardcore_boson_VMC
make all mpi test LDFLAGS=-lstdc++fs
python3 tests/test_mpi.py "$PWD" "$PWD/validation_mpi_20261101/old_${SLURM_JOB_ID}"
python3 tests/test_jastrow_mpi.py "$PWD" "$PWD/validation_mpi_20261101/jastrow_${SLURM_JOB_ID}"
