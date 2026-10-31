#!/bin/bash

#SBATCH --account=<account_or_project_g>
#SBATCH --constraint=gpu&hbm40g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --gpu-bind=none
#SBATCH --time=00:30:00
#SBATCH --job-name=erf
#SBATCH --output=erf_%j.out

# Run with the modules the executable was built with (submitted from install/bin).
# GNU Make builds resolve libmpi_gtl_cuda and libcudart at run time from the loaded
# modules, so a different stack here fails with "libcudart.so.NN: cannot open shared object file".
source ../../Build/machines/perlmutter_erf.profile

# GPU-aware MPI optimizations
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_SUPPORT_ENABLED=1
export SLURM_CPU_BIND=cores

srun -n $((SLURM_NNODES * SLURM_NTASKS_PER_NODE)) --cpus-per-task=${SLURM_CPUS_PER_TASK:-4} --cpu-bind=cores bash -c '
  export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))
  ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most amrex.use_gpu_aware_mpi=1
'
