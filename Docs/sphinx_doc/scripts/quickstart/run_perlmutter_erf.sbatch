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

# GPU-aware MPI optimizations
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_SUPPORT_ENABLED=1
export SLURM_CPU_BIND=cores

srun -n $((SLURM_NNODES * SLURM_NTASKS_PER_NODE)) --cpus-per-task=${SLURM_CPUS_PER_TASK:-4} --cpu-bind=cores bash -c '
  export CUDA_VISIBLE_DEVICES=$((3-SLURM_LOCALID))
  ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most amrex.use_gpu_aware_mpi=1
'
