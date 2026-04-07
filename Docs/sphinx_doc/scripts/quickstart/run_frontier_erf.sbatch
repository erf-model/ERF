#!/usr/bin/env bash

#SBATCH -A <project>
#SBATCH -J erf
#SBATCH -o erf-%j.out
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest

# Frontier networking/runtime settings adapted from OLCF/WarpX examples.
export FI_MR_CACHE_MONITOR=memhooks
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software
export ROCFFT_RTC_CACHE_PATH=/dev/null
export OMP_NUM_THREADS=1

TOTAL_NMPI=$((SLURM_JOB_NUM_NODES * SLURM_NTASKS_PER_NODE))
srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${SLURM_NTASKS_PER_NODE} \
  ./erf_exec ../../Exec/CanonicalTests/ABL/inputs_most
