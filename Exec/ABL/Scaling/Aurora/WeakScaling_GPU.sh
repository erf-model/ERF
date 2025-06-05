#!/bin/bash

export MPICH_GPU_SUPPORT_ENABLED=1

lc=1; nx=512; ny=512; probx=2048; proby=2048;
# the -n argument is (--ntasks-per-node) * (-N) = (number of MPI ranks per node) * (number of nodes)

NRANKS_PER_NODE=12
NDEPTH=4
NTHREADS=1

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=${NTHREADS}
export OMP_PLACES=cores

for n_nodes in 1 2 4 8 16
do
  NTOTRANKS=$(( n_nodes * NRANKS_PER_NODE ))
  echo "Test $lc: $n_nodes nodes, $NRANKS_PER_NODE tasks per node, $NTOTRANKS procs, $probx $proby 1024, $nx $ny 256"

  mpiexec -n $NTOTRANKS --ppn $NRANKS_PER_NODE --depth=${NDEPTH} --cpu-bind depth /soft/tools/mpi_wrapper_utils/gpu_tile_compact.sh\
  bash -c "<exec> inputs_weakscaling geometry.prob_extent=$probx $proby 1024  amr.n_cell=$nx $ny 256 amrex.the_arena_is_managed=0 amrex.use_gpu_aware_mpi=1" \
  > test_$NTOTRANKS.out

  lc=$(( lc+1  ))
  if [ $(( $lc%2 )) -eq 0 ]
  then
    probx=$(( probx*2  )); nx=$(( nx*2  ));
  else
    proby=$(( proby*2  )); ny=$(( ny*2  ));
  fi
done
