rm -rf plt* chk*
make -j
mpirun -np 4 ./ERF3d.gnu.TEST.MPI.ex inputs_moisture
