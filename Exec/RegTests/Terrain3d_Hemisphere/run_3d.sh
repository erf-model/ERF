rm -rf ERF3d.gnu.TEST.MPI.ex
rm -rf plt*
rm -rf Backtrace.*
make -j
mpirun -np 8 ./ERF3d.gnu.TEST.MPI.ex inputs

