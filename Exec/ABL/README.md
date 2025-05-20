# Atmospheric Boundary Layer
This problem setup is for simulation of the Atmospheric Boundary Layer (ABL)
using one of two turbulence schemes (Smagorinsky or Deardorff) and the bottom
boundary condition possibly specified by Monin Obukhov Similarity Theory (MOST).

This version of the ABL problem initializes the data using a hydrostatic profile
with random perturbations in velocity and potential temperature.


## Scaling studies
Scripts have been added to perform scaling studies on GPUs on Perlmutter (NVIDIA A100) 
and Aurora (Intel PVC)

### Perlmutter
1. To compile the code on Perlmutter NVIDIA A100 GPUs
```
cp Scaling/Perlmutter/GNUmakefile_GPU.perlmutter GNUmakefile
make -j8
```
2. Submit the job 
Make sure to put the account id `ACCOUNT_ID` and the executable `<exec>` names in the 
job scripts (.qsub)
```
sbatch WeakScaling_GPU.qsub
sbatch StrongScaling_GPU.qsub
```

### Aurora


