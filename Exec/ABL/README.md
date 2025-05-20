# Atmospheric Boundary Layer
This problem setup is for simulation of the Atmospheric Boundary Layer (ABL)
using one of two turbulence schemes (Smagorinsky or Deardorff) and the bottom
boundary condition possibly specified by Monin Obukhov Similarity Theory (MOST).

This version of the ABL problem initializes the data using a hydrostatic profile
with random perturbations in velocity and potential temperature.


## Scaling studies
Scripts have been added in the `Scaling` directory to perform weak and strong 
scaling studies on GPUs on Perlmutter (NVIDIA A100) and Aurora (Intel PVC). The 
`inputs` files are set to not write output. If output needs to be written, change 
`erf.plot_int_1` parameter in the `inputs` file to the frequency of writing the 
output (ie. number of time steps)

### Perlmutter
1. To compile the code on Perlmutter NVIDIA A100 GPUs
```
cp Scaling/Perlmutter/GNUmakefile_GPU.perlmutter GNUmakefile
make -j8
```
2. Submit the batch job
 
Make sure to put the account id `ACCOUNT_ID` and the executable `<exec>` names in the 
job scripts (.qsub)
```
cd Perlmutter
sbatch WeakScaling_GPU.qsub
sbatch StrongScaling_GPU.qsub
```

3. For small scaling studies, interactive sessions can be requested on Perlmutter for upto 4 nodes (ie. 16 GPUs ( 4 GPUs per node))
```
cd Scaling/Perlmutter
sh get_interactive_nodes.sh
``` 

### Aurora

On Aurora, the scaling studies can be done on an interactive session on the `debug-scaling`
queue that gets upto 31 nodes (ie. 372 Intel GPUs (12 GPUs per node)).

1. To compile the code on Aurora Intel PVC GPUs
```
cp Scaling/Aurora/GNUmakefile_GPU.aurora GNUmakefile
make -j8
```

2. To get an interactive session
```
sh Scaling/Aurora/get_interactive_nodes.sh
```

3. Execute the scripts
```
cd Scaling/Aurora
sh WeakScaling_GPU.sh
sh StrongScaling_GPU.sh
```



