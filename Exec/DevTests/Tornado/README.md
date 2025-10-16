This problem setup is for a tornado simulation from wrfinput and 
wrfbdy files 

# Compilation with GPUs on Perlmutter with radiation (RRTMGP)

```
    cd ERF
    source Build/GNU_Ekat/ekat_build_cuda_commands_Perlmutter.sh
    cd Exec/DevTests/Tornado
    cp GNUmakefile_GPU GNUmakefile
    make -j8
```

