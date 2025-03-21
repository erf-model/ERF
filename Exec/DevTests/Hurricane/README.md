
# Hurricane simulations 

This folder contains examples for hurricane simulations from real weather data.

## Hurricane Henri

1. Download the initial condition file 

```
wget https://zenodo.org/record/15043093/files/ERF_IC_HenriERA5_20210819_VeryLarge.bin
```

2. `make -j8`
3. Run with `inputs_20210819_Henri_NoAMR`

## Hurricane Laura

1. Download the initial condition file 

```
wget https://zenodo.org/record/15043093/files/ERF_IC_Laura_LargeDomain.bin
```

2. `make -j8`.
3. Run with `inputs_20200826_Laura_NoAMR` or `inputs_20200826_Laura_2LevelAMR`.

