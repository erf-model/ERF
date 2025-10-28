# Column diffusion test

This simple problem was used to test implicit vertical diffusion.

Problem: Horizontally homogeneous (periodic boundaries) with a neutrally stratified layer below a
capping inversion and a stably stratified layer aloft; diffusivity is constant throughout. Three
sets of BCs were tested:

* Neumann on zlo, specified theta_grad on zhi
* Specified theta_grad on zlo, zhi
* Surface layer on zlo, specified theta_grad on zhi

Summary of findings

* The original fully explicit setup requires a diffusive CFL < 0.5 for numerical stability
* No-terrain and stretched grid give the same level of error compared to explicit (in terms of MAE
  and RMSE)
* Implicit vertical allows diffusive CFL ~ 50 (dt=60 s) for this test problem, above which the
  solution remains numerically stable but error begins to increase significantly
* To retain full temporal accuracy -- and minimize error -- use explicit diffusion in the final RK
  stage (e.g., with `erf.vert_implicit_fac = 1 1 0`)


## "N" path

Compare with:
```
erf.fixed_dt           = 0.5
erf.fixed_mri_dt_ratio = 4
erf.vert_implicit_fac  = 0.0
```

A diffusive CFL of 1 corresponds to dt = dz*dz/K = 1.33 s. This numerically
stable setup with explicit vertical diffusion corresponds to a CFL of ~0.4.

MAE calculated based on mean profiles of theta, output every 5 min for an hour of runtime

- [x] fully implicit, dt=0.5, `zhi.theta_grad = 0` : MAE = 5e-5
- [x] fully implicit, dt=0.5, `zhi.theta_grad = 0.003` : MAE = 4e-5
- [x] semi-implicit (fac=0.5), dt=0.5, `zhi.theta_grad = 0.003`: MAE = 2e-5
- [x] semi-implicit (fac=0.5), dt=1.0, `zhi.theta_grad = 0.003`: MAE = 4e-5
- [x] semi-implicit (fac=0.5), dt=2.0, `zhi.theta_grad = 0.003`: MAE = 9e-5
- [x] semi-implicit (fac=0.5), dt=4.0, `zhi.theta_grad = 0.003`: MAE = 2e-4
- [x] semi-implicit (fac=0.5), dt=8.0, `zhi.theta_grad = 0.003`: MAE = 3e-4
- [x] semi-implicit (fac=0.5), dt=16.0, `zhi.theta_grad = 0.003`: MAE = 5e-4
- [x] semi-implicit (fac=0.5), dt=30.0, `zhi.theta_grad = 0.003`: MAE = 1e-3, RMSE=2e-3
- [x] semi-implicit (fac=0.5), dt=60.0, `zhi.theta_grad = 0.003`: MAE = 3e-3, RMSE=5e-3 (numerical oscillations)
- [x] semi-implicit (fac=0.5), dt=120.0, `zhi.theta_grad = 0.003`: MAE = 7e-3, RMSE=1e-2 (numerical oscillations)
- [x] semi-implicit (fac=0.5), dt=240.0, `zhi.theta_grad = 0.003`: MAE = 1e-2, RMSE=3e-2 (numerical oscillations)
- [x] fully implicit, dt=60.0, `zhi.theta_grad = 0.003`: MAE = 5e-3, RMSE=7e-3
- [x] fully implicit, dt=120.0, `zhi.theta_grad = 0.003`: MAE = 9e-3, RMSE=1e-2 (deviation from explicit solution)
- [x] fully implicit, dt=240.0, `zhi.theta_grad = 0.003`: MAE = 1e-2, RMSE=2e-2 (deviation from explicit solution)
- [x] fully implicit, dt=300.0, `zhi.theta_grad = 0.003`: MAE = 1e-2, RMSE=2e-2 (deviation from explicit solution)


## "S" path

Same as previous test, with 
```
amr.n_cell           =   4    4   50

erf.initial_dz = 10.0
erf.grid_stretching_ratio = 1.025

zhi.theta_grad = 0.003
```

- [x] semi-implicit (fac=0.5), dt=30, `zhi.theta_grad = 0.003`: MAE = 1e-3, RMSE=2e-3
- [x] fully implicit, dt=60, `zhi.theta_grad = 0.003`: MAE = 5e-3, RMSE=7e-3
- [x] mixed implicit/explicit (fac=1 1 0), dt=60, `zhi.theta_grad = 0.003`: MAE = 2e-4, RMSE=1e-3

