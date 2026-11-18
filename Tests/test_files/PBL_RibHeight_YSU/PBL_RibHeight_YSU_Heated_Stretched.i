# PBL_RibHeight_YSU_Heated on the stretched mesh with ratio 1 (same 32 m levels):
# the PBL height must be the constant-dz one, 230.5 m.  Before erf-model/ERF#5167
# it was 233.6 m, and 229.0 m with only the corrector's Rib measured from the ground.
FILE = PBL_RibHeight_YSU_Heated.i
erf.grid_stretching_ratio = 1.0
erf.initial_dz            = 32.0
