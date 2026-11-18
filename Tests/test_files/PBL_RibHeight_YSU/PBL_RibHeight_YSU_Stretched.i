# PBL_RibHeight_YSU on a "stretched" mesh with a stretching ratio of 1: the same
# 32 m levels, through the terrain-fitted code path (MeshType::StretchedDz).
# The heights above the ground are the same as on the constant-dz mesh, so the
# YSU PBL height must be the same, 146.3 m.  Before erf-model/ERF#5167 this
# path measured the height in the bulk Richardson number from the first cell
# center instead of the ground and gave 153.6 m.
FILE = PBL_RibHeight_YSU.i
erf.grid_stretching_ratio = 1.0
erf.initial_dz            = 32.0
