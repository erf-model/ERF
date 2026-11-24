# PBL_SurfaceCrossing_YSU on a "stretched" mesh with a stretching ratio of 1: the
# same 32 m levels, through the terrain-fitted code path (MeshType::StretchedDz).
# The warm columns take the floor of the YSU PBL height, the top face of the first
# cell, 32 m, as on the constant-dz mesh (erf-model/ERF#5252).
FILE = PBL_SurfaceCrossing_YSU.i
erf.grid_stretching_ratio = 1.0
erf.initial_dz            = 32.0
