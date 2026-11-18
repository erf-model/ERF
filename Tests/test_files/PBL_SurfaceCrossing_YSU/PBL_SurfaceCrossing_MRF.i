# PBL_SurfaceCrossing_YSU with the MRF scheme: the MRF PBL height of a column
# whose surface cell already has a bulk Richardson number at or above Ribcr.
#
# The warm columns of that deck have a surface-cell Rib of about 0.8, above the
# MRF Ribcr = 0.5: they take the floor of the MRF PBL height, the top face of
# the first cell, max(dz, 10 m) = 32 m, the same floor as YSU
# (erf-model/ERF#5252).  Before, the MRF floor was the first cell center,
# max(0.5 dz, 10 m) = 16 m.  The cold columns take 0.9 z(khi) = 907.2 m.  The
# CTest entry requires the pblh of the 2D plotfile to reach both values.
FILE = PBL_SurfaceCrossing_YSU.i
erf.pbl_type            = "MRF"
erf.most.pblh_calc      = "MRF"
