# PBL_SurfaceCrossing_YSU on a graded mesh: a 20 m first cell, and a stretching
# ratio (about 1.03) that ERF solves for so that 32 cells fill the 1024 m domain.
# The warm columns take the floor of the YSU PBL height, the top face of the first
# cell, 20 m (erf-model/ERF#5252).  The mid-point of the first two cell centers
# is 20.14 m, so the test tells the two apart.
FILE = PBL_SurfaceCrossing_YSU.i
erf.initial_dz = 20.0
