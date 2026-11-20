# PBL_SurfaceCrossing_MRF on a graded mesh: a 20 m first cell, and a stretching
# ratio (about 1.03) that ERF solves for so that 32 cells fill the 1024 m domain.
# The warm columns take the floor of the MRF PBL height, the top face of the first
# cell, 20 m.  Before, the floor was the first cell center, 10 m, and these columns
# reported a crossing just above it, at 11.8 m.
FILE = PBL_SurfaceCrossing_MRF.i
erf.initial_dz = 20.0
