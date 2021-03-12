#=============================================================================
# Regression tests
#=============================================================================

# # Run in CI
# if(ERF_DIM GREATER 1)
#   add_test_r(tg-1 TG)
#   add_test_r(tg-2 TG)
#   add_test_r(hit-1 HIT)
#   add_test_r(hit-2 HIT)
#   add_test_r(hit-3 HIT)
#   add_test_r(sod-1 Sod)
# endif()
# if(ERF_ENABLE_MASA)
#   add_test_r(mms-3 MMS)
#   if(ERF_DIM GREATER 1)
#     add_test_r(mms-1 MMS)
#     add_test_r(mms-2 MMS)
#     add_test_r(mms-4 MMS)
#     if(ERF_DIM GREATER 2)
#       add_test_r(mms-5 MMS)
#     endif()
#     if(ERF_ENABLE_AMREX_EB)
#       add_test_r(ebmms-1 EB_MMS)
#     endif()
#   endif()
# endif()
# 
# # Not run in CI
# if(ERF_DIM GREATER 1)
#   add_test_re(pmf-1 PMF)
#   add_test_re(multispecsod-1 MultiSpecSod)
#   if(ERF_DIM GREATER 2)
#     add_test_re(tg-3 TG)
#     add_test_re(tg-4 TG)
#   endif()
# endif()

#=============================================================================
# Verification tests
#=============================================================================
# if(ERF_ENABLE_MASA AND (ERF_DIM GREATER 2))
#   add_test_v1(symmetry MMS)
#   if(ERF_ENABLE_AMREX_EB)
#     add_test_v1(eb-symmetry EB_MMS)
#   endif()
# 
#   # Create list of resolutions we want to test with; make sure to pass it as a string in quotes
#   set(LIST_OF_GRID_SIZES 8 12 16 20)
#   add_test_v2(cns-no-amr MMS "${LIST_OF_GRID_SIZES}")
#   add_test_v2(cns-no-amr-mol MMS "${LIST_OF_GRID_SIZES}")
#   #add_test_v3(cns-amr MMS "${LIST_OF_GRID_SIZES}") # This one takes a while with AMR
# 
#   set(LIST_OF_GRID_SIZES 8 12 16 24)
#   add_test_v2(cns-les-no-amr MMS "${LIST_OF_GRID_SIZES}")
# endif()

#=============================================================================
# Unit tests
#=============================================================================
# add_test_u(unit_tests)

#=============================================================================
# Performance tests
#=============================================================================

