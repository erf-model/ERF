cmake_minimum_required(VERSION 3.24)

# Run one deck twice, once with the CPU default MFIter tiling and once with
# tiling switched off, then require the two plotfiles (3D and 2D) to agree.
#
# PBL schemes build per-tile planar work arrays; a kernel that loops over the
# valid box instead of the tile reads outside them.  Boxes wider than the tile
# size in y are the case that exposes this, and the untiled run is the
# reference, so no gold file is needed.

foreach(_required MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE INPUT
                  WORKING_DIRECTORY FCOMPARE RTOL ATOL PLTFILE PLT2DFILE)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "RunTilingParity.cmake requires ${_required}")
  endif()
endforeach()

set(_mpi_run "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "${NRANKS}")
set(_mpi_one "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" "1")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
  separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
  list(APPEND _mpi_run ${_mpi_preflags})
  list(APPEND _mpi_one ${_mpi_preflags})
endif()

set(_runtime_options)
if(DEFINED RUNTIME_OPTIONS AND NOT "${RUNTIME_OPTIONS}" STREQUAL "")
  separate_arguments(_runtime_options UNIX_COMMAND "${RUNTIME_OPTIONS}")
endif()

# The tiled run spells out the AMReX CPU default (1024000 8 8) so the test does
# not depend on it; the untiled run makes every tile its whole box.
set(_tiled_size   "1024000 8 8")
set(_untiled_size "1024000 1024000 1024000")

foreach(_run IN ITEMS tiled untiled)
  set(_tile_size "${_${_run}_size}")
  file(REMOVE_RECURSE "${WORKING_DIRECTORY}/${_run}_plt" "${WORKING_DIRECTORY}/${_run}_plt2d")
  execute_process(
    COMMAND ${_mpi_run} "${TEST_EXE}" "${INPUT}" ${_runtime_options}
            "fabarray.mfiter_tile_size=${_tile_size}"
            "erf.plot_file_1=${_run}_plt"
            "erf.plot2d_file_1=${_run}_plt2d"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/${_run}.log"
    ERROR_FILE "${WORKING_DIRECTORY}/${_run}.log"
    RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "${_run} run (fabarray.mfiter_tile_size=${_tile_size}) failed "
                        "with exit code ${_result}; see ${WORKING_DIRECTORY}/${_run}.log")
  endif()
endforeach()

foreach(_kind IN ITEMS plt plt2d)
  if(_kind STREQUAL "plt")
    set(_step "${PLTFILE}")
  else()
    set(_step "${PLT2DFILE}")
  endif()
  execute_process(
    COMMAND ${_mpi_one} "${FCOMPARE}" --abort_if_not_all_found -a
            -r "${RTOL}" --abs_tol "${ATOL}"
            "${WORKING_DIRECTORY}/untiled_${_kind}${_step}"
            "${WORKING_DIRECTORY}/tiled_${_kind}${_step}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/fcompare_${_kind}.log"
    ERROR_FILE "${WORKING_DIRECTORY}/fcompare_${_kind}.log"
    RESULT_VARIABLE _result)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "tiled and untiled ${_kind}${_step} differ; "
                        "see ${WORKING_DIRECTORY}/fcompare_${_kind}.log")
  endif()
endforeach()
