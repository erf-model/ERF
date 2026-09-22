cmake_minimum_required(VERSION 3.24)

include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")

# Run one deck twice, once with the CPU default MFIter tiling and once with
# tiling switched off, then require the two plotfiles (3D and 2D) to agree.
#
# PBL schemes build per-tile planar work arrays; a kernel that loops over the
# valid box instead of the tile reads outside them.  Boxes wider than the tile
# size in y are the case that exposes this, and the untiled run is the
# reference, so no gold file is needed.
#
# Two plotfiles that agree prove nothing if the field compared is one constant
# in every column (a PBL height pinned at its floor, or one that the scheme
# never writes), so VARYING_3D and VARYING_2D name fields that must take more
# than one value in the untiled plotfile, checked with fextrema.

# add_test_tiling_parity always passes every -D, so an unset CMake variable
# arrives here defined but empty; reject both.
foreach(_required MPIEXEC MPIEXEC_NUMPROC_FLAG NRANKS TEST_EXE INPUT
                  WORKING_DIRECTORY FCOMPARE RTOL ATOL PLTFILE PLT2DFILE)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "RunTilingParity.cmake requires ${_required}")
  endif()
endforeach()

set(_varying_plt)
set(_varying_plt2d)
if(DEFINED VARYING_3D AND NOT "${VARYING_3D}" STREQUAL "")
  separate_arguments(_varying_plt UNIX_COMMAND "${VARYING_3D}")
endif()
if(DEFINED VARYING_2D AND NOT "${VARYING_2D}" STREQUAL "")
  separate_arguments(_varying_plt2d UNIX_COMMAND "${VARYING_2D}")
endif()
if(_varying_plt OR _varying_plt2d)
  if(NOT DEFINED FEXTREMA OR "${FEXTREMA}" STREQUAL "")
    message(FATAL_ERROR "RunTilingParity.cmake requires FEXTREMA when VARYING_3D or VARYING_2D is set")
  endif()
  # FEXTREMA_EXE is a user-settable cache entry. Pointed at something that is
  # not fextrema (fcompare, say) every check below fails with "could not read
  # <field>", which reads as a broken deck rather than a misconfigured path.
  if(NOT EXISTS "${FEXTREMA}")
    message(FATAL_ERROR "RunTilingParity.cmake: FEXTREMA=${FEXTREMA} does not exist; "
                        "set the FEXTREMA_EXE cache entry to the amrex_fextrema built "
                        "alongside amrex_fcompare")
  endif()
endif()

# MPIEXEC may be a multi-word command such as "flux run", so it is split and
# validated by the shared helper rather than used as a single argv[0].
erf_mpi_launcher_command(_mpi_run
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTilingParity.cmake")
erf_mpi_launcher_command(_mpi_one
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTilingParity.cmake")

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
  # Plotfiles carry the step suffix (tiled_plt00010, tiled_plt2d00010, ...).
  # Remove them all so a stale step-PLTFILE pair from an earlier run cannot be
  # compared if this run no longer writes that step.
  file(GLOB _stale "${WORKING_DIRECTORY}/${_run}_plt*")
  if(_stale)
    file(REMOVE_RECURSE ${_stale})
  endif()
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

# The agreement above is only evidence if the tiled run is in fact tiled.
# FabArrayBase::buildTileArray splits a direction into max(ncells/tilesize,1)
# tiles, so a box splits only once it is twice the tile size: with the 8-cell
# default a 16-cell box gives two tiles and a 15-cell box gives one. A deck
# edit that narrows the boxes -- a smaller amr.max_grid_size_y, a shorter
# domain -- therefore leaves one tile per box, and every test here keeps
# passing while comparing two identical untiled runs.
#
# Only x and y are counted. The schemes this harness covers iterate under
# TileNoZ(), which pins the z tile size to the whole box; counting z would let
# z-tiling alone satisfy the check after the y split that exposes the bug had
# quietly gone away.
function(count_xy_tiles plotfile tile_size out_var)
  set(_cell_h "${plotfile}/Level_0/Cell_H")
  if(NOT EXISTS "${_cell_h}")
    message(FATAL_ERROR "RunTilingParity.cmake: no ${_cell_h}")
  endif()
  separate_arguments(_ts UNIX_COMMAND "${tile_size}")
  list(GET _ts 0 _ts_x)
  list(GET _ts 1 _ts_y)
  # The box count on level 0 is the first number on the fifth line of Cell_H,
  # followed by one "((lo) (hi) (typ))" line per box.
  file(STRINGS "${_cell_h}" _lines)
  list(GET _lines 4 _boxes_line)
  if(NOT _boxes_line MATCHES "^\\(([0-9]+) ")
    message(FATAL_ERROR "RunTilingParity.cmake: cannot read the box count from ${_cell_h}: '${_boxes_line}'")
  endif()
  set(_nboxes "${CMAKE_MATCH_1}")
  set(_seen 0)
  set(_tiles 0)
  foreach(_line IN LISTS _lines)
    if(_line MATCHES "^\\(\\((-?[0-9]+),(-?[0-9]+),(-?[0-9]+)\\) \\((-?[0-9]+),(-?[0-9]+),(-?[0-9]+)\\)")
      math(EXPR _nx "${CMAKE_MATCH_4} - ${CMAKE_MATCH_1} + 1")
      math(EXPR _ny "${CMAKE_MATCH_5} - ${CMAKE_MATCH_2} + 1")
      math(EXPR _tx "${_nx} / ${_ts_x}")
      math(EXPR _ty "${_ny} / ${_ts_y}")
      if(_tx LESS 1)
        set(_tx 1)
      endif()
      if(_ty LESS 1)
        set(_ty 1)
      endif()
      math(EXPR _tiles "${_tiles} + ${_tx} * ${_ty}")
      math(EXPR _seen "${_seen} + 1")
    endif()
  endforeach()
  if(NOT _seen EQUAL _nboxes)
    message(FATAL_ERROR "RunTilingParity.cmake: ${_cell_h} declares ${_nboxes} boxes on level 0 "
                        "but ${_seen} could be parsed")
  endif()
  set(${out_var} ${_tiles} PARENT_SCOPE)
endfunction()
count_xy_tiles("${WORKING_DIRECTORY}/tiled_plt${PLTFILE}"   "${_tiled_size}"   _tiled_tiles)
count_xy_tiles("${WORKING_DIRECTORY}/tiled_plt${PLTFILE}"   "${_untiled_size}" _untiled_tiles)
if(NOT _tiled_tiles GREATER _untiled_tiles)
  message(FATAL_ERROR "RunTilingParity.cmake: fabarray.mfiter_tile_size=${_tiled_size} gives "
                      "${_tiled_tiles} tiles in x-y over the level-0 BoxArray and the untiled run "
                      "${_untiled_tiles}; the tiled run must have more, or the comparison is "
                      "between two untiled runs. Widen the boxes in y (amr.max_grid_size_y) in "
                      "the deck")
endif()
message(STATUS "RunTilingParity: ${_untiled_tiles} boxes on level 0, ${_tiled_tiles} x-y tiles "
               "with fabarray.mfiter_tile_size=${_tiled_size}")

foreach(_kind IN ITEMS plt plt2d)
  if(_kind STREQUAL "plt")
    set(_step "${PLTFILE}")
  else()
    set(_step "${PLT2DFILE}")
  endif()
  # execute_process rewrites the run and fcompare logs, but the fextrema output
  # below is appended one field at a time. Start each run from empty so the
  # ATTACHED_FILES_ON_FAIL artifact holds this run's ranges and not a pile of
  # earlier ones from the same build directory.
  file(WRITE "${WORKING_DIRECTORY}/fextrema_${_kind}.log" "")
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

  # The agreement above is only evidence if the fields vary: require min < max
  # in the untiled reference for every field the test claims to cover.
  foreach(_var IN LISTS _varying_${_kind})
    execute_process(
      COMMAND ${_mpi_one} "${FEXTREMA}" -v "${_var}"
              "${WORKING_DIRECTORY}/untiled_${_kind}${_step}"
      WORKING_DIRECTORY "${WORKING_DIRECTORY}"
      OUTPUT_VARIABLE _extrema
      ERROR_VARIABLE _extrema_err
      RESULT_VARIABLE _result)
    file(APPEND "${WORKING_DIRECTORY}/fextrema_${_kind}.log" "${_extrema}${_extrema_err}")
    # Some launchers label task output with a per-line "<rank>: " prefix (srun
    # --label, flux --label-io). Strip it, normalising to one leading space;
    # unlabelled output is unchanged.
    string(REGEX REPLACE "\n[0-9]+:[ \t]*" "\n " _extrema "\n${_extrema}")
    # fextrema prints one line per variable: " name   min   max"
    if(NOT _result EQUAL 0 OR
       NOT "${_extrema}" MATCHES "\n ${_var}[ \t]+([^ \t\n]+)[ \t]+([^ \t\n]+)[ \t]*\n")
      message(FATAL_ERROR "fextrema could not read ${_var} from untiled_${_kind}${_step}; "
                          "see ${WORKING_DIRECTORY}/fextrema_${_kind}.log")
    endif()
    set(_min "${CMAKE_MATCH_1}")
    set(_max "${CMAKE_MATCH_2}")
    if(NOT _min LESS _max)
      message(FATAL_ERROR "${_var} is uniform (${_min}) in untiled_${_kind}${_step}, so its "
                          "tiled/untiled agreement is vacuous; change the deck so it varies")
    endif()
    message(STATUS "${_kind}${_step} ${_var} ranges ${_min} .. ${_max}")
  endforeach()
endforeach()
