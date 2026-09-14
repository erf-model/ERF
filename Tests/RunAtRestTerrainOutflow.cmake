# Run the at-rest terrain deck twice, with symmetry x boundaries and with outflow x
# boundaries, and require both to stay at rest.
#
# The exact solution of this deck is that nothing happens: a dry, constant-theta atmosphere
# at rest over a hill, with no forcing, no diffusion and no damping.  Any w that appears is
# the scheme's own imbalance.
#
# Under symmetry the extension of the mesh past the domain is an exact mirror, so a ghost
# cell sits at the height of the cell it reflects and the base state there is right by
# construction.  That run is the control: it measures what the interior alone produces, and
# it is what makes a failure of the outflow run attributable to the boundary.
#
# Under outflow the mesh is extrapolated instead of mirrored, so a lateral ghost cell sits
# at a different height than the cell inside it and the base state there has to be built,
# not copied.  This is the run that guards that construction.
if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED FEXTREMA OR NOT DEFINED PLTFILE OR
   NOT DEFINED TOLERANCE OR NOT DEFINED GRADP_TOLERANCE)
    message(FATAL_ERROR "RunAtRestTerrainOutflow.cmake missing required argument")
endif()

# The largest |dp0/dx| the run reported at startup, over every x face of the domain
# including the two on the lateral boundaries.
#
# The base state is one reference atmosphere sampled on this mesh, so its gradient along a
# surface of constant height is a discretization error and nothing more -- on the lateral
# boundary faces as much as anywhere else.  A ghost cell whose base state belongs to some
# other height shows up here, and only here, as an extremum pinned to the boundary face and
# an order of magnitude above what the interior produces.  That is the signature this guards
# against, for both of the ways a lateral boundary extends the mesh past the domain.
function(max_abs_dp0dx RUN_NAME OUT_VAR)
    set(LOG "${WORKING_DIRECTORY}/${RUN_NAME}/simulation.log")
    if(NOT EXISTS "${LOG}")
        message(FATAL_ERROR "${RUN_NAME} wrote no log (${LOG})")
    endif()
    file(READ "${LOG}" log_text)
    if(log_text MATCHES "dp0/dx[ \t]+are[ \t]+zero")
        set(${OUT_VAR} "0" PARENT_SCOPE)
        return()
    endif()
    if(NOT log_text MATCHES "dp0/dx[ \t]+are[ \t]+([-+0-9.eE]+)[ \t]+([-+0-9.eE]+)")
        message(FATAL_ERROR "${RUN_NAME} reported no dp0/dx diagnostic; did erf.v get turned off?")
    endif()
    # Take both captures before the first string(REGEX), which overwrites CMAKE_MATCH_<n>
    set(G_LO "${CMAKE_MATCH_1}")
    set(G_HI "${CMAKE_MATCH_2}")
    string(REGEX REPLACE "^-" "" G_LO "${G_LO}")
    string(REGEX REPLACE "^-" "" G_HI "${G_HI}")
    if(G_LO GREATER G_HI)
        set(${OUT_VAR} "${G_LO}" PARENT_SCOPE)
    else()
        set(${OUT_VAR} "${G_HI}" PARENT_SCOPE)
    endif()
endfunction()

# The largest |w| anywhere in a run's plotfile
function(max_abs_w RUN_NAME OUT_VAR)
    set(PLT "${WORKING_DIRECTORY}/${RUN_NAME}/${PLTFILE}")
    if(NOT EXISTS "${PLT}")
        message(FATAL_ERROR "${RUN_NAME} wrote no plotfile (${PLT})")
    endif()
    execute_process(
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                ${FEXTREMA} -v "z_velocity" "${PLT}"
        OUTPUT_VARIABLE extrema_out
        ERROR_VARIABLE extrema_err
        RESULT_VARIABLE extrema_result)
    if(NOT extrema_result EQUAL 0)
        message(FATAL_ERROR "fextrema failed on ${PLT}: ${extrema_result}\n${extrema_err}")
    endif()
    if(NOT extrema_out MATCHES "z_velocity[ \t]+([-+0-9.eE]+)[ \t]+([-+0-9.eE]+)")
        message(FATAL_ERROR "Cannot read the z_velocity extrema from:\n${extrema_out}")
    endif()
    # Take both captures before the first string(REGEX), which overwrites CMAKE_MATCH_<n>
    set(W_LO "${CMAKE_MATCH_1}")
    set(W_HI "${CMAKE_MATCH_2}")
    string(REGEX REPLACE "^-" "" W_LO "${W_LO}")
    string(REGEX REPLACE "^-" "" W_HI "${W_HI}")
    if(W_LO GREATER W_HI)
        set(${OUT_VAR} "${W_LO}" PARENT_SCOPE)
    else()
        set(${OUT_VAR} "${W_HI}" PARENT_SCOPE)
    endif()
endfunction()

foreach(RUN IN ITEMS "symmetry|Symmetry" "outflow|Outflow")
    string(REPLACE "|" ";" RUN "${RUN}")
    list(GET RUN 0 RUN_NAME)
    list(GET RUN 1 BC_TYPE)

    set(RUN_DIR "${WORKING_DIRECTORY}/${RUN_NAME}")
    file(REMOVE_RECURSE "${RUN_DIR}")
    file(MAKE_DIRECTORY "${RUN_DIR}")

    execute_process(
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS} ${MPIEXEC_PREFLAGS}
                ${TEST_EXE} ${INPUT} "xlo.type=${BC_TYPE}" "xhi.type=${BC_TYPE}"
        WORKING_DIRECTORY "${RUN_DIR}"
        OUTPUT_FILE "${RUN_DIR}/simulation.log"
        ERROR_FILE "${RUN_DIR}/simulation.log"
        RESULT_VARIABLE run_result)
    if(NOT run_result EQUAL 0)
        message(FATAL_ERROR "${RUN_NAME} simulation failed: ${run_result}")
    endif()
endforeach()

max_abs_w(symmetry W_SYMMETRY)
max_abs_w(outflow  W_OUTFLOW)

max_abs_dp0dx(symmetry G_SYMMETRY)
max_abs_dp0dx(outflow  G_OUTFLOW)

file(WRITE "${WORKING_DIRECTORY}/at_rest.log"
     "max|w| with symmetry x boundaries      = ${W_SYMMETRY}\n"
     "max|w| with outflow  x boundaries      = ${W_OUTFLOW}\n"
     "tolerance                              = ${TOLERANCE}\n"
     "max|dp0/dx| with symmetry x boundaries = ${G_SYMMETRY}\n"
     "max|dp0/dx| with outflow  x boundaries = ${G_OUTFLOW}\n"
     "gradp tolerance                        = ${GRADP_TOLERANCE}\n")
message(STATUS "max|w| symmetry = ${W_SYMMETRY}, outflow = ${W_OUTFLOW}, tolerance = ${TOLERANCE}")
message(STATUS "max|dp0/dx| symmetry = ${G_SYMMETRY}, outflow = ${G_OUTFLOW}, tolerance = ${GRADP_TOLERANCE}")

# The control failing means the deck stopped measuring the boundary -- it picked up forcing,
# diffusion or a nonzero initial state -- rather than that the boundary went wrong.
if(W_SYMMETRY GREATER TOLERANCE)
    message(FATAL_ERROR
        "The symmetry control is not at rest: max|w| = ${W_SYMMETRY} exceeds ${TOLERANCE}. "
        "This deck no longer measures the lateral boundary treatment.")
endif()

if(W_OUTFLOW GREATER TOLERANCE)
    message(FATAL_ERROR
        "Outflow x boundaries are not well balanced at rest over terrain: max|w| = ${W_OUTFLOW} "
        "exceeds ${TOLERANCE}, against ${W_SYMMETRY} under symmetry.")
endif()

if(G_SYMMETRY GREATER GRADP_TOLERANCE)
    message(FATAL_ERROR
        "The base state has a horizontal gradient at a symmetry boundary: max|dp0/dx| = "
        "${G_SYMMETRY} exceeds ${GRADP_TOLERANCE}.  A symmetry boundary mirrors the mesh, so "
        "the ghost cell sits at the height of the cell it reflects and the reflected base "
        "state belongs there; a gradient this size means the ghost cell is holding the base "
        "state of some other height.")
endif()

if(G_OUTFLOW GREATER GRADP_TOLERANCE)
    message(FATAL_ERROR
        "The base state has a horizontal gradient at an outflow boundary: max|dp0/dx| = "
        "${G_OUTFLOW} exceeds ${GRADP_TOLERANCE}, against ${G_SYMMETRY} under symmetry.  An "
        "outflow boundary extrapolates the mesh, so the base state in the ghost cells has to "
        "be built at the height the mesh puts them at rather than copied outward.")
endif()
