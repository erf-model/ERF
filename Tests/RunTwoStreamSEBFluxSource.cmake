# Run one two-stream deck three times -- erf.radiation.seb_turbulent_flux_source =
# surface_layer, = defaults, and surface_layer with seb_surface_layer_uses_skin = true
# (two-way) -- then check with check_two_stream_seb_flux_source.py that the surface
# energy balance removed the surface layer's H and LE from the ground (the balance's
# fluxes equal the surface layer's 2D outputs at every step, and the skin ends cooler
# than with the defaults by the energy they carried away), and that in the two-way run
# the surface layer's surface temperature is the balance's skin.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FEXTRACT PYTHON_EXE CHECKER
            STEPS DT HEAT_CAPACITY)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamSEBFluxSource.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamSEBFluxSource.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBFluxSource.cmake: ERF executable")
erf_resolve_executable(FEXTRACT "${FEXTRACT}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBFluxSource.cmake: fextract")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamSEBFluxSource.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)
set(SL_DIR  "${WORKING_DIRECTORY}/surface_layer")
set(DEF_DIR "${WORKING_DIRECTORY}/defaults")
set(TWO_DIR "${WORKING_DIRECTORY}/two_way")
file(REMOVE_RECURSE "${SL_DIR}" "${DEF_DIR}" "${TWO_DIR}")
file(MAKE_DIRECTORY "${SL_DIR}" "${DEF_DIR}" "${TWO_DIR}")

function(run_leg dir source)
    # The deck reads its sounding by a relative name.
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT}
                erf.radiation.seb_turbulent_flux_source=${source}
                max_step=${STEPS} erf.fixed_dt=${DT}
                erf.radiation.seb_surface_heat_capacity=${HEAT_CAPACITY}
                ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamSEBFluxSource.cmake: the ${source} run failed: ${_result} (see ${dir}/simulation.log)")
    endif()
endfunction()

# Extra checker arguments, e.g. --multilevel for a refined deck.
separate_arguments(checker_options UNIX_COMMAND "${CHECKER_OPTIONS}")

run_leg("${SL_DIR}" surface_layer)
run_leg("${DEF_DIR}" defaults)
run_leg("${TWO_DIR}" surface_layer erf.radiation.seb_surface_layer_uses_skin=true)

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}"
            --fextract "${FEXTRACT}"
            --surface-layer-dir "${SL_DIR}"
            --defaults-dir "${DEF_DIR}"
            --two-way-dir "${TWO_DIR}"
            --steps ${STEPS} --dt ${DT} --heat-capacity ${HEAT_CAPACITY}
            ${checker_options}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBFluxSource.cmake: the checker failed: ${check_result}")
endif()
