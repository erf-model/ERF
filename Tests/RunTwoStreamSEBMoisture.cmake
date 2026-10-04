# Run the TwoStream_SEBSurfaceLayerFluxes deck with the balance's skin and moisture both
# coupled to the surface layer (seb_surface_layer_uses_skin, seb_surface_layer_uses_moisture)
# on Noah-MP's silty clay loam (seb_soil_type = 8: wilting point 0.120, field capacity
# 0.387), three times: the soil at the wilting point (dry), at field capacity (wet) and in
# between (mid); and a fourth (veg) like mid with Noah-MP's grassland on 80 % of the surface
# (seb_vegetation_type = 10, LAI 2), so that canopy and soil resistances set the
# availability; then check with check_two_stream_seb_moisture.py.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FEXTRACT PYTHON_EXE CHECKER STEPS DT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTwoStreamSEBMoisture.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunTwoStreamSEBMoisture.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBMoisture.cmake: ERF executable")
erf_resolve_executable(FEXTRACT "${FEXTRACT}" CONFIG "${CONFIG}"
    CONTEXT "RunTwoStreamSEBMoisture.cmake: fextract")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTwoStreamSEBMoisture.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)

function(run_leg name q_s)
    # Extra inputs after q_s.
    set(dir "${WORKING_DIRECTORY}/${name}")
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT}
                max_step=${STEPS} erf.fixed_dt=${DT}
                erf.radiation.seb_surface_layer_uses_skin=true
                erf.radiation.seb_surface_layer_uses_moisture=true
                erf.radiation.seb_soil_type=8
                erf.radiation.seb_q_sfc_default=${q_s}
                erf.radiation.seb_q_deep_default=${q_s}
                "erf.plot2d_vars_1=seb_t_sfc seb_q_sfc seb_lh latent_heat_flux q_surf"
                ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunTwoStreamSEBMoisture.cmake: the ${name} run failed: ${_result} (see ${dir}/simulation.log)")
    endif()
endfunction()

run_leg(dry 0.120)
run_leg(wet 0.387)
run_leg(mid 0.25)
run_leg(veg 0.25 erf.radiation.seb_vegetation_type=10 erf.radiation.seb_vegetation_fraction=0.8
        erf.radiation.seb_leaf_area_index=2.0)

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" --fextract "${FEXTRACT}"
            --dry-dir "${WORKING_DIRECTORY}/dry" --wet-dir "${WORKING_DIRECTORY}/wet"
            --mid-dir "${WORKING_DIRECTORY}/mid" --veg-dir "${WORKING_DIRECTORY}/veg"
            --mid-q 0.25 --wilt 0.120 --fc 0.387
            --steps ${STEPS} --dt ${DT}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunTwoStreamSEBMoisture.cmake: the checker failed: ${check_result}")
endif()
