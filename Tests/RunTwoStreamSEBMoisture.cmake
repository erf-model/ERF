# Run the TwoStream_SEBSurfaceLayerFluxes deck with the balance's skin and moisture both
# coupled to the surface layer (seb_surface_layer_uses_skin, seb_surface_layer_uses_moisture)
# and check the legs with check_two_stream_seb_moisture.py:
#   dry, wet, mid  the linear soil-water factor alone (wilting point 0.120 and field
#                  capacity 0.387 given, no soil type), the soil at the wilting point, at
#                  field capacity and in between;
#   bare           mid's soil as Noah-MP's silty clay loam (seb_soil_type = 8, the same
#                  wilting point and field capacity), so its bare-soil resistance sets beta;
#   bare_f0        bare with grassland at a vegetated fraction of 0, which must give bare's
#                  surface mixing ratio exactly;
#   bare_dry       bare at the wilting point, where the pore air is nearly dry and the soil
#                  must not evaporate (LE <= 0, as in Noah-MP);
#   veg            bare with Noah-MP's grassland on 80 % of the surface
#                  (seb_vegetation_type = 10, LAI 2), canopy and soil resistances;
#   tables         veg for one step from a copy of the deck without erf.most.z0, so that the
#                  surface layer takes its land roughness from Noah-MP's tables and job_info
#                  must record that value.
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
# run_leg reads the deck from leg_input (the caller's, at the call).
set(leg_input "${INPUT}")

function(run_leg name q_s)
    # Extra inputs after q_s.
    set(dir "${WORKING_DIRECTORY}/${name}")
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    file(COPY "${input_dir}/input_sounding" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${leg_input}
                max_step=${STEPS} erf.fixed_dt=${DT}
                erf.radiation.seb_surface_layer_uses_skin=true
                erf.radiation.seb_surface_layer_uses_moisture=true
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

set(linear erf.radiation.seb_soil_moisture_wilt=0.120 erf.radiation.seb_soil_moisture_fc=0.387)
set(soil erf.radiation.seb_soil_type=8)
run_leg(dry 0.120 ${linear})
run_leg(wet 0.387 ${linear})
run_leg(mid 0.25 ${linear})
run_leg(bare 0.25 ${soil})
run_leg(bare_dry 0.120 ${soil})
run_leg(bare_f0 0.25 ${soil} erf.radiation.seb_vegetation_type=10
        erf.radiation.seb_vegetation_fraction=0.0 erf.radiation.seb_leaf_area_index=2.0)
run_leg(veg 0.25 ${soil} erf.radiation.seb_vegetation_type=10 erf.radiation.seb_vegetation_fraction=0.8
        erf.radiation.seb_leaf_area_index=2.0)

# The deck without its erf.most.z0 line: the tables' 0.8 x 0.12 + 0.2 x 0.002 = 0.0964 m.
file(READ "${INPUT}" deck)
string(REGEX REPLACE "\nerf\\.most\\.z0[ \t]*=[^\n]*" "\n" deck_no_z0 "${deck}")
if(deck_no_z0 STREQUAL deck)
    message(FATAL_ERROR "RunTwoStreamSEBMoisture.cmake: no erf.most.z0 line to remove in ${INPUT}")
endif()
set(leg_input "${WORKING_DIRECTORY}/deck_without_z0.i")
file(WRITE "${leg_input}" "${deck_no_z0}")
run_leg(tables 0.25 ${soil} erf.radiation.seb_vegetation_type=10
        erf.radiation.seb_vegetation_fraction=0.8 erf.radiation.seb_leaf_area_index=2.0 max_step=1)

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" --fextract "${FEXTRACT}"
            --dry-dir "${WORKING_DIRECTORY}/dry" --wet-dir "${WORKING_DIRECTORY}/wet"
            --mid-dir "${WORKING_DIRECTORY}/mid" --veg-dir "${WORKING_DIRECTORY}/veg"
            --bare-dir "${WORKING_DIRECTORY}/bare" --bare-f0-dir "${WORKING_DIRECTORY}/bare_f0"
            --bare-dry-dir "${WORKING_DIRECTORY}/bare_dry"
            --tables-dir "${WORKING_DIRECTORY}/tables" --tables-z0 0.0964 --given-z0 0.1
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
