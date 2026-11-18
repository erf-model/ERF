# Run the same flat-terrain RRTMGP/SLM case with terrain correction off and on.
# With a fixed sun and flat ground, factor == 1 must preserve the raw flux
# exactly, including the seven fields passed through lsm_output_ptrs.
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY FCOMPARE FEXTREMA RRTMGP_ROOT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunTopographicRadiationFlat.cmake: ${arg} is required")
    endif()
endforeach()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunTopographicRadiationFlat.cmake: ERF executable")
erf_resolve_executable(FCOMPARE "${FCOMPARE}" CONFIG "${CONFIG}"
    CONTEXT "RunTopographicRadiationFlat.cmake: fcompare")
erf_resolve_executable(FEXTREMA "${FEXTREMA}" CONFIG "${CONFIG}"
    CONTEXT "RunTopographicRadiationFlat.cmake: fextrema")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS} PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTopographicRadiationFlat.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1 PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunTopographicRadiationFlat.cmake")

set(off "${WORKING_DIRECTORY}/topo_off")
set(on  "${WORKING_DIRECTORY}/topo_on")
file(REMOVE_RECURSE "${off}" "${on}")
foreach(dir "${off}" "${on}")
    file(MAKE_DIRECTORY "${dir}/radiation_data")
    file(COPY "${WORKING_DIRECTORY}/input_sounding" DESTINATION "${dir}")
    foreach(name rrtmgp-data-sw-g224-2018-12-04.nc rrtmgp-data-lw-g256-2018-12-04.nc)
        file(CREATE_LINK "${RRTMGP_ROOT}/rrtmgp/data/${name}"
            "${dir}/radiation_data/${name}" SYMBOLIC RESULT link_result)
        if(NOT link_result STREQUAL "0")
            message(FATAL_ERROR "Cannot link RRTMGP data ${name}: ${link_result}")
        endif()
    endforeach()
    foreach(name rrtmgp-cloud-optics-coeffs-sw.nc rrtmgp-cloud-optics-coeffs-lw.nc)
        file(CREATE_LINK "${RRTMGP_ROOT}/extensions/cloud_optics/${name}"
            "${dir}/radiation_data/${name}" SYMBOLIC RESULT link_result)
        if(NOT link_result STREQUAL "0")
            message(FATAL_ERROR "Cannot link RRTMGP data ${name}: ${link_result}")
        endif()
    endforeach()
endforeach()

execute_process(COMMAND ${launch} ${TEST_EXE} ${INPUT}
    erf.topo_shading=false erf.slope_rad=false
    WORKING_DIRECTORY "${off}"
    OUTPUT_FILE "${off}/simulation.log" ERROR_FILE "${off}/simulation.log"
    TIMEOUT 300 RESULT_VARIABLE off_result)
if(NOT off_result EQUAL 0)
    message(FATAL_ERROR "Uncorrected radiation run failed (${off_result}); see ${off}/simulation.log")
endif()
execute_process(COMMAND ${launch} ${TEST_EXE} ${INPUT}
    erf.topo_shading=true erf.slope_rad=true
    WORKING_DIRECTORY "${on}"
    OUTPUT_FILE "${on}/simulation.log" ERROR_FILE "${on}/simulation.log"
    TIMEOUT 300 RESULT_VARIABLE on_result)
if(NOT on_result EQUAL 0)
    message(FATAL_ERROR "Topographic radiation run failed (${on_result}); see ${on}/simulation.log")
endif()

# The first call, the intervening held-forcing step, and the next scheduled
# RRTMGP call must preserve the atmosphere, land forcing, and SLM state exactly.
foreach(step 00001 00002 00003)
    foreach(name "plt${step}" "plt_lsm_${step}" "plt_lsm_2D_${step}")
        foreach(dir "${off}" "${on}")
            if(NOT EXISTS "${dir}/${name}/Header")
                message(FATAL_ERROR "Missing ${name} from ${dir}")
            endif()
        endforeach()
        execute_process(COMMAND ${launch_one} ${FCOMPARE}
            --abort_if_not_all_found --rel_tol 0 --abs_tol 0
            "${off}/${name}" "${on}/${name}"
            OUTPUT_FILE "${WORKING_DIRECTORY}/${name}_parity.log"
            ERROR_FILE "${WORKING_DIRECTORY}/${name}_parity.log"
            RESULT_VARIABLE compare_result)
        if(NOT compare_result EQUAL 0)
            message(FATAL_ERROR "Flat-terrain parity failed for ${name}; see ${name}_parity.log")
        endif()
    endforeach()

    # The active path must actually run. The raw direct-beam factor is exactly
    # one on flat terrain; a missing diagnostic would instead be -999.
    foreach(var topographic_shadow terrain_direct_factor unshaded_sw_flux_dn terrain_sw_flux_dn)
        execute_process(COMMAND ${launch_one} ${FEXTREMA} -v "${var}"
            "${on}/terrain${step}"
            OUTPUT_VARIABLE extrema ERROR_VARIABLE extrema_error
            RESULT_VARIABLE extrema_result)
        if(NOT extrema_result EQUAL 0 OR
           NOT extrema MATCHES "${var}[ \t]+([-+0-9.eE]+)[ \t]+([-+0-9.eE]+)")
            message(FATAL_ERROR "Cannot read ${var} at step ${step}: ${extrema}${extrema_error}")
        endif()
        set(vmin "${CMAKE_MATCH_1}")
        set(vmax "${CMAKE_MATCH_2}")
        if(var STREQUAL "topographic_shadow" AND (NOT vmin EQUAL 0 OR NOT vmax EQUAL 0))
            message(FATAL_ERROR "Flat terrain has a shadow at step ${step}: [${vmin},${vmax}]")
        elseif(var STREQUAL "terrain_direct_factor" AND (NOT vmin EQUAL 1 OR NOT vmax EQUAL 1))
            message(FATAL_ERROR "Flat terrain changed the direct factor at step ${step}: [${vmin},${vmax}]")
        elseif((var STREQUAL "unshaded_sw_flux_dn" OR var STREQUAL "terrain_sw_flux_dn")
               AND (NOT vmin GREATER 0 OR NOT vmax GREATER 0))
            message(FATAL_ERROR "No positive ${var} at step ${step}: [${vmin},${vmax}]")
        endif()
    endforeach()

    # Parity alone could pass if both routes left SLM forcing at zero.
    foreach(var SW_dw_dir_vis LW_dw cos_zenith)
        execute_process(COMMAND ${launch_one} ${FEXTREMA} -v "${var}"
            "${on}/plt_lsm_${step}"
            OUTPUT_VARIABLE extrema ERROR_VARIABLE extrema_error
            RESULT_VARIABLE extrema_result)
        if(NOT extrema_result EQUAL 0 OR
           NOT extrema MATCHES "${var}[ \t]+([-+0-9.eE]+)[ \t]+([-+0-9.eE]+)")
            message(FATAL_ERROR "Cannot read SLM ${var} at step ${step}: ${extrema}${extrema_error}")
        endif()
        set(vmax "${CMAKE_MATCH_2}")
        if(NOT vmax GREATER 0)
            message(FATAL_ERROR "SLM received no positive ${var} at step ${step}")
        endif()
    endforeach()
endforeach()

foreach(dir "${off}" "${on}")
    file(STRINGS "${dir}/simulation.log" updates REGEX "Radiation advancing level 0")
    list(LENGTH updates update_count)
    if(NOT update_count EQUAL 2)
        message(FATAL_ERROR "Expected two RRTMGP calls in ${dir}; found ${update_count}")
    endif()
endforeach()
