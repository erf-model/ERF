include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY CHECKER)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunEBSurfaceLayerActivation.cmake: ${arg} is required")
    endif()
endforeach()
erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}")
erf_resolve_executable(CHECKER "${CHECKER}" CONFIG "${CONFIG}")
erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS} PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunEBSurfaceLayerActivation.cmake")
erf_mpi_launcher_command(launch_one
    LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS 1 PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunEBSurfaceLayerActivation.cmake checker")

function(run_leg name)
    set(dir "${WORKING_DIRECTORY}/${name}")
    file(REMOVE_RECURSE "${dir}")
    file(MAKE_DIRECTORY "${dir}")
    execute_process(COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log" ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 120 RESULT_VARIABLE rc)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "EB SurfaceLayer ${name} failed (${rc}); see ${dir}/simulation.log")
    endif()
endfunction()

function(check name)
    execute_process(COMMAND ${launch_one} ${CHECKER} ${ARGN}
        WORKING_DIRECTORY "${WORKING_DIRECTORY}"
        OUTPUT_FILE "${WORKING_DIRECTORY}/${name}.log"
        ERROR_FILE "${WORKING_DIRECTORY}/${name}.log"
        TIMEOUT 60 RESULT_VARIABLE rc)
    if(NOT rc EQUAL 0)
        message(FATAL_ERROR "EB SurfaceLayer ${name} failed (${rc}); see ${WORKING_DIRECTORY}/${name}.log")
    endif()
endfunction()

run_leg(embedded)
run_leg(physical "zlo.type=surface_layer")
run_leg(inactive "erf.eb_boundary_type=SlipWall")
check(initial state "${WORKING_DIRECTORY}/embedded/chk00001" 0 initial)
check(physical_parity parity "${WORKING_DIRECTORY}/embedded/plt00004"
    "${WORKING_DIRECTORY}/physical/plt00004")
check(drag drag "${WORKING_DIRECTORY}/embedded/plt00004"
    "${WORKING_DIRECTORY}/inactive/plt00004")

# The state filenames must also permit changing the covered physical BC on restart.
run_leg(physical_restart "erf.restart=${WORKING_DIRECTORY}/physical/chk00002")
check(physical_restart_parity parity "${WORKING_DIRECTORY}/embedded/plt00004"
    "${WORKING_DIRECTORY}/physical_restart/plt00004")

# Generate complete node-aligned maps without checking large fixture files in.
# Shift the fine-file transition from x=16 to x=17 to detect wrong file selection.
foreach(lev RANGE 0 1)
    math(EXPR spacing "2 - ${lev}")
    math(EXPR last "16 * (1 + ${lev}) - 1")
    math(EXPR transition "16 + ${lev}")
    set(varying "")
    set(constant "")
    foreach(j RANGE 0 ${last})
        math(EXPR y "${j} * ${spacing}")
        foreach(i RANGE 0 ${last})
            math(EXPR x "${i} * ${spacing}")
            if(x LESS transition)
                set(z0 "0.01")
            else()
                set(z0 "0.1")
            endif()
            string(APPEND varying "${x} ${y} ${z0}\n")
            if(lev EQUAL 0)
                string(APPEND constant "${x} ${y} 0.1\n")
            endif()
        endforeach()
    endforeach()
    set(map${lev} "${WORKING_DIRECTORY}/roughness_lev${lev}.txt")
    file(WRITE "${map${lev}}" "${varying}")
    if(lev EQUAL 0)
        set(constant_map "${WORKING_DIRECTORY}/roughness_constant.txt")
        file(WRITE "${constant_map}" "${constant}")
    endif()
endforeach()

run_leg(constant_file "erf.most.roughness_file_name=${constant_map}")
check(constant_file_parity parity "${WORKING_DIRECTORY}/embedded/plt00004"
    "${WORKING_DIRECTORY}/constant_file/plt00004")

run_leg(varying "erf.most.roughness_file_name=${map0}")
check(varying_initial roughness "${WORKING_DIRECTORY}/varying/chk00001" 0 initial)
check(varying_final roughness "${WORKING_DIRECTORY}/varying/chk00004" 0)
check(varying_drag regional_drag "${WORKING_DIRECTORY}/varying/plt00004"
    "${WORKING_DIRECTORY}/inactive/plt00004")
run_leg(varying_restart "erf.most.roughness_file_name=${map0}"
    "erf.restart=${WORKING_DIRECTORY}/varying/chk00002")
check(varying_restart_state roughness "${WORKING_DIRECTORY}/varying_restart/chk00004" 0)
check(varying_restart_parity parity "${WORKING_DIRECTORY}/varying/plt00004"
    "${WORKING_DIRECTORY}/varying_restart/plt00004")

foreach(mapping per_level interpolated)
    if(mapping STREQUAL "per_level")
        set(files "${map0} ${map1}")
        set(check_option "")
    else()
        set(files "${map0}")
        set(check_option interpolated)
    endif()
    run_leg(${mapping} "amr.max_level=1" "erf.most.roughness_file_name=${files}")
    if(EXISTS "${WORKING_DIRECTORY}/${mapping}/chk00001/Level_1" OR
       NOT EXISTS "${WORKING_DIRECTORY}/${mapping}/chk00003/Level_1/Z0_H")
        message(FATAL_ERROR "The roughness AMR leg must create level 1 after initialization")
    endif()
    check(${mapping}_state roughness "${WORKING_DIRECTORY}/${mapping}/chk00004" 1 ${check_option})
    # Change the tagging box on restart so an existing fine level is remade.
    run_leg(${mapping}_remade "amr.max_level=1" "erf.most.roughness_file_name=${files}"
        "erf.restart=${WORKING_DIRECTORY}/${mapping}/chk00003"
        "erf.wall.in_box_hi=28 24 16" "max_step=5")
    check(${mapping}_remade_state roughness_remade "${WORKING_DIRECTORY}/${mapping}/chk00003"
        "${WORKING_DIRECTORY}/${mapping}_remade/chk00005" ${check_option})
endforeach()

message(STATUS "EB SurfaceLayer activation, varying roughness, drag, checkpoint, and AMR checks passed")
