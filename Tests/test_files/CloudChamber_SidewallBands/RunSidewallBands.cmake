cmake_minimum_required(VERSION 3.20)
include("${MPI_HELPER}")
if(NOT MPIEXEC OR NOT MPIEXEC_NUMPROC_FLAG OR NRANKS LESS 2)
    message(FATAL_ERROR "CloudChamber_SidewallBands requires an MPI launcher and at least two ranks")
endif()
erf_mpi_launcher_command(launch LAUNCHER "${MPIEXEC}" NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS} PREFLAGS "${MPIEXEC_PREFLAGS}" CONTEXT "CloudChamber_SidewallBands")
math(EXPR NZ "8 * ${NRANKS}")

function(run_wall_case leg XLO YHI)
    set(dir "${CMAKE_CURRENT_BINARY_DIR}/${leg}")
    file(MAKE_DIRECTORY "${dir}")
    # ERF appends budget intervals. Never read a previous run's step 1.
    file(REMOVE "${dir}/cloud_chamber_budget.dat")
    configure_file("${CMAKE_CURRENT_LIST_DIR}/CloudChamber_SidewallBands.i.in" "${dir}/inputs" @ONLY)
    execute_process(COMMAND ${launch} "${TEST_EXE}" inputs
        WORKING_DIRECTORY "${dir}" RESULT_VARIABLE rc
        OUTPUT_FILE "${dir}/simulation.log" ERROR_FILE "${dir}/simulation.log")
    if(NOT rc STREQUAL "0")
        message(FATAL_ERROR "${leg}: erf_exec failed (${rc}); see ${dir}/simulation.log")
    endif()
endfunction()

run_wall_case(cold "xlo.temperature = 280" "yhi.temperature = 280")
run_wall_case(hot  "xlo.temperature = 304" "yhi.temperature = 304")
run_wall_case(one_band "xlo.temperature_height_edges = 0 1\nxlo.temperature_values = 280"
                       "yhi.temperature = 280")
run_wall_case(panels "xlo.temperature_height_edges = 0 0.1875 0.90625 1\nxlo.temperature_values = 280 304 280"
                     "yhi.temperature_height_edges = 0 0.265625 1\nyhi.temperature_values = 304 280")

execute_process(COMMAND "${AWK}" -f "${CMAKE_CURRENT_LIST_DIR}/CheckBudgets.awk"
    cold/cloud_chamber_budget.dat hot/cloud_chamber_budget.dat
    one_band/cloud_chamber_budget.dat panels/cloud_chamber_budget.dat
    COMMAND_ERROR_IS_FATAL ANY)

# Same nonuniform sidewalls and MOST options in both runs: a wet horizontal
# plate is valid, but MOST on a banded sidewall must fail before timestepping.
foreach(face IN ITEMS zlo xlo)
    set(dir "${CMAKE_CURRENT_BINARY_DIR}/most_${face}")
    file(MAKE_DIRECTORY "${dir}")
    file(REMOVE "${dir}/cloud_chamber_budget.dat")
    execute_process(COMMAND ${launch} "${TEST_EXE}" "${CMAKE_CURRENT_BINARY_DIR}/panels/inputs"
        zlo.heat=conductive zlo.moisture=wet zlo.temperature=292 prob.U_0=0.1
        ${face}.momentum_transfer_model=bulk_aero ${face}.heat_transfer_model=bulk_aero
        ${face}.vapor_transfer_model=bulk_aero ${face}.coefficient_source=most
        ${face}.z0_m=0.001 ${face}.z0_h=0.001 ${face}.z0_q=0.001 amrex.call_addr2line=0
        WORKING_DIRECTORY "${dir}" RESULT_VARIABLE rc
        OUTPUT_FILE "${dir}/simulation.log" ERROR_FILE "${dir}/simulation.log")
    if(face STREQUAL "zlo")
        if(NOT rc STREQUAL "0")
            message(FATAL_ERROR "MOST plate with banded sidewalls failed (${rc}); see ${dir}/simulation.log")
        endif()
    else()
        file(STRINGS "${dir}/simulation.log" diagnostic
            REGEX "xlo supports coefficient_source = most only on zlo or zhi")
        if(rc STREQUAL "0" OR NOT diagnostic)
            message(FATAL_ERROR "Expected the banded xlo MOST rejection (${rc}); see ${dir}/simulation.log")
        endif()
    endif()
endforeach()
