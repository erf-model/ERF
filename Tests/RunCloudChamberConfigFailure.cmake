if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED TEST_EXE OR NOT DEFINED BASE_INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG OR
   NOT DEFINED CASE OR NOT DEFINED EXPECTED_DIAGNOSTIC)
    message(FATAL_ERROR "RunCloudChamberConfigFailure.cmake missing required argument")
endif()

if(NOT EXISTS "${MPIEXEC}")
    message(FATAL_ERROR "Cloud Chamber configuration test MPI launcher is missing: ${MPIEXEC}")
endif()
if(NOT EXISTS "${TEST_EXE}")
    message(FATAL_ERROR "Cloud Chamber configuration test executable is missing: ${TEST_EXE}")
endif()
if(NOT EXISTS "${BASE_INPUT}")
    message(FATAL_ERROR "Cloud Chamber configuration test fixture is missing: ${BASE_INPUT}")
endif()

file(READ "${BASE_INPUT}" input_text)
set(original_input "${input_text}")

if(CASE STREQUAL "dry_physical_rh")
    string(REPLACE
        "prob.thermodynamic_initialization = physical_temperature_rh"
        "prob.thermodynamic_initialization = physical_temperature_rh\nprob.initial_relative_humidity = 0.95"
        input_text "${input_text}")
elseif(CASE STREQUAL "satadj_physical_missing_rh")
    string(REPLACE "prob.initial_relative_humidity = 0.95\n" ""
        input_text "${input_text}")
elseif(CASE STREQUAL "physical_legacy_theta_perturbation")
    string(REPLACE
        "prob.temperature_perturbation_amplitude = 0.02"
        "prob.temperature_perturbation_amplitude = 0.02\nprob.theta_perturbation_amplitude = 0.02"
        input_text "${input_text}")
elseif(CASE STREQUAL "legacy_physical_temperature")
    string(REPLACE
        "prob.thermodynamic_initialization = physical_temperature_rh"
        "prob.thermodynamic_initialization = legacy_theta_qv\nprob.theta_bottom = 300.0\nprob.theta_top = 284.0\nprob.theta_perturbation_amplitude = 0.02"
        input_text "${input_text}")
else()
    message(FATAL_ERROR "Unknown Cloud Chamber configuration test case: ${CASE}")
endif()

if(input_text STREQUAL original_input)
    message(FATAL_ERROR "Cloud Chamber configuration test case did not modify its fixture: ${CASE}")
endif()

set(generated_input "${WORKING_DIRECTORY}/invalid_${CASE}.i")
file(WRITE "${generated_input}" "${input_text}")

execute_process(
    COMMAND "${MPIEXEC}" "${MPIEXEC_NUMPROC_FLAG}" 1 ${MPIEXEC_PREFLAGS}
            "${TEST_EXE}" "${generated_input}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_VARIABLE simulation_stdout
    ERROR_VARIABLE simulation_stderr
    RESULT_VARIABLE simulation_result)

set(combined_output "${simulation_stdout}\n${simulation_stderr}")
file(WRITE "${LOG}"
    "CASE=${CASE}\nRESULT=${simulation_result}\n\n${combined_output}")

if("${simulation_result}" STREQUAL "0")
    message(FATAL_ERROR
        "Cloud Chamber invalid configuration unexpectedly succeeded (${CASE})")
endif()

string(FIND "${combined_output}" "${EXPECTED_DIAGNOSTIC}" diagnostic_position)
if(diagnostic_position EQUAL -1)
    message(FATAL_ERROR
        "Cloud Chamber invalid configuration stopped without expected diagnostic '${EXPECTED_DIAGNOSTIC}' (${CASE}); see ${LOG}")
endif()
