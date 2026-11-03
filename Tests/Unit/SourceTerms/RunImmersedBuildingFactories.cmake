foreach(name ERF_EXECUTABLE INPUT_FILE EXPECTED_FACTORY TEST_ROOT)
  if(NOT DEFINED ${name} OR "${${name}}" STREQUAL "")
    message(FATAL_ERROR "${name} is required")
  endif()
endforeach()

# Pass only if the run completes and reports the expected factory choice.
string(MAKE_C_IDENTIFIER "${EXPECTED_FACTORY}" run_name)
set(run_dir "${TEST_ROOT}/immersed_building_factories_${run_name}")
file(REMOVE_RECURSE "${run_dir}")
file(MAKE_DIRECTORY "${run_dir}")
separate_arguments(overrides UNIX_COMMAND "${FACTORY_OVERRIDES}")
execute_process(
  COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${overrides}
  WORKING_DIRECTORY "${run_dir}"
  RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "run failed (${result}):\n${stdout}\n${stderr}")
endif()
if(NOT stdout MATCHES "Building EB factories: ${EXPECTED_FACTORY}")
  message(FATAL_ERROR "run did not use ${EXPECTED_FACTORY} factories:\n${stdout}")
endif()
