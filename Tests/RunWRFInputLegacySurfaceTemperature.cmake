cmake_minimum_required(VERSION 3.20)

# Restart from the checkpoints older versions wrote, holding the SST and TSK as potential
# temperature, and check how they are migrated (see legacy_restarts.py).

if(DEFINED ENV{CRAY_MPICH_DIR})
  set(ENV{MPICH_GPU_SUPPORT_ENABLED} 0)
endif()

execute_process(COMMAND "${PYTHON}" "${TEST_DIR}/make_wrfinputs.py"
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/generate.log"
  ERROR_FILE "${RUN_DIR}/generate.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "WRF input generation failed (${status}); inspect ${RUN_DIR}/generate.err")
endif()

file(REMOVE_RECURSE "${RUN_DIR}/chk00000")
execute_process(COMMAND "${ERF_EXECUTABLE}" "${RUN_DIR}/WRFInput_LegacySurfaceTemperature.i"
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/run.log"
  ERROR_FILE "${RUN_DIR}/run.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "WRF input run failed (${status}); inspect ${RUN_DIR}/run.log and run.err")
endif()
if(NOT EXISTS "${RUN_DIR}/chk00000/Level_0/SST_0_H" OR
   NOT EXISTS "${RUN_DIR}/chk00000/Level_0/TSK_0_H")
  message(FATAL_ERROR "The initial run did not write the SST and TSK to chk00000")
endif()

# The same from wrfinput_mixed, part land, whose SST over land is a placeholder
file(REMOVE_RECURSE "${RUN_DIR}/mixed_chk00000")
execute_process(COMMAND "${ERF_EXECUTABLE}" "${RUN_DIR}/WRFInput_LegacySurfaceTemperature.i"
  "erf.nc_init_file_0=wrfinput_mixed" "erf.check_file=mixed_chk"
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/run_mixed.log"
  ERROR_FILE "${RUN_DIR}/run_mixed.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Mixed land/water WRF input run failed (${status}); inspect "
                      "${RUN_DIR}/run_mixed.log and run_mixed.err")
endif()
if(NOT EXISTS "${RUN_DIR}/mixed_chk00000/Level_0/SST_0_H" OR
   NOT EXISTS "${RUN_DIR}/mixed_chk00000/Level_0/TSK_0_H" OR
   NOT EXISTS "${RUN_DIR}/mixed_chk00000/Level_0/LMASK_0_H")
  message(FATAL_ERROR "The mixed run did not write the SST, TSK and LMASK to mixed_chk00000")
endif()

# And from nest_files.i, whose level 1 lists two init files
file(REMOVE_RECURSE "${RUN_DIR}/nest_chk00000")
execute_process(COMMAND "${ERF_EXECUTABLE}" "${RUN_DIR}/nest_files.i"
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/run_nest.log"
  ERROR_FILE "${RUN_DIR}/run_nest.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Nested WRF input run failed (${status}); inspect "
                      "${RUN_DIR}/run_nest.log and run_nest.err")
endif()
if(NOT EXISTS "${RUN_DIR}/nest_chk00000/Level_1/SST_0_H" OR
   NOT EXISTS "${RUN_DIR}/nest_chk00000/Level_1/TSK_0_H")
  message(FATAL_ERROR "The nested run did not write the level-1 SST and TSK to nest_chk00000")
endif()

execute_process(COMMAND "${PYTHON}" "${TEST_DIR}/legacy_restarts.py"
  "${ERF_EXECUTABLE}" "${RUN_DIR}"
  WORKING_DIRECTORY "${RUN_DIR}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Legacy SST/TSK checkpoint migration check failed (${status}); "
                      "inspect the <case>.log files in ${RUN_DIR}")
endif()
