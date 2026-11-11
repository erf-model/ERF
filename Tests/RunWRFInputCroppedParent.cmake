cmake_minimum_required(VERSION 3.20)

if(DEFINED ENV{CRAY_MPICH_DIR})
  set(ENV{MPICH_GPU_SUPPORT_ENABLED} 0)
endif()

execute_process(COMMAND "${PYTHON}" "${TEST_DIR}/make_wrfinputs.py"
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/generate.log"
  ERROR_FILE "${RUN_DIR}/generate.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "WRF input generation failed (${status}); inspect ${RUN_DIR}/generate.err")
endif()

execute_process(COMMAND "${ERF_EXECUTABLE}" "${RUN_DIR}/WRFInput_CroppedParent.i"
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/run.log"
  ERROR_FILE "${RUN_DIR}/run.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "WRF input run failed (${status}); inspect ${RUN_DIR}/run.log and run.err")
endif()

execute_process(COMMAND "${PYTHON}" "${TEST_DIR}/check_result.py"
  "${RUN_DIR}" "${PLOTFILE_READER_DIR}"
  WORKING_DIRECTORY "${RUN_DIR}" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Cropped-parent WRF input check failed (${status})")
endif()

if(NOT EXISTS "${RUN_DIR}/chk00000/Header")
  message(FATAL_ERROR "Initial run did not write chk00000")
endif()
execute_process(COMMAND "${ERF_EXECUTABLE}" "${RUN_DIR}/WRFInput_CroppedParent.i"
  erf.restart=chk00000 erf.check_int=-1 erf.plot_int_1=-1
  WORKING_DIRECTORY "${RUN_DIR}" OUTPUT_FILE "${RUN_DIR}/restart.log"
  ERROR_FILE "${RUN_DIR}/restart.err" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "WRF input restart failed (${status}); inspect ${RUN_DIR}/restart.log and restart.err")
endif()
file(READ "${RUN_DIR}/restart.log" restart_log)
if(NOT restart_log MATCHES "Restart from native checkpoint chk00000" OR
   NOT restart_log MATCHES "Reading  PSFC" OR
   NOT restart_log MATCHES "Reading time slice 0 from wrflowinp_d03")
  message(FATAL_ERROR "Restart did not exercise the level-2 wrflow and PSFC re-read path")
endif()
