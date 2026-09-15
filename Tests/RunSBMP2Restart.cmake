if(NOT DEFINED MPIEXEC OR NOT DEFINED MPIEXEC_NUMPROC_FLAG OR
   NOT DEFINED NRANKS OR NOT DEFINED TEST_EXE OR NOT DEFINED INPUT OR
   NOT DEFINED WORKING_DIRECTORY OR NOT DEFINED LOG)
    message(FATAL_ERROR "RunSBMP2Restart.cmake missing required argument")
endif()

set(_root "${WORKING_DIRECTORY}/restart_equivalence")
file(REMOVE_RECURSE "${_root}")
file(MAKE_DIRECTORY "${_root}")

set(_mpi_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NRANKS})
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(_mpi_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
    list(APPEND _mpi_command ${_mpi_preflags})
endif()

set(_continuous_prefix "${_root}/continuous")
set(_split_prefix "${_root}/split")
set(_restart_prefix "${_root}/restart")
set(_continuous_log "${_root}/continuous.log")
set(_split_log "${_root}/split.log")
set(_restart_log "${_root}/restart.log")

set(_continuous_command ${_mpi_command} ${TEST_EXE} ${INPUT}
    erf.check_int=1 erf.check_file=${_continuous_prefix} max_step=2)
execute_process(
    COMMAND ${_continuous_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${_continuous_log}"
    ERROR_FILE "${_continuous_log}"
    RESULT_VARIABLE _continuous_result)
if(NOT _continuous_result EQUAL 0)
    message(FATAL_ERROR "continuous SBM AMR run failed: ${_continuous_result}")
endif()

set(_split_command ${_mpi_command} ${TEST_EXE} ${INPUT}
    erf.check_int=1 erf.check_file=${_split_prefix} max_step=1)
execute_process(
    COMMAND ${_split_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${_split_log}"
    ERROR_FILE "${_split_log}"
    RESULT_VARIABLE _split_result)
if(NOT _split_result EQUAL 0)
    message(FATAL_ERROR "checkpoint-producing SBM AMR run failed: ${_split_result}")
endif()

set(_split_checkpoint "${_split_prefix}00001")
if(NOT EXISTS "${_split_checkpoint}/Header")
    message(FATAL_ERROR "checkpoint-producing run did not create ${_split_checkpoint}")
endif()

set(_restart_command ${_mpi_command} ${TEST_EXE} ${INPUT}
    amr.restart=${_split_checkpoint}
    erf.check_int=1 erf.check_file=${_restart_prefix} max_step=2)
execute_process(
    COMMAND ${_restart_command}
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${_restart_log}"
    ERROR_FILE "${_restart_log}"
    RESULT_VARIABLE _restart_result)
if(NOT _restart_result EQUAL 0)
    message(FATAL_ERROR "SBM AMR restart run failed: ${_restart_result}")
endif()

set(_continuous_checkpoint "${_continuous_prefix}00002")
set(_restart_checkpoint "${_restart_prefix}00002")
foreach(_checkpoint IN ITEMS "${_continuous_checkpoint}" "${_restart_checkpoint}")
    if(NOT EXISTS "${_checkpoint}/Header")
        message(FATAL_ERROR "missing final checkpoint ${_checkpoint}")
    endif()
endforeach()

file(READ "${_continuous_log}" _continuous_text)
file(READ "${_restart_log}" _restart_text)
foreach(_text IN ITEMS "Coarse STEP 2 ends" "SBM layout identity")
    string(FIND "${_continuous_text}" "${_text}" _continuous_offset)
    string(FIND "${_restart_text}" "${_text}" _restart_offset)
    if(_continuous_offset EQUAL -1 OR _restart_offset EQUAL -1)
        message(FATAL_ERROR "restart equivalence logs are missing expected text: ${_text}")
    endif()
endforeach()

# Compare the authoritative spectral state and compact projection at every
# level, including the per-rank FAB payloads.  Exact hashes are appropriate
# here: both trajectories use the same MPI decomposition and should execute
# the same deterministic AMR lifecycle after the checkpoint boundary.
set(_relative_files)
foreach(_pattern IN ITEMS "Level_*/SBMAux_*" "Level_*/Cell_*")
    file(GLOB _matches RELATIVE "${_continuous_checkpoint}"
         "${_continuous_checkpoint}/${_pattern}")
    list(APPEND _relative_files ${_matches})
endforeach()
list(SORT _relative_files)
if(NOT _relative_files)
    message(FATAL_ERROR "no SBM spectral or compact checkpoint payloads found")
endif()

foreach(_relative IN LISTS _relative_files)
    if(NOT EXISTS "${_restart_checkpoint}/${_relative}")
        message(FATAL_ERROR "restart checkpoint is missing ${_relative}")
    endif()
    file(SHA256 "${_continuous_checkpoint}/${_relative}" _continuous_hash)
    file(SHA256 "${_restart_checkpoint}/${_relative}" _restart_hash)
    if(NOT "${_continuous_hash}" STREQUAL "${_restart_hash}")
        message(FATAL_ERROR "continuous/restart payload mismatch: ${_relative}")
    endif()
endforeach()

file(SHA256 "${_continuous_checkpoint}/SBM_Schema" _continuous_schema_hash)
file(SHA256 "${_restart_checkpoint}/SBM_Schema" _restart_schema_hash)
if(NOT "${_continuous_schema_hash}" STREQUAL "${_restart_schema_hash}")
    message(FATAL_ERROR "continuous/restart SBM schema mismatch")
endif()

file(WRITE "${LOG}" "continuous/restart SBM AMR equivalence passed\n")
