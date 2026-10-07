if(NOT DEFINED ERF_EXECUTABLE OR NOT EXISTS "${ERF_EXECUTABLE}")
  message(FATAL_ERROR "ERF_EXECUTABLE must name a built erf_exec")
endif()
if(NOT DEFINED INPUT_FILE OR NOT EXISTS "${INPUT_FILE}")
  message(FATAL_ERROR "INPUT_FILE must name an SBM M3 periodic input deck")
endif()
if(NOT DEFINED TEST_ROOT)
  message(FATAL_ERROR "TEST_ROOT must name an isolated test output directory")
endif()
if(NOT DEFINED EXPECTED_SPECTRUM_VALUES)
  set(EXPECTED_SPECTRUM_VALUES
    "9.99999999999999955e-07,1.99999999999999991e-06,3.00000000000000008e-06,3.99999999999999982e-06,")
endif()
if(NOT DEFINED EXPECTED_CORE_VALUES)
  set(EXPECTED_CORE_VALUES
    "1.00000000000000000e+00,2.99999999999999886e+02,0.00000000000000000e+00,0.00000000000000000e+00,0.00000000000000000e+00,3.00000000000000008e-06,6.99999999999999989e-06,")
endif()

if(NUMERIC_HEADER_COMPARE)
  # HEADER_COMPARE_ULPS means units of the selected significant decimal digit
  # used by erf_numbers_close; it is not a count of binary IEEE ULPs.
  if(NOT DEFINED HEADER_COMPARE_SIGDIGITS)
    set(HEADER_COMPARE_SIGDIGITS 12)
  endif()
  if(NOT DEFINED HEADER_COMPARE_ULPS)
    set(HEADER_COMPARE_ULPS 8)
  endif()
  if(NOT HEADER_COMPARE_SIGDIGITS MATCHES "^[1-9][0-9]*$" OR
     NOT HEADER_COMPARE_ULPS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR
      "Numeric header comparison requires positive integer significant digits and tolerance units")
  endif()
  include("${CMAKE_CURRENT_LIST_DIR}/../../CompareDataLogs.cmake")
endif()

string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef _run_id)
set(_run_dir "${TEST_ROOT}/sbm_m3_${_run_id}")
file(MAKE_DIRECTORY "${_run_dir}")

function(run_erf description)
  execute_process(
    COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${ARGN}
    WORKING_DIRECTORY "${_run_dir}"
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _stdout
    ERROR_VARIABLE _stderr)
  if(NOT "${_result}" STREQUAL "0")
    message(FATAL_ERROR
      "${description} failed (exit ${_result}).\nstdout:\n${_stdout}\nstderr:\n${_stderr}")
  endif()
endfunction()

function(assert_header_vector header_file expected_csv field_name vector_kind)
  if(NOT EXISTS "${header_file}")
    message(FATAL_ERROR "${field_name} header is missing: ${header_file}")
  endif()

  file(READ "${header_file}" _header)
  if(NOT NUMERIC_HEADER_COMPARE)
    string(FIND "${_header}" "${expected_csv}" _found)
    if(_found EQUAL -1)
      message(FATAL_ERROR
        "${field_name} ${vector_kind} values do not match exactly: ${header_file}")
    endif()
    return()
  endif()

  set(_expected_csv "${expected_csv}")
  string(STRIP "${_expected_csv}" _expected_csv)
  if("${_expected_csv}" MATCHES "^," OR "${_expected_csv}" MATCHES ",," OR
     "${_expected_csv}" STREQUAL "")
    message(FATAL_ERROR
      "${field_name} ${vector_kind} expected vector is malformed: '${expected_csv}'")
  endif()
  string(REGEX REPLACE ",$" "" _expected_csv "${_expected_csv}")
  string(REPLACE "," ";" _expected_values "${_expected_csv}")
  list(LENGTH _expected_values _component_count)
  if(_component_count LESS 1)
    message(FATAL_ERROR
      "${field_name} ${vector_kind} expected vector has no components")
  endif()
  foreach(_expected_value IN LISTS _expected_values)
    string(STRIP "${_expected_value}" _expected_value)
    erf_read_decimal("${_expected_value}" _expected_ok _expected_sign
      _expected_digits _expected_exponent)
    if(NOT _expected_ok)
      message(FATAL_ERROR
        "${field_name} ${vector_kind} expected component is not numeric: '${_expected_value}'")
    endif()
  endforeach()

  if("${vector_kind}" STREQUAL "minimum")
    set(_wanted_section 1)
  elseif("${vector_kind}" STREQUAL "maximum")
    set(_wanted_section 2)
  else()
    message(FATAL_ERROR "Unknown numeric header vector kind '${vector_kind}'")
  endif()

  file(STRINGS "${header_file}" _header_lines)
  list(LENGTH _header_lines _line_count)
  set(_sections_seen 0)
  set(_found_vector_section FALSE)
  set(_first_mismatch "")
  if(_line_count GREATER 0)
    math(EXPR _last_line "${_line_count} - 1")
    foreach(_line_index RANGE 0 ${_last_line})
      list(GET _header_lines ${_line_index} _line)
      if(NOT "${_line}" MATCHES "^([0-9]+),([0-9]+)$")
        continue()
      endif()
      set(_row_count "${CMAKE_MATCH_1}")
      set(_column_count "${CMAKE_MATCH_2}")
      if(NOT _column_count EQUAL _component_count)
        continue()
      endif()
      math(EXPR _sections_seen "${_sections_seen} + 1")
      if(NOT _sections_seen EQUAL _wanted_section)
        continue()
      endif()

      set(_found_vector_section TRUE)
      if(_row_count LESS 1 OR _line_index GREATER_EQUAL _last_line)
        message(FATAL_ERROR
          "${field_name} ${vector_kind} block is malformed in ${header_file}")
      endif()
      math(EXPR _first_row "${_line_index} + 1")
      math(EXPR _last_row "${_first_row} + ${_row_count} - 1")
      if(_last_row GREATER _last_line)
        message(FATAL_ERROR
          "${field_name} ${vector_kind} block is truncated in ${header_file}")
      endif()

      foreach(_row_index RANGE ${_first_row} ${_last_row})
        list(GET _header_lines ${_row_index} _actual_row)
        string(STRIP "${_actual_row}" _actual_row)
        string(REGEX REPLACE ",$" "" _actual_row "${_actual_row}")
        if("${_actual_row}" STREQUAL "" OR "${_actual_row}" MATCHES "(^,|,,|,$)")
          message(FATAL_ERROR
            "${field_name} ${vector_kind} component vector is malformed in ${header_file}: '${_actual_row}'")
        endif()
        string(REPLACE "," ";" _actual_values "${_actual_row}")
        list(LENGTH _actual_values _actual_count)
        if(NOT _actual_count EQUAL _component_count)
          message(FATAL_ERROR
            "${field_name} ${vector_kind} vector in ${header_file} has ${_actual_count} components; expected ${_component_count}")
        endif()

        set(_row_matches TRUE)
        math(EXPR _last_component "${_component_count} - 1")
        foreach(_component_index RANGE 0 ${_last_component})
          list(GET _expected_values ${_component_index} _expected_value)
          list(GET _actual_values ${_component_index} _actual_value)
          string(STRIP "${_expected_value}" _expected_value)
          string(STRIP "${_actual_value}" _actual_value)
          erf_read_decimal("${_actual_value}" _actual_ok _actual_sign
            _actual_digits _actual_exponent)
          if(NOT _actual_ok)
            message(FATAL_ERROR
              "${field_name} ${vector_kind} component ${_component_index} is not numeric in ${header_file}: '${_actual_value}'")
          endif()
          erf_numbers_close("${_actual_value}" "${_expected_value}"
            ${HEADER_COMPARE_SIGDIGITS} ${HEADER_COMPARE_ULPS} _close)
          if(NOT _close)
            set(_row_matches FALSE)
            if("${_first_mismatch}" STREQUAL "")
              set(_first_mismatch
                "component ${_component_index}: expected '${_expected_value}', actual '${_actual_value}'")
            endif()
          endif()
        endforeach()
        if(_row_matches)
          return()
        endif()
      endforeach()
      break()
    endforeach()
  endif()

  if(NOT _found_vector_section)
    message(FATAL_ERROR
      "${field_name} ${vector_kind} ${_component_count}-component vector block was not found in ${header_file}")
  endif()
  message(FATAL_ERROR
    "${field_name} ${vector_kind} vector did not match any ${_component_count}-component per-box row in ${header_file}; ${_first_mismatch}; tolerance is ${HEADER_COMPARE_ULPS} units of the ${HEADER_COMPARE_SIGDIGITS}th significant digit")
endfunction()

function(assert_checkpoint checkpoint)
  set(_checkpoint_dir "${_run_dir}/${checkpoint}")
  foreach(_file IN ITEMS
      "${_checkpoint_dir}/SBM_Schema"
    "${_checkpoint_dir}/Level_0/SBMSpectrum_H"
      "${_checkpoint_dir}/Level_0/Cell_H")
    if(NOT EXISTS "${_file}")
      message(FATAL_ERROR "SBM smoke checkpoint is missing ${_file}")
    endif()
  endforeach()

  file(READ "${_checkpoint_dir}/SBM_Schema" _schema)
  if(NOT _schema MATCHES "constraint-policy=linear-groups-plus-canonical-persisted-v1" OR
     NOT _schema MATCHES "transport=sbm-m3-mapped-group-fct-v1" OR
     NOT _schema MATCHES "high-order=endpoint-number-wenoz3-v1")
    message(FATAL_ERROR "Checkpoint has the wrong SBM M3 schema: ${_schema}")
  endif()

  assert_header_vector("${_checkpoint_dir}/Level_0/SBMSpectrum_H"
    "${EXPECTED_SPECTRUM_VALUES}" "SBM spectrum header" "minimum")
  if(DEFINED EXPECTED_SPECTRUM_MAX_VALUES)
    assert_header_vector("${_checkpoint_dir}/Level_0/SBMSpectrum_H"
      "${EXPECTED_SPECTRUM_MAX_VALUES}" "SBM spectrum header" "maximum")
  endif()

  assert_header_vector("${_checkpoint_dir}/Level_0/Cell_H"
    "${EXPECTED_CORE_VALUES}" "SBM core state header" "minimum")
  if(DEFINED EXPECTED_CORE_MAX_VALUES)
    assert_header_vector("${_checkpoint_dir}/Level_0/Cell_H"
      "${EXPECTED_CORE_MAX_VALUES}" "SBM core state header" "maximum")
  endif()
endfunction()

run_erf("Initial SBM M3 run")
assert_checkpoint("sbm_smoke_chk00001")

if(DEFINED EXPECT_SECOND_CHECKPOINT AND EXPECT_SECOND_CHECKPOINT)
  assert_checkpoint("sbm_smoke_chk00002")
endif()

if(DEFINED SKIP_RESTART AND SKIP_RESTART)
  return()
endif()

run_erf("Exact-schema SBM M3 restart"
  "erf.restart=${_run_dir}/sbm_smoke_chk00001"
  "erf.check_file=sbm_restart_chk"
  "erf.plot_int_1=-1"
  "max_step=2"
  "stop_time=2.e-4")
assert_checkpoint("sbm_restart_chk00002")
