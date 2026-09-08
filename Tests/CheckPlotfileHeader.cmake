cmake_minimum_required(VERSION 3.20)

foreach(required_path HEADER EXPECTED_NAMES_FILE LOG EXPECTED_UNAVAILABLE_FILE)
    if(NOT DEFINED ${required_path} OR "${${required_path}}" STREQUAL "")
        message(FATAL_ERROR "${required_path} is required")
    endif()
    if(NOT EXISTS "${${required_path}}")
        message(FATAL_ERROR "Required plotfile test file does not exist: ${${required_path}}")
    endif()
endforeach()

file(STRINGS "${HEADER}" header_lines)
list(LENGTH header_lines header_line_count)
if(header_line_count LESS 2)
    message(FATAL_ERROR "AMReX Header is too short: ${HEADER}")
endif()

list(GET header_lines 1 header_variable_count)
if(NOT header_variable_count MATCHES "^[0-9]+$")
    message(FATAL_ERROR "AMReX Header variable count is not an integer: ${header_variable_count}")
endif()

set(actual_names)
if(header_variable_count GREATER 0)
    math(EXPR last_header_name_index "1 + ${header_variable_count}")
    if(header_line_count LESS_EQUAL last_header_name_index)
        message(FATAL_ERROR "AMReX Header does not contain ${header_variable_count} variable names")
    endif()
    foreach(index RANGE 2 ${last_header_name_index})
        list(GET header_lines ${index} name)
        list(APPEND actual_names "${name}")
    endforeach()
endif()

function(read_name_file path output_variable)
    file(STRINGS "${path}" raw_lines)
    set(names)
    foreach(line IN LISTS raw_lines)
        string(STRIP "${line}" line)
        if(NOT line STREQUAL "" AND NOT line MATCHES "^#")
            list(APPEND names "${line}")
        endif()
    endforeach()
    set(${output_variable} "${names}" PARENT_SCOPE)
endfunction()

read_name_file("${EXPECTED_NAMES_FILE}" expected_names)
read_name_file("${EXPECTED_UNAVAILABLE_FILE}" expected_unavailable)

list(LENGTH actual_names actual_count)
list(LENGTH expected_names expected_count)
if(NOT header_variable_count EQUAL actual_count OR
   NOT header_variable_count EQUAL expected_count)
    message(FATAL_ERROR
        "Header count mismatch: header=${header_variable_count}, "
        "actual=${actual_count}, expected=${expected_count}")
endif()

foreach(index RANGE 0 ${expected_count})
    if(index LESS expected_count)
        list(GET actual_names ${index} actual_name)
        list(GET expected_names ${index} expected_name)
        if(NOT actual_name STREQUAL expected_name)
            message(FATAL_ERROR
                "Header variable ${index} mismatch: actual='${actual_name}', "
                "expected='${expected_name}'")
        endif()
    endif()
endforeach()

file(READ "${LOG}" log_contents)
set(confirmed_warnings)
foreach(name IN LISTS expected_unavailable)
    string(REGEX MATCHALL
        "Requested to plot variable '${name}' but it is not available"
        name_matches "${log_contents}")
    list(LENGTH name_matches name_match_count)
    if(NOT name_match_count EQUAL 1)
        message(FATAL_ERROR
            "Expected exactly one unavailable warning for '${name}', "
            "found ${name_match_count}")
    endif()
    list(APPEND confirmed_warnings "${name}")
endforeach()

string(REGEX MATCHALL
    "Requested to plot variable '[^']+' but it is not available"
    all_warning_matches "${log_contents}")
list(LENGTH all_warning_matches warning_count)
list(LENGTH expected_unavailable expected_warning_count)
if(NOT warning_count EQUAL expected_warning_count)
    message(FATAL_ERROR
        "Unavailable warning count mismatch: actual=${warning_count}, "
        "expected=${expected_warning_count}")
endif()

message(STATUS "Actual Header names: ${actual_names}")
message(STATUS "Confirmed unavailable warnings: ${confirmed_warnings}")
