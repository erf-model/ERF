# Self-test of Tests/CompareDataLogs.cmake.
#
# A comparison that never fails passes every test it is used in, so the comparator needs its
# own test: it must accept the last-digit disagreement it was written for and still reject a
# real difference, the surface-history factor of the number of stacked boxes above all.
#
# Pure CMake, no ERF run.  Run with
#   cmake -DWORK_DIR=<scratch directory> -P Tests/CompareDataLogsSelfTest.cmake

if("${WORK_DIR}" STREQUAL "")
    message(FATAL_ERROR "CompareDataLogsSelfTest.cmake: WORK_DIR must be given and non-empty")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/CompareDataLogs.cmake")

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

set(selftest_failures 0)
set(selftest_cases 0)

# The data logs are printed with datprecision = 6 significant digits (Source/ERF.H) and
# RunBoxParity.cmake tolerates two units of the last of them
set(SIGDIGITS 6)
set(ULPS 2)

# Compare two logs given as strings and check the verdict against the one expected
function(expect_logs name expected_agree text_a text_b)
    math(EXPR _cases "${selftest_cases} + 1")
    set(selftest_cases ${_cases} PARENT_SCOPE)
    set(_file_a "${WORK_DIR}/${name}_a.dat")
    set(_file_b "${WORK_DIR}/${name}_b.dat")
    file(WRITE "${_file_a}" "${text_a}")
    file(WRITE "${_file_b}" "${text_b}")
    erf_compare_data_logs("${_file_a}" "${_file_b}" ${SIGDIGITS} ${ULPS} _agree _message)
    if(_agree AND NOT expected_agree)
        message(SEND_ERROR "${name}: the logs were accepted but must differ")
        math(EXPR _failures "${selftest_failures} + 1")
        set(selftest_failures ${_failures} PARENT_SCOPE)
    elseif(NOT _agree AND expected_agree)
        message(SEND_ERROR "${name}: the logs were rejected but must agree: ${_message}")
        math(EXPR _failures "${selftest_failures} + 1")
        set(selftest_failures ${_failures} PARENT_SCOPE)
    endif()
endfunction()

set(header "          time        u_star        t_star          olen\n")

# Identical logs agree
expect_logs(identical TRUE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0      0.890592             0         1e+150\n")

# The last printed digit may differ: the two runs decompose the domain differently, so the
# planar sum reduces in a different order and the fields themselves agree only to the
# plotfile tolerance.  This is what the comparison is numeric for.
expect_logs(last_digit TRUE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0      0.890591             0         1e+150\n")

# Two units of the last digit are accepted, three are not
expect_logs(two_ulps TRUE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0      0.890594             0         1e+150\n")
expect_logs(three_ulps FALSE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0      0.890595             0         1e+150\n")

# The bug this comparison was added for: on grids split in z every surface cell was counted
# once per stacked box, so u* came out multiplied by the number of boxes in the column
expect_logs(zsplit_factor FALSE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0       2.67178             0         1e+150\n")

# A sum that cancels to roundoff in one run and to exactly zero in the other is not a
# difference; a value that is small but real is
expect_logs(roundoff_zero TRUE
    "${header}             0             0             0         1e+150\n"
    "${header}             0   1.23457e-16             0         1e+150\n")
expect_logs(small_but_real FALSE
    "${header}             0             0             0         1e+150\n"
    "${header}             0   1.23457e-06             0         1e+150\n")

# The same number in fixed and scientific notation is the same number
expect_logs(mixed_notation TRUE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0   8.90592e-01             0         1e+150\n")

# Values of opposite sign differ however small they are, as long as they are not roundoff
expect_logs(opposite_sign FALSE
    "${header}             0      0.890592     -0.001234         1e+150\n"
    "${header}             0      0.890592      0.001234         1e+150\n")

# Orders of magnitude apart, in either order
expect_logs(orders_apart FALSE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0   8.90592e+09             0         1e+150\n")

# Text fields, the column headers above all, are compared as text
expect_logs(header_renamed FALSE
    "${header}             0      0.890592             0         1e+150\n"
    "          time         u_ast        t_star          olen\n             0      0.890592             0         1e+150\n")

# A run that stopped early, or wrote a different set of columns, is not a match
expect_logs(fewer_lines FALSE
    "${header}             0      0.890592             0         1e+150\n             1      0.890593             0         1e+150\n"
    "${header}             0      0.890592             0         1e+150\n")
expect_logs(fewer_fields FALSE
    "${header}             0      0.890592             0         1e+150\n"
    "${header}             0      0.890592             0\n")

if(selftest_failures GREATER 0)
    message(FATAL_ERROR "CompareDataLogsSelfTest: ${selftest_failures} of ${selftest_cases} cases failed")
endif()
message(STATUS "CompareDataLogsSelfTest: ${selftest_cases} cases passed")
