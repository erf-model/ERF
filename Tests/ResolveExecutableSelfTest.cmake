# Self-test of Tests/ResolveExecutable.cmake.
#
# The resolution runs on Windows, where a wrong answer shows up as a regression test that
# fails to start an hour into a CI job, so it is tested here instead, on whatever platform
# happens to be running: the multi-config layout it has to cope with is just a directory
# tree, and the patterns the ERF build actually hands to erf_resolve_executable expand the
# same way everywhere.  Glob matching is not wholly portable, though -- it is case-folded
# on macOS and Windows and not on Linux -- so the patterns used below stay inside the
# subset that is, which is the `*` the build itself uses.  See the directory_ignored case.
#
# Pure CMake, no ERF run.  Run with
#   cmake -DWORK_DIR=<scratch directory> -P Tests/ResolveExecutableSelfTest.cmake

if("${WORK_DIR}" STREQUAL "")
    message(FATAL_ERROR "ResolveExecutableSelfTest.cmake: WORK_DIR must be given and non-empty")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

set(selftest_failures 0)
set(selftest_cases 0)

function(expect_resolved name expected path)
    math(EXPR _cases "${selftest_cases} + 1")
    set(selftest_cases ${_cases} PARENT_SCOPE)
    erf_resolve_executable(_resolved "${path}" ${ARGN} CONTEXT "${name}")
    if(NOT "${_resolved}" STREQUAL "${expected}")
        message(SEND_ERROR "${name}: ${path} resolved to\n  ${_resolved}\nbut must resolve to\n  ${expected}")
        math(EXPR _failures "${selftest_failures} + 1")
        set(selftest_failures ${_failures} PARENT_SCOPE)
    endif()
endfunction()

# The multi-config layout: the same program built in two configurations, and a directory
# whose name matches the pattern as well, so that a match that cannot be run is seen to be
# skipped rather than merely absent.  That directory is named so that it sorts ahead of
# every config directory: it is then what the first-in-sorted-order fallback would return
# if directories were not filtered out, so a regression in the filtering fails a case
# rather than passing unnoticed behind a config that happened to sort first anyway.
foreach(_config Debug Release)
    file(MAKE_DIRECTORY "${WORK_DIR}/Exec/${_config}")
    file(WRITE "${WORK_DIR}/Exec/${_config}/erf_exec.exe" "not really a program\n")
endforeach()
file(MAKE_DIRECTORY "${WORK_DIR}/Exec/Aborted/erf_exec.exe")
file(MAKE_DIRECTORY "${WORK_DIR}/Single")
file(WRITE "${WORK_DIR}/Single/RelWithDebInfo/amrex_fcompare.exe" "not really a program\n")

# A path without a wildcard is what every single-config generator produces; it is passed
# through untouched, existing or not, so that a missing executable still fails in the run
# that needed it and with the message it always had.
expect_resolved(plain "${WORK_DIR}/Exec/erf_exec" "${WORK_DIR}/Exec/erf_exec")

# One config built: the wildcard has one match and there is nothing to choose between
expect_resolved(single_match
    "${WORK_DIR}/Single/RelWithDebInfo/amrex_fcompare.exe"
    "${WORK_DIR}/Single/*/amrex_fcompare.exe")

# Two configs built: the one CTest was invoked with wins, either way round, rather than
# whichever sorts first
expect_resolved(config_debug
    "${WORK_DIR}/Exec/Debug/erf_exec.exe"
    "${WORK_DIR}/Exec/*/erf_exec.exe" CONFIG "Debug")
expect_resolved(config_release
    "${WORK_DIR}/Exec/Release/erf_exec.exe"
    "${WORK_DIR}/Exec/*/erf_exec.exe" CONFIG "Release")

# A single-config generator expands $<CONFIG> to the empty string when CMAKE_BUILD_TYPE is
# unset, and a config that was never built has no match; either way the choice is the first
# in sorted order, and erf_resolve_executable says which one it took.
expect_resolved(no_config
    "${WORK_DIR}/Exec/Debug/erf_exec.exe"
    "${WORK_DIR}/Exec/*/erf_exec.exe")
expect_resolved(config_not_built
    "${WORK_DIR}/Exec/Debug/erf_exec.exe"
    "${WORK_DIR}/Exec/*/erf_exec.exe" CONFIG "MinSizeRel")

# A directory that matches the pattern is not a program: Aborted/erf_exec.exe is the first
# match in sorted order and must not be the one that comes back.
#
# The pattern deliberately uses only `*`.  A [...] character class would not do: kwsys
# lower-cases each file name before matching it on macOS and Windows (KWSYS_GLOB_CASE_INDEPENDENT
# in Source/kwsys/Glob.cxx) and folds the case of literal pattern characters to suit, but
# copies the contents of a character class through untouched -- so [SD] is tested against
# "debug" and "aborted" there and matches nothing, while matching on Linux.
expect_resolved(directory_ignored
    "${WORK_DIR}/Exec/Debug/erf_exec.exe"
    "${WORK_DIR}/Exec/*/erf_exec.exe")

# Nothing built at all must fail here, with the pattern named, rather than as an exec error
# inside a test log.  The failure is fatal, so it is provoked in a child cmake.
math(EXPR selftest_cases "${selftest_cases} + 1")
file(WRITE "${WORK_DIR}/no_match.cmake"
    "include(\"${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake\")\n"
    "erf_resolve_executable(_exe \"${WORK_DIR}/Exec/*/never_built.exe\" CONTEXT \"no_match\")\n"
    "message(STATUS \"resolved to \${_exe}\")\n")
execute_process(
    COMMAND "${CMAKE_COMMAND}" -P "${WORK_DIR}/no_match.cmake"
    OUTPUT_VARIABLE _no_match_output
    ERROR_VARIABLE _no_match_output
    RESULT_VARIABLE _no_match_result)
if(_no_match_result EQUAL 0)
    message(SEND_ERROR "no_match: a pattern that matches nothing was accepted: ${_no_match_output}")
    math(EXPR selftest_failures "${selftest_failures} + 1")
elseif(NOT _no_match_output MATCHES "never_built")
    message(SEND_ERROR "no_match: the error does not name the pattern: ${_no_match_output}")
    math(EXPR selftest_failures "${selftest_failures} + 1")
endif()

if(selftest_failures GREATER 0)
    message(FATAL_ERROR "ResolveExecutableSelfTest: ${selftest_failures} of ${selftest_cases} cases failed")
endif()
message(STATUS "ResolveExecutableSelfTest: ${selftest_cases} cases passed")
