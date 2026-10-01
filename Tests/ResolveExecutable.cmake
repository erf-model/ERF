#=============================================================================
# Shared executable-path resolution for the checker-driven CTest scripts.
#
# On Windows the build is driven by a multi-config generator, so the config
# subdirectory a binary lands in ("Debug", "Release", ...) is not known when the
# tests are registered.  CTestList.cmake and the top-level CMakeLists.txt
# therefore hand out wildcard paths:
#
#     <build>/Exec/<dir>/*/erf_exec.exe
#     <build>/Submodules/AMReX/Tools/Plotfile/*/amrex_fcompare.exe
#
# The tests that build a shell string run them through `sh -c`, which expands
# the wildcard.  The tests driven through `cmake -P` build an argv list for
# execute_process instead, and execute_process never invokes a shell: the
# wildcard reaches the operating system as a literal file name and the test
# fails before the program starts, with "no such file or directory".
#
# erf_resolve_executable() does the expansion those scripts would otherwise
# have to do without, so a driver script can take the path it was given and use
# it whatever platform produced it.
#=============================================================================

include_guard(GLOBAL)

#-----------------------------------------------------------------------------
# erf_resolve_executable(<out_var> <path>
#                        [CONFIG  <config>]  # build config to prefer, e.g. $<CONFIG>
#                        [CONTEXT <text>])   # prefix for error messages
#
# Sets <out_var> to a program path that execute_process can run.
#
# A path with no wildcard in it is passed through untouched, so the resolution
# is a no-op on the single-config generators used everywhere but Windows, and a
# path that does not exist still produces the same error it always did, from the
# run that needed it.  A path with a wildcard is expanded against the file
# system, and not finding a program is a fatal error here, where the pattern and
# the reason can be named, rather than an exec failure inside a test log.
#
# CONFIG picks between the matches when more than one config has been built: a
# match in a directory of that name wins.  Without it, or with no match under
# it, the first match in lexicographic order is used and the choice is reported,
# since running a binary from a config the user did not ask for is the kind of
# thing that must not happen quietly.
#-----------------------------------------------------------------------------
function(erf_resolve_executable out_var path)
    cmake_parse_arguments(_ARG "" "CONFIG;CONTEXT" "" ${ARGN})

    if(_ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "erf_resolve_executable: unexpected arguments: ${_ARG_UNPARSED_ARGUMENTS}")
    endif()

    set(_context "executable")
    if(NOT "${_ARG_CONTEXT}" STREQUAL "")
        set(_context "${_ARG_CONTEXT}")
    endif()

    if(NOT "${path}" MATCHES "[*?]")
        set(${out_var} "${path}" PARENT_SCOPE)
        return()
    endif()

    file(GLOB _matches "${path}")
    set(_programs "")
    foreach(_match IN LISTS _matches)
        if(NOT IS_DIRECTORY "${_match}")
            list(APPEND _programs "${_match}")
        endif()
    endforeach()

    list(LENGTH _programs _count)
    if(_count EQUAL 0)
        message(FATAL_ERROR
            "${_context}: no program matches ${path}\n"
            "  The pattern is how the multi-config generators name a binary whose\n"
            "  config subdirectory is not known until build time, so an empty match\n"
            "  usually means that config was never built.")
    endif()

    list(SORT _programs)
    set(_chosen "")
    if(NOT "${_ARG_CONFIG}" STREQUAL "")
        foreach(_program IN LISTS _programs)
            get_filename_component(_config_dir "${_program}" DIRECTORY)
            get_filename_component(_config_dir "${_config_dir}" NAME)
            if("${_config_dir}" STREQUAL "${_ARG_CONFIG}")
                set(_chosen "${_program}")
                break()
            endif()
        endforeach()
    endif()

    if("${_chosen}" STREQUAL "")
        list(GET _programs 0 _chosen)
        if(_count GREATER 1)
            string(JOIN "\n    " _all ${_programs})
            message(STATUS
                "${_context}: ${path} matches ${_count} programs, using ${_chosen}\n"
                "    ${_all}")
        endif()
    endif()

    set(${out_var} "${_chosen}" PARENT_SCOPE)
endfunction()
