#=============================================================================
# Shared MPI launcher handling for the checker-driven CTest scripts.
#
# MPIEXEC_EXECUTABLE is not always a single program path. Site configurations
# may point it at a command that carries a subcommand or fixed flags:
#
#     -DMPIEXEC_EXECUTABLE="flux run"
#     -DMPIEXEC_EXECUTABLE="srun --mpibind=off"
#
# The tests that build a shell string (MPI_COMMANDS in CTestList.cmake) cope
# with those values because `sh -c` re-splits the string. The tests driven
# through `cmake -P` build an argv list for execute_process instead, and
# execute_process never invokes a shell: a multi-word value used as argv[0] is
# looked up as one literal file name, so "flux run" fails.
#
# erf_mpi_launcher_command() converts such a value into a validated argv prefix,
# so the driver scripts never touch MPIEXEC directly.
# =============================================================================

include_guard(GLOBAL)

#-----------------------------------------------------------------------------
# erf_mpi_launcher_tokens(<out_var> <launcher>)
#
# Split <launcher> into argv tokens.
#
# A value that names an existing file is kept as a single token, so launcher
# paths that themselves contain spaces keep working (notably on Windows,
# e.g. "C:/Program Files/Microsoft MPI/Bin/mpiexec.exe"). Anything else is
# tokenized with the platform's native quoting rules, which also allows a quoted
# path followed by a subcommand: "\"/opt/my mpi/flux\" run".
# -----------------------------------------------------------------------------
function(erf_mpi_launcher_tokens out_var launcher)
    set(_tokens "")
    if(NOT "${launcher}" STREQUAL "")
        if(EXISTS "${launcher}" AND NOT IS_DIRECTORY "${launcher}")
            set(_tokens "${launcher}")
        else()
            separate_arguments(_tokens NATIVE_COMMAND "${launcher}")
        endif()
    endif()
    set(${out_var} "${_tokens}" PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------------
# erf_mpi_launcher_command(<out_var>
#                          LAUNCHER      <value>   # MPIEXEC / MPIEXEC_EXECUTABLE
#                          [NUMPROC_FLAG <flag>]   # MPIEXEC_NUMPROC_FLAG
#                          [NRANKS       <n>]      # defaults to 1
#                          [PREFLAGS     <string>] # MPIEXEC_PREFLAGS
#                          [CONTEXT      <text>]   # prefix for error messages
#                          [OPTIONAL])             # do not fail on a missing
#                                                  # launcher program
#
# Sets <out_var> to the argv prefix that runs a command under the launcher:
#
#     <launcher tokens...> [<numproc flag> <nranks>] [<preflags...>]
#
# and to an empty list when LAUNCHER is empty (non-MPI builds), so callers can
# unconditionally append the executable and its arguments.
#
# Only the launcher *program* (the first token) is validated; a subcommand
# such as `run` is never mistaken for a file. A bare program name is resolved
# against PATH so that failures are reported here, with the value of MPIEXEC
# that produced them, rather than as a bare exec error inside a test log.
#-----------------------------------------------------------------------------
function(erf_mpi_launcher_command out_var)
    cmake_parse_arguments(_ARG
        "OPTIONAL"
        "LAUNCHER;NUMPROC_FLAG;NRANKS;PREFLAGS;CONTEXT"
        ""
        ${ARGN})

    if(_ARG_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "erf_mpi_launcher_command: unexpected arguments: ${_ARG_UNPARSED_ARGUMENTS}")
    endif()

    set(_context "MPI launcher")
    if(NOT "${_ARG_CONTEXT}" STREQUAL "")
        set(_context "${_ARG_CONTEXT}")
    endif()

    set(_command "")

    if(NOT "${_ARG_LAUNCHER}" STREQUAL "")
        erf_mpi_launcher_tokens(_tokens "${_ARG_LAUNCHER}")
        list(LENGTH _tokens _token_count)
        if(_token_count EQUAL 0)
            message(FATAL_ERROR
                "${_context}: MPI launcher has no program name: \"${_ARG_LAUNCHER}\"")
        endif()

        list(GET _tokens 0 _program)

        # A program given as a path must exist; a bare name is looked up on
        # PATH. Either way the subcommand and flags are left alone.
        string(FIND "${_program}" "/" _slash_index)
        string(FIND "${_program}" "\\" _backslash_index)
        if(IS_ABSOLUTE "${_program}" OR NOT _slash_index EQUAL -1
                OR NOT _backslash_index EQUAL -1)
            if(NOT EXISTS "${_program}" AND NOT _ARG_OPTIONAL)
                message(FATAL_ERROR
                    "${_context}: MPI launcher program is missing: ${_program}\n"
                    "  (MPIEXEC = \"${_ARG_LAUNCHER}\")")
            endif()
        else()
            find_program(_erf_mpi_launcher_program NAMES "${_program}")
            if(_erf_mpi_launcher_program)
                list(REMOVE_AT _tokens 0)
                list(INSERT _tokens 0 "${_erf_mpi_launcher_program}")
            elseif(NOT _ARG_OPTIONAL)
                message(FATAL_ERROR
                    "${_context}: MPI launcher program was not found on PATH: ${_program}\n"
                    "  (MPIEXEC = \"${_ARG_LAUNCHER}\")")
            endif()
            unset(_erf_mpi_launcher_program CACHE)
        endif()

        list(APPEND _command ${_tokens})

        if(NOT "${_ARG_NUMPROC_FLAG}" STREQUAL "")
            set(_nranks 1)
            if(NOT "${_ARG_NRANKS}" STREQUAL "")
                set(_nranks "${_ARG_NRANKS}")
            endif()
            list(APPEND _command "${_ARG_NUMPROC_FLAG}" "${_nranks}")
        endif()

        if(NOT "${_ARG_PREFLAGS}" STREQUAL "")
            # The launcher above is split with NATIVE_COMMAND, these flags with
            # UNIX_COMMAND. The two modes behave identically except on Windows,
            # where NATIVE_COMMAND reads a backslash as a path separator rather
            # than as an escape character. A launcher path needs that; a flag
            # string like "--bind-to core" does not. UNIX_COMMAND is also what
            # every driver script used for preflags before this helper existed,
            # so existing MPIEXEC_PREFLAGS values tokenize exactly as before.
            separate_arguments(_preflags UNIX_COMMAND "${_ARG_PREFLAGS}")
            list(APPEND _command ${_preflags})
        endif()
    endif()

    set(${out_var} "${_command}" PARENT_SCOPE)
endfunction()
