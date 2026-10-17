# Negative tests of check_terrain_poisson_divergence.py on made-up logs: a converged level
# (with --max-mean, a mean above the limit fails and one below it passes),
# passes, a residual just above the limit fails, the mean subtracted from a singular
# problem is allowed to stay (and only it: a level-0 mean is not taken for level 1), the
# largest of several subdomains' divergences is the reference, and a log without solves on
# the level, or with a solve that starts from a negligible divergence or an 'after' line
# without a 'before' line, says nothing (exit 2).
# -DX= defines X as empty, so test for a value, not for DEFINED
foreach(arg PYTHON_EXE CHECKER WORKING_DIRECTORY)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "CheckTerrainPoissonDivergenceSelfTest.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()

file(MAKE_DIRECTORY "${WORKING_DIRECTORY}")
set(b0 "Max/L2 norm of divergence before solve in subdomain 0 at level 0 : 0.5 9.0 and volume-weighted sum 1.0")
set(a0 "Max/L2 norm of divergence after  solve at level 0 : 5.0e-09 1.0e-07 and volume-weighted sum 0.0")
set(b1 "Max/L2 norm of divergence before solve in subdomain 0 at level 1 : 0.4 19.0 and volume-weighted sum 1.0")

# name | log lines (separated by @) | level | expected exit status
foreach(case
        "converged|${b0}@${a0}@${b1}@Max/L2 norm of divergence after  solve at level 1 : 3.9e-06 1.0e-05 and volume-weighted sum 0.0|1|0"
        "just_above|${b1}@Max/L2 norm of divergence after  solve at level 1 : 4.1e-06 1.0e-05 and volume-weighted sum 0.0|1|1"
        "two_subdomains|Max/L2 norm of divergence before solve in subdomain 0 at level 1 : 0.4 1.0 and volume-weighted sum 1.0@Max/L2 norm of divergence before solve in subdomain 1 at level 1 : 0.1 1.0 and volume-weighted sum 1.0@Max/L2 norm of divergence after  solve at level 1 : 3.9e-06 1.0e-05 and volume-weighted sum 0.0|1|0"
        "old_operator|${b1}@Max/L2 norm of divergence after  solve at level 1 : 0.09 15.0 and volume-weighted sum -7.0e4|1|1"
        "mean_left_behind|${b0}@ Subtracting 3.0e-02 from rhs in subdomain 0@${a0}@${b1}@ Subtracting -4.0e-07 from rhs in subdomain 0@Max/L2 norm of divergence after  solve at level 1 : 4.01e-07 3.9e-05 and volume-weighted sum -223.9|1|0"
        "mean_and_more|${b1}@ Subtracting -4.0e-07 from rhs in subdomain 0@Max/L2 norm of divergence after  solve at level 1 : 4.5e-06 3.9e-05 and volume-weighted sum -223.9|1|1"
        "no_level1|${b0}@${a0}|1|2"
        "negligible_before|Max/L2 norm of divergence before solve in subdomain 0 at level 1 : 1.0e-9 1.0e-8 and volume-weighted sum 0.0@Max/L2 norm of divergence after  solve at level 1 : 1.0e-15 1.0e-14 and volume-weighted sum 0.0|1|2"
        "after_without_before|Max/L2 norm of divergence after  solve at level 1 : 1.0e-9 1.0e-8 and volume-weighted sum 0.0|1|2"
        "mean_within_limit|${b1}@ Subtracting -4.0e-10 from rhs in subdomain 0@Max/L2 norm of divergence after  solve at level 1 : 4.01e-10 3.9e-08 and volume-weighted sum -0.2|1|0|--max-mean 1.0e-9"
        "mean_above_limit|${b1}@ Subtracting -4.0e-07 from rhs in subdomain 0@Max/L2 norm of divergence after  solve at level 1 : 4.01e-07 3.9e-05 and volume-weighted sum -223.9|1|1|--max-mean 1.0e-9")
    string(REPLACE "|" ";" fields "${case}")
    list(GET fields 0 name)
    list(GET fields 1 lines)
    list(GET fields 2 level)
    list(GET fields 3 expected)
    set(extra "")
    list(LENGTH fields nfields)
    if(nfields GREATER 4)
        list(GET fields 4 extra_text)
        separate_arguments(extra UNIX_COMMAND "${extra_text}")
    endif()
    string(REPLACE "@" "\n" text "${lines}")
    file(WRITE "${WORKING_DIRECTORY}/${name}.log" "${text}\n")
    execute_process(
        COMMAND "${PYTHON_EXE}" "${CHECKER}" "${WORKING_DIRECTORY}/${name}.log" --level ${level} --max-ratio 1.0e-5 ${extra}
        OUTPUT_VARIABLE out ERROR_VARIABLE out RESULT_VARIABLE rc)
    message("${name}: exit ${rc}\n${out}")
    if(NOT rc EQUAL expected)
        message(FATAL_ERROR "CheckTerrainPoissonDivergenceSelfTest.cmake: ${name} exited ${rc}, expected ${expected}")
    endif()
endforeach()
