# Numeric-aware comparison of two whitespace-separated text tables, e.g. the data logs written
# by erf.data_log.
#
# The logs are printed with std::setprecision(datprecision) (Source/ERF.H), and a planar sum
# reduces in a different order for each domain decomposition, so two runs that agree to roundoff
# can still print a different last digit.  A byte-for-byte comparison of the logs therefore fails
# spuriously on a value that happens to sit near a rounding boundary.  This compares numeric
# fields to within a few units of the last significant digit and everything else, the column
# headers above all, as text.
#
# CMake has no floating-point arithmetic, so a number is read as an integer digit string times a
# power of ten: |value| = digits * 10^exponent.  Two values are brought onto a common power of
# ten and their integer difference is compared with the tolerance expressed in the same units.
# Comparisons alone (if(LESS)) would not do, as the difference itself has to be formed.

# Below 10^this a value counts as zero, so that a sum that cancels to roundoff in one run and to
# exactly zero in the other still compares equal
set(ERF_DATALOG_ZERO_EXPONENT -12)
# A number of more than this many digits once the two are on a common power of ten no longer fits
# in the integer arithmetic; such a pair is orders of magnitude apart and cannot be close anyway
set(ERF_DATALOG_MAX_DIGITS 18)

# Read a decimal number, in fixed or scientific notation, as sign * digits * 10^exponent, with
# the leading zeros of the digit string removed.  Sets OK_VAR to FALSE if the string is not a
# number, which is how a text field, a column header say, is recognized.
function(erf_read_decimal value ok_var sign_var digits_var exponent_var)
    set(${ok_var} FALSE PARENT_SCOPE)
    if(NOT "${value}" MATCHES "^([+-]?)([0-9]*)(\\.([0-9]*))?([eE]([+-]?[0-9]+))?$")
        return()
    endif()
    set(_sign "${CMAKE_MATCH_1}")
    set(_int_part "${CMAKE_MATCH_2}")
    set(_frac_part "${CMAKE_MATCH_4}")
    set(_exponent "${CMAKE_MATCH_6}")
    # The regex also matches "", "+", "." and "e5", none of which is a number
    if("${_int_part}${_frac_part}" STREQUAL "")
        return()
    endif()
    if("${_exponent}" STREQUAL "")
        set(_exponent 0)
    endif()
    string(REGEX REPLACE "^\\+" "" _exponent "${_exponent}")
    string(LENGTH "${_frac_part}" _nfrac)
    math(EXPR _exponent "${_exponent} - ${_nfrac}")
    string(REGEX REPLACE "^0+" "" _digits "${_int_part}${_frac_part}")
    if("${_digits}" STREQUAL "")
        set(_digits 0)
        set(_sign "")
    endif()
    if("${_sign}" STREQUAL "+")
        set(_sign "")
    endif()
    set(${ok_var} TRUE PARENT_SCOPE)
    set(${sign_var} "${_sign}" PARENT_SCOPE)
    set(${digits_var} "${_digits}" PARENT_SCOPE)
    set(${exponent_var} "${_exponent}" PARENT_SCOPE)
endfunction()

# Whether two numbers agree to within ULPS units of their sigdigits-th significant digit
function(erf_numbers_close a b sigdigits ulps out_var)
    set(${out_var} FALSE PARENT_SCOPE)
    if("${a}" STREQUAL "${b}")
        set(${out_var} TRUE PARENT_SCOPE)
        return()
    endif()
    erf_read_decimal("${a}" _a_ok _a_sign _a_digits _a_exp)
    erf_read_decimal("${b}" _b_ok _b_sign _b_digits _b_exp)
    if(NOT _a_ok OR NOT _b_ok)
        return()
    endif()

    # The power of ten of the leading digit, i.e. the exponent of the scientific notation
    string(LENGTH "${_a_digits}" _a_ndig)
    string(LENGTH "${_b_digits}" _b_ndig)
    math(EXPR _a_lead "${_a_exp} + ${_a_ndig} - 1")
    math(EXPR _b_lead "${_b_exp} + ${_b_ndig} - 1")
    set(_a_zero FALSE)
    set(_b_zero FALSE)
    if("${_a_digits}" STREQUAL "0" OR _a_lead LESS ${ERF_DATALOG_ZERO_EXPONENT})
        set(_a_zero TRUE)
    endif()
    if("${_b_digits}" STREQUAL "0" OR _b_lead LESS ${ERF_DATALOG_ZERO_EXPONENT})
        set(_b_zero TRUE)
    endif()
    if(_a_zero OR _b_zero)
        if(_a_zero AND _b_zero)
            set(${out_var} TRUE PARENT_SCOPE)
        endif()
        return()
    endif()
    if(NOT "${_a_sign}" STREQUAL "${_b_sign}")
        return()
    endif()

    # Bring both onto the smaller power of ten by appending that many zeros
    set(_common ${_a_exp})
    if(_b_exp LESS _common)
        set(_common ${_b_exp})
    endif()
    math(EXPR _a_shift "${_a_exp} - ${_common}")
    math(EXPR _b_shift "${_b_exp} - ${_common}")
    if(_a_shift GREATER ${ERF_DATALOG_MAX_DIGITS} OR _b_shift GREATER ${ERF_DATALOG_MAX_DIGITS})
        return()
    endif()
    string(REPEAT "0" ${_a_shift} _a_pad)
    string(REPEAT "0" ${_b_shift} _b_pad)
    set(_a_int "${_a_sign}${_a_digits}${_a_pad}")
    set(_b_int "${_b_sign}${_b_digits}${_b_pad}")
    string(LENGTH "${_a_digits}${_a_pad}" _a_len)
    string(LENGTH "${_b_digits}${_b_pad}" _b_len)
    if(_a_len GREATER ${ERF_DATALOG_MAX_DIGITS} OR _b_len GREATER ${ERF_DATALOG_MAX_DIGITS})
        return()
    endif()

    # ULPS units of the last significant digit of the larger of the two, in those same units
    set(_lead ${_a_lead})
    if(_b_lead GREATER _lead)
        set(_lead ${_b_lead})
    endif()
    math(EXPR _tol_shift "${_lead} - ${sigdigits} + 1 - ${_common}")
    if(_tol_shift LESS 0)
        # Both are printed to more digits than sigdigits: nothing below one unit is tolerated
        set(_tolerance 0)
    else()
        string(REPEAT "0" ${_tol_shift} _tol_pad)
        math(EXPR _tolerance "${ulps}${_tol_pad}")
    endif()

    math(EXPR _difference "(${_a_int}) - (${_b_int})")
    if(_difference LESS 0)
        math(EXPR _difference "0 - ${_difference}")
    endif()
    if(NOT _difference GREATER ${_tolerance})
        set(${out_var} TRUE PARENT_SCOPE)
    endif()
endfunction()

# Compare two files field by field.  Sets OUT_VAR to TRUE if they agree, and MESSAGE_VAR to the
# first disagreement otherwise.
function(erf_compare_data_logs file_a file_b sigdigits ulps out_var message_var)
    set(${out_var} FALSE PARENT_SCOPE)
    file(STRINGS "${file_a}" _lines_a)
    file(STRINGS "${file_b}" _lines_b)
    list(LENGTH _lines_a _na)
    list(LENGTH _lines_b _nb)
    if(NOT _na EQUAL _nb)
        set(${message_var} "the logs have ${_na} and ${_nb} lines" PARENT_SCOPE)
        return()
    endif()
    math(EXPR _last "${_na} - 1")
    foreach(_i RANGE 0 ${_last})
        list(GET _lines_a ${_i} _line_a)
        list(GET _lines_b ${_i} _line_b)
        if("${_line_a}" STREQUAL "${_line_b}")
            continue()
        endif()
        string(REGEX MATCHALL "[^ \t]+" _fields_a "${_line_a}")
        string(REGEX MATCHALL "[^ \t]+" _fields_b "${_line_b}")
        list(LENGTH _fields_a _nfa)
        list(LENGTH _fields_b _nfb)
        math(EXPR _line_number "${_i} + 1")
        if(NOT _nfa EQUAL _nfb)
            set(${message_var}
                "line ${_line_number} has ${_nfa} and ${_nfb} fields\n  ${_line_a}\n  ${_line_b}"
                PARENT_SCOPE)
            return()
        endif()
        math(EXPR _last_field "${_nfa} - 1")
        foreach(_j RANGE 0 ${_last_field})
            list(GET _fields_a ${_j} _field_a)
            list(GET _fields_b ${_j} _field_b)
            erf_numbers_close("${_field_a}" "${_field_b}" ${sigdigits} ${ulps} _close)
            if(NOT _close)
                math(EXPR _column "${_j} + 1")
                set(${message_var}
                    "line ${_line_number}, column ${_column}: '${_field_a}' and '${_field_b}' differ by more than ${ulps} units of the ${sigdigits}th significant digit\n  ${_line_a}\n  ${_line_b}"
                    PARENT_SCOPE)
                return()
            endif()
        endforeach()
    endforeach()
    set(${out_var} TRUE PARENT_SCOPE)
endfunction()
