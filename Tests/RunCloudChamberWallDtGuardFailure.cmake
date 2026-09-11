if(NOT DEFINED TEST_EXE OR NOT DEFINED LOG)
    message(FATAL_ERROR "TEST_EXE and LOG are required")
endif()

execute_process(
    COMMAND "${TEST_EXE}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr)

set(output "${stdout}\n${stderr}")
file(WRITE "${LOG}" "${output}")

if("${result}" STREQUAL "0")
    message(FATAL_ERROR
        "Cloud Chamber fixed-dt guard unexpectedly returned success")
endif()

foreach(expected IN ITEMS
        "Cloud Chamber wall-transfer timestep violation"
        "fixed_dt="
        "wall_dt="
        "max_wall_rate=")
    string(FIND "${output}" "${expected}" position)
    if(position EQUAL -1)
        message(FATAL_ERROR
            "Missing expected fixed-dt guard diagnostic: ${expected}\n${output}")
    endif()
endforeach()
