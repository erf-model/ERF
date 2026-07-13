if(NOT EXISTS "/home/runner/work/ERF/ERF/Build/install_manifest.txt")
    message(FATAL_ERROR "Cannot find install manifest: /home/runner/work/ERF/ERF/Build/install_manifest.txt")
endif()

file(READ "/home/runner/work/ERF/ERF/Build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")

foreach(file ${files})
    message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
    if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
        execute_process(
            COMMAND "/usr/local/bin/cmake" -E remove "$ENV{DESTDIR}${file}"
            RESULT_VARIABLE rm_retval
            OUTPUT_VARIABLE rm_out
            ERROR_VARIABLE rm_err
        )
        if(NOT "${rm_retval}" STREQUAL "0")
            message(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}: ${rm_err}")
        endif()
    else()
        message(STATUS "File $ENV{DESTDIR}${file} does not exist.")
    endif()
endforeach()

message(STATUS "Uninstall complete")
