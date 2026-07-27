if(NOT EXISTS "/g/g10/wise14/compiling/clean/ERF/build/install_manifest.txt")
    message(FATAL_ERROR "Cannot find install manifest: /g/g10/wise14/compiling/clean/ERF/build/install_manifest.txt")
endif()

file(READ "/g/g10/wise14/compiling/clean/ERF/build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")

foreach(file ${files})
    message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
    if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
        execute_process(
            COMMAND "/usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.23.1-mdfqd2l7c33zg7xcvqizwz25vqmp7jfw/bin/cmake" -E remove "$ENV{DESTDIR}${file}"
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
