set (EKAT_MACH_FILES_PATH "${MACH_FILE_DIR}" CACHE PATH "Path to machine files from bash")
message(STATUS "Searching machine files at: ${EKAT_MACH_FILES_PATH}")    
include (${EKAT_MACH_FILES_PATH}/generic.cmake)