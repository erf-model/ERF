# Add uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${PROJECT_SOURCE_DIR}/CMake/cmake_uninstall.cmake.in"
        "${PROJECT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${PROJECT_BINARY_DIR}/cmake_uninstall.cmake
        COMMENT "Uninstalling files listed in install_manifest.txt"
    )
endif()

# Add distclean target
add_custom_target(distclean
    # Header
    COMMAND ${CMAKE_COMMAND} -E echo "=================================================================================="
    COMMAND ${CMAKE_COMMAND} -E echo "Distclean: ${PROJECT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "=================================================================================="

    # CMake configuration files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${PROJECT_BINARY_DIR}/CMakeCache.txt
            ${PROJECT_BINARY_DIR}/cmake_install.cmake
            ${PROJECT_BINARY_DIR}/cmake_uninstall.cmake
            ${PROJECT_BINARY_DIR}/Makefile
            ${PROJECT_BINARY_DIR}/install_manifest.txt
            ${PROJECT_BINARY_DIR}/cray_detected_config.cmake

    # CPack files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${PROJECT_BINARY_DIR}/CPackConfig.cmake
            ${PROJECT_BINARY_DIR}/CPackSourceConfig.cmake

    # CTest files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${PROJECT_BINARY_DIR}/CTestTestfile.cmake
            ${PROJECT_BINARY_DIR}/DartConfiguration.tcl

    # Project-specific generated files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${PROJECT_BINARY_DIR}/ERFConfig.cmake
            ${PROJECT_BINARY_DIR}/compile_commands.json
            ${PROJECT_BINARY_DIR}/git-state.txt

    # CMake-generated directories
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/Testing
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/_deps

    # Build output directories
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/Exec
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/Submodules
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/Tests
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/erf_srclib
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/cmake_packages
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/externals

    # Removing generated files
    COMMAND ${CMAKE_COMMAND} -E echo "Removing generated files from ${PROJECT_BINARY_DIR}..."
    COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/*.pc ${PROJECT_BINARY_DIR}/lib*.a ${PROJECT_BINARY_DIR}/lib*.so ${PROJECT_BINARY_DIR}/build_*.log

    # Summary
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo " DONE: Distclean complete"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Next steps to reconfigure:"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "  If you used a build script:"
    COMMAND ${CMAKE_COMMAND} -E echo "    cd ${PROJECT_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "    ERF_HOME=${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/Build/cmake.sh"
    COMMAND ${CMAKE_COMMAND} -E echo "    or whichever script you used"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "  For manual cmake configuration:"
    COMMAND ${CMAKE_COMMAND} -E echo "    From build directory: cmake ${PROJECT_SOURCE_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Note: Install directories preserved"

    COMMENT "Removing all CMake configuration and build artifacts"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)