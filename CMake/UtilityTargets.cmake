# Add uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_SOURCE_DIR}/CMake/cmake_uninstall.cmake.in"
        "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake
        COMMENT "Uninstalling files listed in install_manifest.txt"
    )
endif()

# Add distclean target
add_custom_target(distclean
    # Header
    COMMAND ${CMAKE_COMMAND} -E echo "=================================================================================="
    COMMAND ${CMAKE_COMMAND} -E echo "Distclean: ${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "=================================================================================="

    # CMake configuration files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${CMAKE_BINARY_DIR}/CMakeCache.txt
            ${CMAKE_BINARY_DIR}/cmake_install.cmake
            ${CMAKE_BINARY_DIR}/cmake_uninstall.cmake
            ${CMAKE_BINARY_DIR}/Makefile
            ${CMAKE_BINARY_DIR}/install_manifest.txt
            ${CMAKE_BINARY_DIR}/cray_detected_config.cmake

    # CPack files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${CMAKE_BINARY_DIR}/CPackConfig.cmake
            ${CMAKE_BINARY_DIR}/CPackSourceConfig.cmake

    # CTest files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${CMAKE_BINARY_DIR}/CTestTestfile.cmake
            ${CMAKE_BINARY_DIR}/DartConfiguration.tcl

    # Project-specific generated files
    COMMAND ${CMAKE_COMMAND} -E remove -f
            ${CMAKE_BINARY_DIR}/ERFConfig.cmake
            ${CMAKE_BINARY_DIR}/compile_commands.json
            ${CMAKE_BINARY_DIR}/git-state.txt

    # CMake-generated directories
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/Testing
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/_deps

    # Build output directories
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/Exec
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/Submodules
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/Tests
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/erf_srclib
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/cmake_packages
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/externals

    # Use shell commands with explicit directory (if needed for glob patterns)
    COMMAND ${CMAKE_COMMAND} -E echo "Removing generated files..."
    COMMAND sh -c "cd ${CMAKE_BINARY_DIR} && rm -f *.pc lib*.a lib*.so build_*.log 2>/dev/null || true"

    # Summary
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo " DONE: Distclean complete"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Next steps to reconfigure:"
    COMMAND ${CMAKE_COMMAND} -E echo "  From build directory: cmake ${CMAKE_SOURCE_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo "  Or more simply:       cmake .."
    COMMAND ${CMAKE_COMMAND} -E echo "  From source directory: cmake -B ${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Note: Install directories preserved"

    COMMENT "Removing all CMake configuration and build artifacts"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)