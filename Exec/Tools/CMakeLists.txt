if(ERF_ENABLE_NETCDF)
    # PlotFile-NetCDF converter
    set(tool_exe erf_plt_to_nc)
    add_executable(${tool_exe}
        ${PROJECT_SOURCE_DIR}/Source/Tools/ERF_PlotfileToNetCDF.cpp
        ${PROJECT_SOURCE_DIR}/Source/IO/ERF_NCPlotFile.cpp
        ${PROJECT_SOURCE_DIR}/Source/IO/ERF_NCInterface.cpp
    )

    target_include_directories(${tool_exe}
        PRIVATE
        ${PROJECT_SOURCE_DIR}/Source/Tools
        ${PROJECT_SOURCE_DIR}/Source/IO
        ${AMREX_SUBMOD_LOCATION}/Src/Base
        ${AMREX_SUBMOD_LOCATION}/Src/EB
        ${AMREX_SUBMOD_LOCATION}/Src/Extern/HDF5
        ${AMReX_BINARY_DIR}
        ${NETCDF_INCLUDE_DIRS}
    )

    target_link_libraries(${tool_exe}
        PUBLIC
        AMReX::amrex
        ${NETCDF_LINK_LIBRARIES}
    )
    if(ERF_ENABLE_MPI)
        target_link_libraries(${tool_exe} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
    endif()

    include(${PROJECT_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
    set_erf_compile_flags(${tool_exe})

    target_compile_features(${tool_exe} PRIVATE cxx_std_17)

    install(TARGETS ${tool_exe}
        RUNTIME DESTINATION bin
    )

else()
    message(WARNING "${tool_exe} was not compiled because ERF_ENABLE_NETCDF was not set")
endif() # ERF_ENABLE_NETCDF
