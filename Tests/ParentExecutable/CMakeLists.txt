# A minimal parent CMake project that uses ERF the way ERF uses AMReX: ERF is a library
# (ERF::erf), and the parent builds its own executable from its own sources and its own
# problem file, with erf_setup_executable. It passes none of ERF's model-directory hooks
# (ERF_STRUCTURE_MODEL_DIRS and the others), and ERF builds only what that executable
# needs (add_subdirectory with EXCLUDE_FROM_ALL). It is a test of ERF's CMake, not part of
# the ERF build: configure this folder, not the ERF root.
#
#   cmake -S Tests/ParentExecutable -B <build dir outside the ERF tree>
#   cmake --build <build dir> --target parent_model fvolumesum fcompare
#   ctest --test-dir <build dir> -R ParentExecutable
#
# The parent's own sources are a copy, outside the ERF tree, of the example structure
# model (Tests/StructureModelExample/ERF_StructureModelExample.cpp) and of ERF's problem
# file (Exec/ERF_Prob.cpp, Exec/ERF_Prob.H), as a repository would keep its own. The
# model registers itself (ERF_REGISTER_STRUCTURE_MODEL) from the parent's executable.
# See README.md.

cmake_minimum_required(VERSION 3.25)
project(ERFParentExecutable CXX C)

get_filename_component(_erf_default "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(ERF_PARENT_ERF_DIR "${_erf_default}" CACHE PATH
    "ERF source tree to add with add_subdirectory (Submodules/ERF in a real parent)")
set(ERF_PARENT_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/OutsideSources" CACHE PATH
    "Folder outside the ERF tree that receives the parent's own sources")

if(NOT EXISTS "${ERF_PARENT_ERF_DIR}/CMake/BuildERFExe.cmake")
  message(FATAL_ERROR "ERF_PARENT_ERF_DIR: ${ERF_PARENT_ERF_DIR} is not an ERF source tree")
endif()
# Compare real paths, so a symlinked prefix (/tmp and /private/tmp on macOS) cannot
# hide a source folder that is inside the ERF tree.
get_filename_component(_erf_real "${ERF_PARENT_ERF_DIR}" REALPATH)
get_filename_component(_src_dir "${ERF_PARENT_SOURCE_DIR}" ABSOLUTE)
file(MAKE_DIRECTORY "${_src_dir}")
get_filename_component(_src_real "${_src_dir}" REALPATH)
file(RELATIVE_PATH _src_rel "${_erf_real}" "${_src_real}")
if(NOT _src_rel MATCHES "^\\.\\./")
  message(FATAL_ERROR "ERF_PARENT_SOURCE_DIR: ${_src_dir} is inside the ERF tree; "
                      "use a build directory outside ${ERF_PARENT_ERF_DIR}")
endif()

# The parent's own sources, as a repository would keep them, outside ERF. The copy is
# remade at each configure, and CMake reconfigures when a source file changes.
set(_parent_files
    "${ERF_PARENT_ERF_DIR}/Tests/StructureModelExample/ERF_StructureModelExample.cpp"
    "${ERF_PARENT_ERF_DIR}/Tests/StructureModelExample/inputs_structure_example"
    "${ERF_PARENT_ERF_DIR}/Exec/ERF_Prob.cpp"
    "${ERF_PARENT_ERF_DIR}/Exec/ERF_Prob.H")
file(REMOVE_RECURSE "${_src_dir}")
file(COPY ${_parent_files} DESTINATION "${_src_dir}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_parent_files})

# ERF options, set as normal variables the way a parent project usually sets them. The
# structure coupling is on; no model directory is handed to ERF.
set(ERF_ENABLE_MOVING_BODIES ON)
set(ERF_ENABLE_FCOMPARE ON)   # builds the AMReX plotfile tools; the test uses fvolumesum

enable_testing()
add_subdirectory("${ERF_PARENT_ERF_DIR}" ERF EXCLUDE_FROM_ALL)

# The parent's executable: its own model and problem file, ERF's main.cpp and the ERF
# library (ERF::erf), with ERF's compile flags.
add_executable(parent_model
               "${_src_dir}/ERF_StructureModelExample.cpp"
               "${_src_dir}/ERF_Prob.cpp")
target_include_directories(parent_model PRIVATE "${_src_dir}")
erf_setup_executable(parent_model)

# The case of the example model (Tests/StructureModelExample/inputs_structure_example):
# one point pushes on still air with 50 N along x. At t = 4 s the volume integral of
# x_velocity is F t / rho_0 = 200 m^4/s, which fvolumesum prints to 6 digits, so the check
# is to 1e-5 relative, as in ERF's Structures_ExampleModel and ParentProject tests.
# "structure_example_fx" in the plotfile header shows the parent's model registered and
# plotted; "t = 4 impulse_x = 200" in the log shows its advance() ran every step.
if(NOT TARGET fvolumesum)
  message(FATAL_ERROR "The parent executable test needs the AMReX tool fvolumesum, which ERF "
                      "builds only with its internal AMReX (ERF_USE_INTERNAL_AMREX=ON)")
endif()
if(ERF_PRECISION STREQUAL "DOUBLE")
  set(_run_dir "${CMAKE_CURRENT_BINARY_DIR}/ParentExecutable_StructureExample")
  file(MAKE_DIRECTORY "${_run_dir}")
  add_test(NAME ParentExecutable_StructureExample
           COMMAND sh -c "rm -rf plt0* chk0* && '$<TARGET_FILE:parent_model>' '${_src_dir}/inputs_structure_example' amrex.call_addr2line=0 > run.log 2>&1 && grep -q structure_example_fx plt00016/Header && grep -q 't = 4 impulse_x = 200 ' run.log && '$<TARGET_FILE:fvolumesum>' -v x_velocity plt00016 > momentum.txt && awk 'BEGIN { ok = 0 } /integral of x_velocity/ { v = $NF; d = v - 200.0; if (d < 0) d = -d; print \"integral of x_velocity = \" v \", expected 200\"; ok = (d < 1.0e-5 * 200.0) } END { exit !ok }' momentum.txt"
           WORKING_DIRECTORY "${_run_dir}")
  set_tests_properties(ParentExecutable_StructureExample PROPERTIES TIMEOUT 600)
else()
  message(STATUS "ParentExecutable_StructureExample needs ERF_PRECISION=DOUBLE; not registered")
endif()

# ERF's general CTest helpers with the parent's executable: ERF_TEST_EXE names it for the
# helpers that take no executable argument. The example model stops a run whose kernel
# width is not positive, which shows the helper ran parent_model (the model is compiled
# into no other executable here).
include("${ERF_PARENT_ERF_DIR}/Tests/ERFTestHelpers.cmake")
set(ERF_TEST_EXE parent_model)
add_test_abort(ParentExecutable_HelperAbort "${_src_dir}" inputs_structure_example
               "structure.example.epsilon must be > 0 and finite, got 0"
               "structure.example.epsilon=0")
