# A minimal parent CMake project that builds ERF with add_subdirectory, the way a
# repository that holds ERF as a git submodule would. It is a test of ERF's CMake,
# not part of the ERF build: configure this folder, not the ERF root.
#
#   cmake -S Tests/ParentProject -B <build dir outside the ERF tree>
#   cmake --build <build dir> --target erf_exec fvolumesum
#   ctest --test-dir <build dir> -R ParentProject
#
# The example structure model (Tests/StructureModelExample) is copied to a folder
# outside the ERF tree and compiled in through ERF_STRUCTURE_MODEL_DIRS. See README.md.

# 3.25, as ERF's internal AMReX. (3.21 is the least for the normal variables set below:
# it has policy CMP0126, without which ERF's set(... CACHE ...) drops
# ERF_STRUCTURE_MODEL_DIRS on the first configure.)
cmake_minimum_required(VERSION 3.25)
project(ERFParentProject CXX C)

get_filename_component(_erf_default "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(ERF_PARENT_ERF_DIR "${_erf_default}" CACHE PATH
    "ERF source tree to add with add_subdirectory (Submodules/ERF in a real parent)")
set(ERF_PARENT_MODEL_DIR "${CMAKE_CURRENT_BINARY_DIR}/OutsideModel" CACHE PATH
    "Folder outside the ERF tree that receives a copy of Tests/StructureModelExample")

if(NOT EXISTS "${ERF_PARENT_ERF_DIR}/CMake/BuildERFExe.cmake")
  message(FATAL_ERROR "ERF_PARENT_ERF_DIR: ${ERF_PARENT_ERF_DIR} is not an ERF source tree")
endif()
# Compare real paths, so a symlinked prefix (/tmp and /private/tmp on macOS) cannot
# hide a model folder that is inside the ERF tree.
get_filename_component(_erf_real "${ERF_PARENT_ERF_DIR}" REALPATH)
get_filename_component(_model_dir "${ERF_PARENT_MODEL_DIR}" ABSOLUTE)
file(MAKE_DIRECTORY "${_model_dir}")
get_filename_component(_model_real "${_model_dir}" REALPATH)
file(RELATIVE_PATH _model_rel "${_erf_real}" "${_model_real}")
if(NOT _model_rel MATCHES "^\\.\\./")
  message(FATAL_ERROR "ERF_PARENT_MODEL_DIR: ${_model_dir} is inside the ERF tree; "
                      "use a build directory outside ${ERF_PARENT_ERF_DIR}")
endif()

# The model folder as a parent repository would keep it, outside ERF. The copy is
# remade at each configure, and CMake reconfigures when a source file changes.
set(_model_src "${ERF_PARENT_ERF_DIR}/Tests/StructureModelExample")
file(GLOB _model_files LIST_DIRECTORIES false "${_model_src}/*")
file(REMOVE_RECURSE "${_model_dir}")
file(COPY ${_model_files} DESTINATION "${_model_dir}")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_model_files} "${_model_src}")

# ERF options, set as normal variables the way a parent project usually sets them.
set(ERF_ENABLE_MOVING_BODIES ON)
set(ERF_STRUCTURE_MODEL_DIRS "${_model_dir}")
set(ERF_ENABLE_FCOMPARE ON)   # builds the AMReX plotfile tools; the test uses fvolumesum

enable_testing()
add_subdirectory("${ERF_PARENT_ERF_DIR}" ERF)

# One point pushes on still air with 50 N along x in a periodic box with no friction
# or gravity (Tests/StructureModelExample/inputs_structure_example). At t = 4 s the
# x-momentum of the air is F t = 200 kg m/s, to round-off, and the volume integral of
# x_velocity is F t / rho_0 = 200 m^4/s to about 1e-6 relative (density changes of
# order u/c). fvolumesum prints 6 significant digits, so the check is to 1e-5 relative,
# as in ERF's own Structures_ExampleModel test. "structure_example_fx" in the plotfile
# header shows the model is compiled in; "t = 4 impulse_x = 200" in the log shows its
# advance() ran every step. Like Structures_ExampleModel, the test needs double precision.
if(NOT TARGET fvolumesum)
  message(FATAL_ERROR "The parent project test needs the AMReX tool fvolumesum, which ERF "
                      "builds only with its internal AMReX (ERF_USE_INTERNAL_AMREX=ON)")
endif()
if(ERF_PRECISION STREQUAL "DOUBLE")
  set(_run_dir "${CMAKE_CURRENT_BINARY_DIR}/ParentProject_StructureExample")
  file(MAKE_DIRECTORY "${_run_dir}")
  add_test(NAME ParentProject_StructureExample
           COMMAND sh -c "rm -rf plt0* chk0* && '$<TARGET_FILE:erf_exec>' '${_model_dir}/inputs_structure_example' amrex.call_addr2line=0 > run.log 2>&1 && grep -q structure_example_fx plt00016/Header && grep -q 't = 4 impulse_x = 200 ' run.log && '$<TARGET_FILE:fvolumesum>' -v x_velocity plt00016 > momentum.txt && awk 'BEGIN { ok = 0 } /integral of x_velocity/ { v = $NF; d = v - 200.0; if (d < 0) d = -d; print \"integral of x_velocity = \" v \", expected 200\"; ok = (d < 1.0e-5 * 200.0) } END { exit !ok }' momentum.txt"
           WORKING_DIRECTORY "${_run_dir}")
  set_tests_properties(ParentProject_StructureExample PROPERTIES TIMEOUT 600)
else()
  message(STATUS "ParentProject_StructureExample needs ERF_PRECISION=DOUBLE; not registered")
endif()
