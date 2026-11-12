# A minimal parent CMake project that builds ERF with add_subdirectory, the way a
# repository that holds ERF as a git submodule would. It is a test of ERF's CMake,
# not part of the ERF build: configure this folder, not the ERF root.
#
#   cmake -S Tests/ParentProject -B <build dir outside the ERF tree>
#   cmake --build <build dir> --target erf_exec fvolumesum fcompare
#   ctest --test-dir <build dir> -R ParentProject
#
# The example structure model (Tests/StructureModelExample) is copied to a folder
# outside the ERF tree and compiled in through ERF_STRUCTURE_MODEL_DIRS. Three more tests
# register through ERF's general CTest helpers (Tests/ERFTestHelpers.cmake), the way a
# parent repository registers its own cases, and one more runs the model's GoogleTest
# when ERF_ENABLE_UNIT_TESTS=ON. See README.md.

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
  # The plotfile that passes is also the gold of the two add_test_r tests below.
  set(_gold_dir "${CMAKE_CURRENT_BINARY_DIR}/GoldFiles")
  add_test(NAME ParentProject_StructureExample
           COMMAND sh -c "rm -rf plt0* chk0* '${_gold_dir}' && '$<TARGET_FILE:erf_exec>' '${_model_dir}/inputs_structure_example' amrex.call_addr2line=0 > run.log 2>&1 && grep -q structure_example_fx plt00016/Header && grep -q 't = 4 impulse_x = 200 ' run.log && '$<TARGET_FILE:fvolumesum>' -v x_velocity plt00016 > momentum.txt && awk 'BEGIN { ok = 0 } /integral of x_velocity/ { v = $NF; d = v - 200.0; if (d < 0) d = -d; print \"integral of x_velocity = \" v \", expected 200\"; ok = (d < 1.0e-5 * 200.0) } END { exit !ok }' momentum.txt && mkdir -p '${_gold_dir}' && cp -R plt00016 '${_gold_dir}/ParentProject_HelperRegression' && cp -R plt00016 '${_gold_dir}/ParentProject_HelperRegressionDetects'"
           WORKING_DIRECTORY "${_run_dir}")
  set_tests_properties(ParentProject_StructureExample PROPERTIES TIMEOUT 600 FIXTURES_SETUP ParentProjectGold)
else()
  message(STATUS "ParentProject_StructureExample needs ERF_PRECISION=DOUBLE; not registered")
endif()

# ERF's general CTest helpers, included in this directory the way a parent repository
# registers its own cases: the inputs in test_files/<test name> of this folder, the gold
# plotfiles in a gold folder of the parent's own.
include("${ERF_PARENT_ERF_DIR}/Tests/ERFTestHelpers.cmake")

# add_test_abort: the example model stops a run whose kernel width is not positive.
add_test_abort(ParentProject_HelperAbort "${_model_dir}" inputs_structure_example
               "structure.example.epsilon must be > 0 and finite, got 0"
               "structure.example.epsilon=0")

# add_test_r: the same case as ParentProject_StructureExample, from this folder's
# test_files, compared by fcompare with the plotfile that test checked against the closed
# form (fixture ParentProjectGold). The normal variable ERF_TEST_GOLD_FILES_DIRECTORY,
# set after add_subdirectory so that ERF's own tests keep ERF's gold folder, points
# add_test_r at the parent's gold folder; the tolerances are ERF's defaults, set here
# only when ERF's tests (ERF_ENABLE_TESTS) have not set them. The gold comes from the
# same executable and deck, so this test checks that add_test_r works from a parent
# (input copy, launcher, gold path, fcompare), not the physics; the closed form above
# checks the physics. ParentProject_HelperRegressionDetects runs the case with 51 N
# instead of 50 N against the same gold: the flow is linear in the force, so fcompare
# must report a relative difference of 1/50 = 0.02 in x_velocity. Matching that number,
# rather than only a failure, keeps a run that fails for another reason (no executable,
# no gold) from passing it. The gold is written on 1 rank, and
# ERF_TEST_NRANKS = 1 runs these tests on 1 rank too with ERF_ENABLE_MPI.
if(ERF_PRECISION STREQUAL "DOUBLE")
  set(ERF_TEST_GOLD_FILES_DIRECTORY "${_gold_dir}")
  foreach(_tol IN ITEMS ERF_TEST_FCOMPARE_RTOL ERF_TEST_FCOMPARE_ATOL)
    if(NOT DEFINED ${_tol})
      set(${_tol} 2.0e-10)
    endif()
  endforeach()
  set(ERF_TEST_NRANKS 1)
  add_test_r(ParentProject_HelperRegression "" erf_exec plt00016)
  set_tests_properties(ParentProject_HelperRegression PROPERTIES FIXTURES_REQUIRED ParentProjectGold)
  add_test_r(ParentProject_HelperRegressionDetects "" erf_exec plt00016
             RUNTIME_OPTIONS "structure.example.force=51 0 0")
  set_tests_properties(ParentProject_HelperRegressionDetects PROPERTIES
                       FIXTURES_REQUIRED ParentProjectGold
                       PASS_REGULAR_EXPRESSION "x_velocity +[0-9.e+-]+ +0\\.0200")
endif()

# The model's GoogleTest, added to erf_unit_tests through its
# erf_structure_model_add_unit_tests hook. The test passes only when exactly one test of
# the suite ran and passed, so a hook that adds nothing fails it ("0 tests").
if(TARGET erf_unit_tests)
  add_test(NAME ParentProject_ModelUnitTest
           COMMAND $<TARGET_FILE:erf_unit_tests> --gtest_filter=StructureModelExample.*)
  set_tests_properties(ParentProject_ModelUnitTest PROPERTIES TIMEOUT 300
                       PASS_REGULAR_EXPRESSION "\\[  PASSED  \\] 1 test\\."
                       FAIL_REGULAR_EXPRESSION "\\[  FAILED  \\]")
endif()
