# Run Tests/test_files/IBSEB_RefinedLevels twice to its step-0 report -- with level 1
# holding the cube only (the deck's box) and holding the cube and the tower -- then with
# level 1 around the tower only and a material per building, and two steps with level 1
# holding no building, and check
# with check_ibseb_refined_levels.py that the cube's faces on level 1 see the tower
# through level 0's column map: the same shadow on every face and view fractions within
# 5 of the 128 hemisphere rays (IBFaceSet::add_outside_occluders()). Before running, the deck's box is checked
# with Exec/CanonicalTests/SEB/ibseb_refinement_box.py, the script users make such boxes
# with: it must accept it and two boxes that split the cube between them, reject a box
# whose edge cuts the cube or a 4 m block only level 1 resolves, and stop on a proposal
# without --fit and on a deck without a height map.
# -DX= defines X as empty, so test for a value, not for DEFINED
include("${CMAKE_CURRENT_LIST_DIR}/MPILauncher.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/ResolveExecutable.cmake")

foreach(arg NRANKS TEST_EXE INPUT WORKING_DIRECTORY PYTHON_EXE CHECKER BOX_SCRIPT)
    if("${${arg}}" STREQUAL "")
        message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: ${arg} must be given and non-empty")
    endif()
endforeach()
if(NOT "${MPIEXEC}" STREQUAL "" AND "${MPIEXEC_NUMPROC_FLAG}" STREQUAL "")
    message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: MPIEXEC_NUMPROC_FLAG must be given with MPIEXEC")
endif()

erf_resolve_executable(TEST_EXE "${TEST_EXE}" CONFIG "${CONFIG}"
    CONTEXT "RunIBSEBRefinedLevels.cmake: ERF executable")

erf_mpi_launcher_command(launch
    LAUNCHER "${MPIEXEC}"
    NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}"
    NRANKS ${NRANKS}
    PREFLAGS "${MPIEXEC_PREFLAGS}"
    CONTEXT "RunIBSEBRefinedLevels.cmake")

get_filename_component(input_dir "${INPUT}" DIRECTORY)
set(ONE_DIR  "${WORKING_DIRECTORY}/cube_only")
set(BOTH_DIR "${WORKING_DIRECTORY}/both")
file(REMOVE_RECURSE "${ONE_DIR}" "${BOTH_DIR}")
file(MAKE_DIRECTORY "${ONE_DIR}/faces" "${BOTH_DIR}/faces")

# The refinement-box script: the deck's box passes (exit 0), a box cut through the cube
# does not (exit 1).
execute_process(COMMAND "${PYTHON_EXE}" "${BOX_SCRIPT}" "${INPUT}"
                OUTPUT_VARIABLE box_out ERROR_VARIABLE box_out RESULT_VARIABLE box_result)
message("${box_out}")
if(NOT box_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: ibseb_refinement_box.py rejects the deck's box (${box_result})")
endif()
file(READ "${INPUT}" deck_text)
string(REPLACE "erf.city.in_box_lo = 100.0 180.0" "erf.city.in_box_lo = 160.0 180.0" cut_text "${deck_text}")
if("${cut_text}" STREQUAL "${deck_text}")
    message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: the deck's erf.city.in_box_lo line changed; update the cut box")
endif()
file(WRITE "${WORKING_DIRECTORY}/cut.i" "${cut_text}")
execute_process(COMMAND "${PYTHON_EXE}" "${BOX_SCRIPT}" "${WORKING_DIRECTORY}/cut.i"
                        --buildings "${input_dir}/cube_and_tower_10m.txt"
                OUTPUT_VARIABLE cut_out ERROR_VARIABLE cut_out RESULT_VARIABLE cut_result)
message("${cut_out}")
if(NOT cut_result EQUAL 1)
    message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: ibseb_refinement_box.py accepts a box whose edge cuts the cube (${cut_result})")
endif()

# Level 1 is the union of its boxes: two boxes splitting the cube between them hold it.
string(REPLACE "erf.refinement_indicators = city" "erf.refinement_indicators = west east" split_text "${deck_text}")
string(REPLACE "erf.city.max_level = 1" "erf.west.max_level = 1\nerf.east.max_level = 1" split_text "${split_text}")
string(REPLACE "erf.city.in_box_lo = 100.0 180.0" "erf.west.in_box_lo = 100.0 180.0\nerf.east.in_box_lo = 160.0 180.0" split_text "${split_text}")
string(REPLACE "erf.city.in_box_hi = 220.0 300.0" "erf.west.in_box_hi = 160.0 300.0\nerf.east.in_box_hi = 220.0 300.0" split_text "${split_text}")
file(WRITE "${WORKING_DIRECTORY}/split.i" "${split_text}")
# The block only level 1 resolves still crosses the box (the script's coarse footprint is
# the most a building can occupy); a proposal needs --fit; a deck without a height map stops.
foreach(case "split.i;0" "cut_low;1" "fit;2" "nomap;1")
    list(GET case 0 name)
    list(GET case 1 expected)
    if(name STREQUAL "split.i")
        set(cmd "${PYTHON_EXE}" "${BOX_SCRIPT}" "${WORKING_DIRECTORY}/split.i" --buildings "${input_dir}/cube_and_tower_10m.txt")
    elseif(name STREQUAL "cut_low")
        set(cmd "${PYTHON_EXE}" "${BOX_SCRIPT}" "${INPUT}" --buildings "${input_dir}/cube_and_low_block_10m.txt")
    elseif(name STREQUAL "fit")
        set(cmd "${PYTHON_EXE}" "${BOX_SCRIPT}" "${INPUT}" --all)
    else()
        file(WRITE "${WORKING_DIRECTORY}/nomap.i" "amr.n_cell = 24 24 16\ngeometry.prob_extent = 480 480 160\n")
        set(cmd "${PYTHON_EXE}" "${BOX_SCRIPT}" "${WORKING_DIRECTORY}/nomap.i")
    endif()
    execute_process(COMMAND ${cmd} OUTPUT_VARIABLE out ERROR_VARIABLE out RESULT_VARIABLE rc)
    message("ibseb_refinement_box.py ${name}: exit ${rc}\n${out}")
    if(NOT rc EQUAL expected)
        message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: ibseb_refinement_box.py ${name} exited ${rc}, expected ${expected}")
    endif()
endforeach()

function(run_leg dir)
    # The deck reads its sounding and height map by relative names.
    file(COPY "${input_dir}/input_sounding" "${input_dir}/cube_and_tower_10m.txt" DESTINATION "${dir}")
    execute_process(
        COMMAND ${launch} ${TEST_EXE} ${INPUT} ${ARGN}
        WORKING_DIRECTORY "${dir}"
        OUTPUT_FILE "${dir}/simulation.log"
        ERROR_FILE "${dir}/simulation.log"
        TIMEOUT 600
        RESULT_VARIABLE _result)
    if(NOT _result EQUAL 0)
        message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: the run in ${dir} failed: ${_result} (see ${dir}/simulation.log)")
    endif()
endfunction()

run_leg("${ONE_DIR}")
run_leg("${BOTH_DIR}" "erf.city.in_box_hi=360.0 300.0")
# Level 1 around the tower only, with a material per building: the tower is building 2 on
# level 0 and building 1 on level 1, and must keep level 0's material for building 2.
set(TOWER_DIR "${WORKING_DIRECTORY}/tower")
file(REMOVE_RECURSE "${TOWER_DIR}")
file(MAKE_DIRECTORY "${TOWER_DIR}/faces")
file(COPY "${input_dir}/materials.csv" DESTINATION "${TOWER_DIR}")
run_leg("${TOWER_DIR}" "erf.city.in_box_lo=220.0 180.0" "erf.city.in_box_hi=340.0 300.0"
        "erf.ibseb.material_file=materials.csv" "erf.ibseb.material_by_building=1 2")
# A refined level that holds no building, with the debug summaries on (its means had
# divided by a face count of zero).
set(EMPTY_DIR "${WORKING_DIRECTORY}/empty")
file(REMOVE_RECURSE "${EMPTY_DIR}")
file(MAKE_DIRECTORY "${EMPTY_DIR}/faces")
run_leg("${EMPTY_DIR}" "erf.city.in_box_lo=320.0 60.0" "erf.city.in_box_hi=420.0 140.0"
        "max_step=2" "erf.ibseb.debug=true")

execute_process(
    COMMAND "${PYTHON_EXE}" "${CHECKER}" --both "${BOTH_DIR}" --one "${ONE_DIR}" --tower "${TOWER_DIR}"
    WORKING_DIRECTORY "${WORKING_DIRECTORY}"
    OUTPUT_FILE "${WORKING_DIRECTORY}/checker.log"
    ERROR_FILE "${WORKING_DIRECTORY}/checker.log"
    RESULT_VARIABLE check_result)
file(READ "${WORKING_DIRECTORY}/checker.log" check_output)
message("${check_output}")
if(NOT check_result EQUAL 0)
    message(FATAL_ERROR "RunIBSEBRefinedLevels.cmake: the checker failed: ${check_result}")
endif()
