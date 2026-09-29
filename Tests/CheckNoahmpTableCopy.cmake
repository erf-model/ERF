# Require Exec/RegTests/NoahMP_Ideal/NoahmpTable.TBL to match the Noah-MP submodule's table.
#
# The nightly NoahMP_Ideal case carries its own copy so that it runs in place, and its stored
# reference is only meaningful for the parameters it was made with. After a submodule bump a
# stale copy would keep testing the old parameters without anyone noticing -- the older copies
# elsewhere in the repository have already drifted -- so a mismatch fails here, in CI, instead.
#
# Inputs: SUBMODULE_TABLE, COPY_TABLE (absolute paths).

foreach(f IN ITEMS "${SUBMODULE_TABLE}" "${COPY_TABLE}")
  if(NOT EXISTS "${f}")
    message(FATAL_ERROR "NoahmpTable.TBL check: '${f}' does not exist.")
  endif()
endforeach()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol "${SUBMODULE_TABLE}" "${COPY_TABLE}"
  RESULT_VARIABLE rc)

if(NOT rc EQUAL 0)
  message(FATAL_ERROR
    "Exec/RegTests/NoahMP_Ideal/NoahmpTable.TBL no longer matches the Noah-MP submodule's "
    "table. Copy the submodule's table over it:\n"
    "  cp Submodules/Noah-MP/parameters/NoahmpTable.TBL Exec/RegTests/NoahMP_Ideal/\n"
    "and regenerate the nightly NoahMP_Ideal reference, since its answer depends on the "
    "table's parameters.")
endif()
message(STATUS "NoahmpTable.TBL copy matches the Noah-MP submodule.")
