foreach(name ERF_EXECUTABLE CHECK_EXECUTABLE INPUT_FILE TEST_ROOT)
  if(NOT DEFINED ${name} OR "${${name}}" STREQUAL "")
    message(FATAL_ERROR "${name} is required")
  endif()
endforeach()

foreach(case "active;800" "coarse;960" "active_native;800" "coarse_native;960" "coarse_ibseb;960")
  list(GET case 0 mode)
  list(GET case 1 extent)
  set(factory_options)
  if(mode MATCHES "_native$")
    list(APPEND factory_options
      "erf.buildings_fc_factories=true" "erf.plot_face_terrain_blanking=true")
  elseif(mode STREQUAL "coarse_ibseb")
    list(APPEND factory_options
      "erf.terrain_type=None" "erf.plot_face_terrain_blanking=true"
      "erf.ibseb.enable=true" "erf.ibseb.prognostic=false"
      "erf.ibseb.radiation=prescribed")
  endif()
  set(run_dir "${TEST_ROOT}/immersed_building_cutoff_${mode}")
  file(REMOVE_RECURSE "${run_dir}")
  file(MAKE_DIRECTORY "${run_dir}")
  execute_process(
    COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
      "erf.prob_name=bubble" "prob.U_0=3" "max_step=0"
      "geometry.prob_extent=${extent} ${extent} 160"
      "eb2.box_lo=200 200 -5" "eb2.box_hi=400 400 50"
      "erf.plot_int_1=1" "erf.plot_file_1=plt"
      "erf.plot_vars_1=x_velocity terrain_IB_mask"
      ${factory_options}
    WORKING_DIRECTORY "${run_dir}"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${mode} run failed (${result}):\n${stdout}\n${stderr}")
  endif()
  if(mode STREQUAL "coarse_ibseb" AND NOT stdout MATCHES "faces=0")
    message(FATAL_ERROR "coarse IBSEB run did not report zero building faces:\n${stdout}")
  endif()

  execute_process(
    COMMAND "${CHECK_EXECUTABLE}" "${run_dir}/plt00000" "${mode}"
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "${mode} plotfile check failed (${result}):\n${stdout}\n${stderr}")
  endif()
endforeach()
