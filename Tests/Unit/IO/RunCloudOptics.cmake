# Radiation-off Kessler writer regression: both streams, reversed/duplicate
# requests, single requests, z-split grids, Headers, metadata, and native values.
file(MAKE_DIRECTORY "${TEST_ROOT}")
# Bubble's existing moist initializer requires full-height boxes. Initialize
# once on supported grids, then use ERF's restart regridding to exercise the
# writer on genuinely z-split grids without changing prognostic initialization.
set(seed "${TEST_ROOT}/seed")
file(MAKE_DIRECTORY "${seed}")
execute_process(COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
  max_step=1 "amr.n_cell=8 4 8" amr.blocking_factor=2
  amr.max_grid_size=8 amr.max_grid_size_z=8 amr.no_box_split_dir=2
  erf.fixed_dt=0.01 erf.fixed_mri_dt_ratio=2 amrex.fpe_trap_invalid=0
  erf.radiation_model=None erf.plot_int_1=-1 erf.plot2d_int_1=-1
  erf.check_int=1 erf.check_file=chk
  WORKING_DIRECTORY "${seed}" OUTPUT_FILE "${seed}/run.log"
  ERROR_FILE "${seed}/error.log" RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT EXISTS "${seed}/chk00001/Header")
  message(FATAL_ERROR "Cloud optics seed checkpoint failed (${status}); see ${seed}")
endif()
function(run_case name vars)
  set(dir "${TEST_ROOT}/${name}")
  file(MAKE_DIRECTORY "${dir}")
  execute_process(COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
    max_step=2 "amr.n_cell=8 4 8" amr.blocking_factor=2
    "erf.restart=${seed}/chk00001" erf.regrid_level_0_on_restart=1
    amr.max_grid_size=4 amr.max_grid_size_z=4 amr.no_box_split_dir=-1
    erf.substepping_type=None erf.fixed_dt=0.01 erf.fixed_mri_dt_ratio=2
    erf.implicit_thermal_diffusion=false erf.implicit_momentum_diffusion=false
    amrex.fpe_trap_invalid=0 erf.radiation_model=None
    erf.plot_int_1=-1 erf.plot2d_int_1=1 erf.plot2d_file_1=plt2d
    "erf.plot2d_vars_1=${vars}"
    erf.plot2d_int_2=1 erf.plot2d_file_2=stream2
    "erf.plot2d_vars_2=cloud_albedo_est cloud_optical_depth_sw"
    erf.station_names=CloudOptics erf.station_output_dir=stations
    erf.station_sampling_interval=1
    "erf.CloudOptics.field=cloud_albedo_est cloud_optical_depth_sw"
    erf.CloudOptics.x=8750 erf.CloudOptics.y=150 ${ARGN}
    WORKING_DIRECTORY "${dir}" OUTPUT_FILE "${dir}/run.log"
    ERROR_FILE "${dir}/error.log" RESULT_VARIABLE status)
  if(name MATCHES "^invalid_")
    if(status EQUAL 0)
      message(FATAL_ERROR "${name}: invalid parameter was accepted")
    endif()
    file(READ "${dir}/error.log" error)
    if(NOT error MATCHES "erf.plot2d.cloud_optics.*must be finite")
      message(FATAL_ERROR "${name}: expected a named cloud optics parameter error: ${error}")
    endif()
  elseif(NOT status EQUAL 0)
    message(FATAL_ERROR "${name} ERF failed (${status}); see ${dir}/run.log and error.log")
  endif()
endfunction()
run_case(both "cloud_albedo_est integrated_qr cloud_optical_depth_sw integrated_qc cloud_albedo_est")
run_case(tau "cloud_optical_depth_sw")
run_case(albedo "cloud_albedo_est")
run_case(invalid_coefficient "cloud_albedo_est" erf.plot2d.cloud_optics.kext_qc=-1)
run_case(invalid_scale "cloud_optical_depth_sw" erf.plot2d.cloud_optics.reflectance_scale=0)
foreach(name IN ITEMS both tau albedo)
  set(plot "${TEST_ROOT}/${name}/plt2d00002")
  file(STRINGS "${plot}/Header" header LIMIT_COUNT 7)
  file(READ "${plot}/2DMetadata.json" metadata)
  if(name STREQUAL "both")
    set(expected "integrated_qc;integrated_qr;cloud_optical_depth_sw;cloud_albedo_est")
    set(count 4)
  elseif(name STREQUAL "tau")
    set(expected cloud_optical_depth_sw)
    set(count 1)
  else()
    set(expected cloud_albedo_est)
    set(count 1)
  endif()
  list(GET header 1 actual_count)
  list(SUBLIST header 2 ${count} actual)
  if(NOT actual_count EQUAL count OR NOT actual STREQUAL expected)
    message(FATAL_ERROR "${name}: wrong Header component order/count: ${header}")
  endif()
  string(JSON metadata_count GET "${metadata}" n_variables)
  if(NOT metadata_count EQUAL count)
    message(FATAL_ERROR "${name}: wrong metadata variable count")
  endif()
  math(EXPR last "${count} - 1")
  foreach(comp RANGE 0 ${last})
    string(JSON var GET "${metadata}" variables ${comp} name)
    string(JSON index GET "${metadata}" variables ${comp} component_index)
    list(GET expected ${comp} expected_var)
    if(NOT var STREQUAL expected_var OR NOT index EQUAL comp)
      message(FATAL_ERROR "${name}: wrong metadata component ${comp}")
    endif()
    if(var MATCHES "^cloud_")
      string(JSON units GET "${metadata}" variables ${comp} units)
      string(JSON category GET "${metadata}" variables ${comp} category)
      string(JSON missing GET "${metadata}" variables ${comp} missing_value)
      if(NOT units STREQUAL "1" OR NOT category STREQUAL "Radiation" OR NOT missing EQUAL -999)
        message(FATAL_ERROR "${name}: incorrect cloud metadata")
      endif()
    endif()
  endforeach()
  file(READ "${TEST_ROOT}/${name}/stations/CloudOptics.dat" station)
  if(NOT station MATCHES "cloud_optical_depth_sw" OR NOT station MATCHES "cloud_albedo_est")
    message(FATAL_ERROR "${name}: missing station diagnostic columns")
  endif()
  file(STRINGS "${TEST_ROOT}/${name}/stations/CloudOptics.dat" rows REGEX "^[ ]*[0-9]")
  if(NOT rows OR rows MATCHES "-999|nan|inf")
    message(FATAL_ERROR "${name}: missing or invalid station values: ${rows}")
  endif()
  foreach(var IN LISTS expected)
    if(NOT metadata MATCHES "\"name\": \"${var}\"")
      message(FATAL_ERROR "${name}: metadata missing ${var}")
    endif()
  endforeach()
  if(name STREQUAL "both")
    execute_process(COMMAND "${CHECK_EXECUTABLE}" "${plot}" RESULT_VARIABLE status)
  else()
    execute_process(COMMAND "${CHECK_EXECUTABLE}" "${plot}"
      "${TEST_ROOT}/both/plt2d00002" RESULT_VARIABLE status)
  endif()
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${name}: numeric Plot2D check failed (${status})")
  endif()
  set(stream2 "${TEST_ROOT}/${name}/stream200002")
  execute_process(COMMAND "${CHECK_EXECUTABLE}" "${stream2}"
    "${TEST_ROOT}/both/plt2d00002" RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${name}: second stream numeric check failed (${status})")
  endif()
endforeach()
message(STATUS "CloudOptics radiation-off z-split writer, both streams and single requests passed")
