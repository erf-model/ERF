cmake_minimum_required(VERSION 3.20)
# One direct radiation-off Kessler writer run; mathematical edge cases live in
# ERF_GTestCloudOptics.cpp. Keep the canonical Bubble input deck unchanged.
string(RANDOM LENGTH 8 ALPHABET 0123456789abcdef run_id)
set(dir "${TEST_ROOT}/run_${run_id}")
file(MAKE_DIRECTORY "${dir}")
execute_process(COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}"
  max_step=2 "amr.n_cell=8 4 8" amr.max_level=0 amr.blocking_factor=2
  amr.max_grid_size=8 amr.max_grid_size_z=8 amr.no_box_split_dir=2
  erf.fixed_dt=0.01 erf.fixed_mri_dt_ratio=2
  erf.moisture_model=Kessler erf.radiation_model=None
  erf.plot_int_1=-1 erf.plot_int_2=-1 erf.check_int=-1 erf.check_per=-1
  erf.plot2d_int_1=2 erf.plot2d_file_1=stream1
  "erf.plot2d_vars_1=cloud_albedo_est integrated_qr cloud_optical_depth_sw integrated_qc"
  erf.plot2d_int_2=2 erf.plot2d_file_2=stream2
  "erf.plot2d_vars_2=cloud_albedo_est"
  erf.plot2d.cloud_optics.kext_qc=50 erf.plot2d.cloud_optics.kext_qr=2
  erf.plot2d.cloud_optics.reflectance_scale=3
  erf.station_names=CloudOptics erf.station_output_dir=stations
  erf.station_sampling_interval=1
  "erf.CloudOptics.field=cloud_albedo_est cloud_optical_depth_sw"
  erf.CloudOptics.x=8250 erf.CloudOptics.y=175
  WORKING_DIRECTORY "${dir}" OUTPUT_FILE "${dir}/run.log"
  ERROR_FILE "${dir}/error.log" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Kessler writer failed (${status}); inspect ${dir}/run.log and error.log")
endif()
file(STRINGS "${dir}/run.log" settings REGEX "^Plot2D cloud optics:")
list(LENGTH settings count)
if(NOT count EQUAL 1 OR NOT settings MATCHES "kext_qc=50 kext_qi=" OR
   NOT settings MATCHES "kext_qr=2 kext_qs=" OR
   NOT settings MATCHES "reflectance_scale=3$")
  message(FATAL_ERROR "Incorrect resolved optical settings: ${settings}")
endif()
file(WRITE "${dir}/cloud_optics_resolved_settings.txt" "${settings}\n")
foreach(stream IN ITEMS stream1 stream2)
  file(GLOB headers "${dir}/${stream}*/Header")
  list(SORT headers)
  list(LENGTH headers count)
  if(count EQUAL 0)
    message(FATAL_ERROR "No ${stream} output in ${dir}")
  endif()
  list(GET headers -1 header_file)
  get_filename_component(${stream} "${header_file}" DIRECTORY)
  file(STRINGS "${header_file}" header LIMIT_COUNT 7)
  file(READ "${${stream}}/2DMetadata.json" metadata)
  if(stream STREQUAL "stream1")
    set(expected integrated_qc integrated_qr cloud_optical_depth_sw cloud_albedo_est)
  else()
    set(expected cloud_albedo_est)
  endif()
  list(LENGTH expected count)
  list(GET header 1 actual_count)
  list(SUBLIST header 2 ${count} names)
  string(JSON metadata_count GET "${metadata}" n_variables)
  if(NOT actual_count EQUAL count OR NOT names STREQUAL expected OR
     NOT metadata_count EQUAL count)
    message(FATAL_ERROR "${stream}: wrong Header/metadata count or order: ${header}")
  endif()
  math(EXPR last "${count} - 1")
  foreach(comp RANGE 0 ${last})
    list(GET expected ${comp} expected_name)
    string(JSON name GET "${metadata}" variables ${comp} name)
    string(JSON index GET "${metadata}" variables ${comp} component_index)
    if(NOT name STREQUAL expected_name OR NOT index EQUAL comp)
      message(FATAL_ERROR "${stream}: incorrect metadata component ${comp}")
    endif()
    if(name MATCHES "^cloud_")
      string(JSON units GET "${metadata}" variables ${comp} units)
      string(JSON category GET "${metadata}" variables ${comp} category)
      string(JSON missing GET "${metadata}" variables ${comp} missing_value)
      if(NOT units STREQUAL "1" OR NOT category STREQUAL "Radiation" OR NOT missing EQUAL -999)
        message(FATAL_ERROR "${stream}: incorrect optical metadata")
      endif()
    endif()
  endforeach()
endforeach()
file(READ "${dir}/stations/CloudOptics.dat" station)
if(NOT station MATCHES "cloud_albedo_est.*cloud_optical_depth_sw")
  message(FATAL_ERROR "Missing or reordered station fields")
endif()
message(STATUS "Checking ${stream1}, ${stream2}; logs retained in ${dir}")
execute_process(COMMAND "${CHECK_EXECUTABLE}" "${stream1}" "${stream2}"
  "${dir}/stations/CloudOptics.dat" RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Cloud optics writer numerical check failed (${status})")
endif()
