# Actual Morrison writer with positive, distinct qc/qi/qr/qs/qg checkpoint
# profiles, followed by a genuine two-level partial-column example.
file(MAKE_DIRECTORY "${TEST_ROOT}/seed")
set(seed "${TEST_ROOT}/seed")
set(common "amr.n_cell=8 4 8" amr.max_level=1 amr.ref_ratio=2
  amr.blocking_factor=2 erf.fixed_dt=0.01 erf.fixed_mri_dt_ratio=2
  erf.moisture_model=Morrison erf.radiation_model=None amrex.fpe_trap_invalid=0
  erf.refinement_indicators=box1 erf.box1.max_level=1 "erf.box1.in_box_lo=0 0 0")
execute_process(COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${common}
  max_step=1 amr.max_grid_size=16 amr.max_grid_size_z=16 amr.no_box_split_dir=2
  "erf.box1.in_box_hi=20000 400 10000" erf.plot_int_1=-1 erf.plot2d_int_1=-1
  erf.check_int=1 erf.check_file=chk
  WORKING_DIRECTORY "${seed}" OUTPUT_FILE "${seed}/run.log"
  ERROR_FILE "${seed}/error.log" RESULT_VARIABLE status)
if(NOT status EQUAL 0 OR NOT EXISTS "${seed}/chk00001/Level_1/Cell_H")
  message(FATAL_ERROR "Two-level Morrison seed failed (${status}); see ${seed}")
endif()
file(MAKE_DIRECTORY "${TEST_ROOT}/fixture")
file(COPY "${seed}/chk00001" DESTINATION "${TEST_ROOT}/fixture")
set(checkpoint "${TEST_ROOT}/fixture/chk00001")
execute_process(COMMAND "${FIXTURE_EXECUTABLE}" "${checkpoint}"
  OUTPUT_FILE "${TEST_ROOT}/fixture/seed.log" ERROR_FILE "${TEST_ROOT}/fixture/error.log"
  RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Synthetic five-species checkpoint fixture failed (${status})")
endif()
foreach(case IN ITEMS mixed partial)
  set(dir "${TEST_ROOT}/${case}")
  file(MAKE_DIRECTORY "${dir}")
  set(options max_step=1 amr.max_grid_size=16 amr.max_grid_size_z=16
    amr.no_box_split_dir=2 "erf.box1.in_box_hi=20000 400 10000")
  if(case STREQUAL "partial")
    # Regrid an already initialized full-height checkpoint: the existing Bubble
    # initializer itself requires full-height boxes and is left untouched.
    set(options max_step=2 amr.max_grid_size=4 amr.max_grid_size_z=4
      amr.no_box_split_dir=-1 amr.n_error_buf=0 erf.regrid_int=1
      "erf.box1.in_box_hi=20000 400 2500" erf.substepping_type=None
      erf.implicit_thermal_diffusion=false erf.implicit_momentum_diffusion=false
      erf.station_names=CloudOptics erf.station_output_dir=stations
      erf.station_sampling_interval=1 erf.CloudOptics.x=8750 erf.CloudOptics.y=150
      "erf.CloudOptics.field=cloud_albedo_est cloud_optical_depth_sw")
  endif()
  execute_process(COMMAND "${ERF_EXECUTABLE}" "${INPUT_FILE}" ${common} ${options}
    "erf.restart=${checkpoint}" erf.plot_int_1=1 erf.plot_file_1=volume
    "erf.plot_vars_1=rhoQ2 rhoQ3 rhoQ4 rhoQ5 rhoQ6"
    erf.plot2d_int_1=1 "erf.plot2d_file_1=${case}"
    "erf.plot2d_vars_1=cloud_albedo_est integrated_qg integrated_qs cloud_optical_depth_sw integrated_qr integrated_qi integrated_qc"
    WORKING_DIRECTORY "${dir}" OUTPUT_FILE "${dir}/run.log"
    ERROR_FILE "${dir}/error.log" RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${case} Morrison writer failed (${status}); see ${dir}")
  endif()
  file(STRINGS "${dir}/run.log" settings REGEX "^Plot2D cloud optics:")
  list(LENGTH settings count)
  if(NOT count EQUAL 1)
    message(FATAL_ERROR "${case}: expected one resolved parameter record")
  endif()
  foreach(key IN ITEMS kext_qc kext_qi kext_qr kext_qs kext_qg reflectance_scale)
    if(NOT settings MATCHES "${key}=[0-9.eE+-]+")
      message(FATAL_ERROR "${case}: missing resolved ${key}")
    endif()
  endforeach()
  file(WRITE "${dir}/cloud_optics_resolved_settings.txt" "${settings}\n")
  set(step 00001)
  if(case STREQUAL "partial")
    set(step 00002)
  endif()
  execute_process(COMMAND "${CHECK_EXECUTABLE}" "--${case}"
    "${dir}/${case}${step}" "${dir}/volume${step}" RESULT_VARIABLE status)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${case}: independent native 3D/2D check failed (${status})")
  endif()
  if(case STREQUAL "partial")
    execute_process(COMMAND "${CHECK_EXECUTABLE}" --station "${dir}/${case}${step}"
      "${dir}/stations/CloudOptics.dat" 8750 150 RESULT_VARIABLE status)
    if(NOT status EQUAL 0)
      message(FATAL_ERROR "partial fine-level station does not match its Plot2D column")
    endif()
  endif()
endforeach()
message(STATUS "Five-species Morrison writer and two-level PARTIAL (not composite) proxies passed")
