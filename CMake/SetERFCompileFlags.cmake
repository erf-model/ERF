function(set_erf_compile_flags target)
  # Logic for handling warnings
  if(ERF_ENABLE_ALL_WARNINGS)
    # GCC, Clang, and Intel seem to accept these
    list(APPEND ERF_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
    if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      # ifort doesn't like -Wall
      list(APPEND ERF_Fortran_FLAGS "-Wall")
    else()
      # Intel always reports some diagnostics we don't necessarily care about
      list(APPEND ERF_CXX_FLAGS "-diag-disable:11074,11076")
    endif()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
      # Avoid notes about -faligned-new with GCC > 7
      list(APPEND ERF_CXX_FLAGS "-faligned-new")
    endif()
  endif()

  # Add our extra flags according to language
  separate_arguments(ERF_CXX_FLAGS)
  target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${ERF_CXX_FLAGS}>)

  separate_arguments(ERF_Fortran_FLAGS)
  target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${ERF_Fortran_FLAGS}>)

  # Check if using both AMReX and Kokkos
  set(USING_KOKKOS FALSE)
  if(ERF_ENABLE_EKAT OR TARGET Kokkos::kokkoscore OR DEFINED Kokkos_DIR)
    set(USING_KOKKOS TRUE)
  endif()

  if(USING_KOKKOS AND TARGET amrex)
    check_gpu_architectures(${target})
  endif()

  # CUDA configuration
  if(ERF_ENABLE_CUDA)
    list(APPEND ERF_CUDA_FLAGS "--expt-relaxed-constexpr")
    list(APPEND ERF_CUDA_FLAGS "--expt-extended-lambda")
    list(APPEND ERF_CUDA_FLAGS "--Wno-deprecated-gpu-targets")
    list(APPEND ERF_CUDA_FLAGS "-m64")
    if(ENABLE_CUDA_FASTMATH)
      list(APPEND ERF_CUDA_FLAGS "--use_fast_math")
    endif()
    separate_arguments(ERF_CUDA_FLAGS)
    target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${ERF_CUDA_FLAGS}>)
    set_cuda_architectures(AMReX_CUDA_ARCH)
    set_target_properties(${target} PROPERTIES
      CUDA_ARCHITECTURES "${AMREX_CUDA_ARCHS}"
      LANGUAGE CUDA
      CUDA_SEPARABLE_COMPILATION ON
      CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )
  endif()

  # HIP configuration - handle duplicate flags
  if(ERF_ENABLE_HIP)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.21 AND DEFINED AMReX_AMD_ARCH)
      set_target_properties(${target} PROPERTIES HIP_ARCHITECTURES "${AMReX_AMD_ARCH}")
    elseif(USING_KOKKOS AND TARGET amrex)
      # Workaround for duplicate HIP flags from both libraries
      target_link_options(${target} PRIVATE "-Wl,--allow-multiple-definition")
    endif()
  endif()
endfunction()

# Check GPU architecture consistency between AMReX and Kokkos
function(check_gpu_architectures target)
  set(amrex_arch "")
  set(kokkos_arch "")

  # Get architectures based on backend
  if(ERF_ENABLE_CUDA)
    if(DEFINED AMReX_CUDA_ARCH)
      set(amrex_arch "${AMReX_CUDA_ARCH}")
    endif()
    # Check common Kokkos CUDA arch flags
    foreach(var Kokkos_ARCH_VOLTA70 Kokkos_ARCH_AMPERE80 Kokkos_ARCH_HOPPER90)
      if(${var})
        string(REGEX MATCH "([0-9]+)$" _ "${var}")
        set(kokkos_arch "sm_${CMAKE_MATCH_1}")
        break()
      endif()
    endforeach()

  elseif(ERF_ENABLE_HIP)
    if(DEFINED AMReX_AMD_ARCH)
      set(amrex_arch "${AMReX_AMD_ARCH}")
    endif()
    # Check Kokkos AMD arch or extract from link flags
    if(Kokkos_ARCH_AMD_GFX90A)
      set(kokkos_arch "gfx90a")
    elseif(TARGET Kokkos::kokkoscore)
      get_target_property(opts Kokkos::kokkoscore INTERFACE_LINK_OPTIONS)
      if(opts)
        string(REGEX MATCH "--offload-arch=([a-z0-9]+)" _ "${opts}")
        if(CMAKE_MATCH_1)
          set(kokkos_arch "${CMAKE_MATCH_1}")
        endif()
      endif()
    endif()

  elseif(ERF_ENABLE_SYCL)
    # SYCL doesn't typically have explicit arch conflicts
    return()
  endif()

  # Warn if architectures don't match
  if(amrex_arch AND kokkos_arch AND NOT "${amrex_arch}" STREQUAL "${kokkos_arch}")
    message(WARNING "GPU architecture mismatch for ${target}: AMReX=${amrex_arch}, Kokkos=${kokkos_arch}")
  endif()
endfunction()