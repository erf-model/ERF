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
  # HIP configuration - deduplicate flags at the target level
  if(ERF_ENABLE_HIP)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.21)
      # Extract architecture from AMReX, CMake, or environment
      set(hip_arch "")
      if(DEFINED AMReX_AMD_ARCH)
        set(hip_arch "${AMReX_AMD_ARCH}")
      elseif(DEFINED CMAKE_HIP_ARCHITECTURES)
        set(hip_arch "${CMAKE_HIP_ARCHITECTURES}")
      elseif(DEFINED ENV{ROCM_GPU})
        set(hip_arch "$ENV{ROCM_GPU}")
      elseif(DEFINED ENV{HIPARCHS})
        set(hip_arch "$ENV{HIPARCHS}")
      elseif(DEFINED ENV{HIP_ARCH})
        set(hip_arch "$ENV{HIP_ARCH}")
      else()
        set(hip_arch "gfx90a")  # Default for Frontier
      endif()
      set_target_properties(${target} PROPERTIES HIP_ARCHITECTURES "${hip_arch}")
      message(VERBOSE "Set HIP_ARCHITECTURES=${hip_arch} for ${target}")
      # If both AMReX and Kokkos are present, remove duplicate link options
      if(TARGET amrex AND (ERF_ENABLE_EKAT OR TARGET Kokkos::kokkoscore))
        # Get current target's link options
        get_target_property(target_link_opts ${target} LINK_OPTIONS)
        if(target_link_opts)
          # Remove duplicate HIP flags
          list(REMOVE_DUPLICATES target_link_opts)
          # Also remove any extra --hip-link and --offload-arch flags
          set(filtered_opts "")
          set(found_hip_link FALSE)
          set(found_offload_arch FALSE)
          foreach(opt ${target_link_opts})
            if(opt MATCHES "^--hip-link$")
              if(NOT found_hip_link)
                list(APPEND filtered_opts ${opt})
                set(found_hip_link TRUE)
              endif()
            elseif(opt MATCHES "^--offload-arch=")
              if(NOT found_offload_arch)
                list(APPEND filtered_opts ${opt})
                set(found_offload_arch TRUE)
              endif()
            else()
              list(APPEND filtered_opts ${opt})
            endif()
          endforeach()
          set_target_properties(${target} PROPERTIES LINK_OPTIONS "${filtered_opts}")
        endif()
      endif()
    else()
      # For older CMake, just add workaround flag
      if(ERF_ENABLE_EKAT OR TARGET Kokkos::kokkoscore)
        target_link_options(${target} PRIVATE "-Wl,--allow-multiple-definition")
      endif()
    endif()
  endif()
endfunction()
