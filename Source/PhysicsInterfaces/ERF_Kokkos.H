#ifndef ERF_KOKKOS_H
#define ERF_KOKKOS_H
using KokkosHostDevice    = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
#if defined(ERF_USE_CUDA)
using KokkosDefaultDevice = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
#elif defined(ERF_USE_HIP)
using KokkosDefaultDevice = Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>;
#elif defined(ERF_USE_SYCL)
using KokkosDefaultDevice = Kokkos::Device<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace>;
#else
using KokkosDefaultDevice = KokkosHostDevice;
#endif
using RealT            = amrex::Real;
using layout_t         = Kokkos::LayoutRight;
using KokkosDefaultMem = Kokkos::DefaultExecutionSpace::memory_space;
using realHost1d_k     = Kokkos::View<RealT* , KokkosHostDevice>;
using realHost2d_k     = Kokkos::View<RealT**, layout_t, KokkosHostDevice>;
using real1d_k         = Kokkos::View<RealT*  , KokkosDefaultDevice>;
using real2d_k         = Kokkos::View<RealT** , layout_t, KokkosDefaultDevice>;
using real3d_k         = Kokkos::View<RealT***, layout_t, KokkosDefaultDevice>;
using int1d_k          = Kokkos::View<int*  , KokkosDefaultDevice>;
using int2d_k          = Kokkos::View<int** , layout_t, KokkosDefaultDevice>;
using int3d_k          = Kokkos::View<int***, layout_t, KokkosDefaultDevice>;
#endif
