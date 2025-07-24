#ifndef ERF_KOKKOS_RADIATION_H
#define ERF_KOKKOS_RADIATION_H
using RadHostDevice    = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
#if defined(ERF_USE_CUDA)
using RadDefaultDevice = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
#elif defined(ERF_USE_HIP)
using RadDefaultDevice = Kokkos::Device<Kokkos::HIP, Kokkos::HIPSpace>;
#elif defined(ERF_USE_SYCL)
using RadDefaultDevice = Kokkos::Device<Kokkos::Experimental::SYCL, Kokkos::Experimental::SYCLDeviceUSMSpace>;
#else
using RadDefaultDevice = RadHostDevice;
#endif
using layout_t         = Kokkos::LayoutRight;
using RealT            = amrex::Real;
using realHost1d_k     = Kokkos::View<RealT* , RadHostDevice>;
using realHost2d_k     = Kokkos::View<RealT**, layout_t, RadHostDevice>;
using real1d_k         = Kokkos::View<RealT*  , RadDefaultDevice>;
using real2d_k         = Kokkos::View<RealT** , layout_t, RadDefaultDevice>;
using real3d_k         = Kokkos::View<RealT***, layout_t, RadDefaultDevice>;
using int1d_k          = Kokkos::View<int*  , RadDefaultDevice>;
using int2d_k          = Kokkos::View<int** , layout_t, RadDefaultDevice>;
using int3d_k          = Kokkos::View<int***, layout_t, RadDefaultDevice>;
#endif
