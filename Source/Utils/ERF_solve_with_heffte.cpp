#include "ERF.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include "ERF_Utils.H"
#ifdef ERF_USE_HEFFTE
#include "heffte.h"
#endif

#ifdef ERF_USE_POISSON_SOLVE
#ifdef ERF_USE_HEFFTE

using namespace amrex;

void ERF::solve_with_heffte (int lev, MultiFab& rhs, MultiFab& phi,
                             Array<MultiFab,AMREX_SPACEDIM>& fluxes)
{
    BoxArray ba(rhs.boxArray());
    DistributionMapping dm(rhs.DistributionMap());

    // The heffte solve uses a solution array with no ghost cells
    MultiFab soln(ba,dm,1,0);

    // Determine the domain length in each direction
    Real L_x = geom[lev].ProbHi(0) - geom[lev].ProbLo(0);
    Real L_y = geom[lev].ProbHi(1) - geom[lev].ProbLo(1);
    Real L_z = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);

    Box domain = geom[lev].Domain();
    auto const& domlo = lbound(domain);
    auto const& domhi = ubound(domain);

    int n_cell_x = domain.length(0);
    int n_cell_y = domain.length(1);
    int n_cell_z = domain.length(2);

    auto dx    = geom[lev].CellSize();
    auto dxinv = geom[lev].InvCellSize();

    // Since there is 1 MPI rank per box, here each MPI rank obtains its local box and the associated boxid
    Box local_box;
    int local_boxid;
    {
        for (int i = 0; i < ba.size(); ++i) {
            Box b = ba[i];
            // each MPI rank has its own local_box Box and local_boxid ID
            if (ParallelDescriptor::MyProc() == dm[i]) {
                local_box = b;
                local_boxid = i;
            }
        }
    }

    // Now each MPI rank works on its own box
    // Ror real->complex fft's, the fft is stored in an (nx/2+1) x ny x nz dataset

    // start by coarsening each box by 2 in the x-direction
    Box c_local_box = amrex::coarsen(local_box, IntVect(AMREX_D_DECL(2,1,1)));

    // If the coarsened box's high-x index is even, we shrink the size in 1 in x
    // this avoids overlap between coarsened boxes
    if (c_local_box.bigEnd(0) * 2 == local_box.bigEnd(0)) {
        c_local_box.setBig(0,c_local_box.bigEnd(0)-1);
    }
    // For any boxes that touch the hi-x domain we increase the size of boxes by 1 in x
    // This makes the overall fft dataset have size (Nx/2+1 x Ny x Nz)
    if (local_box.bigEnd(0) == geom[lev].Domain().bigEnd(0)) {
        c_local_box.growHi(0,1);
    }

    // Each MPI rank gets storage for its piece of the fft
    BaseFab<GpuComplex<Real> > spectral_field(c_local_box, 1, The_Device_Arena());

    // Create real->complex fft objects with the appropriate backend and data about
    // the domain size and its local box size

    bool do_2d_solves = false;

    // ********************************************************************************************
    // ********************************************************************************************
    // ********************************************************************************************

    // ********************************************************************************************
    // NOTE: THIS IS A WIP - IT DOES NOT WORK YET
    // ********************************************************************************************
    if (do_2d_solves) {

#ifdef AMREX_USE_CUDA
        heffte::fft2d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
        heffte::fft2d_r2c<heffte::backend::rocfft> fft
#else
        heffte::fft2d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),0},
          {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,0}},
         {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),0},
          {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,0}},
         0, ParallelDescriptor::Communicator());

        using heffte_complex = typename heffte::fft_output<Real>::type;
        heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();

        // ********************************************************************************************

        for (int k = domlo.z; k <= domhi.z; k++) {
            int offset = k * (n_cell_x*n_cell_y);
            fft.forward(rhs[local_boxid].dataPtr(offset), spectral_data);
        }

        // ********************************************************************************************

        // Now we take the standard FFT and scale it by 1/k^2
        Array4< GpuComplex<Real> > spectral = spectral_field.array();

        ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real a = 2.*M_PI*i / L_x;
            Real b = 2.*M_PI*j / L_y;
            Real c = 2.*M_PI*k / L_z;

            // the values in the upper-half of the spectral array in y and z are here interpreted as negative wavenumbers
            if (j >= n_cell_y/2) b = 2.*M_PI*(n_cell_y-j) / L_y;
            if (k >= n_cell_z/2) c = 2.*M_PI*(n_cell_z-k) / L_z;

            Real k2 =  2.0*(std::cos(a*dx[0])-1.)*(dxinv[0]*dxinv[0]) +
                       2.0*(std::cos(b*dx[1])-1.)*(dxinv[1]*dxinv[1]) ;
            if (k2 != 0.) {
                spectral(i,j,k) /= k2;
            } else {
                spectral(i,j,k) *= 0.; // interpretation here is that the average value of the solution is zero
            }
        });

        // ********************************************************************************************

        for (int k = domlo.z; k <= domhi.z; k++) {
            int offset = k * (n_cell_x*n_cell_y);
            fft.backward(spectral_data, soln[local_boxid].dataPtr(offset));
        }

        // ********************************************************************************************

    } else {

#ifdef AMREX_USE_CUDA
        heffte::fft3d_r2c<heffte::backend::cufft> fft
#elif AMREX_USE_HIP
        heffte::fft3d_r2c<heffte::backend::rocfft> fft
#else
        heffte::fft3d_r2c<heffte::backend::fftw> fft
#endif
        ({{local_box.smallEnd(0),local_box.smallEnd(1),local_box.smallEnd(2)},
          {local_box.bigEnd(0)  ,local_box.bigEnd(1)  ,local_box.bigEnd(2)}},
         {{c_local_box.smallEnd(0),c_local_box.smallEnd(1),c_local_box.smallEnd(2)},
          {c_local_box.bigEnd(0)  ,c_local_box.bigEnd(1)  ,c_local_box.bigEnd(2)}},
         0, ParallelDescriptor::Communicator());

        using heffte_complex = typename heffte::fft_output<Real>::type;
        heffte_complex* spectral_data = (heffte_complex*) spectral_field.dataPtr();

        // ********************************************************************************************

        fft.forward(rhs[local_boxid].dataPtr(), spectral_data);

        // ********************************************************************************************

        // Now we take the standard FFT and scale it by 1/k^2
        Array4< GpuComplex<Real> > spectral = spectral_field.array();

        ParallelFor(c_local_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real a = 2.*M_PI*i / L_x;
            Real b = 2.*M_PI*j / L_y;
            Real c = 2.*M_PI*k / L_z;

            // the values in the upper-half of the spectral array in y and z are here interpreted as negative wavenumbers
            if (j >= n_cell_y/2) b = 2.*M_PI*(n_cell_y-j) / L_y;
            if (k >= n_cell_z/2) c = 2.*M_PI*(n_cell_z-k) / L_z;

            Real k2 =  2.0*(std::cos(a*dx[0])-1.)*(dxinv[0]*dxinv[0]) +
                       2.0*(std::cos(b*dx[1])-1.)*(dxinv[1]*dxinv[1]) +
                       2.0*(std::cos(c*dx[2])-1.)*(dxinv[2]*dxinv[2]);
            if (k2 != 0.) {
                spectral(i,j,k) /= k2;
            } else {
                spectral(i,j,k) *= 0.; // interpretation here is that the average value of the solution is zero
            }
        });

        // ********************************************************************************************

        fft.backward(spectral_data, soln[local_boxid].dataPtr());

        // ********************************************************************************************

    } // 3d solve

    // ********************************************************************************************
    // ********************************************************************************************
    // ********************************************************************************************

    // Scale by 1/npts (both forward and inverse need sqrt(npts) scaling so I am doing it all here)
    Real npts = static_cast<Real>(ba.numPts());
    soln.mult(1./npts);

    // ********************************************************************************************

    phi.copy(soln);
    phi.FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(soln, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const&  p_arr  = phi.array(mfi);

        Box const& xbx = mfi.nodaltilebox(0);
        const Real dx_inv = dxinv[0];
        Array4<Real      > const& fx_arr  = fluxes[0].array(mfi);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i-1,j,k)) * dx_inv;
        });

        Box const& ybx = mfi.nodaltilebox(1);
        const Real dy_inv = dxinv[1];
        Array4<Real      > const& fy_arr  = fluxes[1].array(mfi);
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dy_inv;
        });

        Box const& zbx = mfi.nodaltilebox(2);
        const Real dz_inv = dxinv[2];
        Array4<Real      > const& fz_arr  = fluxes[2].array(mfi);
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dz_inv;
        });
    } // mfi
}

#endif
#endif
