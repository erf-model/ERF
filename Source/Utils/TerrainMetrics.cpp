#include <TerrainMetrics.H>
#include <AMReX_ParmParse.H>
#include <ERF_Constants.H>
#include <cmath>

using namespace amrex;

void
init_zlevels (amrex::Vector<amrex::Real>& zlevels_stag,
              const amrex::Geometry& geom,
              const amrex::Real grid_stretching_ratio,
              const amrex::Real zsurf,
              const amrex::Real dz0)
{
    auto dx = geom.CellSizeArray();
    const amrex::Box& domain = geom.Domain();
    int nz = domain.length(2)+1; // staggered
    zlevels_stag.resize(nz);

    if (grid_stretching_ratio == 0) {
        // This is the default for z_levels
        for (int k = 0; k < nz; k++)
        {
            zlevels_stag[k] = k * dx[2];
        }
    } else {
        // Create stretched grid based on initial dz and stretching ratio
        zlevels_stag[0] = zsurf;
        amrex::Real dz = dz0;
        for (int k = 1; k < nz; k++)
        {
            zlevels_stag[k] = zlevels_stag[k-1] + dz;
            dz *= grid_stretching_ratio;
        }
    }

    ParmParse pp("erf");
    int n_zlevels = pp.countval("terrain_z_levels");
    if (n_zlevels > 0)
    {
        if (n_zlevels != nz) {
            amrex::Print() << "You supplied " << n_zlevels << " staggered terrain_z_levels " << std::endl;
            amrex::Print() << "but n_cell+1 in the z-direction is " <<  nz << std::endl;
            amrex::Abort("You must specify a z_level for every value of k");
        }

        if (grid_stretching_ratio > 0) {
            amrex::Print() << "Note: Found terrain_z_levels, ignoring grid_stretching_ratio" << std::endl;
        }

        pp.queryarr("terrain_z_levels", zlevels_stag, 0, nz);
    }
}

/**
 * Computation of the terrain grid from BTF, STF, or Sullivan TF model
 *
 * NOTE: Multilevel is not yet working for either of these terrain-following coordinates,
 *       but (we think) the issue is deep in ERF and this code will work once the deeper
 *       problem is fixed. For now, make sure to run on a single level. -mmsanders
 */

void
init_terrain_grid (const Geometry& geom, MultiFab& z_phys_nd, amrex::Vector<Real> const& z_levels_h)
{
  auto dx = geom.CellSizeArray();
  auto ProbHiArr = geom.ProbHiArray();

  // z_nd is nodal in all directions
  const amrex::Box& domain = geom.Domain();
  int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
  int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
  int domlo_z = domain.smallEnd(2); int domhi_z = domain.bigEnd(2) + 1;
  int nz = domain.length(2)+1;

  // User-selected method from inputs file (BTF default)
  ParmParse pp("erf");
  int terrain_smoothing = 0;
  pp.query("terrain_smoothing", terrain_smoothing);

  // "custom" refers to terrain analytically specified, such as WitchOfAgnesi
  //std::string terrain_type = "custom";

   amrex::Gpu::DeviceVector<Real> z_levels_d;
   z_levels_d.resize(nz);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy_async(z_levels_d.data(), z_levels_h.data(), sizeof(Real)*nz);
#else
    std::memcpy(z_levels_d.data(), z_levels_h.data(), sizeof(Real)*nz);
#endif

  // Number of ghost cells
  int ngrow = z_phys_nd.nGrow();

  int imin = domlo_x; if (geom.isPeriodic(0)) imin -= z_phys_nd.nGrowVect()[0];
  int jmin = domlo_y; if (geom.isPeriodic(1)) jmin -= z_phys_nd.nGrowVect()[1];

  int imax = domhi_x; if (geom.isPeriodic(0)) imax += z_phys_nd.nGrowVect()[0];
  int jmax = domhi_y; if (geom.isPeriodic(1)) jmax += z_phys_nd.nGrowVect()[1];

  switch(terrain_smoothing) {
    case 0: // BTF Method
    {
      int k0    = 0;
      Real ztop = ProbHiArr[2];

      for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
      {
          // Grown box with corrected ghost cells at top
          Box gbx = mfi.growntilebox(ngrow);
          gbx.setRange(2,domlo_z,domhi_z+1);

          Array4<Real> const& z_arr = z_phys_nd.array(mfi);
          auto const&         z_lev = z_levels_d.data();

          ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
          {
              // Vertical grid stretching
              Real z =  z_lev[k];

              int ii = amrex::max(amrex::min(i,imax),imin);
              int jj = amrex::max(amrex::min(j,jmax),jmin);

              // Fill levels using BTF model from p2163 of Klemp2011
              z_arr(i,j,k) = z + (1. - (z/ztop)) * z_arr(ii,jj,k0);

              // Fill lateral boundaries and below the bottom surface
              if (k == k0)
                  z_arr(i,j,k0  ) = z_arr(ii,jj,k0);
              if (k == 1)
                  z_arr(i,j,k0-1) = 2.0*z_arr(ii,jj,k0) - z_arr(i,j,k);
          });
        }

        break;
    }

    case 1: // STF Method
    {
        // Get Multifab spanning domain with 1 level of ghost cells
        MultiFab h_mf(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, ngrow+1);
        MultiFab h_mf_old(z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, ngrow+1);

        // Save max height for smoothing
        Real max_h;

        // Crete 2D MF without allocation
        MultiFab mf2d;
        {
            BoxList bl2d = h_mf.boxArray().boxList();
            for (auto& b : bl2d) {
                b.setRange(2,0);
            }
            BoxArray ba2d(std::move(bl2d));
            mf2d = MultiFab(ba2d, h_mf.DistributionMap(), 1, ngrow, amrex::MFInfo().SetAlloc(false));
        }

        // Get MultiArray4s from the multifabs
        MultiArray4<Real> const& ma_h_s     = h_mf.arrays();
        MultiArray4<Real> const& ma_h_s_old = h_mf_old.arrays();
        MultiArray4<Real> const& ma_z_phys  = z_phys_nd.arrays();

        // Bottom boundary
        int k0 = domlo_z;

        // Get max value
        max_h = ParReduce(amrex::TypeList<amrex::ReduceOpMax>{}, amrex::TypeList<Real>{}, mf2d, amrex::IntVect(ngrow,ngrow,0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> amrex::GpuTuple<Real>
                {
                  // Get Array4s
                  const auto & h     = ma_h_s[box_no];
                  const auto & z_arr = ma_z_phys[box_no];

                  int ii = amrex::max(amrex::min(i,imax),imin);
                  int jj = amrex::max(amrex::min(j,jmax),jmin);

                  // Fill the lateral boundaries
                  z_arr(i,j,k0) = z_arr(ii,jj,k0);

                  // Populate h with terrain
                  h(i,j,k0) = z_arr(i,j,k0);

                  // Return height for max
                  return { z_arr(i,j,k0) };
                });

        // Fill ghost cells (neglects domain boundary if not periodic)
        h_mf.FillBoundary(geom.periodicity());

        // Make h_mf copy for old values
        MultiFab::Copy(h_mf_old, h_mf,0,0,1,h_mf_old.nGrow());

        // Populate h_mf at k>0 with h_s, solving in ordered 2D slices
        for (int k = domlo_z+1; k <= domhi_z; k++) // skip terrain level
        {
            auto const& z_lev_h = z_levels_h.data();

            Real zz       = z_lev_h[k];
            Real zz_minus = z_lev_h[k-1];

            Real gamma_m = 0.5; // min allowed fractional grid spacing
            Real z_H     = 2.44/(1-gamma_m);

            Real A;
            Real foo = cos((PI/2)*(zz/z_H));
            if(zz < z_H) { A = foo*foo*foo*foo*foo*foo; } // A controls rate of return to atm
            else         { A = 0; }
            Real foo_minus = cos((PI/2)*(zz_minus/z_H));
            Real A_minus;
            if(zz_minus < z_H) { A_minus = foo_minus*foo_minus*foo_minus*foo_minus*foo_minus*foo_minus; } // A controls rate of return to atm
            else               { A_minus = 0; }

            unsigned maxIter = 50; // M_k in paper
            unsigned iter    = 0;
            Real threshold   = gamma_m;
            Real diff        = 1.e20;
            while (iter < maxIter && diff > threshold)
            {

                diff = ParReduce(amrex::TypeList<amrex::ReduceOpMin>{}, amrex::TypeList<Real>{}, mf2d, amrex::IntVect(ngrow,ngrow,0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> amrex::GpuTuple<Real>
                {
                    const auto & h_s     = ma_h_s[box_no];
                    const auto & h_s_old = ma_h_s_old[box_no];

                    Real h_m    = max_h; //high point of hill
                    Real beta_k = 0.2*std::min(zz/(2*h_m),1.0); //smoothing coefficient

                    // Clip indices for ghost-cells
                    int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                    int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                    if (iter == 0) {
                        h_s(i,j,k) = h_s_old(i,j,k-1) + beta_k*(h_s_old(ii+1,jj  ,k-1)
                                                              + h_s_old(ii-1,jj  ,k-1)
                                                              + h_s_old(ii  ,jj+1,k-1)
                                                              + h_s_old(ii  ,jj-1,k-1) - 4*h_s_old(ii,jj,k-1));
                    }
                    else {
                        h_s(i,j,k) = h_s_old(i,j,k  ) + beta_k*(h_s_old(ii+1,jj  ,k  )
                                                              + h_s_old(ii-1,jj  ,k  )
                                                              + h_s_old(ii  ,jj+1,k  )
                                                              + h_s_old(ii  ,jj-1,k  ) - 4*h_s_old(ii,jj,k  ));
                    }

                    return { (zz + A * h_s(i,j,k) - (zz_minus + A_minus * h_s(i,j,k-1))) / (zz - zz_minus) };

                }); //ParReduce

                MultiFab::Copy(h_mf_old, h_mf,0,0,1,h_mf_old.nGrow());

                amrex::ParallelDescriptor::ReduceRealMin(diff);

                iter++;

                //fill ghost points
                h_mf_old.FillBoundary(geom.periodicity());

            } //while

            auto const& z_lev_d = z_levels_d.data();

            //Populate z_phys_nd by solving z_arr(i,j,k) = z + A*h_s(i,j,k)
            for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                // Grown box with no z range
                Box xybx = mfi.growntilebox(ngrow);
                xybx.setRange(2,0);

                Array4<Real> const& h_s   = h_mf_old.array(mfi);
                Array4<Real> const& z_arr = z_phys_nd.array(mfi);

                ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

                    // Location of nodes
                    Real z = z_lev_d[k];

                    // STF model from p2164 of Klemp2011
                    z_arr(i,j,k) = z + A*h_s(i,j,k);

                    // Fill below the bottom surface
                    if (k == 1)
                        z_arr(i,j,k0-1) = 2.0*z_arr(i,j,k0) - z_arr(i,j,k);
                });
            }//mfi
        } //k

        amrex::Gpu::streamSynchronize();

        break;
    }

    case 2: // Sullivan TF Method
    {
        int k0    = 0;
        Real ztop = ProbHiArr[2];

        for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Grown box with corrected ghost cells at top
            Box gbx = mfi.growntilebox(ngrow);
            gbx.setRange(2,domlo_z,domhi_z+1);

            Array4<Real> const& z_arr = z_phys_nd.array(mfi);
            auto const&         z_lev = z_levels_d.data();

            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Vertical grid stretching
                Real z =  z_lev[k];

                int ii = amrex::max(amrex::min(i,imax),imin);
                int jj = amrex::max(amrex::min(j,jmax),jmin);

                // Fill levels using model from Sullivan et. al. 2014
                int omega = 3; //Used to adjust how rapidly grid lines level out. omega=1 is BTF!
                z_arr(i,j,k) = z + (std::pow((1. - (z/ztop)),omega) * z_arr(ii,jj,k0));

                // Fill lateral boundaries and below the bottom surface
                if (k == k0)
                    z_arr(i,j,k0  ) = z_arr(ii,jj,k0);
                if (k == 1)
                    z_arr(i,j,k0-1) = 2.0*z_arr(ii,jj,k0) - z_arr(i,j,k);
            });
        }

        break;
    }
    case 3: // Debugging Test Method -- applies Sullivan TF starting at k = 1 so that domain does not change size
    {
        int k0    = 0;
        Real ztop = ProbHiArr[2];

        for ( amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Grown box with corrected ghost cells at top
            Box gbx = mfi.growntilebox(ngrow);
            gbx.setRange(2,domlo_z,domhi_z+1);

            Array4<Real> const& z_arr = z_phys_nd.array(mfi);
            auto const&         z_lev = z_levels_d.data();

            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Vertical grid stretching
                Real z =  z_lev[k];

                int ii = amrex::max(amrex::min(i,imax),imin);
                int jj = amrex::max(amrex::min(j,jmax),jmin);

                // Fill values outside the lateral boundaries and below the bottom surface (necessary if init_type = "real")
                if (k == k0+1) {
                    z_arr(i,j,k) = z + z_arr(ii,jj,k0);
                } else {
                    // Fill levels using model from Sullivan et. al. 2014
                    int omega = 3; //Used to adjust how rapidly grid lines level out. omega=1 is BTF!
                    z_arr(i,j,k) = z + (std::pow((1. - (z/ztop)),omega) * z_arr(ii,jj,k0));
                }
            });
            gbx.setBig(2,0);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                z_arr(i,j,k  ) = 0.0;
                z_arr(i,j,k-1) = -z_arr(i,j,k+1);
            });
        }

        break;
    }
  } //switch

  // Fill ghost layers and corners (including periodic)
  z_phys_nd.FillBoundary(geom.periodicity());


  //********************************************************************************
  // Populate domain boundary cells in z-direction
  //********************************************************************************
  for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
      // Grown box with no z range
      Box xybx = mfi.growntilebox(ngrow);
      xybx.setRange(2,0);

      Array4<Real> const& z_arr = z_phys_nd.array(mfi);

      // Extrapolate top layer
      ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
          z_arr(i,j,domhi_z+1) = 2.0*z_arr(i,j,domhi_z) - z_arr(i,j,domhi_z-1);
      });
  }

  /*
  // Debug
  for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
       Box gbx = mfi.growntilebox(ngrow);
       gbx.setRange(2,domlo_z,domhi_z+1);

       Array4<Real> z_arr = z_phys_nd.array(mfi);

       int rank = 0;
       amrex::Print(rank) << "Debugging init_terrain_grid" << "\n";
       amrex::Print(rank) << gbx << "\n";

       ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
           amrex::Print(rank) << amrex::IntVect(i,j,k) << "\n";
           amrex::Print(rank) << z_arr(i,j,k) << "\n";
           amrex::Print(rank) << "\n";
       });
       amrex::Print(rank) << "Cleared..." << "\n";
   }
  */
}

/**
 * Computation of detJ at cell-center
 */
void
make_J(const amrex::Geometry& geom,
       amrex::MultiFab& z_phys_nd,
       amrex::MultiFab& detJ_cc)
{
    const auto *dx = geom.CellSize();
    amrex::Real dzInv = 1.0/dx[2];

    // Domain valid box (z_nd is nodal)
    const amrex::Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // Number of ghost cells
    int ngrow = detJ_cc.nGrow();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(detJ_cc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        amrex::Box gbx = mfi.growntilebox(ngrow);
        gbx.setSmall(2,domlo_z);
        amrex::Array4<amrex::Real const> z_nd = z_phys_nd.const_array(mfi);
        amrex::Array4<amrex::Real      > detJ = detJ_cc.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               detJ(i, j, k) = .25 * dzInv * (
                       z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
                      -z_nd(i,j,k  ) - z_nd(i+1,j,k  ) - z_nd(i,j+1,k  ) - z_nd(i+1,j+1,k  ) );
       });
    }
    detJ_cc.FillBoundary(geom.periodicity());
}

/**
 * Computation of z_phys at cell-center
 */
void
make_zcc(const amrex::Geometry& geom,
         amrex::MultiFab& z_phys_nd,
         amrex::MultiFab& z_phys_cc)
{
    // Domain valid box (z_nd is nodal)
    const amrex::Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // Number of ghost cells
    int ngrow = z_phys_cc.nGrow();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(z_phys_cc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        amrex::Box gbx = mfi.growntilebox(ngrow);
        gbx.setSmall(2,domlo_z);
        amrex::Array4<amrex::Real const> z_nd = z_phys_nd.const_array(mfi);
        amrex::Array4<amrex::Real      > z_cc = z_phys_cc.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
           z_cc(i, j, k) = .125 * ( z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  )
                                   +z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1) );
       });
    }
    z_phys_cc.FillBoundary(geom.periodicity());
}
