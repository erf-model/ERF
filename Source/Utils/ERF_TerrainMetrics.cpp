#include <ERF_TerrainMetrics.H>
#include <ERF_Utils.H>
#include <AMReX_ParmParse.H>
#include <ERF_Constants.H>
#include <cmath>

using namespace amrex;

void
init_zlevels (Vector<Real>& zlevels_stag,
              const Geometry& geom,
              const Real grid_stretching_ratio,
              const Real zsurf,
              const Real dz0)
{
    auto dx = geom.CellSizeArray();
    const Box& domain = geom.Domain();
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
        Real dz = dz0;
        Print() << "Stretched grid levels: " << zsurf;
        for (int k = 1; k < nz; k++)
        {
            zlevels_stag[k] = zlevels_stag[k-1] + dz;
            Print() << " " << zlevels_stag[k];
            dz *= grid_stretching_ratio;
        }
        Print() << std::endl;
    }

    // Try reading in terrain_z_levels, which allows arbitrarily spaced grid
    // levels to be specified and will take precedence over the
    // grid_stretching_ratio parameter
    ParmParse pp("erf");
    int n_zlevels = pp.countval("terrain_z_levels");
    if (n_zlevels > 0)
    {
        if (n_zlevels != nz) {
            Print() << "You supplied " << n_zlevels << " staggered terrain_z_levels " << std::endl;
            Print() << "but n_cell+1 in the z-direction is " <<  nz << std::endl;
            Abort("You must specify a z_level for every value of k");
        }

        if (grid_stretching_ratio > 0) {
            Print() << "Note: Found terrain_z_levels, ignoring grid_stretching_ratio" << std::endl;
        }

        pp.queryarr("terrain_z_levels", zlevels_stag, 0, nz);

        // These levels should range from 0 at the surface to the height of the
        // top of model domain (see the coordinate surface height, zeta, in
        // Klemp 2011)
        AMREX_ALWAYS_ASSERT(zlevels_stag[0] == 0);
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
init_terrain_grid (int lev, const Geometry& geom, MultiFab& z_phys_nd,
                   Vector<Real> const& z_levels_h,
                   GpuArray<ERF_BC, AMREX_SPACEDIM*2>& phys_bc_type)
{
    const Box& domain = geom.Domain();

    int domlo_z = domain.smallEnd(2);
    int domhi_z = domain.bigEnd(2) + 1;

    // ****************************************************************************

    if (lev == 0) {
        BoxArray ba(z_phys_nd.boxArray());
        bool all_boxes_touch_bottom = true;
        for (int i = 0; i < ba.size(); i++) {
            if (ba[i].smallEnd(2) != domlo_z) {
                all_boxes_touch_bottom = false;
            }
        }

        if (all_boxes_touch_bottom) {
            init_which_terrain_grid(lev, geom, z_phys_nd, z_levels_h);
        } else {

            BoxArray ba_new(domain);
            ChopGrids2D(ba_new, domain, ParallelDescriptor::NProcs());

            DistributionMapping dm_new(ba_new);
            ba_new.surroundingNodes();

            MultiFab z_phys_nd_new(ba_new, dm_new, 1, z_phys_nd.nGrowVect());

            z_phys_nd_new.ParallelCopy(z_phys_nd,0,0,1,z_phys_nd.nGrowVect(),z_phys_nd.nGrowVect());

            init_which_terrain_grid(lev, geom, z_phys_nd_new, z_levels_h);

            z_phys_nd.ParallelCopy(z_phys_nd_new,0,0,1,z_phys_nd.nGrowVect(),z_phys_nd.nGrowVect());
        }
    } else {
        init_which_terrain_grid(lev, geom, z_phys_nd, z_levels_h);
    }

    // Fill ghost layers and corners (including periodic)
    z_phys_nd.FillBoundary(geom.periodicity());

    if (phys_bc_type[Orientation(0,Orientation::low )] == ERF_BC::symmetry ||
        phys_bc_type[Orientation(0,Orientation::high)] == ERF_BC::symmetry ||
        phys_bc_type[Orientation(1,Orientation::low )] == ERF_BC::symmetry ||
        phys_bc_type[Orientation(1,Orientation::high)] == ERF_BC::symmetry) {

        const auto& dom_lo = lbound(convert(geom.Domain(),IntVect(1,1,1)));
        const auto& dom_hi = ubound(convert(geom.Domain(),IntVect(1,1,1)));

        for (MFIter mfi(z_phys_nd,true); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox();
            const Array4< Real> z_nd_arr = z_phys_nd.array(mfi);
            if (phys_bc_type[Orientation(0,Orientation::low)] == ERF_BC::symmetry && bx.smallEnd(0) == dom_lo.x) {
                ParallelFor(makeSlab(bx,0,1), [=] AMREX_GPU_DEVICE (int , int j, int k)
                {
                    z_nd_arr(dom_lo.x-1,j,k) = z_nd_arr(dom_lo.x+1,j,k);
                });
            }
            if (phys_bc_type[Orientation(0,Orientation::high)] == ERF_BC::symmetry && bx.bigEnd(0) == dom_hi.x) {
                ParallelFor(makeSlab(bx,0,1), [=] AMREX_GPU_DEVICE (int , int j, int k)
                {
                    z_nd_arr(dom_hi.x+1,j,k) = z_nd_arr(dom_hi.x-1,j,k);
                });
            }
            if (phys_bc_type[Orientation(1,Orientation::low)] == ERF_BC::symmetry && bx.smallEnd(1) == dom_lo.y) {
                ParallelFor(makeSlab(bx,1,1), [=] AMREX_GPU_DEVICE (int i, int  , int k)
                {
                    z_nd_arr(i,dom_lo.y-1,k) = z_nd_arr(i,dom_lo.y+1,k);
                });
            }
            if (phys_bc_type[Orientation(1,Orientation::high)] == ERF_BC::symmetry && bx.bigEnd(1) == dom_hi.y) {
                ParallelFor(makeSlab(bx,1,1), [=] AMREX_GPU_DEVICE (int i, int  , int k)
                {
                    z_nd_arr(i,dom_hi.y+1,k) = z_nd_arr(i,dom_hi.y-1,k);
                });
            }
        }
    }

    //********************************************************************************
    // Populate domain boundary cells in z-direction
    //********************************************************************************
    int ngrow = z_phys_nd.nGrow();

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Only set values above top of domain if this box reaches that far
        Box nd_bx = mfi.tilebox();

        // Note that domhi_z is already nodal in the z-direction
        if (nd_bx.bigEnd(2) >= domhi_z) {
            // Grown box with no z range
            Box xybx = mfi.growntilebox(ngrow);
            xybx.setRange(2,0);
            Array4<Real> const& z_arr = z_phys_nd.array(mfi);

            // Extrapolate top layer
            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int ) {
                z_arr(i,j,domhi_z+1) = 2.0*z_arr(i,j,domhi_z) - z_arr(i,j,domhi_z-1);
            });
        }
    }

    /*
    // Debug
    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ngrow);
        gbx.setRange(2,domlo_z,domhi_z+1);

        Array4<Real> z_arr = z_phys_nd.array(mfi);

        int rank = 0;
        Print(rank) << "Debugging init_terrain_grid" << "\n";
        Print(rank) << gbx << "\n";

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Print(rank) << IntVect(i,j,k) << "\n";
            Print(rank) << z_arr(i,j,k) << "\n";
            Print(rank) << "\n";
        });
        Print(rank) << "Cleared..." << "\n";
    }
   */
} // init_terrain_grid

void
init_which_terrain_grid (int lev, const Geometry& geom, MultiFab& z_phys_nd,
                         Vector<Real> const& z_levels_h)
{
    // User-selected method from inputs file (BTF default)
    ParmParse pp("erf");
    int terrain_smoothing = 0;
    pp.query("terrain_smoothing", terrain_smoothing);

    if (lev > 0 && terrain_smoothing != 0) {
        Abort("Must use terrain_smoothing = 0 when doing multilevel");
    }

    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();
    IntVect ngrowVect = z_phys_nd.nGrowVect();

    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2); int domhi_z = domain.bigEnd(2) + 1;

    int imin = domlo_x; if (geom.isPeriodic(0)) imin -= z_phys_nd.nGrowVect()[0];
    int jmin = domlo_y; if (geom.isPeriodic(1)) jmin -= z_phys_nd.nGrowVect()[1];

    int imax = domhi_x; if (geom.isPeriodic(0)) imax += z_phys_nd.nGrowVect()[0];
    int jmax = domhi_y; if (geom.isPeriodic(1)) jmax += z_phys_nd.nGrowVect()[1];

    int nz = domain.length(2)+1; // staggered
    Real z_top = z_levels_h[nz-1];

    Gpu::DeviceVector<Real> z_levels_d;
    z_levels_d.resize(nz);
    Gpu::copy(Gpu::hostToDevice, z_levels_h.begin(), z_levels_h.end(), z_levels_d.begin());

    switch(terrain_smoothing) {
    case 0: // BTF Method
    {
        for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Note that this box is nodal because it is based on z_phys_nd
            const Box&  bx = mfi.validbox();

            int k0 = bx.smallEnd()[2];

            // Grown box with corrected ghost cells at top
            Box gbx = mfi.growntilebox(ngrowVect);

            if (bx.smallEnd(2) == domlo_z) {
                gbx.setSmall(2,domlo_z);
            } else {
                gbx.growLo(2,-1);
            }
            if (bx.bigEnd(2) == domhi_z) {
                gbx.setBig(2,domhi_z);
            } else {
                gbx.growHi(2,-1);
            }

            Array4<Real> const& z_arr = z_phys_nd.array(mfi);
            auto const&         z_lev = z_levels_d.data();

            //
            // Vertical grid stretching using BTF model from p2163 of Klemp2011
            // z_levels are only defined from k = dom_lo to dom_hi (nodal)
            //
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int ii = amrex::max(amrex::min(i,imax),imin);
                int jj = amrex::max(amrex::min(j,jmax),jmin);

                //
                // Start with flat z_lev set either with uniform cell size or specified z_levels
                // If k0 = 0 then z_arr at k0 has already been filled from the terrain data
                // If k0 > 0 then z_arr at k0 has already been filled from interpolation
                //
                Real z         = z_lev[k];
                Real z_sfc     = z_arr(ii,jj,k0);
                Real z_lev_sfc = z_lev[k0];

                z_arr(i,j,k) = ( (z_sfc - z_lev_sfc) * z_top +
                                 (z_top - z_sfc    ) * z     ) / (z_top - z_lev_sfc);
            });

            if (k0 == 0) {
                // Fill lateral boundaries below the bottom surface
                ParallelFor(makeSlab(gbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    z_arr(i,j,-1) = 2.0*z_arr(i,j,0) - z_arr(i,j,1);
                });
            }
        } // mfi
        break;
    } // case 0

    case 1: // STF Method
    {
        // Get Multifab spanning domain with 1 level of ghost cells
        MultiFab h_mf(    z_phys_nd.boxArray(), z_phys_nd.DistributionMap(), 1, ngrow+1);
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
            mf2d = MultiFab(ba2d, h_mf.DistributionMap(), 1, ngrow, MFInfo().SetAlloc(false));
        }

        // Get MultiArray4s from the multifabs
        MultiArray4<Real> const& ma_h_s     = h_mf.arrays();
        MultiArray4<Real> const& ma_h_s_old = h_mf_old.arrays();
        MultiArray4<Real> const& ma_z_phys  = z_phys_nd.arrays();

        // Bottom boundary
        int k0 = domlo_z;

        // Get max value
        max_h = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, mf2d, IntVect(ngrow,ngrow,0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> GpuTuple<Real>
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

                diff = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, mf2d, IntVect(ngrow,ngrow,0),
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int) noexcept
                        -> GpuTuple<Real>
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

                ParallelDescriptor::ReduceRealMin(diff);

                iter++;

                //fill ghost points
                h_mf_old.FillBoundary(geom.periodicity());

            } //while

            auto const& z_lev_d = z_levels_d.data();

            //Populate z_phys_nd by solving z_arr(i,j,k) = z + A*h_s(i,j,k)
            for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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
                    if (k == 1) {
                        z_arr(i,j,k0-1) = 2.0*z_arr(i,j,k0) - z_arr(i,j,k);
                    }
                });
            } // mfi
            } // k

            Gpu::streamSynchronize();

            break;
        } // case 1

        case 2: // Sullivan TF Method
        {
            int k0    = 0;

            for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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
                    z_arr(i,j,k) = z + (std::pow((1. - (z/z_top)),omega) * z_arr(ii,jj,k0));

                    // Fill lateral boundaries and below the bottom surface
                    if (k == k0) {
                        z_arr(i,j,k0  ) = z_arr(ii,jj,k0);
                    }
                    if (k == 1) {
                        z_arr(i,j,k0-1) = 2.0*z_arr(ii,jj,k0) - z_arr(i,j,k);
                    }
                });
            } // mfi
            break;
        } // case 2

        case 3: // Debugging Test Method -- applies Sullivan TF starting at k = 1 so that domain does not change size
        {
            int k0    = 0;

            for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
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
                        z_arr(i,j,k) = z + (std::pow((1. - (z/z_top)),omega) * z_arr(ii,jj,k0));
                    }
                });
                gbx.setBig(2,0);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    z_arr(i,j,k  ) = 0.0;
                    z_arr(i,j,k-1) = -z_arr(i,j,k+1);
                });
            } // mfi
            break;
        } // case 3
    } //switch
}

/**
 * Computation of detJ at cell-center
 */
void
make_J (const Geometry& geom,
        MultiFab& z_phys_nd,
        MultiFab& detJ_cc)
{
    const auto *dx = geom.CellSize();
    Real dzInv = 1.0/dx[2];

    // Domain valid box (z_nd is nodal)
    const Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // Number of ghost cells
    int ngrow= detJ_cc.nGrow();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(detJ_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ngrow);
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }
        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > detJ = detJ_cc.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               detJ(i, j, k) = .25 * dzInv * (
                       z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
                      -z_nd(i,j,k  ) - z_nd(i+1,j,k  ) - z_nd(i,j+1,k  ) - z_nd(i+1,j+1,k  ) );
       });
    }
    detJ_cc.FillBoundary(geom.periodicity());
}

/**
 * Computation of area fractions on faces
 */
void
make_areas (const Geometry& geom,
            MultiFab& z_phys_nd, MultiFab& ax,
            MultiFab& ay, MultiFab& az)
{
    const auto* dx = geom.CellSize();
    Real dzInv = 1.0/dx[2];

    // Domain valid box (z_nd is nodal)
    const Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // The z-faces are always full when using terrain-fitted coordinates
    az.setVal(1.0);

    //
    // x-areas
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(ax, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ax.nGrow());
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }

        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > ax_arr = ax.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               ax_arr(i, j, k) = .5 * dzInv * (
                       z_nd(i,j,k+1) + z_nd(i,j+1,k+1) - z_nd(i,j,k) - z_nd(i,j+1,k));
        });
    }

    //
    // y-areas
    //
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(ay, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ay.nGrow());
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }

        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > ay_arr = ay.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
               ay_arr(i, j, k) = .5 * dzInv * (
                       z_nd(i,j,k+1) + z_nd(i+1,j,k+1) - z_nd(i,j,k) - z_nd(i+1,j,k));
        });
    }

    ax.FillBoundary(geom.periodicity());
    ay.FillBoundary(geom.periodicity());
    az.FillBoundary(geom.periodicity());
}

/**
 * Computation of z_phys at cell-center
 */
void
make_zcc (const Geometry& geom,
          MultiFab& z_phys_nd,
          MultiFab& z_phys_cc)
{
    // Domain valid box (z_nd is nodal)
    const Box& domain = geom.Domain();
    int domlo_z = domain.smallEnd(2);

    // Number of ghost cells
    int ngrow= z_phys_cc.nGrow();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(z_phys_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gbx = mfi.growntilebox(ngrow);
        if (gbx.smallEnd(2) < domlo_z) {
            gbx.setSmall(2,domlo_z);
        }
        Array4<Real const> z_nd = z_phys_nd.const_array(mfi);
        Array4<Real      > z_cc = z_phys_cc.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
           z_cc(i, j, k) = .125 * ( z_nd(i,j,k  ) + z_nd(i+1,j,k  ) + z_nd(i,j+1,k  ) + z_nd(i+1,j+1,k  )
                                   +z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1) );
       });
    }
    z_phys_cc.FillBoundary(geom.periodicity());
}
