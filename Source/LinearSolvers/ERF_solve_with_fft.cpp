#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

#ifdef ERF_USE_FFT
Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM>
ERF::get_fft_bc () const noexcept
{
    //
    // This logic only works for level 0
    // TODO: fix for level > 0
    //
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> r;

    for (int dir = 0; dir <= 1; dir++) {
        if (geom[0].isPeriodic(dir)) {
            r[dir] = std::make_pair(FFT::Boundary::periodic,FFT::Boundary::periodic);
        }
    } // dir

    for (OrientationIter ori; ori != nullptr; ++ori) {
        const int dir  = ori().coordDir();
        if (!geom[0].isPeriodic(dir) && ori().faceDir() == Orientation::low) {
            auto bc_type_lo = domain_bc_type[Orientation(dir,Orientation::low)];
            auto bc_type_hi = domain_bc_type[Orientation(dir,Orientation::high)];
            if ( (bc_type_lo == "Outflow" || bc_type_lo == "Open") &&
                 (bc_type_hi == "Outflow" || bc_type_hi == "Open") ) {
                r[dir] = std::make_pair(FFT::Boundary::odd,FFT::Boundary::odd);
            } else if ( (bc_type_lo != "Outflow" && bc_type_lo != "Open") &&
                        (bc_type_hi == "Outflow" || bc_type_hi == "Outflow") ) {
                r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::odd);
            } else if ( (bc_type_lo == "Outflow" || bc_type_lo == "Open") &&
                        (bc_type_hi != "Outflow" && bc_type_hi != "Outflow") ) {
                r[dir] = std::make_pair(FFT::Boundary::odd,FFT::Boundary::even);
            } else {
                r[dir] = std::make_pair(FFT::Boundary::even,FFT::Boundary::even);
            }
        } // not periodic
    } // ori

    return r;
}

/**
 * Solve the Poisson equation using FFT
 * Note that the level may or may not be level 0.
 */
void ERF::solve_with_fft (int lev, MultiFab& rhs, MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM>& fluxes)
{
    BL_PROFILE("ERF::solve_with_fft()");

    // For now the logic only works correctly for lev == 0
    // TODO: enable FFT for nested domains at lev > 0
    AMREX_ALWAYS_ASSERT(lev == 0);

    bool l_use_terrain = SolverChoice::terrain_type != TerrainType::None;

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);

    // amrex::Print() << "BCLO " << bclo[0] << " " << bclo[1] << " " << bclo[2] << std::endl;
    // amrex::Print() << "BCHI " << bchi[0] << " " << bchi[1] << " " << bchi[2] << std::endl;

    auto dxInv = geom[lev].InvCellSizeArray();

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    // ****************************************************************************
    // FFT solve
    // ****************************************************************************
    AMREX_ALWAYS_ASSERT(lev == 0);
    //
    // No terrain or stretched grids
    // This calls the full 3D FFT solver with bc's set through bc_fft
    //
    if (!l_use_terrain)
    {
        if (mg_verbose > 0) {
            amrex::Print() << "Using the 3D FFT solver..." << std::endl;
        }
        if (!m_3D_poisson) {
            auto bc_fft = get_fft_bc();
            m_3D_poisson = std::make_unique<FFT::Poisson<MultiFab>>(Geom(0),bc_fft);
        }
        m_3D_poisson->solve(phi, rhs);

    //
    // Stretched grids
    // This calls the hybrid 2D FFT solver + tridiagonal in z
    //
    // For right now we can only do this solve for periodic in the x- and y-directions
    // We assume Neumann at top and bottom z-boundaries
    // This will be generalized in future
    //
    //
    } else if (l_use_terrain && SolverChoice::terrain_is_flat)
    {
        if (mg_verbose > 0) {
            amrex::Print() << "Using the hybrid FFT solver..." << std::endl;
        }
        if (!m_2D_poisson) {
            m_2D_poisson = std::make_unique<FFT::PoissonHybrid<MultiFab>>(Geom(0));
        }
        m_2D_poisson->solve(phi, rhs, stretched_dz_d[lev]);

    } else {
        amrex::Abort("FFT isn't appropriate for spatially varying terrain");
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const&  p_arr  = phi.array(mfi);

        Box const& xbx = mfi.nodaltilebox(0);
        const Real dx_inv = dxInv[0];
        Array4<Real> const& fx_arr  = fluxes[0].array(mfi);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i-1,j,k)) * dx_inv;
        });

        Box const& ybx = mfi.nodaltilebox(1);
        const Real dy_inv = dxInv[1];
        Array4<Real> const& fy_arr  = fluxes[1].array(mfi);
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dy_inv;
        });

        Box const& zbx = mfi.nodaltilebox(2);
        Array4<Real> const& fz_arr  = fluxes[2].array(mfi);
        if (l_use_terrain && SolverChoice::terrain_is_flat) {
            Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == dom_lo.z || k == dom_hi.z+1) {
                    fz_arr(i,j,k) = 0.0;
                } else {
                    Real dz = 0.5 * (stretched_dz_d_ptr[k] + stretched_dz_d_ptr[k-1]);
                    fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) / dz;
                }
            });
        } else { // no grid stretching
            const Real dz_inv = dxInv[2];
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == dom_lo.z || k == dom_hi.z+1) {
                    fz_arr(i,j,k) = 0.0;
                } else {
                    fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) * dz_inv;
                }
            });
        }
    } // mfi
}
#endif
