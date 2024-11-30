#include "ERF.H"
#include "ERF_FFT_Utils.H"

using namespace amrex;

#ifdef ERF_USE_FFT
/**
 * Solve the Poisson equation using FFT
 * Note that the level may or may not be level 0.
 */
void ERF::solve_with_fft (int lev, MultiFab& rhs, MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM>& fluxes)
{
    BL_PROFILE("ERF::solve_with_fft()");

    bool l_use_terrain = SolverChoice::terrain_type != TerrainType::None;

    auto const dom_lo = lbound(geom[lev].Domain());
    auto const dom_hi = ubound(geom[lev].Domain());

    auto dxInv = geom[lev].InvCellSizeArray();

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    Box bounding_box(rhs.boxArray().minimalBox());
    auto bc_fft = get_fft_bc(geom[lev],domain_bc_type,bounding_box);

    // ****************************************************************************
    // FFT solve
    // ****************************************************************************
    //
    // No terrain or stretched grids
    // This calls the full 3D FFT solver with bc's set through bc_fft
    //
    if (!l_use_terrain)
    {
        if (mg_verbose > 0) {
            amrex::Print() << "Using the 3D FFT solver..." << std::endl;
        }
        if (m_3D_poisson.size() <= lev) {
            m_3D_poisson.resize(lev+1);
            m_3D_poisson[lev] = std::make_unique<FFT::Poisson<MultiFab>>(Geom(lev),bc_fft);
        }
        m_3D_poisson[lev]->solve(phi, rhs);

    //
    // Stretched grids
    //
    // This calls the hybrid 2D FFT solver + tridiagonal in z with lateral bc's set through bc_fft
    // and Neumann at top and bottom z-boundaries
    //
    } else if (l_use_terrain && SolverChoice::terrain_is_flat)
    {
        if (mg_verbose > 0) {
            amrex::Print() << "Using the hybrid FFT solver..." << std::endl;
        }
        if (m_2D_poisson.size() <= lev) {
            m_2D_poisson.resize(lev+1);
            m_2D_poisson[lev] = std::make_unique<FFT::PoissonHybrid<MultiFab>>(Geom(lev),bc_fft);
        }
        m_2D_poisson[lev]->solve(phi, rhs, stretched_dz_d[lev]);

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
