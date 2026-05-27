#include "ERF.H"
#include "ERF_SolverUtils.H"

using namespace amrex;

#ifdef ERF_USE_FFT
/**
 * Build the FFT solver(s) -- one per subdomain in each level
 */
void ERF::build_fft_solvers (int lev)
{
    BL_PROFILE("ERF::build_fft_solvers()");

    // Clear any FFT solvers that might already be built at this level
    if (m_3D_poisson.size() > lev) {
        for (int n = 0; n < m_3D_poisson[lev].size(); n++) {
            m_3D_poisson[lev][n].reset();
        }
    }
    if (m_2D_poisson.size() > lev) {
        for (int n = 0; n < m_2D_poisson[lev].size(); n++) {
            m_2D_poisson[lev][n].reset();
        }
    }

    if ( (solverChoice.mesh_type == MeshType::ConstantDz )  ||
         (solverChoice.mesh_type == MeshType::StretchedDz) )
    {
        bool will_solve_with_mlmg = (solverChoice.mesh_type == MeshType::ConstantDz);
        bool all_boxes_ok = true;
        for (int isub = 0; isub < subdomains[lev].size(); ++isub) {
            Box my_region(subdomains[lev][isub].minimalBox());
            bool boxes_make_rectangle = (my_region.numPts() == subdomains[lev][isub].numPts());
            if (!boxes_make_rectangle) {
                all_boxes_ok = false;
            }
        } // isub
        if (all_boxes_ok) {
            will_solve_with_mlmg = false;
        }

        for (int isub = 0; isub < subdomains[lev].size(); ++isub)
        {
            Box subdomain(subdomains[lev][isub].minimalBox());

            auto const dom_lo = lbound(Geom(lev).Domain());
            auto const dom_hi = ubound(Geom(lev).Domain());

            auto const sub_lo = lbound(subdomain);
            auto const sub_hi = ubound(subdomain);

            auto dx    = geom[lev].CellSizeArray();
            auto dxInv = geom[lev].InvCellSizeArray();

            auto bc_fft = get_fft_bc(geom[lev],domain_bc_type,subdomain,solverChoice.use_real_bcs);

            Geometry my_geom;

            Array<int,AMREX_SPACEDIM> is_per; is_per[0] = 0; is_per[1] = 0; is_per[2] = 0;
            if (geom[lev].isPeriodic(0) && sub_lo.x == dom_lo.x && sub_hi.x == dom_hi.x) { is_per[0] = 1;}
            if (geom[lev].isPeriodic(1) && sub_lo.y == dom_lo.y && sub_hi.y == dom_hi.y) { is_per[1] = 1;}

            int coord_sys = 0;

            // If subdomain == domain then we pass Geom(lev) to the FFT solver
            if (subdomain == Geom(lev).Domain()) {
                my_geom.define(Geom(lev).Domain(), Geom(lev).ProbDomain(), coord_sys, is_per);
            } else {
                // else we create a new geometry based only on the subdomain
                // The information in my_geom used by the FFT routines is:
                //   1) my_geom.Domain()
                //   2) my_geom.CellSize()
                //   3) my_geom.isAllPeriodic() / my_geom.periodicity()
                RealBox rb( sub_lo.x   *dx[0],  sub_lo.y   *dx[1],  sub_lo.z   *dx[2],
                           (sub_hi.x+1)*dx[0], (sub_hi.y+1)*dx[1], (sub_hi.z+1)*dx[2]);
                my_geom.define(subdomain, rb, coord_sys, is_per);
            }
            bool boxes_make_rectangle = (subdomain.numPts() == subdomains[lev][isub].numPts());
            if (!boxes_make_rectangle && !will_solve_with_mlmg) {
                amrex::Abort("FFT won't work unless the union of boxes is rectangular");
            }

            if (solverChoice.mesh_type == MeshType::ConstantDz && !will_solve_with_mlmg)
            {
                if (mg_verbose > 0) {
                    amrex::Print() << "Building the 3D FFT solver to be used on domain " << my_geom.Domain() << std::endl;
                }
                if (m_3D_poisson.size() <= lev) {
                    m_3D_poisson.resize(lev+1);
                }
                if (m_3D_poisson[lev].size() <= isub) {
                    m_3D_poisson[lev].resize(isub+1);
                }
                m_3D_poisson[lev][isub] = std::make_unique<FFT::Poisson<MultiFab>>(my_geom,bc_fft);

            }
            else if (solverChoice.mesh_type == MeshType::StretchedDz)
            {
                if (mg_verbose > 0) {
                    amrex::Print() << "Building the hybrid FFT solver to be used on domain " << my_geom.Domain() << std::endl;
                }
                if (m_2D_poisson.size() <= lev) {
                    m_2D_poisson.resize(lev+1);
                }
                if (m_2D_poisson[lev].size() <= isub) {
                    m_2D_poisson[lev].resize(isub+1);
                }
                m_2D_poisson[lev][isub] = std::make_unique<FFT::PoissonHybrid<MultiFab>>(my_geom,bc_fft);
            }
        } // isub
    } // constant or stretched dz
}

/**
 * Solve the Poisson equation using FFT
 * Note that the level may or may not be level zero
 */
void ERF::solve_with_fft (int lev, int isub, const Box& subdomain,
                          MultiFab& rhs, MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM>& fluxes)
{
    BL_PROFILE("ERF::solve_with_fft()");

    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);

    Vector<Real> stretched_dz_h_sub;
    Gpu::DeviceVector<Real> stretched_dz_d_sub;

    // ****************************************************************************
    // FFT solve
    // ****************************************************************************
    //
    // This calls the full 3D FFT solver with bc's set through bc_fft
    //
    if (solverChoice.mesh_type == MeshType::ConstantDz)
    {
        if (mg_verbose > 0) {
            amrex::Print() << "Using the 3D FFT solver on domain " << subdomain << std::endl;
        }
        m_3D_poisson[lev][isub]->solve(phi, rhs);

    //
    // Stretched grids
    //
    // This calls the hybrid 2D FFT solver + tridiagonal in z with lateral bc's set through bc_fft
    // and Neumann at top and bottom z-boundaries
    //
    }
    else
    {
        if (mg_verbose > 0) {
            amrex::Print() << "Using the hybrid FFT solver on domain " << subdomain << std::endl;
        }

        // Make new stretched_dz_sub array which only covers my_geom
        stretched_dz_h_sub.resize(subdomain.length(2));

        for (int k = sub_lo.z; k <= sub_hi.z; k++)
        {
            stretched_dz_h_sub[k-sub_lo.z] = stretched_dz_h[lev][k];
        }
        stretched_dz_d_sub.resize(stretched_dz_h_sub.size());
        Gpu::copy(Gpu::hostToDevice, stretched_dz_h_sub.begin(), stretched_dz_h_sub.end(), stretched_dz_d_sub.begin());

        m_2D_poisson[lev][isub]->solve(phi, rhs, stretched_dz_d_sub);
    }

    if (mg_verbose > 0) {
        amrex::Print() << "FFT solve complete on domain " << subdomain << std::endl;
    }

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);

    // ****************************************************************************
    // Compute fluxes which we will subtract from the momenta
    // ****************************************************************************
    auto dxInv = geom[lev].InvCellSizeArray();
    const Real dx_inv = dxInv[0]; const Real dy_inv = dxInv[1];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const&  p_arr  = phi.array(mfi);

        Box const& xbx = mfi.nodaltilebox(0);
        Array4<Real> const& fx_arr  = fluxes[0].array(mfi);
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i-1,j,k)) * dx_inv;
        });

        Box const& ybx = mfi.nodaltilebox(1);
        Array4<Real> const& fy_arr  = fluxes[1].array(mfi);
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j-1,k)) * dy_inv;
        });

        Box const& zbx = mfi.nodaltilebox(2);
        Array4<Real> const& fz_arr  = fluxes[2].array(mfi);
        if (solverChoice.mesh_type != MeshType::ConstantDz) {
            Real* stretched_dz_d_ptr = stretched_dz_d[lev].data();
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == sub_lo.z || k == sub_hi.z+1) {
                    fz_arr(i,j,k) = zero;
                } else {
                    Real dz = myhalf * (stretched_dz_d_ptr[k] + stretched_dz_d_ptr[k-1]);
                    fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) / dz;
                }
            });
        } else { // no grid stretching
            const Real dz_inv = dxInv[2];
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k == sub_lo.z || k == sub_hi.z+1) {
                    fz_arr(i,j,k) = zero;
                } else {
                    fz_arr(i,j,k) = -(p_arr(i,j,k) - p_arr(i,j,k-1)) * dz_inv;
                }
            });
        }
    } // mfi
}
#endif
