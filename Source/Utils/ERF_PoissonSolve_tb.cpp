#include "ERF.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include "Utils.H"

#ifdef ERF_USE_POISSON_SOLVE

using namespace amrex;

/**
 * Project the single-level velocity field to enforce incompressibility with a
 * thin body
 */
void ERF::project_velocities_tb (Vector<MultiFab>& vmf, const Real dt)
{
    Vector<Vector<MultiFab>> tmpmf(1);
    for (auto& mf : vmf) {
        tmpmf[0].emplace_back(mf, amrex::make_alias, 0, mf.nComp());
    }
    project_velocities_tb(tmpmf, dt);
}

/**
 * Project the multi-level velocity field to enforce incompressibility with a
 * thin body; we iterate on the pressure solve to obtain the body force that
 * satisfies no penetration
 */
void
ERF::project_velocities_tb (Vector<Vector<MultiFab>>& vars, const Real dt)
{
    BL_PROFILE("ERF::project_velocities()");

    const int nlevs = geom.size();

    // Use the default settings
    LPInfo info;
    std::unique_ptr<MLPoisson> p_mlpoisson;
#if 0
    if (overset_imask[0]) {
        // Add overset mask around thin body
        p_mlpoisson = std::make_unique<MLPoisson>(geom, grids, dmap, GetVecOfConstPtrs(overset_imask), info);
    }
    else
#endif
    {
        // Use the default settings
        p_mlpoisson = std::make_unique<MLPoisson>(geom, grids, dmap, info);
    }

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    bool need_adjust_rhs = (projection_has_dirichlet(bclo) || projection_has_dirichlet(bchi)) ? false : true;
    p_mlpoisson->setDomainBC(bclo, bchi);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
       p_mlpoisson->setLevelBC(0, nullptr);
    }

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > deltaf; // f^* - f^{n-1}
    Vector<Array<MultiFab,AMREX_SPACEDIM> > u_plus_dtdf; // u + dt*deltaf

    // Used to pass array of const MFs to ComputeDivergence
    Array<MultiFab const*, AMREX_SPACEDIM> u;

    rhs.resize(nlevs);
    phi.resize(nlevs);
    fluxes.resize(nlevs);
    deltaf.resize(nlevs);
    u_plus_dtdf.resize(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev) {
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev].setVal(0.0);
        phi[ilev].setVal(0.0);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            fluxes[ilev][idim].define(
                convert(grids[ilev], IntVect::TheDimensionVector(idim)),
                dmap[ilev], 1, 0);
            u_plus_dtdf[ilev][idim].define(
                convert(grids[ilev], IntVect::TheDimensionVector(idim)),
                dmap[ilev], 1, 0);
            deltaf[ilev][idim].define(
                convert(grids[ilev], IntVect::TheDimensionVector(idim)),
                dmap[ilev], 1, 0);
            deltaf[ilev][idim].setVal(0.0); // start with f^* == f^{n-1}
        }
    }

    // DEBUG
    u[0] = &(vars[0][Vars::xvel]);
    u[1] = &(vars[0][Vars::yvel]);
    u[2] = &(vars[0][Vars::zvel]);
    computeDivergence(rhs[0], u, geom[0]);
    Print() << "Max norm of divergence before solve at level 0 : " << rhs[0].norm0() << std::endl;

    for (int itp = 0; itp < solverChoice.ncorr; ++itp)
    {
        for (int ilev = 0; ilev < nlevs; ++ilev)
        {
            // Calculate u + dt*deltaf
            for (int idim = 0; idim < 3; ++idim) {
                MultiFab::Copy(u_plus_dtdf[ilev][idim], deltaf[ilev][idim], 0, 0, 1, 0);
                u_plus_dtdf[ilev][0].mult(-dt,0,1,0);
            }
            MultiFab::Add(u_plus_dtdf[ilev][0], vars[ilev][Vars::xvel], 0, 0, 1, 0);
            MultiFab::Add(u_plus_dtdf[ilev][1], vars[ilev][Vars::yvel], 0, 0, 1, 0);
            MultiFab::Add(u_plus_dtdf[ilev][2], vars[ilev][Vars::zvel], 0, 0, 1, 0);

            u[0] = &(u_plus_dtdf[ilev][0]);
            u[1] = &(u_plus_dtdf[ilev][1]);
            u[2] = &(u_plus_dtdf[ilev][2]);
            computeDivergence(rhs[ilev], u, geom[ilev]);

            // DEBUG
            if (itp==0) {
                for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real const>& divU = rhs[ilev].const_array(mfi);
                    const Array4<Real const>& uarr = vars[ilev][Vars::xvel].const_array(mfi);
                    const Array4<Real const>& varr = vars[ilev][Vars::yvel].const_array(mfi);
                    const Array4<Real const>& warr = vars[ilev][Vars::zvel].const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        if ((i>=120) && (i<=139) && (j==0) && ((k>=127)&&(k<=128))) {
                            amrex::AllPrint() << "before project div"<<IntVect(i,j,k)<<" = "<< divU(i,j,k)
                                << " u: " << uarr(i,j,k) << " " << uarr(i+1,j,k)
                                << " v: " << varr(i,j,k) << " " << varr(i,j+1,k)
                                << " w: " << warr(i,j,k) << " " << warr(i,j,k+1)
                                << std::endl;
                        }
                    });
                }
            }

            // If all Neumann BCs, adjust RHS to make sure we can converge
            if (need_adjust_rhs) {
                Real offset = volWgtSumMF(ilev, rhs[ilev], 0, *mapfac_m[ilev], false, false);
                amrex::Print() << "Poisson solvability offset = " << offset << std::endl;
                rhs[ilev].plus(-offset, 0, 1);
            }
        }

        // Initialize phi to 0
        for (int ilev = 0; ilev < nlevs; ++ilev)
        {
            phi[ilev].setVal(0.0);
        }

        MLMG mlmg(*p_mlpoisson);
        int max_iter = 100;
        mlmg.setMaxIter(max_iter);

        int verbose = 1;
        mlmg.setVerbose(verbose);
        //mlmg.setBottomVerbose(verbose);

        // solve for dt*p
        mlmg.solve(GetVecOfPtrs(phi),
                   GetVecOfConstPtrs(rhs),
                   solverChoice.poisson_reltol,
                   solverChoice.poisson_abstol);

        mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

        // Calculate new intermediate body force with updated gradp
        for (int ilev = 0; ilev < nlevs; ++ilev)
        {
            if (thin_xforce[ilev]) {
                MultiFab::Copy(   deltaf[ilev][0], fluxes[ilev][0], 0, 0, 1, 0);
                ApplyInvertedMask(deltaf[ilev][0], *xflux_imask[ilev]);
            }
            if (thin_yforce[ilev]) {
                MultiFab::Copy(   deltaf[ilev][1], fluxes[ilev][1], 0, 0, 1, 0);
                ApplyInvertedMask(deltaf[ilev][1], *yflux_imask[ilev]);
            }
            if (thin_zforce[ilev]) {
                MultiFab::Copy(   deltaf[ilev][2], fluxes[ilev][2], 0, 0, 1, 0);
                ApplyInvertedMask(deltaf[ilev][2], *zflux_imask[ilev]);
            }

            // DEBUG
    //        for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    //        {
    //            const Box& bx = mfi.tilebox();
    //            const Array4<Real const>& dfz_arr = deltaf[ilev][2].const_array(mfi);
    //            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    //            {
    //                if ((i>=120) && (i<=139) && (j==0) && (k==128)) {
    //                    amrex::AllPrint()
    //                        << " piter" << itp
    //                        << " dfz"<<IntVect(i,j,k)<<" = "<< dfz_arr(i,j,k)
    //                        << std::endl;
    //                }
    //            });
    //        }
        }

        // Subtract grad(phi) from the velocity components
        Real beta = 1.0;
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            MultiFab::Saxpy(vars[ilev][Vars::xvel], beta, fluxes[ilev][0], 0, 0, 1, 0);
            MultiFab::Saxpy(vars[ilev][Vars::yvel], beta, fluxes[ilev][1], 0, 0, 1, 0);
            MultiFab::Saxpy(vars[ilev][Vars::zvel], beta, fluxes[ilev][2], 0, 0, 1, 0);
            if (thin_xforce[ilev]) {
                ApplyMask(vars[ilev][Vars::xvel], *xflux_imask[ilev]);
            }
            if (thin_yforce[ilev]) {
                ApplyMask(vars[ilev][Vars::yvel], *yflux_imask[ilev]);
            }
            if (thin_zforce[ilev]) {
                ApplyMask(vars[ilev][Vars::zvel], *zflux_imask[ilev]);
            }
        }
    } // pressure-force iterations

    // Subtract grad(phi) from the velocity components
//    Real beta = 1.0;
//    for (int ilev = 0; ilev < nlevs; ++ilev) {
//        MultiFab::Saxpy(vars[ilev][Vars::xvel], beta, fluxes[ilev][0], 0, 0, 1, 0);
//        MultiFab::Saxpy(vars[ilev][Vars::yvel], beta, fluxes[ilev][1], 0, 0, 1, 0);
//        MultiFab::Saxpy(vars[ilev][Vars::zvel], beta, fluxes[ilev][2], 0, 0, 1, 0);
//        if (thin_xforce[ilev]) {
//            ApplyMask(vars[ilev][Vars::xvel], *xflux_imask[ilev]);
//        }
//        if (thin_yforce[ilev]) {
//            ApplyMask(vars[ilev][Vars::yvel], *yflux_imask[ilev]);
//        }
//        if (thin_zforce[ilev]) {
//            ApplyMask(vars[ilev][Vars::zvel], *zflux_imask[ilev]);
//        }
//    }

    // Average down the velocity from finest to coarsest to ensure consistency across levels
    int finest_level = nlevs - 1;
    Array<MultiFab const*, AMREX_SPACEDIM> u_fine;
    Array<MultiFab      *, AMREX_SPACEDIM> u_crse;
    for (int ilev = finest_level; ilev > 0; --ilev)
    {
        IntVect rr  = geom[ilev].Domain().size() / geom[ilev-1].Domain().size();
        u_fine[0] = &(vars[ilev  ][Vars::xvel]);
        u_fine[1] = &(vars[ilev  ][Vars::yvel]);
        u_fine[2] = &(vars[ilev  ][Vars::zvel]);
        u_crse[0] = &(vars[ilev-1][Vars::xvel]);
        u_crse[1] = &(vars[ilev-1][Vars::yvel]);
        u_crse[2] = &(vars[ilev-1][Vars::zvel]);
        average_down_faces(u_fine, u_crse, rr, geom[ilev-1]);
    }

#if 1
    // Confirm that the velocity is now divergence free
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        u[0] = &(vars[ilev][Vars::xvel]);
        u[1] = &(vars[ilev][Vars::yvel]);
        u[2] = &(vars[ilev][Vars::zvel]);
        computeDivergence(rhs[ilev], u, geom[ilev]);
        Print() << "Max norm of divergence after solve at level " << ilev << " : " << rhs[ilev].norm0() << std::endl;

        for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Array4<Real const>& divU = rhs[ilev].const_array(mfi);
            const Array4<Real const>& uarr = u[0]->const_array(mfi);
            const Array4<Real const>& varr = u[1]->const_array(mfi);
            const Array4<Real const>& warr = u[2]->const_array(mfi);
            const Array4<Real const>& fzarr = thin_zforce[ilev]->const_array(mfi);
            const Array4<Real const>& dfzarr = deltaf[ilev][2].const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if ((i>=120) && (i<=139) && (j==0) && ((k>=127)&&(k<=128))) {
                    amrex::AllPrint() << "after project div"<<IntVect(i,j,k)<<" = "<< divU(i,j,k)
                        << " u: " << uarr(i,j,k) << " " << uarr(i+1,j,k)
                        << " v: " << varr(i,j,k) << " " << varr(i,j+1,k)
                        << " w: " << warr(i,j,k) << " " << warr(i,j,k+1)
                        << " fz = " << fzarr(i,j,k) << " + " << dfzarr(i,j,k)
                        << std::endl;
                }
            });
        }
    }
#endif
}
#endif
