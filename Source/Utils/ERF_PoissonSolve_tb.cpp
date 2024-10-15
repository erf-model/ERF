#include "ERF.H"
#include "ERF_Utils.H"

#ifdef ERF_USE_POISSON_SOLVE
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#endif

using namespace amrex;

/**
 * Project the single-level velocity field to enforce incompressibility with a
 * thin body
 */
void ERF::project_velocities_tb (int lev, Real l_dt, Vector<MultiFab>& vmf, MultiFab& pmf)
{
#ifdef ERF_USE_POISSON_SOLVE
    BL_PROFILE("ERF::project_velocities_tb()");
    AMREX_ALWAYS_ASSERT(!solverChoice.use_terrain);

    // Make sure the solver only sees the levels over which we are solving
    Vector<BoxArray>            ba_tmp;   ba_tmp.push_back(vmf[Vars::cons].boxArray());
    Vector<DistributionMapping> dm_tmp;   dm_tmp.push_back(vmf[Vars::cons].DistributionMap());
    Vector<Geometry>          geom_tmp; geom_tmp.push_back(geom[lev]);
    // amrex::Print() << "AT LEVEL " << lev << " BA FOR POISSON SOLVE " << vmf[Vars::cons].boxArray() << std::endl;

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
        p_mlpoisson = std::make_unique<MLPoisson>(geom_tmp, ba_tmp, dm_tmp, info);
    }

    auto bclo = get_projection_bc(Orientation::low);
    auto bchi = get_projection_bc(Orientation::high);
    bool need_adjust_rhs = (projection_has_dirichlet(bclo) || projection_has_dirichlet(bchi)) ? false : true;
    p_mlpoisson->setDomainBC(bclo, bchi);

    if (lev > 0) {
        p_mlpoisson->setCoarseFineBC(nullptr, ref_ratio[lev-1], LinOpBCType::Neumann);
    }

    p_mlpoisson->setLevelBC(0, nullptr);

    Vector<MultiFab> rhs;
    Vector<MultiFab> phi;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes;
    Vector<Array<MultiFab,AMREX_SPACEDIM> > deltaf; // f^* - f^{n-1}
    Vector<Array<MultiFab,AMREX_SPACEDIM> > u_plus_dtdf; // u + dt*deltaf

    // Used to pass array of const MFs to ComputeDivergence
    Array<MultiFab const*, AMREX_SPACEDIM> u;

    rhs.resize(1);
    phi.resize(1);
    fluxes.resize(1);
    deltaf.resize(1);
    u_plus_dtdf.resize(1);

    rhs[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    phi[0].define(ba_tmp[0], dm_tmp[0], 1, 0);
    rhs[0].setVal(0.0);
    phi[0].setVal(0.0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
             fluxes[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
        u_plus_dtdf[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);

             deltaf[0][idim].define(convert(ba_tmp[0], IntVect::TheDimensionVector(idim)), dm_tmp[0], 1, 0);
             deltaf[0][idim].setVal(0.0); // start with f^* == f^{n-1}
    }

#if 0
    // DEBUG
    u[0] = &(vmf[Vars::xvel]);
    u[1] = &(vmf[Vars::yvel]);
    u[2] = &(vmf[Vars::zvel]);
    computeDivergence(rhs[0], u, geom[0]);
    Print() << "Max norm of divergence before solve at level 0 : " << rhs[0].norm0() << std::endl;
#endif

    for (int itp = 0; itp < solverChoice.ncorr; ++itp)
    {
        // Calculate u + dt*deltaf
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Copy(u_plus_dtdf[0][idim], deltaf[0][idim], 0, 0, 1, 0);
            u_plus_dtdf[0][0].mult(-l_dt,0,1,0);
        }
        MultiFab::Add(u_plus_dtdf[0][0], vmf[Vars::xvel], 0, 0, 1, 0);
        MultiFab::Add(u_plus_dtdf[0][1], vmf[Vars::yvel], 0, 0, 1, 0);
        MultiFab::Add(u_plus_dtdf[0][2], vmf[Vars::zvel], 0, 0, 1, 0);

        u[0] = &(u_plus_dtdf[0][0]);
        u[1] = &(u_plus_dtdf[0][1]);
        u[2] = &(u_plus_dtdf[0][2]);
        computeDivergence(rhs[0], u, geom_tmp[0]);

#if 0
        // DEBUG
        if (itp==0) {
            for (MFIter mfi(rhs[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real const>& divU = rhs[0].const_array(mfi);
                const Array4<Real const>& uarr = vmf[Vars::xvel].const_array(mfi);
                const Array4<Real const>& varr = vmf[Vars::yvel].const_array(mfi);
                const Array4<Real const>& warr = vmf[Vars::zvel].const_array(mfi);
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
#endif

        // If all Neumann BCs, adjust RHS to make sure we can converge
        if (need_adjust_rhs) {
            Real offset = volWgtSumMF(lev, rhs[0], 0, *mapfac_m[lev], false, false);
            // amrex::Print() << "Poisson solvability offset = " << offset << std::endl;
            rhs[0].plus(-offset, 0, 1);
        }

        // Initialize phi to 0
        phi[0].setVal(0.0);

        MLMG mlmg(*p_mlpoisson);
        int max_iter = 100;
        mlmg.setMaxIter(max_iter);

        mlmg.setVerbose(mg_verbose);
        //mlmg.setBottomVerbose(mg_verbose);

        // solve for dt*p
        mlmg.solve(GetVecOfPtrs(phi),
                   GetVecOfConstPtrs(rhs),
                   solverChoice.poisson_reltol,
                   solverChoice.poisson_abstol);

        mlmg.getFluxes(GetVecOfArrOfPtrs(fluxes));

        // Calculate new intermediate body force with updated gradp
        if (thin_xforce[lev]) {
            MultiFab::Copy(   deltaf[0][0], fluxes[0][0], 0, 0, 1, 0);
            ApplyInvertedMask(deltaf[0][0], *xflux_imask[0]);
        }
        if (thin_yforce[lev]) {
            MultiFab::Copy(   deltaf[0][1], fluxes[0][1], 0, 0, 1, 0);
            ApplyInvertedMask(deltaf[0][1], *yflux_imask[0]);
        }
        if (thin_zforce[lev]) {
            MultiFab::Copy(   deltaf[0][2], fluxes[0][2], 0, 0, 1, 0);
            ApplyInvertedMask(deltaf[0][2], *zflux_imask[0]);
        }

        // DEBUG
        //        for (MFIter mfi(rhs[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        //        {
        //            const Box& bx = mfi.tilebox();
        //            const Array4<Real const>& dfz_arr = deltaf[0][2].const_array(mfi);
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

        // Update pressure variable with phi -- note that phi is change in pressure, not the full pressure
        MultiFab::Saxpy(pmf, 1.0, phi[0],0,0,1,0);

        // Subtract grad(phi) from the velocity components
        Real beta = 1.0;
        MultiFab::Saxpy(vmf[Vars::xvel], beta, fluxes[0][0], 0, 0, 1, 0);
        MultiFab::Saxpy(vmf[Vars::yvel], beta, fluxes[0][1], 0, 0, 1, 0);
        MultiFab::Saxpy(vmf[Vars::zvel], beta, fluxes[0][2], 0, 0, 1, 0);
        if (thin_xforce[lev]) {
            ApplyMask(vmf[Vars::xvel], *xflux_imask[0]);
        }
        if (thin_yforce[lev]) {
            ApplyMask(vmf[Vars::yvel], *yflux_imask[0]);
        }
        if (thin_zforce[lev]) {
            ApplyMask(vmf[Vars::zvel], *zflux_imask[0]);
        }
    } // itp: pressure-force iterations

    // Subtract grad(phi) from the velocity components
//    Real beta = 1.0;
//    for (int ilev = lev_min; ilev <= lev_max; ++ilev) {
//        MultiFab::Saxpy(vmf[Vars::xvel], beta, fluxes[0][0], 0, 0, 1, 0);
//        MultiFab::Saxpy(vmf[Vars::yvel], beta, fluxes[0][1], 0, 0, 1, 0);
//        MultiFab::Saxpy(vmf[Vars::zvel], beta, fluxes[0][2], 0, 0, 1, 0);
//        if (thin_xforce[lev]) {
//            ApplyMask(vmf[Vars::xvel], *xflux_imask[0]);
//        }
//        if (thin_yforce[lev]) {
//            ApplyMask(vmf[Vars::yvel], *yflux_imask[0]);
//        }
//        if (thin_zforce[lev]) {
//            ApplyMask(vmf[Vars::zvel], *zflux_imask[0]);
//        }
//    }

#if 0
    // Confirm that the velocity is now divergence free
    u[0] = &(vmf[Vars::xvel]);
    u[1] = &(vmf[Vars::yvel]);
    u[2] = &(vmf[Vars::zvel]);
    computeDivergence(rhs[0], u, geom_tmp[0]);
    Print() << "Max norm of divergence after solve at level " << lev << " : " << rhs[0].norm0() << std::endl;

#if 0
    for (MFIter mfi(rhs[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real const>& divU = rhs[0].const_array(mfi);
        const Array4<Real const>& uarr = u[0]->const_array(mfi);
        const Array4<Real const>& varr = u[1]->const_array(mfi);
        const Array4<Real const>& warr = u[2]->const_array(mfi);
        const Array4<Real const>& fzarr = thin_zforce[0]->const_array(mfi);
        const Array4<Real const>& dfzarr = deltaf[0][2].const_array(mfi);
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
    } // mfi
#endif
#endif
#else
    amrex::ignore_unused(lev);
    amrex::ignore_unused(l_dt);
    amrex::ignore_unused(vmf);
    amrex::ignore_unused(pmf);
#endif
}
