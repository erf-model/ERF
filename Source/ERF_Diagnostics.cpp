/**
 * \file ERF_Diagnostics.cpp
 */

#include <limits>

#include "ERF.H"
#include "ERF_SrcHeaders.H"
#include "ERF_Utils.H"

using namespace amrex;

void
ERF::compute_max_pressure_gradient_diagnostic(int lev)
{
    // We don't require HSE when anelastic because the pressure gradient
    //    is computed from the Poisson solve
    if (solverChoice.anelastic[lev]) return;

    auto& lev_new = vars_new[lev];

    int ng = (solverChoice.terrain_type == TerrainType::EB) ? 3 : 1;

    const Real grav = solverChoice.gravity;

    Vector<MultiFab> gradp_temp;  gradp_temp.resize(AMREX_SPACEDIM);
    gradp_temp[0].define(vars_new[lev][Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, 0);
    gradp_temp[0].setVal(0.);
    gradp_temp[1].define(vars_new[lev][Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, 0);
    gradp_temp[1].setVal(0.);
    gradp_temp[2].define(vars_new[lev][Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, 0);
    gradp_temp[2].setVal(0.);

    int comp = 0;

    // Use this region to take max/min of gpx without including xlo,xhi if using real_bcs
    Box xface_domain = surroundingNodes(geom[lev].Domain(), 0);
    int ilo = xface_domain.smallEnd(0);
    int ihi = xface_domain.bigEnd(0);
    if (solverChoice.use_real_bcs) {
        xface_domain.growLo(0,-1);
        xface_domain.growHi(0,-1);
    }

    // Use this region to take max/min of gpy without including ylo,yhi if using real_bcs
    Box yface_domain = surroundingNodes(geom[lev].Domain(), 1);
    int jlo = yface_domain.smallEnd(1);
    int jhi = yface_domain.bigEnd(1);
    if (solverChoice.use_real_bcs) {
        yface_domain.growLo(1,-1);
        yface_domain.growHi(1,-1);
    }


    // Use this region to take max/min of gpz without including top and bottom faces
    Box zface_domain = surroundingNodes(geom[lev].Domain(), 2);
    int klo = zface_domain.smallEnd(2);
    int khi = zface_domain.bigEnd(2);

    zface_domain.growLo(2,-1);
    zface_domain.growHi(2,-1);

    // *******************************************************************************
    // First check that base state satisfies EOS
    // *******************************************************************************

    Print() << " " << std::endl;

    MultiFab  r_hse(base_state[lev], make_alias, BaseState::r0_comp , 1);
    MultiFab  p_hse(base_state[lev], make_alias, BaseState::p0_comp , 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp , 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);

    MultiFab dp(p_hse.boxArray(), p_hse.DistributionMap(), 1, 0);

    // Initialize to zero in case of EB covered cells
    dp.setVal(0.);

    for (MFIter mfi(dp); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        auto const  rhse_arr  =  r_hse.const_array(mfi);
        auto const  phse_arr  =  p_hse.const_array(mfi);
        auto const qvhse_arr  = qv_hse.const_array(mfi);
        auto const thhse_arr  = th_hse.const_array(mfi);
        auto       dpeos_arr  = dp.array(mfi);

        if (solverChoice.terrain_type != TerrainType::EB) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Real rhotheta = rhse_arr(i,j,k) * thhse_arr(i,j,k);
                dpeos_arr(i,j,k) = std::abs(getPgivenRTh(rhotheta, qvhse_arr(i,j,k)) - phse_arr(i,j,k));
            });
        } else {
            Array4<const Real> volfrac = (get_eb(lev).get_const_factory())->getVolFrac().const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (volfrac(i,j,k) > zero) {
                    Real rhotheta = rhse_arr(i,j,k) * thhse_arr(i,j,k);
                    dpeos_arr(i,j,k) = std::abs(getPgivenRTh(rhotheta, qvhse_arr(i,j,k)) - phse_arr(i,j,k));
                }
            });
        }
    }
    // This compares two pressures, so the tolerance has to be scaled to the
    // pressure rather than fixed in Pa.  The base state is O(1e5) Pa, where one
    // single-precision ulp is already ~1e-2 Pa, so the absolute 1e-8 Pa bound
    // used here previously could never be met in single precision even for a
    // perfectly consistent base state.  The double-precision coefficient below
    // reproduces that 1e-8 bound at 1e5 Pa, so the check is unchanged there.
#ifdef AMREX_USE_FLOAT
    const Real eos_rel_tol = Real(1.e-5);
#else
    const Real eos_rel_tol = Real(1.e-13);
#endif
    const Real p_hse_max = p_hse.max(0);
    const Real eos_tol = eos_rel_tol * p_hse_max;

    Real max_diff = dp.max(0);
    if (max_diff > eos_tol) {
        IntVect max_loc = dp.maxIndex(0);
        Print() << "Max value of |p_hse - p_eos| is " << max_diff << std::endl;
        Print() << " with max in cell " << max_loc << std::endl;
        Abort("Base state violates EOS ");
    } else {
        Print() << "Max value of |p_hse - p_eos| is " << max_diff
                << ", less than the tolerance " << eos_tol << std::endl;
    }

    // *******************************************************************************
    // Now compute pressure gradients for base state pressure
    // *******************************************************************************

    compute_gradp(p_hse, geom[lev], *z_phys_nd[lev].get(), *z_phys_cc[lev].get(), mapfac[lev],
                  get_eb(lev), gradp_temp, solverChoice);

    Real min_gpx = gradp_temp[0].min(xface_domain,comp);
    Real max_gpx = gradp_temp[0].max(xface_domain,comp);
    if (max_gpx != zero || min_gpx != zero) {
        Print() << "Min/max value of dp0/dx            are " << min_gpx << " " << max_gpx << std::endl;
        IntVect min_loc = gradp_temp[0].minIndex(comp);
        IntVect max_loc = gradp_temp[0].maxIndex(comp);
        if (min_loc[0] != ilo && min_loc[0] != ihi) amrex::Print() << " with min at face " << min_loc;
        if (max_loc[0] != ilo && max_loc[0] != ihi) amrex::Print() << " with max at face " << max_loc;
        Print() << std::endl;
    } else {
        Print() << "Min/max value of dp0/dx            are zero " << std::endl;
    }

    Real min_gpy = gradp_temp[1].min(yface_domain,comp);
    Real max_gpy = gradp_temp[1].max(yface_domain,comp);
    if (max_gpy != zero || min_gpy != zero) {
        Print() << "Min/max value of dp0/dy            are " << min_gpy << " " << max_gpy << std::endl;
        IntVect min_loc = gradp_temp[1].minIndex(comp);
        IntVect max_loc = gradp_temp[1].maxIndex(comp);
        if (min_loc[1] != jlo && min_loc[1] != jhi) amrex::Print() << " with min at face " << min_loc;
        if (max_loc[1] != jlo && max_loc[1] != jhi) amrex::Print() << " with max at face " << max_loc;
        Print() << std::endl;
    } else {
        Print() << "Min/max value of dp0/dy            are zero " << std::endl;
    }

    if (solverChoice.terrain_type != TerrainType::EB) {
        for (MFIter mfi(gradp_temp[2]); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            if (bx.bigEnd(2)   == khi) bx.growHi(2,-1);
            if (bx.smallEnd(2) == 0  ) bx.growLo(2,-1);
            auto        gpz_arr  = gradp_temp[2].array(mfi);
            auto const  rhse_arr  =  r_hse.const_array(mfi);
            auto const qvhse_arr  = qv_hse.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                gpz_arr(i,j,k) += grav * myhalf * ( rhse_arr(i,j,k  ) * (one + qvhse_arr(i,j,k  ))
                                                   +rhse_arr(i,j,k-1) * (one + qvhse_arr(i,j,k-1)) );
            });
        }
    // EB case: check HSE only for uncovered cells
    } else {
        for (MFIter mfi(gradp_temp[2]); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            if (bx.bigEnd(2)   == khi) bx.growHi(2,-1);
            if (bx.smallEnd(2) == 0  ) bx.growLo(2,-1);
            auto        gpz_arr  = gradp_temp[2].array(mfi);
            auto const  rhse_arr  =  r_hse.const_array(mfi);
            auto const qvhse_arr  = qv_hse.const_array(mfi);
            Array4<const Real> w_volfrac = (get_eb(lev).get_w_const_factory())->getVolFrac().const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                if (w_volfrac(i,j,k) > zero) {
                    gpz_arr(i,j,k) += grav * myhalf * ( rhse_arr(i,j,k  ) * (one + qvhse_arr(i,j,k  ))
                                                       +rhse_arr(i,j,k-1) * (one + qvhse_arr(i,j,k-1)) );
                }
            });
        }
    }

    // The residual tested below, dp0/dz + rho0*|g|, is a cancellation between two
    // terms of size rho0*|g|.  The smallest value it can take is therefore set by
    // the roundoff in p0 spread over one vertical cell, eps*|p0|/dz.  That floor
    // is grid dependent -- in single precision it is ~1e-4 Pa/m at dz = 100 m but
    // ~3e-4 Pa/m at dz = 40 m -- so no fixed tolerance can work across the range
    // of meshes the tests cover.  Use the larger of the historical fixed value and
    // that floor.  In double the floor is O(1e-12), so the fixed value always wins
    // and this check keeps exactly its previous behavior.
    const Real tol_fixed = (solverChoice.terrain_type == TerrainType::EB) ? Real(1.e-4)
                                                                         : Real(1.e-8);

    const Real dz_for_tol = (lev < static_cast<int>(dz_min.size()) && dz_min[lev] > Real(0))
                          ? dz_min[lev] : geom[lev].CellSize(2);

    const Real tol_roundoff = Real(16) * std::numeric_limits<Real>::epsilon()
                            * p_hse_max / dz_for_tol;

    const Real tol = std::max(tol_fixed, tol_roundoff);

    Real min_gpz = gradp_temp[2].min(zface_domain,comp);
    Real max_gpz = gradp_temp[2].max(zface_domain,comp);

    if (std::abs(max_gpz) > tol || std::abs(min_gpz) > tol) {
        IntVect min_loc = gradp_temp[2].minIndex(comp);
        IntVect max_loc = gradp_temp[2].maxIndex(comp);
        Print() << "Min/max value of dp0/dz + rho0*|g| are " << min_gpz << " " << max_gpz;
        if (min_loc[2] != klo && min_loc[2] != khi) amrex::Print() << " with min at face " << min_loc;
        if (max_loc[2] != klo && max_loc[2] != khi) amrex::Print() << " with max at face " << max_loc;
        amrex::Abort("Base state is too far out of HSE");
    } else {
        Print() << "Min/max value of dp0/dz + rho0*|g| are less than " << tol << std::endl;
    }
    Print() << " " << std::endl;

    if (!solverChoice.anelastic[lev]) {

        // *******************************************************************************
        // Now compute for full (moist) pressure
        // *******************************************************************************

        MultiFab p(p_hse.boxArray(), p_hse.DistributionMap(), 1, ng);
        MultiFab rho(lev_new[Vars::cons], make_alias, Rho_comp , 1);

        if (solverChoice.moisture_type != MoistureType::None) {

            for (MFIter mfi(rho); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.tilebox();
                gbx.grow(IntVect(ng,ng,ng));
                if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);

                const Array4<const Real>& cell_data = lev_new[Vars::cons].array(mfi);
                const Array4<const Real>&  r_arr = rho.array(mfi);
                const Array4<      Real>& pp_arr = p.array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real qv_for_p = cell_data(i,j,k,RhoQ1_comp)/r_arr(i,j,k);
                    pp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p);
                });
            }
            compute_gradp(p, geom[lev], *z_phys_nd[lev].get(), *z_phys_cc[lev].get(), mapfac[lev],
                          get_eb(lev), gradp_temp, solverChoice);

            min_gpx = gradp_temp[0].min(xface_domain,comp);
            max_gpx = gradp_temp[0].max(xface_domain,comp);
            if (max_gpx != zero || min_gpx != zero) {
                Print() << "Min/Max value of x-gradient of full (moist) pressure are " << min_gpx << " " << max_gpx;
                IntVect min_loc = gradp_temp[0].minIndex(comp);
                IntVect max_loc = gradp_temp[0].maxIndex(comp);
                if (min_loc[0] != ilo && min_loc[0] != ihi) amrex::Print() << " with min at face " << min_loc;
                if (max_loc[0] != ilo && max_loc[0] != ihi) amrex::Print() << " with max at face " << max_loc;
                Print() << std::endl;
            } else {
                Print() << "Min/max value of x-gradient of full (moist) pressure are zero " << std::endl;
            }

            min_gpy = gradp_temp[1].min(yface_domain,comp);
            max_gpy = gradp_temp[1].max(yface_domain,comp);
            if (max_gpy != zero || min_gpy != zero) {
                Print() << "Min/Max value of y-gradient of full (moist) pressure are " << min_gpy << " " << max_gpy;
                IntVect min_loc = gradp_temp[1].minIndex(comp);
                IntVect max_loc = gradp_temp[1].maxIndex(comp);
                if (min_loc[1] != jlo && min_loc[1] != jhi) amrex::Print() << " with min at face " << min_loc;
                if (max_loc[1] != jlo && max_loc[1] != jhi) amrex::Print() << " with max at face " << max_loc;
                Print() << std::endl;
            } else {
                Print() << "Min/max value of y-gradient of full (moist) pressure are zero " << std::endl;
            }

            MultiFab qt(rho.boxArray(), rho.DistributionMap(), 1, 1);
            int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
            make_qt(lev_new[Vars::cons], qt, n_qstate_into_total);

            for (MFIter mfi(gradp_temp[2]); mfi.isValid(); ++mfi)
            {
                Box bx = mfi.validbox();
                if (bx.bigEnd(2)   == khi) bx.growHi(2,-1);
                if (bx.smallEnd(2) == 0  ) bx.growLo(2,-1);
                auto      gpz_arr   = gradp_temp[2].array(mfi);
                auto const  r_arr   = rho.const_array(mfi);
                auto const qt_arr   =  qt.const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    gpz_arr(i,j,k) += grav * myhalf * (r_arr(i,j,k  )*(one+qt_arr(i,j,k  )) +
                                                       r_arr(i,j,k-1)*(one+qt_arr(i,j,k-1)) );
                });
            }

            min_gpz = gradp_temp[2].min(zface_domain,comp);
            max_gpz = gradp_temp[2].max(zface_domain,comp);
            if (max_gpz != zero || min_gpz != zero) {
                IntVect min_loc = gradp_temp[2].minIndex(comp);
                IntVect max_loc = gradp_temp[2].maxIndex(comp);
                Print() << "Min/max value of moist dp/dz + rho_m*|g|             are " << min_gpz << " " << max_gpz;
                if (min_loc[2] != klo && min_loc[2] != khi) amrex::Print() << " with min at face " << min_loc;
                if (max_loc[2] != klo && max_loc[2] != khi) amrex::Print() << " with max at face " << max_loc;
            } else {
                Print() << "Min/max value of moist dp/dz + rho_m*|g|  are zero " << std::endl;
            }
            Print() << " " << std::endl;
        } // if moist
    } // if !anelastic
}

void
ERF::compute_max_buoyancy_gradp_diagnostic (int lev)
{
    // When anelastic the pressure gradient comes from the projection, not from
    //    a perturbational pressure, so this diagnostic doesn't apply
    if (solverChoice.anelastic[lev]) return;

    auto& lev_new = vars_new[lev];

    const BoxArray&            ba = lev_new[Vars::cons].boxArray();
    const DistributionMapping& dm = lev_new[Vars::cons].DistributionMap();

    // *******************************************************************************
    // Compute the gradient of the perturbational pressure exactly as the dycore does
    // *******************************************************************************

    Vector<MultiFab> lgradp;  lgradp.resize(AMREX_SPACEDIM);
    lgradp[GpVars::gpx].define(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, 0);
    lgradp[GpVars::gpx].setVal(0.);
    lgradp[GpVars::gpy].define(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, 0);
    lgradp[GpVars::gpy].setVal(0.);
    lgradp[GpVars::gpz].define(lev_new[Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, 0);
    lgradp[GpVars::gpz].setVal(0.);

    MultiFab p0(base_state[lev], make_alias, BaseState::p0_comp, 1);

    make_gradp_pert(lev, solverChoice, geom[lev], lev_new, p0,
                    *z_phys_nd[lev].get(), *z_phys_cc[lev].get(), mapfac[lev],
                    get_eb(lev), lgradp);

    // *******************************************************************************
    // Compute the buoyancy term exactly as the dycore does
    // *******************************************************************************

    MultiFab qt(ba, dm, 1, 1);
    qt.setVal(0.);
    int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
    if (solverChoice.moisture_type != MoistureType::None) {
        make_qt(lev_new[Vars::cons], qt, n_qstate_into_total);
    }

    //
    // NOTE: we must fill one ghost cell of S_prim here because make_buoyancy
    //       reads cell_prim(i,j,k-1) at the lower z face of every box
    //
    MultiFab S_prim(ba, dm, lev_new[Vars::cons].nComp()-1, 1);
    cons_to_prim(lev_new[Vars::cons], S_prim, 1);

    // Initialize to zero because buoyancy is not defined on the faces at the top and bottom of the domain
    MultiFab buoyancy(lgradp[GpVars::gpz].boxArray(), lgradp[GpVars::gpz].DistributionMap(), 1, 0);
    buoyancy.setVal(0.);

    make_buoyancy(lev, lev_new, S_prim, qt, buoyancy, geom[lev], solverChoice, base_state[lev],
                  n_qstate_into_total, get_eb(lev), solverChoice.anelastic[lev]);

    // *******************************************************************************
    // Now form the combined term that appears in the z-momentum equation:
    //      -d(p - p0)/dz + buoyancy
    // If the initial state is in HSE and the buoyancy is computed from (rho - rho0)
    //      (i.e. buoyancy_type = 1) then this is identically zero
    // *******************************************************************************

    MultiFab combined(lgradp[GpVars::gpz].boxArray(), lgradp[GpVars::gpz].DistributionMap(), 1, 0);
    combined.setVal(0.);

    MultiFab::Copy    (combined, buoyancy          , 0, 0, 1, 0);
    MultiFab::Subtract(combined, lgradp[GpVars::gpz], 0, 0, 1, 0);

    // Don't include the top and bottom faces of the domain, where neither term is defined
    Box zface_domain = surroundingNodes(geom[lev].Domain(), 2);
    int klo = zface_domain.smallEnd(2);
    int khi = zface_domain.bigEnd(2);
    zface_domain.growLo(2,-1);
    zface_domain.growHi(2,-1);

    int comp = 0;
    Real min_val = combined.min(zface_domain,comp);
    Real max_val = combined.max(zface_domain,comp);

    Print() << " " << std::endl;
    if (max_val != zero || min_val != zero) {
        IntVect min_loc = combined.minIndex(comp);
        IntVect max_loc = combined.maxIndex(comp);
        Print() << "Min/max value of -dp'/dz + buoyancy are " << min_val << " " << max_val;
        if (min_loc[2] != klo && min_loc[2] != khi) Print() << " with min at face " << min_loc;
        if (max_loc[2] != klo && max_loc[2] != khi) Print() << " with max at face " << max_loc;
        Print() << std::endl;
    } else {
        Print() << "Min/max value of -dp'/dz + buoyancy are zero " << std::endl;
    }
    Print() << " " << std::endl;
}
