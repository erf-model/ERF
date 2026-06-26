/**
 * \file ERF_Diagnostics.cpp
 */

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
    Real max_diff = dp.max(0);
    if (max_diff > 1.e-8) {
        IntVect max_loc = dp.maxIndex(0);
        Print() << "Max value of |p_hse - p_eos| is " << max_diff << std::endl;
        Print() << " with max in cell " << max_loc << std::endl;
        Abort("Base state violates EOS ");
    } else {
        Print() << "Max value of |p_hse - p_eos| is less than 1e-8" << std::endl;
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
            Box bx = mfi.validbox(); bx.growHi(2,-1);
            if (bx.smallEnd(2) == 0) bx.growLo(2,-1);
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
            Box bx = mfi.validbox(); bx.growHi(2,-1);
            if (bx.smallEnd(2) == 0) bx.growLo(2,-1);
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

#ifdef AMREX_USE_FLOAT
    Real tol = 1.e-4;
#else
    Real tol = 1.e-8;
#endif

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
                Box bx = mfi.validbox(); bx.growHi(2,-1);
                if (bx.smallEnd(2) == 0) bx.growLo(2,-1);
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
                Print() << std::endl;
            } else {
                Print() << "Min/max value of moist dp/dz + rho_m*|g|  are zero " << std::endl;
            }
            Print() << " " << std::endl;
        } // if moist
    } // if !anelastic
}
