/**
 * \file ERF_Diagnostics.cpp
 */

#include "ERF.H"
#include "ERF_SrcHeaders.H"
#include "ERF_BuoyancyUtils.H"
#include "ERF_Utils.H"

using namespace amrex;

void 
ERF::compute_max_pressure_gradient_diagnostic(int lev)
{
    auto& lev_new = vars_new[lev];

    Vector<MultiFab> gradp_temp;  gradp_temp.resize(AMREX_SPACEDIM);
    gradp_temp[0].define(vars_new[lev][Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, 0);
    gradp_temp[0].setVal(0.);
    gradp_temp[1].define(vars_new[lev][Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, 0);
    gradp_temp[1].setVal(0.);
    gradp_temp[2].define(vars_new[lev][Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, 0);
    gradp_temp[2].setVal(0.);

    int comp = 0;

    // ******************************************************************************* 
    // First compute for base state pressure
    // ******************************************************************************* 

    MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp , 1);
    MultiFab p_hse(base_state[lev], make_alias, BaseState::p0_comp , 1);

    compute_gradp(p_hse, geom[lev], *z_phys_nd[lev].get(), *z_phys_cc[lev].get(), mapfac[lev],
                  get_eb(lev), gradp_temp, solverChoice);

    Real min_gpx = gradp_temp[0].min(comp);
    Real max_gpx = gradp_temp[0].max(comp);
    if (max_gpx != zero || min_gpx != zero) {
        Print() << "Min/Max value of x-gradient of base state pressure are " << min_gpx << " " << max_gpx <<
                    " and occur at faces " << gradp_temp[0].minIndex(comp) << " and " 
                                           << gradp_temp[0].maxIndex(comp) << std::endl;
    } else {
        Print() << "Min/max value of x-gradient of base state pressure are zero " << std::endl;
    }

    Real min_gpy = gradp_temp[1].min(comp);
    Real max_gpy = gradp_temp[1].max(comp);
    if (max_gpy != zero || min_gpy != zero) {
        Print() << "Min/max value of y-gradient of base state pressure are " << min_gpy << " " << max_gpy <<
                    " and occur at faces " << gradp_temp[1].minIndex(comp) << " and " 
                                           << gradp_temp[1].maxIndex(comp) << std::endl;
    } else {
        Print() << "Min/max value of y-gradient of base state pressure are zero " << std::endl;
    }

    const Real grav = solverChoice.gravity;
    for (MFIter mfi(gradp_temp[2]); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox(); bx.growHi(2,-1);
        if (bx.smallEnd(2) == 0) bx.growLo(2,-1);
        auto        gpz_arr  = gradp_temp[2].array(mfi);
        auto const rhse_arr  = r_hse.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            gpz_arr(i,j,k) -= grav * myhalf * (rhse_arr(i,j,k) + rhse_arr(i,j,k-1));
        });
    }

    Real min_gpz = gradp_temp[2].min(comp);
    Real max_gpz = gradp_temp[2].max(comp);
    if (max_gpz > zero) {
        Print() << "Min/max value of dp0/dz - rho0*g  are " << min_gpz << " " << max_gpz <<
                    " and occur at faces " << gradp_temp[2].minIndex(comp) << " and " 
                                           << gradp_temp[2].maxIndex(comp) << std::endl;
    } else {
        Print() << "Min/max value of dp0/dz - rho0*g  are zero " << std::endl;
    }
    Print() << " " << std::endl;

   // ******************************************************************************* 
   // Now compute for full (moist) pressure
   // ******************************************************************************* 

    MultiFab p(p_hse.boxArray(), p_hse.DistributionMap(), 1, 1);
    MultiFab rho(lev_new[Vars::cons], make_alias, Rho_comp , 1);

    if (solverChoice.moisture_type != MoistureType::None) {

        for (MFIter mfi(rho); mfi.isValid(); ++mfi)
        {
            Box gbx = mfi.tilebox();
            gbx.grow(IntVect(1,1,1));
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

        min_gpx = gradp_temp[0].min(comp);
        max_gpx = gradp_temp[0].max(comp);
        if (max_gpx != zero || min_gpx != zero) {
            Print() << "Min/Max value of x-gradient of full (moist) pressure are " << min_gpx << " " << max_gpx <<
                        " and occur at faces " << gradp_temp[0].minIndex(comp) << " and " 
                                               << gradp_temp[0].maxIndex(comp) << std::endl;
        } else {
            Print() << "Min/max value of x-gradient of full (moist) pressure are zero " << std::endl;
        }

        min_gpy = gradp_temp[1].min(comp);
        max_gpy = gradp_temp[1].max(comp);
        if (max_gpy != zero || min_gpy != zero) {
            Print() << "Min/Max value of y-gradient of full (moist) pressure are " << min_gpy << " " << max_gpy <<
                        " and occur at faces " << gradp_temp[1].minIndex(comp) << " and " 
                                               << gradp_temp[1].maxIndex(comp) << std::endl;
        } else {
            Print() << "Min/max value of y-gradient of full (moist) pressure are zero " << std::endl;
        }

        MultiFab qt(rho.boxArray(), rho.DistributionMap(), 1, 1);
        int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
        make_qt(lev_new[IntVars::cons], qt, n_qstate_into_total);

        const Real grav = solverChoice.gravity;
        for (MFIter mfi(gradp_temp[2]); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.validbox(); bx.growHi(2,-1);
            if (bx.smallEnd(2) == 0) bx.growLo(2,-1);
            auto      gpz_arr  = gradp_temp[2].array(mfi);
            auto const  s_arr  = lev_new[Vars::cons].const_array(mfi);
            auto const  r_arr  = rho.const_array(mfi);
            auto const qt_arr  = qt.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                gpz_arr(i,j,k) -= buoyancy_rhopert(i,j,k,grav,r_arr,s_arr,qt_arr);
            });
        }

        min_gpz = gradp_temp[2].min(comp);
        max_gpz = gradp_temp[2].max(comp);
        if (max_gpz != zero || min_gpz != zero) {
            Print() << "Min/max value of moist dp/dz - rho*g  are " << min_gpz << " " << max_gpz <<
                        " and occur at faces " << gradp_temp[2].minIndex(comp) << " and " 
                                               << gradp_temp[2].maxIndex(comp) << std::endl;
        } else {
            Print() << "Min/max value of moist dp/dz - rho*g  are zero " << std::endl;
        }
        Print() << " " << std::endl;
    } // if moist

    // ******************************************************************************* 
    // Now compute for full (dry) pressure
    // ******************************************************************************* 

    for ( MFIter mfi(p); mfi.isValid(); ++mfi)
    {
        Box gbx = mfi.tilebox();
        gbx.grow(IntVect(1,1,1));
        if (gbx.smallEnd(2) < 0) gbx.setSmall(2,0);

        const Array4<const Real>& cell_data = lev_new[Vars::cons].array(mfi);
        const Array4<      Real>& pp_arr = p.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real qv_for_p = zero;
            pp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p);
        });
    }

    compute_gradp(p, geom[lev], *z_phys_nd[lev].get(), *z_phys_cc[lev].get(), mapfac[lev],
                  get_eb(lev), gradp_temp, solverChoice);

    min_gpx = gradp_temp[0].min(comp);
    max_gpx = gradp_temp[0].max(comp);
    if (max_gpx != zero || min_gpx != zero) {
        Print() << "Min/max value of x-gradient of full (dry) pressure are " << min_gpx << " " << max_gpx <<
                    " and occur at faces " << gradp_temp[0].minIndex(comp) << " and " 
                                           << gradp_temp[0].maxIndex(comp) << std::endl;
    } else {
        Print() << "Min/max value of x-gradient of full (dry) pressure are zero " << std::endl;
    }

    min_gpy = gradp_temp[1].min(comp);
    max_gpy = gradp_temp[1].max(comp);
    if (max_gpy != zero || min_gpy != zero) {
        Print() << "Min/max value of y-gradient of full (dry) pressure are " << min_gpy << " " << max_gpy <<
                    " and occur at faces " << gradp_temp[1].minIndex(comp) << " and " 
                                           << gradp_temp[1].maxIndex(comp) << std::endl;
    } else {
        Print() << "Min/max value of y-gradient of full (dry) pressure are zero " << std::endl;
    }


    for (MFIter mfi(gradp_temp[2]); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox(); bx.growHi(2,-1); 
        if (bx.smallEnd(2) == 0) bx.growLo(2,-1);
        auto     gpz_arr  = gradp_temp[2].array(mfi);
        auto const r_arr  = rho.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            gpz_arr(i,j,k) -= grav * myhalf * (r_arr(i,j,k) + r_arr(i,j,k-1));
        });
    }

    min_gpz = gradp_temp[2].min(comp);
    max_gpz = gradp_temp[2].max(comp);
    if (max_gpz != zero || min_gpz != zero) {
        Print() << "Min/max value of dry dp/dz - rho_m*g  are " << min_gpz << " " << max_gpz <<
                    " and occur at faces " << gradp_temp[2].minIndex(comp) << " and " 
                                           << gradp_temp[2].maxIndex(comp) << std::endl;
    } else {
        Print() << "Min/max value of dry dp/dz - rho*g  are zero " << std::endl;
    }
    Print() << " " << std::endl;
}
