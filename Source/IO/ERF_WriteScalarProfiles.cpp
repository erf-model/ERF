#include <iomanip>

#include "ERF.H"
#include "ERF_Derive.H"

using namespace amrex;

/**
 * Computes the integrated quantities on the grid such as the
 * total scalar and total mass quantities. Prints and writes to output file.
 *
 * @param time Current time
 */
void
ERF::sum_integrated_quantities (Real time)
{
    BL_PROFILE("ERF::sum_integrated_quantities()");

    if (verbose <= 0)
      return;

    // Single level sum
    Real mass_sl;

    // Multilevel sums
    Real mass_ml = 0.0;
    Real rhth_ml = 0.0;
    Real scal_ml = 0.0;
    Real mois_ml = 0.0;

    bool local = true;

    auto& mfx0 = *mapfac[0][MapFacType::m_x];
    auto& mfy0 = *mapfac[0][MapFacType::m_x];
    auto&  dJ0 = *detJ_cc[0];

    mass_sl = volWgtSumMF(0,vars_new[0][Vars::cons],Rho_comp,dJ0,mfx0,mfy0,false,local);

    for (int lev = 0; lev <= finest_level; lev++) {
        auto& mfx = *mapfac[lev][MapFacType::m_x];
        auto& mfy = *mapfac[lev][MapFacType::m_x];
        auto&  dJ = *detJ_cc[lev];
        mass_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],Rho_comp,dJ,mfx,mfy,true);
    }

    Real rhth_sl = volWgtSumMF(0,vars_new[0][Vars::cons], RhoTheta_comp,dJ0,mfx0,mfy0,false);
    Real scal_sl = volWgtSumMF(0,vars_new[0][Vars::cons],RhoScalar_comp,dJ0,mfx0,mfy0,false);
    Real mois_sl = 0.0;
    if (solverChoice.moisture_type != MoistureType::None) {
        int n_qstate_moist = micro->Get_Qstate_Moist_Size();
        for (int qoff(0); qoff<n_qstate_moist; ++qoff) {
            mois_sl += volWgtSumMF(0,vars_new[0][Vars::cons],RhoQ1_comp+qoff,dJ0,mfx0,mfy0,false);
        }
    }

    for (int lev = 0; lev <= finest_level; lev++) {
        auto& mfx = *mapfac[lev][MapFacType::m_x];
        auto& mfy = *mapfac[lev][MapFacType::m_x];
        auto&  dJ = *detJ_cc[lev];
        rhth_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons], RhoTheta_comp,dJ,mfx,mfy,true);
        scal_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],RhoScalar_comp,dJ,mfx,mfy,true);
        if (solverChoice.moisture_type != MoistureType::None) {
            int n_qstate_moist = micro->Get_Qstate_Moist_Size();
            for (int qoff(0); qoff<n_qstate_moist; ++qoff) {
                mois_ml += volWgtSumMF(lev,vars_new[lev][Vars::cons],RhoQ1_comp+qoff,dJ,mfx,mfy,false);
            }
        }
    }

    Gpu::HostVector<Real> h_avg_ustar; h_avg_ustar.resize(1);
    Gpu::HostVector<Real> h_avg_tstar; h_avg_tstar.resize(1);
    Gpu::HostVector<Real> h_avg_olen; h_avg_olen.resize(1);
    if ((m_SurfaceLayer != nullptr) && (NumDataLogs() > 0)) {
        Box domain = geom[0].Domain();
        int zdir = 2;
        h_avg_ustar = sumToLine(*m_SurfaceLayer->get_u_star(0),0,1,domain,zdir);
        h_avg_tstar = sumToLine(*m_SurfaceLayer->get_t_star(0),0,1,domain,zdir);
        h_avg_olen  = sumToLine(*m_SurfaceLayer->get_olen(0)  ,0,1,domain,zdir);

        // Divide by the total number of cells we are averaging over
        Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
        h_avg_ustar[0] /= area_z;
        h_avg_tstar[0] /= area_z;
        h_avg_olen[0]  /= area_z;

    } else {
        h_avg_ustar[0] = 0.;
        h_avg_tstar[0] = 0.;
        h_avg_olen[0]  = 0.;
    }

    const int nfoo = 8;
    Real foo[nfoo] = {mass_sl,rhth_sl,scal_sl,mois_sl,mass_ml,rhth_ml,scal_ml,mois_ml};
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(
        foo, nfoo, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
        int i = 0;
        mass_sl = foo[i++];
        rhth_sl = foo[i++];
        scal_sl = foo[i++];
        mois_sl = foo[i++];
        mass_ml = foo[i++];
        rhth_ml = foo[i++];
        scal_ml = foo[i++];
        mois_ml = foo[i++];

        Print() << '\n';
        Print() << "TIME= " << std::setw(datwidth) << std::setprecision(timeprecision) << std::left << time << '\n';
        if (finest_level ==  0) {
#if 1
           Print() << " MASS       = " << mass_sl << '\n';
#else
           Print() << " PERT MASS  = " << mass_sl << '\n';
#endif
           Print() << " RHO THETA  = " << rhth_sl << '\n';
           Print() << " RHO SCALAR = " << scal_sl << '\n';
           Print() << " RHO QTOTAL = " << mois_sl << '\n';
        } else {
#if 1
           Print() << " MASS       SL/ML = " << mass_sl << " " << mass_ml << '\n';
#else
           Print() << " PERT MASS  SL/ML = " << mass_sl << " " << mass_ml << '\n';
#endif
           Print() << " RHO THETA  SL/ML = " << rhth_sl << " " << rhth_ml << '\n';
           Print() << " RHO SCALAR SL/ML = " << scal_sl << " " << scal_ml << '\n';
           Print() << " RHO QTOTAL SL/ML = " << mois_sl << " " << mois_ml << '\n';
        }

        // The first data log only holds scalars
        if (NumDataLogs() > 0)
        {
            int n_d = 0;
            std::ostream& data_log1 = DataLog(n_d);
            if (data_log1.good()) {
                if (time == 0.0) {
                    data_log1 << std::setw(datwidth) << "          time";
                    data_log1 << std::setw(datwidth) << "          u_star";
                    data_log1 << std::setw(datwidth) << "          t_star";
                    data_log1 << std::setw(datwidth) << "          olen";
                    data_log1 << std::endl;
                } // time = 0

              // Write the quantities at this time
              data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time;
              data_log1 << std::setw(datwidth) << std::setprecision(datprecision)  << h_avg_ustar[0];
              data_log1 << std::setw(datwidth) << std::setprecision(datprecision)  << h_avg_tstar[0];
              data_log1 << std::setw(datwidth) << std::setprecision(datprecision)  << h_avg_olen[0];
              data_log1 << std::endl;
            } // if good
        } // loop over i
      } // if IOProcessor
#ifdef AMREX_LAZY
    });
#endif

    // This is just an alias for convenience
    int lev = 0;
    if (NumSamplePointLogs() > 0 && NumSamplePoints() > 0) {
        for (int i = 0; i < NumSamplePoints(); ++i)
        {
            sample_points(lev, time, SamplePoint(i), vars_new[lev][Vars::cons]);
        }
    }
    if (NumSampleLineLogs() > 0 && NumSampleLines() > 0) {
        for (int i = 0; i < NumSampleLines(); ++i)
        {
            sample_lines(lev, time, SampleLine(i), vars_new[lev][Vars::cons]);
        }
    }
}

void
ERF::sum_derived_quantities (Real time)
{
    if (verbose <= 0 || NumDerDataLogs() <= 0) return;

    int lev = 0;

    AMREX_ALWAYS_ASSERT(lev == 0);

    auto& mfx0 = *mapfac[0][MapFacType::m_x];
    auto& mfy0 = *mapfac[0][MapFacType::m_x];
    auto&  dJ0 = *detJ_cc[0];

    // ************************************************************************
    // WARNING: we are not filling ghost cells other than periodic outside the domain
    // ************************************************************************

    MultiFab mf_cc_vel(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
    mf_cc_vel.setVal(0.); // We just do this to avoid uninitialized values

    // Average all three components of velocity (on faces) to the cell center
    average_face_to_cellcenter(mf_cc_vel,0,
                               Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],
                                                        &vars_new[lev][Vars::yvel],
                                                        &vars_new[lev][Vars::zvel]});
    mf_cc_vel.FillBoundary(geom[lev].periodicity());

    if (!geom[lev].isPeriodic(0) || !geom[lev].isPeriodic(1) || !geom[lev].isPeriodic(2)) {
        amrex::Warning("Ghost cells outside non-periodic physical boundaries are not filled -- vel set to 0 there");
    }

    MultiFab r_wted_magvelsq(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));
    MultiFab unwted_magvelsq(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));
    MultiFab     enstrophysq(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
    MultiFab        theta_mf(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(unwted_magvelsq, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto& src_fab = mf_cc_vel[mfi];

        auto& dest1_fab = unwted_magvelsq[mfi];
        derived::erf_dermagvelsq(bx, dest1_fab, 0, 1, src_fab, Geom(lev), t_new[0], nullptr, lev);

        auto& dest2_fab = enstrophysq[mfi];
        derived::erf_derenstrophysq(bx, dest2_fab, 0, 1, src_fab, Geom(lev), t_new[0], nullptr, lev);
    }

    // Copy the MF holding 1/2(u^2 + v^2 + w^2) into the MF that will hold 1/2 rho (u^2 + v^2 + w^2)d
    MultiFab::Copy(r_wted_magvelsq, unwted_magvelsq, 0, 0, 1, 0);

    // Multiply the MF holding 1/2(u^2 + v^2 + w^2) by rho to get  1/2 rho (u^2 + v^2 + w^2)
    MultiFab::Multiply(r_wted_magvelsq, vars_new[lev][Vars::cons], 0, 0, 1, 0);

    // Copy the MF holding (rho theta) into "theta_mf"
    MultiFab::Copy(theta_mf, vars_new[lev][Vars::cons], RhoTheta_comp, 0, 1, 0);

    // Divide (rho theta) by rho to get theta in the MF "theta_mf"
    MultiFab::Divide(theta_mf, vars_new[lev][Vars::cons], Rho_comp, 0, 1, 0);

    Real  unwted_avg = volWgtSumMF(lev, unwted_magvelsq, 0, dJ0, mfx0, mfy0, false);
    Real  r_wted_avg = volWgtSumMF(lev, r_wted_magvelsq, 0, dJ0, mfx0, mfy0, false);
    Real enstrsq_avg = volWgtSumMF(lev, enstrophysq,     0, dJ0, mfx0, mfy0, false);
    Real   theta_avg = volWgtSumMF(lev, theta_mf,        0, dJ0, mfx0, mfy0, false);

    // Get volume including terrain (consistent with volWgtSumMF routine)
    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    auto const& dx = geom[lev].CellSizeArray();
    Real cell_vol  = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        MultiFab::Multiply(volume, *detJ_cc[lev], 0, 0, 1, 0);
    }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(volume, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx  = mfi.tilebox();
        auto dst        = volume.array(mfi);
        const auto& mfx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const auto& mfy = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst(i,j,k) /= (mfx(i,j,0)*mfy(i,j,0));
        });
    }
    Real vol = volume.sum();

     unwted_avg /= vol;
     r_wted_avg /= vol;
    enstrsq_avg /= vol;
      theta_avg /= vol;

    const int nfoo = 4;
    Real foo[nfoo] = {unwted_avg,r_wted_avg,enstrsq_avg,theta_avg};
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(
        foo, nfoo, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
        int i = 0;
        unwted_avg  = foo[i++];
        r_wted_avg  = foo[i++];
        enstrsq_avg = foo[i++];
          theta_avg = foo[i++];

        std::ostream& data_log_der = DerDataLog(0);

        if (time == 0.0) {
            data_log_der << std::setw(datwidth) << "          time";
            data_log_der << std::setw(datwidth) << "        ke_den";
            data_log_der << std::setw(datwidth) << "         velsq";
            data_log_der << std::setw(datwidth) << "     enstrophy";
            data_log_der << std::setw(datwidth) << "    int_energy";
            data_log_der << std::endl;
        }
        data_log_der << std::setw(datwidth) << std::setprecision(timeprecision) << time;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  <<  unwted_avg;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  <<  r_wted_avg;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  << enstrsq_avg;
        data_log_der << std::setw(datwidth) << std::setprecision(datprecision)  <<   theta_avg;
        data_log_der << std::endl;

      } // if IOProcessor
#ifdef AMREX_LAZY
    }
#endif
}

void
ERF::sum_energy_quantities (Real time)
{
    if ( (verbose <= 0) || (tot_e_datalog.size() < 1) ) { return; }

    int lev = 0;

    auto& mfx0 = *mapfac[0][MapFacType::m_x];
    auto& mfy0 = *mapfac[0][MapFacType::m_x];
    auto&  dJ0 = *detJ_cc[0];

    AMREX_ALWAYS_ASSERT(lev == 0);

    bool local = true;

    // ************************************************************************
    // WARNING: we are not filling ghost cells other than periodic outside the domain
    // ************************************************************************

    MultiFab mf_cc_vel(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(1,1,1));
    mf_cc_vel.setVal(0.); // We just do this to avoid uninitialized values

    // Average all three components of velocity (on faces) to the cell center
    average_face_to_cellcenter(mf_cc_vel,0,
                               Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],
                                                        &vars_new[lev][Vars::yvel],
                                                        &vars_new[lev][Vars::zvel]});
    mf_cc_vel.FillBoundary(geom[lev].periodicity());

    if (!geom[lev].isPeriodic(0) || !geom[lev].isPeriodic(1) || !geom[lev].isPeriodic(2)) {
        amrex::Warning("Ghost cells outside non-periodic physical boundaries are not filled -- vel set to 0 there");
    }

    MultiFab tot_mass  (grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));
    MultiFab tot_energy(grids[lev], dmap[lev], AMREX_SPACEDIM, IntVect(0,0,0));

    auto const& dx = geom[lev].CellSizeArray();
    bool is_moist = (solverChoice.moisture_type != MoistureType::None);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(tot_mass, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const Array4<Real>& cc_vel_arr     = mf_cc_vel.array(mfi);
        const Array4<Real>& tot_mass_arr   = tot_mass.array(mfi);
        const Array4<Real>& tot_energy_arr = tot_energy.array(mfi);
        const Array4<const Real>& cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
        const Array4<const Real>& z_arr    = (z_phys_nd[lev]) ? z_phys_nd[lev]->const_array(mfi) :
                                                                Array4<const Real>{};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real Qv   = (is_moist) ? cons_arr(i,j,k,RhoQ1_comp) : 0.0;
            Real Qc   = (is_moist) ? cons_arr(i,j,k,RhoQ2_comp) : 0.0;
            Real Qt   = Qv + Qc;
            Real Rhod = cons_arr(i,j,k,Rho_comp);
            Real Rhot = Rhod * (1.0 + Qt);
            Real Temp = getTgivenRandRTh(Rhod, cons_arr(i,j,k,RhoTheta_comp), Qv);
            Real TKE  = 0.5 * ( cc_vel_arr(i,j,k,0)*cc_vel_arr(i,j,k,0)
                              + cc_vel_arr(i,j,k,1)*cc_vel_arr(i,j,k,1)
                              + cc_vel_arr(i,j,k,2)*cc_vel_arr(i,j,k,2) );
            Real zval = (z_arr) ? z_arr(i,j,k) : Real(k)*dx[2];

            Real Cv   = Cp_d - R_d;
            Real Cvv  = Cp_v - R_v;
            Real Cpv  = Cp_v;

            tot_mass_arr(i,j,k)   = Rhot;
            tot_energy_arr(i,j,k) = Rhod * ( (Cv + Cvv*Qv + Cpv*Qc)*Temp - L_v*Qc
                                           + (1.0 + Qt)*TKE + (1.0 + Qt)*CONST_GRAV*zval );

        });

    }

    Real  tot_mass_avg   = volWgtSumMF(lev, tot_mass  , 0, dJ0, mfx0, mfy0, false, local);
    Real  tot_energy_avg = volWgtSumMF(lev, tot_energy, 0, dJ0, mfx0, mfy0, false, local);

    // Get volume including terrain (consistent with volWgtSumMF routine)
    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    Real cell_vol  = dx[0]*dx[1]*dx[2];
    volume.setVal(cell_vol);
    if (SolverChoice::mesh_type != MeshType::ConstantDz) {
        MultiFab::Multiply(volume, *detJ_cc[lev], 0, 0, 1, 0);
    }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(volume, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx  = mfi.tilebox();
        auto dst        = volume.array(mfi);
        const auto& mfx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const auto& mfy = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst(i,j,k) /= (mfx(i,j,0)*mfy(i,j,0));
        });
    }
    Real vol = volume.sum();

    // Divide by the volume
    tot_mass_avg   /= vol;
    tot_energy_avg /= vol;

    const int nfoo = 2;
    Real foo[nfoo] = {tot_mass_avg,tot_energy_avg};
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
    ParallelDescriptor::ReduceRealSum(
        foo, nfoo, ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor()) {
        int i = 0;
        tot_mass_avg   = foo[i++];
        tot_energy_avg = foo[i++];

        std::ostream& data_log_energy = *tot_e_datalog[0];

        if (time == 0.0) {
            data_log_energy << std::setw(datwidth) << "          time";
            data_log_energy << std::setw(datwidth) << "      tot_mass";
            data_log_energy << std::setw(datwidth) << "    tot_energy";
            data_log_energy << std::endl;
        }
        data_log_energy << std::setw(datwidth) << std::setprecision(timeprecision) << time;
        data_log_energy << std::setw(datwidth) << std::setprecision(datprecision)  << tot_mass_avg;
        data_log_energy << std::setw(datwidth) << std::setprecision(datprecision)  << tot_energy_avg;
        data_log_energy << std::endl;

      } // if IOProcessor
#ifdef AMREX_LAZY
    }
#endif
}

Real
ERF::cloud_fraction (Real /*time*/)
{
    BL_PROFILE("ERF::cloud_fraction()");

    int lev = 0;
    // This holds all of qc
    MultiFab qc(vars_new[lev][Vars::cons],make_alias,RhoQ2_comp,1);

    int direction = 2; // z-direction
    Box const& domain = geom[lev].Domain();

    auto const& qc_arr = qc.const_arrays();

    // qc_2d is an BaseFab<int> holding the max value over the column
    auto qc_2d = ReduceToPlane<ReduceOpMax,int>(direction, domain, qc,
         [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> int
         {
             if (qc_arr[box_no](i,j,k) > 0) {
                 return 1;
             } else {
                 return 0;
             }
         });

    auto* p = qc_2d.dataPtr();

    Long numpts = qc_2d.numPts();

    AMREX_ASSERT(numpts < Long(std::numeric_limits<int>::max));

#if 1
    if (ParallelDescriptor::UseGpuAwareMpi()) {
        ParallelDescriptor::ReduceIntMax(p,static_cast<int>(numpts));
    } else {
        Gpu::PinnedVector<int> hv(numpts);
        Gpu::copyAsync(Gpu::deviceToHost, p, p+numpts, hv.data());
        Gpu::streamSynchronize();
        ParallelDescriptor::ReduceIntMax(hv.data(),static_cast<int>(numpts));
        Gpu::copyAsync(Gpu::hostToDevice, hv.data(), hv.data()+numpts, p);
    }

    // Sum over component 0
    Long num_cloudy = qc_2d.template sum<RunOn::Device>(0);

#else
    //
    // We need this if we allow domain decomposition in the vertical
    //    but for now we leave it commented out
    //
    Long num_cloudy = Reduce::Sum<Long>(numpts,
         [=] AMREX_GPU_DEVICE (Long i) -> Long {
             if (p[i] == 1) {
                 return 1;
             } else {
                 return 0;
             }
         });
    ParallelDescriptor::ReduceLongSum(num_cloudy);
#endif

    Real num_total = qc_2d.box().d_numPts();

    Real cloud_frac = num_cloudy / num_total;

    return cloud_frac;
}

/**
 * Utility function for sampling MultiFab data at a specified cell index.
 *
 * @param lev Level for the associated MultiFab data
 * @param time Current time
 * @param cell IntVect containing the indexes for the cell where we want to sample
 * @param mf MultiFab from which we wish to sample data
 */
void
ERF::sample_points (int /*lev*/, Real time, IntVect cell, MultiFab& mf)
{
    int ifile = 0;

    //
    // Sample the data at a single point in space
    //
    int ncomp = mf.nComp();
    Vector<Real> my_point = get_cell_data(mf, cell);

    if (!my_point.empty()) {

        // HERE DO WHATEVER YOU WANT TO THE DATA BEFORE WRITING

        std::ostream& sample_log = SamplePointLog(ifile);
        if (sample_log.good()) {
          sample_log << std::setw(datwidth) << time;
          for (int i = 0; i < ncomp; ++i)
          {
              sample_log << std::setw(datwidth) << my_point[i];
          }
          sample_log << std::endl;
        } // if good
    } // only write from processor that holds the cell
}

/**
 * Utility function for sampling data along a line along the z-dimension
 * at the (x,y) indices specified and writes it to an output file.
 *
 * @param lev Current level
 * @param time Current time
 * @param cell IntVect containing the x,y-dimension indices to sample along z
 * @param mf MultiFab from which we sample the data
 */
void
ERF::sample_lines (int lev, Real time, IntVect cell, MultiFab& mf)
{
    int ifile = 0;

    const int ncomp = mf.nComp(); // cell-centered state vars

    MultiFab mf_vels(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
    average_face_to_cellcenter(mf_vels, 0,
                               Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});

    //
    // Sample the data at a line (in direction "dir") in space
    // In this case we sample in the vertical direction so dir = 2
    // The "k" value of "cell" is ignored
    //
    int dir = 2;
    MultiFab my_line       = get_line_data(mf,              dir, cell);
    MultiFab my_line_vels  = get_line_data(mf_vels,         dir, cell);
    MultiFab my_line_tau11 = get_line_data(*Tau[lev][TauType::tau11], dir, cell);
    MultiFab my_line_tau12 = get_line_data(*Tau[lev][TauType::tau12], dir, cell);
    MultiFab my_line_tau13 = get_line_data(*Tau[lev][TauType::tau13], dir, cell);
    MultiFab my_line_tau22 = get_line_data(*Tau[lev][TauType::tau22], dir, cell);
    MultiFab my_line_tau23 = get_line_data(*Tau[lev][TauType::tau23], dir, cell);
    MultiFab my_line_tau33 = get_line_data(*Tau[lev][TauType::tau33], dir, cell);

    for (MFIter mfi(my_line, false); mfi.isValid(); ++mfi)
    {
        // HERE DO WHATEVER YOU WANT TO THE DATA BEFORE WRITING

        std::ostream& sample_log = SampleLineLog(ifile);
        if (sample_log.good()) {
          sample_log << std::setw(datwidth) << std::setprecision(datprecision) << time;
          const auto& my_line_arr = my_line[0].const_array();
          const auto& my_line_vels_arr = my_line_vels[0].const_array();
          const auto& my_line_tau11_arr = my_line_tau11[0].const_array();
          const auto& my_line_tau12_arr = my_line_tau12[0].const_array();
          const auto& my_line_tau13_arr = my_line_tau13[0].const_array();
          const auto& my_line_tau22_arr = my_line_tau22[0].const_array();
          const auto& my_line_tau23_arr = my_line_tau23[0].const_array();
          const auto& my_line_tau33_arr = my_line_tau33[0].const_array();
          const Box&  my_box = my_line[0].box();
          const int klo = my_box.smallEnd(2);
          const int khi = my_box.bigEnd(2);
          int i = cell[0];
          int j = cell[1];
          for (int n = 0; n < ncomp; n++) {
              for (int k = klo; k <= khi; k++) {
                  sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_arr(i,j,k,n);
              }
          }
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
              for (int k = klo; k <= khi; k++) {
                  sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_vels_arr(i,j,k,n);
              }
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau11_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau12_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau13_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau22_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau23_arr(i,j,k);
          }
          for (int k = klo; k <= khi; k++) {
              sample_log << std::setw(datwidth) << std::setprecision(datprecision) << my_line_tau33_arr(i,j,k);
          }
          sample_log << std::endl;
        } // if good
    } // mfi
}

/**
 * Helper function which uses the current step number, time, and timestep to
 * determine whether it is time to take an action specified at every interval
 * of timesteps.
 *
 * @param nstep Timestep number
 * @param time Current time
 * @param dtlev Timestep for the current level
 * @param action_interval Interval in number of timesteps for taking action
 * @param action_per Interval in simulation time for taking action
 */
bool
ERF::is_it_time_for_action (int nstep, Real time, Real dtlev, int action_interval, Real action_per)
{
  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > 0.0) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}
