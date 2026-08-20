#if defined(ERF_USE_NETCDF)

#include <ERF_SrcHeaders.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in] S_data  current solution
 * @param[in] source  source terms for moist conserved variables
 * @param[in] dt      slow time step
 * @param[in] old_stage_time_total
 * @param[in] start_bdy_time
 * @param[in] final_bdy_time
 * @param[in] bdy_time_interval
 * @param[in] bdy_factor
 * @param[in] width
 * @param[in] geom   Container for geometric information
 * @param[in] bdy_data_xlo boundary data on interior of low x-face
 * @param[in] bdy_data_xhi boundary data on interior of high x-face
 * @param[in] bdy_data_ylo boundary data on interior of low y-face
 * @param[in] bdy_data_yhi boundary data on interior of high y-face
 * @param[in] m_r2d        boundary data read in using ReadBndryPlane
 */

void add_moist_nudging_terms (const MultiFab& S_data,
                                    MultiFab & source,
                              const int n_qstate,
                              const double& dt,
                              const double& time,
                              const double& start_bdy_time,
                              const double& final_bdy_time,
                              const double& bdy_time_interval,
                              const Real& nudge_factor,
                              const int width,
                              const Geometry& geom,
                              Vector<Vector<FArrayBox>>& bdy_data_xlo,
                              Vector<Vector<FArrayBox>>& bdy_data_xhi,
                              Vector<Vector<FArrayBox>>& bdy_data_ylo,
                              Vector<Vector<FArrayBox>>& bdy_data_yhi,
                              std::unique_ptr<ReadBndryPlanes>& m_r2d,
                              const Real& c_p,
                              const Real& rdOcp,
                              const MoistureComponentIndices& moisture_indices,
                              const bool use_wrf_bdy_density,
                              const int bdy_moist_nudge_type)

{
    BL_PROFILE_REGION("erf_add_moist_nudging_terms()");

    const Box domain = geom.Domain();
    IntVect ng  = S_data.nGrowVect();

    AMREX_ALWAYS_ASSERT(bdy_moist_nudge_type >= 0 && bdy_moist_nudge_type <= 3);
    AMREX_ALWAYS_ASSERT(RhoQ1_comp + n_qstate <= S_data.nComp());
    AMREX_ALWAYS_ASSERT(bdy_moist_nudge_type != 3 || !m_r2d);
    const bool separate_hydrometeors =
        wrf_bdy_has_separate_hydrometeors(moisture_indices);
    const GpuArray<int,6> moisture_comps = {
        moisture_indices.qv, moisture_indices.qc, moisture_indices.qi,
        moisture_indices.qr, moisture_indices.qs, moisture_indices.qg};
    const int n_moisture_targets = separate_hydrometeors
        ? 6 : (moisture_indices.qi >= 0 ? 3 : 2);

    // Temporary MF so we can nudge qv + qc to the bdy data
    //
    // NOTE: realbdy_compute_relaxation reads Rho_comp, RhoTheta_comp (for the Exner
    //       function) and RhoQ1_comp..RhoQ3_comp out of this buffer, so every one of
    //       those components must be filled here.  Rho and RhoTheta are contiguous,
    //       so they come over in a single copy.  The setVal is belt-and-braces so
    //       that any component we do not use is deterministic rather than recycled
    //       arena memory.
    MultiFab S_tmp(S_data.boxArray(), S_data.DistributionMap(), S_data.nComp(), ng);
    S_tmp.setVal(zero);
    MultiFab::Copy(S_tmp, S_data, Rho_comp, Rho_comp, 2, ng);

    if (bdy_moist_nudge_type == 0) {

        // This is our original approach in which we only nudge qv
        // based on the difference between qv and the values of qv in the bdy files
        MultiFab::Copy(S_tmp, S_data, RhoQ1_comp, RhoQ1_comp, 1, ng);

    } else if (bdy_moist_nudge_type == 1) {

        // This is a new approach in which we only nudge qv but
        // the amount by which it is nudged is  based on the difference
        // between (qv+qc+qi) and the values of qv in the bdy files
        // There is no associated heating/cooling.

        MultiFab::Copy(S_tmp, S_data, RhoQ1_comp, RhoQ1_comp, 1, ng);
        MultiFab::Add (S_tmp, S_data, RhoQ2_comp, RhoQ1_comp, 1, ng);
        if (n_qstate > 3) { // n_qstate > 3 guarantees that RhoQ3 is ice
            MultiFab::Add (S_tmp, S_data, RhoQ3_comp, RhoQ1_comp, 1, ng);
        }

    } else if (bdy_moist_nudge_type == 2) {

       // This is a new approach in which we nudge qv
       // based on the difference between qv and the values of qv in the bdy files
       // and also nudge qc (and qi) to 0.   The removal of qc (and qi) is treated
       // as if it is converted into qv to generate a source term for (rho theta),
       // but qv itself is not increased by the change in qc (and qi).
       //
       // NOTE: the relaxation reads qv and qc, plus qi when there is an ice
       //       species; it does not touch the precipitating species.
        int ncomp_q = (n_qstate > 3) ? 3 : 2;
        MultiFab::Copy(S_tmp, S_data, RhoQ1_comp, RhoQ1_comp, ncomp_q, ng);

    } else if (bdy_moist_nudge_type == 3) {

        for (int n = 0; n < n_moisture_targets; ++n) {
            const int comp = moisture_comps[n];
            if (comp >= 0) {
                AMREX_ALWAYS_ASSERT(comp < S_data.nComp());
                MultiFab::Copy(S_tmp, S_data, comp, comp, 1, ng);
            }
        }

    }

    int bdy_comp = BCVars::RhoQ1_bc_comp;
    Array4<Real> bdatxlo, bdatxhi, bdatylo, bdatyhi;

    // Relaxation constants
    Real F1 = one/(nudge_factor*static_cast<Real>(dt));

    // Domain bounds
    const auto& dom_hi = ubound(domain);
    const auto& dom_lo = lbound(domain);
    auto dx = geom.CellSizeArray();
    auto ProbHi = geom.ProbHiArray();
    auto ProbLo = geom.ProbLoArray();

    // Time interpolation
    double dT_d = bdy_time_interval;
    Real   dT   = static_cast<Real>(dT_d);

    //
    // Note that time (= start_time+old_stage_time)  is measured as total time
    //           start_bdy_time and final_bdy_time are also measured as total time
    //

    int n_time    = static_cast<int>( (time-start_bdy_time) /  dT);
    int n_time_p1 = n_time + 1;
    Real alpha    = static_cast<Real>(((time-start_bdy_time) - n_time * dT) / dT);

    // Do not over run the last bdy file
    if (time >= final_bdy_time) {
      n_time    = static_cast<int>( (final_bdy_time - start_bdy_time)/ dT);
      n_time_p1 = n_time;
      alpha     = zero;
    }

    AMREX_ALWAYS_ASSERT( alpha >= zero && alpha <= one);
    Real oma   = one - alpha;

    // Limiting offset
    int offset = width - 1;

    for (MFIter mfi(S_tmp,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box tbx  = mfi.tilebox();

        const Array4<const Real>& cons_arr = S_tmp.const_array(mfi);

        // We will add to source terms for moist variables and
        //    to source term for (rho theta)
        const Array4<      Real>& src_arr        = source.array(mfi);

        //
        // Note that old_stage_time_total = start_time+old_stage_time is total time
        //           start_bdy_time and final_bdy_time are total time
        //
        // moist_set_rhs(geom, tbx, new_cons_const, src_arr,
        //               old_stage_time_total, dt,
        //               start_bdy_time, final_bdy_time, bdy_time_interval,
        //               bdy_factor, width, domain,
        //               bdy_data_xlo, bdy_data_xhi,
        //               bdy_data_ylo, bdy_data_yhi,
        //               m_r2d);

        // Get bndry data
        if (m_r2d) {
            Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);
            bdatxlo = (*bndry_data[0])[0].array();
            bdatylo = (*bndry_data[1])[0].array();
            bdatxhi = (*bndry_data[3])[0].array();
            bdatyhi = (*bndry_data[4])[0].array();
        }

        // NOTE: The sizing of the temporary BDY FABS is
        //       GLOBAL and occurs over the entire BDY region.

        // Size the FABs
        //==========================================================
        // NOTE: No ghost cells, we force mask to be idx type (0,0,0)
        IntVect ng_vect(0);
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        realbdy_interior_bxs_xy(gdom, domain, width,
                                bx_xlo, bx_xhi,
                                bx_ylo, bx_yhi,
                                ng_vect, true);

        // Temporary FABs for storage (owned/filled on all ranks)
        const int ntargets = (bdy_moist_nudge_type == 3) ? n_moisture_targets : 1;
        FArrayBox QV_xlo, QV_xhi, QV_ylo, QV_yhi;
        QV_xlo.resize(bx_xlo,ntargets,The_Async_Arena()); QV_xhi.resize(bx_xhi,ntargets,The_Async_Arena());
        QV_ylo.resize(bx_ylo,ntargets,The_Async_Arena()); QV_yhi.resize(bx_yhi,ntargets,The_Async_Arena());

        // Populate FABs from bdy interpolation (primitive vars)
        //==========================================================
        const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][WRFBdyVars::QV].const_array();
        const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][WRFBdyVars::QV].const_array();
        const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][WRFBdyVars::QV].const_array();
        const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][WRFBdyVars::QV].const_array();
        const auto& bdatylo_n   = bdy_data_ylo[n_time   ][WRFBdyVars::QV].const_array();
        const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][WRFBdyVars::QV].const_array();
        const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][WRFBdyVars::QV].const_array();
        const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][WRFBdyVars::QV].const_array();
        Array4<const Real> qcdatxlo_n, qcdatxlo_np1, qcdatxhi_n, qcdatxhi_np1;
        Array4<const Real> qcdatylo_n, qcdatylo_np1, qcdatyhi_n, qcdatyhi_np1;
        Array4<const Real> qidatxlo_n, qidatxlo_np1, qidatxhi_n, qidatxhi_np1;
        Array4<const Real> qidatylo_n, qidatylo_np1, qidatyhi_n, qidatyhi_np1;
        Array4<const Real> qrdatxlo_n, qrdatxlo_np1, qrdatxhi_n, qrdatxhi_np1;
        Array4<const Real> qrdatylo_n, qrdatylo_np1, qrdatyhi_n, qrdatyhi_np1;
        Array4<const Real> qsdatxlo_n, qsdatxlo_np1, qsdatxhi_n, qsdatxhi_np1;
        Array4<const Real> qsdatylo_n, qsdatylo_np1, qsdatyhi_n, qsdatyhi_np1;
        Array4<const Real> qgdatxlo_n, qgdatxlo_np1, qgdatxhi_n, qgdatxhi_np1;
        Array4<const Real> qgdatylo_n, qgdatylo_np1, qgdatyhi_n, qgdatyhi_np1;
        if (bdy_moist_nudge_type == 3) {
            qcdatxlo_n   = bdy_data_xlo[n_time   ][RealBdyVars::QC].const_array();
            qcdatxlo_np1 = bdy_data_xlo[n_time_p1][RealBdyVars::QC].const_array();
            qcdatxhi_n   = bdy_data_xhi[n_time   ][RealBdyVars::QC].const_array();
            qcdatxhi_np1 = bdy_data_xhi[n_time_p1][RealBdyVars::QC].const_array();
            qcdatylo_n   = bdy_data_ylo[n_time   ][RealBdyVars::QC].const_array();
            qcdatylo_np1 = bdy_data_ylo[n_time_p1][RealBdyVars::QC].const_array();
            qcdatyhi_n   = bdy_data_yhi[n_time   ][RealBdyVars::QC].const_array();
            qcdatyhi_np1 = bdy_data_yhi[n_time_p1][RealBdyVars::QC].const_array();
            qidatxlo_n   = bdy_data_xlo[n_time   ][RealBdyVars::QI].const_array();
            qidatxlo_np1 = bdy_data_xlo[n_time_p1][RealBdyVars::QI].const_array();
            qidatxhi_n   = bdy_data_xhi[n_time   ][RealBdyVars::QI].const_array();
            qidatxhi_np1 = bdy_data_xhi[n_time_p1][RealBdyVars::QI].const_array();
            qidatylo_n   = bdy_data_ylo[n_time   ][RealBdyVars::QI].const_array();
            qidatylo_np1 = bdy_data_ylo[n_time_p1][RealBdyVars::QI].const_array();
            qidatyhi_n   = bdy_data_yhi[n_time   ][RealBdyVars::QI].const_array();
            qidatyhi_np1 = bdy_data_yhi[n_time_p1][RealBdyVars::QI].const_array();
            if (separate_hydrometeors) {
                qrdatxlo_n   = bdy_data_xlo[n_time   ][RealBdyHydrometeorVars::QR].const_array();
                qrdatxlo_np1 = bdy_data_xlo[n_time_p1][RealBdyHydrometeorVars::QR].const_array();
                qrdatxhi_n   = bdy_data_xhi[n_time   ][RealBdyHydrometeorVars::QR].const_array();
                qrdatxhi_np1 = bdy_data_xhi[n_time_p1][RealBdyHydrometeorVars::QR].const_array();
                qrdatylo_n   = bdy_data_ylo[n_time   ][RealBdyHydrometeorVars::QR].const_array();
                qrdatylo_np1 = bdy_data_ylo[n_time_p1][RealBdyHydrometeorVars::QR].const_array();
                qrdatyhi_n   = bdy_data_yhi[n_time   ][RealBdyHydrometeorVars::QR].const_array();
                qrdatyhi_np1 = bdy_data_yhi[n_time_p1][RealBdyHydrometeorVars::QR].const_array();
                qsdatxlo_n   = bdy_data_xlo[n_time   ][RealBdyHydrometeorVars::QS].const_array();
                qsdatxlo_np1 = bdy_data_xlo[n_time_p1][RealBdyHydrometeorVars::QS].const_array();
                qsdatxhi_n   = bdy_data_xhi[n_time   ][RealBdyHydrometeorVars::QS].const_array();
                qsdatxhi_np1 = bdy_data_xhi[n_time_p1][RealBdyHydrometeorVars::QS].const_array();
                qsdatylo_n   = bdy_data_ylo[n_time   ][RealBdyHydrometeorVars::QS].const_array();
                qsdatylo_np1 = bdy_data_ylo[n_time_p1][RealBdyHydrometeorVars::QS].const_array();
                qsdatyhi_n   = bdy_data_yhi[n_time   ][RealBdyHydrometeorVars::QS].const_array();
                qsdatyhi_np1 = bdy_data_yhi[n_time_p1][RealBdyHydrometeorVars::QS].const_array();
                qgdatxlo_n   = bdy_data_xlo[n_time   ][RealBdyHydrometeorVars::QG].const_array();
                qgdatxlo_np1 = bdy_data_xlo[n_time_p1][RealBdyHydrometeorVars::QG].const_array();
                qgdatxhi_n   = bdy_data_xhi[n_time   ][RealBdyHydrometeorVars::QG].const_array();
                qgdatxhi_np1 = bdy_data_xhi[n_time_p1][RealBdyHydrometeorVars::QG].const_array();
                qgdatylo_n   = bdy_data_ylo[n_time   ][RealBdyHydrometeorVars::QG].const_array();
                qgdatylo_np1 = bdy_data_ylo[n_time_p1][RealBdyHydrometeorVars::QG].const_array();
                qgdatyhi_n   = bdy_data_yhi[n_time   ][RealBdyHydrometeorVars::QG].const_array();
                qgdatyhi_np1 = bdy_data_yhi[n_time_p1][RealBdyHydrometeorVars::QG].const_array();
            }
        }
        Array4<const Real> rdatxlo_n, rdatxlo_np1, rdatxhi_n, rdatxhi_np1;
        Array4<const Real> rdatylo_n, rdatylo_np1, rdatyhi_n, rdatyhi_np1;
        if (use_wrf_bdy_density) {
            rdatxlo_n   = bdy_data_xlo[n_time   ][WRFBdyVars::R].const_array();
            rdatxlo_np1 = bdy_data_xlo[n_time_p1][WRFBdyVars::R].const_array();
            rdatxhi_n   = bdy_data_xhi[n_time   ][WRFBdyVars::R].const_array();
            rdatxhi_np1 = bdy_data_xhi[n_time_p1][WRFBdyVars::R].const_array();
            rdatylo_n   = bdy_data_ylo[n_time   ][WRFBdyVars::R].const_array();
            rdatylo_np1 = bdy_data_ylo[n_time_p1][WRFBdyVars::R].const_array();
            rdatyhi_n   = bdy_data_yhi[n_time   ][WRFBdyVars::R].const_array();
            rdatyhi_np1 = bdy_data_yhi[n_time_p1][WRFBdyVars::R].const_array();
        }

        // Get Array4 of interpolated values
        Array4<Real> arr_xlo = QV_xlo.array();  Array4<Real> arr_xhi = QV_xhi.array();
        Array4<Real> arr_ylo = QV_ylo.array();  Array4<Real> arr_yhi = QV_yhi.array();

        Box gtbx = grow(tbx,ng_vect);
        Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
        realbdy_interior_bxs_xy(gtbx, domain, width,
                                tbx_xlo, tbx_xhi,
                                tbx_ylo, tbx_yhi,
                                ng_vect, true);

        // Populate with interpolation (protect from ghost cells)
        const bool has_cloud_ice = moisture_indices.qi >= 0;
        ParallelFor(tbx_xlo, tbx_xhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::min(std::max(i , dom_lo.x), dom_lo.x+offset);
            int jj = std::min(std::max(j , dom_lo.y), dom_hi.y       );
            Real rho = cons_arr(i,j,k,Rho_comp);
            if (use_wrf_bdy_density) {
                rho = oma*rdatxlo_n(ii,jj,k) + alpha*rdatxlo_np1(ii,jj,k);
            }
            if (bdy_moist_nudge_type == 3) {
                arr_xlo(i,j,k,0) = wrf_moisture_target(rho, bdatxlo_n(ii,jj,k),
                                                       bdatxlo_np1(ii,jj,k), alpha);
                arr_xlo(i,j,k,1) = wrf_moisture_target(rho, qcdatxlo_n(ii,jj,k),
                                                       qcdatxlo_np1(ii,jj,k), alpha);
                if (has_cloud_ice) {
                    arr_xlo(i,j,k,2) = wrf_moisture_target(rho, qidatxlo_n(ii,jj,k),
                                                          qidatxlo_np1(ii,jj,k), alpha);
                }
                if (separate_hydrometeors) {
                    arr_xlo(i,j,k,3) = wrf_moisture_target(rho, qrdatxlo_n(ii,jj,k), qrdatxlo_np1(ii,jj,k), alpha);
                    arr_xlo(i,j,k,4) = wrf_moisture_target(rho, qsdatxlo_n(ii,jj,k), qsdatxlo_np1(ii,jj,k), alpha);
                    arr_xlo(i,j,k,5) = wrf_moisture_target(rho, qgdatxlo_n(ii,jj,k), qgdatxlo_np1(ii,jj,k), alpha);
                }
            } else {
                arr_xlo(i,j,k) = (bdatxlo) ? rho * bdatxlo(ii,jj,k,bdy_comp) :
                    rho * ( oma   * bdatxlo_n  (ii,jj,k)
                          + alpha * bdatxlo_np1(ii,jj,k) );
            }
        }    ,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::min(std::max(i , dom_hi.x-offset), dom_hi.x);
            int jj = std::min(std::max(j , dom_lo.y       ), dom_hi.y);
            Real rho = cons_arr(i,j,k,Rho_comp);
            if (use_wrf_bdy_density) {
                rho = oma*rdatxhi_n(ii,jj,k) + alpha*rdatxhi_np1(ii,jj,k);
            }
            if (bdy_moist_nudge_type == 3) {
                arr_xhi(i,j,k,0) = wrf_moisture_target(rho, bdatxhi_n(ii,jj,k),
                                                       bdatxhi_np1(ii,jj,k), alpha);
                arr_xhi(i,j,k,1) = wrf_moisture_target(rho, qcdatxhi_n(ii,jj,k),
                                                       qcdatxhi_np1(ii,jj,k), alpha);
                if (has_cloud_ice) {
                    arr_xhi(i,j,k,2) = wrf_moisture_target(rho, qidatxhi_n(ii,jj,k),
                                                          qidatxhi_np1(ii,jj,k), alpha);
                }
                if (separate_hydrometeors) {
                    arr_xhi(i,j,k,3) = wrf_moisture_target(rho, qrdatxhi_n(ii,jj,k), qrdatxhi_np1(ii,jj,k), alpha);
                    arr_xhi(i,j,k,4) = wrf_moisture_target(rho, qsdatxhi_n(ii,jj,k), qsdatxhi_np1(ii,jj,k), alpha);
                    arr_xhi(i,j,k,5) = wrf_moisture_target(rho, qgdatxhi_n(ii,jj,k), qgdatxhi_np1(ii,jj,k), alpha);
                }
            } else {
                arr_xhi(i,j,k) = (bdatxhi) ? rho * bdatxhi(ii,jj,k,bdy_comp) :
                    rho * ( oma   * bdatxhi_n  (ii,jj,k)
                          + alpha * bdatxhi_np1(ii,jj,k) );
            }
        });

        ParallelFor(tbx_ylo, tbx_yhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::min(std::max(i , dom_lo.x), dom_hi.x       );
            int jj = std::min(std::max(j , dom_lo.y), dom_lo.y+offset);
            Real rho = cons_arr(i,j,k,Rho_comp);
            if (use_wrf_bdy_density) {
                rho = oma*rdatylo_n(ii,jj,k) + alpha*rdatylo_np1(ii,jj,k);
            }
            if (bdy_moist_nudge_type == 3) {
                arr_ylo(i,j,k,0) = wrf_moisture_target(rho, bdatylo_n(ii,jj,k),
                                                       bdatylo_np1(ii,jj,k), alpha);
                arr_ylo(i,j,k,1) = wrf_moisture_target(rho, qcdatylo_n(ii,jj,k),
                                                       qcdatylo_np1(ii,jj,k), alpha);
                if (has_cloud_ice) {
                    arr_ylo(i,j,k,2) = wrf_moisture_target(rho, qidatylo_n(ii,jj,k),
                                                          qidatylo_np1(ii,jj,k), alpha);
                }
                if (separate_hydrometeors) {
                    arr_ylo(i,j,k,3) = wrf_moisture_target(rho, qrdatylo_n(ii,jj,k), qrdatylo_np1(ii,jj,k), alpha);
                    arr_ylo(i,j,k,4) = wrf_moisture_target(rho, qsdatylo_n(ii,jj,k), qsdatylo_np1(ii,jj,k), alpha);
                    arr_ylo(i,j,k,5) = wrf_moisture_target(rho, qgdatylo_n(ii,jj,k), qgdatylo_np1(ii,jj,k), alpha);
                }
            } else {
                arr_ylo(i,j,k) = (bdatylo) ? rho * bdatylo(ii,jj,k,bdy_comp) :
                    rho * ( oma   * bdatylo_n  (ii,jj,k)
                          + alpha * bdatylo_np1(ii,jj,k) );
            }
        },
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           int ii = std::min(std::max(i , dom_lo.x       ), dom_hi.x);
           int jj = std::min(std::max(j , dom_hi.y-offset), dom_hi.y);
               jj = std::min(jj, dom_hi.y);
           Real rho = cons_arr(i,j,k,Rho_comp);
           if (use_wrf_bdy_density) {
               rho = oma*rdatyhi_n(ii,jj,k) + alpha*rdatyhi_np1(ii,jj,k);
           }
           if (bdy_moist_nudge_type == 3) {
               arr_yhi(i,j,k,0) = wrf_moisture_target(rho, bdatyhi_n(ii,jj,k),
                                                      bdatyhi_np1(ii,jj,k), alpha);
               arr_yhi(i,j,k,1) = wrf_moisture_target(rho, qcdatyhi_n(ii,jj,k),
                                                      qcdatyhi_np1(ii,jj,k), alpha);
               if (has_cloud_ice) {
                   arr_yhi(i,j,k,2) = wrf_moisture_target(rho, qidatyhi_n(ii,jj,k),
                                                         qidatyhi_np1(ii,jj,k), alpha);
               }
               if (separate_hydrometeors) {
                   arr_yhi(i,j,k,3) = wrf_moisture_target(rho, qrdatyhi_n(ii,jj,k), qrdatyhi_np1(ii,jj,k), alpha);
                   arr_yhi(i,j,k,4) = wrf_moisture_target(rho, qsdatyhi_n(ii,jj,k), qsdatyhi_np1(ii,jj,k), alpha);
                   arr_yhi(i,j,k,5) = wrf_moisture_target(rho, qgdatyhi_n(ii,jj,k), qgdatyhi_np1(ii,jj,k), alpha);
               }
           } else {
               arr_yhi(i,j,k) = (bdatyhi) ? rho * bdatyhi(ii,jj,k,bdy_comp) :
                   rho * ( oma   * bdatyhi_n  (ii,jj,k)
                         + alpha * bdatyhi_np1(ii,jj,k) );
           }
       });

       realbdy_interior_bxs_xy(tbx, domain, width,
                               tbx_xlo, tbx_xhi,
                               tbx_ylo, tbx_yhi,
                               ng_vect);

       //
       // Add relaxation terms for moist variables and (rho theta) to existing source terms
       //
       if (bdy_moist_nudge_type == 3) {
           realbdy_compute_moisture_relaxation(
               moisture_comps, n_moisture_targets, width, dx, ProbLo, ProbHi, F1,
               tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi,
               arr_xlo, arr_xhi, arr_ylo, arr_yhi, cons_arr, src_arr);
       } else {
           realbdy_compute_relaxation(RhoQ1_comp, n_qstate,
                                      width, dx, ProbLo, ProbHi, F1,
                                      tbx_xlo , tbx_xhi , tbx_ylo , tbx_yhi ,
                                      arr_xlo , arr_xhi , arr_ylo , arr_yhi ,
                                      cons_arr, src_arr, c_p, rdOcp,
                                      bdy_moist_nudge_type);
       }
   } // mfi
}
#endif
