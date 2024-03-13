#include <iomanip>

#include "ERF.H"
#include "EOS.H"

using namespace amrex;

/**
 * Writes 1-dimensional averaged quantities as profiles to output log files
 * at the given time. Quantities are output at their native grid locations,
 * therefore W and associated quantities <(*)'w'>, tau13, and tau23 (where *
 * includes u, v, p, theta, and k) will be output at staggered heights (i.e., z
 * faces) rather than cell-center heights to avoid performing additional
 * averaging. The unstaggered quantities are associated with the left cell face,
 * i.e., they share the same k index. These quantities will have a zero value at
 * the max height of Nz+1.
 *
 * @param time Current time
 */
void
ERF::write_1D_profiles_stag(Real time)
{
    BL_PROFILE("ERF::write_1D_profiles()");

    if (verbose <= 0)
      return;

    int datwidth = 14;
    int datprecision = 6;
    int timeprecision = 13; // e.g., 1-yr LES: 31,536,000 s with dt ~ 0.01 ==> min prec = 10

    if (verbose > 0 && NumDataLogs() > 1)
    {
        // Define the 1d arrays we will need
        Gpu::HostVector<Real> h_avg_u, h_avg_v, h_avg_w;
        Gpu::HostVector<Real> h_avg_rho, h_avg_th, h_avg_ksgs;
        Gpu::HostVector<Real> h_avg_uth, h_avg_vth, h_avg_wth, h_avg_thth;
        Gpu::HostVector<Real> h_avg_uu, h_avg_uv, h_avg_uw, h_avg_vv, h_avg_vw, h_avg_ww;
        Gpu::HostVector<Real> h_avg_uiuiu, h_avg_uiuiv, h_avg_uiuiw;
        Gpu::HostVector<Real> h_avg_p, h_avg_pu, h_avg_pv, h_avg_pw;
        Gpu::HostVector<Real> h_avg_tau11, h_avg_tau12, h_avg_tau13, h_avg_tau22, h_avg_tau23, h_avg_tau33;
        Gpu::HostVector<Real> h_avg_sgshfx, h_avg_sgsdiss; // only output tau_{theta,w} and epsilon for now

        if (NumDataLogs() > 1) {
            derive_diag_profiles_stag(h_avg_u, h_avg_v, h_avg_w,
                                      h_avg_rho, h_avg_th, h_avg_ksgs,
                                      h_avg_uu, h_avg_uv, h_avg_uw, h_avg_vv, h_avg_vw, h_avg_ww,
                                      h_avg_uth, h_avg_vth, h_avg_wth, h_avg_thth,
                                      h_avg_uiuiu, h_avg_uiuiv, h_avg_uiuiw,
                                      h_avg_p, h_avg_pu, h_avg_pv, h_avg_pw);
        }

        if (NumDataLogs() > 3) {
            derive_stress_profiles_stag(h_avg_tau11, h_avg_tau12, h_avg_tau13,
                                        h_avg_tau22, h_avg_tau23, h_avg_tau33,
                                        h_avg_sgshfx,
                                        h_avg_sgsdiss);
        }

        int unstag_size =  h_avg_w.size() - 1; // _un_staggered heights

        auto const& dx = geom[0].CellSizeArray();
        if (ParallelDescriptor::IOProcessor()) {
            if (NumDataLogs() > 1) {
                std::ostream& data_log1 = DataLog(1);
                if (data_log1.good()) {
                  // Write the quantities at this time
                  for (int k = 0; k < unstag_size; k++) {
                      Real z = k * dx[2];
                      data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_u[k]  << " " << h_avg_v[k] << " " << h_avg_w[k] << " "
                                << h_avg_rho[k] << " " << h_avg_th[k] << " " << h_avg_ksgs[k]
                                << std::endl;
                  } // loop over z
                  // Write top face values
                  data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << unstag_size * dx[2] << " "
                            << 0  << " " << 0 << " " << h_avg_w[unstag_size+1] << " "
                            << 0 << " " << 0 << " " << 0
                            << std::endl;
                } // if good
            } // NumDataLogs

            if (NumDataLogs() > 2) {
                std::ostream& data_log2 = DataLog(2);
                if (data_log2.good()) {
                  // Write the perturbational quantities at this time
                  // For surface values (k=0), assume w = uw = vw = ww = 0
                  Real w_cc  = h_avg_w[1] / 2;  // w at first cell center
                  Real uw_cc = h_avg_uw[1] / 2; // u*w at first cell center
                  Real vw_cc = h_avg_vw[1] / 2; // v*w at first cell center
                  Real ww_cc = h_avg_ww[1] / 2; // w*w at first cell center
                  data_log2 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << 0 << " "
                            << h_avg_uu[0]   - h_avg_u[0]*h_avg_u[0]   << " " // u'u'
                            << h_avg_uv[0]   - h_avg_u[0]*h_avg_v[0]   << " " // u'v'
                            << 0                                       << " " // u'w'
                            << h_avg_vv[0]   - h_avg_v[0]*h_avg_v[0]   << " " // v'v'
                            << 0                                       << " " // v'w'
                            << 0                                       << " " // w'w'
                            << h_avg_uth[0]  - h_avg_u[0]*h_avg_th[0]  << " " // u'th'
                            << h_avg_vth[0]  - h_avg_v[0]*h_avg_th[0]  << " " // v'th'
                            << 0                                       << " " // w'th'
                            << h_avg_thth[0] - h_avg_th[0]*h_avg_th[0] << " " // th'th'
                            << h_avg_uiuiu[0]
                             - (h_avg_uu[0] + h_avg_vv[0] + ww_cc)*h_avg_u[0]
                             - 2*(h_avg_u[0]*h_avg_uu[0] + h_avg_v[0]*h_avg_uv[0] + w_cc*uw_cc)
                             + 2*(h_avg_u[0]*h_avg_u[0] + h_avg_v[0]*h_avg_v[0] + w_cc*w_cc)*h_avg_u[0]
                               << " " // (u'_i u'_i)u'
                            << h_avg_uiuiv[0]
                             - (h_avg_uu[0] + h_avg_vv[0] + ww_cc)*h_avg_v[0]
                             - 2*(h_avg_u[0]*h_avg_uv[0] + h_avg_v[0]*h_avg_vv[0] + w_cc*vw_cc)
                             + 2*(h_avg_u[0]*h_avg_u[0] + h_avg_v[0]*h_avg_v[0] + w_cc*w_cc)*h_avg_v[0]
                               << " " // (u'_i u'_i)v'
                            << 0 << " " // (u'_i u'_i)w'
                            << h_avg_pu[0]   - h_avg_p[0]*h_avg_u[0]   << " " // pu'
                            << h_avg_pv[0]   - h_avg_p[0]*h_avg_v[0]   << " " // pv'
                            << 0                                              // pw'
                            << std::endl;

                  // For internal values, interpolate scalar quantities to faces
                  for (int k = 1; k < unstag_size; k++) {
                      Real z = k * dx[2];
                      Real uface  = 0.5*(h_avg_u[k]  + h_avg_u[k-1]);
                      Real vface  = 0.5*(h_avg_v[k]  + h_avg_v[k-1]);
                      Real thface = 0.5*(h_avg_th[k] + h_avg_th[k-1]);
                      Real pface  = 0.5*(h_avg_p[k]  + h_avg_p[k-1]);
                      Real uuface = 0.5*(h_avg_uu[k] + h_avg_uu[k-1]);
                      Real uvface = 0.5*(h_avg_uv[k] + h_avg_uv[k-1]);
                      Real vvface = 0.5*(h_avg_vv[k] + h_avg_vv[k-1]);
                      w_cc   = 0.5*(h_avg_w[k-1]  + h_avg_w[k]);
                      uw_cc  = 0.5*(h_avg_uw[k-1] + h_avg_uw[k]);
                      vw_cc  = 0.5*(h_avg_vw[k-1] + h_avg_vw[k]);
                      ww_cc  = 0.5*(h_avg_ww[k-1] + h_avg_ww[k]);
                      amrex::ignore_unused(uvface);
                      data_log2 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_uu[k]   - h_avg_u[k]*h_avg_u[k]   << " " // u'u'
                                << h_avg_uv[k]   - h_avg_u[k]*h_avg_v[k]   << " " // u'v'
                                << h_avg_uw[k]   -      uface*h_avg_w[k]   << " " // u'w'
                                << h_avg_vv[k]   - h_avg_v[k]*h_avg_v[k]   << " " // v'v'
                                << h_avg_vw[k]   -      vface*h_avg_w[k]   << " " // v'w'
                                << h_avg_ww[k]   - h_avg_w[k]*h_avg_w[k]   << " " // w'w'
                                << h_avg_uth[k]  - h_avg_u[k]*h_avg_th[k]  << " " // u'th'
                                << h_avg_vth[k]  - h_avg_v[k]*h_avg_th[k]  << " " // v'th'
                                << h_avg_wth[k]  - h_avg_w[k]*thface       << " " // w'th'
                                << h_avg_thth[k] - h_avg_th[k]*h_avg_th[k] << " " // th'th'
                                // Note: <u'_i u'_i u'_j> =   <u_i u_i u_j>
                                //                        -   <u_i u_i> * <u_j>
                                //                        - 2*<u_i> * <u_i u_j>
                                //                        + 2*<u_i>*<u_i> * <u_j>
                                << h_avg_uiuiu[k]
                                 - (h_avg_uu[k] + h_avg_vv[k] + ww_cc)*h_avg_u[k]
                                 - 2*(h_avg_u[k]*h_avg_uu[k] + h_avg_v[k]*h_avg_uv[k] + w_cc*uw_cc)
                                 + 2*(h_avg_u[k]*h_avg_u[k] + h_avg_v[k]*h_avg_v[k] + w_cc*w_cc)*h_avg_u[k]
                                   << " " // cell-centered (u'_i u'_i)u'
                                << h_avg_uiuiv[k]
                                 - (h_avg_uu[k] + h_avg_vv[k] + ww_cc)*h_avg_v[k]
                                 - 2*(h_avg_u[k]*h_avg_uv[k] + h_avg_v[k]*h_avg_vv[k] + w_cc*vw_cc)
                                 + 2*(h_avg_u[k]*h_avg_u[k] + h_avg_v[k]*h_avg_v[k] + w_cc*w_cc)*h_avg_v[k]
                                   << " " // cell-centered (u'_i u'_i)v'
                                << h_avg_uiuiw[k]
                                 - (uuface + vvface + h_avg_ww[k])*h_avg_w[k]
                                 - 2*(uface*h_avg_uw[k] + vface*h_avg_vw[k] + h_avg_w[k]*h_avg_ww[k])
                                 + 2*(uface*uface + vface*vface + h_avg_w[k]*h_avg_w[k])*h_avg_w[k]
                                   << " " // face-centered (u'_i u'_i)w'
                                << h_avg_pu[k]   - h_avg_p[k]*h_avg_u[k]   << " " // pu'
                                << h_avg_pv[k]   - h_avg_p[k]*h_avg_v[k]   << " " // pv'
                                << h_avg_pw[k]   -      pface*h_avg_w[k]          // pw'
                                << std::endl;
                  } // loop over z

                  // Write top face values, extrapolating scalar quantities
                  const int k = unstag_size;
                  Real uface  = 1.5*h_avg_u[k-1]  - 0.5*h_avg_u[k-2];
                  Real vface  = 1.5*h_avg_v[k-1]  - 0.5*h_avg_v[k-2];
                  Real thface = 1.5*h_avg_th[k-1] - 0.5*h_avg_th[k-2];
                  Real pface  = 1.5*h_avg_p[k-1]  - 0.5*h_avg_p[k-2];
                  Real uuface = 1.5*h_avg_uu[k-1] - 0.5*h_avg_uu[k-2];
                  Real uvface = 1.5*h_avg_uv[k-1] - 0.5*h_avg_uv[k-2];
                  Real vvface = 1.5*h_avg_vv[k-1] - 0.5*h_avg_vv[k-2];
                  amrex::ignore_unused(uvface);
                  data_log2 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << k * dx[2] << " "
                            << 0                                     << " " // u'u'
                            << 0                                     << " " // u'v'
                            << h_avg_uw[k]   -      uface*h_avg_w[k] << " " // u'w'
                            << 0                                     << " " // v'v'
                            << h_avg_vw[k]   -      vface*h_avg_w[k] << " " // v'w'
                            << h_avg_ww[k]   - h_avg_w[k]*h_avg_w[k] << " " // w'w'
                            << 0                                     << " " // u'th'
                            << 0                                     << " " // v'th'
                            << h_avg_wth[k]  -     thface*h_avg_w[k] << " " // w'th'
                            << 0                                     << " " // th'th'
                            << 0                                     << " " // (u'_i u'_i)u'
                            << 0                                     << " " // (u'_i u'_i)v'
                            << h_avg_uiuiw[k]
                             - (uuface + vvface + h_avg_ww[k])*h_avg_w[k]
                             - 2*(uface*h_avg_uw[k] + vface*h_avg_vw[k] + h_avg_w[k]*h_avg_ww[k])
                             + 2*(uface*uface + vface*vface + h_avg_w[k]*h_avg_w[k])*h_avg_w[k]
                               << " " // (u'_i u'_i)w'
                            << 0                                     << " " // pu'
                            << 0                                     << " " // pv'
                            << h_avg_pw[k]   -      pface*h_avg_w[k]        // pw'
                            << std::endl;
                } // if good
            } // NumDataLogs

            if (NumDataLogs() > 3) {
                std::ostream& data_log3 = DataLog(3);
                if (data_log3.good()) {
                  // Write the average stresses
                  for (int k = 0; k < unstag_size; k++) {
                      Real z = k * dx[2];
                      data_log3 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_tau11[k] << " " << h_avg_tau12[k] << " " << h_avg_tau13[k] << " "
                                << h_avg_tau22[k] << " " << h_avg_tau23[k] << " " << h_avg_tau33[k] << " "
                                << h_avg_sgshfx[k] << " " << h_avg_sgsdiss[k]
                                << std::endl;
                  } // loop over z
                  // Write top face values
                  data_log3 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << unstag_size * dx[2] << " "
                            << 0 << " " << 0 << " " << h_avg_tau13[unstag_size] << " "
                            << 0 << " " << h_avg_tau23[unstag_size] << " " << 0 << " "
                            << 0 << " " << 0
                            << std::endl;
                } // if good
            } // NumDataLogs
        } // if IOProcessor
    } // if verbose
}

/**
 * Computes the profiles for diagnostic quantities at _staggered_ heights.
 *
 * @param h_avg_u Profile for x-velocity on Host
 * @param h_avg_v Profile for y-velocity on Host
 * @param h_avg_w Profile for z-velocity on Host
 * @param h_avg_rho Profile for density on Host
 * @param h_avg_th Profile for potential temperature on Host
 * @param h_avg_ksgs Profile for Kinetic Energy on Host
 * @param h_avg_uu Profile for x-velocity squared on Host
 * @param h_avg_uv Profile for x-velocity * y-velocity on Host
 * @param h_avg_uw Profile for x-velocity * z-velocity on Host
 * @param h_avg_vv Profile for y-velocity squared on Host
 * @param h_avg_vw Profile for y-velocity * z-velocity on Host
 * @param h_avg_ww Profile for z-velocity squared on Host
 * @param h_avg_uth Profile for x-velocity * potential temperature on Host
 * @param h_avg_uiuiu Profile for u_i*u_i*u triple product on Host
 * @param h_avg_uiuiv Profile for u_i*u_i*v triple product on Host
 * @param h_avg_uiuiw Profile for u_i*u_i*w triple product on Host
 * @param h_avg_p Profile for pressure perturbation on Host
 * @param h_avg_pu Profile for pressure perturbation * x-velocity on Host
 * @param h_avg_pv Profile for pressure perturbation * y-velocity on Host
 * @param h_avg_pw Profile for pressure perturbation * z-velocity on Host
 */
void
ERF::derive_diag_profiles_stag(Gpu::HostVector<Real>& h_avg_u   , Gpu::HostVector<Real>& h_avg_v  , Gpu::HostVector<Real>& h_avg_w,
                               Gpu::HostVector<Real>& h_avg_rho , Gpu::HostVector<Real>& h_avg_th , Gpu::HostVector<Real>& h_avg_ksgs,
                               Gpu::HostVector<Real>& h_avg_uu  , Gpu::HostVector<Real>& h_avg_uv , Gpu::HostVector<Real>& h_avg_uw,
                               Gpu::HostVector<Real>& h_avg_vv  , Gpu::HostVector<Real>& h_avg_vw , Gpu::HostVector<Real>& h_avg_ww,
                               Gpu::HostVector<Real>& h_avg_uth , Gpu::HostVector<Real>& h_avg_vth, Gpu::HostVector<Real>& h_avg_wth,
                               Gpu::HostVector<Real>& h_avg_thth,
                               Gpu::HostVector<Real>& h_avg_uiuiu  , Gpu::HostVector<Real>& h_avg_uiuiv , Gpu::HostVector<Real>& h_avg_uiuiw,
                               Gpu::HostVector<Real>& h_avg_p,
                               Gpu::HostVector<Real>& h_avg_pu  , Gpu::HostVector<Real>& h_avg_pv , Gpu::HostVector<Real>& h_avg_pw)
{
    // We assume that this is always called at level 0
    int lev = 0;

    bool l_use_KE  = (solverChoice.turbChoice[lev].les_type == LESType::Deardorff);
    bool l_use_QKE = solverChoice.turbChoice[lev].use_QKE && solverChoice.turbChoice[lev].advect_QKE;

    // This will hold rho, theta, ksgs, uu, uv, vv, uth, vth, thth, uiuiu, uiuiv, p, pu, pv
    //         index:   0      1     2   3   4   5    6    7     8      9     10  11 12  13
    MultiFab mf_out(grids[lev], dmap[lev], 14, 0);

    // This will hold uw, vw, ww, wth, uiuiw, pw (note: uiui == u_i*u_i = u*u + v*v + w*w)
    //         index:  0   1   2    3      4   5
    MultiFab mf_out_stag(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 6, 0);

    // This is only used to average u and v
    MultiFab mf_vels(grids[lev], dmap[lev], 2, 0);

    MultiFab  u_cc(mf_vels, make_alias, 0, 1); // u at cell centers
    MultiFab  v_cc(mf_vels, make_alias, 1, 1); // v at cell centers
    MultiFab  w_fc(vars_new[lev][Vars::zvel], make_alias, 0, 1); // w at face centers

    int zdir = 2;
    auto domain = geom[0].Domain();
    Box stag_domain = domain;
    stag_domain.convert(IntVect(0,0,1));

    int nvars = vars_new[lev][Vars::cons].nComp();
    MultiFab mf_cons(vars_new[lev][Vars::cons], make_alias, 0, nvars);

    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component

    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);

    for ( MFIter mfi(mf_cons,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& fab_arr = mf_out.array(mfi);
        const Array4<Real>& fab_arr_stag = mf_out_stag.array(mfi);
        const Array4<Real>& u_arr = vars_new[lev][Vars::xvel].array(mfi);
        const Array4<Real>& v_arr = vars_new[lev][Vars::yvel].array(mfi);
        const Array4<Real>& u_cc_arr =  u_cc.array(mfi);
        const Array4<Real>& v_cc_arr =  v_cc.array(mfi);
        const Array4<Real>& w_fc_arr =  w_fc.array(mfi);
        const Array4<Real>& cons_arr = mf_cons.array(mfi);
        const Array4<Real>&   p0_arr = p_hse.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_cc_arr(i,j,k) = 0.5 * (u_arr(i,j,k) + u_arr(i+1,j  ,k));
            v_cc_arr(i,j,k) = 0.5 * (v_arr(i,j,k) + v_arr(i  ,j+1,k));

            Real theta = cons_arr(i,j,k,RhoTheta_comp) / cons_arr(i,j,k,Rho_comp);
            fab_arr(i, j, k, 0) = cons_arr(i,j,k,Rho_comp);
            fab_arr(i, j, k, 1) = theta;
            Real ksgs = 0.0;
            if (l_use_KE)
                ksgs = cons_arr(i,j,k,RhoKE_comp) / cons_arr(i,j,k,Rho_comp);
            else if (l_use_QKE)
                ksgs = cons_arr(i,j,k,RhoQKE_comp) / cons_arr(i,j,k,Rho_comp);
            fab_arr(i, j, k, 2) = ksgs;
            fab_arr(i, j, k, 3) = u_cc_arr(i,j,k) * u_cc_arr(i,j,k);   // u*u
            fab_arr(i, j, k, 4) = u_cc_arr(i,j,k) * v_cc_arr(i,j,k);   // u*v
            fab_arr(i, j, k, 5) = v_cc_arr(i,j,k) * v_cc_arr(i,j,k);   // v*v
            fab_arr(i, j, k, 6) = u_cc_arr(i,j,k) * theta;             // u*th
            fab_arr(i, j, k, 7) = v_cc_arr(i,j,k) * theta;             // v*th
            fab_arr(i, j, k, 8) = theta * theta;                       // th*th

            Real wcc = 0.5 * (w_fc_arr(i,j,k) + w_fc_arr(i,j,k+1));
            Real uiui = fab_arr(i,j,k,3) + fab_arr(i,j,k,5) + wcc*wcc;
            fab_arr(i, j, k,9 ) = uiui * u_cc_arr(i,j,k);           // (ui*ui)*u
            fab_arr(i, j, k,10) = uiui * v_cc_arr(i,j,k);           // (ui*ui)*v

            if (!use_moisture) {
                Real p = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp));
                p -= p0_arr(i,j,k);
                fab_arr(i, j, k,11) = p;                       // p'
                fab_arr(i, j, k,12) = p * u_cc_arr(i,j,k);     // p'*u
                fab_arr(i, j, k,13) = p * v_cc_arr(i,j,k);     // p'*v
            }
        });

        const Box& zbx = mfi.tilebox(IntVect(0,0,1));
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // average to z faces (first to cell centers, then in z)
            Real uface = 0.25 * ( u_arr(i  ,j,k) + u_arr(i  ,j,k-1)
                                + u_arr(i+1,j,k) + u_arr(i+1,j,k-1));
            Real vface = 0.25 * ( v_arr(i,j  ,k) + v_arr(i,j  ,k-1)
                                + v_arr(i,j+1,k) + v_arr(i,j+1,k-1));
            Real theta0 = cons_arr(i,j,k  ,RhoTheta_comp) / cons_arr(i,j,k  ,Rho_comp);
            Real theta1 = cons_arr(i,j,k-1,RhoTheta_comp) / cons_arr(i,j,k-1,Rho_comp);
            Real thface = 0.5*(theta0 + theta1);
            fab_arr_stag(i,j,k,0) =           uface * w_fc_arr(i,j,k); // u*w
            fab_arr_stag(i,j,k,1) =           vface * w_fc_arr(i,j,k); // v*w
            fab_arr_stag(i,j,k,2) = w_fc_arr(i,j,k) * w_fc_arr(i,j,k); // w*w
            fab_arr_stag(i,j,k,3) =          thface * w_fc_arr(i,j,k); // th*w
            Real uiui = uface*uface + vface*vface + fab_arr_stag(i,j,k,2);
            fab_arr_stag(i,j,k,4) = uiui * w_fc_arr(i,j,k); // (ui*ui)*w
            if (!use_moisture) {
                Real p0 = getPgivenRTh(cons_arr(i, j, k  , RhoTheta_comp)) - p0_arr(i,j,k  );
                Real p1 = getPgivenRTh(cons_arr(i, j, k-1, RhoTheta_comp)) - p0_arr(i,j,k-1);
                Real pface = 0.5 * (p0 + p1);
                fab_arr_stag(i,j,k,5) = pface * w_fc_arr(i,j,k);       // p'*w
            }
        });

    } // mfi

    if (use_moisture)
    {
        for ( MFIter mfi(mf_cons,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Array4<Real>& fab_arr  = mf_out.array(mfi);
            const Array4<Real>& fab_arr_stag  = mf_out_stag.array(mfi);
            const Array4<Real>& cons_arr = mf_cons.array(mfi);
            const Array4<Real>& u_cc_arr =  u_cc.array(mfi);
            const Array4<Real>& v_cc_arr =  v_cc.array(mfi);
            const Array4<Real>& w_fc_arr =  w_fc.array(mfi);
            const Array4<Real>&   p0_arr = p_hse.array(mfi);
            const Array4<Real>&   qv_arr = qmoist[0][0]->array(mfi); // TODO: Is this written only on lev 0?

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real p = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp), qv_arr(i,j,k));
                p -= p0_arr(i,j,k);
                fab_arr(i, j, k,11) = p;                       // p'
                fab_arr(i, j, k,12) = p * u_cc_arr(i,j,k);     // p'*u
                fab_arr(i, j, k,13) = p * v_cc_arr(i,j,k);     // p'*v
            });

            const Box& zbx = mfi.tilebox(IntVect(0,0,1));
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real p0 = getPgivenRTh(cons_arr(i, j, k  , RhoTheta_comp), qv_arr(i,j,k  )) - p0_arr(i,j,k  );
                Real p1 = getPgivenRTh(cons_arr(i, j, k-1, RhoTheta_comp), qv_arr(i,j,k-1)) - p0_arr(i,j,k-1);
                Real pface = 0.5 * (p0 + p1);
                fab_arr_stag(i,j,k,5) = pface * w_fc_arr(i,j,k); // p'*w
            });
        } // mfi
    } // use_moisture

    // Sum in the horizontal plane
    h_avg_u  = sumToLine(u_cc,0,1,     domain,zdir);
    h_avg_v  = sumToLine(v_cc,0,1,     domain,zdir);
    h_avg_w  = sumToLine(w_fc,0,1,stag_domain,zdir);

    h_avg_rho  = sumToLine(mf_out, 0,1,domain,zdir);
    h_avg_th   = sumToLine(mf_out, 1,1,domain,zdir);
    h_avg_ksgs = sumToLine(mf_out, 2,1,domain,zdir);
    h_avg_uu   = sumToLine(mf_out, 3,1,domain,zdir);
    h_avg_uv   = sumToLine(mf_out, 4,1,domain,zdir);
    h_avg_vv   = sumToLine(mf_out, 5,1,domain,zdir);
    h_avg_uth  = sumToLine(mf_out, 6,1,domain,zdir);
    h_avg_vth  = sumToLine(mf_out, 7,1,domain,zdir);
    h_avg_thth = sumToLine(mf_out, 8,1,domain,zdir);
    h_avg_uiuiu= sumToLine(mf_out, 9,1,stag_domain,zdir);
    h_avg_uiuiv= sumToLine(mf_out,10,1,stag_domain,zdir);
    h_avg_p    = sumToLine(mf_out,11,1,domain,zdir);
    h_avg_pu   = sumToLine(mf_out,12,1,domain,zdir);
    h_avg_pv   = sumToLine(mf_out,13,1,domain,zdir);

    h_avg_uw   = sumToLine(mf_out_stag,0,1,stag_domain,zdir);
    h_avg_vw   = sumToLine(mf_out_stag,1,1,stag_domain,zdir);
    h_avg_ww   = sumToLine(mf_out_stag,2,1,stag_domain,zdir);
    h_avg_wth  = sumToLine(mf_out_stag,3,1,stag_domain,zdir);
    h_avg_uiuiw= sumToLine(mf_out_stag,4,1,stag_domain,zdir);
    h_avg_pw   = sumToLine(mf_out_stag,5,1,stag_domain,zdir);

    // Divide by the total number of cells we are averaging over
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    int unstag_size = h_avg_w.size() - 1; // _un_staggered heights
    for (int k = 0; k < unstag_size; ++k) {
        h_avg_u[k]     /= area_z;
        h_avg_v[k]     /= area_z;
        h_avg_rho[k]   /= area_z;
        h_avg_ksgs[k]  /= area_z;
        h_avg_th[k]    /= area_z;
        h_avg_thth[k]  /= area_z;
        h_avg_uu[k]    /= area_z;
        h_avg_uv[k]    /= area_z;
        h_avg_vv[k]    /= area_z;
        h_avg_uth[k]   /= area_z;
        h_avg_vth[k]   /= area_z;
        h_avg_uiuiu[k] /= area_z;
        h_avg_uiuiv[k] /= area_z;
        h_avg_p[k]     /= area_z;
        h_avg_pu[k]    /= area_z;
        h_avg_pv[k]    /= area_z;
    }

    for (int k = 0; k < unstag_size+1; ++k) { // staggered heights
        h_avg_w[k]     /= area_z;
        h_avg_uw[k]    /= area_z;
        h_avg_vw[k]    /= area_z;
        h_avg_ww[k]    /= area_z;
        h_avg_wth[k]   /= area_z;
        h_avg_uiuiw[k] /= area_z;
        h_avg_pw[k]    /= area_z;
    }
}

void
ERF::derive_stress_profiles_stag(Gpu::HostVector<Real>& h_avg_tau11, Gpu::HostVector<Real>& h_avg_tau12,
                                 Gpu::HostVector<Real>& h_avg_tau13, Gpu::HostVector<Real>& h_avg_tau22,
                                 Gpu::HostVector<Real>& h_avg_tau23, Gpu::HostVector<Real>& h_avg_tau33,
                                 Gpu::HostVector<Real>& h_avg_hfx3,  Gpu::HostVector<Real>& h_avg_diss)
{
    int lev = 0;

    // This will hold the stress tensor components
    MultiFab mf_out(grids[lev], dmap[lev], 8, 0);

    // This will hold Tau13 and Tau23
    MultiFab mf_out_stag(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 2, 0);

    for ( MFIter mfi(mf_out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& fab_arr = mf_out.array(mfi);
        const Array4<Real>& fab_arr_stag = mf_out_stag.array(mfi);

        // NOTE: These are from the last RK stage...
        const Array4<const Real>& tau11_arr = Tau11_lev[lev]->const_array(mfi);
        const Array4<const Real>& tau12_arr = Tau12_lev[lev]->const_array(mfi);
        const Array4<const Real>& tau13_arr = Tau13_lev[lev]->const_array(mfi);
        const Array4<const Real>& tau22_arr = Tau22_lev[lev]->const_array(mfi);
        const Array4<const Real>& tau23_arr = Tau23_lev[lev]->const_array(mfi);
        const Array4<const Real>& tau33_arr = Tau33_lev[lev]->const_array(mfi);

        // These should be re-calculated during ERF_slow_rhs_post
        // -- just vertical SFS kinematic heat flux for now
        //const Array4<const Real>& hfx1_arr = SFS_hfx1_lev[lev]->const_array(mfi);
        //const Array4<const Real>& hfx2_arr = SFS_hfx2_lev[lev]->const_array(mfi);
        const Array4<const Real>& hfx3_arr = SFS_hfx3_lev[lev]->const_array(mfi);
        const Array4<const Real>& diss_arr = SFS_diss_lev[lev]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            fab_arr(i, j, k, 0) = tau11_arr(i,j,k);
            fab_arr(i, j, k, 1) = 0.25 * ( tau12_arr(i,j  ,k) + tau12_arr(i+1,j  ,k)
                                         + tau12_arr(i,j+1,k) + tau12_arr(i+1,j+1,k) );
//          fab_arr(i, j, k, 2) = 0.25 * ( tau13_arr(i,j,k  ) + tau13_arr(i+1,j,k)
//                                       + tau13_arr(i,j,k+1) + tau13_arr(i+1,j,k+1) );
            fab_arr(i, j, k, 3) = tau22_arr(i,j,k);
//          fab_arr(i, j, k, 4) = 0.25 * ( tau23_arr(i,j,k  ) + tau23_arr(i,j+1,k)
//                                       + tau23_arr(i,j,k+1) + tau23_arr(i,j+1,k+1) );
            fab_arr(i, j, k, 5) = tau33_arr(i,j,k);
            fab_arr(i, j, k, 6) =  hfx3_arr(i,j,k);
            fab_arr(i, j, k, 7) =  diss_arr(i,j,k);
        });

        const Box& zbx = mfi.tilebox(IntVect(0,0,1));
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            fab_arr_stag(i,j,k,0) = 0.5*(tau13_arr(i,j,k) + tau13_arr(i+1,j  ,k));
            fab_arr_stag(i,j,k,1) = 0.5*(tau23_arr(i,j,k) + tau23_arr(i  ,j+1,k));
        });
    }

    int zdir = 2;
    auto domain = geom[0].Domain();
    Box stag_domain = domain;
    stag_domain.convert(IntVect(0,0,1));

    h_avg_tau11 = sumToLine(mf_out,0,1,domain,zdir);
    h_avg_tau12 = sumToLine(mf_out,1,1,domain,zdir);
//  h_avg_tau13 = sumToLine(mf_out,2,1,domain,zdir);
    h_avg_tau22 = sumToLine(mf_out,3,1,domain,zdir);
//  h_avg_tau23 = sumToLine(mf_out,4,1,domain,zdir);
    h_avg_tau33 = sumToLine(mf_out,5,1,domain,zdir);
    h_avg_hfx3  = sumToLine(mf_out,6,1,domain,zdir);
    h_avg_diss  = sumToLine(mf_out,7,1,domain,zdir);

    h_avg_tau13 = sumToLine(mf_out_stag,0,1,stag_domain,zdir);
    h_avg_tau23 = sumToLine(mf_out_stag,1,1,stag_domain,zdir);

    int ht_size =  h_avg_tau11.size(); // _un_staggered

    // Divide by the total number of cells we are averaging over
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    for (int k = 0; k < ht_size; ++k) {
        h_avg_tau11[k] /= area_z; h_avg_tau12[k] /= area_z;
        h_avg_tau22[k] /= area_z; h_avg_tau33[k] /= area_z;
        h_avg_hfx3[k] /= area_z;
        h_avg_diss[k] /= area_z;
    }
    for (int k = 0; k < ht_size+1; ++k) { // staggered heights
        h_avg_tau13[k] /= area_z;
        h_avg_tau23[k] /= area_z;
    }
}
