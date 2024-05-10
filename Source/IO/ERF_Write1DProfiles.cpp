#include <iomanip>

#include "ERF.H"
#include "EOS.H"

using namespace amrex;

/**
 * Writes 1-dimensional averaged quantities as profiles to output log files
 * at the given time.
 *
 * @param time Current time
 */
void
ERF::write_1D_profiles (Real time)
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
        Gpu::HostVector<Real> h_avg_qv, h_avg_qc, h_avg_qr, h_avg_wqv, h_avg_wqc, h_avg_wqr, h_avg_qi, h_avg_qs, h_avg_qg;
        Gpu::HostVector<Real> h_avg_uth, h_avg_vth, h_avg_wth, h_avg_thth;
        Gpu::HostVector<Real> h_avg_uu, h_avg_uv, h_avg_uw, h_avg_vv, h_avg_vw, h_avg_ww;
        Gpu::HostVector<Real> h_avg_uiuiu, h_avg_uiuiv, h_avg_uiuiw;
        Gpu::HostVector<Real> h_avg_p, h_avg_pu, h_avg_pv, h_avg_pw;
        Gpu::HostVector<Real> h_avg_tau11, h_avg_tau12, h_avg_tau13, h_avg_tau22, h_avg_tau23, h_avg_tau33;
        Gpu::HostVector<Real> h_avg_sgshfx, h_avg_sgsdiss; // only output tau_{theta,w} and epsilon for now

        if (NumDataLogs() > 1) {
            derive_diag_profiles(h_avg_u, h_avg_v, h_avg_w,
                                 h_avg_rho, h_avg_th, h_avg_ksgs,
                                 h_avg_qv, h_avg_qc, h_avg_qr, h_avg_wqv, h_avg_wqc, h_avg_wqr, h_avg_qi, h_avg_qs, h_avg_qg,
                                 h_avg_uu, h_avg_uv, h_avg_uw, h_avg_vv, h_avg_vw, h_avg_ww,
                                 h_avg_uth, h_avg_vth, h_avg_wth, h_avg_thth,
                                 h_avg_uiuiu, h_avg_uiuiv, h_avg_uiuiw,
                                 h_avg_p, h_avg_pu, h_avg_pv, h_avg_pw);
        }

        if (NumDataLogs() > 3) {
            derive_stress_profiles(h_avg_tau11, h_avg_tau12, h_avg_tau13,
                                   h_avg_tau22, h_avg_tau23, h_avg_tau33,
                                   h_avg_sgshfx,
                                   h_avg_sgsdiss);
        }

        int hu_size =  h_avg_u.size();

        auto const& dx = geom[0].CellSizeArray();
        if (ParallelDescriptor::IOProcessor()) {
            if (NumDataLogs() > 1) {
                std::ostream& data_log1 = DataLog(1);
                if (data_log1.good()) {
                  // Write the quantities at this time
                  for (int k = 0; k < hu_size; k++) {
                      Real z;
                      if (zlevels_stag.size() > 1) {
                          z = 0.5 * (zlevels_stag[k] + zlevels_stag[k+1]);
                      } else {
                          z = (k + 0.5)* dx[2];
                      }
                      data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_u[k]  << " " << h_avg_v[k] << " " << h_avg_w[k] << " "
                                << h_avg_rho[k] << " " << h_avg_th[k] << " " << h_avg_ksgs[k] << " "
                                << h_avg_qv[k] << " " << h_avg_qc[k] << " " << h_avg_qr[k] << " "
                                << h_avg_qi[k] << " " << h_avg_qs[k] << " " << h_avg_qg[k]
                                << std::endl;
                  } // loop over z
                } // if good
            } // NumDataLogs

            if (NumDataLogs() > 2) {
                std::ostream& data_log2 = DataLog(2);
                if (data_log2.good()) {
                  // Write the perturbational quantities at this time
                  for (int k = 0; k < hu_size; k++) {
                      Real z;
                      if (zlevels_stag.size() > 1) {
                          z = 0.5 * (zlevels_stag[k] + zlevels_stag[k+1]);
                      } else {
                          z = (k + 0.5)* dx[2];
                      }
                      data_log2 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_uu[k]   - h_avg_u[k]*h_avg_u[k]  << " "
                                << h_avg_uv[k]   - h_avg_u[k]*h_avg_v[k]  << " "
                                << h_avg_uw[k]   - h_avg_u[k]*h_avg_w[k]  << " "
                                << h_avg_vv[k]   - h_avg_v[k]*h_avg_v[k]  << " "
                                << h_avg_vw[k]   - h_avg_v[k]*h_avg_w[k]  << " "
                                << h_avg_ww[k]   - h_avg_w[k]*h_avg_w[k]  << " "
                                << h_avg_uth[k]  - h_avg_u[k]*h_avg_th[k] << " "
                                << h_avg_vth[k]  - h_avg_v[k]*h_avg_th[k] << " "
                                << h_avg_wth[k]  - h_avg_w[k]*h_avg_th[k] << " "
                                << h_avg_thth[k] - h_avg_th[k]*h_avg_th[k] << " "
                                // Note: <u'_i u'_i u'_j> =   <u_i u_i u_j>
                                //                        -   <u_i u_i> * <u_j>
                                //                        - 2*<u_i> * <u_i u_j>
                                //                        + 2*<u_i>*<u_i> * <u_j>
                                << h_avg_uiuiu[k]
                                 - (h_avg_uu[k] + h_avg_vv[k] + h_avg_ww[k])*h_avg_u[k]
                                 - 2*(h_avg_u[k]*h_avg_uu[k] + h_avg_v[k]*h_avg_uv[k] + h_avg_w[k]*h_avg_uw[k])
                                 + 2*(h_avg_u[k]*h_avg_u[k] + h_avg_v[k]*h_avg_v[k] + h_avg_w[k]*h_avg_w[k])*h_avg_u[k]
                                   << " " // (u'_i u'_i)u'
                                << h_avg_uiuiv[k]
                                 - (h_avg_uu[k] + h_avg_vv[k] + h_avg_ww[k])*h_avg_v[k]
                                 - 2*(h_avg_u[k]*h_avg_uv[k] + h_avg_v[k]*h_avg_vv[k] + h_avg_w[k]*h_avg_vw[k])
                                 + 2*(h_avg_u[k]*h_avg_u[k] + h_avg_v[k]*h_avg_v[k] + h_avg_w[k]*h_avg_w[k])*h_avg_v[k]
                                   << " " // (u'_i u'_i)v'
                                << h_avg_uiuiw[k]
                                 - (h_avg_uu[k] + h_avg_vv[k] + h_avg_ww[k])*h_avg_w[k]
                                 - 2*(h_avg_u[k]*h_avg_uw[k] + h_avg_v[k]*h_avg_vw[k] + h_avg_w[k]*h_avg_ww[k])
                                 + 2*(h_avg_u[k]*h_avg_u[k] + h_avg_v[k]*h_avg_v[k] + h_avg_w[k]*h_avg_w[k])*h_avg_w[k]
                                   << " " // (u'_i u'_i)w'
                                << h_avg_pu[k]   - h_avg_p[k]*h_avg_u[k] << " "
                                << h_avg_pv[k]   - h_avg_p[k]*h_avg_v[k] << " "
                                << h_avg_pw[k]   - h_avg_p[k]*h_avg_w[k] << " "
                                << h_avg_wqv[k]   - h_avg_qv[k]*h_avg_w[k] << " "
                                << h_avg_wqc[k]   - h_avg_qc[k]*h_avg_w[k] << " "
                                << h_avg_wqr[k]   - h_avg_qr[k]*h_avg_w[k]
                                << std::endl;
                  } // loop over z
                } // if good
            } // NumDataLogs

            if (NumDataLogs() > 3) {
                std::ostream& data_log3 = DataLog(3);
                if (data_log3.good()) {
                  // Write the average stresses
                  for (int k = 0; k < hu_size; k++) {
                      Real z;
                      if (zlevels_stag.size() > 1) {
                          z = 0.5 * (zlevels_stag[k] + zlevels_stag[k+1]);
                      } else {
                          z = (k + 0.5)* dx[2];
                      }
                      data_log3 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_tau11[k] << " " << h_avg_tau12[k] << " " << h_avg_tau13[k] << " "
                                << h_avg_tau22[k] << " " << h_avg_tau23[k] << " " << h_avg_tau33[k] << " "
                                << h_avg_sgshfx[k] << " " << h_avg_sgsdiss[k]
                                << std::endl;
                  } // loop over z
                } // if good
            } // NumDataLogs
        } // if IOProcessor
    } // if verbose
}

/**
 * Computes the profiles for diagnostic quantities.
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

void ERF::derive_diag_profiles(Gpu::HostVector<Real>& h_avg_u   , Gpu::HostVector<Real>& h_avg_v  , Gpu::HostVector<Real>& h_avg_w,
                          Gpu::HostVector<Real>& h_avg_rho , Gpu::HostVector<Real>& h_avg_th , Gpu::HostVector<Real>& h_avg_ksgs,
                          Gpu::HostVector<Real>& h_avg_qv  , Gpu::HostVector<Real>& h_avg_qc , Gpu::HostVector<Real>& h_avg_qr,
                          Gpu::HostVector<Real>& h_avg_wqv , Gpu::HostVector<Real>& h_avg_wqc, Gpu::HostVector<Real>& h_avg_wqr,
                          Gpu::HostVector<Real>& h_avg_qi  , Gpu::HostVector<Real>& h_avg_qs , Gpu::HostVector<Real>& h_avg_qg,
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

    // This will hold rho, theta, ksgs, uu, uv, uw, vv, vw, ww, uth, vth, wth,
    //                  0      1     2   3   4   5   6   7   8    9   10   11
    //                thth, uiuiu, uiuiv, uiuiw, p, pu, pv, pw, qv, qc, qr, wqv, wqc, wqr, qi, qs, qg
    //                  12     13     14     15 16  17  18  19  20  21  22  23   24   25   26  27  28
    MultiFab mf_out(grids[lev], dmap[lev], 29, 0);

    MultiFab mf_vels(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);

    MultiFab  u_cc(mf_vels, make_alias, 0, 1); // u at cell centers
    MultiFab  v_cc(mf_vels, make_alias, 1, 1); // v at cell centers
    MultiFab  w_cc(mf_vels, make_alias, 2, 1); // w at cell centers

    average_face_to_cellcenter(mf_vels,0,
            Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});

    int zdir = 2;
    auto domain = geom[0].Domain();

    // Sum in the horizontal plane
    h_avg_u  = sumToLine(mf_vels ,0,1,domain,zdir);
    h_avg_v  = sumToLine(mf_vels ,1,1,domain,zdir);
    h_avg_w  = sumToLine(mf_vels ,2,1,domain,zdir);

    int hu_size =  h_avg_u.size();

    // Divide by the total number of cells we are averaging over
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    for (int k = 0; k < hu_size; ++k) {
        h_avg_u[k] /= area_z; h_avg_v[k] /= area_z; h_avg_w[k] /= area_z;
    }

    int nvars = vars_new[lev][Vars::cons].nComp();
    MultiFab mf_cons(vars_new[lev][Vars::cons], make_alias, 0, nvars);

    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component

    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);

    for ( MFIter mfi(mf_cons,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& fab_arr = mf_out.array(mfi);
        const Array4<Real>& u_cc_arr =  u_cc.array(mfi);
        const Array4<Real>& v_cc_arr =  v_cc.array(mfi);
        const Array4<Real>& w_cc_arr =  w_cc.array(mfi);
        const Array4<Real>& cons_arr = mf_cons.array(mfi);
        const Array4<Real>&   p0_arr = p_hse.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
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
            fab_arr(i, j, k, 5) = u_cc_arr(i,j,k) * w_cc_arr(i,j,k);   // u*w
            fab_arr(i, j, k, 6) = v_cc_arr(i,j,k) * v_cc_arr(i,j,k);   // v*v
            fab_arr(i, j, k, 7) = v_cc_arr(i,j,k) * w_cc_arr(i,j,k);   // v*w
            fab_arr(i, j, k, 8) = w_cc_arr(i,j,k) * w_cc_arr(i,j,k);   // w*w
            fab_arr(i, j, k, 9) = u_cc_arr(i,j,k) * theta; // u*th
            fab_arr(i, j, k,10) = v_cc_arr(i,j,k) * theta; // v*th
            fab_arr(i, j, k,11) = w_cc_arr(i,j,k) * theta; // w*th
            fab_arr(i, j, k,12) = theta * theta;           // th*th
            Real uiui = fab_arr(i,j,k,3) + fab_arr(i,j,k,6) + fab_arr(i,j,k,8);
            fab_arr(i, j, k,13) = uiui * u_cc_arr(i,j,k);   // (ui*ui)*u
            fab_arr(i, j, k,14) = uiui * v_cc_arr(i,j,k);   // (ui*ui)*v
            fab_arr(i, j, k,15) = uiui * w_cc_arr(i,j,k);   // (ui*ui)*w

            if (!use_moisture) {
                Real p = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp));
                p -= p0_arr(i,j,k);
                fab_arr(i, j, k,16) = p;                       // p'
                fab_arr(i, j, k,17) = p * u_cc_arr(i,j,k);     // p'*u
                fab_arr(i, j, k,18) = p * v_cc_arr(i,j,k);     // p'*v
                fab_arr(i, j, k,19) = p * w_cc_arr(i,j,k);     // p'*w
                fab_arr(i, j, k,20) = 0.;  // qv
                fab_arr(i, j, k,21) = 0.;  // qc
                fab_arr(i, j, k,22) = 0.;  // qr
                fab_arr(i, j, k,23) = 0.;  // w*qv
                fab_arr(i, j, k,24) = 0.;  // w*qc
                fab_arr(i, j, k,25) = 0.;  // w*qr
                fab_arr(i, j, k,26) = 0.;  // qi
                fab_arr(i, j, k,27) = 0.;  // qs
                fab_arr(i, j, k,28) = 0.;  // qg
            }
        });
    } // mfi

    if (use_moisture)
    {
        int RhoQr_comp;
        int n_qstate = micro->Get_Qstate_Size();
        if (n_qstate > 3) {
            RhoQr_comp = RhoQ4_comp;
        } else {
            RhoQr_comp = RhoQ3_comp;
        }

        for ( MFIter mfi(mf_cons,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Array4<Real>& fab_arr  = mf_out.array(mfi);
            const Array4<Real>& cons_arr = mf_cons.array(mfi);
            const Array4<Real>& u_cc_arr =  u_cc.array(mfi);
            const Array4<Real>& v_cc_arr =  v_cc.array(mfi);
            const Array4<Real>& w_cc_arr =  w_cc.array(mfi);
            const Array4<Real>&   p0_arr = p_hse.array(mfi);
            const Array4<Real>&   qv_arr = qmoist[0][0]->array(mfi); // TODO: Is this written only on lev 0?

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real p = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp), qv_arr(i,j,k));

                p -= p0_arr(i,j,k);
                fab_arr(i, j, k,16) = p;                       // p'
                fab_arr(i, j, k,17) = p * u_cc_arr(i,j,k);     // p'*u
                fab_arr(i, j, k,18) = p * v_cc_arr(i,j,k);     // p'*v
                fab_arr(i, j, k,19) = p * w_cc_arr(i,j,k);     // p'*w
                fab_arr(i, j, k,20) = cons_arr(i,j,k,RhoQ1_comp) / cons_arr(i,j,k,Rho_comp);  // qv
                fab_arr(i, j, k,21) = cons_arr(i,j,k,RhoQ2_comp) / cons_arr(i,j,k,Rho_comp);  // qc
                fab_arr(i, j, k,22) = cons_arr(i,j,k,RhoQr_comp) / cons_arr(i,j,k,Rho_comp);  // qr
                fab_arr(i, j, k,23) = w_cc_arr(i,j,k) * cons_arr(i,j,k,RhoQ1_comp) / cons_arr(i,j,k,Rho_comp);  // w*qv
                fab_arr(i, j, k,24) = w_cc_arr(i,j,k) * cons_arr(i,j,k,RhoQ2_comp) / cons_arr(i,j,k,Rho_comp);  // w*qc
                fab_arr(i, j, k,25) = w_cc_arr(i,j,k) * cons_arr(i,j,k,RhoQr_comp) / cons_arr(i,j,k,Rho_comp);  // w*qr
                if (n_qstate > 3) {
                    fab_arr(i, j, k,26) = cons_arr(i,j,k,RhoQ3_comp) / cons_arr(i,j,k,Rho_comp);  // qi
                    fab_arr(i, j, k,27) = cons_arr(i,j,k,RhoQ5_comp) / cons_arr(i,j,k,Rho_comp);  // qs
                    fab_arr(i, j, k,28) = cons_arr(i,j,k,RhoQ6_comp) / cons_arr(i,j,k,Rho_comp);  // qg
                } else {
                    fab_arr(i, j, k,26) = 0.0;  // qi
                    fab_arr(i, j, k,27) = 0.0;  // qs
                    fab_arr(i, j, k,28) = 0.0;  // qg
                }
            });
        } // mfi
    } // use_moisture

    h_avg_rho  = sumToLine(mf_out, 0,1,domain,zdir);
    h_avg_th   = sumToLine(mf_out, 1,1,domain,zdir);
    h_avg_ksgs = sumToLine(mf_out, 2,1,domain,zdir);
    h_avg_uu   = sumToLine(mf_out, 3,1,domain,zdir);
    h_avg_uv   = sumToLine(mf_out, 4,1,domain,zdir);
    h_avg_uw   = sumToLine(mf_out, 5,1,domain,zdir);
    h_avg_vv   = sumToLine(mf_out, 6,1,domain,zdir);
    h_avg_vw   = sumToLine(mf_out, 7,1,domain,zdir);
    h_avg_ww   = sumToLine(mf_out, 8,1,domain,zdir);
    h_avg_uth  = sumToLine(mf_out, 9,1,domain,zdir);
    h_avg_vth  = sumToLine(mf_out,10,1,domain,zdir);
    h_avg_wth  = sumToLine(mf_out,11,1,domain,zdir);
    h_avg_thth = sumToLine(mf_out,12,1,domain,zdir);
    h_avg_uiuiu= sumToLine(mf_out,13,1,domain,zdir);
    h_avg_uiuiv= sumToLine(mf_out,14,1,domain,zdir);
    h_avg_uiuiw= sumToLine(mf_out,15,1,domain,zdir);
    h_avg_p    = sumToLine(mf_out,16,1,domain,zdir);
    h_avg_pu   = sumToLine(mf_out,17,1,domain,zdir);
    h_avg_pv   = sumToLine(mf_out,18,1,domain,zdir);
    h_avg_pw   = sumToLine(mf_out,19,1,domain,zdir);
    h_avg_qv   = sumToLine(mf_out,20,1,domain,zdir);
    h_avg_qc   = sumToLine(mf_out,21,1,domain,zdir);
    h_avg_qr   = sumToLine(mf_out,22,1,domain,zdir);
    h_avg_wqv  = sumToLine(mf_out,23,1,domain,zdir);
    h_avg_wqc  = sumToLine(mf_out,24,1,domain,zdir);
    h_avg_wqr  = sumToLine(mf_out,25,1,domain,zdir);
    h_avg_qi   = sumToLine(mf_out,26,1,domain,zdir);
    h_avg_qs   = sumToLine(mf_out,27,1,domain,zdir);
    h_avg_qg   = sumToLine(mf_out,28,1,domain,zdir);

    // Divide by the total number of cells we are averaging over
    int h_avg_u_size = static_cast<int>(h_avg_u.size());
    for (int k = 0; k < h_avg_u_size; ++k) {
        h_avg_rho[k] /= area_z;  h_avg_ksgs[k] /= area_z;
        h_avg_th[k]  /= area_z;  h_avg_thth[k] /= area_z;
        h_avg_uu[k]  /= area_z;  h_avg_uv[k]   /= area_z;  h_avg_uw[k]  /= area_z;
        h_avg_vv[k]  /= area_z;  h_avg_vw[k]   /= area_z;  h_avg_ww[k]  /= area_z;
        h_avg_uth[k] /= area_z;  h_avg_vth[k]  /= area_z;  h_avg_wth[k] /= area_z;
        h_avg_uiuiu[k] /= area_z;
        h_avg_uiuiv[k] /= area_z;
        h_avg_uiuiw[k] /= area_z;
        h_avg_p[k]   /= area_z;
        h_avg_pu[k]  /= area_z;  h_avg_pv[k]   /= area_z;  h_avg_pw[k]  /= area_z;
        h_avg_qv[k]  /= area_z;  h_avg_qc[k]   /= area_z;  h_avg_qr[k]  /= area_z;
        h_avg_wqv[k] /= area_z;  h_avg_wqc[k]  /= area_z;  h_avg_wqr[k] /= area_z;
        h_avg_qi[k]  /= area_z;  h_avg_qs[k]   /= area_z;  h_avg_qg[k]  /= area_z;
    }
}
void
ERF::derive_stress_profiles (Gpu::HostVector<Real>& h_avg_tau11, Gpu::HostVector<Real>& h_avg_tau12,
                             Gpu::HostVector<Real>& h_avg_tau13, Gpu::HostVector<Real>& h_avg_tau22,
                             Gpu::HostVector<Real>& h_avg_tau23, Gpu::HostVector<Real>& h_avg_tau33,
                             Gpu::HostVector<Real>& h_avg_hfx3,  Gpu::HostVector<Real>& h_avg_diss)
{
    int lev = 0;

    // This will hold the stress tensor components
    MultiFab mf_out(grids[lev], dmap[lev], 8, 0);

    for ( MFIter mfi(mf_out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& fab_arr = mf_out.array(mfi);

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
            fab_arr(i, j, k, 2) = 0.25 * ( tau13_arr(i,j,k  ) + tau13_arr(i+1,j,k)
                                         + tau13_arr(i,j,k+1) + tau13_arr(i+1,j,k+1) );
            fab_arr(i, j, k, 3) = tau22_arr(i,j,k);
            fab_arr(i, j, k, 4) = 0.25 * ( tau23_arr(i,j,k  ) + tau23_arr(i,j+1,k)
                                         + tau23_arr(i,j,k+1) + tau23_arr(i,j+1,k+1) );
            fab_arr(i, j, k, 5) = tau33_arr(i,j,k);
            fab_arr(i, j, k, 6) =  hfx3_arr(i,j,k);
            fab_arr(i, j, k, 7) =  diss_arr(i,j,k);
        });
    }

    int zdir = 2;
    auto domain = geom[0].Domain();

    h_avg_tau11 = sumToLine(mf_out,0,1,domain,zdir);
    h_avg_tau12 = sumToLine(mf_out,1,1,domain,zdir);
    h_avg_tau13 = sumToLine(mf_out,2,1,domain,zdir);
    h_avg_tau22 = sumToLine(mf_out,3,1,domain,zdir);
    h_avg_tau23 = sumToLine(mf_out,4,1,domain,zdir);
    h_avg_tau33 = sumToLine(mf_out,5,1,domain,zdir);
    h_avg_hfx3  = sumToLine(mf_out,6,1,domain,zdir);
    h_avg_diss  = sumToLine(mf_out,7,1,domain,zdir);

    int ht_size =  h_avg_tau11.size();

    // Divide by the total number of cells we are averaging over
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    for (int k = 0; k < ht_size; ++k) {
        h_avg_tau11[k] /= area_z; h_avg_tau12[k] /= area_z; h_avg_tau13[k] /= area_z;
        h_avg_tau22[k] /= area_z; h_avg_tau23[k] /= area_z; h_avg_tau33[k] /= area_z;
        h_avg_hfx3[k] /= area_z;
        h_avg_diss[k] /= area_z;
    }
}
