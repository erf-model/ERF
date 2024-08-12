#include <iomanip>

#include "ERF.H"
#include "EOS.H"

using namespace amrex;

/**
 * Writes 1-dimensional averaged quantities as profiles to output log files
 * at the given time.
 *
 * Quantities are output at their native grid locations. Therefore, w and
 * associated flux quantities <(•)'w'>, tau13, and tau23 (where '•' includes
 * u, v, p, theta, ...) will be output at staggered heights (i.e., coincident
 * with z faces) rather than cell-center heights to avoid performing additional
 * averaging. Unstaggered (i.e., cell-centered) quantities are output alongside
 * staggered quantities at the lower cell faces in the log file; these
 * quantities will have a zero value at the big end, corresponding to k=Nz+1.
 *
 * @param time Current time
 */
void
ERF::write_1D_profiles_stag (Real time)
{
    BL_PROFILE("ERF::write_1D_profiles()");

    int datwidth = 14;
    int datprecision = 6;
    int timeprecision = 13; // e.g., 1-yr LES: 31,536,000 s with dt ~ 0.01 ==> min prec = 10

    if (NumDataLogs() > 1)
    {
        // Define the 1d arrays we will need
        Gpu::HostVector<Real> h_avg_u, h_avg_v, h_avg_w;
        Gpu::HostVector<Real> h_avg_rho, h_avg_th, h_avg_ksgs, h_avg_Kmv, h_avg_Khv;
        Gpu::HostVector<Real> h_avg_qv, h_avg_qc, h_avg_qr, h_avg_wqv, h_avg_wqc, h_avg_wqr, h_avg_qi, h_avg_qs, h_avg_qg;
        Gpu::HostVector<Real> h_avg_wthv;
        Gpu::HostVector<Real> h_avg_uth, h_avg_vth, h_avg_wth, h_avg_thth;
        Gpu::HostVector<Real> h_avg_uu, h_avg_uv, h_avg_uw, h_avg_vv, h_avg_vw, h_avg_ww;
        Gpu::HostVector<Real> h_avg_uiuiu, h_avg_uiuiv, h_avg_uiuiw;
        Gpu::HostVector<Real> h_avg_p, h_avg_pu, h_avg_pv, h_avg_pw;
        Gpu::HostVector<Real> h_avg_tau11, h_avg_tau12, h_avg_tau13, h_avg_tau22, h_avg_tau23, h_avg_tau33;
        Gpu::HostVector<Real> h_avg_sgshfx, h_avg_sgsq1fx, h_avg_sgsq2fx, h_avg_sgsdiss; // only output tau_{theta,w} and epsilon for now

        if (NumDataLogs() > 1) {
            derive_diag_profiles_stag(time,
                                      h_avg_u, h_avg_v, h_avg_w,
                                      h_avg_rho, h_avg_th, h_avg_ksgs,
                                      h_avg_Kmv, h_avg_Khv,
                                      h_avg_qv, h_avg_qc, h_avg_qr,
                                      h_avg_wqv, h_avg_wqc, h_avg_wqr,
                                      h_avg_qi, h_avg_qs, h_avg_qg,
                                      h_avg_uu, h_avg_uv, h_avg_uw, h_avg_vv, h_avg_vw, h_avg_ww,
                                      h_avg_uth, h_avg_vth, h_avg_wth, h_avg_thth,
                                      h_avg_uiuiu, h_avg_uiuiv, h_avg_uiuiw,
                                      h_avg_p, h_avg_pu, h_avg_pv, h_avg_pw,
                                      h_avg_wthv);
        }

        if (NumDataLogs() > 3 && time > 0.) {
            derive_stress_profiles_stag(h_avg_tau11, h_avg_tau12, h_avg_tau13,
                                        h_avg_tau22, h_avg_tau23, h_avg_tau33,
                                        h_avg_sgshfx, h_avg_sgsq1fx, h_avg_sgsq2fx,
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
                      Real z = (zlevels_stag.size() > 1) ? zlevels_stag[k] : k * dx[2];
                      data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_u[k]   << " " << h_avg_v[k]   << " " << h_avg_w[k]     << " "
                                << h_avg_rho[k] << " " << h_avg_th[k]  << " " << h_avg_ksgs[k] << " "
                                << h_avg_Kmv[k] << " " << h_avg_Khv[k] << " "
                                << h_avg_qv[k]  << " " << h_avg_qc[k]  << " " << h_avg_qr[k]    << " "
                                << h_avg_qi[k]  << " " << h_avg_qs[k]  << " " << h_avg_qg[k]
                                << std::endl;
                  } // loop over z
                  // Write top face values
                  Real z = (zlevels_stag.size() > 1) ? zlevels_stag[unstag_size] : unstag_size * dx[2];
                  data_log1 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                            << 0 << " " << 0 << " " << h_avg_w[unstag_size+1] << " "
                            << 0 << " " << 0 << " " << 0 << " " // rho, theta, ksgs
                            << 0 << " " << 0 << " "             // Kmv, Khv
                            << 0 << " " << 0 << " " << 0 << " " // qv, qc, qr
                            << 0 << " " << 0 << " " << 0        // qi, qs, qg
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
                            << h_avg_pu[0]   - h_avg_p[0]*h_avg_u[0]   << " " // p'u'
                            << h_avg_pv[0]   - h_avg_p[0]*h_avg_v[0]   << " " // p'v'
                            << 0                                       << " " // p'w'
                            << 0                                       << " " // qv'w'
                            << 0                                       << " " // qc'w'
                            << 0                                       << " " // qr'w'
                            << 0                                              // thv'w'
                            << std::endl;

                  // For internal values, interpolate scalar quantities to faces
                  for (int k = 1; k < unstag_size; k++) {
                      Real z = (zlevels_stag.size() > 1) ? zlevels_stag[k] : k * dx[2];
                      Real uface  = 0.5*(h_avg_u[k]  + h_avg_u[k-1]);
                      Real vface  = 0.5*(h_avg_v[k]  + h_avg_v[k-1]);
                      Real thface = 0.5*(h_avg_th[k] + h_avg_th[k-1]);
                      Real pface  = 0.5*(h_avg_p[k]  + h_avg_p[k-1]);
                      Real qvface = 0.5*(h_avg_qv[k] + h_avg_qv[k-1]);
                      Real qcface = 0.5*(h_avg_qc[k] + h_avg_qc[k-1]);
                      Real qrface = 0.5*(h_avg_qr[k] + h_avg_qr[k-1]);
                      Real uuface = 0.5*(h_avg_uu[k] + h_avg_uu[k-1]);
                      Real vvface = 0.5*(h_avg_vv[k] + h_avg_vv[k-1]);
                      Real thvface = thface * (1 + 0.61*qvface - qcface - qrface);
                      w_cc   = 0.5*(h_avg_w[k-1]  + h_avg_w[k]);
                      uw_cc  = 0.5*(h_avg_uw[k-1] + h_avg_uw[k]);
                      vw_cc  = 0.5*(h_avg_vw[k-1] + h_avg_vw[k]);
                      ww_cc  = 0.5*(h_avg_ww[k-1] + h_avg_ww[k]);
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
                                << h_avg_pu[k]   - h_avg_p[k]*h_avg_u[k]   << " " // cell-centered p'u'
                                << h_avg_pv[k]   - h_avg_p[k]*h_avg_v[k]   << " " // cell-centered p'v'
                                << h_avg_pw[k]   -      pface*h_avg_w[k]   << " " // face-centered p'w'
                                << h_avg_wqv[k]  -     qvface*h_avg_w[k]   << " "
                                << h_avg_wqc[k]  -     qcface*h_avg_w[k]   << " "
                                << h_avg_wqr[k]  -     qrface*h_avg_w[k]   << " "
                                << h_avg_wthv[k] -    thvface*h_avg_w[k]
                                << std::endl;
                  } // loop over z

                  // Write top face values, extrapolating scalar quantities
                  const int k = unstag_size;
                  Real uface  = 1.5*h_avg_u[k-1]  - 0.5*h_avg_u[k-2];
                  Real vface  = 1.5*h_avg_v[k-1]  - 0.5*h_avg_v[k-2];
                  Real thface = 1.5*h_avg_th[k-1] - 0.5*h_avg_th[k-2];
                  Real pface  = 1.5*h_avg_p[k-1]  - 0.5*h_avg_p[k-2];
                  Real qvface = 1.5*h_avg_qv[k-1] - 0.5*h_avg_qv[k-2];
                  Real qcface = 1.5*h_avg_qc[k-1] - 0.5*h_avg_qc[k-2];
                  Real qrface = 1.5*h_avg_qr[k-1] - 0.5*h_avg_qr[k-2];
                  Real uuface = 1.5*h_avg_uu[k-1] - 0.5*h_avg_uu[k-2];
                  Real vvface = 1.5*h_avg_vv[k-1] - 0.5*h_avg_vv[k-2];
                  Real thvface = thface * (1 + 0.61*qvface - qcface - qrface);
                  Real z = (zlevels_stag.size() > 1) ? zlevels_stag[unstag_size] : unstag_size * dx[2];
                  data_log2 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
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
                            << h_avg_pw[k]   -      pface*h_avg_w[k] << " " // pw'
                            << h_avg_wqv[k]  -     qvface*h_avg_w[k] << " "
                            << h_avg_wqc[k]  -     qcface*h_avg_w[k] << " "
                            << h_avg_wqr[k]  -     qrface*h_avg_w[k] << " "
                            << h_avg_wthv[k] -    thvface*h_avg_w[k]
                            << std::endl;
                } // if good
            } // NumDataLogs

            if (NumDataLogs() > 3 && time > 0.) {
                std::ostream& data_log3 = DataLog(3);
                if (data_log3.good()) {
                  // Write the average stresses
                  for (int k = 0; k < unstag_size; k++) {
                      Real z = (zlevels_stag.size() > 1) ? zlevels_stag[k] : k * dx[2];
                      data_log3 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                                << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                                << h_avg_tau11[k]  << " " << h_avg_tau12[k] << " " << h_avg_tau13[k] << " "
                                << h_avg_tau22[k]  << " " << h_avg_tau23[k] << " " << h_avg_tau33[k] << " "
                                << h_avg_sgshfx[k] << " "
                                << h_avg_sgsq1fx[k] << " " << h_avg_sgsq2fx[k] << " "
                                << h_avg_sgsdiss[k]
                                << std::endl;
                  } // loop over z
                  // Write top face values
                  Real NANval = 0.0;
                  Real z = (zlevels_stag.size() > 1) ? zlevels_stag[unstag_size] : unstag_size * dx[2];
                  data_log3 << std::setw(datwidth) << std::setprecision(timeprecision) << time << " "
                            << std::setw(datwidth) << std::setprecision(datprecision) << z << " "
                            << NANval << " " << NANval << " " << h_avg_tau13[unstag_size] << " "
                            << NANval << " " << h_avg_tau23[unstag_size] << " " << NANval << " "
                            << h_avg_sgshfx[unstag_size] << " "
                            << h_avg_sgsq1fx[unstag_size] << " " << h_avg_sgsq2fx[unstag_size] << " "
                            << NANval
                            << std::endl;
                } // if good
            } // if (NumDataLogs() > 3)
        } // if IOProcessor
    } // if (NumDataLogs() > 1)
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
ERF::derive_diag_profiles_stag (Real /*time*/,
                                Gpu::HostVector<Real>& h_avg_u   , Gpu::HostVector<Real>& h_avg_v  , Gpu::HostVector<Real>& h_avg_w,
                                Gpu::HostVector<Real>& h_avg_rho , Gpu::HostVector<Real>& h_avg_th , Gpu::HostVector<Real>& h_avg_ksgs,
                                Gpu::HostVector<Real>& h_avg_Kmv , Gpu::HostVector<Real>& h_avg_Khv,
                                Gpu::HostVector<Real>& h_avg_qv  , Gpu::HostVector<Real>& h_avg_qc , Gpu::HostVector<Real>& h_avg_qr,
                                Gpu::HostVector<Real>& h_avg_wqv , Gpu::HostVector<Real>& h_avg_wqc, Gpu::HostVector<Real>& h_avg_wqr,
                                Gpu::HostVector<Real>& h_avg_qi  , Gpu::HostVector<Real>& h_avg_qs , Gpu::HostVector<Real>& h_avg_qg,
                                Gpu::HostVector<Real>& h_avg_uu  , Gpu::HostVector<Real>& h_avg_uv , Gpu::HostVector<Real>& h_avg_uw,
                                Gpu::HostVector<Real>& h_avg_vv  , Gpu::HostVector<Real>& h_avg_vw , Gpu::HostVector<Real>& h_avg_ww,
                                Gpu::HostVector<Real>& h_avg_uth , Gpu::HostVector<Real>& h_avg_vth, Gpu::HostVector<Real>& h_avg_wth,
                                Gpu::HostVector<Real>& h_avg_thth,
                                Gpu::HostVector<Real>& h_avg_uiuiu, Gpu::HostVector<Real>& h_avg_uiuiv, Gpu::HostVector<Real>& h_avg_uiuiw,
                                Gpu::HostVector<Real>& h_avg_p,
                                Gpu::HostVector<Real>& h_avg_pu  , Gpu::HostVector<Real>& h_avg_pv , Gpu::HostVector<Real>& h_avg_pw,
                                Gpu::HostVector<Real>& h_avg_wthv)
{
    // We assume that this is always called at level 0
    int lev = 0;

    bool l_use_kturb = ((solverChoice.turbChoice[lev].les_type != LESType::None) ||
                        (solverChoice.turbChoice[lev].pbl_type != PBLType::None));
    bool l_use_KE   = (solverChoice.turbChoice[lev].les_type == LESType::Deardorff);
    bool l_use_QKE  = solverChoice.turbChoice[lev].use_QKE && solverChoice.turbChoice[lev].advect_QKE;

    // Note: "uiui" == u_i*u_i = u*u + v*v + w*w
    // This will hold rho, theta, ksgs, Kmh, Kmv, uu, uv, vv, uth, vth,
    //       indices:   0      1     2    3    4   5   6   7    8    9
    //                thth, uiuiu, uiuiv, p, pu, pv, qv, qc, qr, qi, qs, qg
    //                  10     11     12 13  14  15  16  17  18  19  20  21
    MultiFab mf_out(grids[lev], dmap[lev], 22, 0);

    // This will hold uw, vw, ww, wth, uiuiw, pw, wqv, wqc, wqr, wthv
    //       indices:  0   1   2    3      4   5    6    7    8     9
    MultiFab mf_out_stag(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 10, 0);

    // This is only used to average u and v; w is not averaged to cell centers
    MultiFab mf_vels(grids[lev], dmap[lev], 2, 0);

    MultiFab  u_cc(mf_vels, make_alias, 0, 1); // u at cell centers
    MultiFab  v_cc(mf_vels, make_alias, 1, 1); // v at cell centers
    MultiFab  w_fc(vars_new[lev][Vars::zvel], make_alias, 0, 1); // w at face centers (staggered)

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
        const Array4<const Real>& eta_arr = (l_use_kturb) ? eddyDiffs_lev[lev]->const_array(mfi) :
                                                            Array4<const Real>{};

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_cc_arr(i,j,k) = 0.5 * (u_arr(i,j,k) + u_arr(i+1,j  ,k));
            v_cc_arr(i,j,k) = 0.5 * (v_arr(i,j,k) + v_arr(i  ,j+1,k));

            Real theta = cons_arr(i,j,k,RhoTheta_comp) / cons_arr(i,j,k,Rho_comp);
            fab_arr(i, j, k, 0) = cons_arr(i,j,k,Rho_comp);
            fab_arr(i, j, k, 1) = theta;
            Real ksgs = 0.0;
            if (l_use_KE) {
                ksgs = cons_arr(i,j,k,RhoKE_comp) / cons_arr(i,j,k,Rho_comp);
            } else if (l_use_QKE) {
                ksgs = cons_arr(i,j,k,RhoQKE_comp) / cons_arr(i,j,k,Rho_comp);
            }
            fab_arr(i, j, k, 2) = ksgs;
            if (l_use_kturb) {
                fab_arr(i, j, k, 3) = eta_arr(i,j,k,EddyDiff::Mom_v); // Kmv
                fab_arr(i, j, k, 4) = eta_arr(i,j,k,EddyDiff::Theta_v); // Khv
            } else {
                fab_arr(i, j, k, 3) = 0.0;
                fab_arr(i, j, k, 4) = 0.0;
            }
            fab_arr(i, j, k, 5) = u_cc_arr(i,j,k) * u_cc_arr(i,j,k);   // u*u
            fab_arr(i, j, k, 6) = u_cc_arr(i,j,k) * v_cc_arr(i,j,k);   // u*v
            fab_arr(i, j, k, 7) = v_cc_arr(i,j,k) * v_cc_arr(i,j,k);   // v*v
            fab_arr(i, j, k, 8) = u_cc_arr(i,j,k) * theta;             // u*th
            fab_arr(i, j, k, 9) = v_cc_arr(i,j,k) * theta;             // v*th
            fab_arr(i, j, k,10) = theta * theta;                       // th*th

            Real wcc = 0.5 * (w_fc_arr(i,j,k) + w_fc_arr(i,j,k+1));

            // if the number of fields is changed above, then be sure to update
            // the following def!
            Real uiui = fab_arr(i,j,k,5) + fab_arr(i,j,k,7) + wcc*wcc;
            fab_arr(i, j, k,11) = uiui * u_cc_arr(i,j,k);           // (ui*ui)*u
            fab_arr(i, j, k,12) = uiui * v_cc_arr(i,j,k);           // (ui*ui)*v

            if (!use_moisture) {
                Real p = getPgivenRTh(cons_arr(i, j, k, RhoTheta_comp));
                p -= p0_arr(i,j,k);
                fab_arr(i, j, k,13) = p;                       // p
                fab_arr(i, j, k,14) = p * u_cc_arr(i,j,k);     // p*u
                fab_arr(i, j, k,15) = p * v_cc_arr(i,j,k);     // p*v
                fab_arr(i, j, k,16) = 0.;  // qv
                fab_arr(i, j, k,17) = 0.;  // qc
                fab_arr(i, j, k,18) = 0.;  // qr
                fab_arr(i, j, k,19) = 0.;  // qi
                fab_arr(i, j, k,20) = 0.;  // qs
                fab_arr(i, j, k,21) = 0.;  // qg
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
                fab_arr_stag(i,j,k,5) = pface * w_fc_arr(i,j,k);       // p*w
                fab_arr_stag(i,j,k,6) = 0.;  // w*qv
                fab_arr_stag(i,j,k,7) = 0.;  // w*qc
                fab_arr_stag(i,j,k,8) = 0.;  // w*qr
                fab_arr_stag(i,j,k,9) = 0.;  // w*thv
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
                fab_arr(i, j, k,13) = p;                       // p
                fab_arr(i, j, k,14) = p * u_cc_arr(i,j,k);     // p*u
                fab_arr(i, j, k,15) = p * v_cc_arr(i,j,k);     // p*v
                fab_arr(i, j, k,16) = cons_arr(i,j,k,RhoQ1_comp) / cons_arr(i,j,k,Rho_comp);  // qv
                fab_arr(i, j, k,17) = cons_arr(i,j,k,RhoQ2_comp) / cons_arr(i,j,k,Rho_comp);  // qc
                fab_arr(i, j, k,18) = cons_arr(i,j,k,RhoQr_comp) / cons_arr(i,j,k,Rho_comp);  // qr
                if (n_qstate > 3) {
                    fab_arr(i, j, k,19) = cons_arr(i,j,k,RhoQ3_comp) / cons_arr(i,j,k,Rho_comp);  // qi
                    fab_arr(i, j, k,20) = cons_arr(i,j,k,RhoQ5_comp) / cons_arr(i,j,k,Rho_comp);  // qs
                    fab_arr(i, j, k,21) = cons_arr(i,j,k,RhoQ6_comp) / cons_arr(i,j,k,Rho_comp);  // qg
                } else {
                    fab_arr(i, j, k,19) = 0.0;  // qi
                    fab_arr(i, j, k,20) = 0.0;  // qs
                    fab_arr(i, j, k,21) = 0.0;  // qg
                }
            });

            const Box& zbx = mfi.tilebox(IntVect(0,0,1));
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real p0 = getPgivenRTh(cons_arr(i, j, k  , RhoTheta_comp), qv_arr(i,j,k  )) - p0_arr(i,j,k  );
                Real p1 = getPgivenRTh(cons_arr(i, j, k-1, RhoTheta_comp), qv_arr(i,j,k-1)) - p0_arr(i,j,k-1);
                Real pface = 0.5 * (p0 + p1);

                Real qv0 = cons_arr(i,j,k  ,RhoQ1_comp) / cons_arr(i,j,k  ,Rho_comp);
                Real qv1 = cons_arr(i,j,k-1,RhoQ1_comp) / cons_arr(i,j,k-1,Rho_comp);
                Real qc0 = cons_arr(i,j,k  ,RhoQ2_comp) / cons_arr(i,j,k  ,Rho_comp);
                Real qc1 = cons_arr(i,j,k-1,RhoQ2_comp) / cons_arr(i,j,k-1,Rho_comp);
                Real qr0 = cons_arr(i,j,k  ,RhoQ3_comp) / cons_arr(i,j,k  ,Rho_comp);
                Real qr1 = cons_arr(i,j,k-1,RhoQ3_comp) / cons_arr(i,j,k-1,Rho_comp);
                Real qvface = 0.5 * (qv0 + qv1);
                Real qcface = 0.5 * (qc0 + qc1);
                Real qrface = 0.5 * (qr0 + qr1);

                Real theta0 = cons_arr(i,j,k  ,RhoTheta_comp) / cons_arr(i,j,k  ,Rho_comp);
                Real theta1 = cons_arr(i,j,k-1,RhoTheta_comp) / cons_arr(i,j,k-1,Rho_comp);
                Real thface = 0.5*(theta0 + theta1);
                Real ql = qcface + qrface;
                Real thv = thface * (1 + 0.61*qvface - ql);

                fab_arr_stag(i,j,k,5) = pface  * w_fc_arr(i,j,k); // p*w
                fab_arr_stag(i,j,k,6) = qvface * w_fc_arr(i,j,k); // w*qv
                fab_arr_stag(i,j,k,7) = qcface * w_fc_arr(i,j,k); // w*qc
                fab_arr_stag(i,j,k,8) = qrface * w_fc_arr(i,j,k); // w*qr
                fab_arr_stag(i,j,k,9) = thv    * w_fc_arr(i,j,k); // w*thv
            });
        } // mfi
    } // use_moisture

    // Sum in the horizontal plane
    h_avg_u = sumToLine(u_cc,0,1,     domain,zdir);
    h_avg_v = sumToLine(v_cc,0,1,     domain,zdir);
    h_avg_w = sumToLine(w_fc,0,1,stag_domain,zdir);

    h_avg_rho   = sumToLine(mf_out, 0,1,domain,zdir);
    h_avg_th    = sumToLine(mf_out, 1,1,domain,zdir);
    h_avg_ksgs  = sumToLine(mf_out, 2,1,domain,zdir);
    h_avg_Kmv   = sumToLine(mf_out, 3,1,domain,zdir);
    h_avg_Khv   = sumToLine(mf_out, 4,1,domain,zdir);
    h_avg_uu    = sumToLine(mf_out, 5,1,domain,zdir);
    h_avg_uv    = sumToLine(mf_out, 6,1,domain,zdir);
    h_avg_vv    = sumToLine(mf_out, 7,1,domain,zdir);
    h_avg_uth   = sumToLine(mf_out, 8,1,domain,zdir);
    h_avg_vth   = sumToLine(mf_out, 9,1,domain,zdir);
    h_avg_thth  = sumToLine(mf_out,10,1,domain,zdir);
    h_avg_uiuiu = sumToLine(mf_out,11,1,domain,zdir);
    h_avg_uiuiv = sumToLine(mf_out,12,1,domain,zdir);
    h_avg_p     = sumToLine(mf_out,13,1,domain,zdir);
    h_avg_pu    = sumToLine(mf_out,14,1,domain,zdir);
    h_avg_pv    = sumToLine(mf_out,15,1,domain,zdir);
    h_avg_qv    = sumToLine(mf_out,16,1,domain,zdir);
    h_avg_qc    = sumToLine(mf_out,17,1,domain,zdir);
    h_avg_qr    = sumToLine(mf_out,18,1,domain,zdir);
    h_avg_qi    = sumToLine(mf_out,19,1,domain,zdir);
    h_avg_qs    = sumToLine(mf_out,20,1,domain,zdir);
    h_avg_qg    = sumToLine(mf_out,21,1,domain,zdir);

    h_avg_uw    = sumToLine(mf_out_stag,0,1,stag_domain,zdir);
    h_avg_vw    = sumToLine(mf_out_stag,1,1,stag_domain,zdir);
    h_avg_ww    = sumToLine(mf_out_stag,2,1,stag_domain,zdir);
    h_avg_wth   = sumToLine(mf_out_stag,3,1,stag_domain,zdir);
    h_avg_uiuiw = sumToLine(mf_out_stag,4,1,stag_domain,zdir);
    h_avg_pw    = sumToLine(mf_out_stag,5,1,stag_domain,zdir);
    h_avg_wqv   = sumToLine(mf_out_stag,6,1,stag_domain,zdir);
    h_avg_wqc   = sumToLine(mf_out_stag,7,1,stag_domain,zdir);
    h_avg_wqr   = sumToLine(mf_out_stag,8,1,stag_domain,zdir);
    h_avg_wthv  = sumToLine(mf_out_stag,9,1,stag_domain,zdir);

    // Divide by the total number of cells we are averaging over
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    int unstag_size = h_avg_w.size() - 1; // _un_staggered heights
    for (int k = 0; k < unstag_size; ++k) {
        h_avg_u[k]     /= area_z;
        h_avg_v[k]     /= area_z;
        h_avg_rho[k]   /= area_z;
        h_avg_ksgs[k]  /= area_z;
        h_avg_Kmv[k]   /= area_z;
        h_avg_Khv[k]   /= area_z;
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
        h_avg_qv[k]    /= area_z;
        h_avg_qc[k]    /= area_z;
        h_avg_qr[k]    /= area_z;
        h_avg_qi[k]    /= area_z;
        h_avg_qs[k]    /= area_z;
        h_avg_qg[k]    /= area_z;
    }

    for (int k = 0; k < unstag_size+1; ++k) { // staggered heights
        h_avg_w[k]     /= area_z;
        h_avg_uw[k]    /= area_z;
        h_avg_vw[k]    /= area_z;
        h_avg_ww[k]    /= area_z;
        h_avg_wth[k]   /= area_z;
        h_avg_uiuiw[k] /= area_z;
        h_avg_pw[k]    /= area_z;
        h_avg_wqv[k]   /= area_z;
        h_avg_wqc[k]   /= area_z;
        h_avg_wqr[k]   /= area_z;
        h_avg_wthv[k]  /= area_z;
    }
}

void
ERF::derive_stress_profiles_stag (Gpu::HostVector<Real>& h_avg_tau11, Gpu::HostVector<Real>& h_avg_tau12,
                                  Gpu::HostVector<Real>& h_avg_tau13, Gpu::HostVector<Real>& h_avg_tau22,
                                  Gpu::HostVector<Real>& h_avg_tau23, Gpu::HostVector<Real>& h_avg_tau33,
                                  Gpu::HostVector<Real>& h_avg_hfx3,  Gpu::HostVector<Real>& h_avg_q1fx3,
                                  Gpu::HostVector<Real>& h_avg_q2fx3, Gpu::HostVector<Real>& h_avg_diss)
{
    int lev = 0;

    // This will hold the stress tensor components
    MultiFab mf_out(grids[lev], dmap[lev], 10, 0);

    // This will hold Tau13 and Tau23
    MultiFab mf_out_stag(convert(grids[lev], IntVect(0,0,1)), dmap[lev], 5, 0);

    bool l_use_moist   = ( solverChoice.moisture_type != MoistureType::None );

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
        const Array4<const Real>& q1fx3_arr = (l_use_moist) ? SFS_q1fx3_lev[lev]->const_array(mfi) :
                                                              Array4<const Real>{};
        const Array4<const Real>& q2fx3_arr = (l_use_moist) ? SFS_q2fx3_lev[lev]->const_array(mfi) :
                                                              Array4<const Real>{};
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
//          fab_arr(i, j, k, 6) =  hfx3_arr(i,j,k);
//          fab_arr(i, j, k, 7) =  (l_use_moist) ? q1fx3_arr(i,j,k) : 0.0;
//          fab_arr(i, j, k, 8) =  (l_use_moist) ? q2fx3_arr(i,j,k) : 0.0;
            fab_arr(i, j, k, 9) =  diss_arr(i,j,k);
        });

        const Box& zbx = mfi.tilebox(IntVect(0,0,1));
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // average from edge to face center
            fab_arr_stag(i,j,k,0) = 0.5*(tau13_arr(i,j,k) + tau13_arr(i+1,j  ,k));
            fab_arr_stag(i,j,k,1) = 0.5*(tau23_arr(i,j,k) + tau23_arr(i  ,j+1,k));

            fab_arr_stag(i,j,k,2) =  hfx3_arr(i,j,k);
            fab_arr_stag(i,j,k,3) =  (l_use_moist) ? q1fx3_arr(i,j,k) : 0.0;
            fab_arr_stag(i,j,k,4) =  (l_use_moist) ? q2fx3_arr(i,j,k) : 0.0;
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
//  h_avg_hfx3  = sumToLine(mf_out,6,1,domain,zdir);
//  h_avg_q1fx3 = sumToLine(mf_out,7,1,domain,zdir);
//  h_avg_q2fx3 = sumToLine(mf_out,8,1,domain,zdir);
    h_avg_diss  = sumToLine(mf_out,9,1,domain,zdir);

    h_avg_tau13 = sumToLine(mf_out_stag,0,1,stag_domain,zdir);
    h_avg_tau23 = sumToLine(mf_out_stag,1,1,stag_domain,zdir);
    h_avg_hfx3  = sumToLine(mf_out_stag,2,1,stag_domain,zdir);
    h_avg_q1fx3 = sumToLine(mf_out_stag,3,1,stag_domain,zdir);
    h_avg_q2fx3 = sumToLine(mf_out_stag,4,1,stag_domain,zdir);

    int ht_size =  h_avg_tau11.size(); // _un_staggered

    // Divide by the total number of cells we are averaging over
    Real area_z = static_cast<Real>(domain.length(0)*domain.length(1));
    for (int k = 0; k < ht_size; ++k) {
        h_avg_tau11[k] /= area_z; h_avg_tau12[k] /= area_z;
        h_avg_tau22[k] /= area_z;
        h_avg_tau33[k] /= area_z;
        h_avg_diss[k] /= area_z;
    }
    for (int k = 0; k < ht_size+1; ++k) { // staggered heights
        h_avg_tau13[k] /= area_z;
        h_avg_tau23[k] /= area_z;
        h_avg_hfx3[k] /= area_z;
        h_avg_q1fx3[k] /= area_z; h_avg_q2fx3[k] /= area_z;
    }
}
