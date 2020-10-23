#include "MOL.H"

void
pc_compute_hyp_mol_flux(
  const amrex::Box& cbox,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const int plm_iorder
)
{
  const int R_RHO = 0;
  const int R_UN = 1;
  const int R_UT1 = 2;
  const int R_UT2 = 3;
  const int R_P = 4;
  const int R_Y = 5;
  const int bc_test_val = 1;

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::FArrayBox dq_fab(cbox, QVAR);
    amrex::Elixir dq_fab_eli = dq_fab.elixir();
    auto const& dq = dq_fab.array();
    setV(cbox, QVAR, dq, 0.0);

    // dimensional indexing
    const amrex::GpuArray<const int, 3> bdim{dir == 0, dir == 1, dir == 2};
    const amrex::GpuArray<const int, 3> q_idx{
      bdim[0] * QU + bdim[1] * QV + bdim[2] * QW,
      bdim[0] * QV + bdim[1] * QU + bdim[2] * QU,
      bdim[0] * QW + bdim[1] * QW + bdim[2] * QV};
    const amrex::GpuArray<const int, 3> f_idx{
      bdim[0] * UMX + bdim[1] * UMY + bdim[2] * UMZ,
      bdim[0] * UMY + bdim[1] * UMX + bdim[2] * UMX,
      bdim[0] * UMZ + bdim[1] * UMZ + bdim[2] * UMY};

    if (plm_iorder != 1) {
      amrex::ParallelFor(
        cbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          mol_slope(
            i, j, k, bdim, q_idx, q, qaux, dq
          );
        });
    }
    const amrex::Box tbox = amrex::grow(cbox, dir, -1);
    const amrex::Box ebox = amrex::surroundingNodes(tbox, dir);
    amrex::ParallelFor(
      ebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const int ii = i - bdim[0];
        const int jj = j - bdim[1];
        const int kk = k - bdim[2];

        amrex::Real qtempl[5 + NUM_SPECIES] = {0.0};
        qtempl[R_UN] =
          q(ii, jj, kk, q_idx[0]) +
          0.5 * ((dq(ii, jj, kk, 1) - dq(ii, jj, kk, 0)) / q(ii, jj, kk, QRHO));
        qtempl[R_P] =
          q(ii, jj, kk, QPRES) +
          0.5 * (dq(ii, jj, kk, 0) + dq(ii, jj, kk, 1)) * qaux(ii, jj, kk, QC);
        qtempl[R_UT1] = q(ii, jj, kk, q_idx[1]) + 0.5 * dq(ii, jj, kk, 2);
        qtempl[R_UT2] = q(ii, jj, kk, q_idx[2]) + 0.5 * dq(ii, jj, kk, 3);
        qtempl[R_RHO] = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempl[R_Y + n] = q(ii, jj, kk, QFS + n) * q(ii, jj, kk, QRHO) +
                            0.5 * (dq(ii, jj, kk, 4 + n) +
                                   q(ii, jj, kk, QFS + n) *
                                     (dq(ii, jj, kk, 0) + dq(ii, jj, kk, 1)) /
                                     qaux(ii, jj, kk, QC));
          qtempl[R_RHO] += qtempl[R_Y + n];
        }

        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempl[R_Y + n] = qtempl[R_Y + n] / qtempl[R_RHO];
        }

        amrex::Real qtempr[5 + NUM_SPECIES] = {0.0};
        qtempr[R_UN] =
          q(i, j, k, q_idx[0]) -
          0.5 * ((dq(i, j, k, 1) - dq(i, j, k, 0)) / q(i, j, k, QRHO));
        qtempr[R_P] = q(i, j, k, QPRES) - 0.5 *
                                            (dq(i, j, k, 0) + dq(i, j, k, 1)) *
                                            qaux(i, j, k, QC);
        qtempr[R_UT1] = q(i, j, k, q_idx[1]) - 0.5 * dq(i, j, k, 2);
        qtempr[R_UT2] = q(i, j, k, q_idx[2]) - 0.5 * dq(i, j, k, 3);
        qtempr[R_RHO] = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempr[R_Y + n] =
            q(i, j, k, QFS + n) * q(i, j, k, QRHO) -
            0.5 * (dq(i, j, k, 4 + n) + q(i, j, k, QFS + n) *
                                          (dq(i, j, k, 0) + dq(i, j, k, 1)) /
                                          qaux(i, j, k, QC));
          qtempr[R_RHO] += qtempr[R_Y + n];
        }
        for (int n = 0; n < NUM_SPECIES; n++) {
          qtempr[R_Y + n] = qtempr[R_Y + n] / qtempr[R_RHO];
        }

        const amrex::Real cavg =
          0.5 * (qaux(i, j, k, QC) + qaux(ii, jj, kk, QC));
        const amrex::Real csmall =
          amrex::min(qaux(i, j, k, QCSML), qaux(ii, jj, kk, QCSML));

        amrex::Real eos_state_rho, eos_state_p, eos_state_e, eos_state_cs,
          eos_state_gamma, eos_state_T;

        eos_state_rho = qtempl[R_RHO];
        eos_state_p = qtempl[R_P];
        amrex::Real spl[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
          spl[n] = qtempl[R_Y + n];
        }
        EOS::RYP2T(eos_state_rho, spl, eos_state_p, eos_state_T);
        EOS::RYP2E(eos_state_rho, spl, eos_state_p, eos_state_e);
        EOS::TY2G(eos_state_T, spl, eos_state_gamma);
        EOS::RPY2Cs(eos_state_rho, eos_state_p, spl, eos_state_cs);
        const amrex::Real rhoe_l = eos_state_rho * eos_state_e;
        const amrex::Real gamc_l = eos_state_gamma;

        eos_state_rho = qtempr[R_RHO];
        eos_state_p = qtempr[R_P];
        amrex::Real spr[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
          spr[n] = qtempr[R_Y + n];
        }
        EOS::RYP2T(eos_state_rho, spr, eos_state_p, eos_state_T);
        EOS::RYP2E(eos_state_rho, spr, eos_state_p, eos_state_e);
        EOS::TY2G(eos_state_T, spr, eos_state_gamma);
        EOS::RPY2Cs(eos_state_rho, eos_state_p, spr, eos_state_cs);
        const amrex::Real rhoe_r = eos_state_rho * eos_state_e;
        const amrex::Real gamc_r = eos_state_gamma;

        amrex::Real flux_tmp[NVAR] = {0.0};
        amrex::Real ustar = 0.0;

        amrex::Real tmp0, tmp1, tmp2, tmp3, tmp4;
        riemann(
          qtempl[R_RHO], qtempl[R_UN], qtempl[R_UT1], qtempl[R_UT2],
          qtempl[R_P], rhoe_l, spl, gamc_l, qtempr[R_RHO], qtempr[R_UN],
          qtempr[R_UT1], qtempr[R_UT2], qtempr[R_P], rhoe_r, spr, gamc_r,
          bc_test_val, csmall, cavg, ustar, flux_tmp[URHO], flux_tmp[f_idx[0]],
          flux_tmp[f_idx[1]], flux_tmp[f_idx[2]], flux_tmp[UEDEN],
          flux_tmp[UEINT], tmp0, tmp1, tmp2, tmp3, tmp4);

        for (int n = 0; n < NUM_SPECIES; n++) {
          flux_tmp[UFS + n] = (ustar > 0.0) ? flux_tmp[URHO] * qtempl[R_Y + n]
                                            : flux_tmp[URHO] * qtempr[R_Y + n];
          flux_tmp[UFS + n] =
            (ustar == 0.0)
              ? flux_tmp[URHO] * 0.5 * (qtempl[R_Y + n] + qtempr[R_Y + n])
              : flux_tmp[UFS + n];
        }

        flux_tmp[UTEMP] = 0.0;
        for (int n = UFX; n < UFX + NUM_AUX; n++) {
          flux_tmp[n] = (NUM_AUX > 0) ? 0.0 : flux_tmp[n];
        }
        for (int n = UFA; n < UFA + NUM_ADV; n++) {
          flux_tmp[n] = (NUM_ADV > 0) ? 0.0 : flux_tmp[n];
        }

        for (int ivar = 0; ivar < NVAR; ivar++) {
          flx[dir](i, j, k, ivar) += flux_tmp[ivar] * a[dir](i, j, k);
        }
      });
  }
}
