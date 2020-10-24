#ifndef _MOL_H_
#define _MOL_H_

#include "Diffterm.H"
#include "IndexDefines.H"
#include "PeleC.H"
#include "EOS.H"
#include "Riemann.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
mol_slope(
  const int i,
  const int j,
  const int k,
  const amrex::GpuArray<const int, 3> bdim,
  const amrex::GpuArray<const int, 3> q_idx,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<amrex::Real>& dq
)
{
  bool flagArrayL = true;
  bool flagArrayR = true;

  amrex::Real dlft[QVAR] = {0.0};
  amrex::Real drgt[QVAR] = {0.0};
  const int im = i - bdim[0];
  const int jm = j - bdim[1];
  const int km = k - bdim[2];
  const int ip = i + bdim[0];
  const int jp = j + bdim[1];
  const int kp = k + bdim[2];

  const amrex::Real gamma = 1.4;
  const amrex::Real c = std::sqrt(gamma * q(i,j,k,QPRES)/q(i,j,k,QRHO));

  dlft[0] = flagArrayL ? 0.5 * (q(i, j, k, QPRES) - q(im, jm, km, QPRES)) / c -
                         0.5 *  q(i, j, k, QRHO) *
                               (q(i, j, k, q_idx[0]) - q(im, jm, km, q_idx[0]))
                       : 0.0;
  dlft[1] = flagArrayL ? 0.5 * (q(i, j, k, QPRES) - q(im, jm, km, QPRES)) / c + 
                         0.5 *  q(i, j, k, QRHO) *
                               (q(i, j, k, q_idx[0]) - q(im, jm, km, q_idx[0]))
                       : 0.0;
  dlft[2] = flagArrayL ? q(i, j, k, q_idx[1]) - q(im, jm, km, q_idx[1]) : 0.0;
  dlft[3] = flagArrayL ? q(i, j, k, q_idx[2]) - q(im, jm, km, q_idx[2]) : 0.0;

  drgt[0] = flagArrayR ? 0.5 * (q(ip, jp, kp, QPRES) - q(i, j, k, QPRES)) / c -
                           0.5 * q(i, j, k, QRHO) *
                             (q(ip, jp, kp, q_idx[0]) - q(i, j, k, q_idx[0]))
                       : 0.0;
  drgt[1] = flagArrayR ? 0.5 * (q(ip, jp, kp, QPRES) - q(i, j, k, QPRES)) / c + 
                           0.5 * q(i, j, k, QRHO) *
                             (q(ip, jp, kp, q_idx[0]) - q(i, j, k, q_idx[0]))
                       : 0.0;
  drgt[2] = flagArrayR ? q(ip, jp, kp, q_idx[1]) - q(i, j, k, q_idx[1]) : 0.0;
  drgt[3] = flagArrayR ? q(ip, jp, kp, q_idx[2]) - q(i, j, k, q_idx[2]) : 0.0;

  for (int n = 0; n < QVAR; n++) {
    const amrex::Real dcen = 0.5 * (dlft[n] + drgt[n]);
    const amrex::Real dlim =
      dlft[n] * drgt[n] >= 0.0
        ? 2.0 * amrex::min(amrex::Math::abs(dlft[n]), amrex::Math::abs(drgt[n]))
        : 0.0;
    dq(i, j, k, n) = amrex::Math::copysign(1.0, dcen) *
                     amrex::min(dlim, amrex::Math::abs(dcen));
  }
}

void pc_compute_hyp_mol_flux(
  const amrex::Box& cbox,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    a,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> del,
  const int plm_iorder
);

#endif
