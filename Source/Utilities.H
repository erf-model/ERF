#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <AMReX_FArrayBox.H>
#include "Constants.H"
#include "IndexDefines.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
erf_move_transcoefs_to_ec(
  const int i,
  const int j,
  const int k,
  const int n,
  const amrex::Array4<const amrex::Real>& carr,
  amrex::Real* earr,
  const int dir,
  const int do_harmonic)

{
  if (do_harmonic == 0) {
    if (dir == 0)
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i - 1, j, k, n));
    else if (dir == 1)
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j - 1, k, n));
    else
      earr[n] = 0.5 * (carr(i, j, k, n) + carr(i, j, k - 1, n));
  } else {
    amrex::Real a = carr(i, j, k, n);
    if (dir == 0) {
      amrex::Real b = carr(i - 1, j, k, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    } else if (dir == 1) {
      amrex::Real b = carr(i, j - 1, k, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    } else {
      amrex::Real b = carr(i, j, k - 1, n);
      if (a * b > 0.0)
        earr[n] = 2.0 * (a * b) / (a + b);
      else
        earr[n] = 0.0;
    }
  }
}

AMREX_FORCE_INLINE
void
setV(
  const amrex::Box box,
  const int Ncomp,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      array(i, j, k, n) = val;
    });
}

AMREX_FORCE_INLINE
void
setC(
  const amrex::Box box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<amrex::Real>& array,
  const amrex::Real val)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++)
      array(i, j, k, n) = val;
  });
}

AMREX_FORCE_INLINE
void
copy_array4(
  const amrex::Box box,
  const int Ncomp,
  const amrex::Array4<const amrex::Real>& in,
  const amrex::Array4<amrex::Real>& out)
{
  amrex::ParallelFor(
    box, Ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      out(i, j, k, n) = in(i, j, k, n);
    });
}

AMREX_FORCE_INLINE
void
lincomb_array4(
  const amrex::Box box,
  const int Nstart,
  const int Nend,
  const amrex::Array4<const amrex::Real>& f1,
  const amrex::Array4<const amrex::Real>& f2,
  const amrex::Real alpha,
  const amrex::Real beta,
  const amrex::Array4<amrex::Real>& d)
{
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    for (int n = Nstart; n < Nend; n++)
      d(i, j, k, n) = alpha * f1(i, j, k, n) + beta * f2(i, j, k, n);
  });
}

#endif
