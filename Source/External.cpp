#include "ERF.H"
#include "IndexDefines.H"

void
ERF::construct_old_ext_source(amrex::Real time, amrex::Real dt)
{
  amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[ext_src]->setVal(0.0);

  if (!add_ext_src)
    return;

  fill_ext_source(time, dt, S_old, S_old, *old_sources[ext_src], ng);

  old_sources[ext_src]->FillBoundary(geom.periodicity());
}

void
ERF::construct_new_ext_source(amrex::Real time, amrex::Real dt)
{
  amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[ext_src]->setVal(0.0);

  if (!add_ext_src)
    return;

  fill_ext_source(time, dt, S_old, S_new, *new_sources[ext_src], ng);
}

void
ERF::fill_ext_source(
  amrex::Real time,
  amrex::Real dt,
  const amrex::MultiFab& state_old,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& ext_src,
  int ng)
{
//  const amrex::Real* dx = geom.CellSize();
//  const amrex::Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(ext_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

//    auto const& So = state_old.array(mfi);
//    auto const& Sn = state_new.array(mfi);
    auto const& Farr = ext_src.array(mfi);

    // Evaluate the external source
    amrex::ParallelFor(
      bx, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        Farr(i, j, k, n) = 0.0;
      });
  }
}
