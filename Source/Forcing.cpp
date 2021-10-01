#include "ERF.H"
#include "Forcing.H"

forcing_params fparms;

void
ERF::construct_old_forcing_source(amrex::Real time, amrex::Real dt)
{
  amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[forcing_src]->setVal(0.0);

  if (!add_forcing_src)
    return;

  fill_forcing_source(time, dt, S_old, S_old, *old_sources[forcing_src], ng);

  old_sources[forcing_src]->FillBoundary(geom.periodicity());
}

void
ERF::construct_new_forcing_source(amrex::Real time, amrex::Real dt)
{
  amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[forcing_src]->setVal(0.0);

  if (!add_forcing_src)
    return;

  fill_forcing_source(time, dt, S_old, S_new, *new_sources[forcing_src], ng);
}

void
ERF::fill_forcing_source(
  amrex::Real time,
  amrex::Real dt,
  const amrex::MultiFab& state_old,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& forcing_src,
  int ng)
{
//  const amrex::Real* dx = geom.CellSize();
//  const amrex::Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(forcing_src, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);
    amrex::RealBox gridloc =
      amrex::RealBox(grids[mfi.index()], geom.CellSize(), geom.ProbLo());

//    auto const& sarr = state_new.array(mfi);
//    auto const& src = forcing_src.array(mfi);

    // Evaluate the linear forcing term
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    //amrex::ParallelFor(bx, [=,
    //                        forcing=fparms.forcing,
    //                        u0 = fparms.u0,
    //                        v0 = fparms.v0,
    //                        w0 = fparms.w0] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
//    src(i, j, k, UMX) = forcing * sarr(i, j, k, Rho_comp) *
//                        (sarr(i, j, k, UMX) - u0);
//    src(i, j, k, UMY) = forcing * sarr(i, j, k, Rho_comp) *
//                        (sarr(i, j, k, UMY) - v0);
//    src(i, j, k, UMZ) = forcing * sarr(i, j, k, Rho_comp) *
//                        (sarr(i, j, k, UMZ) - w0);
    });
  }
}
