#include <iomanip>

#include "ERF.H"

amrex::Real
ERF::sumDerive(const std::string& name, amrex::Real time, bool local)
{
  amrex::Real sum = 0.0;
  auto mf = derive(name, time, 0);

  AMREX_ASSERT(!(mf == 0));

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion()) reduction(+ : sum)
#endif
  {
    for (amrex::MFIter mfi(*mf, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      sum += (*mf)[mfi].sum<amrex::RunOn::Device>(mfi.tilebox(), 0);
    }
  }

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

amrex::Real
ERF::volWgtSum(
  const std::string& name, amrex::Real time, bool local, bool finemask)
{
  BL_PROFILE("ERF::volWgtSum()");

  amrex::Real sum = 0.0;
  auto mf = derive(name, time, 0);

  AMREX_ASSERT(mf != 0);

  if (level < parent->finestLevel() && finemask) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

  sum = amrex::MultiFab::Dot(*mf, 0, volume, 0, 1, 0, local);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

amrex::Real
ERF::volWgtSquaredSum(const std::string& name, amrex::Real time, bool local)
{
  BL_PROFILE("ERF::volWgtSquaredSum()");

  amrex::Real sum = 0.0;
  auto mf = derive(name, time, 0);

  AMREX_ASSERT(mf != 0);

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

  amrex::MultiFab::Multiply(*mf, *mf, 0, 0, 1, 0);

  amrex::MultiFab vol(grids, dmap, 1, 0);
  amrex::MultiFab::Copy(vol, volume, 0, 0, 1, 0);
  sum = amrex::MultiFab::Dot(*mf, 0, vol, 0, 1, 0, local);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

amrex::Real
ERF ::volWgtSquaredSumDiff(int comp, amrex::Real /*time*/, bool local)
{

  // Calculate volume weighted sum of the square of the difference
  // between the old and new quantity

  amrex::Real sum = 0.0;
  amrex::MultiFab& S_old = get_old_data(State_Type);
  amrex::MultiFab& S_new = get_new_data(State_Type);
  amrex::MultiFab diff(grids, dmap, 1, 0);

  // Calculate the difference between the states
  amrex::MultiFab::Copy(diff, S_old, comp, 0, 1, 0);
  amrex::MultiFab::Subtract(diff, S_new, comp, 0, 1, 0);

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(diff, mask, 0, 0, 1, 0);
  }

  amrex::MultiFab::Multiply(diff, diff, 0, 0, 1, 0);

  amrex::MultiFab vol(grids, dmap, 1, 0);
  amrex::MultiFab::Copy(vol, volume, 0, 0, 1, 0);
  sum = amrex::MultiFab::Dot(diff, 0, vol, 0, 1, 0, local);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

amrex::Real
ERF::volWgtSumMF(
  const amrex::MultiFab& mf, int comp, bool local, bool finemask)
{
  BL_PROFILE("ERF::volWgtSumMF()");

  amrex::Real sum = 0.0;
  amrex::MultiFab vol(grids, dmap, 1, 0);
  amrex::MultiFab::Copy(vol, mf, comp, 0, 1, 0);

  if (level < parent->finestLevel() && finemask) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(vol, mask, 0, 0, 1, 0);
  }

  sum = amrex::MultiFab::Dot(vol, 0, volume, 0, 1, 0, local);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

amrex::Real
ERF::maxDerive(const std::string& name, amrex::Real time, bool local)
{
  auto mf = derive(name, time, 0);

  BL_ASSERT(!(mf == 0));

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

  return mf->max(0, 0, local);
}

bool
ERF::is_it_time_for_action(int action_interval, amrex::Real action_per)
{
  int nstep = parent->levelSteps(0);
  amrex::Real dtlev = parent->dtLevel(0);
  amrex::Real cumtime = parent->cumTime();
  if (cumtime != 0) {
    cumtime += dtlev;
  }

  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > 0.0) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((cumtime - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((cumtime) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}
