#include <iomanip>

#include "ERF.H"

using namespace amrex;

void
ERF::sum_integrated_quantities(Real time)
{
    BL_PROFILE("ERF::sum_integrated_quantities()");

    if (verbose <= 0)
      return;

    int datwidth = 14;
    int datprecision = 6;

    amrex::Real scalar = 0.0;
    amrex::Real mass   = 0.0;

    for (int lev = 0; lev <= finest_level; lev++) {
        mass   += volWgtSumMF(lev,vars_new[lev][Vars::cons],Rho_comp,true,true);
        scalar += volWgtSumMF(lev,vars_new[lev][Vars::cons],RhoScalar_comp,true,true);
    }

    if (verbose > 0) {
        const int nfoo = 2;
        amrex::Real foo[nfoo] = {mass,scalar};
#ifdef AMREX_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
        amrex::ParallelDescriptor::ReduceRealSum(
            foo, nfoo, amrex::ParallelDescriptor::IOProcessorNumber());

          if (amrex::ParallelDescriptor::IOProcessor()) {
            int i = 0;
            mass   = foo[i++];
            scalar = foo[i++];

            amrex::Print() << '\n';
            amrex::Print() << "TIME= " << time << " MASS        = " << mass   << '\n';
            amrex::Print() << "TIME= " << time << " SCALAR      = " << scalar << '\n';

            if (NumDataLogs() > 0) {
                std::ostream& data_log1 = DataLog(0);
                if (data_log1.good()) {
                    if (time == 0.0) {
                        data_log1 << std::setw(datwidth) << "          time";
                        data_log1 << std::setw(datwidth) << "          mass";
                        data_log1 << std::setw(datwidth) << "        scalar";
                        data_log1 << std::endl;
                    }

                  // Write the quantities at this time
                  data_log1 << std::setw(datwidth) << time;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << mass;
                  data_log1 << std::setw(datwidth) << std::setprecision(datprecision)
                            << scalar;
                  data_log1 << std::endl;
              }
            }
          }
#ifdef AMREX_LAZY
        });
#endif
    }
}

amrex::Real
ERF::volWgtSumMF(int lev,
  const amrex::MultiFab& mf, int comp, bool local, bool finemask)
{
  BL_PROFILE("ERF::volWgtSumMF()");

  amrex::Real sum = 0.0;
  amrex::MultiFab tmp(grids[lev], dmap[lev], 1, 0);
  amrex::MultiFab::Copy(tmp, mf, comp, 0, 1, 0);

  if (lev < finest_level && finemask) {
    const amrex::MultiFab& mask = build_fine_mask(lev+1);
    amrex::MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
  }

  amrex::MultiFab volume(grids[lev], dmap[lev], 1, 0);
  auto const& dx = geom[lev].CellSizeArray();
  Real cell_vol = dx[0]*dx[1]*dx[2];
  volume.setVal(cell_vol);
  sum = amrex::MultiFab::Dot(tmp, 0, volume, 0, 1, 0, local);

  if (!local)
    amrex::ParallelDescriptor::ReduceRealSum(sum);

  return sum;
}

amrex::MultiFab&
ERF::build_fine_mask(int level)
{
  // Mask for zeroing covered cells
  AMREX_ASSERT(level > 0);

  const amrex::BoxArray& cba = grids[level-1];
  const amrex::DistributionMapping& cdm = dmap[level-1];

  // TODO -- we should make a vector of these a member of ERF class
  fine_mask.define(cba, cdm, 1, 0, amrex::MFInfo());
  fine_mask.setVal(1.0);

  amrex::BoxArray fba = grids[level];
  amrex::iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, ref_ratio[level-1], 1, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(fine_mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto& fab = fine_mask[mfi];
        auto& ifab = ifine_mask[mfi];
        const auto arr = fab.array();
        const auto iarr = ifab.array();
        amrex::ParallelFor(
          fab.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef _OPENMP
#pragma omp atomic write
#endif
            arr(i, j, k) = iarr(i, j, k);
          });
    }

    return fine_mask;
}

bool
ERF::is_it_time_for_action(int nstep, Real time, Real dtlev, int action_interval, amrex::Real action_per)
{
  bool int_test = (action_interval > 0 && nstep % action_interval == 0);

  bool per_test = false;
  if (action_per > 0.0) {
    const int num_per_old = static_cast<int>(amrex::Math::floor((time - dtlev) / action_per));
    const int num_per_new = static_cast<int>(amrex::Math::floor((time) / action_per));

    if (num_per_old != num_per_new) {
      per_test = true;
    }
  }

  return int_test || per_test;
}
