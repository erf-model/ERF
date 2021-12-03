#include <iomanip>

#include "ERF.H"

void
ERF::sum_integrated_quantities()
{
  BL_PROFILE("ERF::sum_integrated_quantities()");

  if (verbose <= 0)
    return;

  bool local_flag = true;

  int finest_level = parent->finestLevel();
  amrex::Real time = state[State_Type].curTime();

  amrex::Real scalar = 0.0;
  amrex::Real mass   = 0.0;

  int datwidth = 14;
  int datprecision = 6;

  for (int lev = 0; lev <= finest_level; lev++) {
    ERF& erf_lev = getLevel(lev);
    mass   += erf_lev.volWgtSum("density", time, local_flag);
    scalar += erf_lev.volWgtSum("adv_0"  , time, local_flag);
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

        if (parent->NumDataLogs() > 0) {
          std::ostream& data_log1 = parent->DataLog(0);
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
