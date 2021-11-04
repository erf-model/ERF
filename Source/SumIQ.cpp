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

  amrex::Real mass      = 0.0;
  int datwidth = 14;
  int datprecision = 6;

  // Compute sum on all levels -- this uses the data at
  //    the finest level available
  for (int lev = 0; lev <= finest_level; lev++) {
      ERF& erf_lev = getLevel(lev);
      mass += erf_lev.volWgtSum("density", time, local_flag);
  }

  // Compute sum on coarsest level only -- this definition of
  //    mass is conserved with OneWay coupling.
  amrex::Real mass_crse;
  {
      ERF& erf_lev = getLevel(0);
      mass_crse = erf_lev.volWgtSum("density", time, local_flag, false);
  }

  if (verbose > 0) {
    const int nfoo = 2;
    amrex::Real foo[nfoo] = {mass,mass_crse};
#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealSum(
        foo, nfoo, amrex::ParallelDescriptor::IOProcessorNumber());

      if (amrex::ParallelDescriptor::IOProcessor()) {
        mass      = foo[0];
        mass_crse = foo[1];

        amrex::Print() << '\n';
        amrex::Print() << "TIME= " << time << " MASS (all levels / level 0) = " << mass << " " << mass_crse << '\n';

        if (parent->NumDataLogs() > 0) {
          std::ostream& data_log1 = parent->DataLog(0);
          if (data_log1.good()) {
            if (time == 0.0) {
              data_log1 << std::setw(datwidth) << "          time";
              data_log1 << std::setw(datwidth) << "          mass";
              data_log1 << std::setw(datwidth) << "     mass_crse";
            }

            // Write the quantities at this time
            data_log1 << std::setw(datwidth) << time;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision) << mass;
            data_log1 << std::setw(datwidth) << std::setprecision(datprecision) << mass_crse;
            data_log1 << std::endl;
          }
        }
      }
#ifdef AMREX_LAZY
    });
#endif
  }
}
