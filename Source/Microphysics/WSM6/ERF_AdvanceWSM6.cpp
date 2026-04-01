#include "ERF_WSM6.H"

using namespace amrex;

void
WSM6::Advance(const Real& dt_advance,
              const SolverChoice&)
{
    dt = dt_advance;

    // Phase-1 scaffold: Fortran WSM6 bridge is added next.
    // No-op here to keep build and dispatch integration incremental.
}
