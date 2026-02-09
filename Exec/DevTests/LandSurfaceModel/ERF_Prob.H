#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
}; // namespace ProbParm
class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/);

protected:

private:
    ProbParm parms;
};

#endif
