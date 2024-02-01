#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "prob_common.H"

struct ProbParm : ProbParmDefaults {
}; // namespace ProbParm
class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/);

protected:
    std::string name() override { return "LandSurfaceModel"; }

private:
    ProbParm parms;
};

#endif
