#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_prob_common.H"

struct ProbParm : ProbParmDefaults {
}; // namespace ProbParm
class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* problo, const amrex::Real* probhi);

protected:
    std::string name() override { return "Metgrid"; }

private:
    ProbParm parms;
};

#endif
