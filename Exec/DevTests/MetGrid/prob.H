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
    Problem();

protected:
    std::string name() override { return "Metgrid"; }

private:
    ProbParm parms;
};

#endif
