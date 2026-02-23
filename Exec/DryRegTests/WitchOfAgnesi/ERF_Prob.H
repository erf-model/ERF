#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "ERF_ProbCommon.H"

struct ProbParmCore : ProbParmDefaults {};

struct ProbParm : ProbParmCore {
  std::string custom_terrain_type;
};

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitCustomPert.H"
};

#endif
