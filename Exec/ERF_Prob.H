#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "ERF_ProbCommon.H"
#include "Prob/ERF_CloudChamber.H"

class Problem : public ProblemBase
{
public:
    Problem (const amrex::Real* problo, const amrex::Real* probhi);

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitCustomPert.H"
#include "Prob/ERF_InitRayleighDamping.H"
#include "Prob/ERF_UpdateSources.H"

private:
    erf_cloud_chamber::Config m_cloud_chamber_config;
};

#endif
