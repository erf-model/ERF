#include <cstdio>

#include <AMReX_LevelBld.H>
#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>

#include "ERF.H"
#include "Derive.H"
#include "IndexDefines.H"
#include "prob.H"

void
ERF::clear_prob()
{
  erf_prob_close();
}

void
ERF::variableCleanUp()
{
  desc_lst.clear();

  clear_prob();
}

