#include <chrono>
#include <ctime>

#include "IOManager.H"

#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

IOManager::IOManager(ERF& erf)
    : erf(erf) {}

IOManager::~IOManager() = default;

