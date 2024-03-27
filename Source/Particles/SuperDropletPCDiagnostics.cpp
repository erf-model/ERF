#include <limits>
#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Compute diagnostics (max, min, avg radius, mass, etc) */
void SuperDropletPC::Diagnostics()
{
    BL_PROFILE("SuperDropletPC::Diagnostics()");
    using PTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;

    Long num_total_particles = TotalNumberOfParticles();

    auto min_par_radius = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i]; } );

    auto max_par_radius = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i]; } );

    auto avg_par_radius = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     {
                                         auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto r = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                                         return n*r;
                                     } );

    auto min_par_mass   = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i]; } );

    auto max_par_mass   = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     { return ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i]; } );

    auto avg_par_mass   = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                                     {
                                         auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                                         return n*m;
                                     } );


    ParallelDescriptor::ReduceRealMin(&min_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMin(&min_par_radius,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealMax(&max_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(&max_par_radius,1,ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::ReduceRealSum(&avg_par_mass,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealSum(&avg_par_radius,1,ParallelDescriptor::IOProcessorNumber());
    avg_par_mass /= static_cast<Real>(num_total_particles);
    avg_par_radius /= static_cast<Real>(num_total_particles);

    Print() << "SuperDropletPC(" << m_name << ") diagnostics:\n"
            << "    particle mass (min,max,avg): "
            << min_par_mass << ", " << max_par_mass << ", " << avg_par_mass << "\n"
            << "    particle radius (min,max,avg): "
            << min_par_radius << ", " << max_par_radius << ", " << avg_par_radius << "\n";

}

#endif

