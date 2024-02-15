#include <AMReX_ParticleInterpolators.H>
#include <ERF_Constants.H>
#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! This returns the total number of particles that all the super-droplets represent */
Real SuperDropletPC::TotalNumberOfParticles ()
{
    BL_PROFILE("SuperDropletPC::TotalNumberOfParticles()");
    using PTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;
    Real count = ReduceSum(*this,
                           [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real
                           { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; });
    ParallelDescriptor::ReduceRealSum(&count, 1, ParallelDescriptor::IOProcessorNumber());
    return count;
}

/*! Computes the number density of the particles over a mesh */
void SuperDropletPC::numberDensity ( MultiFab&  a_mf,  /*!< Number density multifab */
                                     const int& a_comp /*!< Multifab component to fill with number density */) const
{
    BL_PROFILE("SuperDropletPC::numberDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    a_mf.setVal(0.0);

    ParticleToMesh( *this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (  const SuperDropletPC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType& part, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    return num_par*inv_cell_volume;
                });
        });

    return;
}

/*! Computes the mass density of the particles over a mesh: this does
    include the aerosol mass*/
void SuperDropletPC::massDensity ( MultiFab&  a_mf,  /*!< Mass density multifab */
                                   const int& a_comp /*!< Multifav component to fill with mass density */) const
{
    BL_PROFILE("SuperDropletPC::numberDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(0.0);

    ParticleToMesh( *this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (  const SuperDropletPC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType& part, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                    return num_par*par_mass*inv_cell_volume;
                });
        });

    return;
}

#endif


