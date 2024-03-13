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
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
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
                                   const int& a_comp /*!< Multifab component to fill with mass density */) const
{
    BL_PROFILE("SuperDropletPC::massDensity()");

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
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                    return num_par*par_mass*inv_cell_volume;
                });
        });

    return;
}

/*! Computes the mass density of the condensate over a mesh: this does
    include the aerosol mass*/
void SuperDropletPC::massDensityCondensate ( MultiFab&  a_mf,  /*!< Mass density multifab */
                                             const int& a_comp /*!< Multifab component to fill with mass density */) const
{
    BL_PROFILE("SuperDropletPC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(0.0);
    const int num_aerosols = m_num_aerosols;

    ParticleToMesh( *this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (  const SuperDropletPC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                    ParticleReal solute_mass = 0.0;
                    for (int ai = 0; ai < num_aerosols; ai++) {
                        solute_mass += (ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::ncomps+ai][i]);
                    }
                    return num_par*(par_mass-solute_mass)*inv_cell_volume;
                });
        });

    return;
}

/*! Computes the particle velocity components over a mesh */
void SuperDropletPC::massFlux ( MultiFab&  a_mf,  /*!< Mass flux multifab */
                                const int& a_dim, /*!< Flux component */
                                const int& a_comp /*!< Multifab component to fill with mass density */) const
{
    BL_PROFILE("SuperDropletPC::velocityComp()");

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
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                    auto par_velocity = ptd.m_rdata[SuperDropletsRealIdxSoA::vx+a_dim][i];
                    return num_par * par_mass * par_velocity * inv_cell_volume;
                });
        });

    return;
}

/*! Computes the particle velocity components over a mesh */
void SuperDropletPC::massFluxCondensate ( MultiFab&  a_mf,  /*!< Condensate Mass flux multifab */
                                          const int& a_dim, /*!< Flux component */
                                          const int& a_comp /*!< Multifab component to fill with mass density */) const
{
    BL_PROFILE("SuperDropletPC::velocityComp()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(0.0);
    const int num_aerosols = m_num_aerosols;

    ParticleToMesh( *this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (  const SuperDropletPC::ParticleTileType::ConstParticleTileDataType& ptd,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                    ParticleReal solute_mass = 0.0;
                    for (int ai = 0; ai < num_aerosols; ai++) {
                        solute_mass += (ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::ncomps+ai][i]);
                    }
                    auto par_velocity = ptd.m_rdata[SuperDropletsRealIdxSoA::vx+a_dim][i];
                    return num_par * (par_mass-solute_mass) * par_velocity * inv_cell_volume;
                });
        });

    return;
}

/*! Computes the number density of the particles over a mesh */
void SuperDropletPC::aerosolMassDensity ( MultiFab&  a_mf,  /*!< Aerosol mass density multifab */
                                          const int& a_idx, /*!< Aerosol index */
                                          const int& a_comp /*!< Multifab component to fill with number density */) const
{
    BL_PROFILE("SuperDropletPC::aerosolMassDensity()");

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
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto aero_mass = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::ncomps+a_idx][i];
                    return num_par*aero_mass*inv_cell_volume;
                });
        });

    return;
}

/*! Computes the number density of the particles over a mesh */
void SuperDropletPC::aerosolMassFlux ( MultiFab&  a_mf,  /*!< Aerosol mass density multifab */
                                       const int& a_idx, /*!< Aerosol index */
                                       const int& a_dim, /*!< Flux component */
                                       const int& a_comp /*!< Multifab component to fill with number density */) const
{
    BL_PROFILE("SuperDropletPC::aerosolMassDensity()");

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
                [=] AMREX_GPU_DEVICE ( const SuperDropletPC::ParticleType&, int)
                {
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto aero_mass = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::ncomps+a_idx][i];
                    auto par_velocity = ptd.m_rdata[SuperDropletsRealIdxSoA::vx+a_dim][i];
                    return num_par*aero_mass*par_velocity*inv_cell_volume;
                });
        });

    return;
}

#endif


