#include <AMReX_ParticleInterpolators.H>
#include <ERF_Constants.H>
#include "SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! This returns the total number of particles that all the super-droplets represent */
Long SuperDropletPC::TotalNumberOfParticles ()
{
    BL_PROFILE("SuperDropletPC::TotalNumberOfParticles()");

    Long count = 0;
    for (ParIterType pti(*this, m_lev); pti.isValid(); ++pti) {
        const auto& particle_tile = ParticlesAt(m_lev, pti);
        auto& soa = particle_tile.GetStructOfArrays();
        auto& aos = particle_tile.GetArrayOfStructs();

        auto* mult_ptr = soa.GetIntData(SuperDropletsIntIdxSoA::multiplicity).data();
        const int n = aos.numParticles();

        count += Reduce::Sum(n, mult_ptr);
    }

    ParallelDescriptor::ReduceLongSum(&count, 1, ParallelDescriptor::IOProcessorNumber());
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
                    Real num_par = (Real) ptd.m_runtime_idata[SuperDropletsIntIdxSoA::multiplicity][i];
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
                    int num_par = ptd.m_runtime_idata[SuperDropletsIntIdxSoA::multiplicity][i];
                    Real par_mass = p.rdata(SuperDropletsRealIdxAoS::mass);
                    return num_par*par_mass*inv_cell_volume;
                });
        });

    return;
}

#endif


