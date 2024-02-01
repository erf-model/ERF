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
                                     const int& a_comp /*!< Multifav component to fill with number density */) const
{
    BL_PROFILE("SuperDropletPC::numberDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (  const SuperParticleType& p,
                                Array4<Real> const& fab_arr )
        {
            auto iv = getParticleCell(p, plo, dxi, domain);
            int num_par = p.idata(   SuperDropletsIntIdxAoS::ncomps
                                   + SuperDropletsIntIdxSoA::multiplicity );
            Real num_dens = num_par * inv_cell_volume;
            Gpu::Atomic::AddNoRet( &fab_arr(iv, a_comp), num_dens );
        }, false);

    return;
}

/*! Computes the mass density of the particles over a mesh: this does not
    include the aerosol mass */
void SuperDropletPC::massDensity ( MultiFab&  a_mf,  /*!< Mass density multifab */
                                   const int& a_comp /*!< Multifav component to fill with mass density */) const
{
    BL_PROFILE("SuperDropletPC::numberDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const Real mat_density = 1.0; /* TODO: material object */

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (  const SuperParticleType& p,
                                Array4<Real> const& fab_arr )
        {
            auto iv = getParticleCell(p, plo, dxi, domain);
            int num_par = p.idata(   SuperDropletsIntIdxAoS::ncomps
                                   + SuperDropletsIntIdxSoA::multiplicity );
            Real radius = p.rdata(   SuperDropletsRealIdxAoS::ncomps
                                   + SuperDropletsRealIdxSoA::radius );
            Real par_volume = (4.0/3.0)*PI*(radius*radius*radius);
            Real mass = mat_density*par_volume*num_par;
            Real mass_dens = mass * inv_cell_volume;
            Gpu::Atomic::AddNoRet( &fab_arr(iv, a_comp), mass_dens );
        }, false);

    return;
}

#endif


