#include <AMReX_ParticleInterpolators.H>
#include <ERF_Constants.H>
#include <ERFPC.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;

void ERFPC::massDensity ( MultiFab&  a_mf,
                          const int& a_lev,
                          const int& a_comp ) const
{
    BL_PROFILE("ERFPC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const Real inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    a_mf.setVal(0.0);

    ParticleToMesh( *this, a_mf, a_lev,
        [=] AMREX_GPU_DEVICE (  const ERFPC::ParticleTileType::ConstParticleTileDataType& ptr,
                                int i, Array4<Real> const& rho)
        {
            auto p = ptr.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh ( p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE ( const ERFPC::ParticleType&, int)
                {
                    auto mass = ptr.m_rdata[ERFParticlesRealIdxSoA::mass][i];
                    return mass*inv_cell_volume;
                });
        });

    return;
}

#endif
