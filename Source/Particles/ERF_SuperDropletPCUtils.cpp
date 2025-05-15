#include <AMReX_ParticleInterpolators.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;

/*! Compute mesh variable from particles */
void SuperDropletPC::computeMeshVar( const std::string&  a_var_name,
                                     MultiFab&           a_mf,
                                     const int           a_lev) const
{
    a_mf.setVal(0.0);
    if (a_lev == 0) {
        if (a_var_name == "number_density") {
            numberDensity( a_mf );
        } else if (a_var_name == "mass_density") {
            massDensity( a_mf );
        } else if (a_var_name == ("mass_density_"+getEnumNameString(m_vapour_mat->m_name))) {
            massDensityCondensate( a_mf );
        } else if (a_var_name == ("mass_flux_x_"+getEnumNameString(m_vapour_mat->m_name))) {
            massFluxCondensate( a_mf, 0 );
        } else if (a_var_name == ("mass_flux_y_"+getEnumNameString(m_vapour_mat->m_name))) {
            massFluxCondensate( a_mf, 1 );
        } else if (a_var_name == ("mass_flux_z_"+getEnumNameString(m_vapour_mat->m_name))) {
            massFluxCondensate( a_mf, 2 );
        } else if (a_var_name == "mass_flux_x") {
            massFlux( a_mf, 0 );
        } else if (a_var_name == "mass_flux_y") {
            massFlux( a_mf, 1 );
        } else if (a_var_name == "mass_flux_z") {
            massFlux( a_mf, 2 );
        } else if (a_var_name == "radius") {
            effectiveRadius( a_mf );
        } else {
            for (int i = 0; i < m_num_aerosols; i++) {
                std::string var_name = "aerosol_mass_density_"+getEnumNameString(m_aerosol_mat[i]->m_name);
                if (a_var_name == var_name) {
                    aerosolMassDensity( a_mf, i );
                    break;
                }
            }
            for (int i = 0; i < m_num_aerosols; i++) {
                std::string var_name = "aerosol_mass_flux_x_"+getEnumNameString(m_aerosol_mat[i]->m_name);
                if (a_var_name == var_name) {
                    aerosolMassFlux( a_mf, i, 0 );
                    break;
                }
            }
            for (int i = 0; i < m_num_aerosols; i++) {
                std::string var_name = "aerosol_mass_flux_y_"+getEnumNameString(m_aerosol_mat[i]->m_name);
                if (a_var_name == var_name) {
                    aerosolMassFlux( a_mf, i, 1 );
                    break;
                }
            }
            for (int i = 0; i < m_num_aerosols; i++) {
                std::string var_name = "aerosol_mass_flux_z_"+getEnumNameString(m_aerosol_mat[i]->m_name);
                if (a_var_name == var_name) {
                    aerosolMassFlux( a_mf, i, 2 );
                    break;
                }
            }
        }
    }
}

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

/*! This returns the total number of deactivated superdroplets (multiplicity zero) */
Long SuperDropletPC::NumSDDeactivated ()
{
    BL_PROFILE("SuperDropletPC::NumSDDeactivated()");
    using PTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;
    auto count = ReduceSum( *this,
                            [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Long
                            {
                                auto mult = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                if (mult == 0) { return Long(1); }
                                else           { return Long(0); }
                            } );
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

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

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

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
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
void SuperDropletPC::massDensityCondensate ( MultiFab&  a_mf, /*!< Mass density multifab */
                                             const Real a_rmin, /*!< minimum radius */
                                             const Real a_rmax, /*!< maximum radius */
                                             const int& a_comp /*!< Multifab component to fill with mass density */ ) const
{
    BL_PROFILE("SuperDropletPC::massDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
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
                    auto radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                    if ((radius < a_rmin) || (radius >= a_rmax)) {
                        return 0.0;
                    } else {
                        auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                        auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                        ParticleReal solute_mass = 0.0;
                        for (int ai = 0; ai < num_aerosols; ai++) {
                            solute_mass += (ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::ncomps+ai][i]);
                        }
                        ParticleReal condensate_mass = std::max(0.0, par_mass-solute_mass);
                        return num_par*condensate_mass*inv_cell_volume;
                    }
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

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
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
                    if (a_dim == 2) {
                        auto term_vel = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                        par_velocity -= term_vel;
                    }
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

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
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
                    if (a_dim == 2) {
                        auto term_vel = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                        par_velocity -= term_vel;
                    }
                    ParticleReal condensate_mass = std::max(0.0, par_mass-solute_mass);
                    return num_par * condensate_mass * par_velocity * inv_cell_volume;
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

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

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

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

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
                    if (a_dim == 2) {
                        auto term_vel = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                        par_velocity -= term_vel;
                    }
                    return num_par*aero_mass*par_velocity*inv_cell_volume;
                });
        });

    return;
}

/*! Computes the effective radius of the particles over a mesh */
void SuperDropletPC::effectiveRadius (  MultiFab&  a_mf,  /*!< Effective radius multifab */
                                        const int& a_comp /*!< Multifab component to fill with number density */) const
{
    BL_PROFILE("SuperDropletPC::effectiveRadius()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    a_mf.setVal(0.0);

    MultiFab number_density( a_mf.boxArray(), a_mf.DistributionMap(), 1, a_mf.nGrowVect() );
    numberDensity(number_density);

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
                    auto radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                    return num_par*radius*inv_cell_volume;
                });
        });

    for ( MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto mf_arr = a_mf.array(mfi);
        const auto nd_arr = number_density.const_array(mfi);
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                          {
                              if (nd_arr(i,j,k,0) > 0) {
                                  mf_arr(i,j,k,a_comp) /= nd_arr(i,j,k,0);
                              } else {
                                  mf_arr(i,j,k,a_comp) = 0.0;
                              }
                          } );
    }


    return;
}

#endif


