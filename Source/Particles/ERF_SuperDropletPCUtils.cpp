#include <AMReX_ParticleInterpolators.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;
using SDTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;

/*! Initialize device property arrays for species and aerosol materials */
void SuperDropletPC::initializeDeviceProperties()
{
    if (m_device_props_initialized) return;

    const int num_sp = m_species_mat.size();
    const int num_ae = m_aerosol_mat.size();

    // Resize device vectors for species properties
    m_sp_density.resize(num_sp);
    m_sp_solubility.resize(num_sp);
    m_sp_ionization.resize(num_sp);
    m_sp_mol_weight.resize(num_sp);
    m_sp_is_INP.resize(num_sp);

    // Resize device vectors for aerosol properties
    m_ae_density.resize(num_ae);
    m_ae_solubility.resize(num_ae);
    m_ae_ionization.resize(num_ae);
    m_ae_mol_weight.resize(num_ae);
    m_ae_is_INP.resize(num_ae);

    // Create host vectors for species properties
    amrex::Vector<ParticleReal> sp_density_h(num_sp);
    amrex::Vector<int> sp_solubility_h(num_sp);
    amrex::Vector<ParticleReal> sp_ionization_h(num_sp);
    amrex::Vector<ParticleReal> sp_mol_weight_h(num_sp);
    amrex::Vector<int> sp_is_INP_h(num_sp);

    // Create host vectors for aerosol properties
    amrex::Vector<ParticleReal> ae_density_h(num_ae);
    amrex::Vector<int> ae_solubility_h(num_ae);
    amrex::Vector<ParticleReal> ae_ionization_h(num_ae);
    amrex::Vector<ParticleReal> ae_mol_weight_h(num_ae);
    amrex::Vector<int> ae_is_INP_h(num_ae);

    // Fill host vectors with species material properties
    for (int i = 0; i < num_sp; i++) {
        sp_density_h[i] = m_species_mat[i]->m_density;
        sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
        sp_ionization_h[i] = m_species_mat[i]->m_ionization;
        sp_mol_weight_h[i] = m_species_mat[i]->m_mol_weight;
        sp_is_INP_h[i] = static_cast<int>(m_species_mat[i]->m_is_INP);
    }

    // Fill host vectors with aerosol material properties
    for (int i = 0; i < num_ae; i++) {
        ae_density_h[i] = m_aerosol_mat[i]->m_density;
        ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
        ae_ionization_h[i] = m_aerosol_mat[i]->m_ionization;
        ae_mol_weight_h[i] = m_aerosol_mat[i]->m_mol_weight;
        ae_is_INP_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_INP);
    }

    // Transfer species data to device
    Gpu::copy(Gpu::hostToDevice, sp_density_h.begin(), sp_density_h.end(), m_sp_density.begin());
    Gpu::copy(Gpu::hostToDevice, sp_solubility_h.begin(), sp_solubility_h.end(), m_sp_solubility.begin());
    Gpu::copy(Gpu::hostToDevice, sp_ionization_h.begin(), sp_ionization_h.end(), m_sp_ionization.begin());
    Gpu::copy(Gpu::hostToDevice, sp_mol_weight_h.begin(), sp_mol_weight_h.end(), m_sp_mol_weight.begin());
    Gpu::copy(Gpu::hostToDevice, sp_is_INP_h.begin(), sp_is_INP_h.end(), m_sp_is_INP.begin());

    // Transfer aerosol data to device
    Gpu::copy(Gpu::hostToDevice, ae_density_h.begin(), ae_density_h.end(), m_ae_density.begin());
    Gpu::copy(Gpu::hostToDevice, ae_solubility_h.begin(), ae_solubility_h.end(), m_ae_solubility.begin());
    Gpu::copy(Gpu::hostToDevice, ae_ionization_h.begin(), ae_ionization_h.end(), m_ae_ionization.begin());
    Gpu::copy(Gpu::hostToDevice, ae_mol_weight_h.begin(), ae_mol_weight_h.end(), m_ae_mol_weight.begin());
    Gpu::copy(Gpu::hostToDevice, ae_is_INP_h.begin(), ae_is_INP_h.end(), m_ae_is_INP.begin());

    m_device_props_initialized = true;
}

/*! Update device properties if material properties change */
void SuperDropletPC::updateDeviceProperties()
{
    m_device_props_initialized = false;
    initializeDeviceProperties();
}

/*! Compute mesh variable from particles */
void SuperDropletPC::computeMeshVar( const std::string&  a_var_name,
                                     MultiFab&           a_mf,
                                     const int           a_lev) const
{
    BL_PROFILE("SuperDropletPC::computeMeshVar()");
    a_mf.setVal(0.0);

    if (a_lev != 0) { return; }

    // Check basic variables
    if (a_var_name == "number_density") {
        numberDensity(a_mf); return;
    }
    if (a_var_name == "sd_number_density") {
        SDNumberDensity(a_mf); return;
    }
    if (a_var_name == "mass_density") {
        massDensity(a_mf); return;
    }
    if (a_var_name == "radius") {
        effectiveRadius(a_mf); return;
    }

    // Check water species shortcuts
    const std::string water_name = getEnumNameString(m_species_mat[m_idx_w]->m_name);
    if (a_var_name == "mass_density_" + water_name) {
        speciesMassDensity(a_mf, m_idx_w); return;
    }
    if (a_var_name == "mass_flux_x_" + water_name) {
        speciesMassFlux(a_mf, m_idx_w, 0); return;
    }
    if (a_var_name == "mass_flux_y_" + water_name) {
        speciesMassFlux(a_mf, m_idx_w, 1); return;
    }
    if (a_var_name == "mass_flux_z_" + water_name) {
        speciesMassFlux(a_mf, m_idx_w, 2); return;
    }

    // Check total mass flux
    if (a_var_name == "mass_flux_x") { massFlux(a_mf, 0); return; }
    if (a_var_name == "mass_flux_y") { massFlux(a_mf, 1); return; }
    if (a_var_name == "mass_flux_z") { massFlux(a_mf, 2); return; }

    // Check species variables
    for (int i = 0; i < m_num_species; i++) {
        const std::string name = getEnumNameString(m_species_mat[i]->m_name);
        if (a_var_name == "species_mass_density_" + name) {
            speciesMassDensity(a_mf, i); return;
        }
        if (a_var_name == "species_mass_flux_x_" + name) {
            speciesMassFlux(a_mf, i, 0); return;
        }
        if (a_var_name == "species_mass_flux_y_" + name) {
            speciesMassFlux(a_mf, i, 1); return;
        }
        if (a_var_name == "species_mass_flux_z_" + name) {
            speciesMassFlux(a_mf, i, 2); return;
        }
    }

    // Check aerosol variables
    for (int i = 0; i < m_num_aerosols; i++) {
        const std::string name = getEnumNameString(m_aerosol_mat[i]->m_name);
        if (a_var_name == "aerosol_mass_density_" + name) {
            aerosolMassDensity(a_mf, i); return;
        }
        if (a_var_name == "aerosol_mass_flux_x_" + name) {
            aerosolMassFlux(a_mf, i, 0); return;
        }
        if (a_var_name == "aerosol_mass_flux_y_" + name) {
            aerosolMassFlux(a_mf, i, 1); return;
        }
        if (a_var_name == "aerosol_mass_flux_z_" + name) {
            aerosolMassFlux(a_mf, i, 2); return;
        }
    }
}

/*! This returns the total number of particles that all the super-droplets represent */
Real SuperDropletPC::TotalNumberOfParticles ()
{
    BL_PROFILE("SuperDropletPC::TotalNumberOfParticles()");
    Real count = ReduceSum(*this,
                           [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                           {
                                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                return ai*ni;
                           });
    ParallelDescriptor::ReduceRealSum(&count, 1);
    return count;
}

/*! This returns the total number of deactivated superdroplets (multiplicity zero) */
Long SuperDropletPC::NumSDDeactivated ()
{
    BL_PROFILE("SuperDropletPC::NumSDDeactivated()");
    auto count = ReduceSum( *this,
                            [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Long
                            {
                                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                if (ai == 0) { return Long(1); }
                                else         { return Long(0); }
                            } );
    ParallelDescriptor::ReduceLongSum(&count, 1);
    return count;
}

/*! Helper function for ParticleToMesh operations
 *
 *  This template consolidates the common boilerplate for all particle-to-mesh
 *  deposition functions. The ValueFunc should return the value to deposit
 *  (before multiplication by inv_cell_volume).
 */
template<typename ValueFunc>
void SuperDropletPC::particleToMeshHelper(MultiFab& a_mf, int a_comp, ValueFunc&& value_func) const
{
    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    a_mf.setVal(0.0);

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleTileType::ConstParticleTileDataType& ptd,
                              int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    return value_func(ptd, i) * inv_cell_volume;
                });
        });
}

/*! Computes the number density of the SDs over a mesh */
void SuperDropletPC::SDNumberDensity(MultiFab& a_mf, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::SDNumberDensity()");
    particleToMeshHelper(a_mf, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            return ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
        });
}

/*! Computes the number density of the particles over a mesh */
void SuperDropletPC::numberDensity(MultiFab& a_mf, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::numberDensity()");
    particleToMeshHelper(a_mf, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            return ai * num_par;
        });
}

/*! Computes the mass density of the particles over a mesh */
void SuperDropletPC::massDensity(MultiFab& a_mf, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::massDensity()");
    particleToMeshHelper(a_mf, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
            return ai * num_par * par_mass;
        });
}

/*! Computes the particle mass flux over a mesh */
void SuperDropletPC::massFlux(MultiFab& a_mf, const int a_dim, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::massFlux()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    a_mf.setVal(0.0);

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto par_mass = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                    auto par_velocity = ptd.m_rdata[SuperDropletsRealIdxSoA::vx+a_dim][i];
                    if (a_dim == 2) {
                        par_velocity -= ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                    }
                    return ai * num_par * par_mass * par_velocity * inv_cell_volume;
                });
        });
}

/*! Computes the aerosol mass density of the particles over a mesh */
void SuperDropletPC::aerosolMassDensity(MultiFab& a_mf, const int a_idx, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::aerosolMassDensity()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    particleToMeshHelper(a_mf, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto aero_mass = ptd.m_runtime_rdata[ridx_a(a_idx,na,ns)][i];
            return ai * num_par * aero_mass;
        });
}

/*! Computes the aerosol mass flux of the particles over a mesh */
void SuperDropletPC::aerosolMassFlux(MultiFab& a_mf, const int a_idx, const int a_dim, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::aerosolMassFlux()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;

    a_mf.setVal(0.0);

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto aero_mass = ptd.m_runtime_rdata[ridx_a(a_idx,na,ns)][i];
                    auto par_velocity = ptd.m_rdata[SuperDropletsRealIdxSoA::vx+a_dim][i];
                    if (a_dim == 2) {
                        par_velocity -= ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                    }
                    return ai * num_par * aero_mass * par_velocity * inv_cell_volume;
                });
        });
}

/*! Computes the species mass density of the particles over a mesh */
void SuperDropletPC::speciesMassDensity(MultiFab& a_mf, const int a_idx, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::speciesMassDensity()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    particleToMeshHelper(a_mf, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto species_mass = ptd.m_runtime_rdata[ridx_s(a_idx,na,ns)][i];
            return ai * num_par * species_mass;
        });
}

/*! Computes the cloud/rain mass density of the particles over a mesh */
void SuperDropletPC::cloudRainDensity(MultiFab& a_mf, const Real a_rmin, const Real a_rmax, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::cloudRainDensity()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx = m_idx_w;

    a_mf.setVal(0.0);

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    auto radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                    if ((radius < a_rmin) || (radius >= a_rmax)) {
                        return 0.0;
                    }
                    auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto species_mass = ptd.m_runtime_rdata[ridx_s(idx,na,ns)][i];
                    return ai * num_par * species_mass * inv_cell_volume;
                });
        });
}

/*! Computes ice category mass density (ice, snow, graupel, or total) over a mesh
 *
 *  Ice categories are distinguished by rime mass fraction and number of monomers:
 *  - Ice: frac < threshold AND nmono == 1
 *  - Snow: frac < threshold AND nmono > 1
 *  - Graupel: frac >= threshold
 *  - Total: all ice particles
 */
void SuperDropletPC::iceCategoryDensity(MultiFab& a_mf, IceCategory a_category,
                                        const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::iceCategoryDensity()");

    a_mf.setVal(0.0);
    if (m_idx_i < 0) { return; }

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx = m_idx_i;
    const auto category = a_category;

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    auto mass = ptd.m_runtime_rdata[ridx_s(idx,na,ns)][i];
                    auto mrime = ptd.m_runtime_rdata[ridx_ice_mrime(na,ns)][i];
                    auto nmono = ptd.m_runtime_rdata[ridx_ice_nmono(na,ns)][i];
                    auto frac = (mass > 0.0 ? mrime / mass : 0.0);

                    bool include = false;
                    switch (category) {
                        case IceCategory::Ice:
                            include = (frac < a_mrime_frac) && (nmono == 1);
                            break;
                        case IceCategory::Snow:
                            include = (frac < a_mrime_frac) && (nmono > 1);
                            break;
                        case IceCategory::Graupel:
                            include = (frac >= a_mrime_frac);
                            break;
                        case IceCategory::Total:
                            include = true;
                            break;
                    }

                    if (include) {
                        auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                        auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                        return ai * num_par * mass * inv_cell_volume;
                    }
                    return 0.0;
                });
        });
}

/*! Computes the ice mass density of the particles over a mesh */
void SuperDropletPC::iceDensity(MultiFab& a_mf, const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::iceDensity()");
    iceCategoryDensity(a_mf, IceCategory::Ice, a_mrime_frac, a_comp);
}

/*! Computes the snow mass density of the particles over a mesh */
void SuperDropletPC::snowDensity(MultiFab& a_mf, const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::snowDensity()");
    iceCategoryDensity(a_mf, IceCategory::Snow, a_mrime_frac, a_comp);
}

/*! Computes the graupel mass density of the particles over a mesh */
void SuperDropletPC::graupelDensity(MultiFab& a_mf, const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::graupelDensity()");
    iceCategoryDensity(a_mf, IceCategory::Graupel, a_mrime_frac, a_comp);
}

/*! Computes the total frozen water mass density of the particles over a mesh */
void SuperDropletPC::totalIceDensity(MultiFab& a_mf, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::totalIceDensity()");
    iceCategoryDensity(a_mf, IceCategory::Total, 0.0, a_comp);
}

/*! Computes the species mass flux of the particles over a mesh */
void SuperDropletPC::speciesMassFlux(MultiFab& a_mf, const int a_idx, const int a_dim, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::speciesMassFlux()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;

    a_mf.setVal(0.0);

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto species_mass = ptd.m_runtime_rdata[ridx_s(a_idx,na,ns)][i];
                    auto par_velocity = ptd.m_rdata[SuperDropletsRealIdxSoA::vx+a_dim][i];
                    if (a_dim == 2) {
                        par_velocity -= ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                    }
                    return ai * num_par * species_mass * par_velocity * inv_cell_volume;
                });
        });
}

/*! Computes the effective radius of the particles over a mesh */
void SuperDropletPC::effectiveRadius(MultiFab& a_mf, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::effectiveRadius()");

    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    const auto& geom = Geom(m_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];

    a_mf.setVal(0.0);

    MultiFab number_density(a_mf.boxArray(), a_mf.DistributionMap(), 1, a_mf.nGrowVect());
    numberDensity(number_density);

    ParticleToMesh(*this, a_mf, m_lev,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i, Array4<Real> const& rho)
        {
            auto p = ptd.m_aos[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho, 0, a_comp, 1,
                [=] AMREX_GPU_DEVICE (const SuperDropletPC::ParticleType&, int)
                {
                    auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                    auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    auto radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                    return ai * num_par * radius * inv_cell_volume;
                });
        });

    for (MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto mf_arr = a_mf.array(mfi);
        const auto nd_arr = number_density.const_array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (nd_arr(i,j,k,0) > 0) {
                mf_arr(i,j,k,a_comp) /= nd_arr(i,j,k,0);
            } else {
                mf_arr(i,j,k,a_comp) = 0.0;
            }
        });
    }
}

#endif
