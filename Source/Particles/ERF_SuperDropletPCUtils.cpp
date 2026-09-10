#include "ERF_SuperDropletPC.H"
#include <ERFPCParticleToMesh.H>
#include <AMReX_DenseBins.H>
#include <AMReX_Scan.H>
#include <AMReX_Reduce.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleUtil.H>

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;
using SDTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;

namespace {
//! Map a droplet radius to a log-spaced bin index over a per-cell [rmin,rmax]
//! range.  inv_lnrange = nbins / (ln rmax - ln rmin), or <= 0 when the cell is
//! monodisperse (all SDs collapse to bin 0).  Used to keep AMR merges within a
//! single size class so the droplet spectrum (and its rain-seeding tail) is
//! preserved.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int sdmRadiusBin (amrex::ParticleReal r, amrex::ParticleReal lnrmin,
                  amrex::ParticleReal inv_lnrange, int nbins)
{
    if (inv_lnrange <= amrex::ParticleReal(0) || r <= amrex::ParticleReal(0)) { return 0; }
    int b = static_cast<int>((std::log(r) - lnrmin) * inv_lnrange);
    return amrex::min(amrex::max(b, 0), nbins-1);
}
}

/*! Initialize device property arrays for species and aerosol materials */
void SuperDropletPC::initializeDeviceProperties()
{
    if (m_device_props_initialized) return;

    const int num_sp = static_cast<int>(m_species_mat.size());
    const int num_ae = static_cast<int>(m_aerosol_mat.size());

    m_sp_density.resize(num_sp);
    m_sp_solubility.resize(num_sp);
    m_sp_ionization.resize(num_sp);
    m_sp_mol_weight.resize(num_sp);
    m_sp_is_INP.resize(num_sp);
    m_ae_density.resize(num_ae);
    m_ae_solubility.resize(num_ae);
    m_ae_ionization.resize(num_ae);
    m_ae_mol_weight.resize(num_ae);
    m_ae_is_INP.resize(num_ae);

    amrex::Vector<ParticleReal> sp_density_h(num_sp), sp_ionization_h(num_sp), sp_mol_weight_h(num_sp);
    amrex::Vector<int> sp_solubility_h(num_sp), sp_is_INP_h(num_sp);
    amrex::Vector<ParticleReal> ae_density_h(num_ae), ae_ionization_h(num_ae), ae_mol_weight_h(num_ae);
    amrex::Vector<int> ae_solubility_h(num_ae), ae_is_INP_h(num_ae);

    for (int i = 0; i < num_sp; i++) {
        sp_density_h[i] = m_species_mat[i]->m_density;
        sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
        sp_ionization_h[i] = m_species_mat[i]->m_ionization;
        sp_mol_weight_h[i] = m_species_mat[i]->m_mol_weight;
        sp_is_INP_h[i] = static_cast<int>(m_species_mat[i]->m_is_INP);
    }
    for (int i = 0; i < num_ae; i++) {
        ae_density_h[i] = m_aerosol_mat[i]->m_density;
        ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
        ae_ionization_h[i] = m_aerosol_mat[i]->m_ionization;
        ae_mol_weight_h[i] = m_aerosol_mat[i]->m_mol_weight;
        ae_is_INP_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_INP);
    }

    Gpu::copy(Gpu::hostToDevice, sp_density_h.begin(), sp_density_h.end(), m_sp_density.begin());
    Gpu::copy(Gpu::hostToDevice, sp_solubility_h.begin(), sp_solubility_h.end(), m_sp_solubility.begin());
    Gpu::copy(Gpu::hostToDevice, sp_ionization_h.begin(), sp_ionization_h.end(), m_sp_ionization.begin());
    Gpu::copy(Gpu::hostToDevice, sp_mol_weight_h.begin(), sp_mol_weight_h.end(), m_sp_mol_weight.begin());
    Gpu::copy(Gpu::hostToDevice, sp_is_INP_h.begin(), sp_is_INP_h.end(), m_sp_is_INP.begin());
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
void SuperDropletPC::computeMeshVar( const std::string& a_var_name,
                                     MultiFab&          a_mf,
                                     const MultiFab&    a_z_phys_nd,
                                     const int          a_lev) const
{
    BL_PROFILE("SuperDropletPC::computeMeshVar()");
    a_mf.setVal(0.0);

    // Check basic variables
    if (a_var_name == "number_density") {
        numberDensity(a_mf, a_z_phys_nd, a_lev); return;
    }
    if (a_var_name == "sd_number_density") {
        SDNumberDensity(a_mf, a_z_phys_nd, a_lev); return;
    }
    if (a_var_name == "mass_density") {
        massDensity(a_mf, a_z_phys_nd, a_lev); return;
    }
    if (a_var_name == "radius") {
        effectiveRadius(a_mf, a_z_phys_nd, a_lev); return;
    }

    // Check water species shortcuts
    const std::string water_name = getEnumNameString(m_species_mat[m_idx_w]->m_name);
    if (a_var_name == "mass_density_" + water_name) {
        speciesMassDensity(a_mf, a_z_phys_nd, a_lev, m_idx_w); return;
    }
    if (a_var_name == "mass_flux_x_" + water_name) {
        speciesMassFlux(a_mf, a_z_phys_nd, a_lev, m_idx_w, 0); return;
    }
    if (a_var_name == "mass_flux_y_" + water_name) {
        speciesMassFlux(a_mf, a_z_phys_nd, a_lev, m_idx_w, 1); return;
    }
    if (a_var_name == "mass_flux_z_" + water_name) {
        speciesMassFlux(a_mf, a_z_phys_nd, a_lev, m_idx_w, 2); return;
    }

    // Check total mass flux
    if (a_var_name == "mass_flux_x") { massFlux(a_mf, a_z_phys_nd, a_lev, 0); return; }
    if (a_var_name == "mass_flux_y") { massFlux(a_mf, a_z_phys_nd, a_lev, 1); return; }
    if (a_var_name == "mass_flux_z") { massFlux(a_mf, a_z_phys_nd, a_lev, 2); return; }

    // Check species variables
    for (int i = 0; i < m_num_species; i++) {
        const std::string name = getEnumNameString(m_species_mat[i]->m_name);
        if (a_var_name == "species_mass_density_" + name) {
            speciesMassDensity(a_mf, a_z_phys_nd, a_lev, i); return;
        }
        if (a_var_name == "species_mass_flux_x_" + name) {
            speciesMassFlux(a_mf, a_z_phys_nd, a_lev, i, 0); return;
        }
        if (a_var_name == "species_mass_flux_y_" + name) {
            speciesMassFlux(a_mf, a_z_phys_nd, a_lev, i, 1); return;
        }
        if (a_var_name == "species_mass_flux_z_" + name) {
            speciesMassFlux(a_mf, a_z_phys_nd, a_lev, i, 2); return;
        }
    }

    // Check aerosol variables
    for (int i = 0; i < m_num_aerosols; i++) {
        const std::string name = getEnumNameString(m_aerosol_mat[i]->m_name);
        if (a_var_name == "aerosol_mass_density_" + name) {
            aerosolMassDensity(a_mf, a_z_phys_nd, a_lev, i); return;
        }
        if (a_var_name == "aerosol_mass_flux_x_" + name) {
            aerosolMassFlux(a_mf, a_z_phys_nd, a_lev, i, 0); return;
        }
        if (a_var_name == "aerosol_mass_flux_y_" + name) {
            aerosolMassFlux(a_mf, a_z_phys_nd, a_lev, i, 1); return;
        }
        if (a_var_name == "aerosol_mass_flux_z_" + name) {
            aerosolMassFlux(a_mf, a_z_phys_nd, a_lev, i, 2); return;
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
                                int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
                                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                return static_cast<Real>(ai*ni);
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
                                int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
                                if (ai == 0) { return Long(1); }
                                else         { return Long(0); }
                            } );
    ParallelDescriptor::ReduceLongSum(&count, 1);
    return count;
}

/*! Computes the number density of the SDs over a mesh */
void SuperDropletPC::SDNumberDensity ( MultiFab& a_mf,
                                       const MultiFab& a_z_phys_nd,
                                       int a_lev,
                                       const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::SDNumberDensity()");
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            return (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
        });
}

/*! Computes the number density of the particles over a mesh */
void SuperDropletPC::numberDensity ( MultiFab& a_mf,
                                     const MultiFab& a_z_phys_nd,
                                     int a_lev,
                                     const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::numberDensity()");
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            return ai * num_par;
        });
}

/*! Computes the mass density of the particles over a mesh */
void SuperDropletPC::massDensity ( MultiFab& a_mf,
                                   const MultiFab& a_z_phys_nd,
                                   const int& a_lev,
                                   const int& a_comp ) const
{
    BL_PROFILE("SuperDropletPC::massDensity()");
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto par_mass = ptd.m_rdata[SuperDropletsRealIdx::mass][i];
            return ai * num_par * par_mass;
        });
}

/*! Computes the particle velocity components over a mesh */
void SuperDropletPC::massFlux ( MultiFab& a_mf,
                                const MultiFab& a_z_phys_nd,
                                int a_lev,
                                const int a_dim,
                                const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::massFlux()");
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto par_mass = ptd.m_rdata[SuperDropletsRealIdx::mass][i];
            auto par_velocity = ptd.m_rdata[SuperDropletsRealIdx::vx+a_dim][i];
            if (a_dim == 2) {
                par_velocity -= ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
            }
            return ai * num_par * par_mass * par_velocity;
        });
}

/*! Computes the aerosol mass density of the particles over a mesh */
void SuperDropletPC::aerosolMassDensity ( MultiFab& a_mf,
                                          const MultiFab& a_z_phys_nd,
                                          int a_lev,
                                          const int a_idx,
                                          const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::aerosolMassDensity()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto aero_mass = ptd.m_runtime_rdata[ridx_a(a_idx,na,ns)][i];
            return ai * num_par * aero_mass;
        });
}

/*! Computes the aerosol mass flux of the particles over a mesh */
void SuperDropletPC::aerosolMassFlux ( MultiFab& a_mf,
                                       const MultiFab& a_z_phys_nd,
                                       int a_lev,
                                       const int a_idx,
                                       const int a_dim,
                                       const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::aerosolMassFlux()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto aero_mass = ptd.m_runtime_rdata[ridx_a(a_idx,na,ns)][i];
            auto par_velocity = ptd.m_rdata[SuperDropletsRealIdx::vx+a_dim][i];
            if (a_dim == 2) {
                par_velocity -= ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
            }
            return ai * num_par * aero_mass * par_velocity;
        });
}

/*! Computes the species mass density of the particles over a mesh */
void SuperDropletPC::speciesMassDensity ( MultiFab&  a_mf,
                                          const MultiFab& a_z_phys_nd,
                                          int        a_lev,
                                          int        a_idx,
                                          const int  a_comp ) const
{
    BL_PROFILE("SuperDropletPC::speciesMassDensity()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto species_mass = ptd.m_runtime_rdata[ridx_s(a_idx,na,ns)][i];
            return ai * num_par * species_mass;
        });
}

/*! Computes the cloud/rain mass density of the particles over a mesh */
void SuperDropletPC::cloudRainDensity(MultiFab& a_mf, const MultiFab& a_z_phys_nd, int a_lev, const Real a_rmin, const Real a_rmax, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::cloudRainDensity()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx = m_idx_w;
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            auto radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
            if ((radius < a_rmin) || (radius >= a_rmax)) {
                return ParticleReal(0);
            }
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto species_mass = ptd.m_runtime_rdata[ridx_s(idx,na,ns)][i];
            return ai * num_par * species_mass;
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
void SuperDropletPC::iceCategoryDensity(MultiFab& a_mf, const MultiFab& a_z_phys_nd,
                                        int a_lev, IceCategory a_category,
                                        const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::iceCategoryDensity()");

    if (m_idx_i < 0) { a_mf.setVal(0.0); return; }

    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx = m_idx_i;
    const auto category = a_category;
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            auto mass = ptd.m_runtime_rdata[ridx_s(idx,na,ns)][i];
            auto mrime = ptd.m_runtime_rdata[ridx_ice_mrime(na,ns)][i];
            auto nmono = ptd.m_runtime_rdata[ridx_ice_nmono(na,ns)][i];
            auto frac = (mass > amrex::ParticleReal(0) ? mrime / mass : amrex::ParticleReal(0));

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

            if (!include) { return amrex::ParticleReal(0); }

            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            return ai * num_par * mass;
        });
}

/*! Computes the ice mass density of the particles over a mesh */
void SuperDropletPC::iceDensity(MultiFab& a_mf, const MultiFab& a_z_phys_nd,
                                int a_lev, const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::iceDensity()");
    iceCategoryDensity(a_mf, a_z_phys_nd, a_lev, IceCategory::Ice, a_mrime_frac, a_comp);
}

/*! Computes the snow mass density of the particles over a mesh */
void SuperDropletPC::snowDensity(MultiFab& a_mf, const MultiFab& a_z_phys_nd,
                                 int a_lev, const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::snowDensity()");
    iceCategoryDensity(a_mf, a_z_phys_nd, a_lev, IceCategory::Snow, a_mrime_frac, a_comp);
}

/*! Computes the graupel mass density of the particles over a mesh */
void SuperDropletPC::graupelDensity(MultiFab& a_mf, const MultiFab& a_z_phys_nd,
                                    int a_lev, const Real a_mrime_frac, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::graupelDensity()");
    iceCategoryDensity(a_mf, a_z_phys_nd, a_lev, IceCategory::Graupel, a_mrime_frac, a_comp);
}

/*! Computes the total frozen water mass density of the particles over a mesh */
void SuperDropletPC::totalIceDensity(MultiFab& a_mf, const MultiFab& a_z_phys_nd,
                                     int a_lev, const int a_comp) const
{
    BL_PROFILE("SuperDropletPC::totalIceDensity()");
    iceCategoryDensity(a_mf, a_z_phys_nd, a_lev, IceCategory::Total, 0.0, a_comp);
}

/*! Computes the species mass flux of the particles over a mesh */
void SuperDropletPC::speciesMassFlux ( MultiFab& a_mf,
                                       const MultiFab& a_z_phys_nd,
                                       int a_lev,
                                       const int a_idx,
                                       const int a_dim,
                                       const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::speciesMassFlux()");
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto species_mass = ptd.m_runtime_rdata[ridx_s(a_idx,na,ns)][i];
            auto par_velocity = ptd.m_rdata[SuperDropletsRealIdx::vx+a_dim][i];
            if (a_dim == 2) {
                par_velocity -= ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
            }
            return ai * num_par * species_mass * par_velocity;
        });
}

/*! Computes the effective radius of the particles over a mesh */
void SuperDropletPC::effectiveRadius (  MultiFab& a_mf,
                                        const MultiFab& a_z_phys_nd,
                                        int a_lev,
                                        const int a_comp ) const
{
    BL_PROFILE("SuperDropletPC::effectiveRadius()");

    MultiFab number_density(a_mf.boxArray(), a_mf.DistributionMap(), 1, a_mf.nGrowVect());
    numberDensity(number_density, a_z_phys_nd, a_lev);

    ERFPCParticleToMesh(a_mf, a_z_phys_nd, a_lev, a_comp,
        [=] AMREX_GPU_DEVICE (const SDTDType& ptd, int i) {
            int ai = (ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i] > 0) ? 1 : 0;
            auto num_par = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
            auto radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
            return ai * num_par * radius;
        });

    for (MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();
        auto mf_arr = a_mf.array(mfi);
        const auto nd_arr = number_density.const_array(mfi);
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                          {
                              if (nd_arr(i,j,k,0) > 0) {
                                  mf_arr(i,j,k,a_comp) /= nd_arr(i,j,k,0);
                              } else {
                                  mf_arr(i,j,k,a_comp) = zero;
                              }
                          } );
    }
}

/*! \brief Split super-droplets that ended up on a level deeper than their
 *  tag indicates (called once after the post-regrid Redistribute moves SDs
 *  to their leaf level).  For every (target_lev, source_tag) pair where
 *  source_tag < target_lev + 1, the SDs with `tag == source_tag` on
 *  `target_lev` are split by the cumulative refinement factor
 *  prod_{j=source_tag-1}^{target_lev-1} refRatio(j) and tagged with
 *  `target_lev + 1` (the destination level's native tag).  Daughter
 *  positions are jittered across the source-level parent cell so the
 *  following Redistribute can spread them into the appropriate finer
 *  sub-cells.
 *
 *  This single call handles cascading multi-level refinement (e.g. a fresh
 *  regrid creating both L1 and L2 from L0): an L0-native SD that
 *  Redistribute placed on L2 is split by refRatio(0)*refRatio(1), so the
 *  next Redistribute lands the right SD density on every level.
 *
 *  Only particles whose multiplicity is >= the relevant split factor are
 *  split; lower-multiplicity particles are skipped. */
void SuperDropletPC::SplitParticlesForRefinement (int finest_level)
{
    BL_PROFILE("SuperDropletPC::SplitParticlesForRefinement()");

    if (!m_split_merge_amr) { return; }
    if (finest_level < 1) { return; }

    constexpr int rt_off_r = SuperDropletsRealIdx::ncomps;
    constexpr int rt_off_i = SuperDropletsIntIdx::ncomps;
    const int mult_idx   = rt_off_r + SuperDropletsRealIdxSoA_RT::multiplicity;
    const int active_idx = rt_off_i + SuperDropletsIntIdxSoA_RT::active;
    const int my_proc = ParallelDescriptor::MyProc();

    for (int lev = 1; lev <= finest_level; lev++) {

        const int lev_native = lev + 1;

        for (int source_tag = 1; source_tag < lev_native; source_tag++) {
            const int source_lev = source_tag - 1;

            IntVect cum_ref = IntVect::TheUnitVector();
            for (int j = source_lev; j < lev; j++) {
                cum_ref *= m_gdb->refRatio(j);
            }
            const int split_factor = AMREX_D_TERM(cum_ref[0],
                                                  *cum_ref[1],
                                                  *cum_ref[2]);
            if (split_factor <= 1) { continue; }
            const int daughter_tag = lev_native;
            const int sf_count = split_factor;
            const int src_tag = source_tag;

            // Parent-cell size at the source level; daughter positions get
            // jittered across this range so the next Redistribute spreads
            // them into the appropriate finer sub-cells.
            const auto src_dx = m_gdb->Geom(source_lev).CellSizeArray();
            const auto src_plo = m_gdb->Geom(source_lev).ProbLoArray();
            const auto src_dxi = m_gdb->Geom(source_lev).InvCellSizeArray();

            Long n_to_split         = 0;
            Long n_skipped_low_mult = 0;

            for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti) {
                const auto& ptile = ParticlesAt(lev, pti);
                const auto& aos   = ptile.GetArrayOfStructs();
                int np = aos.numParticles();
                if (np == 0) { continue; }

                const auto* p_pbox = aos().data();
                const auto& soa    = ptile.GetStructOfArrays();
                const auto* mult_ptr   = soa.GetRealData(mult_idx).data();
                const auto* active_ptr = soa.GetIntData(active_idx).data();

                n_to_split += static_cast<Long>(amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int {
                        return (p_pbox[i].id() > 0
                                && active_ptr[i] == src_tag
                                && mult_ptr[i] >= static_cast<ParticleReal>(sf_count)) ? 1 : 0;
                    }));
                n_skipped_low_mult += static_cast<Long>(amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int {
                        return (p_pbox[i].id() > 0
                                && active_ptr[i] == src_tag
                                && mult_ptr[i] < static_cast<ParticleReal>(sf_count)) ? 1 : 0;
                    }));
            }
            ParallelDescriptor::ReduceLongSum(n_to_split);
            ParallelDescriptor::ReduceLongSum(n_skipped_low_mult);

            Long n_new_particles = n_to_split * (static_cast<Long>(split_factor) - 1);
            if (n_new_particles == 0) {
                if (n_skipped_low_mult > 0) {
                    Print() << "SplitParticlesForRefinement: L" << source_lev
                            << " -> L" << lev << " -- skipped "
                            << n_skipped_low_mult
                            << " particles with multiplicity < "
                            << split_factor << "\n";
                }
                continue;
            }

            {
                size_t bytes_per_p = sizeof(ParticleType)
                                   + static_cast<size_t>(NumRealComps()
                                       + NumRuntimeRealComps()) * sizeof(ParticleReal)
                                   + static_cast<size_t>(NumIntComps()
                                       + NumRuntimeIntComps()) * sizeof(int);
                size_t bytes_needed = static_cast<size_t>(n_new_particles) * bytes_per_p;
#ifdef AMREX_USE_GPU
                size_t free_mem = Gpu::Device::freeMemAvailable();
                if (bytes_needed > static_cast<size_t>(static_cast<double>(free_mem) * 0.8)) {
                    Print() << "WARNING: SplitParticlesForRefinement: not enough GPU memory for L"
                            << source_lev << " -> L" << lev << ".  Need "
                            << bytes_needed/(1024*1024) << " MB, available "
                            << free_mem/(1024*1024) << " MB.  Skipping.\n";
                    continue;
                }
#endif
                amrex::ignore_unused(bytes_needed);
            }

            const int sf = split_factor;

            for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
                auto& ptile = ParticlesAt(lev, pti);
                auto& aos   = ptile.GetArrayOfStructs();
                int np = aos.numParticles();
                if (np == 0) { continue; }

                Gpu::DeviceVector<int> split_mask_d(np);
                auto* mask_d_ptr = split_mask_d.data();
                auto* p_pbox = aos().data();
                auto& soa = ptile.GetStructOfArrays();
                const auto* mult_data   = soa.GetRealData(mult_idx).data();
                const auto* active_data = soa.GetIntData(active_idx).data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    mask_d_ptr[i] = (p_pbox[i].id() > 0
                                     && active_data[i] == src_tag
                                     && mult_data[i] >= static_cast<ParticleReal>(sf)) ? 1 : 0;
                });

                int n_split_tile = amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int { return mask_d_ptr[i]; });
                if (n_split_tile == 0) { continue; }

                ParticleTileType src_tile;
                src_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                src_tile.resize(n_split_tile);
                amrex::filterParticles(src_tile, ptile, mask_d_ptr, 0, 0, np);

                int n_copies = n_split_tile * (sf - 1);
                ParticleTileType copy_tile;
                copy_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                copy_tile.resize(n_copies);
                for (int c = 0; c < sf - 1; c++) {
                    amrex::copyParticles(copy_tile, src_tile, 0, c * n_split_tile, n_split_tile);
                }

                Long pid = ParticleType::NextID();
                ParticleType::NextID(pid + static_cast<Long>(n_copies));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    static_cast<Long>(pid + n_copies) < LastParticleID,
                    "particle id overflow in SplitParticlesForRefinement");

                auto& copy_aos = copy_tile.GetArrayOfStructs();
                auto* copy_p = copy_aos().data();
                const Long pid_start = pid;
                const int proc = my_proc;
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    copy_p[i].id()  = static_cast<int>(pid_start + i);
                    copy_p[i].cpu() = proc;
                });

                auto& copy_soa = copy_tile.GetStructOfArrays();
                auto* copy_mult = copy_soa.GetRealData(mult_idx).data();

                const int nst = n_split_tile;
                const int kf  = split_factor;
                const auto local_src_dx  = src_dx;
                const auto local_src_plo = src_plo;
                const auto local_src_dxi = src_dxi;

                // Jitter daughter positions uniformly across the source-level
                // parent cell.  The next Redistribute spreads them into the
                // appropriate finer sub-cells of that parent.
                amrex::ParallelForRNG(n_copies,
                    [=] AMREX_GPU_DEVICE (int i, const amrex::RandomEngine& rng) {
                    for (int d = 0; d < AMREX_SPACEDIM; d++) {
                        int src_cell = static_cast<int>(amrex::Math::floor(
                            (copy_p[i].pos(d) - local_src_plo[d]) * local_src_dxi[d]));
                        ParticleReal cell_lo = local_src_plo[d] + src_cell * local_src_dx[d];
                        ParticleReal cell_hi = cell_lo + local_src_dx[d];
                        ParticleReal u = amrex::Random(rng);
                        ParticleReal pos_new = cell_lo + u * local_src_dx[d];
                        constexpr ParticleReal eps = ParticleReal(1.0e-10);
                        pos_new = amrex::max(pos_new, cell_lo + eps);
                        pos_new = amrex::min(pos_new, cell_hi - eps);
                        copy_p[i].pos(d) = pos_new;
                    }
                });

                // Integer multiplicity split.  Each of the kf pieces (the
                // original + kf-1 copies) gets base = floor(mult/kf) real
                // droplets; the remainder mult - kf*base is handed out one
                // droplet each to the first `rem` pieces (piece 0 = original).
                // Every piece is an integer >= 1 and the parent multiplicity is
                // conserved, so coalescence cannot drive super-droplets below
                // one real droplet.
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    long xil  = static_cast<long>(copy_mult[i] + ParticleReal(0.5));
                    long base = xil / kf;
                    long rem  = xil - base * kf;
                    int  p    = (i / nst) + 1;            // copy piece index (1..kf-1)
                    copy_mult[i] = static_cast<ParticleReal>(base + (p < rem ? 1 : 0));
                });

                auto* orig_mult   = soa.GetRealData(mult_idx).data();
                auto* orig_active = soa.GetIntData(active_idx).data();
                const int dtag = daughter_tag;
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    if (mask_d_ptr[i]) {
                        long xil  = static_cast<long>(orig_mult[i] + ParticleReal(0.5));
                        long base = xil / kf;
                        long rem  = xil - base * kf;
                        orig_mult[i]   = static_cast<ParticleReal>(base + (rem > 0 ? 1 : 0));
                        orig_active[i]  = dtag;
                    }
                });
                Gpu::synchronize();

                auto old_np = ptile.numParticles();
                ptile.resize(old_np + n_copies);
                amrex::copyParticles(ptile, copy_tile, 0, old_np, n_copies);

                auto* copy_active_data = ptile.GetStructOfArrays()
                    .GetIntData(active_idx).data();
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    copy_active_data[old_np + i] = dtag;
                });
            }

            Print() << "SplitParticlesForRefinement: L" << source_lev
                    << " -> L" << lev << " split " << n_to_split
                    << " tag=" << source_tag << " particles by factor "
                    << split_factor
                    << " (created " << n_new_particles << " new SDs)\n";
            if (n_skipped_low_mult > 0) {
                Print() << "  (skipped " << n_skipped_low_mult
                        << " particles with multiplicity < "
                        << split_factor << ")\n";
            }
        }
    }
}

/*! \brief Regrid-time merge in cells that lost fine-level coverage.  Builds
 *  the de-refinement mask on `a_lev - 1` as old_fine_BA minus current_fine_BA
 *  (refined to clev) and, for masked cells, reduces `np_bin` super-droplets
 *  to `np_bin / merge_factor` via multiplicity-weighted host-cycling merge.
 *  Excess SDs are deactivated (id = -1); the surviving heads retain summed
 *  multiplicity and multiplicity-weighted averaged attributes.  Without this
 *  pass, oscillating L_{n} -> L_{n+1} -> L_{n} cycles in moist runs leak SDs
 *  by `merge_factor` per cycle. */
void SuperDropletPC::MergeParticlesAtDerefinement (
    int a_lev,
    const amrex::BoxArray& a_old_fine_ba,
    const amrex::IntVect& a_ref_ratio)
{
    BL_PROFILE("SuperDropletPC::MergeParticlesAtDerefinement()");

    if (!m_split_merge_amr) { return; }

    const int merge_factor = AMREX_D_TERM(a_ref_ratio[0], *a_ref_ratio[1], *a_ref_ratio[2]);
    if (merge_factor <= 1) { return; }
    if (a_old_fine_ba.empty()) { return; }

    const int clev = a_lev - 1;
    AMREX_ALWAYS_ASSERT(clev >= 0);

    const auto& cba = ParticleBoxArray(clev);
    const auto& cdm = ParticleDistributionMap(clev);
    iMultiFab mask = makeFineMask(cba, cdm, a_old_fine_ba, a_ref_ratio,
                                  /*crse_value=*/0, /*fine_value=*/1);

    // Subtract cells still covered by the current fine level: only cells that
    // genuinely lost coverage need to be merged.
    if (a_lev <= finestLevel()) {
        const BoxArray& cur_fine_ba = ParticleBoxArray(a_lev);
        iMultiFab cur_mask = makeFineMask(cba, cdm, cur_fine_ba, a_ref_ratio,
                                          /*crse_value=*/0, /*fine_value=*/1);
        for (MFIter mfi(mask); mfi.isValid(); ++mfi) {
            auto old_arr = mask[mfi].array();
            auto cur_arr = cur_mask[mfi].const_array();
            const Box& bx = mfi.validbox();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (cur_arr(i,j,k)) { old_arr(i,j,k) = 0; }
            });
        }
        Gpu::synchronize();
    }

    const auto& geom = m_gdb->Geom(clev);
    const auto plo    = geom.ProbLoArray();
    const auto dxi    = geom.InvCellSizeArray();
    const auto domain = geom.Domain();

    const int num_sp = m_num_species;
    const int num_ae = m_num_aerosols;
    const int mf     = merge_factor;
    const int n_rbin = (m_num_sd_per_cell >= 2) ? amrex::min(m_num_sd_per_cell, 256) : 60;
    const auto ctx   = buildProcessContext(clev);
    const IntVect bin_size = {AMREX_D_DECL(1,1,1)};
    Long n_merged = 0;

    forEachParticleTileSerial(clev, ctx,
        [&](ParIterType& /*pti*/, int grid, ParticleType* pstruct_ptr,
            const SDProcess::ParticlePointers& ptrs,
            const SDProcess::ProcessContext& /*ctx_unused*/)
    {
        const size_t np = static_cast<size_t>(ptrs.num_particles);
        if (np == 0) { return; }

        Box box = mask[grid].box();
        int ntiles = numTilesInBox(box, true, bin_size);
        auto binner = GetParticleBinERF{plo, dxi, domain, bin_size, box};
        DenseBins<ParticleType> bins;
        bins.build(np, pstruct_ptr, ntiles, binner);
        auto inds    = bins.permutationPtr();
        auto offsets = bins.offsetsPtr();

        auto mask_arr = mask[grid].const_array();

        auto* mult_p   = ptrs.mult_ptr;
        auto* v_p0     = ptrs.v_ptr[0];
        auto* v_p1     = ptrs.v_ptr[1];
        auto* v_p2     = ptrs.v_ptr[2];
        auto* radius_p = ptrs.radius_ptr;
        auto* mass_p   = ptrs.mass_ptr;
        auto sp_mass_p = ptrs.sp_mass_ptrs;
        auto ae_mass_p = ptrs.ae_mass_ptrs;
        auto* sp_rho_p = ptrs.sp_rho_arr;
        auto* sp_sol_p = ptrs.sp_sol_arr;
        auto* ae_rho_p = ptrs.ae_rho_arr;
        auto* ae_sol_p = ptrs.ae_sol_arr;
        auto idx_w     = ctx.idx_water;
        auto rho_w     = ctx.rho_water;

        int num_bins = bins.numBins();

        auto tile_merged = amrex::Reduce::Sum<Long>(num_bins,
            [=] AMREX_GPU_DEVICE (int i_bin) -> Long
        {
            auto bin_start = offsets[i_bin];
            auto bin_stop  = offsets[i_bin+1];
            int np_bin = static_cast<int>(bin_stop - bin_start);
            if (np_bin < 2) { return 0; }

            unsigned int first_idx = inds[bin_start];
            IntVect iv = amrex::getParticleCell(pstruct_ptr[first_idx], plo, dxi, domain);
            if (!mask_arr(iv)) { return 0; }

            // Per-cell adaptive radius range over live SDs.
            ParticleReal rmin_c = ParticleReal(-1);
            ParticleReal rmax_c = ParticleReal(-1);
            for (int k = 0; k < np_bin; k++) {
                unsigned int idx = inds[bin_start + k];
                if (pstruct_ptr[idx].id() <= 0) { continue; }
                ParticleReal r = radius_p[idx];
                if (r > ParticleReal(0)) {
                    if (rmin_c < ParticleReal(0) || r < rmin_c) { rmin_c = r; }
                    if (r > rmax_c) { rmax_c = r; }
                }
            }
            if (rmax_c <= ParticleReal(0)) { return 0; }
            const ParticleReal lnrmin  = std::log(rmin_c);
            const ParticleReal lnrange = std::log(rmax_c) - lnrmin;
            const ParticleReal inv_lnr = (lnrange > ParticleReal(1.0e-12))
                                       ? (static_cast<ParticleReal>(n_rbin) / lnrange)
                                       : ParticleReal(0);

            constexpr int KEEP_MAX = 16;
            Long bin_merged = 0;

            // Reduce each occupied radius bin to max(1, round(count/mf)),
            // merging only same-size super-droplets.  Singletons and the
            // sparsely-populated large-droplet tail are left untouched.
            for (int rb = 0; rb < n_rbin; rb++) {
                int cnt = 0;
                for (int k = 0; k < np_bin; k++) {
                    unsigned int idx = inds[bin_start + k];
                    if (pstruct_ptr[idx].id() <= 0) { continue; }
                    if (sdmRadiusBin(radius_p[idx], lnrmin, inv_lnr, n_rbin) == rb) { cnt++; }
                }
                if (cnt < 2) { continue; }
                int keep = (cnt + mf/2) / mf;
                if (keep < 1) { keep = 1; }
                if (keep > KEEP_MAX) { keep = KEEP_MAX; }
                if (cnt <= keep) { continue; }

                int surv[KEEP_MAX];
                int ns = 0, cyc = 0;
                for (int k = 0; k < np_bin; k++) {
                    unsigned int idx = inds[bin_start + k];
                    if (pstruct_ptr[idx].id() <= 0) { continue; }
                    if (sdmRadiusBin(radius_p[idx], lnrmin, inv_lnr, n_rbin) != rb) { continue; }
                    if (ns < keep) { surv[ns++] = static_cast<int>(idx); continue; }
                    unsigned int h = static_cast<unsigned int>(surv[cyc % keep]); cyc++;

                    ParticleReal xi_e = mult_p[idx];
                    ParticleReal xi_s = mult_p[h];
                    ParticleReal xi_new = xi_e + xi_s;
                    if (xi_new <= ParticleReal(0)) { continue; }
                    ParticleReal inv_xi = ParticleReal(1.0) / xi_new;

                    mult_p[h] = xi_new;
                    v_p0[h] = (xi_e * v_p0[idx] + xi_s * v_p0[h]) * inv_xi;
                    v_p1[h] = (xi_e * v_p1[idx] + xi_s * v_p1[h]) * inv_xi;
                    v_p2[h] = (xi_e * v_p2[idx] + xi_s * v_p2[h]) * inv_xi;
                    for (int s = 0; s < num_sp; s++) {
                        sp_mass_p[s][h] = (xi_e * sp_mass_p[s][idx]
                                         + xi_s * sp_mass_p[s][h]) * inv_xi;
                    }
                    for (int a = 0; a < num_ae; a++) {
                        ae_mass_p[a][h] = (xi_e * ae_mass_p[a][idx]
                                         + xi_s * ae_mass_p[a][h]) * inv_xi;
                    }
                    pstruct_ptr[idx].id() = -1;
                    bin_merged++;
                }
            }

            for (int k = 0; k < np_bin; k++) {
                unsigned int idx = inds[bin_start + k];
                if (pstruct_ptr[idx].id() > 0) {
                    updateParticleAttributes(
                        idx, radius_p, mass_p,
                        idx_w, rho_w, num_sp, num_ae,
                        sp_sol_p, ae_sol_p,
                        sp_mass_p, ae_mass_p,
                        sp_rho_p, ae_rho_p);
                }
            }

            return bin_merged;
        });
        Gpu::synchronize();
        n_merged += tile_merged;
    });

    ParallelDescriptor::ReduceLongSum(n_merged);
    if (n_merged > 0) {
        Redistribute(clev, clev);
        Print() << "MergeParticlesAtDerefinement: merged " << n_merged
                << " excess particles on L" << clev
                << " (de-refined from level " << a_lev << ")\n";
    }
}

/*! \brief Continuously split new entrants on fine levels, merge departees on
 *  coarse levels, and clean up leftover tags from disappeared levels.  Uses
 *  the tag = level + 1 convention (0 deactivated).
 *
 *  Two loops:
 *    - Per-(clev, flev) pair (flev = 1..finest, clev = flev-1):
 *        Part 1: split new entrants on flev (active particles whose tag is
 *                below flev's native value `flev + 1`) by refRatio(clev) and
 *                tag them native to flev.  Entrants with multiplicity below
 *                split_factor are retagged native without splitting.
 *        Part 3: repopulate empty fine cells on flev by cloning a face-
 *                adjacent neighbor inside the same coarse parent.
 *    - Per-level (lev = 0..finest):
 *        Part 2: in each cell on `lev`, merge any super-droplets with tag
 *                > lev+1 (leftovers from a finer level that has since drifted
 *                away or vanished) into a native (tag == lev+1) host via
 *                multiplicity-weighted host cycling.  If no native host is
 *                present, the over-tag SDs are retagged native without
 *                merging; count reduction at de-refinement is handled by
 *                MergeParticlesAtDerefinement at regrid time.
 */
void SuperDropletPC::SplitMergeAtLevelBoundary ()
{
    BL_PROFILE("SuperDropletPC::SplitMergeAtLevelBoundary()");

    if (!m_split_merge_amr) { return; }
    int finest = finestLevel();

    constexpr int rtoff_i = SuperDropletsIntIdx::ncomps;
    constexpr int rtoff_r = SuperDropletsRealIdx::ncomps;
    const int mult_idx   = rtoff_r + SuperDropletsRealIdxSoA_RT::multiplicity;
    const int active_idx = rtoff_i + SuperDropletsIntIdxSoA_RT::active;
    const int my_proc = ParallelDescriptor::MyProc();

    for (int flev = 1; flev <= finest; flev++) {
        const int clev = flev - 1;
        const int flev_native = flev + 1;
        const IntVect ref_ratio = m_gdb->refRatio(clev);
        const int split_factor = AMREX_D_TERM(ref_ratio[0], *ref_ratio[1], *ref_ratio[2]);
        if (split_factor <= 1) { continue; }

        // ---- Part 1: split new entrants on the fine level ----
        //
        // A new entrant on `flev` is any active particle whose tag is below
        // this level's native value `flev + 1` -- it has not yet been split
        // to match flev's density.  Common case (per-step drift L_{flev-1} ->
        // L_flev): tag == clev + 1 == flev.
        {
            Long n_new_entrants = 0;

            for (ParIterType pti(*this, flev); pti.isValid(); ++pti) {
                auto& ptile = ParticlesAt(flev, pti);
                auto& aos = ptile.GetArrayOfStructs();
                int np = aos.numParticles();
                if (np == 0) { continue; }

                auto& soa = ptile.GetStructOfArrays();
                auto* p_pbox = aos().data();
                auto* active_data = soa.GetIntData(active_idx).data();
                const auto* mult_data = soa.GetRealData(mult_idx).data();
                const ParticleReal sf_floor = static_cast<ParticleReal>(split_factor);

                // Multiplicity floor: an entrant with multiplicity < split_factor
                // cannot be split into split_factor daughters without driving
                // each below one real droplet.  Retag those native (they already
                // represent ~the coarse density) so they are not re-flagged.
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    if (p_pbox[i].id() > 0
                        && active_data[i] > 0
                        && active_data[i] < flev_native
                        && mult_data[i] < sf_floor) {
                        active_data[i] = flev_native;
                    }
                });

                int n_ent_tile = amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int {
                        return (p_pbox[i].id() > 0
                                && active_data[i] > 0
                                && active_data[i] < flev_native
                                && mult_data[i] >= sf_floor) ? 1 : 0;
                    });
                if (n_ent_tile == 0) { continue; }

                Gpu::DeviceVector<int> ent_mask_d(np);
                auto* ent_mask = ent_mask_d.data();
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    ent_mask[i] = (p_pbox[i].id() > 0
                                   && active_data[i] > 0
                                   && active_data[i] < flev_native
                                   && mult_data[i] >= sf_floor) ? 1 : 0;
                });

                ParticleTileType src_tile;
                src_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                src_tile.resize(n_ent_tile);
                amrex::filterParticles(src_tile, ptile, ent_mask, 0, 0, np);

                int n_copies = n_ent_tile * (split_factor - 1);
                ParticleTileType copy_tile;
                copy_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                copy_tile.resize(n_copies);
                for (int c = 0; c < split_factor - 1; c++) {
                    amrex::copyParticles(copy_tile, src_tile, 0, c * n_ent_tile, n_ent_tile);
                }

                Long pid = ParticleType::NextID();
                ParticleType::NextID(pid + static_cast<Long>(n_copies));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    static_cast<Long>(pid + n_copies) < LastParticleID,
                    "particle id overflow in SplitMergeAtLevelBoundary split");

                auto& copy_aos = copy_tile.GetArrayOfStructs();
                auto* copy_p = copy_aos().data();
                const Long pid_start = pid;
                const int proc = my_proc;
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    copy_p[i].id()  = static_cast<int>(pid_start + i);
                    copy_p[i].cpu() = proc;
                });

                // Integer multiplicity split (see SplitParticlesForRefinement):
                // each of the split_factor pieces gets base = floor(mult/kf)
                // real droplets, with the remainder handed out one droplet each
                // to the first `rem` pieces (piece 0 = original).  Keeps every
                // piece an integer >= 1 and conserves the parent multiplicity.
                auto& copy_soa = copy_tile.GetStructOfArrays();
                auto* copy_mult = copy_soa.GetRealData(mult_idx).data();

                const int nst = n_ent_tile;
                const int kf  = split_factor;
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    long xil  = static_cast<long>(copy_mult[i] + ParticleReal(0.5));
                    long base = xil / kf;
                    long rem  = xil - base * kf;
                    int  p    = (i / nst) + 1;            // copy piece index (1..kf-1)
                    copy_mult[i] = static_cast<ParticleReal>(base + (p < rem ? 1 : 0));
                });

                auto* orig_mult   = soa.GetRealData(mult_idx).data();
                auto* orig_active = soa.GetIntData(active_idx).data();
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    if (ent_mask[i]) {
                        long xil  = static_cast<long>(orig_mult[i] + ParticleReal(0.5));
                        long base = xil / kf;
                        long rem  = xil - base * kf;
                        orig_mult[i]   = static_cast<ParticleReal>(base + (rem > 0 ? 1 : 0));
                        orig_active[i]  = flev_native;
                    }
                });
                Gpu::synchronize();

                auto old_np = ptile.numParticles();
                ptile.resize(old_np + n_copies);
                amrex::copyParticles(ptile, copy_tile, 0, old_np, n_copies);

                auto* new_active = ptile.GetStructOfArrays()
                    .GetIntData(active_idx).data();
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    new_active[old_np + i] = flev_native;
                });

                n_new_entrants += n_ent_tile;
            }

            ParallelDescriptor::ReduceLongSum(n_new_entrants);
            if (n_new_entrants > 0) {
                Redistribute(flev, flev);
                Print() << "SplitMergeAtLevelBoundary: split " << n_new_entrants
                        << " new entrants on L" << flev << "\n";
            }
        }

        // ---- Part 3: repopulate empty fine cells by cloning a neighbor ----
        {
            const auto& geom_f = m_gdb->Geom(flev);
            const auto plo_f = geom_f.ProbLoArray();
            const auto dxi_f = geom_f.InvCellSizeArray();
            const auto dx_f  = geom_f.CellSizeArray();
            const auto domain_f = geom_f.Domain();
            Long n_injected = 0;

            for (ParIterType pti(*this, flev); pti.isValid(); ++pti) {
                auto& ptile = ParticlesAt(flev, pti);
                auto& aos = ptile.GetArrayOfStructs();
                auto& soa = ptile.GetStructOfArrays();
                int np = aos.numParticles();
                if (np == 0) { continue; }

                Box box = ParticleBoxArray(flev)[pti.index()];
                auto* pstruct = aos().data();

                const IntVect bin_sz{AMREX_D_DECL(1,1,1)};
                auto binner = GetParticleBinERF{plo_f, dxi_f, domain_f, bin_sz, box};
                DenseBins<ParticleType> bins;
                bins.build(np, pstruct, numTilesInBox(box, true, bin_sz), binner);
                int num_bins = bins.numBins();
                auto* d_offsets = bins.offsetsPtr();
                auto* d_inds    = bins.permutationPtr();

                Gpu::DeviceVector<int> active_count_d(num_bins, 0);
                Gpu::DeviceVector<int> first_active_d(num_bins, -1);
                auto* active_count = active_count_d.data();
                auto* first_active = first_active_d.data();

                amrex::ParallelFor(num_bins, [=] AMREX_GPU_DEVICE (int b) {
                    auto bstart = d_offsets[b];
                    auto bstop  = d_offsets[b+1];
                    int count = 0;
                    int first = -1;
                    for (auto k = bstart; k < bstop; k++) {
                        unsigned int idx = d_inds[k];
                        if (pstruct[idx].id() > 0) {
                            count++;
                            if (first < 0) { first = static_cast<int>(idx); }
                        }
                    }
                    active_count[b] = count;
                    first_active[b] = first;
                });

                Gpu::DeviceVector<int> donor_pidx_d(num_bins, -1);
                auto* donor_pidx = donor_pidx_d.data();
                const auto blen = box.length();
                const auto blo  = box.smallEnd();
                const auto rr   = ref_ratio;

                amrex::ParallelFor(num_bins, [=] AMREX_GPU_DEVICE (int b) {
                    if (active_count[b] > 0) { return; }
                    int bx = b % blen[0];
                    int by = (b / blen[0]) % blen[1];
                    int bz = b / (blen[0] * blen[1]);

                    // Global fine-level cell index and the corresponding
                    // coarse parent index.  Cloning is restricted to face-
                    // adjacent neighbours inside the same coarse parent so
                    // SDs cannot leak out of their parent's footprint --
                    // unrestricted cloning + count-based refinement creates
                    // a positive feedback loop that grows the fine BoxArray
                    // to cover the whole domain.
                    const int gi = blo[0] + bx;
                    const int gj = blo[1] + by;
                    const int gk = blo[2] + bz;
                    const int pi = gi / rr[0];
                    const int pj = gj / rr[1];
                    const int pk = gk / rr[2];

                    const int nbr_off[6][3] = {
                        {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
                    };
                    // Prefer neighbors with >=2 active particles (steal one)
                    for (int n = 0; n < 6; n++) {
                        int nx = bx + nbr_off[n][0];
                        int ny = by + nbr_off[n][1];
                        int nz = bz + nbr_off[n][2];
                        if (nx < 0 || nx >= blen[0] || ny < 0 || ny >= blen[1] ||
                            nz < 0 || nz >= blen[2]) { continue; }
                        const int ngi = blo[0] + nx;
                        const int ngj = blo[1] + ny;
                        const int ngk = blo[2] + nz;
                        if (ngi/rr[0] != pi || ngj/rr[1] != pj || ngk/rr[2] != pk) { continue; }
                        int nb = nx + blen[0] * (ny + blen[1] * nz);
                        if (active_count[nb] >= 2 && first_active[nb] >= 0) {
                            donor_pidx[b] = first_active[nb];
                            return;
                        }
                    }
                    // Fall back: clone from any same-parent neighbor with >=1 active
                    for (int n = 0; n < 6; n++) {
                        int nx = bx + nbr_off[n][0];
                        int ny = by + nbr_off[n][1];
                        int nz = bz + nbr_off[n][2];
                        if (nx < 0 || nx >= blen[0] || ny < 0 || ny >= blen[1] ||
                            nz < 0 || nz >= blen[2]) { continue; }
                        const int ngi = blo[0] + nx;
                        const int ngj = blo[1] + ny;
                        const int ngk = blo[2] + nz;
                        if (ngi/rr[0] != pi || ngj/rr[1] != pj || ngk/rr[2] != pk) { continue; }
                        int nb = nx + blen[0] * (ny + blen[1] * nz);
                        if (first_active[nb] >= 0) {
                            donor_pidx[b] = first_active[nb];
                            return;
                        }
                    }
                });

                int n_inject = amrex::Reduce::Sum<int>(num_bins,
                    [=] AMREX_GPU_DEVICE (int b) -> int {
                        return (donor_pidx[b] >= 0) ? 1 : 0;
                    });
                if (n_inject == 0) { continue; }

                Gpu::DeviceVector<int> inject_flag_d(num_bins);
                auto* inject_flag = inject_flag_d.data();
                amrex::ParallelFor(num_bins, [=] AMREX_GPU_DEVICE (int b) {
                    inject_flag[b] = (donor_pidx[b] >= 0) ? 1 : 0;
                });

                Gpu::DeviceVector<int> inject_scan_d(num_bins);
                auto* inject_scan = inject_scan_d.data();
                Gpu::exclusive_scan(inject_flag, inject_flag + num_bins, inject_scan);

                Gpu::DeviceVector<int> compact_donor_d(n_inject);
                Gpu::DeviceVector<int> compact_bin_d(n_inject);
                auto* compact_donor = compact_donor_d.data();
                auto* compact_bin   = compact_bin_d.data();
                amrex::ParallelFor(num_bins, [=] AMREX_GPU_DEVICE (int b) {
                    if (donor_pidx[b] >= 0) {
                        int pos = inject_scan[b];
                        compact_donor[pos] = donor_pidx[b];
                        compact_bin[pos]   = b;
                    }
                });

                ParticleTileType inject_tile;
                inject_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                inject_tile.resize(n_inject);

                auto& inject_aos = inject_tile.GetArrayOfStructs();
                auto* inject_p = inject_aos().data();
                amrex::ParallelFor(n_inject, [=] AMREX_GPU_DEVICE (int i) {
                    inject_p[i] = pstruct[compact_donor[i]];
                });

                auto& inject_soa = inject_tile.GetStructOfArrays();
                for (int c = 0; c < soa.NumRealComps(); c++) {
                    auto* src = soa.GetRealData(c).data();
                    auto* dst = inject_soa.GetRealData(c).data();
                    amrex::ParallelFor(n_inject, [=] AMREX_GPU_DEVICE (int i) {
                        dst[i] = src[compact_donor[i]];
                    });
                }
                for (int c = 0; c < soa.NumIntComps(); c++) {
                    auto* src = soa.GetIntData(c).data();
                    auto* dst = inject_soa.GetIntData(c).data();
                    amrex::ParallelFor(n_inject, [=] AMREX_GPU_DEVICE (int i) {
                        dst[i] = src[compact_donor[i]];
                    });
                }

                Long pid = ParticleType::NextID();
                ParticleType::NextID(pid + static_cast<Long>(n_inject));
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    static_cast<Long>(pid + n_inject) < LastParticleID,
                    "particle id overflow during empty-cell repopulation");

                const Long pid_start = pid;
                const int proc = my_proc;
                const auto local_plo  = plo_f;
                const auto local_dx   = dx_f;
                const auto local_blen = blen;
                const auto local_blo  = box.smallEnd();
                auto* inject_active = inject_soa.GetIntData(active_idx).data();

                amrex::ParallelFor(n_inject, [=] AMREX_GPU_DEVICE (int i) {
                    inject_p[i].id()  = static_cast<int>(pid_start + i);
                    inject_p[i].cpu() = proc;

                    int b = compact_bin[i];
                    int cell_i = (b % local_blen[0]) + local_blo[0];
                    int cell_j = ((b / local_blen[0]) % local_blen[1]) + local_blo[1];
                    int cell_k = (b / (local_blen[0] * local_blen[1])) + local_blo[2];

                    inject_p[i].pos(0) = local_plo[0] + (cell_i + Real(0.5)) * local_dx[0];
                    inject_p[i].pos(1) = local_plo[1] + (cell_j + Real(0.5)) * local_dx[1];
                    inject_p[i].pos(2) = local_plo[2] + (cell_k + Real(0.5)) * local_dx[2];

                    inject_active[i] = flev_native;
                });
                Gpu::synchronize();

                auto old_np = ptile.numParticles();
                ptile.resize(old_np + n_inject);
                amrex::copyParticles(ptile, inject_tile, 0, old_np, n_inject);

                n_injected += n_inject;
            }

            ParallelDescriptor::ReduceLongSum(n_injected);
            if (n_injected > 0) {
                Redistribute(flev, flev);
                Print() << "SplitMergeAtLevelBoundary: repopulated " << n_injected
                        << " empty cells on L" << flev << "\n";
            }
        }
    }

    // ----- Per-level tag-normalizing merge sweep -----
    //
    // For every level lev in [0, finest], merge any super-droplets whose tag
    // exceeds lev+1 (departees from a finer level) into a native (tag == lev+1)
    // host in the same cell AND the same log-radius bin, so only same-size
    // droplets combine and the rain-seeding large-droplet tail is preserved.
    // A size bin with no native promotes its first departee to native.
    for (int lev = 0; lev <= finest; lev++) {
        const int lev_native = lev + 1;
        const auto& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto domain = geom.Domain();

        const int num_sp = m_num_species;
        const int num_ae = m_num_aerosols;
        const int n_rbin = (m_num_sd_per_cell >= 2) ? amrex::min(m_num_sd_per_cell, 256) : 60;
        const auto ctx = buildProcessContext(lev);
        const IntVect bin_size = {AMREX_D_DECL(1,1,1)};
        Long n_merged = 0;

        forEachParticleTileSerial(lev, ctx,
            [&](ParIterType& /*pti*/, int grid, ParticleType* pstruct_ptr,
                const SDProcess::ParticlePointers& ptrs,
                const SDProcess::ProcessContext& /*ctx_unused*/)
        {
            const size_t np = static_cast<size_t>(ptrs.num_particles);
            if (np == 0) { return; }

            int n_excess = amrex::Reduce::Sum<int>(static_cast<int>(np),
                [=] AMREX_GPU_DEVICE (int i) -> int {
                    return (pstruct_ptr[i].id() > 0
                            && ptrs.active_ptr[i] > lev_native) ? 1 : 0;
                });
            if (n_excess == 0) { return; }

            Box box = ParticleBoxArray(lev)[grid];
            int ntiles = numTilesInBox(box, true, bin_size);
            auto binner = GetParticleBinERF{plo, dxi, domain, bin_size, box};
            DenseBins<ParticleType> bins;
            bins.build(np, pstruct_ptr, ntiles, binner);
            auto inds    = bins.permutationPtr();
            auto offsets = bins.offsetsPtr();

            auto* mult_p   = ptrs.mult_ptr;
            auto* v_p0     = ptrs.v_ptr[0];
            auto* v_p1     = ptrs.v_ptr[1];
            auto* v_p2     = ptrs.v_ptr[2];
            auto* radius_p = ptrs.radius_ptr;
            auto* mass_p   = ptrs.mass_ptr;
            auto sp_mass_p = ptrs.sp_mass_ptrs;
            auto ae_mass_p = ptrs.ae_mass_ptrs;
            auto* sp_rho_p = ptrs.sp_rho_arr;
            auto* sp_sol_p = ptrs.sp_sol_arr;
            auto* ae_rho_p = ptrs.ae_rho_arr;
            auto* ae_sol_p = ptrs.ae_sol_arr;
            auto idx_w     = ctx.idx_water;
            auto rho_w     = ctx.rho_water;
            auto* active_int_p = ptrs.active_ptr;
            const int nat_tag = lev_native;

            int num_bins = bins.numBins();

            auto tile_merged = amrex::Reduce::Sum<Long>(num_bins,
                [=] AMREX_GPU_DEVICE (int i_bin) -> Long
            {
                auto bin_start = offsets[i_bin];
                auto bin_stop  = offsets[i_bin+1];
                int np_bin = static_cast<int>(bin_stop - bin_start);
                if (np_bin < 2) {
                    if (np_bin == 1) {
                        unsigned int idx = inds[bin_start];
                        if (pstruct_ptr[idx].id() > 0 && active_int_p[idx] > nat_tag) {
                            active_int_p[idx] = nat_tag;
                        }
                    }
                    return 0;
                }

                // Per-cell adaptive radius range over live SDs.
                ParticleReal rmin_c = ParticleReal(-1);
                ParticleReal rmax_c = ParticleReal(-1);
                for (int k = 0; k < np_bin; k++) {
                    unsigned int idx = inds[bin_start + k];
                    if (pstruct_ptr[idx].id() <= 0) { continue; }
                    ParticleReal r = radius_p[idx];
                    if (r > ParticleReal(0)) {
                        if (rmin_c < ParticleReal(0) || r < rmin_c) { rmin_c = r; }
                        if (r > rmax_c) { rmax_c = r; }
                    }
                }
                if (rmax_c <= ParticleReal(0)) { return 0; }
                const ParticleReal lnrmin  = std::log(rmin_c);
                const ParticleReal lnrange = std::log(rmax_c) - lnrmin;
                const ParticleReal inv_lnr = (lnrange > ParticleReal(1.0e-12))
                                           ? (static_cast<ParticleReal>(n_rbin) / lnrange)
                                           : ParticleReal(0);

                Long bin_merged = 0;

                // Absorb departees (tag > native) into a same-size native host;
                // if a size bin has no native, promote its first departee.
                // Natives are never merged with one another, so the boundary
                // cell's resolved spectrum is preserved.
                for (int rb = 0; rb < n_rbin; rb++) {
                    int host = -1;
                    for (int k = 0; k < np_bin; k++) {
                        unsigned int idx = inds[bin_start + k];
                        if (pstruct_ptr[idx].id() <= 0 || active_int_p[idx] != nat_tag) { continue; }
                        if (sdmRadiusBin(radius_p[idx], lnrmin, inv_lnr, n_rbin) == rb) {
                            host = static_cast<int>(idx); break;
                        }
                    }
                    if (host < 0) {
                        for (int k = 0; k < np_bin; k++) {
                            unsigned int idx = inds[bin_start + k];
                            if (pstruct_ptr[idx].id() <= 0 || active_int_p[idx] <= nat_tag) { continue; }
                            if (sdmRadiusBin(radius_p[idx], lnrmin, inv_lnr, n_rbin) == rb) {
                                active_int_p[idx] = nat_tag;
                                host = static_cast<int>(idx); break;
                            }
                        }
                    }
                    if (host < 0) { continue; }

                    unsigned int h = static_cast<unsigned int>(host);
                    for (int k = 0; k < np_bin; k++) {
                        unsigned int idx = inds[bin_start + k];
                        if (static_cast<int>(idx) == host) { continue; }
                        if (pstruct_ptr[idx].id() <= 0 || active_int_p[idx] <= nat_tag) { continue; }
                        if (sdmRadiusBin(radius_p[idx], lnrmin, inv_lnr, n_rbin) != rb) { continue; }

                        ParticleReal xi_e = mult_p[idx];
                        ParticleReal xi_s = mult_p[h];
                        ParticleReal xi_new = xi_e + xi_s;
                        if (xi_new <= ParticleReal(0)) { continue; }
                        ParticleReal inv_xi = ParticleReal(1.0) / xi_new;

                        mult_p[h] = xi_new;
                        v_p0[h] = (xi_e * v_p0[idx] + xi_s * v_p0[h]) * inv_xi;
                        v_p1[h] = (xi_e * v_p1[idx] + xi_s * v_p1[h]) * inv_xi;
                        v_p2[h] = (xi_e * v_p2[idx] + xi_s * v_p2[h]) * inv_xi;
                        for (int s = 0; s < num_sp; s++) {
                            sp_mass_p[s][h] = (xi_e * sp_mass_p[s][idx]
                                             + xi_s * sp_mass_p[s][h]) * inv_xi;
                        }
                        for (int a = 0; a < num_ae; a++) {
                            ae_mass_p[a][h] = (xi_e * ae_mass_p[a][idx]
                                             + xi_s * ae_mass_p[a][h]) * inv_xi;
                        }
                        pstruct_ptr[idx].id() = -1;
                        bin_merged++;
                    }
                }

                for (int k = 0; k < np_bin; k++) {
                    unsigned int idx = inds[bin_start + k];
                    if (pstruct_ptr[idx].id() > 0 && active_int_p[idx] == nat_tag) {
                        updateParticleAttributes(
                            idx, radius_p, mass_p,
                            idx_w, rho_w, num_sp, num_ae,
                            sp_sol_p, ae_sol_p,
                            sp_mass_p, ae_mass_p,
                            sp_rho_p, ae_rho_p);
                    }
                }
                return bin_merged;
            });
            Gpu::synchronize();
            n_merged += tile_merged;
        });

        ParallelDescriptor::ReduceLongSum(n_merged);
        if (n_merged > 0) {
            Redistribute(lev, lev);
            Print() << "SplitMergeAtLevelBoundary: merged " << n_merged
                    << " excess-tag particles on L" << lev << "\n";
        }
    }
}

#endif
