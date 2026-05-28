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
    m_ae_density.resize(num_ae);
    m_ae_solubility.resize(num_ae);
    m_ae_ionization.resize(num_ae);
    m_ae_mol_weight.resize(num_ae);

    amrex::Vector<ParticleReal> sp_density_h(num_sp), sp_ionization_h(num_sp), sp_mol_weight_h(num_sp);
    amrex::Vector<int> sp_solubility_h(num_sp);
    amrex::Vector<ParticleReal> ae_density_h(num_ae), ae_ionization_h(num_ae), ae_mol_weight_h(num_ae);
    amrex::Vector<int> ae_solubility_h(num_ae);

    for (int i = 0; i < num_sp; i++) {
        sp_density_h[i] = m_species_mat[i]->m_density;
        sp_solubility_h[i] = static_cast<int>(m_species_mat[i]->m_is_soluble);
        sp_ionization_h[i] = m_species_mat[i]->m_ionization;
        sp_mol_weight_h[i] = m_species_mat[i]->m_mol_weight;
    }
    for (int i = 0; i < num_ae; i++) {
        ae_density_h[i] = m_aerosol_mat[i]->m_density;
        ae_solubility_h[i] = static_cast<int>(m_aerosol_mat[i]->m_is_soluble);
        ae_ionization_h[i] = m_aerosol_mat[i]->m_ionization;
        ae_mol_weight_h[i] = m_aerosol_mat[i]->m_mol_weight;
    }

    Gpu::copy(Gpu::hostToDevice, sp_density_h.begin(), sp_density_h.end(), m_sp_density.begin());
    Gpu::copy(Gpu::hostToDevice, sp_solubility_h.begin(), sp_solubility_h.end(), m_sp_solubility.begin());
    Gpu::copy(Gpu::hostToDevice, sp_ionization_h.begin(), sp_ionization_h.end(), m_sp_ionization.begin());
    Gpu::copy(Gpu::hostToDevice, sp_mol_weight_h.begin(), sp_mol_weight_h.end(), m_sp_mol_weight.begin());
    Gpu::copy(Gpu::hostToDevice, ae_density_h.begin(), ae_density_h.end(), m_ae_density.begin());
    Gpu::copy(Gpu::hostToDevice, ae_solubility_h.begin(), ae_solubility_h.end(), m_ae_solubility.begin());
    Gpu::copy(Gpu::hostToDevice, ae_ionization_h.begin(), ae_ionization_h.end(), m_ae_ionization.begin());
    Gpu::copy(Gpu::hostToDevice, ae_mol_weight_h.begin(), ae_mol_weight_h.end(), m_ae_mol_weight.begin());

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

/*! \brief Split coarse-level particles in newly-refined cells.  Driven by
 *  per-source-level "finest covering level" masks: for source level `clev`,
 *  every cell maps to the highest level k in [clev, finest_level] whose
 *  ParticleBoxArray covers it.  Particles on `clev` in cells with k > clev
 *  are split by the cumulative refinement factor product_{j=clev}^{k-1}
 *  refRatio(j) and tagged with `k + 1` (the destination level's native tag).
 *
 *  This single call handles cascading multi-level refinement (e.g. a fresh
 *  regrid creating both L1 and L2 from L0): coarse particles under L2 are
 *  split by R(0->2) = refRatio(0) * refRatio(1) in one shot, so the next
 *  Redistribute lands the right SD density on every level.
 *
 *  Only particles whose multiplicity is >= the relevant split factor are
 *  split; lower-multiplicity particles are skipped (the count is reported). */
void SuperDropletPC::SplitParticlesForRefinement (int finest_level)
{
    BL_PROFILE("SuperDropletPC::SplitParticlesForRefinement()");

    if (!m_split_merge_amr) { return; }
    if (finest_level < 1) { return; }

    // Resync the particle container's internal dummy MultiFabs with the
    // current mesh BoxArrays/DistributionMaps.  Required when this function
    // is called immediately after AMR regrid: the ParConstIter and the
    // iMultiFab `finest_cov` below must agree on level layout, otherwise the
    // FabArray::operator[](MFIter) assert fires (LocalIndex mismatch).
    for (int lev = 0; lev <= finest_level; lev++) {
        RedefineDummyMF(lev);
    }

    constexpr int rt_off_r = SuperDropletsRealIdx::ncomps;
    constexpr int rt_off_i = SuperDropletsIntIdx::ncomps;
    const int mult_idx   = rt_off_r + SuperDropletsRealIdxSoA_RT::multiplicity;
    const int active_idx = rt_off_i + SuperDropletsIntIdxSoA_RT::active;
    const int my_proc = ParallelDescriptor::MyProc();

    for (int clev = 0; clev < finest_level; clev++) {

        const auto& cba = ParticleBoxArray(clev);
        const auto& cdm = ParticleDistributionMap(clev);
        const auto& geom = m_gdb->Geom(clev);
        const auto plo    = geom.ProbLoArray();
        const auto dxi    = geom.InvCellSizeArray();
        const auto dx     = geom.CellSizeArray();
        const auto domain = geom.Domain();

        // Per-cell finest covering level on the clev grid: for cell (i,j,k),
        // the highest level m in [clev, finest_level] whose ParticleBoxArray
        // covers it.  Iterates m = clev+1 .. finest_level, OR-ing the per-m
        // mask into finest_cov via a max kernel.  Proper nesting (L_m always
        // covered by L_{m-1}) ensures the deepest covering level wins.
        iMultiFab finest_cov(cba, cdm, 1, 0);
        finest_cov.setVal(clev);
        {
            IntVect cum_ref = IntVect::TheUnitVector();
            for (int m = clev+1; m <= finest_level; m++) {
                cum_ref *= m_gdb->refRatio(m-1);
                iMultiFab cov_m = makeFineMask(cba, cdm,
                                               ParticleBoxArray(m), cum_ref,
                                               /*crse_value=*/clev,
                                               /*fine_value=*/m);
                for (MFIter mfi(finest_cov); mfi.isValid(); ++mfi) {
                    auto dst_arr = finest_cov.array(mfi);
                    auto src_arr = cov_m.const_array(mfi);
                    const Box& bx = mfi.validbox();
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) {
                        dst_arr(i,j,kk) = amrex::max(dst_arr(i,j,kk), src_arr(i,j,kk));
                    });
                }
            }
        }

        // One split pass per destination level k > clev.  Each pass uses the
        // cumulative refinement factor from clev down to k and selects only
        // those clev cells whose finest covering level is exactly k.
        IntVect cum_ref_to_k = IntVect::TheUnitVector();
        for (int k = clev+1; k <= finest_level; k++) {
            cum_ref_to_k *= m_gdb->refRatio(k-1);
            const int split_factor = AMREX_D_TERM(cum_ref_to_k[0],
                                                  *cum_ref_to_k[1],
                                                  *cum_ref_to_k[2]);
            if (split_factor <= 1) { continue; }
            const ParticleReal inv_split = ParticleReal(1.0)
                / static_cast<ParticleReal>(split_factor);
            const int daughter_tag = k + 1;
            const int sf_count = split_factor;
            const int target_lev = k;

            Long n_to_split         = 0;
            Long n_skipped_low_mult = 0;

            for (ParConstIterType pti(*this, clev); pti.isValid(); ++pti) {
                const auto& ptile = ParticlesAt(clev, pti);
                const auto& aos   = ptile.GetArrayOfStructs();
                int np = aos.numParticles();
                if (np == 0) { continue; }

                auto fc_arr        = finest_cov[pti].const_array();
                const auto* p_pbox = aos().data();
                const auto& soa    = ptile.GetStructOfArrays();
                const auto* mult_ptr = soa.GetRealData(mult_idx).data();

                n_to_split += static_cast<Long>(amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int {
                        const auto& p = p_pbox[i];
                        if (p.id() <= 0) { return 0; }
                        IntVect iv = amrex::getParticleCell(p, plo, dxi, domain);
                        return (fc_arr(iv) == target_lev
                                && mult_ptr[i] >= static_cast<ParticleReal>(sf_count)) ? 1 : 0;
                    }));
                n_skipped_low_mult += static_cast<Long>(amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int {
                        const auto& p = p_pbox[i];
                        if (p.id() <= 0) { return 0; }
                        IntVect iv = amrex::getParticleCell(p, plo, dxi, domain);
                        return (fc_arr(iv) == target_lev
                                && mult_ptr[i] < static_cast<ParticleReal>(sf_count)) ? 1 : 0;
                    }));
            }
            ParallelDescriptor::ReduceLongSum(n_to_split);
            ParallelDescriptor::ReduceLongSum(n_skipped_low_mult);

            Long n_new_particles = n_to_split * (static_cast<Long>(split_factor) - 1);
            if (n_new_particles == 0) {
                Print() << "SplitParticlesForRefinement: L" << clev
                        << " -> L" << k << " -- no particles to split (factor "
                        << split_factor << ")\n";
                if (n_skipped_low_mult > 0) {
                    Print() << "  (skipped " << n_skipped_low_mult
                            << " particles with multiplicity < "
                            << split_factor << ")\n";
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
                            << clev << " -> L" << k << ".  Need "
                            << bytes_needed/(1024*1024) << " MB, available "
                            << free_mem/(1024*1024) << " MB.  Skipping.\n";
                    continue;
                }
#endif
                amrex::ignore_unused(bytes_needed);
            }

            const int sf = split_factor;

            for (ParIterType pti(*this, clev); pti.isValid(); ++pti) {
                auto& ptile = ParticlesAt(clev, pti);
                auto& aos   = ptile.GetArrayOfStructs();
                int np = aos.numParticles();
                if (np == 0) { continue; }

                auto fc_arr = finest_cov[pti].const_array();

                Gpu::DeviceVector<int> split_mask_d(np);
                auto* mask_d_ptr = split_mask_d.data();
                auto* p_pbox = aos().data();
                auto& soa = ptile.GetStructOfArrays();
                const auto* mult_data = soa.GetRealData(mult_idx).data();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    const auto& p = p_pbox[i];
                    if (p.id() <= 0) { mask_d_ptr[i] = 0; return; }
                    IntVect iv = amrex::getParticleCell(p, plo, dxi, domain);
                    mask_d_ptr[i] = (fc_arr(iv) == target_lev
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

                Gpu::DeviceVector<ParticleReal> copy_weight_d(n_copies, ParticleReal(0.0));
                auto* copy_weight = copy_weight_d.data();
                Gpu::DeviceVector<ParticleReal> weight_sum_d(n_split_tile, ParticleReal(0.0));
                auto* weight_sum = weight_sum_d.data();

                const ParticleReal alpha = ParticleReal(0.4);
                const int nst = n_split_tile;
                const int n_copies_per = split_factor - 1;
                const auto local_dx  = dx;
                const auto local_plo = plo;
                const auto local_dxi = dxi;

                amrex::ParallelForRNG(n_copies,
                    [=] AMREX_GPU_DEVICE (int i, const amrex::RandomEngine& rng) {
                    ParticleReal w = ParticleReal(1.0)
                                   + alpha * (amrex::Random(rng) - ParticleReal(0.5));
                    copy_weight[i] = w;
                    int j = i % nst;
                    Gpu::Atomic::AddNoRet(&weight_sum[j], w);

                    for (int d = 0; d < AMREX_SPACEDIM; d++) {
                        ParticleReal pert = local_dx[d] * ParticleReal(0.01)
                            * (ParticleReal(2.0) * amrex::Random(rng) - ParticleReal(1.0));
                        ParticleReal pos_new = copy_p[i].pos(d) + pert;
                        int cell_idx = static_cast<int>(amrex::Math::floor(
                            (copy_p[i].pos(d) - local_plo[d]) * local_dxi[d]));
                        ParticleReal cell_lo = local_plo[d] + cell_idx * local_dx[d];
                        ParticleReal cell_hi = cell_lo + local_dx[d];
                        constexpr ParticleReal eps = ParticleReal(1.0e-10);
                        pos_new = amrex::max(pos_new, cell_lo + eps);
                        pos_new = amrex::min(pos_new, cell_hi - eps);
                        copy_p[i].pos(d) = pos_new;
                    }
                });

                const ParticleReal copy_share
                    = ParticleReal(n_copies_per) / ParticleReal(split_factor);
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    int j = i % nst;
                    ParticleReal orig = copy_mult[i];
                    copy_mult[i] = orig * copy_share * copy_weight[i] / weight_sum[j];
                });

                Gpu::DeviceVector<int> prefix_d(np);
                amrex::Gpu::exclusive_scan(mask_d_ptr, mask_d_ptr + np, prefix_d.data());

                auto* orig_mult   = soa.GetRealData(mult_idx).data();
                auto* orig_active = soa.GetIntData(active_idx).data();
                const int dtag = daughter_tag;
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    if (mask_d_ptr[i]) {
                        orig_mult[i]   *= inv_split;
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

            Print() << "SplitParticlesForRefinement: L" << clev
                    << " -> L" << k << " split " << n_to_split
                    << " particles by factor " << split_factor
                    << " (created " << n_new_particles << " new SDs)\n";
            if (n_skipped_low_mult > 0) {
                Print() << "  (skipped " << n_skipped_low_mult
                        << " particles with multiplicity < "
                        << split_factor << ")\n";
            }
        }
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
 *                tag them native to flev.
 *        Part 3: repopulate empty fine cells on flev by cloning a face-
 *                adjacent neighbor.
 *    - Per-level (lev = 0..finest):
 *        Part 2: in each cell on `lev`, merge any super-droplets with tag
 *                > lev+1 (leftovers from a finer level that has since drifted
 *                away or vanished) into a native (tag == lev+1) host via
 *                multiplicity-weighted host cycling.  If no native host is
 *                present, pick the first surviving SD as the host, merge all
 *                others in, and retag it native.  This single sweep handles
 *                per-step departees, regrid-time leftovers from disappeared
 *                finer levels, and cascading de-refinement in one pass.
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
        const ParticleReal inv_split = ParticleReal(1.0) / static_cast<ParticleReal>(split_factor);

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
                const auto* active_data = soa.GetIntData(active_idx).data();

                int n_ent_tile = amrex::Reduce::Sum<int>(np,
                    [=] AMREX_GPU_DEVICE (int i) -> int {
                        return (p_pbox[i].id() > 0
                                && active_data[i] > 0
                                && active_data[i] < flev_native) ? 1 : 0;
                    });
                if (n_ent_tile == 0) { continue; }

                Gpu::DeviceVector<int> ent_mask_d(np);
                auto* ent_mask = ent_mask_d.data();
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    ent_mask[i] = (p_pbox[i].id() > 0
                                   && active_data[i] > 0
                                   && active_data[i] < flev_native) ? 1 : 0;
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

                // Multiplicity: original keeps 1/split_factor; copies share
                // the remaining (split_factor-1)/split_factor, distributed by
                // randomized weights (alpha = 0.4) normalized per source.
                auto& copy_soa = copy_tile.GetStructOfArrays();
                auto* copy_mult = copy_soa.GetRealData(mult_idx).data();

                Gpu::DeviceVector<ParticleReal> copy_weight_d(n_copies, ParticleReal(0.0));
                auto* copy_weight = copy_weight_d.data();
                Gpu::DeviceVector<ParticleReal> weight_sum_d(n_ent_tile, ParticleReal(0.0));
                auto* weight_sum = weight_sum_d.data();

                const ParticleReal alpha = ParticleReal(0.4);
                const int nst = n_ent_tile;
                amrex::ParallelForRNG(n_copies,
                    [=] AMREX_GPU_DEVICE (int i, const amrex::RandomEngine& rng) {
                    ParticleReal w = ParticleReal(1.0)
                                   + alpha * (amrex::Random(rng) - ParticleReal(0.5));
                    copy_weight[i] = w;
                    int j = i % nst;
                    Gpu::Atomic::AddNoRet(&weight_sum[j], w);
                });

                const int n_copies_per = split_factor - 1;
                const ParticleReal copy_share
                    = ParticleReal(n_copies_per) / ParticleReal(split_factor);
                amrex::ParallelFor(n_copies, [=] AMREX_GPU_DEVICE (int i) {
                    int j = i % nst;
                    ParticleReal orig = copy_mult[i];
                    copy_mult[i] = orig * copy_share * copy_weight[i] / weight_sum[j];
                });

                // Reduce original new-entrant multiplicity and tag with this
                // level's native value (flev + 1).
                auto* orig_mult   = soa.GetRealData(mult_idx).data();
                auto* orig_active = soa.GetIntData(active_idx).data();
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) {
                    if (ent_mask[i]) {
                        orig_mult[i]   *= inv_split;
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

                amrex::ParallelFor(num_bins, [=] AMREX_GPU_DEVICE (int b) {
                    if (active_count[b] > 0) { return; }
                    int bx = b % blen[0];
                    int by = (b / blen[0]) % blen[1];
                    int bz = b / (blen[0] * blen[1]);

                    const int fo[6][3] = {
                        {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
                    };
                    // Prefer neighbors with >=2 active particles (steal one)
                    for (int n = 0; n < 6; n++) {
                        int nx = bx + fo[n][0];
                        int ny = by + fo[n][1];
                        int nz = bz + fo[n][2];
                        if (nx < 0 || nx >= blen[0] || ny < 0 || ny >= blen[1] ||
                            nz < 0 || nz >= blen[2]) { continue; }
                        int nb = nx + blen[0] * (ny + blen[1] * nz);
                        if (active_count[nb] >= 2 && first_active[nb] >= 0) {
                            donor_pidx[b] = first_active[nb];
                            return;
                        }
                    }
                    // Fall back: clone from any neighbor with >=1 active
                    for (int n = 0; n < 6; n++) {
                        int nx = bx + fo[n][0];
                        int ny = by + fo[n][1];
                        int nz = bz + fo[n][2];
                        if (nx < 0 || nx >= blen[0] || ny < 0 || ny >= blen[1] ||
                            nz < 0 || nz >= blen[2]) { continue; }
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
    // exceeds lev+1 (leftovers from a finer level that has since drifted away
    // or vanished entirely) into a native (tag == lev+1) host within the same
    // cell.  If no native host is present in a cell, the first surviving SD
    // is promoted to host: all other SDs in that cell are merged into it via
    // multiplicity-weighted averaging and the host is retagged native.  This
    // single sweep handles per-step departees, regrid-time leftovers from
    // disappeared finer levels, and cascading de-refinement (e.g. simultaneous
    // loss of L1 and L2) in one pass.
    for (int lev = 0; lev <= finest; lev++) {
        const int lev_native = lev + 1;
        const auto& geom = m_gdb->Geom(lev);
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();
        const auto domain = geom.Domain();

        const int num_sp = m_num_species;
        const int num_ae = m_num_aerosols;
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

                int n_hosts = 0;
                for (int k = 0; k < np_bin; k++) {
                    unsigned int idx = inds[bin_start + k];
                    if (pstruct_ptr[idx].id() > 0 && active_int_p[idx] == nat_tag) {
                        n_hosts++;
                    }
                }

                Long bin_merged = 0;
                if (n_hosts > 0) {
                    int host_cycle = 0;
                    for (int k = 0; k < np_bin; k++) {
                        unsigned int idx = inds[bin_start + k];
                        if (pstruct_ptr[idx].id() <= 0 || active_int_p[idx] <= nat_tag) { continue; }

                        int target_h = host_cycle % n_hosts;
                        unsigned int h_idx = 0;
                        int h_count = 0;
                        for (int h = 0; h < np_bin; h++) {
                            unsigned int candidate = inds[bin_start + h];
                            if (pstruct_ptr[candidate].id() > 0
                                && active_int_p[candidate] == nat_tag) {
                                if (h_count == target_h) { h_idx = candidate; break; }
                                h_count++;
                            }
                        }
                        host_cycle++;

                        ParticleReal xi_e = mult_p[idx];
                        ParticleReal xi_s = mult_p[h_idx];
                        ParticleReal xi_new = xi_e + xi_s;
                        if (xi_new <= ParticleReal(0)) { continue; }
                        ParticleReal inv_xi = ParticleReal(1.0) / xi_new;

                        mult_p[h_idx] = xi_new;
                        v_p0[h_idx] = (xi_e * v_p0[idx] + xi_s * v_p0[h_idx]) * inv_xi;
                        v_p1[h_idx] = (xi_e * v_p1[idx] + xi_s * v_p1[h_idx]) * inv_xi;
                        v_p2[h_idx] = (xi_e * v_p2[idx] + xi_s * v_p2[h_idx]) * inv_xi;

                        for (int s = 0; s < num_sp; s++) {
                            sp_mass_p[s][h_idx] = (xi_e * sp_mass_p[s][idx]
                                                 + xi_s * sp_mass_p[s][h_idx]) * inv_xi;
                        }
                        for (int a = 0; a < num_ae; a++) {
                            ae_mass_p[a][h_idx] = (xi_e * ae_mass_p[a][idx]
                                                 + xi_s * ae_mass_p[a][h_idx]) * inv_xi;
                        }

                        pstruct_ptr[idx].id() = -1;
                        bin_merged++;
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
                } else {
                    // No-host branch: every SD in this cell has tag > nat_tag
                    // (cascading de-refinement leftovers).  Pick the first
                    // surviving SD as the host, merge all others into it, and
                    // retag the host as nat_tag.
                    int host_k = -1;
                    for (int k = 0; k < np_bin; k++) {
                        unsigned int idx = inds[bin_start + k];
                        if (pstruct_ptr[idx].id() > 0) { host_k = k; break; }
                    }
                    if (host_k < 0) { return 0; }
                    unsigned int h_idx = inds[bin_start + host_k];

                    for (int k = host_k + 1; k < np_bin; k++) {
                        unsigned int idx = inds[bin_start + k];
                        if (pstruct_ptr[idx].id() <= 0) { continue; }

                        ParticleReal xi_e = mult_p[idx];
                        ParticleReal xi_s = mult_p[h_idx];
                        ParticleReal xi_new = xi_e + xi_s;
                        if (xi_new <= ParticleReal(0)) { continue; }
                        ParticleReal inv_xi = ParticleReal(1.0) / xi_new;

                        mult_p[h_idx] = xi_new;
                        v_p0[h_idx] = (xi_e * v_p0[idx] + xi_s * v_p0[h_idx]) * inv_xi;
                        v_p1[h_idx] = (xi_e * v_p1[idx] + xi_s * v_p1[h_idx]) * inv_xi;
                        v_p2[h_idx] = (xi_e * v_p2[idx] + xi_s * v_p2[h_idx]) * inv_xi;

                        for (int s = 0; s < num_sp; s++) {
                            sp_mass_p[s][h_idx] = (xi_e * sp_mass_p[s][idx]
                                                 + xi_s * sp_mass_p[s][h_idx]) * inv_xi;
                        }
                        for (int a = 0; a < num_ae; a++) {
                            ae_mass_p[a][h_idx] = (xi_e * ae_mass_p[a][idx]
                                                 + xi_s * ae_mass_p[a][h_idx]) * inv_xi;
                        }

                        pstruct_ptr[idx].id() = -1;
                        bin_merged++;
                    }

                    active_int_p[h_idx] = nat_tag;
                    updateParticleAttributes(
                        h_idx, radius_p, mass_p,
                        idx_w, rho_w, num_sp, num_ae,
                        sp_sol_p, ae_sol_p,
                        sp_mass_p, ae_mass_p,
                        sp_rho_p, ae_rho_p);
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
