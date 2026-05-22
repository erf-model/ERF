#include <stdio.h>
#include <limits>
#include <fstream>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include "ERF_Constants.H"
#include "ERF_SuperDropletPC.H"

#ifdef ERF_USE_PARTICLES

using namespace amrex;
using namespace SDPCDefn;
using SDTDType = typename SuperDropletPC::ParticleTileType::ConstParticleTileDataType;

/*! Get real-type particle attribute names */
Vector<std::string> SuperDropletPC::varNames () const
{
    BL_PROFILE("ERFPCPC::varNames()");
    amrex::Vector<std::string> retval = {   AMREX_D_DECL("xvel","yvel","zvel"),
                                            "particle_mass",
                                            "temperature",
                                            "radius",
                                            "multiplicity",
                                            "terminal_velocity",
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                                            "condensation_tendency",
#endif
                                            "uid" };
    for (int i = 0; i < m_num_aerosols; i++) {
        retval.push_back(std::string("aerosol_mass_"+getEnumNameString(m_aerosol_mat[i]->m_name)));
    }
    for (int i = 0; i < m_num_species; i++) {
        retval.push_back(std::string("species_mass_"+getEnumNameString(m_species_mat[i]->m_name)));
    }
    return retval;
}

/*! Get Eulerian plot variable names */
Vector<std::string> SuperDropletPC::meshPlotVarNames () const
{
    BL_PROFILE("ERFPCPC::varNames()");
    amrex::Vector<std::string> retval = { AMREX_D_DECL("mass_flux_x",
                                                       "mass_flux_y",
                                                       "mass_flux_z"),
                                          "number_density",
                                          "sd_number_density",
                                          "mass_density",
                                          "radius",
                                          ("mass_density_"+getEnumNameString(m_species_mat[m_idx_w]->m_name)),
                                          AMREX_D_DECL (
                                              ("mass_flux_x_"+getEnumNameString(m_species_mat[m_idx_w]->m_name)),
                                              ("mass_flux_y_"+getEnumNameString(m_species_mat[m_idx_w]->m_name)),
                                              ("mass_flux_z_"+getEnumNameString(m_species_mat[m_idx_w]->m_name)) ) };
    for (int i = 0; i < m_num_aerosols; i++) {
        retval.push_back(std::string("aerosol_mass_density_"+getEnumNameString(m_aerosol_mat[i]->m_name)));
        retval.push_back(std::string("aerosol_mass_flux_x_"+getEnumNameString(m_aerosol_mat[i]->m_name)));
        retval.push_back(std::string("aerosol_mass_flux_y_"+getEnumNameString(m_aerosol_mat[i]->m_name)));
        retval.push_back(std::string("aerosol_mass_flux_z_"+getEnumNameString(m_aerosol_mat[i]->m_name)));
    }
    for (int i = 0; i < m_num_species; i++) {
        retval.push_back(std::string("species_mass_density_"+getEnumNameString(m_species_mat[i]->m_name)));
        retval.push_back(std::string("species_mass_flux_x_" +getEnumNameString(m_species_mat[i]->m_name)));
        retval.push_back(std::string("species_mass_flux_y_" +getEnumNameString(m_species_mat[i]->m_name)));
        retval.push_back(std::string("species_mass_flux_z_" +getEnumNameString(m_species_mat[i]->m_name)));
    }
    return retval;
}

/*! Compute diagnostics (max, min, avg radius, mass, etc) */
void SuperDropletPC::Diagnostics( const int a_iter,
                                  const int a_lev,
                                  const Real a_time,
                                  const bool a_flag )
{
    BL_PROFILE("SuperDropletPC::Diagnostics()");

    auto num_total_particles = TotalNumberOfParticles();
    auto num_superdroplets = static_cast<Real>(NumSuperDroplets());

    // Number of base attributes to reduce (radius, multiplicity, mass, vx, vy, vz, term_vel)
    constexpr int NUM_BASE_ATTRS = 7;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    constexpr int NUM_ATTRS = NUM_BASE_ATTRS + 1;  // + cond_tendency
#else
    constexpr int NUM_ATTRS = NUM_BASE_ATTRS;
#endif

    // Fused reduction: compute min, max, and weighted sum for all attributes in a single pass
    using ReduceTuple = GpuTuple<Real, Real, Real, Real, Real, Real, Real,  // min values (7)
                                 Real, Real, Real, Real, Real, Real, Real,  // max values (7)
                                 Real, Real, Real, Real, Real, Real, Real,  // sum values (7)
                                 Real                                        // multiplicity sum
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                                 , Real, Real, Real  // min, max, sum for cond_tendency
#endif
                                 >;

    ReduceOps<ReduceOpMin, ReduceOpMin, ReduceOpMin, ReduceOpMin, ReduceOpMin, ReduceOpMin, ReduceOpMin,
              ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax, ReduceOpMax,
              ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,
              ReduceOpSum
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
              , ReduceOpMin, ReduceOpMax, ReduceOpSum
#endif
              > reduce_ops;

    auto r = ReduceData<Real, Real, Real, Real, Real, Real, Real,
                        Real, Real, Real, Real, Real, Real, Real,
                        Real, Real, Real, Real, Real, Real, Real,
                        Real
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                        , Real, Real, Real
#endif
                        >(reduce_ops);

    for (int lev = 0; lev <= finestLevel(); ++lev) {
        for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti) {
            const auto& ptd = pti.GetParticleTile().getConstParticleTileData();
            const int np = pti.numParticles();

            reduce_ops.eval(np, r,
                [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                {
                    const Real n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                    const Real radius = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                    const Real mass = ptd.m_rdata[SuperDropletsRealIdx::mass][i];
                    const Real vx = ptd.m_rdata[SuperDropletsRealIdx::vx][i];
                    const Real vy = ptd.m_rdata[SuperDropletsRealIdx::vy][i];
                    const Real vz = ptd.m_rdata[SuperDropletsRealIdx::vz][i];
                    const Real term_vel = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                    const Real cond_t = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::cond_tendency][i];
#endif

                    return {
                        // Min values
                        radius, n, mass, vx, vy, vz, term_vel,
                        // Max values
                        radius, n, mass, vx, vy, vz, term_vel,
                        // Weighted sums (n * value)
                        n * radius, n, n * mass, n * vx, n * vy, n * vz, n * term_vel,
                        // Multiplicity sum
                        n
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                        , cond_t, cond_t, n * cond_t  // min, max, sum for cond_tendency
#endif
                    };
                });
        }
    }

    auto hv = r.value(reduce_ops);

    // Extract results from tuple
    Real min_par_radius = amrex::get<0>(hv);
    Real min_multiplic  = amrex::get<1>(hv);
    Real min_par_mass   = amrex::get<2>(hv);
    Real min_par_vx     = amrex::get<3>(hv);
    Real min_par_vy     = amrex::get<4>(hv);
    Real min_par_vz     = amrex::get<5>(hv);
    Real min_term_v     = amrex::get<6>(hv);

    Real max_par_radius = amrex::get<7>(hv);
    Real max_multiplic  = amrex::get<8>(hv);
    Real max_par_mass   = amrex::get<9>(hv);
    Real max_par_vx     = amrex::get<10>(hv);
    Real max_par_vy     = amrex::get<11>(hv);
    Real max_par_vz     = amrex::get<12>(hv);
    Real max_term_v     = amrex::get<13>(hv);

    Real avg_par_radius = amrex::get<14>(hv);
    Real avg_multiplic  = amrex::get<15>(hv);
    Real avg_par_mass   = amrex::get<16>(hv);
    Real avg_par_vx     = amrex::get<17>(hv);
    Real avg_par_vy     = amrex::get<18>(hv);
    Real avg_par_vz     = amrex::get<19>(hv);
    Real avg_term_v     = amrex::get<20>(hv);

    // Note: amrex::get<21>(hv) is the multiplicity sum, used for avg_multiplic normalization

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    Real min_cond_t = amrex::get<22>(hv);
    Real max_cond_t = amrex::get<23>(hv);
    Real avg_cond_t = amrex::get<24>(hv);
#endif

    // MPI reductions - batch into arrays for efficiency
    Real min_vals[NUM_ATTRS] = {min_par_radius, min_multiplic, min_par_mass,
                                min_par_vx, min_par_vy, min_par_vz, min_term_v
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                                , min_cond_t
#endif
                               };
    Real max_vals[NUM_ATTRS] = {max_par_radius, max_multiplic, max_par_mass,
                                max_par_vx, max_par_vy, max_par_vz, max_term_v
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                                , max_cond_t
#endif
                               };
    Real sum_vals[NUM_ATTRS] = {avg_par_radius, avg_multiplic, avg_par_mass,
                                avg_par_vx, avg_par_vy, avg_par_vz, avg_term_v
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                                , avg_cond_t
#endif
                               };

    ParallelDescriptor::ReduceRealMin(min_vals, NUM_ATTRS);
    ParallelDescriptor::ReduceRealMax(max_vals, NUM_ATTRS);
    ParallelDescriptor::ReduceRealSum(sum_vals, NUM_ATTRS);

    // Unpack back to named variables
    min_par_radius = min_vals[0]; min_multiplic = min_vals[1]; min_par_mass = min_vals[2];
    min_par_vx = min_vals[3]; min_par_vy = min_vals[4]; min_par_vz = min_vals[5];
    min_term_v = min_vals[6];

    max_par_radius = max_vals[0]; max_multiplic = max_vals[1]; max_par_mass = max_vals[2];
    max_par_vx = max_vals[3]; max_par_vy = max_vals[4]; max_par_vz = max_vals[5];
    max_term_v = max_vals[6];

    avg_par_radius = sum_vals[0]; avg_multiplic = sum_vals[1]; avg_par_mass = sum_vals[2];
    avg_par_vx = sum_vals[3]; avg_par_vy = sum_vals[4]; avg_par_vz = sum_vals[5];
    avg_term_v = sum_vals[6];

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    min_cond_t = min_vals[7];
    max_cond_t = max_vals[7];
    avg_cond_t = sum_vals[7];
#endif

    if (num_total_particles > 0) {
        avg_par_mass /= num_total_particles;
        avg_par_radius /= num_total_particles;
        avg_multiplic /= num_superdroplets;
        avg_par_vx /= num_total_particles;
        avg_par_vy /= num_total_particles;
        avg_par_vz /= num_total_particles;
        avg_term_v /= num_total_particles;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        avg_cond_t /= num_total_particles;
#endif
    } else {
        min_par_mass = 0;
        min_par_radius = 0;
        min_multiplic = 0;
        min_par_vx = 0;
        min_par_vy = 0;
        min_par_vz = 0;
        min_term_v = 0;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        min_cond_t = 0;
#endif

        max_par_mass = 0;
        max_par_radius = 0;
        max_multiplic = 0;
        max_par_vx = 0;
        max_par_vy = 0;
        max_par_vz = 0;
        max_term_v = 0;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        max_cond_t = 0;
#endif

        avg_par_mass = 0;
        avg_par_radius = 0;
        avg_multiplic = 0;
        avg_par_vx = 0;
        avg_par_vy = 0;
        avg_par_vz = 0;
        avg_term_v = 0;
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        avg_cond_t = 0;
#endif
    }

    std::vector<Real> min_mass_species(m_num_species);
    std::vector<Real> max_mass_species(m_num_species);
    std::vector<Real> avg_mass_species(m_num_species);

    auto na = m_num_aerosols;
    auto ns = m_num_species;

    for (int is = 0; is < m_num_species; is++) {
        min_mass_species[is] = ReduceMin( *this,
                                          [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                          { return ptd.m_runtime_rdata[ridx_s(is,na,ns)][i]; } );
        max_mass_species[is] = ReduceMax( *this,
                                          [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                          { return ptd.m_runtime_rdata[ridx_s(is,na,ns)][i]; } );
        avg_mass_species[is] = ReduceSum( *this,
                                          [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                          {
                                              auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                              auto m = ptd.m_runtime_rdata[ridx_s(is,na,ns)][i];
                                              return n*m;
                                          } );
    }
    ParallelDescriptor::ReduceRealMin(min_mass_species.data(),m_num_species);
    ParallelDescriptor::ReduceRealMax(max_mass_species.data(),m_num_species);
    ParallelDescriptor::ReduceRealSum(avg_mass_species.data(),m_num_species);
    for (int is = 0; is < m_num_species; is++) {
        if (num_total_particles > 0) {
            avg_mass_species[is] /= num_total_particles;
        } else {
            min_mass_species[is] = zero;
            max_mass_species[is] = zero;
            avg_mass_species[is] = zero;
        }
    }

    std::vector<Real> min_mass_aerosols(m_num_aerosols);
    std::vector<Real> max_mass_aerosols(m_num_aerosols);
    std::vector<Real> avg_mass_aerosols(m_num_aerosols);

    for (int ia = 0; ia < m_num_aerosols; ia++) {
        min_mass_aerosols[ia] = ReduceMin( *this,
                                           [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                           { return ptd.m_runtime_rdata[ridx_a(ia,na,ns)][i]; } );
        max_mass_aerosols[ia] = ReduceMax( *this,
                                           [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                           { return ptd.m_runtime_rdata[ridx_a(ia,na,ns)][i]; } );
        avg_mass_aerosols[ia] = ReduceSum( *this,
                                           [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                           {
                                               auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                               auto m = ptd.m_runtime_rdata[ridx_a(ia,na,ns)][i];
                                               return n*m;
                                           } );
    }
    ParallelDescriptor::ReduceRealMin(min_mass_aerosols.data(),m_num_aerosols);
    ParallelDescriptor::ReduceRealMax(max_mass_aerosols.data(),m_num_aerosols);
    ParallelDescriptor::ReduceRealSum(avg_mass_aerosols.data(),m_num_aerosols);
    for (int ia = 0; ia < m_num_aerosols; ia++) {
        if (num_total_particles > 0) {
            avg_mass_aerosols[ia] /= num_total_particles;
        } else {
            min_mass_aerosols[ia] = zero;
            max_mass_aerosols[ia] = zero;
            avg_mass_aerosols[ia] = zero;
        }
    }

    if (num_total_particles > 0) {
        Print() << "SuperDropletPC(" << m_name << ") attributes (min, max, avg):\n"
                << "    mass [kg]: "
                << min_par_mass << ", " << max_par_mass << ", " << avg_par_mass << "\n"
                << "    radius [m]: "
                << min_par_radius << ", " << max_par_radius << ", " << avg_par_radius << "\n"
                << "    multiplicity: "
                << min_multiplic << ", " << max_multiplic << ", " << avg_multiplic << "\n"
                << "    velocity components [m/s]:\n"
                << "        x: "
                << min_par_vx << ", " << max_par_vx << ", " << avg_par_vx << "\n"
                << "        y: "
                << min_par_vy << ", " << max_par_vy << ", " << avg_par_vy << "\n"
                << "        z: "
                << min_par_vz << ", " << max_par_vz << ", " << avg_par_vz << "\n"
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
                << "    condensation/evaporation tendency [m/s]: "
                << min_cond_t << ", " << max_cond_t << ", " << avg_cond_t << "\n"
#endif
                << "    terminal velocity [m/s]: "
                << min_term_v << ", " << max_term_v << ", " << avg_term_v << "\n";

        Print() << "    species masses [kg]:\n";
        for (int is = 0; is < m_num_species; is++) {
            Print() << "        " << getEnumNameString(m_species_mat[is]->m_name)
                    << ": "
                    << min_mass_species[is] << ", " << max_mass_species[is] << ", " << avg_mass_species[is] << "\n";
        }
        Print() << "    aerosol masses [kg]:\n";
        for (int ia = 0; ia < m_num_aerosols; ia++) {
            Print() << "        " << getEnumNameString(m_aerosol_mat[ia]->m_name)
                    << ": "
                    << min_mass_aerosols[ia] << ", " << max_mass_aerosols[ia] << ", " << avg_mass_aerosols[ia] << "\n";
        }
    }


    Long num_unconverged_particles = m_num_unconverged_particles;
    m_num_unconverged_particles = 0;
    ParallelDescriptor::ReduceLongSum( &num_unconverged_particles, 1 );

    if (num_unconverged_particles > 0) {
        Print() << "SuperDropletPC::MassChange(): Warning - "
                << num_unconverged_particles
                << " particles did not converge.\n";
    }

    if (a_flag && (num_total_particles > 0)) {
        auto r_eff_min = min_par_radius;
        auto r_eff_max = max_par_radius;
        for (int ia = 0; ia < m_num_species; ia++) {
            const auto rho = m_species_mat[ia]->m_density;
            auto r_eff_min_species = std::cbrt( min_mass_species[ia] / (four_thirds_pi*rho) );
            if ((r_eff_min_species < r_eff_min) && (r_eff_min_species > Real(1.0e-10))) { r_eff_min = r_eff_min_species; }
        }
        for (int ia = 0; ia < m_num_aerosols; ia++) {
            const auto rho = m_aerosol_mat[ia]->m_density;
            auto r_eff_min_aero = std::cbrt( min_mass_aerosols[ia] / (four_thirds_pi*rho) );
            if ((r_eff_min_aero < r_eff_min) && (r_eff_min_aero > Real(1.0e-10))) { r_eff_min = r_eff_min_aero; }
        }
        ComputeDistributions( a_iter, a_lev, r_eff_min, r_eff_max );
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        ComputeBinnedDistributions( a_iter, a_lev);
        ComputeBinnedDistributionsCell( a_iter, a_lev, a_time);
#else
        amrex::ignore_unused(a_time);
#endif
    }
}

/*! Compute and write the distributions (as a function of the log of
    the droplet radius. The file written is a text file with multiple columns:
    R, g_mass(ln R), g_n(ln R) */
void SuperDropletPC::ComputeDistributions( const int a_iter,
                                           const int a_lev,
                                           const ParticleReal a_r_min,
                                           const ParticleReal a_r_max )
{
    BL_PROFILE("SuperDropletPC::ComputeDistributions()");
    int Nr = m_distribution_grid_size;

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto dxi = geom.InvCellSizeArray();

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = one / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[1])
                 * static_cast<ParticleReal>(m_coalescence_bin_size[2]) );
    const ParticleReal inv_bin_volume = inv_cell_volume*inv_bin_size;

    Vector<Real> ln_R, g_mass_ln_R;
    ln_R.resize(Nr);
    g_mass_ln_R.resize(Nr);

    Vector<Vector<Real>> g_amass_ln_R;
    g_amass_ln_R.resize(m_num_aerosols+m_num_species);
    for (int ia = 0; ia < g_amass_ln_R.size(); ia++) {
        g_amass_ln_R[ia].resize(Nr);
    }

    // Set ln R grid
    for (int n = 0; n < Nr; n++) {
        ln_R[n] = std::log(a_r_min) + n*(std::log(a_r_max)-std::log(a_r_min))/(Nr-1);
    }

    const auto np = NumSuperDroplets();
    const ParticleReal sigma = m_sigma0 * std::pow(static_cast<ParticleReal>(np), Real(-0.2));
    const ParticleReal lambda = one / (two*sigma*sigma);
    const ParticleReal gamma = one/(std::sqrt(two*PI)*sigma) * inv_bin_volume;

    const auto rho_w = m_species_mat[m_idx_w]->m_density;
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx_w = m_idx_w;

    // compute g(ln R)
    for (int n = 0; n < Nr; n++) {
        const auto lnR = ln_R[n];
        g_mass_ln_R[n] = ReduceSum(  *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     {
                                         auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                         auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                                         auto ri = std::cbrt( mi / (four_thirds_pi*rho_w) );
                                         auto lnRi = std::log(ri);
                                         return gamma*ai*ni*mi*std::exp(-lambda*(lnR-lnRi)*(lnR-lnRi));
                                     } );
        for (int ia = 0; ia < m_num_aerosols; ia++) {
            const auto rho = m_aerosol_mat[ia]->m_density;
            g_amass_ln_R[ia][n] = ReduceSum(*this,
                                            [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                            {
                                                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                                auto mi = ptd.m_runtime_rdata[ridx_a(ia,na,ns)][i];
                                                auto ri = std::cbrt( mi / (four_thirds_pi*rho) );
                                                auto lnRi = std::log(ri);
                                                return gamma*ai*ni*mi*std::exp(-lambda*(lnR-lnRi)*(lnR-lnRi));
                                            } );
        }
        for (int ia = m_num_aerosols; ia < m_num_aerosols+m_num_species; ia++) {
            const int is = ia - m_num_aerosols;
            const auto rho = m_species_mat[is]->m_density;
            g_amass_ln_R[ia][n] = ReduceSum(*this,
                                            [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                            {
                                                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                                auto mi = ptd.m_runtime_rdata[ridx_s(is,na,ns)][i];
                                                auto ri = std::cbrt( mi / (four_thirds_pi*rho) );
                                                auto lnRi = std::log(ri);
                                                return gamma*ai*ni*mi*std::exp(-lambda*(lnR-lnRi)*(lnR-lnRi));
                                            } );
        }
    }

    // Sum g(ln R) over MPI subdomains
    ParallelDescriptor::ReduceRealSum(g_mass_ln_R.dataPtr(),Nr);
    for (int ia = 0; ia < g_amass_ln_R.size(); ia++) {
        ParallelDescriptor::ReduceRealSum(g_amass_ln_R[ia].dataPtr(),Nr);
    }

    // Write to file
    char iter_str[12]; snprintf(iter_str, sizeof(iter_str), "%05d", a_iter+1);
    std::string output_filename =   m_name
                                    + "_g_lnR_"
                                    + std::string(iter_str) + ".txt";
    Print() << "Writing " << output_filename << "\n";
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream outfile;
        outfile.open(output_filename.c_str(), std::ios::out|std::ios::trunc);
        if (!outfile.good()) { amrex::FileOpenFailed(output_filename); }

        for (int n = 0; n < Nr; n++) {
            outfile << std::exp(ln_R[n])
                    << " " << g_mass_ln_R[n];
            for (int ia = 0; ia < g_amass_ln_R.size(); ia++) {
                outfile << " " << g_amass_ln_R[ia][n];
            }
            outfile << "\n";
        }

        outfile.flush();
        outfile.close();
        if (!outfile.good()) { amrex::Abort("problem writing output file"); }
    }
}

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
/*! Compute and write the distributions (as a function of the log of
    the droplet radius. The file written is a text file with multiple columns:
    R, g_mass(ln R), g_n(ln R) */
void SuperDropletPC::ComputeBinnedDistributions( const int a_iter, const int a_lev)
{
    BL_PROFILE("SuperDropletPC::ComputeBinnedDistributions()");
    int Nbin = m_distribution_grid_size;
    auto r_min = m_bindist_rmin;
    auto r_max = m_bindist_rmax;

    const Geometry& geom = m_gdb->Geom(a_lev);
    const auto dxi = geom.InvCellSizeArray();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2]; // divide by cell volume

    Vector<Real> ln_R, g_num_ln_R, g_mass_ln_R;
    ln_R.resize(Nbin+1); // one extra: denotes left and right bounds of each DSD bin
    g_num_ln_R.resize(Nbin);
    g_mass_ln_R.resize(Nbin);

    // Set ln R grid
    auto dln_R = (std::log(r_max)-std::log(r_min))/(Nbin);
    for (int n = 0; n < Nbin+1; n++) {
        ln_R[n] = std::log(r_min) + n*dln_R;
    }

    const auto density = m_species_mat[m_idx_w]->m_density;
    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx_w = m_idx_w;

    // compute g(ln R)
    for (int n = 0; n < Nbin; n++) {
        auto r_l = std::exp(ln_R[n]);
        auto r_r = std::exp(ln_R[n+1]);

        g_mass_ln_R[n] = ReduceSum(  *this,
                                    [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                    {
                                        auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                        auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                        auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                                        auto ri = std::cbrt( mi / (four_thirds_pi*density) );
                                        auto inbin = (r_l <= ri && ri < r_r) ? one : zero;
                                        return ai*ni*mi*inbin * inv_cell_volume / dln_R;
                                    } );
        g_num_ln_R[n] = ReduceSum(  *this,
                                    [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                    {
                                        auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                        auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                        auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                                        auto ri = std::cbrt( mi / (four_thirds_pi*density) );
                                        auto inbin = (r_l <= ri && ri < r_r) ? one : zero;
                                        return ai*ni*inbin * inv_cell_volume / dln_R;
                                    } );
    }

    // Sum g(ln R) over MPI subdomains
    ParallelDescriptor::ReduceRealSum(g_mass_ln_R.dataPtr(), Nbin);
    ParallelDescriptor::ReduceRealSum(g_num_ln_R.dataPtr(), Nbin);

    // Write to file
    char iter_str[12]; snprintf(iter_str, sizeof(iter_str), "%05d", a_iter+1);
    std::string output_filename =   m_name
                                    + "_binned_dsd_"
                                    + std::string(iter_str) + ".txt";
                                    Print() << "Writing " << output_filename << "\n";
    if (ParallelDescriptor::IOProcessor()) {
    std::ofstream outfile;
    outfile.open(output_filename.c_str(), std::ios::out|std::ios::trunc);
    if (!outfile.good()) { amrex::FileOpenFailed(output_filename); }

    for (int n = 0; n < Nbin; n++) {
        outfile << std::exp(ln_R[n])
                << " " << g_mass_ln_R[n]
                << " " << g_num_ln_R[n];
        outfile << "\n";
    }
    outfile << std::exp(ln_R[Nbin]) << "\n";

    outfile.flush();
    outfile.close();
    if (!outfile.good()) { amrex::Abort("problem writing output file"); }
    }
}

/*! Compute and write the cell-wise distributions (as a function of the log of
    the droplet radius. The file written is a text file with multiple columns:
    R, g_mass(ln R), g_n(ln R) */
void SuperDropletPC::ComputeBinnedDistributionsCell( const int a_iter,
                                                     const int a_lev,
                                                     const Real a_time )
{
    BL_PROFILE("SuperDropletPC::ComputeBinnedDistributionsCell()");
    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    int Nbin = m_distribution_grid_size;
    auto r_min = m_bindist_rmin;
    auto r_max = m_bindist_rmax;

    const auto& geom = Geom(a_lev);
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    const auto domain = geom.Domain();
    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2]; // divide by cell volume

    // Set ln R grid
    Vector<Real> ln_R;
    ln_R.resize(Nbin+1); // one extra: denotes left and right bounds of each DSD bin
    auto dln_R = (std::log(r_max)-std::log(r_min))/(Nbin);
    for (int n = 0; n < Nbin+1; n++) {
        ln_R[n] = std::log(r_min) + n*dln_R;
    }

    auto density = m_species_mat[m_idx_w]->m_density;
    Vector<std::string> varnames(Nbin,"");
    m_mass_ln_R_mf.setVal(0.0);
    m_num_ln_R_mf.setVal(0.0);

    const auto na = m_num_aerosols;
    const auto ns = m_num_species;
    const auto idx_w = m_idx_w;

    for (int n = 0; n < Nbin; n++) {
        auto r_l = std::exp(ln_R[n]);
        auto r_r = std::exp(ln_R[n+1]);

        char r_str[12]; snprintf(r_str, sizeof(r_str), "%1.4e", r_l);
        varnames[n] = std::string(r_str);

        ParticleToMesh( *this, m_mass_ln_R_mf, a_lev,
            [=] AMREX_GPU_DEVICE (  const SDTDType& ptd, int i, Array4<Real> const& mf_arr)
            {
                auto p = ptd.m_aos[i];
                auto iv = getParticleCell(p, plo, dxi, domain);

                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                auto ri = std::cbrt( mi / (four_thirds_pi*density) );
                auto inbin = (r_l <= ri && ri < r_r) ? one : zero;

                Gpu::Atomic::AddNoRet(&mf_arr(iv, n), (ai*ni*mi*inbin * inv_cell_volume / dln_R));
            }, false);

        ParticleToMesh( *this, m_num_ln_R_mf, a_lev,
            [=] AMREX_GPU_DEVICE (  const SDTDType& ptd, int i, Array4<Real> const& mf_arr)
            {
                auto p = ptd.m_aos[i];
                auto iv = getParticleCell(p, plo, dxi, domain);

                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                auto ri = std::cbrt( mi / (four_thirds_pi*density) );
                auto inbin = (r_l <= ri && ri < r_r) ? one : zero;

                Gpu::Atomic::AddNoRet(&mf_arr(iv, n), (ai*ni*inbin * inv_cell_volume / dln_R) );
            }, false);

    }


    // Write to file
    {
        char iter_str[12]; snprintf(iter_str, sizeof(iter_str), "%05d", a_iter+1);
        std::string op_fname =   m_name
                               + "_binned_dsd_mass_"
                               + std::string(iter_str);
        Print() << "Writing " << op_fname << "\n";
        WriteSingleLevelPlotfile ( op_fname,
                                   m_mass_ln_R_mf,
                                   varnames,
                                   geom,
                                   a_time,
                                   a_iter );
    }
    {
        char iter_str[12]; snprintf(iter_str, sizeof(iter_str), "%05d", a_iter+1);
        std::string op_fname =   m_name
                               + "_binned_dsd_number_"
                               + std::string(iter_str);
        Print() << "Writing " << op_fname << "\n";
        WriteSingleLevelPlotfile ( op_fname,
                                   m_num_ln_R_mf,
                                   varnames,
                                   geom,
                                   a_time,
                                   a_iter );
    }
}
#endif

#endif
