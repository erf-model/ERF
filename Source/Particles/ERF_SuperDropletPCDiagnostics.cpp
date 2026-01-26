#include <stdio.h>
#include <limits>
#include <fstream>
#include <AMReX_PlotFileUtil.H>
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
                                  const Real a_time,
                                  const bool a_flag )
{
    BL_PROFILE("SuperDropletPC::Diagnostics()");

    auto num_total_particles = TotalNumberOfParticles();
    auto num_superdroplets = static_cast<Real>(NumSuperDroplets());

    auto min_par_radius = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i]; } );

    auto max_par_radius = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i]; } );

    auto avg_par_radius = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     {
                                         auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto r = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::radius][i];
                                         return n*r;
                                     } );

    auto min_multiplic  = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; } );

    auto max_multiplic  = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; } );

    auto avg_multiplic  = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i]; } );

    auto min_par_mass   = ReduceMin( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i]; } );

    auto max_par_mass   = ReduceMax( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     { return ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i]; } );

    auto avg_par_mass   = ReduceSum( *this,
                                     [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                     {
                                         auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                         auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::mass][i];
                                         return n*m;
                                     } );

    auto min_par_vx = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vx][i]; } );

    auto max_par_vx = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vx][i]; } );

    auto avg_par_vx = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::vx][i];
                                     return n*m;
                                 } );

    auto min_par_vy = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vy][i]; } );

    auto max_par_vy = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vy][i]; } );

    auto avg_par_vy = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::vy][i];
                                     return n*m;
                                 } );

    auto min_par_vz = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vz][i]; } );

    auto max_par_vz = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_rdata[SuperDropletsRealIdxSoA::vz][i]; } );

    auto avg_par_vz = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_rdata[SuperDropletsRealIdxSoA::vz][i];
                                     return n*m;
                                 } );

    auto min_term_v = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i]; } );

    auto max_term_v = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i]; } );

    auto avg_term_v = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::term_vel][i];
                                     return n*m;
                                 } );

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    auto min_cond_t = ReduceMin( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::cond_tendency][i]; } );

    auto max_cond_t = ReduceMax( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 { return ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::cond_tendency][i]; } );

    auto avg_cond_t = ReduceSum( *this,
                                 [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                 {
                                     auto n = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                     auto m = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::cond_tendency][i];
                                     return n*m;
                                 } );
#endif

    ParallelDescriptor::ReduceRealMin(&min_par_mass,1);
    ParallelDescriptor::ReduceRealMin(&min_par_radius,1);
    ParallelDescriptor::ReduceRealMin(&min_multiplic,1);
    ParallelDescriptor::ReduceRealMin(&min_par_vx,1);
    ParallelDescriptor::ReduceRealMin(&min_par_vy,1);
    ParallelDescriptor::ReduceRealMin(&min_par_vz,1);
    ParallelDescriptor::ReduceRealMin(&min_term_v,1);
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    ParallelDescriptor::ReduceRealMin(&min_cond_t,1);
#endif

    ParallelDescriptor::ReduceRealMax(&max_par_mass,1);
    ParallelDescriptor::ReduceRealMax(&max_par_radius,1);
    ParallelDescriptor::ReduceRealMax(&max_multiplic,1);
    ParallelDescriptor::ReduceRealMax(&max_par_vx,1);
    ParallelDescriptor::ReduceRealMax(&max_par_vy,1);
    ParallelDescriptor::ReduceRealMax(&max_par_vz,1);
    ParallelDescriptor::ReduceRealMax(&max_term_v,1);
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    ParallelDescriptor::ReduceRealMax(&max_cond_t,1);
#endif

    ParallelDescriptor::ReduceRealSum(&avg_par_mass,1);
    ParallelDescriptor::ReduceRealSum(&avg_par_radius,1);
    ParallelDescriptor::ReduceRealSum(&avg_multiplic,1);
    ParallelDescriptor::ReduceRealSum(&avg_par_vx,1);
    ParallelDescriptor::ReduceRealSum(&avg_par_vy,1);
    ParallelDescriptor::ReduceRealSum(&avg_par_vz,1);
    ParallelDescriptor::ReduceRealSum(&avg_term_v,1);
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
    ParallelDescriptor::ReduceRealSum(&avg_cond_t,1);
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
            min_mass_species[is] = 0.0;
            max_mass_species[is] = 0.0;
            avg_mass_species[is] = 0.0;
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
            min_mass_aerosols[ia] = 0.0;
            max_mass_aerosols[ia] = 0.0;
            avg_mass_aerosols[ia] = 0.0;
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
            auto r_eff_min_species = std::cbrt( min_mass_species[ia] / ((4.0/3.0)*PI*rho) );
            if ((r_eff_min_species < r_eff_min) && (r_eff_min_species > 1.0e-10)) { r_eff_min = r_eff_min_species; }
        }
        for (int ia = 0; ia < m_num_aerosols; ia++) {
            const auto rho = m_aerosol_mat[ia]->m_density;
            auto r_eff_min_aero = std::cbrt( min_mass_aerosols[ia] / ((4.0/3.0)*PI*rho) );
            if ((r_eff_min_aero < r_eff_min) && (r_eff_min_aero > 1.0e-10)) { r_eff_min = r_eff_min_aero; }
        }
        ComputeDistributions( a_iter, r_eff_min, r_eff_max );
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        ComputeBinnedDistributions( a_iter);
        ComputeBinnedDistributionsCell( a_iter, a_time);
#else
        amrex::ignore_unused(a_time);
#endif
    }
}

/*! Compute and write the distributions (as a function of the log of
    the droplet radius. The file written is a text file with multiple columns:
    R, g_mass(ln R), g_n(ln R) */
void SuperDropletPC::ComputeDistributions( const int a_iter,
                                           const ParticleReal a_r_min,
                                           const ParticleReal a_r_max )
{
    BL_PROFILE("SuperDropletPC::ComputeDistributions()");
    int Nr = m_distribution_grid_size;

    const Geometry& geom = m_gdb->Geom(m_lev);
    const auto dxi = geom.InvCellSizeArray();

    const ParticleReal inv_cell_volume = dxi[0]*dxi[1]*dxi[2];
    const ParticleReal inv_bin_size
        = 1.0 / (  static_cast<ParticleReal>(m_coalescence_bin_size[0])
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
    const ParticleReal sigma = m_sigma0 * std::exp(-0.2*std::log(static_cast<ParticleReal>(np)));
    const ParticleReal lambda = 1.0 / (2.0*sigma*sigma);
    const ParticleReal gamma = 1.0/(std::sqrt(2.0*PI)*sigma) * inv_bin_volume;

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
                                         auto ri = std::cbrt( mi / ((4.0/3.0)*PI*rho_w) );
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
                                                auto ri = std::cbrt( mi / ((4.0/3.0)*PI*rho) );
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
                                                auto ri = std::cbrt( mi / ((4.0/3.0)*PI*rho) );
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
void SuperDropletPC::ComputeBinnedDistributions( const int a_iter)
{
    BL_PROFILE("SuperDropletPC::ComputeBinnedDistributions()");
    int Nbin = m_distribution_grid_size;
    auto r_min = m_bindist_rmin;
    auto r_max = m_bindist_rmax;

    const Geometry& geom = m_gdb->Geom(m_lev);
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
                                        auto ri = std::cbrt( mi / ((4.0/3.0)*PI*density) );
                                        auto inbin = (r_l <= ri && ri < r_r) ? 1.0 : 0.0;
                                        return ai*ni*mi*inbin * inv_cell_volume / dln_R;
                                    } );
        g_num_ln_R[n] = ReduceSum(  *this,
                                    [=] AMREX_GPU_HOST_DEVICE (const SDTDType& ptd, const int i) -> Real
                                    {
                                        auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                                        auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                                        auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                                        auto ri = std::cbrt( mi / ((4.0/3.0)*PI*density) );
                                        auto inbin = (r_l <= ri && ri < r_r) ? 1.0 : 0.0;
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
                                                     const Real a_time )
{
    BL_PROFILE("SuperDropletPC::ComputeBinnedDistributionsCell()");
    AMREX_ASSERT(OK());
    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    int Nbin = m_distribution_grid_size;
    auto r_min = m_bindist_rmin;
    auto r_max = m_bindist_rmax;

    const auto& geom = Geom(m_lev);
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

        ParticleToMesh( *this, m_mass_ln_R_mf, m_lev,
            [=] AMREX_GPU_DEVICE (  const SDTDType& ptd, int i, Array4<Real> const& mf_arr)
            {
                auto p = ptd.m_aos[i];
                auto iv = getParticleCell(p, plo, dxi, domain);

                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                auto ri = std::cbrt( mi / ((4.0/3.0)*PI*density) );
                auto inbin = (r_l <= ri && ri < r_r) ? 1.0 : 0.0;

                Gpu::Atomic::AddNoRet(&mf_arr(iv, n), (ai*ni*mi*inbin * inv_cell_volume / dln_R));
            }, false);

        ParticleToMesh( *this, m_num_ln_R_mf, m_lev,
            [=] AMREX_GPU_DEVICE (  const SDTDType& ptd, int i, Array4<Real> const& mf_arr)
            {
                auto p = ptd.m_aos[i];
                auto iv = getParticleCell(p, plo, dxi, domain);

                auto ai = ptd.m_runtime_idata[SuperDropletsIntIdxSoA_RT::active][i];
                auto ni = ptd.m_runtime_rdata[SuperDropletsRealIdxSoA_RT::multiplicity][i];
                auto mi = ptd.m_runtime_rdata[ridx_s(idx_w,na,ns)][i];
                auto ri = std::cbrt( mi / ((4.0/3.0)*PI*density) );
                auto inbin = (r_l <= ri && ri < r_r) ? 1.0 : 0.0;

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
