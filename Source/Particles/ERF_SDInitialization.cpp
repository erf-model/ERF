#include <cmath>
#include "ERF_SDInitialization.H"

using MatVec = std::vector<std::unique_ptr<MaterialProperties>>;

void SDInitProperties::setDefaults ( const amrex::Geometry& a_geom,
                                     const MatVec& a_species_mat,
                                     const MatVec& a_aerosol_mat )
{
    BL_PROFILE("SDInitProperties::setDefaults");

    // Default
    m_init_particle_p1.resize(AMREX_SPACEDIM);
    m_init_particle_p2.resize(AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        m_init_particle_p1[i] = a_geom.ProbLo(i);
        m_init_particle_p2[i] = a_geom.ProbHi(i);
    }

    m_num_species = a_species_mat.size();
    m_mass_species_min.resize(m_num_species);
    m_mass_species_max.resize(m_num_species);
    m_mass_species_mean.resize(m_num_species);
    m_radius_species_min.resize(m_num_species);
    m_radius_species_max.resize(m_num_species);
    m_radius_species_mean.resize(m_num_species);
    m_radius_species_geom_std.resize(m_num_species);
    m_species_init_type.resize(m_num_species);

    for (int i = 0; i < m_num_species; i++) {
        // default values
        m_species_init_type[i] = SupDropInit::attrib_init_const;
        m_mass_species_min[i]   = 4.1887902e-42;
        m_mass_species_max[i]   = 4.1887902e-42;
        m_mass_species_mean[i]  = 4.1887902e-42;
        m_radius_species_min[i] = 1.0e-15;
        m_radius_species_max[i] = 1.0e-15;
        m_radius_species_mean[i] = 1.0e-15;
        m_radius_species_geom_std[i] = 2.0;
    }

    m_num_aerosols = a_aerosol_mat.size();
    m_mass_aerosol_min.resize(m_num_aerosols);
    m_mass_aerosol_max.resize(m_num_aerosols);
    m_mass_aerosol_mean.resize(m_num_aerosols);
    m_radius_aerosol_min.resize(m_num_aerosols);
    m_radius_aerosol_max.resize(m_num_aerosols);
    m_radius_aerosol_mean.resize(m_num_aerosols);
    m_radius_aerosol_geom_std.resize(m_num_aerosols);
    m_aerosol_init_type.resize(m_num_aerosols);

    for (int i = 0; i < m_num_aerosols; i++) {
        // default values
        m_aerosol_init_type[i] = SupDropInit::attrib_init_const;
        m_mass_aerosol_min[i]   = 0.0;
        m_mass_aerosol_max[i]   = 0.0;
        m_mass_aerosol_mean[i]  = 0.0;
        m_radius_aerosol_min[i] = 1.0e-9;
        m_radius_aerosol_max[i] = 1.0e-6;
        m_radius_aerosol_mean[i] = 1.0e-40;
        m_radius_aerosol_geom_std[i] = 2.0;
    }

    m_mult_type = SDMultiplicityType::sampled;
}

void SDInitProperties::readInputs ( const std::string& a_prefix,
                                    const std::string& a_key,
                                    const amrex::Geometry& a_geom,
                                    const MatVec& a_species_mat,
                                    const MatVec& a_aerosol_mat )
{
    BL_PROFILE("SDInitProperties::readInputs");
    amrex::ignore_unused(a_geom);
    using namespace amrex;

    amrex::ParmParse pp(a_prefix);
    pp.query(std::string(a_key+"distribution_type").c_str(), m_type);
    pp.query("maximum_multiplicity", m_max_multiplicity);
    pp.query("multiplicity_type", m_mult_type);

    pp.query(std::string(a_key+"particles_per_cell").c_str(), m_ppc);

    if (m_type == SDInitShape::uniform) {

        pp.queryAdd("particle_box_lo", m_init_particle_p1, AMREX_SPACEDIM);
        AMREX_ASSERT(m_init_particle_p1.size() == AMREX_SPACEDIM);

        pp.queryAdd("particle_box_hi", m_init_particle_p2, AMREX_SPACEDIM);
        AMREX_ASSERT(m_init_particle_p2.size() == AMREX_SPACEDIM);

        m_particle_domain.setLo(m_init_particle_p1);
        m_particle_domain.setHi(m_init_particle_p2);
    } else if (m_type == SDInitShape::bubble){

        pp.queryAdd("particle_bubble_center", m_init_particle_p1, AMREX_SPACEDIM);
        AMREX_ASSERT(m_init_particle_p1.size() == AMREX_SPACEDIM);

        pp.queryAdd("particle_bubble_radius", m_init_particle_p2, AMREX_SPACEDIM);
        AMREX_ASSERT(m_init_particle_p2.size() == AMREX_SPACEDIM);

        m_particle_domain.setLo(m_init_particle_p1);
        m_particle_domain.setHi(m_init_particle_p2);
    }

    // Backward compatibility
    for (int i = 0; i < m_num_species; i++) {
        pp.query(std::string(a_key+"condensate_distribution_type").c_str(), m_species_init_type[i]);
        pp.query(std::string(a_key+"condensate_mass_min").c_str(), m_mass_species_min[i]);
        pp.query(std::string(a_key+"condensate_mass_mean").c_str(), m_mass_species_mean[i]);
        pp.query(std::string(a_key+"condensate_min_radius").c_str(), m_radius_species_min[i]);
        pp.query(std::string(a_key+"condensate_max_radius").c_str(), m_radius_species_max[i]);
    }
    for (int i = 0; i < m_num_species; i++) {
        {
            std::string key = a_key+"species_distribution_type_"+getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_species_init_type[i]);
        }
        {
            std::string key = a_key+"species_min_mass_" + getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_mass_species_min[i]);
        }
        {
            std::string key = a_key+"species_mean_mass_" + getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_mass_species_mean[i]);
        }
        {
            m_mass_species_max[i] = 5 * m_mass_species_mean[i]; // default
            std::string key = a_key+"species_max_mass_" + getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_mass_species_max[i]);
        }
        {
            std::string key = a_key+"species_min_radius_" + getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_radius_species_min[i]);
        }
        {
            std::string key = a_key+"species_max_radius_" + getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_radius_species_max[i]);
        }
        {
            std::string key = a_key+"species_mean_radius_" + getEnumNameString(a_species_mat[i]->m_name);
            pp.query(key.c_str(), m_radius_species_mean[i]);
        }
        {
            std::string key_std = a_key+"species_std_radius_" + getEnumNameString(a_species_mat[i]->m_name);
            std::string key_gstd = a_key+"species_geomstd_radius_" + getEnumNameString(a_species_mat[i]->m_name);
            if (pp.contains(key_std.c_str()) && pp.contains(key_gstd.c_str())) {
                amrex::Abort("Cannot specify BOTH initial_species_std_radius and initial_species_geomstd_radius");
            }
            if (pp.contains(key_std.c_str())) {
                pp.get(key_std.c_str(), m_radius_species_geom_std[i]);
                m_radius_species_geom_std[i] = std::exp(m_radius_species_geom_std[i]);
            } else {
                pp.query(key_gstd.c_str(), m_radius_species_geom_std[i]);
            }
        }
    }

    for (int i = 0; i < m_num_aerosols; i++) {
        {
            std::string key = a_key+"aerosol_distribution_type_"+getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_aerosol_init_type[i]);
        }
        {
            std::string key = a_key+"aerosol_min_mass_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_mass_aerosol_min[i]);
        }
        {
            std::string key = a_key+"aerosol_mean_mass_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_mass_aerosol_mean[i]);
        }
        {
            m_mass_aerosol_max[i] = 5 * m_mass_aerosol_mean[i]; // default
            std::string key = a_key+"aerosol_max_mass_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_mass_aerosol_max[i]);
        }
        {
            std::string key = a_key+"aerosol_min_radius_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_radius_aerosol_min[i]);
        }
        {
            std::string key = a_key+"aerosol_max_radius_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_radius_aerosol_max[i]);
        }
        {
            std::string key = a_key+"aerosol_mean_radius_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            pp.query(key.c_str(), m_radius_aerosol_mean[i]);
        }
        {
            std::string key_std = a_key+"aerosol_std_radius_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            std::string key_gstd = a_key+"aerosol_geomstd_radius_" + getEnumNameString(a_aerosol_mat[i]->m_name);
            if (pp.contains(key_std.c_str()) && pp.contains(key_gstd.c_str())) {
                amrex::Abort("Cannot specify BOTH initial_species_std_radius and initial_species_geomstd_radius");
            }
            if (pp.contains(key_std.c_str())) {
                pp.get(key_std.c_str(), m_radius_aerosol_geom_std[i]);
                m_radius_aerosol_geom_std[i] = std::exp(m_radius_aerosol_geom_std[i]);
            } else {
                pp.query(key_gstd.c_str(), m_radius_aerosol_geom_std[i]);
            }
        }
    }

}

void SDInitialization::readInputs ( const std::string& a_prefix,
                                    const amrex::Geometry& a_geom,
                                    const MatVec& a_species_mat,
                                    const MatVec& a_aerosol_mat )
{
    BL_PROFILE("SDInitialization::readInputs");

    SDInitProperties::readInputs( a_prefix, "initial_", a_geom, a_species_mat, a_aerosol_mat);

    amrex::ignore_unused(a_geom);
    using namespace amrex;

    amrex::ParmParse pp(a_prefix);
    pp.query("initial_number_density", this->m_numdens);
    pp.query("initial_super_droplet_density", m_numdens_sd_init);
}

void SDInjection::readInputs ( const std::string& a_prefix,
                               const amrex::Geometry& a_geom,
                               const MatVec& a_species_mat,
                               const MatVec& a_aerosol_mat,
                               const amrex::Real a_dt )
{
    BL_PROFILE("SDInjection::readInputs");

    SDInitProperties::readInputs( a_prefix, "", a_geom, a_species_mat, a_aerosol_mat);

    amrex::ignore_unused(a_geom);
    using namespace amrex;

    amrex::ParmParse pp(a_prefix);
    pp.query("rate", m_inj_rate);
    pp.query("sd_rate", m_sd_inj_rate);
    pp.query("t_start", m_tstart);
    pp.query("t_stop", m_tstop);
    pp.queryarr("domain_velocity", m_domain_vel);

    this->m_numdens = m_inj_rate * a_dt;
    m_numdens_sd = (m_sd_inj_rate > 0 ? std::max(m_sd_inj_rate*a_dt, 1.0) : -1);
}

void SDInitProperties::printParameters ( const MatVec& a_species_mat,
                                         const MatVec& a_aerosol_mat ) const
{
    using namespace amrex;
    if (m_type == SDInitShape::uniform) {
        Print() << "    Particle box: " << m_particle_domain << "\n";
    } else if (m_type == SDInitShape::bubble) {
        Print() << "    Particle bubble (radius, center): " << m_particle_domain << "\n";
    }
    Print() << "    Multiplicity type: " << amrex::getEnumNameString(m_mult_type) << "\n";
    Print() << "    Particles per cell: " << m_ppc << "\n";

    Print() << "    Vapour/Condensate Species material:\n";
    for (unsigned long i=0; i < a_species_mat.size(); i++) {
        Print() << "        "
                << getEnumNameString(a_species_mat[i]->m_name)
                << " (distribution: " << m_species_init_type[i];
        if (m_species_init_type[i] == SupDropInit::attrib_init_const) {
            Print() << ", value=" << m_mass_species_mean[i];
        } else if (m_species_init_type[i] == SupDropInit::attrib_init_exp) {
            Print() << ", min=" << m_mass_species_min[i]
                    << ", mean=" << m_mass_species_mean[i]
                    << ", max=" << m_mass_species_max[i];
            AMREX_ALWAYS_ASSERT(m_mass_species_min[i] > 0.0);
            AMREX_ALWAYS_ASSERT(m_mass_species_max[i] >= m_mass_species_min[i]);
            AMREX_ALWAYS_ASSERT(    (m_mass_species_mean[i] >= m_mass_species_min[i])
                                 && (m_mass_species_mean[i] <= m_mass_species_max[i]) );
        } else if (m_species_init_type[i] == SupDropInit::attrib_init_lnr) {
            Print() << ", min=" << m_radius_species_min[i]
                    << ", max=" << m_radius_species_max[i]
                    << ", mean=" << m_radius_species_mean[i]
                    << ", std=" << m_radius_species_geom_std[i];
            AMREX_ALWAYS_ASSERT(m_radius_species_min[i] > 0.0);
            AMREX_ALWAYS_ASSERT(m_radius_species_max[i] >= m_radius_species_min[i]);
            AMREX_ALWAYS_ASSERT(    (m_radius_species_mean[i] >= m_radius_species_min[i])
                                 && (m_radius_species_mean[i] <= m_radius_species_max[i]) );
        } else if (m_species_init_type[i] == SupDropInit::attrib_init_lnr_auto) {
            Print() << ", mean=" << m_radius_species_mean[i]
                    << ", std=" << m_radius_species_geom_std[i];
            AMREX_ALWAYS_ASSERT(m_radius_species_mean[i] > 0.0);
        }
        Print() << ")" << "\n";
    }

    if (a_aerosol_mat.size() > 0) {
        Print() << "    Aerosols material:\n";
        for (unsigned long i=0; i < a_aerosol_mat.size(); i++) {
            Print() << "        "
                    << getEnumNameString(a_aerosol_mat[i]->m_name)
                    << " (distribution: " << m_aerosol_init_type[i];
            if (m_aerosol_init_type[i] == SupDropInit::attrib_init_const) {
                Print() << ", value=" << m_mass_aerosol_mean[i];
                AMREX_ALWAYS_ASSERT(m_mass_aerosol_mean[i] > 0.0);
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_exp) {
                Print() << ", min=" << m_mass_aerosol_min[i]
                        << ", mean=" << m_mass_aerosol_mean[i]
                        << ", max=" << m_mass_aerosol_max[i];
                AMREX_ALWAYS_ASSERT(m_mass_aerosol_min[i] > 0.0);
                AMREX_ALWAYS_ASSERT(m_mass_aerosol_max[i] >= m_mass_aerosol_min[i]);
                AMREX_ALWAYS_ASSERT(    (m_mass_aerosol_mean[i] >= m_mass_aerosol_min[i])
                                     && (m_mass_aerosol_mean[i] <= m_mass_aerosol_max[i]) );
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr) {
                Print() << ", min=" << m_radius_aerosol_min[i]
                        << ", max=" << m_radius_aerosol_max[i]
                        << ", mean=" << m_radius_aerosol_mean[i]
                        << ", std=" << m_radius_aerosol_geom_std[i];
                AMREX_ALWAYS_ASSERT(m_radius_aerosol_min[i] > 0.0);
                AMREX_ALWAYS_ASSERT(m_radius_aerosol_max[i] >= m_radius_aerosol_min[i]);
                AMREX_ALWAYS_ASSERT(    (m_radius_aerosol_mean[i] >= m_radius_aerosol_min[i])
                                     && (m_radius_aerosol_mean[i] <= m_radius_aerosol_max[i]) );
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr_auto) {
                Print() << ", mean=" << m_radius_aerosol_mean[i]
                        << ", std=" << m_radius_aerosol_geom_std[i];
                AMREX_ALWAYS_ASSERT(m_radius_aerosol_mean[i] > 0.0);
            }
            Print() << ")" << "\n";
        }
    }
}

void SDInitialization::printParameters ( const MatVec& a_species_mat,
                                         const MatVec& a_aerosol_mat ) const
{
    using namespace amrex;
    Print() << "    Initial number density: " << this->m_numdens << "\n"
            << "    Initial super-droplets number density: " << m_numdens_sd_init << "\n";
    SDInitProperties::printParameters(a_species_mat, a_aerosol_mat);
}

void SDInjection::printParameters ( const MatVec& a_species_mat,
                                    const MatVec& a_aerosol_mat ) const
{
    using namespace amrex;
    Print() << "    Injection rate (# m^{-3} s^{-1}): " << m_inj_rate << "\n"
            << "    Injection domain velocity [m/s]: "
            << m_domain_vel[0] << ","
            << m_domain_vel[1] << ","
            << m_domain_vel[2] << "\n"
            << "    Time (start, stop) [s]: " << m_tstart << ", " << m_tstop << "\n"
            << "    SD injection rate: " << m_sd_inj_rate << "\n";
    SDInitProperties::printParameters(a_species_mat, a_aerosol_mat);
}

void SDInitProperties::getDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                         const int a_np,
                                         const amrex::Real a_density,
                                         const std::string& a_init_type,
                                         const amrex::Real a_mass_min,
                                         const amrex::Real /*a_mass_max*/,
                                         const amrex::Real a_mass_mean,
                                         const amrex::Real a_radius_min,
                                         const amrex::Real a_radius_max,
                                         const amrex::Real a_radius_mean,
                                         const amrex::Real a_radius_gstd,
                                         std::mt19937& a_rng ) const
{
    a_mass.resize(a_np);
    if (a_init_type == SupDropInit::attrib_init_const) {
        for (int n = 0; n < a_np; n++) {
            a_mass[n] = a_mass_mean;
        }
    } else if (a_init_type == SupDropInit::attrib_init_exp) {
        auto delta = a_mass_mean - a_mass_min;
        std::exponential_distribution<amrex::Real> ed(1.0/delta);
        for (int n = 0; n < a_np; n++) {
            a_mass[n] = ed(a_rng) + a_mass_min;
        }
    } else if (a_init_type == SupDropInit::attrib_init_lnr) {
        std::normal_distribution<amrex::Real> nrd(std::log(a_radius_mean),
                                                  std::log(a_radius_gstd));
        for (int n = 0; n < a_np; n++) {
            auto dry_r = std::exp(nrd(a_rng));
            int count = 0;
            while ((dry_r < a_radius_min) || (dry_r > a_radius_max)) {
                dry_r = std::exp(nrd(a_rng));
                count++;
                if (count > 100) { break; }
            }
            a_mass[n] = (4.0/3.0) * PI
                                * dry_r * dry_r * dry_r
                                * a_density;
        }
    } else {
        amrex::Abort("Unknown a_init_type!");
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static amrex::Real SD_erfinv(const amrex::Real x) {
    amrex::Real a = 0.147;
    amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real term = std::log(1 - x * x + eps);
    amrex::Real p1 = 2 / (PI * a) + term / 2.0;
    amrex::Real p2 = term / a;
    return std::sqrt(std::sqrt(p1 * p1 - p2) - p1);
}

void SDInitProperties::getDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                         amrex::Vector<amrex::Real>& a_mult,
                                         const amrex::Real a_dV,
                                         const int a_np,
                                         const amrex::Real a_density,
                                         const std::string& a_init_type,
                                         const amrex::Real a_mass_min,
                                         const amrex::Real a_mass_max,
                                         const amrex::Real a_mass_mean,
                                         const amrex::Real a_radius_min,
                                         const amrex::Real a_radius_max,
                                         const amrex::Real a_radius_mean,
                                         const amrex::Real a_radius_gstd,
                                         std::mt19937& a_rng ) const
{
    a_mass.resize(a_np);
    AMREX_ALWAYS_ASSERT(a_mult.size() == a_np);
    if (a_init_type == SupDropInit::attrib_init_const) {
        std::uniform_real_distribution<> urd(0.0, 1.0);
        for (int n = 0; n < a_np; n++) {
            a_mass[n] = a_mass_mean;
            a_mult[n] += urd(a_rng); // initially this will be a non-integer; later we will rescale to an integer.
        }
    } else if (a_init_type == SupDropInit::attrib_init_exp) {
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto delta = a_mass_mean - a_mass_min;
        auto lnrng = std::log(a_mass_max) - std::log(a_mass_min);
        auto lnmin = std::log(a_mass_min);
        for (int n = 0; n < a_np; n++) {
            auto tmp = lnmin + urd(a_rng) * lnrng;
            a_mass[n] = std::exp(tmp);
            a_mult[n] += (m_numdens * a_dV) * std::exp(-a_mass[n] / delta);
        }
    } else if (a_init_type == SupDropInit::attrib_init_lnr) {
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto sigma = std::log(a_radius_gstd);
        auto mu = a_radius_mean;
        auto lnrng = std::log(a_radius_max) - std::log(a_radius_min);
        auto lnmin = std::log(a_radius_min);
        for (int n = 0; n < a_np; n++) {
            auto tmp = lnmin + urd(a_rng) * lnrng;
            auto dry_r = std::exp(tmp);
            a_mass[n] = (4.0/3.0) * PI * dry_r * dry_r * dry_r * a_density;
            auto term = std::exp(-std::log(dry_r/mu)*std::log(dry_r/mu)/(2.0*sigma*sigma));
            a_mult[n] += ( m_numdens * a_dV ) / (sigma*std::sqrt(2*PI)) * term;
        }
    } else if (a_init_type == SupDropInit::attrib_init_lnr_auto) {
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto sigma = std::log(a_radius_gstd);
        auto mu = a_radius_mean;
        // automatically find the min and max radius of superdroplets, using Dziekan & Pawlowska 2017
        auto rmin = 1e-9;
        auto rmax = 1.0;
        auto dlnr = (std::log(rmax) - std::log(rmin)) / a_np;
        auto P_min = 0.0;
        auto P_max = 1.0;
        auto tol = 1.0 / (m_numdens * a_dV);
        int a_np_tail = static_cast<int>(std::ceil(0.01*a_np)); // this is an approximation for now; saves 1% of SDs for the tail
        amrex::Vector<amrex::Real> tmp_mass(a_np);
        amrex::Vector<amrex::Real> tmp_mult(a_np);
        amrex::Print() << "Finding aerosol radius sampling range\n";
        while ((P_max >= 1.0 - tol) || (P_min <= tol)) {
            if (P_max >= 1.0 - tol) {
                rmax = rmax * 0.99;
            }
            if (P_min <= tol) {
                rmin = rmin * 1.01;
            }
            P_min = (1 + std::erf((std::log(rmin / mu)) / sigma / std::sqrt(2))) / 2;
            P_max = (1 + std::erf((std::log(rmax / mu)) / sigma / std::sqrt(2))) / 2;
        }
        dlnr = (std::log(rmax) - std::log(rmin));
        amrex::Print() << "Range: rmin =" << rmin << ", rmax = " << rmax << ", dlnr = " << dlnr << "\n";

        // initialize the main distribution
        amrex::Print() << "Initializing radii\n";
        auto lnrmin = std::log(rmin);
        for (int n = 0; n < a_np; n++) {
            auto tmp = lnrmin + urd(a_rng)*dlnr;
            auto dry_r = std::exp(tmp);
            tmp_mass[n] = (4.0/3.0) * PI * dry_r * dry_r * dry_r * a_density;
            auto term = std::exp(-std::log(dry_r/mu)*std::log(dry_r/mu)/(2.0*sigma*sigma));
            tmp_mult[n] =  (m_numdens * a_dV)/ (sigma*std::sqrt(2*PI)) * term;
        }

        // initialize the tail using approximate erfinv
        amrex::Print() << "Initializing tail: " << a_np_tail << " particles\n";
        auto tail_mult = std::exp(-std::log(rmax/mu)*std::log(rmax/mu)/(2.0*sigma*sigma)) / (sigma*std::sqrt(2*PI));
        for (int n = 0; n < a_np_tail; n++) {
            int sd_id = static_cast<int>(std::round(urd(a_rng) * a_np));
            auto tmp = P_max + (1.0 - P_max) * urd(a_rng);
            auto tmp2 = SD_erfinv(2 * tmp - 1);
            auto dry_r = mu * std::exp(sigma * std::sqrt(2) * tmp2);
            tmp_mass[sd_id] = (4.0/3.0) * PI * dry_r * dry_r * dry_r * a_density;
            // set the multiplicity to the same as for the 99th percentile aerosol
            tmp_mult[sd_id] = (m_numdens * a_dV) * tail_mult;
        }
        // Update SD multiplicity and mass with the initialized main + tail distribution
        for (int n = 0; n < a_np; n++) {
            a_mult[n] += tmp_mult[n];
            a_mass[n] += tmp_mass[n];
        }
        amrex::Print() << "Done sampling\n";
    } else {
        amrex::Abort("Unknown m_init_type!");
    }
}
