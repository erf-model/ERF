#include <cmath>
#include "ERF_SDInitialization.H"

void SDInitialization::setDefaults ( const MatVec& a_species_mat,
                                     const MatVec& a_aerosol_mat )
{
    BL_PROFILE("SDInitialization::setDefaults");

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
        m_radius_species_geom_std[i] = 1.0;
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
        m_radius_aerosol_geom_std[i] = 1.0;
    }

    m_mult_type = SDMultiplicityType::sampled;
}

void SDInitialization::readInputs ( const std::string& a_prefix,
                                    const amrex::Geometry& a_geom,
                                    const MatVec& a_species_mat,
                                    const MatVec& a_aerosol_mat )
{
    BL_PROFILE("SDInitialization::readInputs");

    amrex::ParmParse pp(a_prefix);
    pp.query("initial_distribution_type", m_type);
    pp.query("initial_particles_per_cell", m_ppc_init);
    pp.query("initial_number_density", m_numdens_init);
    pp.query("initial_super_droplet_density", m_numdens_sd_init);
    pp.query("maximum_multiplicity", m_max_multiplicity);

    pp.query("multiplicity_type", m_mult_type);

    if (m_type == SupDropInit::init_uniform) {
        amrex::Vector<amrex::Real> particle_box_lo(AMREX_SPACEDIM);
        amrex::Vector<amrex::Real> particle_box_hi(AMREX_SPACEDIM);

        // Defaults
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            particle_box_lo[i] = a_geom.ProbLo(i);
            particle_box_hi[i] = a_geom.ProbHi(i);
        }

        pp.queryAdd("particle_box_lo", particle_box_lo, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_box_lo.size() == AMREX_SPACEDIM);

        pp.queryAdd("particle_box_hi", particle_box_hi, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_box_hi.size() == AMREX_SPACEDIM);

        m_init_particle_box.setLo(particle_box_lo);
        m_init_particle_box.setHi(particle_box_hi);
    } else if (m_type == SupDropInit::init_bubble){
        amrex::Vector<amrex::Real> particle_bubble_center(AMREX_SPACEDIM);
        amrex::Vector<amrex::Real> particle_bubble_radius(AMREX_SPACEDIM);

        // Defaults
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            particle_bubble_center[i] = a_geom.ProbHi(i)/2;
            particle_bubble_radius[i] = a_geom.ProbHi(i)/2;
        }

        pp.queryAdd("particle_bubble_center", particle_bubble_center, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_bubble_center.size() == AMREX_SPACEDIM);

        pp.queryAdd("particle_bubble_radius", particle_bubble_radius, AMREX_SPACEDIM);
        AMREX_ASSERT(particle_bubble_radius.size() == AMREX_SPACEDIM);

        m_init_particle_box.setLo(particle_bubble_radius);
        m_init_particle_box.setHi(particle_bubble_center);
    }

    // Backward compatibility
    for (int i = 0; i < m_num_species; i++) {
        pp.query("initial_condensate_distribution_type", m_species_init_type[i]);
        pp.query("initial_condensate_mass_min", m_mass_species_min[i]);
        pp.query("initial_condensate_mass_mean", m_mass_species_mean[i]);
        pp.query("initial_condensate_min_radius", m_radius_species_min[i]);
        pp.query("initial_condensate_max_radius", m_radius_species_max[i]);
    }
    for (int i = 0; i < m_num_species; i++) {
        {
            std::string key = "initial_species_distribution_type_"+a_species_mat[i]->name();
            pp.query(key.c_str(), m_species_init_type[i]);
        }
        {
            std::string key = "initial_species_min_mass_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_mass_species_min[i]);
        }
        {
            std::string key = "initial_species_mean_mass_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_mass_species_mean[i]);
        }
        {
            m_mass_species_max[i] = 5 * m_mass_species_mean[i]; // default
            std::string key = "initial_species_max_mass_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_mass_species_max[i]);
        }
        {
            std::string key = "initial_species_min_radius_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_radius_species_min[i]);
        }
        {
            std::string key = "initial_species_max_radius_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_radius_species_max[i]);
        }
        {
            m_radius_species_mean[i] = std::exp(0.5*(std::log(m_radius_species_min[i])+std::log(m_radius_species_max[i])));
            std::string key = "initial_species_mean_radius_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_radius_species_mean[i]);
        }
        {
            m_radius_species_geom_std[i] = 2.0;
            std::string key = "initial_species_std_radius_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_radius_species_geom_std[i]);
            m_radius_species_geom_std[i] = std::exp(m_radius_species_geom_std[i]);
        }
        {
            std::string key = "initial_species_geomstd_radius_" + a_species_mat[i]->name();
            pp.query(key.c_str(), m_radius_species_geom_std[i]);
        }
    }

    for (int i = 0; i < m_num_aerosols; i++) {
        {
            std::string key = "initial_aerosol_distribution_type_"+a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_aerosol_init_type[i]);
        }
        {
            std::string key = "initial_aerosol_min_mass_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_mass_aerosol_min[i]);
        }
        {
            std::string key = "initial_aerosol_mean_mass_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_mass_aerosol_mean[i]);
        }
        {
            m_mass_aerosol_max[i] = 5 * m_mass_aerosol_mean[i]; // default
            std::string key = "initial_aerosol_max_mass_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_mass_aerosol_max[i]);
        }
        {
            std::string key = "initial_aerosol_min_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_min[i]);
        }
        {
            std::string key = "initial_aerosol_max_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_max[i]);
        }
        {
            m_radius_aerosol_mean[i] = std::exp(0.5*(std::log(m_radius_aerosol_min[i])+std::log(m_radius_aerosol_max[i])));
            std::string key = "initial_aerosol_mean_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_mean[i]);
        }
        {
            m_radius_aerosol_geom_std[i] = 2.0;
            std::string key = "initial_aerosol_std_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_geom_std[i]);
            m_radius_aerosol_geom_std[i] = std::exp(m_radius_aerosol_geom_std[i]);
        }
        {
            std::string key = "initial_aerosol_geomstd_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_geom_std[i]);
        }
    }

}

void SDInitialization::printParameters ( const MatVec& a_species_mat,
                                         const MatVec& a_aerosol_mat ) const
{
    using namespace amrex;
    Print() << "    Initial particle box: " << m_init_particle_box << "\n"
            << "    Initial number density: " << m_numdens_init << "\n"
            << "    Inital super-droplets number density: " << m_numdens_sd_init << "\n"
            << "    Initial multiplicity type: " << amrex::getEnumNameString(m_mult_type) << "\n";

    Print() << "    Vapour/Condensate Species material:\n";
    for (unsigned long i=0; i < a_species_mat.size(); i++) {
        Print() << "        "
                << a_species_mat[i]->name()
                << " (Initial distribution: " << m_species_init_type[i];
        if (m_species_init_type[i] == SupDropInit::attrib_init_const) {
            Print() << ", value=" << m_mass_species_mean[i];
        } else if (m_species_init_type[i] == SupDropInit::attrib_init_exp) {
            Print() << ", min=" << m_mass_species_min[i]
                    << ", mean=" << m_mass_species_mean[i]
                    << ", max=" << m_mass_species_max[i];
        } else if (m_species_init_type[i] == SupDropInit::attrib_init_lnr) {
            Print() << ", min=" << m_radius_species_min[i]
                    << ", max=" << m_radius_species_max[i]
                    << ", mean=" << m_radius_species_mean[i]
                    << ", std=" << m_radius_species_geom_std[i];
        } else if (m_species_init_type[i] == SupDropInit::attrib_init_lnr_auto) {
            Print() << ", mean=" << m_radius_species_mean[i]
                    << ", std=" << m_radius_species_geom_std[i];
        }
        Print() << ")" << "\n";
    }

    if (a_aerosol_mat.size() > 0) {
        Print() << "    Aerosols material:\n";
        for (unsigned long i=0; i < a_aerosol_mat.size(); i++) {
            Print() << "        "
                    << a_aerosol_mat[i]->name()
                    << " (Initial distribution: " << m_aerosol_init_type[i];
            if (m_aerosol_init_type[i] == SupDropInit::attrib_init_const) {
                Print() << ", value=" << m_mass_aerosol_mean[i];
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_exp) {
                Print() << ", min=" << m_mass_aerosol_min[i]
                        << ", mean=" << m_mass_aerosol_mean[i]
                        << ", max=" << m_mass_aerosol_max[i];
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr) {
                Print() << ", min=" << m_radius_aerosol_min[i]
                        << ", max=" << m_radius_aerosol_max[i]
                        << ", mean=" << m_radius_aerosol_mean[i]
                        << ", std=" << m_radius_aerosol_geom_std[i];
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr_auto) {
                Print() << ", mean=" << m_radius_aerosol_mean[i]
                        << ", std=" << m_radius_aerosol_geom_std[i];
            }
            Print() << ")" << "\n";
        }
    }

    Print() << "    Number of seed particles per cell: " << m_ppc_seed << "\n"
            << "    Seed condensate mass: " << m_seed_mass << "\n";

}

void SDInitialization::getDistribution ( amrex::Vector<amrex::Real>& a_mass,
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
    amrex::Real term = std::log(1 - x * x);
    amrex::Real p1 = 2 / (PI * a) + term / 2.0;
    amrex::Real p2 = term / a;
    return std::sqrt(std::sqrt(p1 * p1 - p2) - p1);
}

void SDInitialization::getDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                         amrex::Vector<amrex::Real>& a_mult,
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
            a_mult[n] += std::exp(-a_mass[n] / delta);
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
            a_mult[n] += 1.0 / (sigma*std::sqrt(2*PI)*dry_r) * term;
        }
    } else if (a_init_type == SupDropInit::attrib_init_lnr_auto) {
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto sigma = std::log(a_radius_gstd);
        auto mu = a_radius_mean;
        // automatically find the min and max radius of superdroplets, using Dziekan & Pawlowska 2017
        auto rmin = 1e-9;
        auto rmax = 1e-3;
        auto dlnr = (std::log(rmax) - std::log(rmin)) / a_np;
        auto P_min = 0.0;
        auto P_max = 1.0;
        auto tol = 1.0 / m_numdens_init;
        int a_np_tail = 0.01 * a_np; // this is an approximation for now; saves 1% of SDs for the tail
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
            a_mass[n] = (4.0/3.0) * PI * dry_r * dry_r * dry_r * a_density;
            auto term = std::exp(-std::log(dry_r/mu)*std::log(dry_r/mu)/(2.0*sigma*sigma));
            a_mult[n] = 1.0 / (sigma*std::sqrt(2*PI)*dry_r) * term;
        }

        // initialize the tail using approximate erfinv
        amrex::Print() << "Initializing tail: " << a_np_tail << " particles\n";
        auto tail_mult = std::exp(-std::log(rmax/mu)*std::log(rmax/mu)/(2.0*sigma*sigma)) / (sigma*std::sqrt(2*PI)*rmax);
        for (int n = 0; n < a_np_tail; n++) {
            int sd_id = urd(a_rng) * a_np;
            auto tmp = P_max + (1.0 - P_max) * urd(a_rng);
            auto tmp2 = SD_erfinv(2 * tmp - 1);
            auto dry_r = mu * std::exp(sigma * std::sqrt(2) * tmp2);
            a_mass[sd_id] = (4.0/3.0) * PI * dry_r * dry_r * dry_r * a_density;
            // set the multiplicity to the same as for the 99th percentile aerosol
            a_mult[sd_id] = tail_mult;
        }
        amrex::Print() << "Done sampling\n";
    } else {
        amrex::Abort("Unknown m_init_type!");
    }
}
