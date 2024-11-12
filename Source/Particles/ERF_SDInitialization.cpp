#include "ERF_SDInitialization.H"

void SDInitialization::setDefaults ( const std::vector<std::unique_ptr<MaterialProperties>>& a_aerosol_mat )
{
    BL_PROFILE("SDInitialization::setDefaults");

    m_num_aerosols = a_aerosol_mat.size();
    m_mass_aerosol_min.resize(m_num_aerosols);
    m_mass_aerosol_max.resize(m_num_aerosols);
    m_mass_aerosol_mean.resize(m_num_aerosols);
    m_radius_aerosol_min.resize(m_num_aerosols);
    m_radius_aerosol_max.resize(m_num_aerosols);
    m_radius_aerosol_mean.resize(m_num_aerosols);
    m_radius_aerosol_std.resize(m_num_aerosols);
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
        m_radius_aerosol_std[i] = 1.0;
    }
}

void SDInitialization::readInputs ( const std::string& a_prefix,
                                    const amrex::Geometry& a_geom,
                                    const std::vector<std::unique_ptr<MaterialProperties>>& a_aerosol_mat )
{
    BL_PROFILE("SDInitialization::readInputs");

    amrex::ParmParse pp(a_prefix);
    pp.query("initial_distribution_type", m_type);
    pp.query("initial_particles_per_cell", m_ppc_init);
    pp.query("initial_number_density", m_numdens_init);
    pp.query("initial_super_droplet_density", m_numdens_sd_init);
    pp.query("maximum_multiplicity", m_max_multiplicity);
    pp.query("initial_condensate_distribution_type", m_condensate_init_type);
    pp.query("initial_condensate_mass", m_mass_condensate_mean);
    pp.query("initial_condensate_min_radius", m_radius_condensate_min);
    pp.query("initial_condensate_max_radius", m_radius_condensate_max);

    pp.query("initial_seeds_per_cell", m_ppc_seed);
    pp.query("seed_condensate_mass", m_seed_mass);

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
            m_radius_aerosol_std[i] = 2.0;
            std::string key = "initial_aerosol_std_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_std[i]);
        }
    }

}

void SDInitialization::printParameters ( const std::vector<std::unique_ptr<MaterialProperties>>& a_aerosol_mat ) const
{
    using namespace amrex;
    Print() << "    Initial particle box: " << m_init_particle_box << "\n"
            << "    Initial number density: " << m_numdens_init << "\n"
            << "    Iniital super-droplets number density: " << m_numdens_sd_init << "\n";

    Print() << "    Condensate initial distribution: " << m_condensate_init_type << " (";
    if (m_condensate_init_type == SupDropInit::attrib_init_const) {
        Print() << "value=" << m_mass_condensate_mean;
    } else if (m_condensate_init_type == SupDropInit::attrib_init_exp) {
        Print() << "mean=" << m_mass_condensate_mean;
    } else if (m_condensate_init_type == SupDropInit::attrib_init_lnr) {
        Print() << "min=" << m_radius_condensate_min
                << ", max=" << m_radius_condensate_max;
    }
    Print() << ")\n";

    if (a_aerosol_mat.size() > 0) {
        Print() << "    Aerosol materials:\n";
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
                        << ", std=" << m_radius_aerosol_std[i];
            }
            Print() << ")" << "\n";
        }
    }

    Print() << "    Number of seed particles per cell: " << m_ppc_seed << "\n"
            << "    Seed condensate mass: " << m_seed_mass << "\n";

}

void SDInitialization::getAerosolDistribution ( amrex::Vector<amrex::Real>& a_aerosol_mass,
                                                amrex::Vector<amrex::Real>& a_multiplicity,
                                                const int a_idx,
                                                const int a_np,
                                                const amrex::Real a_density ) const
{
    a_aerosol_mass.resize(a_np);
    AMREX_ALWAYS_ASSERT(a_multiplicity.size() == a_np);
    if (m_aerosol_init_type[a_idx] == SupDropInit::attrib_init_const) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<> urd(0.0, 1.0);
        for (int n = 0; n < a_np; n++) {
            a_aerosol_mass[n] = m_mass_aerosol_mean[a_idx];
            a_multiplicity[n] += urd(rng); // initially this will be a non-integer; later we will rescale to an integer.
        }
    } else if (m_aerosol_init_type[a_idx] == SupDropInit::attrib_init_exp) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto delta = m_mass_aerosol_mean[a_idx] - m_mass_aerosol_min[a_idx];
        auto lnrng = std::log(m_mass_aerosol_max[a_idx]) - std::log(m_mass_aerosol_min[a_idx]);
        auto lnmin = std::log(m_mass_aerosol_min[a_idx]);
        for (int n = 0; n < a_np; n++) {
            auto tmp = lnmin + urd(rng) * lnrng;
            a_aerosol_mass[n] = std::exp(tmp);
            a_multiplicity[n] += std::exp(-a_aerosol_mass[n] / delta);
        }
    } else if (m_aerosol_init_type[a_idx] == SupDropInit::attrib_init_lnr) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto sigma = m_radius_aerosol_std[a_idx];
        auto mu = m_radius_aerosol_mean[a_idx];
        auto lnrng = std::log(m_radius_aerosol_max[a_idx]) - std::log(m_radius_aerosol_min[a_idx]);
        auto lnmin = std::log(m_radius_aerosol_min[a_idx]);
        for (int n = 0; n < a_np; n++) {
            auto tmp = lnmin + urd(rng) * lnrng;
            auto dry_r = std::exp(tmp);
            a_aerosol_mass[n] = (4.0/3.0) * PI * dry_r * dry_r * dry_r * a_density;
            auto term = std::exp(-std::log(dry_r/mu)*std::log(dry_r/mu)/(2.0*sigma*sigma));
            a_multiplicity[n] += 1.0 / (sigma*std::sqrt(2*PI)*dry_r) * term;
        }
    } else {
        amrex::Abort("Unknown m_aerosol_init_type!");
    }
}

void SDInitialization::getCondensateDistribution ( amrex::Vector<amrex::Real>& a_condensate_mass,
                                                   const int a_np,
                                                   const amrex::Real a_density ) const
{
    a_condensate_mass.resize(a_np);
    if (m_condensate_init_type == SupDropInit::attrib_init_exp) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::exponential_distribution<amrex::Real> ed(1.0/m_mass_condensate_mean);
        for (int n = 0; n < a_np; n++) {
            a_condensate_mass[n] = ed(rng);
        }
    } else if (m_condensate_init_type == SupDropInit::attrib_init_const) {
        for (int n = 0; n < a_np; n++) {
            a_condensate_mass[n] = m_mass_condensate_mean;
        }
    } else if (m_condensate_init_type == SupDropInit::attrib_init_lnr) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto delta =    std::log(m_radius_condensate_max)
                     -  std::log(m_radius_condensate_min);
        for (int n = 0; n < a_np; n++) {
            auto term = std::log(m_radius_condensate_min) + urd(rng)*delta;
            auto radius = std::exp(term);
            a_condensate_mass[n] =    (4.0/3.0) * PI
                                    * radius * radius * radius
                                    * a_density;
        }
    } else {
        amrex::Abort("Unknown m_condensate_init_type!");
    }
}

