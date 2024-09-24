#include "ERF_SDInitialization.H"

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

    m_num_aerosols = a_aerosol_mat.size();
    m_mass_aerosol_min.resize(m_num_aerosols);
    m_mass_aerosol_mean.resize(m_num_aerosols);
    m_radius_aerosol_min.resize(m_num_aerosols);
    m_radius_aerosol_max.resize(m_num_aerosols);
    m_aerosol_init_type.resize(m_num_aerosols);

    for (int i = 0; i < m_num_aerosols; i++) {
        // default values
        m_aerosol_init_type[i] = SupDropInit::attrib_init_const;
        m_mass_aerosol_min[i]   = 0.0;
        m_mass_aerosol_mean[i]  = 0.0;
        m_radius_aerosol_min[i] = 1.0e-40;
        m_radius_aerosol_max[i] = 1.0e-40;

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
            std::string key = "initial_aerosol_min_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_min[i]);
        }
        {
            std::string key = "initial_aerosol_max_radius_" + a_aerosol_mat[i]->name();
            pp.query(key.c_str(), m_radius_aerosol_max[i]);
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
                        << ", mean=" << m_mass_aerosol_mean[i];
            } else if (m_aerosol_init_type[i] == SupDropInit::attrib_init_lnr) {
                Print() << ", min=" << m_radius_aerosol_min[i]
                        << ", max=" << m_radius_aerosol_max[i];
            }
            Print() << ")" << "\n";
        }
    }

    Print() << "    Number of seed particles per cell: " << m_ppc_seed << "\n"
            << "    Seed condensate mass: " << m_seed_mass << "\n";

}

void SDInitialization::getAerosolDistribution ( amrex::Vector<amrex::Real>& a_aerosol_mass,
                                                const int a_idx,
                                                const int a_np,
                                                const amrex::Real a_density ) const
{
    a_aerosol_mass.resize(a_np);
    if (m_aerosol_init_type[a_idx] == SupDropInit::attrib_init_const) {
        for (int n = 0; n < a_np; n++) {
            a_aerosol_mass[n] = m_mass_aerosol_mean[a_idx];
        }
    } else if (m_aerosol_init_type[a_idx] == SupDropInit::attrib_init_exp) {
        std::random_device rd;
        std::mt19937 rng(rd());
        auto delta = m_mass_aerosol_mean[a_idx] - m_mass_aerosol_min[a_idx];
        std::exponential_distribution<amrex::Real> ed(1.0/delta);
        for (int n = 0; n < a_np; n++) {
            a_aerosol_mass[n] = ed(rng) + m_mass_aerosol_min[a_idx];
        }
    } else if (m_aerosol_init_type[a_idx] == SupDropInit::attrib_init_lnr) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<> urd(0.0, 1.0);
        auto delta =   std::log(m_radius_aerosol_max[a_idx])
                     - std::log(m_radius_aerosol_min[a_idx]);
        for (int n = 0; n < a_np; n++) {
            auto term = std::log(m_radius_aerosol_min[a_idx]) + urd(rng)*delta;
            auto dry_r = std::exp(term);
            a_aerosol_mass[n] = (4.0/3.0) * PI
                                * dry_r * dry_r * dry_r
                                * a_density;
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

