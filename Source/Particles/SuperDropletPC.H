#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include "MaterialProperties.H"
#include "ERFPC.H"

/*! \brief Real-type array-of-struct attributes of a super-droplet */
using SuperDropletsRealIdxAoS = ERFParticlesRealIdxAoS;

/*! \brief Int-type array-of-struct attributes of a super-droplet */
using SuperDropletsIntIdxAoS = ERFParticlesIntIdxAoS;

/*! \brief Real-type struct-of-array attributes of a super-droplet */
using SuperDropletsRealIdxSoA = ERFParticlesRealIdxSoA;

/*! \brief Real-type struct-of-array attributes added during runtime */
struct SuperDropletsRealIdxSoA_RT
{
    enum {
        radius = 0,     /*!< Radius of physical particles */
        multiplicity,   /*!< Number of droplets that this super-droplet represents
                             (Real type to handle large values) */
        sd_mass,        /*!< Superdroplet mass */
        ncomps
    };
};

/*! \brief Int-type struct-of-array attributes of a super-droplet */
using SuperDropletsIntIdxSoA = ERFParticlesIntIdxSoA;

/*! \brief Utilities for super-droplet PC */
namespace SuperDropletsUtils
{
    using namespace amrex;

    /*! \brief RHS function for super-droplet radius ODE */
    struct dRsqdt_RHSFunc
    {
        ParticleReal a;
        ParticleReal b;
        ParticleReal L;
        ParticleReal K;
        ParticleReal Rv;
        ParticleReal rho_l;
        ParticleReal D;

        AMREX_GPU_HOST_DEVICE
        ParticleReal operator() (   const ParticleReal a_R_sq,
                                    const ParticleReal a_S,
                                    const ParticleReal a_T,
                                    const ParticleReal a_e_s,
                                    const ParticleReal a_M_s ) const noexcept
        {
            ParticleReal F_k = ( L/(Rv*a_T) - 1.0) * ((L*rho_l) / (K*a_T));
            ParticleReal F_d = (rho_l*Rv*a_T) / (D*a_e_s);

            ParticleReal alpha = 2.0 * (a_S-1.0) / (F_k + F_d);
            ParticleReal beta = -2.0 * (a/a_T) / (F_k + F_d);
            ParticleReal gamma = 2.0 * b * a_M_s / (F_k + F_d);

            ParticleReal R_inv = 1.0/std::sqrt(a_R_sq);
            ParticleReal R_inv_cubed = R_inv*R_inv*R_inv;

            ParticleReal retval = alpha + beta*R_inv + gamma*R_inv_cubed;
            return retval;
        }

    };

    /*! \brief Jacobian for super-droplet radius ODE */
    struct dRsqdt_RHSJac
    {
        ParticleReal a;
        ParticleReal b;
        ParticleReal L;
        ParticleReal K;
        ParticleReal Rv;
        ParticleReal rho_l;
        ParticleReal D;

        AMREX_GPU_HOST_DEVICE
        ParticleReal operator() (   const ParticleReal a_R_sq,
                                    const ParticleReal a_T,
                                    const ParticleReal a_e_s,
                                    const ParticleReal a_M_s ) const noexcept
        {
            ParticleReal F_k = ( L/(Rv*a_T) - 1.0) * ((L*rho_l) / (K*a_T));
            ParticleReal F_d = (rho_l*Rv*a_T) / (D*a_e_s);

            ParticleReal beta = -2.0 * (a/a_T) / (F_k + F_d);
            ParticleReal gamma = 2.0 * b * a_M_s / (F_k + F_d);

            ParticleReal R_inv = 1.0/std::sqrt(a_R_sq);
            ParticleReal R_inv_3 = R_inv*R_inv*R_inv;
            ParticleReal R_inv_5 = R_inv_3*R_inv*R_inv;

            ParticleReal retval = -0.5 * (beta*R_inv_3 + 3.0*gamma*R_inv_5);
            return retval;
        }

    };

    /*! \brief Newton solver for super-droplet radius ODE */
    template<typename RHSFunc, typename JacFunc, typename RT>
    struct NewtonSolver
    {
        RHSFunc m_rhs;
        JacFunc m_jac;

        RT  m_rtol;
        RT  m_atol;
        RT  m_stol;
        int m_maxits;

        AMREX_GPU_HOST_DEVICE
        void operator()  (  RT&             a_u,
                            RT&             a_r,
                            const RT&       a_dt,
                            const RT&       a_S,
                            const RT&       a_T,
                            const RT&       a_e_s,
                            const RT&       a_M_s,
                            bool&           a_converged ) const
        {
            a_converged = false;
            RT res_norm0 = 0.0;

            for (int k = 0; k < m_maxits; k++) {
                RT residual = a_u - (a_r + a_dt * m_rhs( a_u, a_S, a_T, a_e_s, a_M_s ));
                RT res_norm_a = std::sqrt(residual*residual);

                if (k == 0) { res_norm0 = res_norm_a; }
                RT res_norm_r = res_norm_a / res_norm0;

                if (res_norm_a <= m_atol) {
                    a_converged = true;
                    break;
                }
                if (res_norm_r <= m_rtol) {
                    a_converged = true;
                    break;
                }

                RT slope = 1.0 - a_dt * m_jac( a_u, a_T, a_e_s, a_M_s );
                RT du = 0.0;
                if (slope*slope <  1.0e-40) {
                    du = 0.0;
                } else {
                    du = - residual / slope;
                }

                RT du_norm = std::sqrt(du*du);
                if (du_norm <= m_stol) {
                    a_converged = true;
                    break;
                }

                if (a_u + du < 1e-16) {
                    a_u += 0.01*a_u;
                } else {
                    a_u += du;
                }
            }
        }
    };
}

/*! \brief List of super-droplet initializations */
namespace SupDropInit
{
    /*! Uniform distribution throughout domain */
    const std::string init_uniform = "uniform";
    /*! No particles are initialized; nucleation generates them */
    const std::string init_null = "null";

    /*! Maximum number of aerosols */
    const int num_aerosols_max = 8;

    const std::string attrib_init_const = "constant";
    const std::string attrib_init_exp = "exponential";

    const std::string alg_coalescence_dsmc = "DSMC";
    const std::string alg_coalescence_cudsmc = "cuDSMC";
}

using SDAerosolMassArr = amrex::GpuArray<amrex::ParticleReal*,SupDropInit::num_aerosols_max>;
using MFPtr = std::unique_ptr<amrex::MultiFab>;

class SuperDropletPC : public ERFPC
{
    public:

        /*! \brief Constructor */
        SuperDropletPC (  amrex::ParGDBBase* a_gdb,
                          const std::string& a_vap_mat,
                          const std::vector<std::string>& a_aerosol_mat,
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_gdb, a_name)
        {
            define(a_vap_mat, a_aerosol_mat);
        }

        /*! \brief Constructor */
        SuperDropletPC (  const amrex::Geometry& a_geom,
                          const amrex::DistributionMapping& a_dmap,
                          const amrex::BoxArray& a_ba,
                          const std::string& a_vap_mat,
                          const std::vector<std::string>& a_aerosol_mat,
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_geom, a_dmap, a_ba, a_name)
        {
            define(a_vap_mat, a_aerosol_mat);
        }

        /*! \brief Set vapour material */
        virtual inline void setVapourMaterial ( const std::string& a_name )
        {
            m_vapour_mat = std::make_unique<MaterialProperties>(a_name);
        }

        /*! \brief Get vapour material */
        virtual inline MaterialProperties& getVapourMaterial() const { return (*m_vapour_mat); }

        /*! \brief Set aerosol material */
        virtual inline void setAerosolMaterial ( const std::string& a_name )
        {
            m_aerosol_mat.push_back(std::make_unique<MaterialProperties>(a_name));
        }

        /*! \brief Set aerosol material */
        virtual inline void setAerosolMaterial ( const std::vector<std::string>& a_names )
        {
            for (unsigned long i = 0; i < a_names.size(); i++) {
                setAerosolMaterial( a_names[i] );
            }
        }

        /*! \brief Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const override
        {
            BL_PROFILE("ERFPCPC::varNames()");
            amrex::Vector<std::string> retval = {   AMREX_D_DECL("xvel","yvel","zvel"),
                                                    "particle_mass",
                                                    "radius",
                                                    "multiplicity",
                                                    "superdroplet_mass"};
            for (int i = 0; i < m_num_aerosols; i++) {
                retval.push_back(std::string("aerosol_mass_"+m_aerosol_mat[i]->name()));
            }
            return retval;
        }

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> meshPlotVarNames () const override
        {
            BL_PROFILE("ERFPCPC::varNames()");
            amrex::Vector<std::string> retval = { AMREX_D_DECL("mass_flux_x",
                                                               "mass_flux_y",
                                                               "mass_flux_z"),
                                                  "number_density",
                                                  "mass_density",
                                                  ("mass_density_"+m_vapour_mat->name()),
                                                  AMREX_D_DECL (
                                                      ("mass_flux_x_"+m_vapour_mat->name()),
                                                      ("mass_flux_y_"+m_vapour_mat->name()),
                                                      ("mass_flux_z_"+m_vapour_mat->name()) ) };
            for (int i = 0; i < m_num_aerosols; i++) {
                retval.push_back(std::string("aerosol_mass_density_"+m_aerosol_mat[i]->name()));
                retval.push_back(std::string("aerosol_mass_flux_x_"+m_aerosol_mat[i]->name()));
                retval.push_back(std::string("aerosol_mass_flux_y_"+m_aerosol_mat[i]->name()));
                retval.push_back(std::string("aerosol_mass_flux_z_"+m_aerosol_mat[i]->name()));
            }
            return retval;
        }

        /*! Compute mesh variable from particles */
        virtual void computeMeshVar(    const std::string&  a_var_name,
                                        amrex::MultiFab&    a_mf,
                                        const int           a_lev) const override
        {
            a_mf.setVal(0.0);
            if (a_lev == 0) {
                if (a_var_name == "number_density") {
                    numberDensity( a_mf );
                } else if (a_var_name == "mass_density") {
                    massDensity( a_mf );
                } else if (a_var_name == ("mass_density_"+m_vapour_mat->name())) {
                    massDensityCondensate( a_mf );
                } else if (a_var_name == ("mass_flux_x_"+m_vapour_mat->name())) {
                    massFluxCondensate( a_mf, 0 );
                } else if (a_var_name == ("mass_flux_y_"+m_vapour_mat->name())) {
                    massFluxCondensate( a_mf, 1 );
                } else if (a_var_name == ("mass_flux_z_"+m_vapour_mat->name())) {
                    massFluxCondensate( a_mf, 2 );
                } else if (a_var_name == "mass_flux_x") {
                    massFlux( a_mf, 0 );
                } else if (a_var_name == "mass_flux_y") {
                    massFlux( a_mf, 1 );
                } else if (a_var_name == "mass_flux_z") {
                    massFlux( a_mf, 2 );
                } else {
                    for (int i = 0; i < m_num_aerosols; i++) {
                        std::string var_name = "aerosol_mass_density_"+m_aerosol_mat[i]->name();
                        if (a_var_name == var_name) {
                            aerosolMassDensity( a_mf, i );
                            break;
                        }
                    }
                    for (int i = 0; i < m_num_aerosols; i++) {
                        std::string var_name = "aerosol_mass_flux_x_"+m_aerosol_mat[i]->name();
                        if (a_var_name == var_name) {
                            aerosolMassFlux( a_mf, i, 0 );
                            break;
                        }
                    }
                    for (int i = 0; i < m_num_aerosols; i++) {
                        std::string var_name = "aerosol_mass_flux_y_"+m_aerosol_mat[i]->name();
                        if (a_var_name == var_name) {
                            aerosolMassFlux( a_mf, i, 1 );
                            break;
                        }
                    }
                    for (int i = 0; i < m_num_aerosols; i++) {
                        std::string var_name = "aerosol_mass_flux_z_"+m_aerosol_mat[i]->name();
                        if (a_var_name == var_name) {
                            aerosolMassFlux( a_mf, i, 2 );
                            break;
                        }
                    }
                }
            }
        }

        /*! \brief Initialize super-droplets in domain */
        inline virtual void InitializeParticles (const MFPtr& a_ptr = nullptr) override
        {
            InitializeParticles( m_initialization_type, a_ptr );
        }

        /*! \brief Initialize super-droplets in domain */
        virtual void InitializeParticles (  const std::string&,
                                            const MFPtr& a_ptr = nullptr);

        /*! \brief Particle initialization - uniform distribution */
        inline void initializeParticlesUniformDistribution (const MFPtr& a_height = nullptr)
        {
            initializeParticlesUniformDistribution(a_height, m_init_particle_box);
        }
        /*! \brief Particle initialization - uniform distribution */
        void initializeParticlesUniformDistribution (const MFPtr&,
                                                     const amrex::RealBox& );

        /*! \brief Set super-droplets attributes from a given condensate mass density */
        virtual void SetAttributes (  amrex::MultiFab& );

        /*! \brief Evolve particles for one time step */
        virtual void EvolveParticles ( int,
                                       amrex::Real,
                                       amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                       const amrex::Vector<MFPtr>& ) override
        {
            amrex::Abort("SuperDropletPC::EvolveParticles() is intentionally disabled.");
        }

        /*! \brief Advect particles for one time step */
        virtual void AdvectParticles (   int,
                                         amrex::Real,
                                         amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                         const amrex::Vector<MFPtr>& );

        /*! \brief Grow/shrink particles for one time step */
        virtual void MassChange (   int,
                                    amrex::Real,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::Vector<MFPtr>& );

        /*! \brief Coalescence of super-droplets */
        virtual void Coalescence ( int, amrex::Real, const amrex::MultiFab& );

        /*! \brief Return the number of super-droplets
            This returns the number of AMReX particles being used in the simulation. */
        inline virtual amrex::Long NumSuperDroplets ()
        {
            return ERFPC::TotalNumberOfParticles();
        }

        /*! \brief Return the total number of physical particles */
        virtual amrex::Real TotalNumberOfParticles ();

        /*! \brief Compute and print diagnostics */
        virtual void Diagnostics (const int&, const bool);
        /*! \brief Compute mass density distribution */
        virtual void MassDensityDistribution (  const int&,
                                                const amrex::ParticleReal&,
                                                const amrex::ParticleReal& );

        /*! \brief Compute the particle number density on a mesh */
        virtual void numberDensity ( amrex::MultiFab&, const int& a_comp=0 ) const;
        /*! \brief Compute the particle mass density on a mesh */
        virtual void massDensity ( amrex::MultiFab&, const int& a_comp=0 ) const;
        /*! \brief Compute the condensate mass density on a mesh */
        virtual void massDensityCondensate ( amrex::MultiFab&, const int& a_comp=0 ) const;
        /*! \brief Compute the mass flux component on a mesh */
        virtual void massFlux ( amrex::MultiFab&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the condensate mass flux component on a mesh */
        virtual void massFluxCondensate ( amrex::MultiFab&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the aerosol mass density on a mesh */
        virtual void aerosolMassDensity ( amrex::MultiFab&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the aerosol mass flux component on a mesh */
        virtual void aerosolMassFlux ( amrex::MultiFab&, const int&, const int&, const int& a_comp=0 ) const;


    protected:

        const int m_lev = 0; /*!< AMR Level: always 0 for super-droplets PC */

        amrex::Real m_max_multiplicity;     /*!< Maximum multiplicity */
        amrex::Real m_numdens_init;         /*!< Initial number density of physical particles */
        amrex::Real m_numdens_sd_init;      /*!< Initial number density of super-droplets */

        std::string m_mass_condensate_init_type; /*!< type of initial distribution */
        amrex::Real m_mass_condensate_init; /*!< Initial condensate mass of each physical particle */
        std::string m_coalescence_alg; /*!< Choice of algorithm for coalescence */

        int m_ppc_seed; /*!< Number of seed particles per cell */
        amrex::Real m_seed_mass; /*!< Condensate mass of seed */

        /*! Initial aerosol masses of each physical particle */
        std::vector<amrex::Real> m_mass_aerosol_init;
        /*! Initial distribution type for aerosol masses */
        std::vector<std::string> m_mass_aerosol_init_type;
        int m_num_aerosols; /*!< Number of aerosols */

        int m_num_sd_per_cell;  /*!< Number of super-droplets per cell */
        bool m_nucleate_particles;  /*!< nucleate new superdroplets from vapour */
        amrex::RealBox m_init_particle_box; /*!<box within which to place particles*/

        std::unique_ptr<MaterialProperties> m_vapour_mat; /*!< Vapour material */
        std::vector<std::unique_ptr<MaterialProperties>> m_aerosol_mat; /*!< Aerosol materials */

        /* Newton solver parameters */
        amrex::Real m_newton_rtol;  /*!< Newton solver - relative tolerance */
        amrex::Real m_newton_atol;  /*!< Newton solver - absolute tolerance */
        amrex::Real m_newton_stol;  /*!< Newton solver - step tolerance */
        int m_newton_maxits;        /*!< Newton solver - maximum iterations */

        int m_mass_change_max_substeps; /*!< Max substeps for mass change update */
        int m_max_substeps_actual; /*!< Actual max substeps taken by mass change function */
        long m_num_unconverged_particles; /*!< total number of unconverged particles */

        amrex::IntVect m_coalescence_bin_size; /*!< bin size for coalescence */

        int m_distribution_grid_size; /*!< Size of 1D grid to compute distributions on */

        /*! \brief read inputs from file */
        virtual void readInputs ();

        /*! \brief Particle initialization - null (no particles are initialized) */
        void initializeParticlesNull (const MFPtr&) { }

    private:

        /*! \brief define super-droplets */
        inline void define ( const std::string& a_vap_mat,
                             const std::vector<std::string>& a_aerosol_mat)
        {
            m_num_unconverged_particles = 0;
            m_max_substeps_actual = 1;

            m_aerosol_mat.clear();

            setVapourMaterial( a_vap_mat );
            setAerosolMaterial( a_aerosol_mat );
            m_num_aerosols = m_aerosol_mat.size();
            m_mass_aerosol_init.resize(m_num_aerosols);
            m_mass_aerosol_init_type.resize(m_num_aerosols);
            AMREX_ASSERT(m_num_aerosols <= SupDropInit::num_aerosols_max);

            add_superdroplet_attributes();
            readInputs();
        }

        /*! \brief add super-droplet method-specific particle attributes */
        void add_superdroplet_attributes();

        /* do not use */
        using ERFPC::massDensity;

};

#endif
#endif
