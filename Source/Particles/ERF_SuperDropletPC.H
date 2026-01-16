#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include <random>
#include "ERF_Constants.H"
#include "ERF_MaterialProperties.H"
#include "ERF_SuperDropletPCDefinitions.H"

class SuperDropletPC : public ERFPC
{
    using MFPtr = std::unique_ptr<amrex::MultiFab>;
    using BCTypeArr = amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>;

    public:

        /*! \brief Constructor */
        SuperDropletPC (  amrex::ParGDBBase* a_gdb,
                          const std::vector<Species::Name>& a_species_mat,
                          const std::vector<Species::Name>& a_aerosol_mat,
                          const amrex::Real a_dt,
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_gdb, a_name)
        {
            define( a_species_mat,
                    a_aerosol_mat,
                    a_gdb->ParticleBoxArray(m_lev),
                    a_gdb->ParticleDistributionMap(m_lev),
                    a_dt );
        }

        /*! \brief Constructor */
        SuperDropletPC (  const amrex::Geometry& a_geom,
                          const amrex::DistributionMapping& a_dmap,
                          const amrex::BoxArray& a_ba,
                          const std::vector<Species::Name>& a_species_mat,
                          const std::vector<Species::Name>& a_aerosol_mat,
                          const amrex::Real a_dt,
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_geom, a_dmap, a_ba, a_name)
        {
            define(a_species_mat, a_aerosol_mat, a_ba, a_dmap, a_dt);
        }

        /*! \brief Destructor */
        ~SuperDropletPC()
        {
            if (m_mass_change_logging) {
                fclose(m_mass_change_log);
            }
        }

        /*! \brief Set vapour species material */
        AMREX_FORCE_INLINE
        virtual void setSpeciesMaterial ( const Species::Name& a_name )
        {
            if (a_name == Species::Name::H2O) {
                m_idx_w = m_species_mat.size();
            }
            m_species_mat.push_back(std::make_unique<MaterialProperties>(a_name));
        }

        /*! \brief Set species material */
        AMREX_FORCE_INLINE
        virtual void setSpeciesMaterial ( const std::vector<Species::Name>& a_names )
        {
            for (auto& name : a_names) { setSpeciesMaterial(name); }
        }

        /*! \brief Get vapour material */
        AMREX_FORCE_INLINE
        virtual const MaterialProperties& getSpeciesMaterial(const Species::Name& a_name) const
        {
            for (auto& species : m_species_mat) {
                if (species->m_name == a_name) { return *species; }
            }
            amrex::Abort("SuperDropletPC::getSpeciesMaterial() - species not found");
            return *m_species_mat[0];
        }

        /*! \brief Set aerosol material */
        AMREX_FORCE_INLINE
        virtual void setAerosolMaterial ( const Species::Name& a_name )
        {
            m_aerosol_mat.push_back(std::make_unique<MaterialProperties>(a_name));
        }

        /*! \brief Set aerosol material */
        AMREX_FORCE_INLINE
        virtual void setAerosolMaterial ( const std::vector<Species::Name>& a_names )
        {
            for (auto& name : a_names) { setAerosolMaterial(name); }
        }

        /*! \brief Get real-type particle attribute names */
        [[nodiscard]] virtual amrex::Vector<std::string> varNames () const override;

        /*! \brief Get Eulerian plot variable names */
        [[nodiscard]] virtual amrex::Vector<std::string> meshPlotVarNames () const override;

        /*! \brief Compute mesh variable from particles */
        virtual void computeMeshVar( const std::string&,
                                     amrex::MultiFab&,
                                     const int ) const override;

        /*! \brief Initialize super-droplets in domain */
        virtual void InitializeParticles (const amrex::Real, const MFPtr& a_ptr) override;

        /*! \brief Inject super-droplets in domain */
        virtual void InjectParticles (const amrex::Real, const MFPtr& a_ptr, const amrex::Real);

        /*! \brief add particles for a given set of initialization parameters */
        void addParticles ( const MFPtr&, const SDInitProperties& );

        /*! \brief Set initial number of superdroplets per cell as a box with uniform density */
        void setNumSDBoxDistribution( amrex::iMultiFab&,
                                      const int,
                                      const MFPtr&,
                                      const amrex::RealBox&,
                                      const bool );

        /*! \brief Set initial number of superdroplets per cell as a bubble with uniform density */
        void setNumSDBubbleDistribution( amrex::iMultiFab&,
                                      const int,
                                      const MFPtr&,
                                      const amrex::RealBox&,
                                      const bool );

        /*! \brief Set super-droplets attributes from a given condensate mass density
         * \param[in] a_mf Multifab containing the initial condensate mass density
         */
        virtual void SetAttributes (  amrex::MultiFab& a_mf );

        /*! \brief Scale the SDM number density with air density
         * \param[in] a_mf Multifab containing the air density field
         */
        virtual void DensityScaling ( const amrex::MultiFab& a_mf );

        /*! \brief Evolve particles for one time step */
        virtual void EvolveParticles ( int,
                                       amrex::Real,
                                       amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                       const amrex::Vector<MFPtr>& ) override
        {
            amrex::Abort("SuperDropletPC::EvolveParticles() is intentionally disabled.");
        }

        /*! \brief Advect particles for one time step
         * \param[in] iteration Current iteration number
         * \param[in] dt Timestep for advection
         * \param[in] time Current simulation time
         * \param[in] moving_umac Array of face-based velocities for moving reference frame
         * \param[in] umac_x X-direction face-based velocity
         * \param[in] umac_y Y-direction face-based velocity
         * \param[in] umac_z Z-direction face-based velocity
         * \param[in] mf_array Array of additional multifabs
         * \param[in] bc_arr Array of boundary conditions
         * \param[in] verbose Flag to enable verbose output during advection
         */
        virtual void AdvectParticles (   int iteration,
                                         amrex::Real dt,
                                         amrex::Real time,
                                         const amrex::MultiFab* const moving_umac,
                                         const amrex::MultiFab& umac_x,
                                         const amrex::MultiFab& umac_y,
                                         const amrex::MultiFab& umac_z,
                                         const amrex::Vector<MFPtr>& mf_array,
                                         const BCTypeArr& bc_arr,
                                         const bool verbose);

        /*! \brief Grow/shrink particles for one time step
         * \param[in] iteration Current iteration number
         * \param[in] dt Timestep for mass change
         * \param[in] species_name Name of the species undergoing phase change
         * \param[in] qv Vapour mixing ratio field
         * \param[in] temperature Temperature field
         * \param[in] pressure Pressure field
         * \param[in] density Density field
         * \param[in] mf_array Array of additional multifabs
         * \param[in] verbose Flag to enable verbose output during mass change
         */
        virtual void MassChange (   int iteration,
                                    amrex::Real dt,
                                    const Species::Name& species_name,
                                    const amrex::MultiFab& qv,
                                    const amrex::MultiFab& temperature,
                                    const amrex::MultiFab& pressure,
                                    const amrex::MultiFab& density,
                                    const amrex::Vector<MFPtr>& mf_array,
                                    const bool verbose);

        /*! \brief Coalescence of super-droplets
         * \param[in] iteration Current iteration number
         * \param[in] dt Timestep for coalescence
         * \param[in] temperature Temperature field used in coalescence calculation
         * \param[in] density Density field used in coalescence calculation
         */
        virtual void Coalescence (  int iteration,
                                    amrex::Real dt,
                                    const amrex::MultiFab& temperature,
                                    const amrex::MultiFab& density);

        /*! \brief Recycle super-droplets
         * \param[in] iteration Current iteration number
         * \param[in] mf_array Array of additional multifabs
         * \param[in] species_idx Index of the species to recycle
         * \param[in] time Current simulation time
         * \param[in] verbose Flag to enable verbose output during recycling
         */
        virtual void Recycle ( const int iteration,
                               const amrex::Vector<MFPtr>& mf_array,
                               const int species_idx,
                               const amrex::Real time,
                               const bool verbose);

        /*! \brief Return the number of super-droplets
            This returns the number of AMReX particles being used in the simulation. */
        [[nodiscard]] inline virtual amrex::Long NumSuperDroplets ()
        {
            return ERFPC::TotalNumberOfParticles();
        }

        /*! \brief Return the total number of physical particles */
        [[nodiscard]] virtual amrex::Real TotalNumberOfParticles ();

        /*! \brief Return the number of deactivated super-droplets */
        [[nodiscard]] virtual amrex::Long NumSDDeactivated ();

        /*! \brief Compute and print diagnostics */
        virtual void Diagnostics (int, amrex::Real, bool);
        /*! \brief Compute mass density distribution */
        virtual void ComputeDistributions (  int,
                                             amrex::ParticleReal,
                                             amrex::ParticleReal );
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        /*! \brief Compute binned size distributions */
        virtual void ComputeBinnedDistributions ( int );
        /*! \brief Compute cell-wise binned size distributions */
        virtual void ComputeBinnedDistributionsCell ( int, amrex::Real );
#endif

        /*! \brief Compute the SD number density on a mesh */
        virtual void SDNumberDensity ( amrex::MultiFab&, const int a_comp=0 ) const;
        /*! \brief Compute the particle number density on a mesh */
        virtual void numberDensity ( amrex::MultiFab&, const int a_comp=0 ) const;
        /*! \brief Compute the particle mass density on a mesh */
        virtual void massDensity ( amrex::MultiFab&, const int a_comp=0 ) const;
        /*! \brief Compute the mass flux component on a mesh */
        virtual void massFlux ( amrex::MultiFab&, const int, const int a_comp=0 ) const;

        /*! \brief Compute the condensate mass density on a mesh */
        virtual void speciesMassDensity ( amrex::MultiFab& a_mf, const int a_i ) const
        {
            speciesMassDensity(a_mf, a_i, 0.0, DBL_MAX);
        }
        /*! \brief Compute the condensate mass density on a mesh */
        virtual void speciesMassDensity ( amrex::MultiFab&,
                                          const int,
                                          const amrex::Real,
                                          const amrex::Real,
                                          const int a_comp = 0) const;
        /*! \brief Compute the condensate mass flux component on a mesh */
        virtual void speciesMassFlux ( amrex::MultiFab&, const int, const int, const int a_comp=0 ) const;

        /*! \brief Compute the aerosol mass density on a mesh */
        virtual void aerosolMassDensity ( amrex::MultiFab&, const int, const int a_comp=0 ) const;
        /*! \brief Compute the aerosol mass flux component on a mesh */
        virtual void aerosolMassFlux ( amrex::MultiFab&, const int, const int, const int a_comp=0 ) const;

        /*! \brief Compute the particle effective radius on a mesh */
        virtual void effectiveRadius ( amrex::MultiFab&, const int a_comp=0 ) const;

        /*! \brief Applies boundary treatment to the particle container */
        virtual void applyBoundaryTreatment ( int, const amrex::Vector<MFPtr>&, const BCTypeArr&, const bool );

    protected:

        const int m_lev = 0; /*!< AMR Level: always 0 for super-droplets PC */

        SDCoalescenceKernelType m_coalescence_kernel; /*!< Coalescence kernel */
        bool m_include_brownian_coalescence; /*!< Include Brownian coalescence? */

        SDTerminalVelocityType m_term_vel_type; /*!< Terminal velocity model */

        int m_num_sd_per_cell;  /*!< Number of super-droplets per cell */
        bool m_density_scaling; /*!< Scale initial number density with air density */
        bool m_nucleate_particles;  /*!< nucleate new superdroplets from vapour */
        bool m_prescribed_advection; /*!< move particles with prescribed vertical velocity */

        int m_num_species; /*!< Number of vapour/condensate species */
        std::vector<std::unique_ptr<MaterialProperties>> m_species_mat; /*!< Vapour/condensate material */
        int m_num_aerosols; /*!< Number of aerosols */
        std::vector<std::unique_ptr<MaterialProperties>> m_aerosol_mat; /*!< Aerosol materials */

        /* Mass change equation Newton solver parameters */
        amrex::Real m_newton_rtol;  /*!< Newton solver - relative tolerance */
        amrex::Real m_newton_atol;  /*!< Newton solver - absolute tolerance */
        amrex::Real m_newton_stol;  /*!< Newton solver - step tolerance */
        int m_newton_maxits;        /*!< Newton solver - maximum iterations */

        /* Mass change equation (un)convergence logs */
        bool m_mass_change_logging; /*!< Whether to log unconverged particles */
        FILE* m_mass_change_log; /*!< File handle for log with unconverged info */
        std::string m_mass_change_log_fname; /*!< Name of unconverged info log file */

        amrex::Real m_mass_change_cfl; /*!< CFL for phase change equation */
        SDMassChangeTIMethod m_mass_change_ti; /*!< time integrator for mass change ODE */
        long m_num_unconverged_particles; /*!< total number of unconverged particles */

        amrex::IntVect m_coalescence_bin_size; /*!< bin size for coalescence */

        int m_distribution_grid_size; /*!< Size of 1D grid to compute distributions on */

        amrex::MultiFab m_mf_buf;

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        amrex::Real m_bindist_rmin; /*!< min radius for binned distribution */
        amrex::Real m_bindist_rmax; /*!< max radius for binned distribution */

        /*! Cell-wise mass distribution */
        amrex::MultiFab m_mass_ln_R_mf;
        /*! Cell-wise number distribution */
        amrex::MultiFab m_num_ln_R_mf;
#endif

        /*! maximum coalescence time scale */
        amrex::Real m_t_coalescence;

        /*! Initialization parameters objects */
        int m_num_initializations = 1;
        /*! vector of initializations */
        std::vector< std::unique_ptr<SDInitialization> > m_initializations;

        /*! Injection parameters objects */
        int m_num_injections = 0;
        /*! vector of injections */
        std::vector< std::unique_ptr<SDInjection> > m_injections;

        /*! sigma_0 parameter for mass distribution calculation */
        amrex::Real m_sigma0;

        /*! place the initial SDs randomly within cells */
        bool m_place_randomly_in_cells;

        /*! random engine */
        std::mt19937 m_rndeng;

        /*! species index of water */
        int m_idx_w;

        /*! inactive particles threshold */
        amrex::Real m_deac_threshold;

        /*! save inactive particles to file */
        bool m_save_inactive;

        /* recycled particle position bounds */
        amrex::Real m_recyc_xmin;
        amrex::Real m_recyc_xmax;
        amrex::Real m_recyc_ymin;
        amrex::Real m_recyc_ymax;
        amrex::Real m_recyc_zmin;
        amrex::Real m_recyc_zmax;

        /*! \brief read inputs from file */
        virtual void readInputs () override
        {
            amrex::Abort("SuperDropletPC::readInputs(): Do not use this interface.");
        }

        /*! \brief read inputs from file */
        virtual void readInputs (const amrex::Real);

        /*! \brief Particle initialization - null (no particles are initialized) */
        void initializeParticlesNull (const MFPtr&) { }

    private:

        /*! \brief define super-droplets */
        void define ( const std::vector<Species::Name>&,
                      const std::vector<Species::Name>&,
                      const amrex::BoxArray&,
                      const amrex::DistributionMapping&,
                      const amrex::Real );

        /*! \brief add super-droplet method-specific particle attributes */
        void add_superdroplet_attributes();

        /* do not use */
        using ERFPC::massDensity;

};

#endif
#endif
