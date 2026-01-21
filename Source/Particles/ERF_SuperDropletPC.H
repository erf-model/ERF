#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include <random>
#include "ERF_Constants.H"
#include "ERF_MaterialProperties.H"
#include "ERF_SuperDropletPCDefinitions.H"
#include <AMReX_StructOfArrays.H>

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
            if (a_name == Species::Name::H2O) { m_idx_w = m_species_mat.size(); }
            if (a_name == Species::Name::ice) { m_idx_i = m_species_mat.size(); }
            m_species_mat.push_back(std::make_unique<MaterialProperties>(a_name));
            m_device_props_initialized = false;  // Properties have changed, need to reinitialize
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
            m_device_props_initialized = false;  // Properties have changed, need to reinitialize
        }

        /*! \brief Set aerosol material */
        AMREX_FORCE_INLINE
        virtual void setAerosolMaterial ( const std::vector<Species::Name>& a_names )
        {
            for (auto& name : a_names) { setAerosolMaterial(name); }
        }

        /*! \brief Get device pointer to species densities */
        const amrex::ParticleReal* getSpeciesDensitiesDevice() const;

        /*! \brief Get device pointer to species solubilities */
        const int* getSpeciesSolubilitiesDevice() const;

        /*! \brief Get device pointer to species ionization values */
        const amrex::ParticleReal* getSpeciesIonizationDevice() const;

        /*! \brief Get device pointer to species molecular weight values */
        const amrex::ParticleReal* getSpeciesMolWeightDevice() const;

        /*! \brief Get device pointer to aerosol densities */
        const amrex::ParticleReal* getAerosolDensitiesDevice() const;

        /*! \brief Get device pointer to aerosol solubilities */
        const int* getAerosolSolubilitiesDevice() const;

        /*! \brief Get device pointer to aerosol ionization values */
        const amrex::ParticleReal* getAerosolIonizationDevice() const;

        /*! \brief Get device pointer to aerosol molecular weight values */
        const amrex::ParticleReal* getAerosolMolWeightDevice() const;

        /*! \brief Get device pointer to species ice nucleating particle flags */
        const int* getSpeciesINPFlagsDevice() const;

        /*! \brief Get device pointer to aerosol ice nucleating particle flags */
        const int* getAerosolINPFlagsDevice() const;

        /*! \brief Update device properties if material properties change */
        void updateDeviceProperties();

        /*! \brief Setup species and aerosol mass pointers from SOA
         *
         * This helper method sets up the species and aerosol mass pointers from the
         * Structure of Arrays (SOA) data. It handles the common pattern of setting up
         * these pointers that appears in multiple places in the code.
         *
         * @param soa The Structure of Arrays containing particle data
         * @param sp_mass_ptrs Output array for species mass pointers
         * @param ae_mass_ptrs Output array for aerosol mass pointers
         */
        template<typename SOAType>
        void setupMassPointers(
            SOAType& soa,
            SDPCDefn::SDSpeciesMassArr& sp_mass_ptrs,
            SDPCDefn::SDAerosolMassArr& ae_mass_ptrs) const
        {
            // Setup species mass pointers
            for (int i = 0; i < m_num_species; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i, m_num_aerosols, m_num_species)).data();
            }

            // Setup aerosol mass pointers
            for (int i = 0; i < m_num_aerosols; i++) {
                ae_mass_ptrs[i] = soa.GetRealData(idx_a(i, m_num_aerosols, m_num_species)).data();
            }
        }

        /*! \brief Helper method to update common particle attributes (radius and mass)
         *
         * @param particle_idx Index of the particle to update
         * @param radius_ptr Pointer to the radius array
         * @param mass_ptr Pointer to the mass array
         * @param idx_w Index of water species
         * @param rho_w Density of water
         * @param num_sp Number of species
         * @param num_ae Number of aerosols
         * @param sp_sol_arr Species solubility array
         * @param ae_sol_arr Aerosol solubility array
         * @param sp_mass_ptrs Species mass pointers
         * @param ae_mass_ptrs Aerosol mass pointers
         * @param sp_rho_arr Species density array
         * @param ae_rho_arr Aerosol density array
         */
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        static void updateParticleAttributes(
            int particle_idx,
            amrex::ParticleReal* radius_ptr,
            amrex::ParticleReal* mass_ptr,
            int idx_w,
            amrex::ParticleReal rho_w,
            int num_sp, int num_ae,
            const int* sp_sol_arr,
            const int* ae_sol_arr,
            const SDPCDefn::SDSpeciesMassArr sp_mass_ptrs,
            const SDPCDefn::SDAerosolMassArr ae_mass_ptrs,
            const amrex::ParticleReal* sp_rho_arr,
            const amrex::ParticleReal* ae_rho_arr)
        {
            // Update effective radius
            radius_ptr[particle_idx] = SD_effective_radius(
                particle_idx, idx_w, rho_w,
                num_sp, num_ae,
                sp_sol_arr, ae_sol_arr,
                sp_mass_ptrs, ae_mass_ptrs,
                sp_rho_arr, ae_rho_arr);

            // Update total mass
            mass_ptr[particle_idx] = SD_total_mass(
                particle_idx, num_sp, num_ae,
                sp_mass_ptrs, ae_mass_ptrs);
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
         * \param[in] a_lev AMR level
         * \param[in] a_time Current simulation time
         * \param[in] a_dt Timestep for advection
         * \param[in] a_flow_vel Array of face-based velocities
         * \param[in] a_density Density field
         * \param[in] a_pressure Pressure field
         * \param[in] a_temperature Temperature field
         * \param[in] a_z_phys_nd Array of terrain heights
         * \param[in] a_bctypes Array of boundary condition types
         * \param[in] a_recycle Flag to enable particle recycling
         */
        virtual void AdvectParticles (   int a_lev,
                                         amrex::Real a_time,
                                         amrex::Real a_dt,
                                         const amrex::MultiFab* const a_flow_vel,
                                         const amrex::MultiFab& a_density,
                                         const amrex::MultiFab& a_pressure,
                                         const amrex::MultiFab& a_temperature,
                                         const amrex::Vector<MFPtr>& a_z_phys_nd,
                                         const BCTypeArr& a_bctypes,
                                         const bool a_recycle);

        /*! \brief Condensation and evaporation growth/shrinking for one time step */
        virtual void MassChange_LV (int,
                                    amrex::Real,
                                    const Species::Name&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::Vector<MFPtr>&,
                                    const bool );

        /*! \brief Freezing/melting for one time step */
        virtual void MassChange_SL (int,
                                    amrex::Real,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::Vector<MFPtr>& );

        /*! \brief Deposition and sublimation growth/shrinking for one time step */
        virtual void MassChange_SV (int,
                                    amrex::Real,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::Vector<MFPtr>& );

        /*! \brief Coalescence of super-droplets */
        virtual void Coalescence (  int,
                                    amrex::Real,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::Vector<MFPtr>& );

        /*! \brief Recycle super-droplets
         * \param[in] a_lev AMR level
         * \param[in] a_z_phys_nd Array of terrain heights
         * \param[in] a_iter Current iteration number
         * \param[in] a_dt Timestep size
         * \param[in] a_recycle Flag to enable recycling
         */
        virtual void Recycle ( const int a_lev,
                               const amrex::Vector<MFPtr>& a_z_phys_nd,
                               const int a_iter,
                               const amrex::Real a_dt,
                               const bool a_recycle);

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

        /*! \brief Compute the species mass density on a mesh */
        virtual void speciesMassDensity ( amrex::MultiFab&,
                                          int,
                                          const int a_comp = 0) const;

        /*! \brief Compute the cloud/rain mass density on a mesh */
        virtual void cloudRainDensity ( amrex::MultiFab&,
                                        amrex::Real,
                                        amrex::Real,
                                        const int a_comp = 0) const;

        /*! \brief Compute the ice mass density on a mesh */
        virtual void iceDensity ( amrex::MultiFab&,
                                  amrex::Real,
                                  const int a_comp = 0) const;

        /*! \brief Compute the snow mass density on a mesh */
        virtual void snowDensity ( amrex::MultiFab&,
                                   amrex::Real,
                                   const int a_comp = 0) const;

        /*! \brief Compute the graupel mass density on a mesh */
        virtual void graupelDensity ( amrex::MultiFab&,
                                      amrex::Real,
                                      const int a_comp = 0) const;

        /*! \brief Compute the frozen water mass density on a mesh */
        virtual void totalIceDensity ( amrex::MultiFab&,
                                       const int a_comp = 0) const;

        /*! \brief Compute the species mass flux component on a mesh */
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

        SDTerminalVelocityType m_term_vel_type_w; /*!< Terminal velocity model for water droplets */
        SDTerminalVelocityType m_term_vel_type_i; /*!< Terminal velocity model for ice particles */

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
        int m_idx_w = -1;

        /*! species index of ice */
        int m_idx_i = -1;

        /*! inactive particles threshold */
        amrex::Real m_deac_threshold;

        /*! save inactive particles to file */
        bool m_save_inactive;

        /*! Persistent device vectors for material properties */
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sp_density;    // Species densities
        amrex::Gpu::DeviceVector<int> m_sp_solubility;                 // Species solubilities
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sp_ionization; // Species ionization
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sp_mol_weight; // Species molecular weight
        amrex::Gpu::DeviceVector<int> m_sp_is_INP;                     // Species ice nucleating particle flags

        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_ae_density;    // Aerosol densities
        amrex::Gpu::DeviceVector<int> m_ae_solubility;                 // Aerosol solubilities
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_ae_ionization; // Aerosol ionization
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_ae_mol_weight; // Aerosol molecular weight
        amrex::Gpu::DeviceVector<int> m_ae_is_INP;                     // Aerosol ice nucleating particle flags

        /*! Flag to track if device properties are initialized */
        bool m_device_props_initialized = false;

        /*! Method to initialize device properties */
        void initializeDeviceProperties();

        /*! \brief Helper method to get species and aerosol material properties
         *
         * This helper method provides pointers to species and aerosol material properties
         * from the persistent device vectors, ensuring they are initialized first.
         * It avoids redundant host-to-device transfers when accessing material properties.
         *
         * @param sp_density_ptr Output pointer for species densities
         * @param sp_solubility_ptr Output pointer for species solubilities
         * @param ae_density_ptr Output pointer for aerosol densities
         * @param ae_solubility_ptr Output pointer for aerosol solubilities
         */
        AMREX_FORCE_INLINE
        void getMaterialPropertiesDevice(
            const amrex::ParticleReal*& sp_density_ptr,
            const int*& sp_solubility_ptr,
            const amrex::ParticleReal*& ae_density_ptr,
            const int*& ae_solubility_ptr) const
        {
            // Ensure device properties are initialized
            if (!m_device_props_initialized) {
                const_cast<SuperDropletPC*>(this)->initializeDeviceProperties();
            }

            // Get pointers to persistent device data
            sp_density_ptr = getSpeciesDensitiesDevice();
            sp_solubility_ptr = getSpeciesSolubilitiesDevice();
            ae_density_ptr = getAerosolDensitiesDevice();
            ae_solubility_ptr = getAerosolSolubilitiesDevice();
        }

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
