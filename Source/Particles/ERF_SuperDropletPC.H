#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include <random>
#include "ERF_Constants.H"
#include "ERF_MaterialProperties.H"
#include "ERF_SuperDropletPCDefinitions.H"
#include "ERF_SuperDropletPCProcess.H"
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
            m_species_mat.push_back(std::make_unique<MaterialProperties>(a_name));
            m_device_props_initialized = false;
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
            m_device_props_initialized = false;
        }

        /*! \brief Set aerosol material */
        AMREX_FORCE_INLINE
        virtual void setAerosolMaterial ( const std::vector<Species::Name>& a_names )
        {
            for (auto& name : a_names) { setAerosolMaterial(name); }
        }

        /*! \brief Update device properties if material properties change */
        void updateDeviceProperties();

        /*! \brief Setup species and aerosol mass pointers from SOA */
        template<typename SOAType>
        void setupMassPointers(SOAType& soa, SDPCDefn::SDSpeciesMassArr& sp_mass_ptrs,
                               SDPCDefn::SDAerosolMassArr& ae_mass_ptrs) const
        {
            for (int i = 0; i < m_num_species; i++) {
                sp_mass_ptrs[i] = soa.GetRealData(idx_s(i, m_num_aerosols, m_num_species)).data();
            }
            for (int i = 0; i < m_num_aerosols; i++) {
                ae_mass_ptrs[i] = soa.GetRealData(idx_a(i, m_num_aerosols, m_num_species)).data();
            }
        }

        /*! \brief Update particle radius and total mass from species/aerosol masses */
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
                                     const amrex::MultiFab&,
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

        /*! \brief Advect particles for one time step */
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

        /*! \brief Grow/shrink particles for one time step */
        virtual void MassChange (   int a_lev,
                                    amrex::Real a_dt,
                                    const Species::Name& a_vap_name,
                                    const amrex::MultiFab& a_temperature,
                                    const amrex::MultiFab& a_pressure,
                                    const amrex::MultiFab& a_sat_pressure,
                                    const amrex::MultiFab& a_sat_ratio,
                                    const amrex::Vector<MFPtr>& a_z_phys_nd,
                                    const bool a_is_water);

        /*! \brief Coalescence of super-droplets */
        virtual void Coalescence (  int a_lev,
                                    amrex::Real a_dt,
                                    const amrex::MultiFab& a_pressure,
                                    const amrex::MultiFab& a_temperature);

        /*! \brief Recycle super-droplets */
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
        virtual void SDNumberDensity ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int a_comp=0 ) const;
        /*! \brief Compute the particle number density on a mesh */
        virtual void numberDensity ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int a_comp=0 ) const;
        /*! \brief Compute the particle mass density on a mesh */
        virtual void massDensity ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int a_comp=0 ) const;
        /*! \brief Compute the mass flux component on a mesh */
        virtual void massFlux ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int, const int a_comp=0 ) const;

        /*! \brief Compute the species mass density on a mesh */
        virtual void speciesMassDensity ( amrex::MultiFab&,
                                          const amrex::MultiFab& a_z_phys_nd,
                                          int,
                                          const int a_comp = 0) const;

        /*! \brief Compute the cloud/rain mass density on a mesh */
        virtual void cloudRainDensity ( amrex::MultiFab&,
                                        const amrex::MultiFab& a_z_phys_nd,
                                        amrex::Real,
                                        amrex::Real,
                                        const int a_comp = 0) const;

        /*! \brief Compute the species mass flux component on a mesh */
        virtual void speciesMassFlux ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int, const int, const int a_comp=0 ) const;

        /*! \brief Compute the aerosol mass density on a mesh */
        virtual void aerosolMassDensity ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int, const int a_comp=0 ) const;
        /*! \brief Compute the aerosol mass flux component on a mesh */
        virtual void aerosolMassFlux ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int, const int, const int a_comp=0 ) const;

        /*! \brief Compute the particle effective radius on a mesh */
        virtual void effectiveRadius ( amrex::MultiFab&, const amrex::MultiFab& a_z_phys_nd, const int a_comp=0 ) const;

        /*! \brief Applies boundary treatment to the particle container */
        virtual void applyBoundaryTreatment ( int, const amrex::Vector<MFPtr>&, const BCTypeArr&, const bool );

    protected:

        const int m_lev = 0; /*!< AMR Level: always 0 for super-droplets PC */

        SDCoalescenceKernelType m_coalescence_kernel; /*!< Coalescence kernel */
        bool m_include_brownian_coalescence; /*!< Include Brownian coalescence? */

        SDTerminalVelocityType m_term_vel_type_w; /*!< Terminal velocity model for water droplets */

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

        /*! inactive particles threshold */
        amrex::Real m_deac_threshold;

        /*! save inactive particles to file */
        bool m_save_inactive;

        /*! Persistent device vectors for material properties */
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sp_density;
        amrex::Gpu::DeviceVector<int> m_sp_solubility;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sp_ionization;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_sp_mol_weight;

        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_ae_density;
        amrex::Gpu::DeviceVector<int> m_ae_solubility;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_ae_ionization;
        amrex::Gpu::DeviceVector<amrex::ParticleReal> m_ae_mol_weight;

        /*! Flag to track if device properties are initialized */
        bool m_device_props_initialized = false;

        /* recycled particle position bounds */
        amrex::Real m_recyc_xmin;
        amrex::Real m_recyc_xmax;
        amrex::Real m_recyc_ymin;
        amrex::Real m_recyc_ymax;
        amrex::Real m_recyc_zmin;
        amrex::Real m_recyc_zmax;

        /*! Method to initialize device properties */
        void initializeDeviceProperties();

        /*! \brief Build ProcessContext with geometry and species info */
        SDProcess::ProcessContext buildProcessContext(int a_lev) const
        {
            SDProcess::ProcessContext ctx;
            const amrex::Geometry& geom = m_gdb->Geom(a_lev);
            ctx.plo = geom.ProbLoArray();
            ctx.phi = geom.ProbHiArray();
            ctx.dxi = geom.InvCellSizeArray();
            ctx.dx = geom.CellSizeArray();
            ctx.domain = geom.Domain();
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                ctx.is_periodic[d] = geom.isPeriodic(d) ? 1 : 0;
            }
            const auto cell_size = geom.CellSize();
            ctx.cell_volume = AMREX_D_TERM(cell_size[0], *cell_size[1], *cell_size[2]);
            ctx.num_species = m_num_species;
            ctx.num_aerosols = m_num_aerosols;
            ctx.idx_water = m_idx_w;
            ctx.rho_water = m_species_mat[m_idx_w]->m_density;
            return ctx;
        }

        /*! \brief Setup particle attribute pointers for a tile */
        template<typename SOAType, typename AOSType>
        void setupParticlePointers(
            SOAType& soa,
            AOSType& aos,
            SDProcess::ParticlePointers& ptrs) const
        {
            using namespace SDPCDefn;

            constexpr int rtoff_i = SuperDropletsIntIdx::ncomps;
            constexpr int rtoff_r = SuperDropletsRealIdx::ncomps;

            ptrs.num_particles = aos.numParticles();
            ptrs.mass_ptr = soa.GetRealData(SuperDropletsRealIdx::mass).data();
            ptrs.radius_ptr = soa.GetRealData(rtoff_r + SuperDropletsRealIdxSoA_RT::radius).data();
            ptrs.active_ptr = soa.GetIntData(rtoff_i + SuperDropletsIntIdxSoA_RT::active).data();
            ptrs.v_ptr[0] = soa.GetRealData(SuperDropletsRealIdx::vx).data();
            ptrs.v_ptr[1] = soa.GetRealData(SuperDropletsRealIdx::vy).data();
            ptrs.v_ptr[2] = soa.GetRealData(SuperDropletsRealIdx::vz).data();
            ptrs.vterm_ptr = soa.GetRealData(rtoff_r + SuperDropletsRealIdxSoA_RT::term_vel).data();
            ptrs.mult_ptr = soa.GetRealData(rtoff_r + SuperDropletsRealIdxSoA_RT::multiplicity).data();
            setupMassPointers(soa, ptrs.sp_mass_ptrs, ptrs.ae_mass_ptrs);
            if (!m_device_props_initialized) {
                const_cast<SuperDropletPC*>(this)->initializeDeviceProperties();
            }
            ptrs.sp_rho_arr = m_sp_density.data();
            ptrs.sp_sol_arr = m_sp_solubility.data();
            ptrs.sp_ion_arr = m_sp_ionization.data();
            ptrs.sp_mw_arr = m_sp_mol_weight.data();
            ptrs.ae_rho_arr = m_ae_density.data();
            ptrs.ae_sol_arr = m_ae_solubility.data();
            ptrs.ae_ion_arr = m_ae_ionization.data();
            ptrs.ae_mw_arr = m_ae_mol_weight.data();
        }

        /*! \brief Common body for forEachParticleTile iterations */
        template<typename TileFunc>
        void forEachParticleTileBody(
            ParIterType& pti, int a_lev,
            const SDProcess::ProcessContext& ctx,
            TileFunc&& func)
        {
            int grid = pti.index();
            auto& ptile = ParticlesAt(a_lev, pti);
            auto& aos = ptile.GetArrayOfStructs();
            auto& soa = ptile.GetStructOfArrays();
            if (aos.numParticles() == 0) { return; }
            auto* p_pbox = aos().data();
            SDProcess::ParticlePointers ptrs;
            setupParticlePointers(soa, aos, ptrs);
            func(pti, grid, p_pbox, ptrs, ctx);
        }

        /*! \brief Iterate over particle tiles with common setup (parallel) */
        template<typename TileFunc>
        AMREX_FORCE_INLINE
        void forEachParticleTile(
            int a_lev,
            const SDProcess::ProcessContext& ctx,
            TileFunc&& func)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
                forEachParticleTileBody(pti, a_lev, ctx, std::forward<TileFunc>(func));
            }
        }

        /*! \brief Iterate over particle tiles WITHOUT OpenMP (serial)
         *  Use when per-tile operations are not thread-safe (e.g., DenseBins).
         */
        template<typename TileFunc>
        AMREX_FORCE_INLINE
        void forEachParticleTileSerial(int a_lev, const SDProcess::ProcessContext& ctx, TileFunc&& func)
        {
            for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
                forEachParticleTileBody(pti, a_lev, ctx, std::forward<TileFunc>(func));
            }
        }

        /*! \brief Lightweight iteration (no ProcessContext) */
        template<typename TileFunc>
        void forEachParticleTile(int a_lev, TileFunc&& func)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (ParIterType pti(*this, a_lev); pti.isValid(); ++pti) {
                int grid = pti.index();
                auto& ptile = ParticlesAt(a_lev, pti);
                auto& aos = ptile.GetArrayOfStructs();
                auto& soa = ptile.GetStructOfArrays();
                const int num_particles = aos.numParticles();
                if (num_particles == 0) { continue; }
                auto* p_pbox = aos().data();
                SDProcess::ParticlePointers ptrs;
                setupParticlePointers(soa, aos, ptrs);
                func(pti, grid, p_pbox, ptrs, num_particles);
            }
        }

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
