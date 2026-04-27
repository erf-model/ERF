#ifndef SUPERDROPLETSMOIST_H
#define SUPERDROPLETSMOIST_H

#ifdef ERF_USE_PARTICLES

#include <limits>
#include <string>
#include <AMReX_Enum.H>
#include <AMReX_Geometry.H>
#include "ERF_NullMoistLagrangian.H"
#include "ERF_SuperDropletPC.H"

namespace MicVar_SD {
    enum {
        rho = 0,     /*!< Density */
        theta,       /*!< Potential temperature */
        temperature, /*!< Temperature */
        pressure,    /*!< Pressure */
        q_t,         /*!< Total (vapour + cloud) */
        q_v,         /*!< Water vapour */
        q_c,         /*!< Cloud condensate */
        dqcdt,       /*!< Condensation rate */
        q_r,         /*!< Rain */
        rh,          /*!< relative humidity */
        rain_accum,  /*!< rain accumulation (mm) */
        NumVars
    };
}

namespace MicVar_SD_Species {
    enum {
        q_t = 0,    /*!< Total density ratio */
        q_v,        /*!< Vapour density ratio */
        q_c,        /*!< Condensate density ratio */
        sr,         /*!< Saturation ratio */
        accum,      /*!< Ground accumulation */
        NumVars
    };
}

namespace MicVar_SD_Aerosols {
    enum {
        accum = 0,  /*!< Ground accumulation */
        NumVars
    };
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_qt_idx (const int a_i /*!< species index */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-1)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::q_t;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_qv_idx (const int a_i /*!< species index */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-1)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::q_v;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_qc_idx (const int a_i /*!< species index */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-1)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::q_c;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_sr_idx (const int a_i /*!< species index */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-1)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::sr;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_accum_idx (const int a_i /*!< species index */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-1)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::accum;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int a_accum_idx (const int a_ns, /*!< number of species */
                        const int a_i /*!< aerosol index */ )
{
    AMREX_ALWAYS_ASSERT(a_i >= 0);
    return      MicVar_SD::NumVars
            +   (a_ns-1)*MicVar_SD_Species::NumVars
            +   a_i*MicVar_SD_Aerosols::NumVars
            +   MicVar_SD_Aerosols::accum;
}

/*! \brief List of super-droplet moisture model initializations */
AMREX_ENUM(SDMoistInit,
    uniform,
    condensate_density
);

AMREX_ENUM(SDMSimulationDim,
    one_d_z, two_d_xz, two_d_yz, three_d
);

class SuperDropletsMoist : public NullMoistLagrangian {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;
    using MFPtr = std::unique_ptr<amrex::MultiFab>;
    using BCTypeArr = amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>;

    public:

        /*! \brief Constructor */
        SuperDropletsMoist()
        {
            m_qmoist_size = 7;
            m_qstate_moist_size = 3; // qv, qc, qr
            m_qstate_nonmoist_size = -1; // intentional
            m_name = "super_droplets_moisture";
            m_idx_w = -1;
            m_num_species = 0;
            m_num_aerosols = 0;
            readInputs();
        }

        /*! \brief Destructor */
        virtual ~SuperDropletsMoist() = default;

        /*! \brief Define super-droplet moisture model parameters */
        virtual void Define (SolverChoice&) override;

        /*! \brief Initialize super-droplet moisture model data structures */
        virtual void Init ( const amrex::MultiFab&,
                            const amrex::BoxArray&,
                            const amrex::Geometry&,
                            const amrex::Real&,
                            MFPtr&,
                            MFPtr& ) override;

        /*! \brief Initialize particles */
        virtual void InitParticles ( MFPtr& ) override;

        /*! \brief Restart particles */
        virtual void RestartParticles ( amrex::ParGDBBase*, const std::string& ) override;

        /*! \brief finish initializations */
        virtual void FinishInit(const int&,
                                amrex::MultiFab&,
                                const amrex::Vector<MFPtr>&) override;

        void Advance ( const amrex::Real&, const SolverChoice& ) override
        {
            amrex::Abort("Do not use this advance() for Lagrangian microphysics.");
        }

        /*! \brief Advance for one timestep
         * \param[in] dt Timestep
         * \param[in] iteration Current iteration number
         * \param[in] time Current simulation time
         * \param[in,out] cons_vars Array of conserved variables
         * \param[in] mf_array Array of additional multifabs
         * \param[in] bc_arr Array of boundary conditions
         */
        virtual void Advance (const amrex::Real& dt,
                              const int& iteration,
                              const amrex::Real& time,
                              amrex::Vector<amrex::Vector<amrex::MultiFab>>& cons_vars,
                              const amrex::Vector<MFPtr>& mf_array,
                              const BCTypeArr& bc_arr) override;

        /*! \brief Update microphysics variables
         * \param[in] a_cons_vars Conservative variables multifab
         */
        virtual void Update_Micro_Vars (amrex::MultiFab& a_cons_vars) override;

        /*! \brief Update state variables from microphysics variables
         * \param[in,out] a_cons_vars Conservative variables multifab to update
         * \param[in] a_z_phys_nd Terrain heights (nodal)
         */
        virtual void Update_State_Vars (amrex::MultiFab& a_cons_vars,
                                        const amrex::MultiFab& a_z_phys_nd) override;

        /*! \brief Copy moisture model vars from state to member multifabs
         * \param[in] a_state State multifab to copy from
         */
        virtual void Copy_State_to_Micro (const amrex::MultiFab& a_state) override;

        /*! \brief Copy moisture model vars to state from member multifabs
         * \param[out] a_state State multifab to copy to
         */
        virtual void Copy_Micro_to_State (amrex::MultiFab& a_state) override;

        /*! \brief returns the moisture variable asked for */
        inline virtual amrex::MultiFab* Qmoist_Ptr (const int& a_idx) override
        {
            AMREX_ALWAYS_ASSERT( a_idx < m_qmoist_size );
            return m_mic_fab_vars[m_mic_var_map[a_idx]].get();
        }

        /*! \brief returns number of variables in this moisture model */
        inline virtual int Qmoist_Size () override
        {
            return m_qmoist_size;
        }

        /*! \brief returns number of moist conserved variables */
        inline virtual int Qstate_Moist_Size () override
        {
            return m_qstate_moist_size;
        }

        /*! \brief returns number of non-water species conserved variables */
        inline virtual int Qstate_NonMoist_Size () override
        {
            AMREX_ALWAYS_ASSERT(m_qstate_nonmoist_size >= 0);
            return m_qstate_nonmoist_size;
        }

        /*! \brief returns pointer to super-droplets particle container */
        inline virtual ERFPC* getParticleContainer() override
        {
            return m_super_droplets;
        }

        inline virtual const std::string& getName() const override
        {
            return m_name;
        }

        /*! \brief compute cloud/rain mixing ratio */
        virtual void computeQcQrWater (const amrex::MultiFab& a_z_phys_nd);

        /*! \brief compute condensate for all non-water species */
        virtual void computeQcSpecies (const amrex::MultiFab& a_z_phys_nd)
        {
            for (int is = 1; is < m_num_species; is++) { computeQcSpecies(is, a_z_phys_nd); }
        }

        /*! \brief compute condensate mixing ratio for a non-water species */
        virtual void computeQcSpecies (const int a_i, const amrex::MultiFab& a_z_phys_nd);

        /*! \brief compute condensate mixing ratio for a species */
        virtual void computeQc (const int a_i, const amrex::MultiFab& a_z_phys_nd)
        {
            if (a_i == m_idx_w) { computeQcQrWater(a_z_phys_nd); }
            else                { computeQcSpecies(a_i, a_z_phys_nd); }
        }

        /*! \brief compute total water mixing ratio */
        virtual void computeQtWater ();

        /*! \brief compute qt (total) for all non-water species */
        virtual void computeQtSpecies ()
        {
            for (int is = 1; is < m_num_species; is++) { computeQtSpecies(is); }
        }

        /*! \brief compute qt (total) for a non-water species */
        virtual void computeQtSpecies (const  int a_i);

        /*! \brief compute qt (total) for a species */
        virtual void computeQt (const  int a_i)
        {
            if (a_i == m_idx_w) { computeQtWater(); }
            else                { computeQtSpecies(a_i); }
        }

        /*! \brief Compute rain accumulation */
        virtual void rainAccumulation (const amrex::MultiFab& a_z_phys_nd);

        /*! \brief Compute non-water species accumulation */
        virtual void speciesAccumulation (const amrex::MultiFab& a_z_phys_nd);

        /*! \brief Compute aerosol accumulation */
        virtual void aerosolAccumulation (const amrex::MultiFab& a_z_phys_nd);

        /*! \brief Convert a multifab containing density of something to its mixing ratio */
        void densityToRatio ( amrex::MultiFab&, const int a_comp = 0 );
        /*! \brief Convert a multifab containing the mixing ratio of something to its density */
        void ratioToDensity ( amrex::MultiFab&, const int a_comp = 0 );

        /*! \brief Compute phase changes (vapour <--> cloud water)
         * \param[in] a_dt Timestep for phase change calculation
         * \param[in] a_z Array of terrain heights
         * \param[in] a_update_qv Flag to enable updating vapour mixing ratio
         */
        virtual void phaseChange (  const amrex::Real& a_dt,
                                    const amrex::Vector<MFPtr>& a_z,
                                    const bool a_update_qv);

        virtual void GetPlotVarNames (amrex::Vector<std::string>& a_names) const override
        {
            for (int v = 1; v < m_num_species; v++) {
                a_names.push_back("qv_"+amrex::getEnumNameString(m_species[v]));
                a_names.push_back("qc_"+amrex::getEnumNameString(m_species[v]));
                a_names.push_back("qt_"+amrex::getEnumNameString(m_species[v]));
                a_names.push_back("sat_ratio_"+amrex::getEnumNameString(m_species[v]));
                a_names.push_back("accum_"+amrex::getEnumNameString(m_species[v]));
            }
            for (int v = 0; v < m_num_aerosols; v++) {
                a_names.push_back("accum_"+amrex::getEnumNameString(m_aerosols[v]));
            }
        }

        virtual void GetPlotVar (const std::string& a_name,
                                 amrex::MultiFab& a_mf )  const override
        {
            AMREX_ASSERT(a_mf.nComp() >= 1);

            // Helper to copy MultiFab if name matches
            auto try_copy = [&](const std::string& prefix, int idx) -> bool {
                if (a_name == prefix) {
                    amrex::MultiFab::Copy(a_mf, *m_mic_fab_vars[idx],
                                          0, 0, 1, amrex::IntVect::TheZeroVector());
                    return true;
                }
                return false;
            };

            // Species variables
            for (int v = 1; v < m_num_species; v++) {
                std::string sp_name = amrex::getEnumNameString(m_species[v]);
                if (try_copy("qv_" + sp_name, s_qv_idx(v)) ||
                    try_copy("qc_" + sp_name, s_qc_idx(v)) ||
                    try_copy("qt_" + sp_name, s_qt_idx(v)) ||
                    try_copy("sat_ratio_" + sp_name, s_sr_idx(v)) ||
                    try_copy("accum_" + sp_name, s_accum_idx(v))) {
                    return;
                }
            }

            // Aerosol variables
            for (int v = 0; v < m_num_aerosols; v++) {
                std::string ae_name = amrex::getEnumNameString(m_aerosols[v]);
                if (try_copy("accum_" + ae_name, a_accum_idx(m_num_species, v))) { return; }
            }

            amrex::Abort("SuperDropletsMoist::GetPlotVar() called with invalid name");
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int q_qv_idx (const int a_i /*!< species index */ )
        {
            AMREX_ALWAYS_ASSERT(a_i > 0);
            return RhoQ1_comp + m_qstate_moist_size + 2*(a_i-1) + 0;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int q_qc_idx (const int a_i /*!< species index */ )
        {
            AMREX_ALWAYS_ASSERT(a_i > 0);
            return RhoQ1_comp + m_qstate_moist_size + 2*(a_i-1) + 1;
        }

    protected:

        bool m_flag_phase_change;   /*!< Enable/disable phase change */
        bool m_flag_advection;      /*!< Enable/disable advection */
        bool m_flag_coalescence;    /*!< Enable/disable coalescence */

        amrex::Real m_Cp; /*!< specific heat at constant pressure for dry air */

        amrex::Real m_r_rain; /*!< minimum radius to be considered rain */

        /*! let initial superdroplets "relax" to a physically-appropriate size? */
        bool m_init_phase_change;
        /*! time (in seconds) of initial relaxation */
        amrex::Real m_init_phase_change_time;

        int m_diagnostics_iter; /*!< number of iterations between computing diagnostics */

        /*! initialization type */
        SDMoistInit m_init_type;

        /*! name of model and its particle container */
        std::string m_name;
        /*! Geometry object */
        amrex::Geometry m_geom;

        /*! number of microphysics variables */
        int m_qmoist_size;
        /*! number of water-related state variables for this moisture model */
        int m_qstate_moist_size;
        /*! number of non-water-related state variables for this moisture model */
        int m_qstate_nonmoist_size;

        amrex::Real m_dt; /*!< timestep */
        amrex::Vector<int> m_mic_var_map; /*!< moisture model variables map */

        /*! moisture model variables */
        amrex::Vector<FabPtr> m_mic_fab_vars;

        /*! names of vapour/condensate species */
        std::vector<Species::Name> m_species;
        int m_num_species;

        /*! names of aerosols */
        std::vector<Species::Name> m_aerosols;
        int m_num_aerosols;

        /*! number of substeps for phase change */
        int m_num_substeps_phase_change;

        /*! kinematic mode? */
        bool m_kinematic_mode;

        /*! Dimensionality of simulation*/
        SDMSimulationDim m_dimensionality;

        /*! recycle particles */
        bool m_recycle_particles;

        /*! species index of water */
        int m_idx_w;

        /*! particle container for super-droplets
         *  Owned by ERF::particleData which deletes it during teardown;
         *  raw pointer here to avoid double-free. */
        SuperDropletPC* m_super_droplets;

        /*! \brief read inputs */
        virtual void readInputs();

    private:

};

#endif
#endif
