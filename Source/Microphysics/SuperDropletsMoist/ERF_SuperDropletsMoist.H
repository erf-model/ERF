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
        q_c,         /*!< Liquid water */
        q_i,         /*!< Ice */
        q_r,         /*!< Rain */
        q_s,         /*!< Snow */
        q_g,         /*!< Graupel */
        dqcdt,       /*!< Condensation rate */
        rh,          /*!< relative humidity */
        rain_accum,  /*!< rain accumulation (mm) */
        graup_accum,  /*!< rain accumulation (mm) */
        snow_accum,  /*!< snow accumulation (mm) */
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_qt_idx (const int a_i, /*!< species index */
                     const int a_is /*!< Index of first non-water/ice species */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-a_is)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::q_t;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_qv_idx (const int a_i, /*!< Species index */
                     const int a_is /*!< Index of first non-water/ice species */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-a_is)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::q_v;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_qc_idx (const int a_i, /*!< Species index */
                     const int a_is /*!< Index of first non-water/ice species */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-a_is)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::q_c;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_sr_idx (const int a_i, /*!< Species index */
                     const int a_is /*!< Index of first non-water/ice species */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-a_is)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::sr;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static int s_accum_idx (const int a_i, /*!< Species index */
                     const int a_is /*!< Index of first non-water/ice species */ )
{
    AMREX_ALWAYS_ASSERT(a_i > 0);
    return      MicVar_SD::NumVars
            +   (a_i-a_is)*MicVar_SD_Species::NumVars
            +   MicVar_SD_Species::accum;
}

/*! \brief List of super-droplet moisture model initializations */
namespace SuperDropletsMoistInitializations
{
    /*! Create super-droplets based on given condensate mass density */
    const std::string init_rhoc = "condensate_density";

}

AMREX_ENUM(SDMSimulationDim,
    one_d_z, two_d_xz, two_d_yz, three_d
);

using namespace amrex;

class SuperDropletsMoist : public NullMoistLagrangian {

    using FabPtr = std::shared_ptr<MultiFab>;
    using MFPtr = std::unique_ptr<MultiFab>;
    using BCTypeArr = amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>;

    public:

        /*! \brief Constructor */
        SuperDropletsMoist()
        {
            m_qmoist_size = 11;
            m_qstate_moist_size = 6; // qv, qc, qi, qrain, qsnow, qgraup
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
        virtual void Init ( const MultiFab&,
                            const BoxArray&,
                            const Geometry&,
                            const Real&,
                            MFPtr&,
                            MFPtr& ) override;

        /*! \brief Initialize particles */
        virtual void InitParticles ( MFPtr& ) override;

        /*! \brief Restart particles */
        virtual void RestartParticles ( ParGDBBase*, const std::string& ) override;

        /*! \brief finish initializations */
        virtual void FinishInit(const int&,
                                MultiFab&,
                                const Vector<MFPtr>&) override;

        void Advance ( const amrex::Real&, const SolverChoice& ) override
        {
            amrex::Abort("Do not use this advance() for Lagrangian microphysics.");
        }

        /*! \brief Advance for one timestep */
        virtual void Advance (const Real&,
                              const int&,
                              const Real&,
                              Vector<Vector<MultiFab>>&,
                              const Vector<MFPtr>&,
                              const BCTypeArr& ) override;

        /*! \brief Update microphysics variables */
        virtual void Update_Micro_Vars (MultiFab& a_cons_vars) override;

        /*! \brief Update state variables from microphysics variables */
        virtual void Update_State_Vars (MultiFab& a_cons_vars) override;

        /*! \brief Copy moisture model vars from state to member multifabs */
        virtual void Copy_State_to_Micro (const MultiFab&) override;

        /*! \brief Copy moisture model vars to state from member multifabs */
        virtual void Copy_Micro_to_State (MultiFab&) override;

        /*! \brief returns the moisture variable asked for */
        inline virtual MultiFab* Qmoist_Ptr (const int& a_idx) override
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

        /*! \brief compute cloud/rain mixing ratios */
        virtual void computeQcQrWater ();

        /*! \brief compute ice/graupel/snow mixing ratios */
        virtual void computeQiQgQsWater ();

        /*! \brief compute condensate for all non-water species */
        virtual void computeQcSpecies ()
        {
            for (int is = m_istart_sp; is < m_num_species; is++) { computeQcSpecies(is); }
        }

        /*! \brief compute condensate mixing ratio for a non-water species */
        virtual void computeQcSpecies (const  int a_i);

        /*! \brief compute condensate mixing ratio for a species */
        virtual void computeQc (const  int a_i)
        {
            if (a_i == m_idx_w) {
                computeQcQrWater();
                computeQiQgQsWater();
            } else {
                computeQcSpecies(a_i);
            }
        }

        /*! \brief compute total water mixing ratio */
        virtual void computeQtWater ();

        /*! \brief compute qt (total) for all non-water species */
        virtual void computeQtSpecies ()
        {
            for (int is = m_istart_sp; is < m_num_species; is++) { computeQtSpecies(is); }
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
        virtual void rainAccumulation ();

        /*! \brief Compute rain accumulation */
        virtual void snowAccumulation ();

        /*! \brief Compute rain accumulation */
        virtual void speciesAccumulation ();

        /*! \brief Convert a multifab containing density of something to its mixing ratio */
        void densityToRatio ( MultiFab&, const int a_comp = 0 );
        /*! \brief Convert a multifab containing the mixing ratio of something to its density */
        void ratioToDensity ( MultiFab&, const int a_comp = 0 );

        /*! \brief Compute phase changes (vapour <--> cloud water) */
        virtual void phaseChange (  const Real&,
                                    const Vector<MFPtr>& );

        virtual void GetPlotVarNames (Vector<std::string>& a_names) const override
        {
            for (int v = m_istart_sp; v < m_num_species; v++) {
                a_names.push_back("qv_"+getEnumNameString(m_species[v]));
                a_names.push_back("qc_"+getEnumNameString(m_species[v]));
                a_names.push_back("qt_"+getEnumNameString(m_species[v]));
                a_names.push_back("sat_ratio_"+getEnumNameString(m_species[v]));
                a_names.push_back("accum_"+getEnumNameString(m_species[v]));
            }
        }

        virtual void GetPlotVar (const std::string& a_name,
                                 MultiFab& a_mf )  const override
        {
            AMREX_ASSERT(a_mf.nComp() >= 1);
            for (int v = m_istart_sp; v < m_num_species; v++) {
                {
                    std::string varname = "qv_" + getEnumNameString(m_species[v]);
                    if (a_name == varname) {
                        MultiFab::Copy( a_mf,
                                        *m_mic_fab_vars[s_qv_idx(v,m_istart_sp)],
                                        0, 0, 1, IntVect::TheZeroVector() );
                        return;
                    }
                }
                {
                    std::string varname = "qc_" + getEnumNameString(m_species[v]);
                    if (a_name == varname) {
                        MultiFab::Copy( a_mf,
                                        *m_mic_fab_vars[s_qc_idx(v,m_istart_sp)],
                                        0, 0, 1, IntVect::TheZeroVector() );
                        return;
                    }
                }
                {
                    std::string varname = "qt_" + getEnumNameString(m_species[v]);
                    if (a_name == varname) {
                        MultiFab::Copy( a_mf,
                                        *m_mic_fab_vars[s_qt_idx(v,m_istart_sp)],
                                        0, 0, 1, IntVect::TheZeroVector() );
                        return;
                    }
                }
                {
                    std::string varname = "sat_ratio_" + getEnumNameString(m_species[v]);
                    if (a_name == varname) {
                        MultiFab::Copy( a_mf,
                                        *m_mic_fab_vars[s_sr_idx(v,m_istart_sp)],
                                        0, 0, 1, IntVect::TheZeroVector() );
                        return;
                    }
                }
                {
                    std::string varname = "accum_" + getEnumNameString(m_species[v]);
                    if (a_name == varname) {
                        MultiFab::Copy( a_mf,
                                        *m_mic_fab_vars[s_accum_idx(v,m_istart_sp)],
                                        0, 0, 1, IntVect::TheZeroVector() );
                        return;
                    }
                }
            }
            amrex::Abort("SuperDropletsMoist::GetPlotVar() called with invalid name");
            return;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int q_qv_idx (const int a_i, /*!< species index */
                      const int a_is /*!< index of first non-water/ice species */ )
        {
            AMREX_ALWAYS_ASSERT(a_i > 0);
            return RhoQ1_comp + m_qstate_moist_size + 2*(a_i-a_is) + 0;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int q_qc_idx (const int a_i, /*!< species index */
                      const int a_is /*!< index of first non-water/ice species */ )
        {
            AMREX_ALWAYS_ASSERT(a_i > 0);
            return RhoQ1_comp + m_qstate_moist_size + 2*(a_i-a_is) + 1;
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
        std::string m_init_type;

        /*! name of model and its particle container */
        std::string m_name;
        /*! particle container for super-droplets */
        SuperDropletPC* m_super_droplets;

        /*! Geometry object */
        amrex::Geometry m_geom;

        /*! number of microphysics variables */
        int m_qmoist_size;
        /*! number of water-related state variables for this moisture model */
        int m_qstate_moist_size;
        /*! number of non-water-related state variables for this moisture model */
        int m_qstate_nonmoist_size;

        Real m_dt; /*!< timestep */
        Vector<int> m_mic_var_map; /*!< moisture model variables map */

        /*! moisture model variables */
        Vector<FabPtr> m_mic_fab_vars;

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

        /*! species index of ice */
        int m_idx_i;

        /*! starting index of non-moisture species */
        int m_istart_sp;

        /*! number of non-moisture species */
        int m_num_nonmoist_sp;

        /*! include cold processes? */
        bool m_with_ice;

        /*! \brief read inputs */
        virtual void readInputs();

    private:

};

#endif
#endif
