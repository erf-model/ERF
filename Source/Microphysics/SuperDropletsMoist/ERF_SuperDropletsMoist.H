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
            m_qmoist_size = 7;
            m_qstate_moist_size = 3; // qv, qc, qr
            m_name = "super_droplets_moisture";
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

        /*! \brief returns number of conserved variables in this moisture model */
        inline virtual int Qstate_Moist_Size () override
        {
            return m_qstate_moist_size;
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

        /*! \brief compute condensate mixing ratio */
        virtual void computeQcQr ();

        /*! \brief compute qt (total) */
        virtual void computeQt ();

        /*! \brief Compute rain accumulation */
        virtual void rainAccumulation ();

        /*! \brief Convert a multifab containing density of something to its mixing ratio */
        void densityToRatio ( MultiFab&, const int a_comp = 0 );
        /*! \brief Convert a multifab containing the mixing ratio of something to its density */
        void ratioToDensity ( MultiFab&, const int a_comp = 0 );

        /*! \brief Compute phase changes (vapour <--> cloud water) */
        virtual void phaseChange (  const Real&,
                                    const Vector<MFPtr>&,
                                    const bool );

    protected:

        bool m_flag_phase_change;   /*!< Enable/disable phase change */
        bool m_flag_advection;      /*!< Enable/disable advection */
        bool m_flag_coalescence;    /*!< Enable/disable coalescence */

        amrex::Real m_fac_cond; /*!< latent heat factor */

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
        /*! number of state variables for this moisture model */
        int m_qstate_moist_size;

        Real m_dt; /*!< timestep */
        Vector<int> m_mic_var_map; /*!< moisture model variables map */

        /*! moisture model variables */
        Array<FabPtr, MicVar_SD::NumVars> m_mic_fab_vars;

        /*! names of aerosols */
        std::vector<Species::Name> m_aerosols;

        /*! number of substeps for phase change */
        int m_num_substeps_phase_change;

        /*! kinematic mode? */
        bool m_kinematic_mode;

        /*! update qv from dycore */
        bool m_update_qv;

        /*! Dimensionality of simulation*/
        SDMSimulationDim m_dimensionality;

        /*! recycle particles */
        bool m_recycle_particles;

        /*! \brief read inputs */
        virtual void readInputs();

    private:

};

#endif
#endif
