#ifndef SUPERDROPLETSMOIST_H
#define SUPERDROPLETSMOIST_H

#ifdef ERF_USE_PARTICLES

#include <limits>
#include <string>
#include <AMReX_Geometry.H>
#include "NullMoistLagrangian.H"
#include "SuperDropletPC.H"

namespace MicVar_SD {
    enum {
        rho = 0,     /*!< Density */
        temperature, /*!< Temperature */
        pressure,    /*!< Pressure */
        q_t,         /*!< Total (vapour + cloud) */
        q_v,         /*!< Water vapour */
        q_c,         /*!< Cloud condensate */
        NumVars
    };
}

/*! \brief List of super-droplet moisture model initializations */
namespace SuperDropletsMoistInitializations
{
    /*! Create super-droplets based on given condensate mass density */
    const std::string init_rhoc = "condensate_density";

}

using namespace amrex;

class SuperDropletsMoist : public NullMoistLagrangian {

    using FabPtr = std::shared_ptr<MultiFab>;
    using MFPtr = std::unique_ptr<MultiFab>;

    public:

        /*! \brief Constructor */
        SuperDropletsMoist() { }

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

        void Advance ( const amrex::Real&, const SolverChoice& ) override
        {
            amrex::Abort("Do not use this advance() for Lagrangian microphysics.");
        }

        /*! \brief Advance for one timestep */
        virtual void Advance (const Real&,
                              const int&,
                              const Real&,
                              Vector<Vector<MultiFab>>&,
                              const Vector<MFPtr>& ) override;

        /*! \brief Update microphysics variables */
        inline virtual void Update_Micro_Vars (MultiFab& a_cons_vars) override
        {
            Copy_State_to_Micro(a_cons_vars);

            static bool first_call(true);
            if (first_call) { FinishInit(a_cons_vars); }
            first_call = false;
        }

        /*! \brief state variables from microphysics variables */
        inline virtual void Update_State_Vars (MultiFab& a_cons_vars) override
        {
            computeQc();
            computeQt();
            Copy_Micro_to_State(a_cons_vars);
        }

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
        inline virtual int Qmoist_Size () override { return m_qmoist_size; }

        /*! \brief returns number of conserved variables in this moisture model */
        inline virtual int Qstate_Size () override { return m_qstate_size; }

        /*! \brief returns pointer to super-droplets particle container */
        inline virtual ERFPC* getParticleContainer() override
        {
            return m_super_droplets;
        }

        inline virtual const std::string& getName() const override
        {
            return m_name;
        }

    //protected: ideally, but can't due to extended lambda capture rules (CUDA)

        /*! \brief compute condensate mixing ratio */
        inline virtual void computeQc ()
        {
            m_super_droplets->massDensityCondensate(*(m_mic_fab_vars[MicVar_SD::q_c]));
            densityToRatio(*(m_mic_fab_vars[MicVar_SD::q_c]));
        }

        /*! \brief compute qt (total) */
        inline virtual void computeQt ()
        {
            for ( MFIter mfi(*m_mic_fab_vars[MicVar_SD::q_t]); mfi.isValid(); ++mfi) {

                Box bx = mfi.tilebox();
                bx.grow(m_mic_fab_vars[MicVar_SD::q_t]->nGrowVect());

                auto q_c_arr = m_mic_fab_vars[MicVar_SD::q_c]->const_array(mfi);
                auto q_v_arr = m_mic_fab_vars[MicVar_SD::q_v]->const_array(mfi);
                auto q_t_arr = m_mic_fab_vars[MicVar_SD::q_t]->array(mfi);

                ParallelFor( bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                { q_t_arr(i,j,k) = q_v_arr(i,j,k) + q_c_arr(i,j,k); });
            }
        }

        /*! \brief finish initializations */
        void FinishInit(MultiFab&);

        /*! \brief Convert a multifab containing density of something to its mixing ratio */
        void densityToRatio ( MultiFab&, const int a_comp = 0 );
        /*! \brief Convert a multifab containing the mixing ratio of something to its density */
        void ratioToDensity ( MultiFab&, const int a_comp = 0 );

        /*! \brief Compute phase changes (vapour <--> cloud water) */
        virtual void phaseChange (  const Real&,
                                    const Vector<MFPtr>& );

    protected:

        bool m_flag_phase_change;   /*!< Enable/disable phase change */
        bool m_flag_advection;      /*!< Enable/disable advection */
        bool m_flag_coalescence;    /*!< Enable/disable coalescence */

        int m_diagnostics_iter; /*!< number of iterations between computing diagnostics */

        /*! initialization type */
        std::string m_init_type;

        /*! name of model and its particle container */
        const std::string m_name = "super_droplets_moisture";
        /*! particle container for super-droplets */
        SuperDropletPC* m_super_droplets;

        /*! Geometry object */
        amrex::Geometry m_geom;

        /*! vapour and condensed water density */
        const int m_qmoist_size = 3;
        /*! vapour and condensed water density */
        const int m_qstate_size = 2;

        Real m_dt; /*!< timestep */
        Vector<int> m_mic_var_map; /*!< moisture model variables map */

        /*! moisture model variables */
        Array<FabPtr, MicVar_SD::NumVars> m_mic_fab_vars;

        /*! names of aerosols */
        std::vector<std::string> m_aerosols;

        /*! number of substeps for phase change */
        int m_num_substeps_phase_change;

        /*! \brief read inputs */
        virtual void readInputs();

    private:

};

#endif
#endif
