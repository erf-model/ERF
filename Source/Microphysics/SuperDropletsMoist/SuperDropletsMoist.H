#ifndef SUPERDROPLETSMOIST_H
#define SUPERDROPLETSMOIST_H

#ifdef ERF_USE_PARTICLES

#include <string>

#include "NullMoist.H"
#include "SuperDropletPC.H"

namespace MicVar_SD {
    enum {
        rho_v, /*!< Water vapour */
        rho_c, /*!< Cloud condensate */
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

    public:

        /*! \brief Constructor */
        SuperDropletsMoist() { }

        /*! \brief Destructor */
        virtual ~SuperDropletsMoist() = default;

        /*! \brief Define super-droplet moisture model parameters */
        virtual void Define (SolverChoice&);

        /*! \brief Initialize super-droplet moisture model data structures */
        virtual void Init ( const MultiFab&,
                            const BoxArray&,
                            const Geometry&,
                            const Real& );

        /* This is here just to avoid CUDA compiler warnings about virtual
           function override. Don't ever call this! */
        using NullMoist::Advance;

        /*! \brief Advance for one timestep */
        virtual void Advance (const Real&,
                              Vector<Vector<MultiFab>>&,
                              const Vector<std::unique_ptr<MultiFab>>& );

        /*! \brief Update microphysics variables */
        inline virtual void Update_Micro_Vars (MultiFab& a_cons_vars)
        {
            Copy_State_to_Micro(a_cons_vars);

            static bool first_call(true);
            if (first_call) { FinishInit(); }
            first_call = false;
        }

        /*! \brief state variables from microphysics variables */
        inline virtual void Update_State_Vars (MultiFab& a_cons_vars)
        {
            Copy_Micro_to_State(a_cons_vars);
        }

        /*! \brief TODO */
        virtual void Diagnose () { }

        /*! \brief Copy moisture model vars from state to member multifabs */
        virtual void Copy_State_to_Micro (const MultiFab&);

        /*! \brief Copy moisture model vars to state from member multifabs */
        virtual void Copy_Micro_to_State (MultiFab&);

        /*! \brief returns the moisture variable asked for */
        inline virtual MultiFab* Qmoist_Ptr (const int& a_idx)
        {
            AMREX_ALWAYS_ASSERT( a_idx < m_qmoist_size );
            return m_mic_fab_vars[m_mic_var_map[a_idx]].get();
        }

        /*! \brief returns number of variables in this moisture model */
        inline virtual int Qmoist_Size () { return m_qmoist_size; }

        /*! \brief returns number of conserved variables in this moisture model */
        inline virtual int Qstate_Size () { return m_qstate_size; }

        /*! \brief returns pointer to super-droplets particle container */
        inline virtual ERFPC* getParticleContainer()
        {
            return m_super_droplets;
        }

        inline virtual const std::string& getName() const
        {
            return m_name;
        }

    protected:

        /*! initialization type */
        std::string m_init_type;

        /*! name of model and its particle container */
        const std::string m_name = "super_droplets_moisture";
        /*! particle container for super-droplets */
        SuperDropletPC* m_super_droplets;

        /*! vapour and condensed water density */
        const int m_qmoist_size = 2;
        /*! vapour and condensed water density */
        const int m_qstate_size = 2;

        Real m_dt; /*!< timestep */
        Vector<int> m_mic_var_map; /*!< moisture model variables map */

        /*! moisture model variables */
        Array<FabPtr, MicVar_SD::NumVars> m_mic_fab_vars;

        /*! \brief read inputs */
        virtual void readInputs();

        /*! \brief finish initializations */
        void FinishInit();

    private:

};

#endif
#endif
