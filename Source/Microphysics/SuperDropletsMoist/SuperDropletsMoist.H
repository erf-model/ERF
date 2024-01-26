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

        /*! \brief Advance for one timestep */
        virtual void Advance (const Real&,
                              Vector<Vector<MultiFab>>&,
                              const Vector<std::unique_ptr<MultiFab>>& );

        /*! \brief TODO */
        virtual void Update_Micro_Vars (MultiFab&) { }

        /*! \brief TODO */
        virtual void Update_State_Vars (MultiFab&) { }

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

    private:

};

#endif
#endif
