/*! @file ERF_LagrangianMicrophysics.H
 *  \brief Contains the Lagrangian microphysics class
 */

#ifndef ERF_LAGRANGIANMICROPHYSICS_H
#define ERF_LAGRANGIANMICROPHYSICS_H

#ifdef ERF_USE_PARTICLES

#include <utility>
#include <string>

#include "ERF_NullMoistLagrangian.H"
#include "ERF_SuperDropletsMoist.H"
#include "ERF_Microphysics.H"

/* forward declaration */
class ERFPC;

/*! \brief Eulerian microphysics interface
 *
 * One key difference from #EulerianMicrophysics is that only the base AMR
 * level has the moisture model. Thus, for higher AMR levels, the microphysics
 * interface need not do anything. The number of conserved state variables for
 * moisture (RhoQ1, RhoQ2, ...) are the same for all AMR levels; however, for
 * levels other than the base, they just contain and evolve zeros. */
class LagrangianMicrophysics : public Microphysics {

public:

    /*! \brief Null constructor */
    LagrangianMicrophysics () { }

    /*! \brief default destructor */
    ~LagrangianMicrophysics () = default;

    /*! \brief Constructor: create the moisture model */
    LagrangianMicrophysics (const int& /* nlev */, /*!< Number of AMR levels */
                            const MoistureType& a_model_type /*!< moisture model */ )
    {
        AMREX_ASSERT( Microphysics::modelType(a_model_type) == MoistureModelType::Lagrangian );
        if (a_model_type == MoistureType::SuperDroplets) {
            SetModel<SuperDropletsMoist>();
            amrex::Print() << "Super-Droplets moisture model!\n";
        } else {
            amrex::Abort("Unknown Lagrangian moisture model!") ;
        }
    }

    /*! \brief Define the moisture model */
    void Define (const int& lev, /*!< AMR level */
                 SolverChoice& sc /*!< Solver choice object */) override
    {
        if (lev > 0) return;
        m_moist_model->Define(sc);
    }

    /*! \brief Initialize the moisture model */
    void Init (const int& lev, /*!< AMR level */
               const amrex::MultiFab& cons_in, /*!< Conserved state variables */
               const amrex::BoxArray& grids,   /*!< Grids */
               const amrex::Geometry& geom,    /*!< Geometry */
               const amrex::Real& dt_advance,  /*!< Time step */
               std::unique_ptr<amrex::MultiFab>& z_phys_nd,/*< Nodal z heights */
               std::unique_ptr<amrex::MultiFab>& detJ_cc /*< CC Jacobian determinants */) override
    {
        if (lev == 0) {
            m_moist_model->Init(cons_in, grids, geom, dt_advance,
                                z_phys_nd, detJ_cc);
        } else {
            m_moist_model->InitLevel(lev, cons_in);
        }
    }

    /*! \brief finish initializations steps that require flow variables */
    void FinishInit (const int& lev, /*!< AMR level */
                     amrex::MultiFab& cons_in, /*!< Conserved state variable */
                     const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd /*!< terrain */) override
    {
        if (lev > 0) return;
        m_moist_model->FinishInit( lev, cons_in, z_phys_nd );
    }

    /*! \brief Advance the moisture model for one time step */
    void Advance (const int& lev, /*!< AMR level */
                  const amrex::Real& dt_advance, /*!< Time step */
                  const int& iter, /*!< iteration number */
                  const amrex::Real& time, /*!< current time */
                  const SolverChoice& /*solverChoice*/, /*!< Solver choice object */
                  amrex::Vector<amrex::Vector<amrex::MultiFab>>& a_vars, /*!< Dycore state variables */
                  const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z/*!< terrain */,
                  const amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>& a_phys_bc_types ) override
    {
        m_moist_model->SetCurrentLevel(lev);
        m_moist_model->Advance(dt_advance, iter, time, a_vars, a_z, a_phys_bc_types);
    }

    /*! \brief update microphysics variables from ERF state variables */
    void Update_Micro_Vars_Lev (const int& lev, /*! AMR level */
                                amrex::MultiFab& cons_in /*!< Conserved state variables */) override
    {
        m_moist_model->SetCurrentLevel(lev);
        m_moist_model->Update_Micro_Vars(cons_in);
    }

    /*! \brief update ERF state variables from microphysics variables */
    void Update_State_Vars_Lev (const int& lev, /*!< AMR level */
                                amrex::MultiFab& cons_in, /*!< Conserved state variables */
                                const amrex::MultiFab& z_phys_nd /*!< Nodal terrain heights */) override
    {
        m_moist_model->SetCurrentLevel(lev);
        m_moist_model->Update_State_Vars(cons_in, z_phys_nd);
    }

    /*! \brief average down moisture multifabs from finest_level to 0 */
    void AverageDownMicroVars (const int finest_level) override
    {
        m_moist_model->AverageDownMicroVars(finest_level);
    }

    /*! \brief get pointer to a moisture variable */
    amrex::MultiFab* Get_Qmoist_Ptr (const int& lev, /*!< AMR level */
                                     const int& varIdx /*!< moisture variable index */) override
    {
        m_moist_model->SetCurrentLevel(lev);
        return m_moist_model->Qmoist_Ptr(varIdx);
    }

    /*! \brief get the number of moisture model variables */
    int Get_Qmoist_Size (const int& /*a_lev*/ /*!< AMR level */) override
    {
        return m_moist_model->Qmoist_Size();
    }

    /*! \brief get the number of  microphysics conserved moist (water-related) state variables */
    int Get_Qstate_Moist_Size () override
    {
        return m_moist_model->Qstate_Moist_Size();
    }

    /*! \brief get the number of  microphysics conserved non-moist (non-water, i.e., other vapor/condensed species) state variables */
    int Get_Qstate_NonMoist_Size () override
    {
        return m_moist_model->Qstate_NonMoist_Size();
    }

    /*! \brief initialize particles in particle container at level a_lev */
    inline void initParticles ( const int a_lev, /*!< AMR level */
                                std::unique_ptr<amrex::MultiFab>& z_phys_nd /*!< Nodal z heights at level a_lev */ )
    {
        m_moist_model->InitParticles(a_lev, z_phys_nd);
    }

    /*! \brief restart particles in particle container */
    inline void restartParticles ( amrex::ParGDBBase* a_gdb, const std::string& a_fname)
    {
        m_moist_model->RestartParticles(a_gdb, a_fname);
    }

    /*! \brief get the particle container from the moisture model */
    inline ERFPC* getParticleContainer () const
    {
        return m_moist_model->getParticleContainer();
    }

    /*! \brief get the name of the moisture model's particle container */
    inline const std::string& getName () const
    {
        return m_moist_model->getName();
    }

    /*! \brief get the indices and names of moisture model variables for restart
               at a given level */
    void Get_Qmoist_Restart_Vars ( const int a_lev, /*!< level */
                                   const SolverChoice& a_sc, /*!< Solver choice object */
                                   std::vector<int>& a_idx, /*!< indices */
                                   std::vector<std::string>& a_names /*!< names */ ) const override
    {
        if (a_lev == 0) {
            m_moist_model->Qmoist_Restart_Vars( a_sc, a_idx, a_names );
        } else {
            a_idx.clear();
            a_names.clear();
        }
    }

    /*! \brief Returns a list of additional plot variable names */
    virtual void GetPlotVarNames (amrex::Vector<std::string>& a_vec) const override
    {
        m_moist_model->GetPlotVarNames(a_vec);
    }

    /*! \brief Fills in a MultiFab for plotting */
    virtual void GetPlotVar (const std::string& a_name,
                             amrex::MultiFab& a_mf,
                             const int a_lev) const override
    {
        m_moist_model->SetCurrentLevel(a_lev);
        m_moist_model->GetPlotVar(a_name, a_mf, a_lev);
    }

    /*! \brief get the diagnostics interval */
    inline int getDiagnosticsInterval () const
    {
        return m_moist_model->getDiagnosticsInterval();
    }

     /*! \brief Populates dz_min in micro model for precipitation */
    virtual void Set_dzmin (const int /*a_lev*/,
                            const amrex::Real dz_min) const override
    {
        m_moist_model->Set_dzmin(dz_min);
    }

protected:

    /*! \brief Create and set the specified moisture model */
    template<class NewMoistModel>
    void SetModel ()
    {
        m_moist_model = std::make_unique<NewMoistModel>();
    }

    std::unique_ptr<NullMoistLagrangian> m_moist_model; /*!< moisture model */
};

#endif
#endif
