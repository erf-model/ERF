/*! @file ERF_NullMoistLagrangian.H
 *  \brief Contains the Lagrangian moisture model base class */

#ifndef ERF_NULLMOISTLAGRANGIAN_H
#define ERF_NULLMOISTLAGRANGIAN_H

#ifdef ERF_USE_PARTICLES

#include <AMReX_AmrParGDB.H>
#include "ERF_IndexDefines.H"
#include "ERF_NullMoist.H"

/* forward declaration */
class ERFPC;

/*! \brief Base class for a Lagrangian moisture model
 *
 * Extends #NullMoist for a Lagrangian model */
class NullMoistLagrangian : public NullMoist {

public:

    /*! \brief Null constructor */
    NullMoistLagrangian () {}

    /*! \brief Default destructor */
    virtual ~NullMoistLagrangian () = default;

    /*! \brief Initialize particles at AMR level */
    virtual
    void InitParticles ( const int, std::unique_ptr<amrex::MultiFab>& ) { }

    /*! \brief Restart particles */
    virtual
    void RestartParticles ( amrex::ParGDBBase*, const std::string& ) { }

    /*! \brief get the particle container */
    virtual ERFPC* getParticleContainer ()
    {
        return nullptr;
    }

    /*! \brief get the name */
    virtual const std::string& getName () const
    {
        return m_name;
    }

    /*! \brief finish initializations steps that require flow variables */
    virtual void FinishInit (const int&, /* AMR level */
                             amrex::MultiFab&, /* Conserved state variable */
                             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& /*terrain*/) {}

    /*! \brief average down moisture multifabs from finest_level to 0 */
    virtual void AverageDownMicroVars (const int /*finest_level*/) {}

    void Advance ( const amrex::Real&, const SolverChoice& ) override
    {
        amrex::Abort("Do not use this advance() for Lagrangian microphysics.");
    }

    void SetCurrentLevel (const int&) override {}
    void InitLevel (const int, const amrex::MultiFab&) override {}
    int getDiagnosticsInterval () const override { return 0; }

    /*! \brief advance the moisture model by one time step */
    virtual void
    Advance (const amrex::Real&, /* dt */
             const int&, /* iter */
             const amrex::Real&, /* time */
             amrex::Vector<amrex::Vector<amrex::MultiFab>>&, /* state variables */
             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>&, /* terrain */
             const amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>& ) { }

protected:

private:
    const std::string m_name = "null"; /*!< name of the moisture model */
};

#endif
#endif
