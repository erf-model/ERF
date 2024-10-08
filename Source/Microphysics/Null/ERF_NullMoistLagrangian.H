/*! @file ERF_NullMoistLagrangian.H
 *  \brief Contains the Lagrangian moisture model base class */

#ifndef ERF_NULLMOISTLAGRANGIAN_H
#define ERF_NULLMOISTLAGRANGIAN_H

#ifdef ERF_USE_PARTICLES

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

    /*! \brief Initialize particles */
    virtual
    void InitParticles ( std::unique_ptr<amrex::MultiFab>& ) { }

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

    using NullMoist::Advance;

    /*! \brief advance the moisture model by one time step */
    virtual void
    Advance (const amrex::Real&, /* dt */
             const int&, /* iter */
             const amrex::Real&, /* time */
             amrex::Vector<amrex::Vector<amrex::MultiFab>>&, /* state variables */
             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& /* terrain */) { }

protected:

private:
    const std::string m_name = "null"; /*!< name of the moisture model */
};

#endif
#endif
