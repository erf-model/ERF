#ifndef NULLMOISTLAGRANGIAN_H
#define NULLMOISTLAGRANGIAN_H

#ifdef ERF_USE_PARTICLES

#include "NullMoist.H"

/* forward declaration */
class ERFPC;

class NullMoistLagrangian : public NullMoist {

public:
    NullMoistLagrangian() {}

    virtual ~NullMoistLagrangian() = default;

    virtual ERFPC* getParticleContainer()
    {
        return nullptr;
    }

    virtual const std::string& getName() const
    {
        return m_name;
    }

    using NullMoist::Advance;

    virtual void
    Advance (const amrex::Real&,
             amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
             const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& ) { }

protected:

private:
    const std::string m_name = "null";
};

#endif
#endif
