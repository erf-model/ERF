#ifndef ERF_OCEANSURF_H
#define ERF_OCEANSURF_H

#include <ERF_NullSurf.H>

#include <AMReX_Vector.H>

class OceanSurf : public NullSurf {
public:
    OceanSurf () = default;
    ~OceanSurf () override = default;

    void
    SetOceanState (const amrex::Vector<amrex::MultiFab*>* ocean_state)
    {
        m_ocean_state = ocean_state;
    }

    [[nodiscard]]
    const amrex::Vector<amrex::MultiFab*>* OceanState () const
    {
        return m_ocean_state;
    }

private:
    const amrex::Vector<amrex::MultiFab*>* m_ocean_state = nullptr;
};
#endif
