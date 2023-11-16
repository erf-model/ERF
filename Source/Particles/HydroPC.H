#ifndef HYDRO_PC_H_
#define HYDRO_PC_H_

#include <AMReX_Particles.H>

struct HydroRealIdx
{
    enum {
        vx = 0,
        vy, vz,
        mass,
        ncomps
    };
};

struct HydroIntIdx
{
    enum {
        k = 0,
        ncomps
    };
};

struct HydroAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (P const& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               const amrex::Box& domain) const noexcept
    {
        amrex::IntVect iv(
            AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                         int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                         p.idata(0)));
        iv[0] += domain.smallEnd()[0];
        iv[1] += domain.smallEnd()[1];
        return iv;
    }
};

class HydroPC
    : public amrex::ParticleContainer<HydroRealIdx::ncomps, HydroIntIdx::ncomps, 0, 0,
                                      amrex::DefaultAllocator, HydroAssignor>
{

public:

    HydroPC (amrex::ParGDBBase* gdb)
        : amrex::ParticleContainer<HydroRealIdx::ncomps, HydroIntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, HydroAssignor>(gdb)
        {}

    HydroPC (const amrex::Geometry            & geom,
              const amrex::DistributionMapping & dmap,
              const amrex::BoxArray            & ba)
        : amrex::ParticleContainer<HydroRealIdx::ncomps, HydroIntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, HydroAssignor>(geom, dmap, ba)
        {}

    void InitParticles (const amrex::MultiFab& a_z_height);

    void AdvectWithGravity (int level, amrex::Real dt, const amrex::MultiFab& a_z_height);
};

#endif
