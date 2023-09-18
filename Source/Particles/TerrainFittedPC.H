#ifndef TERRAIN_FITTED_PC_H_
#define TERRAIN_FITTED_PC_H_

#include <AMReX_Particles.H>

struct RealIdx
{
    enum {
        vx = 0,
        vy, vz,
        ncomps
    };
};

struct IntIdx
{
    enum {
        k = 0,
        ncomps
    };
};

struct TerrainFittedAssignor
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

class TerrainFittedPC
    : public amrex::ParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0,
                                      amrex::DefaultAllocator, TerrainFittedAssignor>
{

public:

    TerrainFittedPC (amrex::ParGDBBase* gdb)
        : amrex::ParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, TerrainFittedAssignor>(gdb)
        {}

    TerrainFittedPC (const amrex::Geometry            & geom,
                     const amrex::DistributionMapping & dmap,
                     const amrex::BoxArray            & ba)
        : amrex::ParticleContainer<RealIdx::ncomps, IntIdx::ncomps, 0, 0,
                                   amrex::DefaultAllocator, TerrainFittedAssignor>(geom, dmap, ba)
        {}

    void InitParticles (const amrex::MultiFab& a_z_height);

    void AdvectWithUmac (amrex::MultiFab* umac, int level, amrex::Real dt,
                         const amrex::MultiFab& a_z_height);

    void Increment (amrex::MultiFab& mf, int level);
};

#endif
