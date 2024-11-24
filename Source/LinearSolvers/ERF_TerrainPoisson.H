#ifndef ERF_TERRAIN_POISSON_H_
#define ERF_TERRAIN_POISSON_H_

#include "ERF_TerrainPoisson_3D_K.H"
#ifdef ERF_USE_FFT
#include "ERF_FFT_TerrainPrecond.H"
#endif

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

class TerrainPoisson
{
public:

    using RT = amrex::Real;

    TerrainPoisson (amrex::Geometry const& geom, amrex::BoxArray const& ba,
                    amrex::DistributionMapping const& dm,
                    amrex::MultiFab const* z_phys_nd);

    void apply(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    void usePrecond(bool precond_flag);

    void getFluxes(amrex::MultiFab& phi, amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fluxes);

    void assign(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    void scale(amrex::MultiFab& lhs, amrex::Real fac);

    amrex::Real dotProduct(amrex::MultiFab const& v1, amrex::MultiFab const& v2);

    void increment(amrex::MultiFab& lhs, amrex::MultiFab const& rhs, amrex::Real a);

    void linComb(amrex::MultiFab& lhs, amrex::Real a, amrex::MultiFab const& rhs_a,
                 amrex::Real b, amrex::MultiFab const& rhs_b);

    amrex::MultiFab makeVecRHS();

    amrex::MultiFab makeVecLHS();

    amrex::Real norm2(amrex::MultiFab const& v);

    void precond(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    void setToZero(amrex::MultiFab& v);

private:
    bool m_use_precond = false;
    amrex::Geometry m_geom;
    amrex::BoxArray m_grids;
    amrex::DistributionMapping m_dmap;
    const amrex::MultiFab* m_zphys;
#ifdef ERF_USE_FFT
    std::unique_ptr<amrex::FFT::PoissonTerrainPrecond<amrex::MultiFab>> m_2D_fft_precond;
#endif
};

#endif
