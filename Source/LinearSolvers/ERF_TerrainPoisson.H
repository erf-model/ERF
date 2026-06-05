#ifndef ERF_TERRAIN_POISSON_H_
#define ERF_TERRAIN_POISSON_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#ifdef ERF_USE_FFT

#include "ERF_TerrainPoisson_3D_K.H"
#include <AMReX_FFT_Poisson.H>

class TerrainPoisson
{
public:

    using RT = amrex::Real;

    TerrainPoisson (amrex::Geometry const& geom, amrex::BoxArray const& ba,
                    amrex::DistributionMapping const& dm,
                    amrex::Array<std::string,2*AMREX_SPACEDIM>& domain_bcs_type,
                    amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_lev_d,
                    const amrex::MultiFab& ax, const amrex::MultiFab& ay,
                    const amrex::MultiFab& az, const amrex::MultiFab& dJ,
                    amrex::MultiFab const* z_phys_nd,
                    bool use_real_bcs);

    void apply(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    void apply_bcs(amrex::MultiFab& phi);

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
    amrex::Array<std::string,2*AMREX_SPACEDIM> m_domain_bcs_type;
    amrex::Gpu::DeviceVector<amrex::Real> m_stretched_dz_d;
    const amrex::MultiFab& m_ax;
    const amrex::MultiFab& m_ay;
    const amrex::MultiFab& m_az;
    const amrex::MultiFab& m_dJ;
    const amrex::MultiFab* m_zphys;

    std::unique_ptr<amrex::FFT::PoissonHybrid<amrex::MultiFab>> m_2D_fft_precond;
    amrex::Array<std::pair<amrex::FFT::Boundary,amrex::FFT::Boundary>,AMREX_SPACEDIM> bc_fft;
};

#endif
#endif
