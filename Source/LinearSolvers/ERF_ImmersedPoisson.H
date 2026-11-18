/**
 * \file ERF_ImmersedPoisson.H
 */
#ifndef ERF_IMMERSED_POISSON_H_
#define ERF_IMMERSED_POISSON_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#ifdef ERF_USE_FFT
#include <AMReX_FFT_Poisson.H>

/**
 * Operator of the anelastic projection with the implicit immersed-forcing drag
 * (erf.if_implicit_projection), for GMRES on a constant-dz mesh:
 *
 *     L(phi) = div( sigma grad phi ),   sigma = 1 / (1 + dt * rate)  on each face,
 *
 * where rate is the drag rate of the terrain immersed forcing on that face. The constant
 * coefficient FFT Poisson solve (sigma = 1) is the preconditioner. The ghost cells of phi
 * follow the FFT boundary conditions of the subdomain, as in TerrainPoisson.
 */
class ImmersedPoisson
{
public:
    using RT = amrex::Real;

    /**
     * @param geom Geometry of the subdomain being solved
     * @param lev_geom Geometry of the level, used for the FFT boundary conditions
     * @param ba BoxArray of the subdomain
     * @param dm DistributionMapping of the subdomain
     * @param domain_bcs_type Domain boundary-condition names
     * @param sigma Face coefficients on the x-, y- and z-faces of the subdomain
     * @param use_real_bcs Whether real lateral boundary conditions are active
     */
    ImmersedPoisson (amrex::Geometry const& geom, amrex::Geometry const& lev_geom,
                     amrex::BoxArray const& ba, amrex::DistributionMapping const& dm,
                     amrex::Array<std::string,2*AMREX_SPACEDIM> const& domain_bcs_type,
                     amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> const& sigma,
                     bool use_real_bcs);

    //! lhs = div(sigma grad rhs)
    void apply (amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    //! Fill the ghost cells of phi from the FFT boundary conditions and the periodicity
    void apply_bcs (amrex::MultiFab& phi);

    //! Enable or disable the FFT preconditioner
    void usePrecond (bool use_precond_in) { m_use_precond = use_precond_in; }

    //! fluxes = -grad phi on each face (without sigma)
    void getFluxes (amrex::MultiFab& phi, amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fluxes);

    void assign (amrex::MultiFab& lhs, amrex::MultiFab const& rhs);
    void scale (amrex::MultiFab& lhs, amrex::Real fac);
    amrex::Real dotProduct (amrex::MultiFab const& v1, amrex::MultiFab const& v2);
    void increment (amrex::MultiFab& lhs, amrex::MultiFab const& rhs, amrex::Real a);
    void linComb (amrex::MultiFab& lhs, amrex::Real a, amrex::MultiFab const& rhs_a,
                  amrex::Real b, amrex::MultiFab const& rhs_b);
    amrex::MultiFab makeVecRHS ();
    amrex::MultiFab makeVecLHS ();
    amrex::Real norm2 (amrex::MultiFab const& v);
    //! lhs = the constant-coefficient FFT Poisson solve of rhs (or a copy without the preconditioner)
    void precond (amrex::MultiFab& lhs, amrex::MultiFab const& rhs);
    void setToZero (amrex::MultiFab& v);

private:
    bool m_use_precond = false;
    amrex::Geometry m_geom;
    amrex::BoxArray m_grids;
    amrex::DistributionMapping m_dmap;
    amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> m_sigma;
    std::unique_ptr<amrex::FFT::Poisson<amrex::MultiFab>> m_fft_precond;
    amrex::Array<std::pair<amrex::FFT::Boundary,amrex::FFT::Boundary>,AMREX_SPACEDIM> m_bc_fft;
};

#endif
#endif
