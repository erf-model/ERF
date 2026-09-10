/**
 * \file ERF_TerrainPoisson.H
 */
#ifndef ERF_TERRAIN_POISSON_H_
#define ERF_TERRAIN_POISSON_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

#ifdef ERF_USE_FFT

#include "ERF_TerrainPoisson_3D_K.H"
#include <AMReX_FFT_Poisson.H>

/**
 * Matrix-vector interface for the terrain Poisson operator used by GMRES.
 */
class TerrainPoisson
{
public:

    using RT = amrex::Real;

    /**
     * Construct a terrain Poisson operator and its FFT preconditioner metadata.
     *
     * @param geom Geometry for the subdomain being solved
     * @param ba BoxArray for operator and vector allocation
     * @param dm DistributionMapping for operator and vector allocation
     * @param domain_bcs_type Domain boundary-condition names
     * @param stretched_dz_lev_d Device stretched-dz spacing for the level
     * @param ax Terrain metric coefficient on x-faces
     * @param ay Terrain metric coefficient on y-faces
     * @param az Terrain metric coefficient on z-faces
     * @param dJ Cell-centered Jacobian determinant
     * @param z_phys_nd Node-centered physical height field
     * @param use_real_bcs Whether real lateral boundary conditions are active
     */
    TerrainPoisson (amrex::Geometry const& geom, amrex::Geometry const& lev_geom,
                    amrex::BoxArray const& ba,
                    amrex::DistributionMapping const& dm,
                    amrex::Array<std::string,2*AMREX_SPACEDIM>& domain_bcs_type,
                    amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_lev_d,
                    const amrex::MultiFab& ax, const amrex::MultiFab& ay,
                    const amrex::MultiFab& az, const amrex::MultiFab& dJ,
                    amrex::MultiFab const* z_phys_nd,
                    bool use_real_bcs);

    /**
     * Apply the terrain Poisson operator.
     *
     * @param lhs Destination for the operator result
     * @param rhs Input vector with ghost cells available for boundary fills
     */
    void apply(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    /**
     * Fill ghost cells according to the configured Poisson boundary conditions.
     *
     * @param phi Field whose ghost cells are updated
     */
    void apply_bcs(amrex::MultiFab& phi);

    /**
     * Enable or disable the FFT preconditioner.
     *
     * @param precond_flag True to use the FFT preconditioner
     */
    void usePrecond(bool precond_flag);

    /**
     * Compute face-centered terrain Poisson fluxes from a solution field.
     *
     * @param phi Solution field used to compute gradients
     * @param fluxes Face-centered fluxes to fill
     */
    void getFluxes(amrex::MultiFab& phi, amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& fluxes);

    /**
     * Copy one vector into another.
     *
     * @param lhs Destination vector
     * @param rhs Source vector
     */
    void assign(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    /**
     * Scale a vector in place.
     *
     * @param lhs Vector to scale
     * @param fac Multiplicative scale factor
     */
    void scale(amrex::MultiFab& lhs, amrex::Real fac);

    /**
     * Compute the vector dot product.
     *
     * @param v1 First input vector
     * @param v2 Second input vector
     * @return Dot product of v1 and v2
     */
    amrex::Real dotProduct(amrex::MultiFab const& v1, amrex::MultiFab const& v2);

    /**
     * Add a scaled vector into another vector.
     *
     * @param lhs Destination vector to increment
     * @param rhs Source vector
     * @param a Scale factor applied to rhs
     */
    void increment(amrex::MultiFab& lhs, amrex::MultiFab const& rhs, amrex::Real a);

    /**
     * Compute a linear combination of two source vectors.
     *
     * @param lhs Destination vector
     * @param a Scale factor applied to rhs_a
     * @param rhs_a First source vector
     * @param b Scale factor applied to rhs_b
     * @param rhs_b Second source vector
     */
    void linComb(amrex::MultiFab& lhs, amrex::Real a, amrex::MultiFab const& rhs_a,
                 amrex::Real b, amrex::MultiFab const& rhs_b);

    /**
     * Allocate a vector with the right-hand-side ghost-cell layout.
     *
     * @return Newly allocated right-hand-side vector
     */
    amrex::MultiFab makeVecRHS();

    /**
     * Allocate a vector with the solution ghost-cell layout.
     *
     * @return Newly allocated solution vector
     */
    amrex::MultiFab makeVecLHS();

    /**
     * Compute the L2 norm of a vector.
     *
     * @param v Input vector
     * @return L2 norm of v
     */
    amrex::Real norm2(amrex::MultiFab const& v);

    /**
     * Apply the configured preconditioner.
     *
     * @param lhs Destination for the preconditioned result
     * @param rhs Input vector to precondition
     */
    void precond(amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    /**
     * Set a vector to zero.
     *
     * @param v Vector to zero
     */
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
