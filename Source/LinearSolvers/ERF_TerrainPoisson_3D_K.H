/**
 * \file ERF_TerrainPoisson_3D_K.H
 */
#ifndef ERF_TERRAINPOISSON_3D_K_H_
#define ERF_TERRAINPOISSON_3D_K_H_

#include <AMReX_FArrayBox.H>
#include "ERF_Constants.H"

/**
 * Compute the x-face terrain Poisson flux.
 *
 * @tparam T Real-valued array and spacing type
 * @param i x-index of the x-face
 * @param j y-index of the x-face
 * @param k z-index of the x-face
 * @param sol Cell-centered solution field
 * @param zp Node-centered physical height field
 * @param dxinv Inverse x grid spacing
 * @return x-face flux for the terrain Poisson operator
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T terrpoisson_flux_x (int i, int j, int k,
                      amrex::Array4<T const> const& sol,
                      amrex::Array4<T const> const& zp,
                      T dxinv) noexcept
{
    using amrex::Real;

    Real h_xi, h_zeta;

    // On x-face
    Real px_lo = (sol(i,j,k) - sol(i-1,j,k)) * dxinv;

    // On y-edges
    Real pz_lo_md_hi  = myhalf * ( sol(i,j,k+1) + sol(i-1,j,k+1)
                                     -sol(i,j,k  ) - sol(i-1,j,k  ) );
         h_xi         = fourth * ( zp(i+1,j  ,k+1) - zp(i-1,j  ,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i-1,j+1,k+1) ) * dxinv;
         h_zeta       = fourth * ( zp(i,j  ,k+2) - zp(i  ,j  ,k  )
                                      +zp(i,j+1,k+2) - zp(i  ,j+1,k  ) );
         pz_lo_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = myhalf * ( sol(i,j,k  ) + sol(i-1,j,k  )
                                     -sol(i,j,k-1) - sol(i-1,j,k-1) );

         h_xi         = fourth * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-face
    px_lo -= myhalf * ( pz_lo_md_hi + pz_lo_md_lo );

    return -px_lo;
}

/**
 * Compute the y-face terrain Poisson flux.
 *
 * @tparam T Real-valued array and spacing type
 * @param i x-index of the y-face
 * @param j y-index of the y-face
 * @param k z-index of the y-face
 * @param sol Cell-centered solution field
 * @param zp Node-centered physical height field
 * @param dyinv Inverse y grid spacing
 * @return y-face flux for the terrain Poisson operator
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T terrpoisson_flux_y (int i, int j, int k,
                      amrex::Array4<T const> const& sol,
                      amrex::Array4<T const> const& zp,
                      T dyinv) noexcept
{
    using amrex::Real;

    Real h_eta, h_zeta;

    Real py_lo = (sol(i,j,k) - sol(i,j-1,k)) * dyinv;

    // On x-edges
    Real pz_md_lo_hi  = myhalf * ( sol(i,j,k+1) + sol(i,j-1,k+1)
                                     -sol(i,j,k  ) - sol(i,j-1,k  ) );
         h_eta        = fourth * ( zp(i  ,j+1,k+1) - zp(i  ,j-1,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j-1,k+1) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j,k+2) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+2) - zp(i+1,j  ,k  ) );
         pz_md_lo_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = myhalf * ( sol(i,j,k  ) + sol(i,j-1,k  )
                                     -sol(i,j,k-1) - sol(i,j-1,k-1) );

         h_eta        = fourth * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k  )
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On y-face
    py_lo -= myhalf * ( pz_md_lo_hi + pz_md_lo_lo );
    return -py_lo;
}

/**
 * Compute the z-face terrain Poisson flux.
 *
 * @tparam T Real-valued array and spacing type
 * @param i x-index of the z-face
 * @param j y-index of the z-face
 * @param k z-index of the z-face
 * @param sol Cell-centered solution field
 * @param zp Node-centered physical height field
 * @param dxinv Inverse x grid spacing
 * @param dyinv Inverse y grid spacing
 * @return z-face flux for the terrain Poisson operator
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T terrpoisson_flux_z (int i, int j, int k,
                      amrex::Array4<T const> const& sol,
                      amrex::Array4<T const> const& zp,
                      T dxinv, T dyinv) noexcept
{
    using amrex::Real;

    Real h_xi, h_eta, h_zeta;

    // On z-face
    Real pz_lo = (sol(i,j,k  ) - sol(i,j,k-1));
    Real hzeta_inv_on_zlo = amrex::Real(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
    pz_lo *= hzeta_inv_on_zlo;

    // On corners
    Real px_hi_md_lo = myhalf * ( sol(i+1,j  ,k  ) - sol(i  ,j  ,k  )
                                    +sol(i+1,j  ,k-1) - sol(i  ,j  ,k-1)) * dxinv;
    Real px_lo_md_lo = myhalf * ( sol(i  ,j  ,k  ) - sol(i-1,j  ,k  )
                                    +sol(i  ,j  ,k-1) - sol(i-1,j  ,k-1)) * dxinv;
    Real py_md_hi_lo = myhalf * ( sol(i  ,j+1,k  ) - sol(i  ,j  ,k  )
                                    +sol(i  ,j+1,k-1) - sol(i  ,j  ,k-1)) * dyinv;
    Real py_md_lo_lo = myhalf * ( sol(i  ,j  ,k  ) - sol(i  ,j-1,k  )
                                    +sol(i  ,j  ,k-1) - sol(i  ,j-1,k-1)) * dyinv;

    // On y-edges
    Real pz_hi_md_lo  = myhalf * ( sol(i+1,j,k  ) + sol(i,j,k  )
                                     -sol(i+1,j,k-1) - sol(i,j,k-1) );
         h_xi         = fourth * ( zp(i+2,j  ,k  ) - zp(i  ,j  ,k  )
                                      +zp(i+2,j+1,k  ) - zp(i  ,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = myhalf * ( sol(i,j,k  ) + sol(i-1,j,k  )
                                     -sol(i,j,k-1) - sol(i-1,j,k-1) );
         h_xi         = fourth * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-edges
    Real pz_md_hi_lo  = myhalf * ( sol(i,j+1,k  ) + sol(i,j,k  )
                                     -sol(i,j+1,k-1) - sol(i,j,k-1) );
         h_eta        = fourth * ( zp(i  ,j+2,k  ) - zp(i  ,j  ,k  )
                                      +zp(i+1,j+2,k  ) - zp(i+1,j  ,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = myhalf * ( sol(i,j,k  ) + sol(i,j-1,k  )
                                     -sol(i,j,k-1) - sol(i,j-1,k-1) );
         h_eta        = fourth * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k  )
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On z-face
    Real h_xi_on_zlo  = myhalf * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k)) * dxinv;
    Real h_eta_on_zlo = myhalf * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k)) * dyinv;

    pz_lo -= myhalf * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= myhalf * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );

    if (k == 0) pz_lo = zero;

    return -pz_lo;
}

/**
 * Simplified gradient calc with dirichlet BC (phi_zlo==0)
 * - compare with terrpoisson_flux_z
 *
 * @tparam T Real-valued array and spacing type
 * @param i x-index of the z-face
 * @param j y-index of the z-face
 * @param k z-index of the z-face
 * @param sol Cell-centered solution field
 * @param zp Node-centered physical height field
 * @param dxinv Inverse x grid spacing
 * @param dyinv Inverse y grid spacing
 * @return z-face flux with a low-z Dirichlet value
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T terrpoisson_flux_zlo_dir (int i, int j, int k,
                            amrex::Array4<T const> const& sol,
                            amrex::Array4<T const> const& zp,
                            T dxinv, T dyinv) noexcept
{
    using amrex::Real;

    Real h_xi, h_eta, h_zeta;

    // On z-face, one-sided
    Real pz_lo = sol(i,j,k);
    Real hzeta_inv_on_zlo = amrex::Real(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
    pz_lo *= hzeta_inv_on_zlo;

    // On corners
    constexpr Real px_hi_md_lo = zero;
    constexpr Real px_lo_md_lo = zero;
    constexpr Real py_md_hi_lo = zero;
    constexpr Real py_md_lo_lo = zero;

    // On y-edges, one-sided
    Real pz_hi_md_lo  = myhalf * ( sol(i+1,j,k  ) + sol(i,j,k  ));
         h_xi         = fourth * ( zp(i+2,j  ,k  ) - zp(i  ,j  ,k  )
                                      +zp(i+2,j+1,k  ) - zp(i  ,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k  )
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k  ) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On y-edges, one-sided
    Real pz_lo_md_lo  = myhalf * ( sol(i,j,k  ) + sol(i-1,j,k  ));
         h_xi         = fourth * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i,j  ,k+1) - zp(i  ,j  ,k  )
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k  ) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-edges, one-sided
    Real pz_md_hi_lo  = myhalf * ( sol(i,j+1,k  ) + sol(i,j,k  ));
         h_eta        = fourth * ( zp(i  ,j+2,k  ) - zp(i  ,j  ,k  )
                                      +zp(i+1,j+2,k  ) - zp(i+1,j  ,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k  )
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k  ) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On x-edges, one-sided
    Real pz_md_lo_lo  = myhalf * ( sol(i,j,k  ) + sol(i,j-1,k  ));
         h_eta        = fourth * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k  )
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j,k+1) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k  ) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On z-face
    Real h_xi_on_zlo  = myhalf * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k)) * dxinv;
    Real h_eta_on_zlo = myhalf * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k)) * dyinv;

    pz_lo -= myhalf * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= myhalf * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );

    return -pz_lo;
}

/**
 * Apply the terrain Poisson matrix to a cell-centered field.
 *
 * @tparam T Real-valued array and spacing type
 * @param i x-index of the cell
 * @param j y-index of the cell
 * @param k z-index of the cell
 * @param y Destination operator result
 * @param x Source solution field
 * @param ax Terrain metric coefficient on x-faces
 * @param ay Terrain metric coefficient on y-faces
 * @param az Terrain metric coefficient on z-faces
 * @param dJ Cell-centered Jacobian determinant
 * @param zp Node-centered physical height field
 * @param dxinv Inverse x grid spacing
 * @param dyinv Inverse y grid spacing
 * @param dzinv Inverse z grid spacing
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_adotx (int i, int j, int k, amrex::Array4<T> const& y,
                        amrex::Array4<T const> const& x,
                        amrex::Array4<T const> const& ax,
                        amrex::Array4<T const> const& ay,
                        amrex::Array4<T const> const& az,
                        amrex::Array4<T const> const& dJ,
                        amrex::Array4<T const> const& zp,
                        T dxinv, T dyinv, T dzinv) noexcept
{
    using amrex::Real;
    Real h_xi, h_eta, h_zeta;

#if 1
    // *********************************************************
    // Hi x-face
    // *********************************************************
    // On x-face
    Real px_hi = (x(i+1,j,k) - x(i,j,k)) * dxinv;

    // On y-edges
    Real pz_hi_md_hi  = myhalf  * (  x(i+1,j  ,k+1) +  x(i  ,j  ,k+1)
                                    -x(i+1,j  ,k  ) -  x(i  ,j  ,k  ) );
         h_xi         = fourth * ( zp(i+2,j  ,k+1) - zp(i  ,j  ,k+1)
                                  +zp(i+2,j+1,k+1) - zp(i  ,j+1,k+1) ) * dxinv;
         h_zeta       = fourth * ( zp(i+1,j  ,k+2) - zp(i+1,j  ,k  )
                                  +zp(i+1,j+1,k+2) - zp(i+1,j+1,k  ) );
         pz_hi_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_hi_md_lo  = myhalf  *  ( x(i+1,j  ,k  ) +  x(i  ,j  ,k  )
                                    -x(i+1,j  ,k-1) -  x(i  ,j  ,k-1) );
         h_xi         = fourth * ( zp(i+2,j  ,k  ) - zp(i  ,j  ,k  )
                                  +zp(i+2,j+1,k  ) - zp(i  ,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k-1)
                                  +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On x-face
    px_hi -= myhalf * ( pz_hi_md_hi + pz_hi_md_lo );

    // *********************************************************
    // Lo x-face
    // ********************************************************* // On x-face
    Real px_lo = (x(i,j,k) - x(i-1,j,k)) * dxinv;

    // On y-edges
    Real pz_lo_md_hi  = myhalf  * (  x(i,j,k+1) + x(i-1,j,k+1)
                                    -x(i,j,k  ) - x(i-1,j,k  ) );
         h_xi         = fourth * ( zp(i+1,j  ,k+1) - zp(i-1,j  ,k+1)
                                  +zp(i+1,j+1,k+1) - zp(i-1,j+1,k+1) ) * dxinv;
         h_zeta       = fourth * ( zp(i,j  ,k+2) - zp(i  ,j  ,k  )
                                  +zp(i,j+1,k+2) - zp(i  ,j+1,k  ) );
         pz_lo_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = myhalf  * (  x(i,j,k  ) + x(i-1,j,k  )
                                    -x(i,j,k-1) - x(i-1,j,k-1) );
         h_xi         = fourth * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                  +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) ) * dxinv;
         h_zeta       = fourth * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                  +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-face
    px_lo -= myhalf * ( pz_lo_md_hi + pz_lo_md_lo );

    // *********************************************************
    // Hi y-face
    // *********************************************************
    // On y-face
    Real py_hi = (x(i,j+1,k) - x(i,j,k)) * dyinv;

    // On x-edges
    Real pz_md_hi_hi  = myhalf  * (  x(i,j+1,k+1) + x(i,j,k+1)
                                       -x(i,j+1,k  ) - x(i,j,k  ) );
         h_eta        = fourth * ( zp(i  ,j+2,k+1) - zp(i  ,j  ,k+1)
                                      +zp(i+1,j+2,k+1) - zp(i+1,j  ,k+1) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j+1,k+2) - zp(i  ,j+1,k  )
                                      +zp(i+1,j+1,k+2) - zp(i+1,j+1,k  ) );
         pz_md_hi_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_hi_lo  = myhalf  * (  x(i,j+1,k  ) + x(i,j,k  )
                                       -x(i,j+1,k-1) - x(i,j,k-1) );
         h_eta        = fourth * ( zp(i  ,j+2,k  ) - zp(i  ,j  ,k  )
                                      +zp(i+1,j+2,k  ) - zp(i+1,j  ,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On y-face
    py_hi -= myhalf * ( pz_md_hi_hi + pz_md_hi_lo );

    // *********************************************************
    // Lo y-face
    // *********************************************************
    // On y-face
    Real py_lo = (x(i,j,k) - x(i,j-1,k)) * dyinv;

    // On x-edges
    Real pz_md_lo_hi  = myhalf  * (  x(i  ,j,k+1) +  x(i  ,j-1,k+1)
                                       -x(i  ,j,k  ) -  x(i  ,j-1,k  ) );
         h_eta        = fourth * ( zp(i  ,j+1,k+1) - zp(i  ,j-1,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j-1,k+1) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j,k+2) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+2) - zp(i+1,j  ,k  ) );
         pz_md_lo_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = myhalf  * (  x(i  ,j,k  ) +  x(i  ,j-1,k  )
                                       -x(i  ,j,k-1) -  x(i  ,j-1,k-1) );
         h_eta        = fourth * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k  )
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k  ) ) * dyinv;
         h_zeta       = fourth * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On y-face
    py_lo -= myhalf * ( pz_md_lo_hi + pz_md_lo_lo );

    // *********************************************************
    // Hi z-face
    // *********************************************************
    // On z-face
    Real pz_hi    = x(i,j,k+1) - x(i,j,k  );
    Real hzeta_inv_on_zhi = amrex::Real(8.0) / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
    pz_hi *= hzeta_inv_on_zhi;

    // On corners
    Real px_hi_md_hi = myhalf * ( x(i+1,j,k+1) - x(i  ,j  ,k+1)
                                    +x(i+1,j,k  ) - x(i  ,j  ,k  )) * dxinv;
    Real px_lo_md_hi = myhalf * ( x(i  ,j,k+1) - x(i-1,j  ,k+1)
                                    +x(i  ,j,k  ) - x(i-1,j  ,k  )) * dxinv;
    Real py_md_hi_hi = myhalf * ( x(i,j+1,k+1) - x(i  ,j  ,k+1)
                                    +x(i,j+1,k  ) - x(i  ,j  ,k  )) * dyinv;
    Real py_md_lo_hi = myhalf * ( x(i,j  ,k+1) - x(i  ,j-1,k+1)
                                    +x(i,j  ,k  ) - x(i  ,j-1,k  )) * dyinv;

    // On z-face
    Real h_xi_on_zhi  = myhalf * ( zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1) ) * dxinv;
    Real h_eta_on_zhi = myhalf * ( zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1) ) * dyinv;
    //
    // Note we do not need to recalculate pz_...hi here
    //
    pz_hi -= myhalf * h_xi_on_zhi  * ( (px_hi_md_hi + px_lo_md_hi) - (pz_hi_md_hi + pz_lo_md_hi) );
    pz_hi -= myhalf * h_eta_on_zhi * ( (py_md_hi_hi + py_md_lo_hi) - (pz_md_hi_hi + pz_md_lo_hi) );

    // *********************************************************
    // Lo z-face
    // *********************************************************
    // On z-face
    Real pz_lo = x(i,j,k  ) - x(i,j,k-1);
    Real hzeta_inv_on_zlo = amrex::Real(8.0) / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
    pz_lo *= hzeta_inv_on_zlo;

    // On corners
    Real px_hi_md_lo = myhalf * ( x(i+1,j  ,k  ) - x(i  ,j  ,k  )
                                    +x(i+1,j  ,k-1) - x(i  ,j  ,k-1)) * dxinv;
    Real px_lo_md_lo = myhalf * ( x(i  ,j  ,k  ) - x(i-1,j  ,k  )
                                    +x(i  ,j  ,k-1) - x(i-1,j  ,k-1)) * dxinv;
    Real py_md_hi_lo = myhalf * ( x(i  ,j+1,k  ) - x(i  ,j  ,k  )
                                    +x(i  ,j+1,k-1) - x(i  ,j  ,k-1)) * dyinv;
    Real py_md_lo_lo = myhalf * ( x(i  ,j  ,k  ) - x(i  ,j-1,k  )
                                    +x(i  ,j  ,k-1) - x(i  ,j-1,k-1)) * dyinv;

    // On z-face
    Real h_xi_on_zlo  = myhalf * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k)) * dxinv;
    Real h_eta_on_zlo = myhalf * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k)) * dyinv;
    //
    // Note we do not need to recalculate pz_...lo here
    //
    pz_lo -= myhalf * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= myhalf * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );

#else

    // *********************************************************
    // Version which calls flux routines
    // *********************************************************
    //
    // This option uses calls to the flux routines so there is
    //      some duplicated computation
    // This option should give the same answer as above
    //
        Real px_lo = -terrpoisson_flux_x(i  ,j,k,x,zp,dxinv);
        Real px_hi = -terrpoisson_flux_x(i+1,j,k,x,zp,dxinv);
        Real py_lo = -terrpoisson_flux_y(i,j  ,k,x,zp,dyinv);
        Real py_hi = -terrpoisson_flux_y(i,j+1,k,x,zp,dyinv);
        Real pz_lo = -terrpoisson_flux_z(i,j,k  ,x,zp,dxinv,dyinv);
        Real pz_hi = -terrpoisson_flux_z(i,j,k+1,x,zp,dxinv,dyinv);
#endif
    //
    // *********************************************************
    // Adotx
    // *********************************************************
    if (k == 0) pz_lo = zero;

    y(i,j,k) =  (ax(i+1,j,k)*px_hi - ax(i,j,k)*px_lo) * dxinv
               +(ay(i,j+1,k)*py_hi - ay(i,j,k)*py_lo) * dyinv
               +(az(i,j,k+1)*pz_hi - az(i,j,k)*pz_lo) * dzinv;
    y(i,j,k) /=  dJ(i,j,k);
}
#endif
