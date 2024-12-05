#ifndef ERF_TERRAINPOISSON_3D_K_H_
#define ERF_TERRAINPOISSON_3D_K_H_

#include <AMReX_FArrayBox.H>

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
    Real pz_lo_md_hi  = Real(0.5) * ( sol(i,j,k+1) + sol(i-1,j,k+1)
                                     -sol(i,j,k  ) - sol(i-1,j,k  ) );
         h_xi         = Real(0.25) * ( zp(i,j  ,k+1) - zp(i-2,j  ,k+1)
                                      +zp(i,j+1,k+1) - zp(i-2,j+1,k+1) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+2) - zp(i  ,j  ,k  )
                                      +zp(i,j+1,k+2) - zp(i  ,j+1,k  ) );
         pz_lo_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i-1,j,k  )
                                     -sol(i,j,k-1) - sol(i-1,j,k-1) );

         h_xi         = Real(0.25) * ( zp(i,j  ,k  ) - zp(i-2,j  ,k  )
                                      +zp(i,j+1,k  ) - zp(i-2,j+1,k  ) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-face
    px_lo -= Real(0.5) * ( pz_lo_md_hi + pz_lo_md_lo );

    return -px_lo;
}
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
    Real pz_md_lo_hi  = Real(0.5) * ( sol(i,j,k+1) + sol(i,j-1,k+1)
                                     -sol(i,j,k  ) - sol(i,j-1,k  ) );
         h_eta        = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j-2,k+1)
                                      +zp(i+1,j,k+1) - zp(i+1,j-2,k+1) ) * dyinv;
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+2) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+2) - zp(i+1,j  ,k  ) );
         pz_md_lo_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i,j-1,k  )
                                     -sol(i,j,k-1) - sol(i,j-1,k-1) );

         h_eta        = Real(0.25) * ( zp(i  ,j,k  ) - zp(i  ,j-2,k  )
                                      +zp(i+1,j,k  ) - zp(i+1,j-2,k  ) ) * dyinv;
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On y-face
    py_lo -= Real(0.5) * ( pz_md_lo_hi + pz_md_lo_lo );

    return -py_lo;
}

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
    Real hzeta_inv_on_zlo = 8.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
    pz_lo *= hzeta_inv_on_zlo;

    // On corners
    Real px_hi_md_lo = Real(0.5) * ( sol(i+1,j  ,k  ) - sol(i  ,j  ,k  )
                                    +sol(i+1,j  ,k-1) - sol(i  ,j  ,k-1)) * dxinv;
    Real px_lo_md_lo = Real(0.5) * ( sol(i  ,j  ,k  ) - sol(i-1,j  ,k  )
                                    +sol(i  ,j  ,k-1) - sol(i-1,j  ,k-1)) * dxinv;
    Real py_md_hi_lo = Real(0.5) * ( sol(i  ,j+1,k  ) - sol(i  ,j  ,k  )
                                    +sol(i  ,j+1,k-1) - sol(i  ,j  ,k-1)) * dyinv;
    Real py_md_lo_lo = Real(0.5) * ( sol(i  ,j  ,k  ) - sol(i  ,j-1,k  )
                                    +sol(i  ,j  ,k-1) - sol(i  ,j-1,k-1)) * dyinv;

    // On y-edges
    Real pz_hi_md_lo  = Real(0.5) * ( sol(i+1,j,k  ) + sol(i,j,k  )
                                     -sol(i+1,j,k-1) - sol(i,j,k-1) );
         h_xi         = Real(0.25) * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i-1,j,k  )
                                     -sol(i,j,k-1) - sol(i-1,j,k-1) );
         h_xi         = Real(0.25) * ( zp(i,j  ,k  ) - zp(i-2,j  ,k  )
                                      +zp(i,j+1,k  ) - zp(i-2,j+1,k  ) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-edges
    Real pz_md_hi_lo  = Real(0.5) * ( sol(i,j+1,k  ) + sol(i,j,k  )
                                     -sol(i,j+1,k-1) - sol(i,j,k-1) );
         h_eta        = Real(0.25) * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k)
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k) ) * dyinv;
         h_zeta       = Real(0.25) * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i,j-1,k  )
                                     -sol(i,j,k-1) - sol(i,j-1,k-1) );
         h_eta        = Real(0.25) * ( zp(i  ,j,k  ) - zp(i  ,j-2,k  )
                                      +zp(i+1,j,k  ) - zp(i+1,j-2,k  ) ) * dyinv;
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On z-face
    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k)) * dxinv;
    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k)) * dyinv;

    pz_lo -= 0.5 * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= 0.5 * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );

    return -pz_lo;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_adotx (int i, int j, int k, amrex::Array4<T> const& y,
                        amrex::Array4<T const> const& x,
                        amrex::Array4<T const> const& zp,
                        T dxinv, T dyinv) noexcept
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
    Real pz_hi_md_hi  = Real(0.25) * ( x(i+1,j,k+1) + x(i,j,k+1)
                                      -x(i+1,j,k  ) - x(i,j,k  ) );
         h_xi         = Real(0.25) * ( zp(i+1,j  ,k+1) - zp(i-1,j  ,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i-1,j+1,k+1) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i+1,j  ,k+2) - zp(i+1,j  ,k  )
                                      +zp(i+1,j+1,k+2) - zp(i+1,j+1,k  ) );
         pz_hi_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_hi_md_lo  = Real(0.5) * ( x(i+1,j,k  ) + x(i,j,k  )
                                     -x(i+1,j,k-1) - x(i,j,k-1) );
         h_xi         = Real(0.25) * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On x-face
    px_hi -= Real(0.5) * ( pz_hi_md_hi + pz_hi_md_lo );

    // *********************************************************
    // Lo x-face
    // ********************************************************* // On x-face
    Real px_lo = (x(i,j,k) - x(i-1,j,k)) * dxinv;

    // On y-edges
    Real pz_lo_md_hi  = Real(0.5) * ( x(i,j,k+1) + x(i-1,j,k+1)
                                     -x(i,j,k  ) - x(i-1,j,k  ) );
         h_xi         = Real(0.25) * ( zp(i,j  ,k+1) - zp(i-2,j  ,k+1)
                                      +zp(i,j+1,k+1) - zp(i-2,j+1,k+1) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+2) - zp(i  ,j  ,k  )
                                      +zp(i,j+1,k+2) - zp(i  ,j+1,k  ) );
         pz_lo_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = Real(0.5) * ( x(i,j,k  ) + x(i-1,j,k  )
                                     -x(i,j,k-1) - x(i-1,j,k-1) );

         h_xi         = Real(0.25) * ( zp(i,j  ,k  ) - zp(i-2,j  ,k  )
                                      +zp(i,j+1,k  ) - zp(i-2,j+1,k  ) ) * dxinv;
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-face
    px_lo -= Real(0.5) * ( pz_lo_md_hi + pz_lo_md_lo );

    // *********************************************************
    // Hi y-face
    // *********************************************************
    // On y-face
    Real py_hi = (x(i,j+1,k) - x(i,j,k)) * dyinv;

    // On x-edges
    Real pz_md_hi_hi  = Real(0.25) * ( x(i,j+1,k+1) + x(i,j,k+1)
                                      -x(i,j+1,k  ) - x(i,j,k  ) );
         h_eta        = Real(0.25) * ( zp(i  ,j+1,k+1) - zp(i  ,j-1,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j-1,k+1) ) * dyinv;

         h_zeta       = Real(0.25) * ( zp(i  ,j+1,k+2) - zp(i  ,j+1,k  )
                                      +zp(i+1,j+1,k+2) - zp(i+1,j+1,k  ) );
         pz_md_hi_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_hi_lo  = Real(0.5) * ( x(i,j+1,k  ) + x(i,j,k  )
                                     -x(i,j+1,k-1) - x(i,j,k-1) );

         h_eta        = Real(0.25) * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k)
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k) ) * dyinv;

         h_zeta       = Real(0.25) * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On y-face
    py_hi -= Real(0.5) * ( pz_md_hi_hi + pz_md_hi_lo );

    // *********************************************************
    // Lo y-face
    // *********************************************************
    // On y-face
    Real py_lo = (x(i,j,k) - x(i,j-1,k)) * dyinv;

    // On x-edges
    Real pz_md_lo_hi  = Real(0.5) * ( x(i,j,k+1) + x(i,j-1,k+1)
                                     -x(i,j,k  ) - x(i,j-1,k  ) );

         h_eta        = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j-2,k+1)
                                      +zp(i+1,j,k+1) - zp(i+1,j-2,k+1) ) * dyinv;

         h_zeta       = Real(0.25) * ( zp(i  ,j,k+2) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+2) - zp(i+1,j  ,k  ) );
         pz_md_lo_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = Real(0.5) * ( x(i,j,k  ) + x(i,j-1,k  )
                                     -x(i,j,k-1) - x(i,j-1,k-1) );

         h_eta        = Real(0.25) * ( zp(i  ,j,k  ) - zp(i  ,j-2,k  )
                                      +zp(i+1,j,k  ) - zp(i+1,j-2,k  ) ) * dyinv;

         h_zeta       = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On y-face
    py_lo -= Real(0.5) * ( pz_md_lo_hi + pz_md_lo_lo );

    // *********************************************************
    // Hi z-face
    // *********************************************************
    // On z-face
    Real pz_hi    = x(i,j,k+1) - x(i,j,k  );

    // On corners
    Real px_hi_md_hi = Real(0.5) * ( x(i+1,j,k+1) - x(i  ,j  ,k+1)
                                    +x(i+1,j,k  ) - x(i  ,j  ,k  )) * dxinv;
    Real px_lo_md_hi = Real(0.5) * ( x(i  ,j,k+1) - x(i-1,j  ,k+1)
                                    +x(i  ,j,k  ) - x(i-1,j  ,k  )) * dxinv;
    Real py_md_hi_hi = Real(0.5) * ( x(i,j+1,k+1) - x(i  ,j  ,k+1)
                                    +x(i,j+1,k  ) - x(i  ,j  ,k  )) * dyinv;
    Real py_md_lo_hi = Real(0.5) * ( x(i,j  ,k+1) - x(i  ,j-1,k+1)
                                    +x(i,j  ,k  ) - x(i  ,j-1,k  )) * dyinv;

    // On z-face
    Real h_xi_on_zhi  = 0.5 * ( zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1) ) * dxinv;
    Real h_eta_on_zhi = 0.5 * ( zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1) ) * dyinv;

    Real hzeta_inv_on_zhi = 8.0 / ( (zp(i,j,k+2) + zp(i+1,j,k+2) + zp(i,j+1,k+2) + zp(i+1,j+1,k+2))
                                   -(zp(i,j,k  ) + zp(i+1,j,k  ) + zp(i,j+1,k  ) + zp(i+1,j+1,k  )) );
    // pz_hi *= hzeta_inv_on_zhi * (1.0 + h_xi_on_zhi*h_xi_on_zhi + h_eta_on_zhi*h_eta_on_zhi);
    pz_hi *= hzeta_inv_on_zhi;

    pz_hi -= 0.5 * h_xi_on_zhi  * ( (px_hi_md_hi + px_lo_md_hi) - (pz_hi_md_hi + pz_lo_md_hi) );
    pz_hi -= 0.5 * h_eta_on_zhi * ( (py_md_hi_hi + py_md_lo_hi) - (pz_md_hi_hi + pz_md_lo_hi) );

    // *********************************************************
    // Lo z-face
    // *********************************************************
    // On z-face
    Real pz_lo = x(i,j,k  ) - x(i,j,k-1);

    // On corners
    Real px_hi_md_lo = Real(0.5) * ( x(i+1,j  ,k  ) - x(i  ,j  ,k  )
                                    +x(i+1,j  ,k-1) - x(i  ,j  ,k-1)) * dxinv;
    Real px_lo_md_lo = Real(0.5) * ( x(i  ,j  ,k  ) - x(i-1,j  ,k  )
                                    +x(i  ,j  ,k-1) - x(i-1,j  ,k-1)) * dxinv;
    Real py_md_hi_lo = Real(0.5) * ( x(i  ,j+1,k  ) - x(i  ,j  ,k  )
                                    +x(i  ,j+1,k-1) - x(i  ,j  ,k-1)) * dyinv;
    Real py_md_lo_lo = Real(0.5) * ( x(i  ,j  ,k  ) - x(i  ,j-1,k  )
                                    +x(i  ,j  ,k-1) - x(i  ,j-1,k-1)) * dyinv;

    // On z-face
    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k)) * dxinv;
    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k)) * dyinv;

    Real hzeta_inv_on_zlo = 8.0 / ( (zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1))
                                   -(zp(i,j,k-1) + zp(i+1,j,k-1) + zp(i,j+1,k-1) + zp(i+1,j+1,k-1)) );
    // pz_lo *= hzeta_inv_on_zlo * (1.0 + h_xi_on_zlo*h_xi_on_zlo + h_eta_on_zlo*h_eta_on_zlo);
    pz_lo *= hzeta_inv_on_zlo;

    pz_lo -= 0.5 * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= 0.5 * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );

#else
    //
    // This option uses calls to the flux routines so there is
    //      some duplicated computation
    // This option should give the same answer as above
    //
    Real px_lo = -terrpoisson_flux_x(i  ,j,k,x,zp,dxinv);
    Real px_hi = -terrpoisson_flux_x(i+1,j,k,x,zp,dxinv);
    Real py_lo = -terrpoisson_flux_y(i,j  ,k,x,zp,dxinv);
    Real py_hi = -terrpoisson_flux_y(i,j+1,k,x,zp,dyinv);
    Real pz_lo = -terrpoisson_flux_z(i,j,k  ,x,zp,dxinv,dyinv);
    Real pz_hi = -terrpoisson_flux_z(i,j,k+1,x,zp,dxinv,dyinv);
#endif

    // *********************************************************
    // Adotx
    // *********************************************************
    Real invdJ = 4.0 / ( zp(i,j,k+1) + zp(i+1,j,k+1) + zp(i,j+1,k+1) + zp(i+1,j+1,k+1)
                        -zp(i,j,k  ) - zp(i+1,j,k  ) - zp(i,j+1,k  ) - zp(i+1,j+1,k  ) );

    Real ax_lo = .5 * (zp(i  ,j,k+1) + zp(i  ,j+1,k+1) - zp(i  ,j,k) - zp(i  ,j+1,k));
    Real ax_hi = .5 * (zp(i+1,j,k+1) + zp(i+1,j+1,k+1) - zp(i+1,j,k) - zp(i+1,j+1,k));
    Real ay_lo = .5 * (zp(i,j  ,k+1) + zp(i+1,j  ,k+1) - zp(i,j  ,k) - zp(i+1,j  ,k));
    Real ay_hi = .5 * (zp(i,j+1,k+1) + zp(i+1,j+1,k+1) - zp(i,j+1,k) - zp(i+1,j+1,k));

    y(i,j,k) =  ( (ax_hi*px_hi - ax_lo*px_lo) * dxinv
                 +(ay_hi*py_hi - ay_lo*py_lo) * dyinv
                 +(      pz_hi -       pz_lo)        ) * invdJ;
}
#endif
