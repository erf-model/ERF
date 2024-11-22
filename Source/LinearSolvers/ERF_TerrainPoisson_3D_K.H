#ifndef ERF_TERRAINPOISSON_3D_K_H_
#define ERF_TERRAINPOISSON_3D_K_H_

#include <AMReX_FArrayBox.H>

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_adotx (int i, int j, int k, amrex::Array4<T> const& y,
                        amrex::Array4<T const> const& x,
                        amrex::Array4<T const> const& zp,
                        T dxinv, T invdy, T dzinv) noexcept
{
    using amrex::Real;
#if 1
    Real h_xi, h_eta, h_zeta;

    // *********************************************************
    // Hi x-face
    // *********************************************************
    // On x-face
    Real px_hi    = x(i+1,j,k) - x(i  ,j,k);

    // On y-edges
    Real pz_hi_md_hi  = Real(0.25) * ( x(i+1,j,k+1) + x(i,j,k+1)
                                      -x(i+1,j,k  ) - x(i,j,k  ) );
         h_xi         = Real(0.25) * ( zp(i+1,j  ,k+1) - zp(i-1,j  ,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i-1,j+1,k+1) );
         h_zeta       = Real(0.25) * ( zp(i+1,j  ,k+2) - zp(i+1,j  ,k  )
                                      +zp(i+1,j+1,k+2) - zp(i+1,j+1,k  ) );
         pz_hi_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_hi_md_lo  = Real(0.5) * ( x(i+1,j,k  ) + x(i,j,k  )
                                     -x(i+1,j,k-1) - x(i,j,k-1) );
         h_xi         = Real(0.25) * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) );
         h_zeta       = Real(0.25) * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On x-face
    px_hi -= Real(0.5) * ( pz_hi_md_hi + pz_hi_md_lo );
    px_hi *= Real(0.5) * ( zp(i+1,j,k+1) + zp(i+1,j+1,k+1) - zp(i+1,j,k) - zp(i+1,j+1,k) ) * dzinv;

    // *********************************************************
    // Lo x-face
    // *********************************************************
    // On x-face
    Real px_lo = x(i  ,j,k) - x(i-1,j,k);

    // On y-edges
    Real pz_lo_md_hi  = Real(0.5) * ( x(i,j,k+1) + x(i-1,j,k+1)
                                     -x(i,j,k  ) - x(i-1,j,k  ) );
         h_xi         = Real(0.25) * ( zp(i,j  ,k+1) - zp(i-2,j  ,k+1)
                                      +zp(i,j+1,k+1) - zp(i-2,j+1,k+1) );
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+2) - zp(i  ,j  ,k  )
                                      +zp(i,j+1,k+2) - zp(i  ,j+1,k  ) );
         pz_lo_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = Real(0.5) * ( x(i,j,k  ) + x(i-1,j,k  )
                                     -x(i,j,k-1) - x(i-1,j,k-1) );

         h_xi         = Real(0.25) * ( zp(i,j  ,k  ) - zp(i-2,j  ,k  )
                                      +zp(i,j+1,k  ) - zp(i-2,j+1,k  ) );
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-face
    px_lo -= Real(0.5) * ( pz_lo_md_hi + pz_lo_md_lo );
    px_lo *= Real(0.5) * ( zp(i,j,k+1) + zp(i,j+1,k+1) - zp(i,j,k) - zp(i,j+1,k) ) * dzinv;

    // *********************************************************
    // Hi y-face
    // *********************************************************
    // On y-face
    Real py_hi    = x(i,j+1,k) - x(i,j,k);

    // On x-edges
    Real pz_md_hi_hi  = Real(0.25) * ( x(i,j+1,k+1) + x(i,j,k+1)
                                      -x(i,j+1,k  ) - x(i,j,k  ) );
         h_eta        = Real(0.25) * ( zp(i  ,j+1,k+1) - zp(i  ,j-1,k+1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j-1,k+1) );
         h_zeta       = Real(0.25) * ( zp(i  ,j+1,k+2) - zp(i  ,j+1,k  )
                                      +zp(i+1,j+1,k+2) - zp(i+1,j+1,k  ) );
         pz_md_hi_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_hi_lo  = Real(0.5) * ( x(i,j+1,k  ) + x(i,j,k  )
                                     -x(i,j+1,k-1) - x(i,j,k-1) );
         h_eta        = Real(0.25) * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k)
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k) );
         h_zeta       = Real(0.25) * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On y-face
    py_hi -= Real(0.5) * ( pz_md_hi_hi + pz_md_hi_lo );
    py_hi *= Real(0.5) * ( zp(i+1,j,k+1) + zp(i+1,j+1,k+1) - zp(i+1,j,k) - zp(i+1,j+1,k) ) * dzinv;

    // *********************************************************
    // Lo y-face
    // *********************************************************
    // On y-face
    Real py_lo = x(i,j,k) - x(i,j-1,k);

    // On x-edges
    Real pz_md_lo_hi  = Real(0.5) * ( x(i,j,k+1) + x(i,j-1,k+1)
                                     -x(i,j,k  ) - x(i,j-1,k  ) );
         h_eta        = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j-2,k+1)
                                      +zp(i+1,j,k+1) - zp(i+1,j-2,k+1) );
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+2) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+2) - zp(i+1,j  ,k  ) );
         pz_md_lo_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = Real(0.5) * ( x(i,j,k  ) + x(i,j-1,k  )
                                     -x(i,j,k-1) - x(i,j-1,k-1) );

         h_eta        = Real(0.25) * ( zp(i  ,j,k  ) - zp(i  ,j-2,k  )
                                      +zp(i+1,j,k  ) - zp(i+1,j-2,k  ) );
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On y-face
    py_lo -= Real(0.5) * ( pz_md_lo_hi + pz_md_lo_lo );
    py_lo *= Real(0.5) * ( zp(i,j,k+1) + zp(i+1,j,k+1) - zp(i,j,k) - zp(i+1,j,k) ) * dzinv;

    // *********************************************************
    // Hi z-face
    // *********************************************************
    // On z-face
    Real pz_hi    = x(i,j,k+1) - x(i,j,k  );

    // On corners
    Real px_hi_md_hi = Real(0.5) * ( x(i+1,j,k+1) - x(i  ,j,k+1)
                                    +x(i+1,j,k  ) - x(i  ,j,k  ));
    Real px_lo_md_hi = Real(0.5) * ( x(i  ,j,k+1) - x(i-1,j,k+1)
                                    +x(i  ,j,k  ) - x(i-1,j,k  ));
    Real py_md_hi_hi = Real(0.5) * ( x(i,j+1,k+1) - x(i,j  ,k+1)
                                    +x(i,j+1,k  ) - x(i,j  ,k  ));
    Real py_md_lo_hi = Real(0.5) * ( x(i,j  ,k+1) - x(i,j-1,k+1)
                                    +x(i,j  ,k  ) - x(i-1,  j,k  ));

    // On z-face
    Real h_xi_on_zhi  = 0.5 * ( zp(i+1,j+1,k+1) + zp(i+1,j,k+1) - zp(i,j+1,k+1) - zp(i,j,k+1) );
    Real h_eta_on_zhi = 0.5 * ( zp(i+1,j+1,k+1) + zp(i,j+1,k+1) - zp(i+1,j,k+1) - zp(i,j,k+1) );
    pz_hi -= 0.5 * h_xi_on_zhi  * ( (px_hi_md_hi + px_lo_md_hi) - (pz_hi_md_hi + pz_lo_md_hi) );
    pz_hi -= 0.5 * h_eta_on_zhi * ( (py_md_hi_hi + py_md_lo_hi) - (pz_md_hi_hi + pz_md_lo_hi) );

    // *********************************************************
    // Lo z-face
    // *********************************************************
    // On z-face
    Real pz_lo = x(i,j,k  ) - x(i,j,k-1);

    // On corners
    Real px_hi_md_lo = Real(0.5) * ( x(i+1,j  ,k  ) - x(i  ,j  ,k  )
                                    +x(i+1,j  ,k-1) - x(i  ,j  ,k-1));
    Real px_lo_md_lo = Real(0.5) * ( x(i  ,j  ,k  ) - x(i-1,j  ,k  )
                                    +x(i  ,j  ,k-1) - x(i-1,j  ,k-1));
    Real py_md_hi_lo = Real(0.5) * ( x(i  ,j+1,k  ) - x(i  ,j  ,k  )
                                    +x(i  ,j+1,k-1) - x(i  ,j  ,k-1));
    Real py_md_lo_lo = Real(0.5) * ( x(i  ,j  ,k  ) - x(i  ,j-1,k  )
                                    +x(i  ,j  ,k-1) - x(i  ,j-1,k-1));

    // On z-face
    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k));
    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k));
    pz_lo -= 0.5 * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= 0.5 * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );
#endif

    // *********************************************************
    // Adotx
    // *********************************************************
    px_hi    = x(i+1,j,k) - x(i  ,j,k);
    px_lo    = x(i  ,j,k) - x(i-1,j,k);
    py_hi    = x(i,j+1,k) - x(i,j  ,k);
    py_lo    = x(i,j  ,k) - x(i,j-1,k);
    pz_hi    = x(i,j,k+1) - x(i,j,k  );
    pz_lo    = x(i,j,k  ) - x(i,j,k-1);

    y(i,j,k) = dxinv * dxinv * (px_hi - px_lo)
        +      invdy * invdy * (py_hi - py_lo)
        +      dzinv * dzinv * (pz_hi - pz_lo);
    amrex::Print() << "Y IS " << amrex::IntVect(i,j,k) << "  " << y(i,j,k) << std::endl;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_flux_x (int i, int j, int k, amrex::Array4<T> const& fx,
                         amrex::Array4<T const> const& sol,
                         amrex::Array4<T const> const& zp,
                         T dxinv, T dzinv) noexcept
{
    using amrex::Real;

    Real h_xi, h_zeta;

    // On x-face
    Real px_lo = (sol(i,j,k) - sol(i-1,j,k)) * dxinv;

    // On y-edges
    Real pz_lo_md_hi  = Real(0.5) * ( sol(i,j,k+1) + sol(i-1,j,k+1)
                                     -sol(i,j,k  ) - sol(i-1,j,k  ) );
         h_xi         = Real(0.25) * ( zp(i,j  ,k+1) - zp(i-2,j  ,k+1)
                                      +zp(i,j+1,k+1) - zp(i-2,j+1,k+1) );
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+2) - zp(i  ,j  ,k  )
                                      +zp(i,j+1,k+2) - zp(i  ,j+1,k  ) );
         pz_lo_md_hi *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i-1,j,k  )
                                     -sol(i,j,k-1) - sol(i-1,j,k-1) );

         h_xi         = Real(0.25) * ( zp(i,j  ,k  ) - zp(i-2,j  ,k  )
                                      +zp(i,j+1,k  ) - zp(i-2,j+1,k  ) );
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-face
    px_lo -= Real(0.5) * ( pz_lo_md_hi + pz_lo_md_lo );
    px_lo *= Real(0.5) * ( zp(i,j,k+1) + zp(i,j+1,k+1) - zp(i,j,k) - zp(i,j+1,k) ) * dzinv;

    // On y-edges
    fx(i,j,k) = -px_lo;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_flux_y (int i, int j, int k, amrex::Array4<T> const& fy,
                         amrex::Array4<T const> const& sol,
                         amrex::Array4<T const> const& zp,
                         T dyinv, T dzinv) noexcept
{
    using amrex::Real;

    Real h_eta, h_zeta;

    Real py_lo = (sol(i,j,k) - sol(i,j-1,k)) * dyinv;

    // On x-edges
    Real pz_md_lo_hi  = Real(0.5) * ( sol(i,j,k+1) + sol(i,j-1,k+1)
                                     -sol(i,j,k  ) - sol(i,j-1,k  ) );
         h_eta        = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j-2,k+1)
                                      +zp(i+1,j,k+1) - zp(i+1,j-2,k+1) );
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+2) - zp(i  ,j  ,k  )
                                      +zp(i+1,j,k+2) - zp(i+1,j  ,k  ) );
         pz_md_lo_hi *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i,j-1,k  )
                                     -sol(i,j,k-1) - sol(i,j-1,k-1) );

         h_eta        = Real(0.25) * ( zp(i  ,j,k  ) - zp(i  ,j-2,k  )
                                      +zp(i+1,j,k  ) - zp(i+1,j-2,k  ) );
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On y-face
    py_lo -= Real(0.5) * ( pz_md_lo_hi + pz_md_lo_lo );
    py_lo *= Real(0.5) * ( zp(i,j,k+1) + zp(i+1,j,k+1) - zp(i,j,k) - zp(i+1,j,k) ) * dzinv;

    fy(i,j,k) = -py_lo;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_flux_z (int i, int j, int k, amrex::Array4<T> const& fz,
                         amrex::Array4<T const> const& sol,
                         amrex::Array4<T const> const& zp,
                         T dxinv, T dyinv, T dzinv) noexcept
{
    using amrex::Real;

    Real h_xi, h_eta, h_zeta;

    // On z-face
    Real pz_lo = (sol(i,j,k  ) - sol(i,j,k-1)) * dzinv;

    // On corners
    Real px_hi_md_lo = Real(0.5) * ( sol(i+1,j  ,k  ) - sol(i  ,j  ,k  )
                                    +sol(i+1,j  ,k-1) - sol(i  ,j  ,k-1));
    Real px_lo_md_lo = Real(0.5) * ( sol(i  ,j  ,k  ) - sol(i-1,j  ,k  )
                                    +sol(i  ,j  ,k-1) - sol(i-1,j  ,k-1));
    Real py_md_hi_lo = Real(0.5) * ( sol(i  ,j+1,k  ) - sol(i  ,j  ,k  )
                                    +sol(i  ,j+1,k-1) - sol(i  ,j  ,k-1));
    Real py_md_lo_lo = Real(0.5) * ( sol(i  ,j  ,k  ) - sol(i  ,j-1,k  )
                                    +sol(i  ,j  ,k-1) - sol(i  ,j-1,k-1));

    // On y-edges
    Real pz_hi_md_lo  = Real(0.5) * ( sol(i+1,j,k  ) + sol(i,j,k  )
                                     -sol(i+1,j,k-1) - sol(i,j,k-1) );
         h_xi         = Real(0.25) * ( zp(i+1,j  ,k  ) - zp(i-1,j  ,k  )
                                      +zp(i+1,j+1,k  ) - zp(i-1,j+1,k  ) );
         h_zeta       = Real(0.25) * ( zp(i+1,j  ,k+1) - zp(i+1,j  ,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_hi_md_lo *= h_xi / h_zeta;

    // On y-edges
    Real pz_lo_md_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i-1,j,k  )
                                     -sol(i,j,k-1) - sol(i-1,j,k-1) );
         h_xi         = Real(0.25) * ( zp(i,j  ,k  ) - zp(i-2,j  ,k  )
                                      +zp(i,j+1,k  ) - zp(i-2,j+1,k  ) );
         h_zeta       = Real(0.25) * ( zp(i,j  ,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i,j+1,k+1) - zp(i  ,j+1,k-1) );
         pz_lo_md_lo *= h_xi / h_zeta;

    // On x-edges
    Real pz_md_hi_lo  = Real(0.5) * ( sol(i,j+1,k  ) + sol(i,j,k  )
                                     -sol(i,j+1,k-1) - sol(i,j,k-1) );
         h_eta        = Real(0.25) * ( zp(i  ,j+1,k  ) - zp(i  ,j-1,k)
                                      +zp(i+1,j+1,k  ) - zp(i+1,j-1,k) );
         h_zeta       = Real(0.25) * ( zp(i  ,j+1,k+1) - zp(i  ,j+1,k-1)
                                      +zp(i+1,j+1,k+1) - zp(i+1,j+1,k-1) );
         pz_md_hi_lo *= h_eta/ h_zeta;

    // On x-edges
    Real pz_md_lo_lo  = Real(0.5) * ( sol(i,j,k  ) + sol(i,j-1,k  )
                                     -sol(i,j,k-1) - sol(i,j-1,k-1) );
         h_eta        = Real(0.25) * ( zp(i  ,j,k  ) - zp(i  ,j-2,k  )
                                      +zp(i+1,j,k  ) - zp(i+1,j-2,k  ) );
         h_zeta       = Real(0.25) * ( zp(i  ,j,k+1) - zp(i  ,j  ,k-1)
                                      +zp(i+1,j,k+1) - zp(i+1,j  ,k-1) );
         pz_md_lo_lo *= h_eta/ h_zeta;

    // On z-face
    Real h_xi_on_zlo  = 0.5 * (zp(i+1,j+1,k) + zp(i+1,j,k) - zp(i,j+1,k) - zp(i,j,k)) * dxinv;
    Real h_eta_on_zlo = 0.5 * (zp(i+1,j+1,k) + zp(i,j+1,k) - zp(i+1,j,k) - zp(i,j,k)) * dyinv;

    pz_lo -= 0.5 * h_xi_on_zlo  * ( (px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo) );
    pz_lo -= 0.5 * h_eta_on_zlo * ( (py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo) );

    fz(i,j,k) = -pz_lo;
}

#endif
