#ifndef ERF_TERRAINPOISSON_3D_K_H_
#define ERF_TERRAINPOISSON_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_adotx (int i, int j, int k, Array4<T> const& y,
                        Array4<T const> const& x,
                        Array4<T const> const& zp,
                        T odx, T ody, T odz) noexcept
{
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
    px_hi *= Real(0.5) * ( zp(i+1,j,k+1) + zp(i+1,j+1,k+1) - zp(i+1,j,k) - zp(i+1,j+1,k) ) * odz;

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
    px_lo *= Real(0.5) * ( zp(i,j,k+1) + zp(i,j+1,k+1) - zp(i,j,k) - zp(i,j+1,k) ) * odz;

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
    py_hi *= Real(0.5) * ( zp(i+1,j,k+1) + zp(i+1,j+1,k+1) - zp(i+1,j,k) - zp(i+1,j+1,k) ) * odz;

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
    py_lo *= Real(0.5) * ( zp(i,j,k+1) + zp(i+1,j,k+1) - zp(i,j,k) - zp(i+1,j,k) ) * odz;

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

    y(i,j,k) = odx * odx * (px_hi - px_lo)
        +      ody * ody * (py_hi - py_lo)
        +      odz * odz * (pz_hi - pz_lo);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_flux_x (Box const& box, Array4<T> const& fx,
                         Array4<T const> const& sol,
                         Array4<T const> const& zp,
                         T dxinv, T dzinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    Real h_xi, h_zeta;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                // On x-face
                Real px_lo = (sol(i,j,k) - sol(i-1,j,k)) * dxinv;

#if 1
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
#endif
                // On y-edges
                fx(i,j,k) = -px_lo;
            }
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_flux_y (Box const& box, Array4<T> const& fy,
                         Array4<T const> const& sol,
                         Array4<T const> const& zp,
                         T dyinv, T dzinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    Real h_eta, h_zeta;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {

                Real py_lo = (sol(i,j,k) - sol(i,j-1,k)) * dyinv;

#if 1
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
#endif

                fy(i,j,k) = -py_lo;
            }
        }
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrpoisson_flux_z (Box const& box, Array4<T> const& fz,
                         Array4<T const> const& sol,
                         Array4<T const> const& zp,
                         T dxinv, T dyinv, T dzinv) noexcept
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    Real h_xi, h_eta, h_zeta;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                // On z-face
                Real pz_lo = (sol(i,j,k  ) - sol(i,j,k-1)) * dzinv;

#if 1
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

#endif
                fz(i,j,k) = -pz_lo;
            }
        }
    }
}

} // amrex
#endif
