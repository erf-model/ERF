#include "ERF_HSEUtils.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Rebalance density and potential temperature columns to satisfy hydrostatic equilibrium.
 *
 * @param[in,out] rho Density field
 * @param[in,out] theta Potential temperature field
 * @param[in] qv Water vapor mixing ratio
 * @param[in] qt Total water mixing ratio
 * @param[in] z_phys Physical height field
 * @param[in] geom Grid geometry
 * @param[in] maintain_Th Whether to maintain the existing potential temperature profile
 * @param[in] use_sfc Whether to use a surface boundary condition for initialization
 */
void
rebalance_columns (MultiFab& rho,
                   MultiFab& theta,
                   const MultiFab& qv,
                   const MultiFab& qt,
                   const MultiFab* z_phys,
                   const Geometry& geom,
                   const bool& maintain_Th,
                   bool use_sfc)
{

#ifdef AMREX_USE_FLOAT
    Real tol  = Real(1.0e-6);
#else
    Real tol  = Real(1.0e-10);
#endif
    Real grav = CONST_GRAV;

    // int ncomp    = cons.nComp();
    int k_dom_lo = geom.Domain().smallEnd(2);
    int k_dom_hi = geom.Domain().bigEnd(2);

    for (MFIter mfi(rho,TileNoZ()); mfi.isValid(); ++mfi) {
        Box bx  = mfi.tilebox();
        int klo = bx.smallEnd(2);
        int khi = bx.bigEnd(2);
        AMREX_ALWAYS_ASSERT((klo == k_dom_lo) && (khi == k_dom_hi));
        bx.makeSlab(2,klo);

        const Array4<      Real>& rho_arr = rho.array(mfi);
        const Array4<      Real>&  th_arr = theta.array(mfi);
        const Array4<const Real>&  qv_arr =    qv.const_array(mfi);
        const Array4<const Real>&  qt_arr =    qt.const_array(mfi);

        const Array4<const Real>&   z_arr = z_phys->const_array(mfi);

        ParallelFor(bx, [=,RdoCp_d=RdoCp] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
        {
            // integrate from surface to domain top
            Real dz, F, C;
            Real rho_tot_hi, rho_tot_lo;
            Real z_lo, z_hi;
            Real R_lo, R_hi;
            Real qv_lo, qv_hi;
            Real qt_lo, qt_hi;
            Real Th_lo, Th_hi;
            Real T_hi;
            Real P_lo, P_hi;

            // Integrate from z=0
            if (use_sfc) {
                z_lo = zero; // corresponding to p_0
                z_hi = Real(0.125) * (z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                     +z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(i,j+1,klo+1) + z_arr(i+1,j+1,klo+1));
                dz = z_hi - z_lo;

                // Establish known constant
                qt_lo = qt_arr(i,j,klo);
                qv_lo = qv_arr(i,j,klo);
                Th_lo = th_arr(i,j,klo);
                P_lo  = p_0;
                R_lo  = getRhogivenThetaPress(Th_lo, P_lo, RdoCp_d, qv_lo);
                rho_tot_lo = R_lo * (one + qt_lo);
                C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

                // Initial guess and residual
                qt_hi = qt_arr(i,j,klo);
                qv_hi = qv_arr(i,j,klo);
                Th_hi = th_arr(i,j,klo);
                P_hi  = p_0;
                T_hi  = getTgivenPandTh(P_hi, Th_hi, RdoCp_d);
                R_hi  = getRhogivenThetaPress(Th_hi, P_hi, RdoCp_d, qv_hi);
                rho_tot_hi = R_hi * (one + qt_hi);
                F = P_hi + myhalf*rho_tot_hi*grav*dz + C;

                // Do iterations
                HSEutils::Newton_Raphson_hse(tol, RdoCp_d, dz,
                                             grav, C, Th_hi, T_hi,
                                             qt_hi, qv_hi,
                                             P_hi, R_hi, F, maintain_Th);

                // Assign data
                rho_arr(i,j,klo) = R_hi;
                if (!maintain_Th) { th_arr(i,j,klo)  = getThgivenTandP(T_hi, P_hi, RdoCp_d); }
                P_lo = P_hi;
                z_lo = z_hi;

            // Use SFC state at first CC
            } else {
                z_lo = Real(0.125) * (z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                     +z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(i,j+1,klo+1) + z_arr(i+1,j+1,klo+1));
                P_lo = getPgivenRTh(rho_arr(i,j,klo)*th_arr(i,j,klo),qv_arr(i,j,klo));
                P_hi = P_lo;
            }

            for (int k(klo+1); k<=khi; ++k)
            {
              z_hi = Real(0.125) * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                   +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1));
              dz   = z_hi - z_lo;

              // Establish known constant
              qt_lo = qt_arr(i,j,k-1);
              qv_lo = qv_arr(i,j,k-1);
              Th_lo = th_arr(i,j,k-1);
              R_lo  = getRhogivenThetaPress(Th_lo, P_lo, RdoCp_d, qv_lo);
              rho_tot_lo = R_lo * (one + qt_lo);
              C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

              // Initial guess and residual
              qt_hi = qt_arr(i,j,k);
              qv_hi = qv_arr(i,j,k);
              Th_hi = th_arr(i,j,k);
              T_hi  = getTgivenPandTh(P_hi, Th_hi, RdoCp_d);
              R_hi  = getRhogivenThetaPress(Th_hi, P_hi, RdoCp_d, qv_hi);
              rho_tot_hi = R_hi * (one + qt_hi);
              F = P_hi + myhalf*rho_tot_hi*grav*dz + C;

              // Do iterations
              HSEutils::Newton_Raphson_hse(tol, RdoCp_d, dz,
                                           grav, C, Th_hi, T_hi,
                                           qt_hi, qv_hi,
                                           P_hi, R_hi, F, maintain_Th);

              // Assign data
              rho_arr(i,j,k) = R_hi;
              if (!maintain_Th) { th_arr(i,j,k)  = getThgivenTandP(T_hi, P_hi, RdoCp_d); }
              P_lo = P_hi;
              z_lo = z_hi;
            }
        });
    } // mfi
} // rebalance_columns
