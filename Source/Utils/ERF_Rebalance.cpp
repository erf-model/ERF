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

    for (MFIter mfi(rho,TileNoZ()); mfi.isValid(); ++mfi) {
        Box bx  = mfi.tilebox();
        int klo = bx.smallEnd(2);
        int khi = bx.bigEnd(2);

        //
        // This is a bottom-up integration: the value in cell k depends only on cells at
        // or below k.  A box that stops below the top of the domain is therefore perfectly
        // well defined -- it produces exactly the values the full-height column would have
        // had in the cells it does contain -- so we deliberately do NOT require
        // khi == geom.Domain().bigEnd(2) here.  That matters for a refined level whose
        // patch covers only the lower part of the domain, which is the normal way to nest
        // an LES region inside a mesoscale parent.
        //
        // What the integration does require is a valid starting value in its lowest cell.
        // The use_sfc seeding below marches up from p_0 at z = 0, so it is only meaningful
        // for a box that reaches the ground; a box whose klo is in the interior has no
        // surface to start from.  (With amr.refine_grid_layout_z = 0, which is the ERF
        // default set in main.cpp, grids are never chopped in z and this always holds.)
        //
        // NOTE: TileNoZ() above guarantees that klo/khi are the *box's* z extent rather
        //       than a tile's, so each column is integrated exactly once.
        //
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!use_sfc || (klo == k_dom_lo),
                                         "rebalance_columns with use_sfc requires boxes that "
                                         "reach the bottom of the domain: the integration is "
                                         "seeded from p_0 at the surface. Set erf.max_grid_size_z "
                                         "large enough that the grids are not decomposed in z.");
        bx.makeSlab(2,klo);

        const Array4<      Real>& rho_arr = rho.array(mfi);
        const Array4<      Real>&  th_arr = theta.array(mfi);
        const Array4<const Real>&  qv_arr =    qv.const_array(mfi);
        const Array4<const Real>&  qt_arr =    qt.const_array(mfi);

        const Array4<const Real>&   z_arr = z_phys->const_array(mfi);

        ParallelFor(bx, [=,RdoCp_d=RdoCp] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
        {
            // Integrate upward from the bottom of this box to its top
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
