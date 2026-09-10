#include <set>

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

    const BoxArray& ba = rho.boxArray();

    //
    // This is a bottom-up integration: the value in cell k depends only on cells at or below
    // k.  A box that stops below the top of the domain therefore produces exactly the values
    // the full-height column would have had in the cells it does contain, so we do NOT require
    // khi == geom.Domain().bigEnd(2).  That matters for a refined level whose patch covers only
    // the lower part of the domain, which is the normal way to nest an LES region inside a
    // mesoscale parent.
    //
    // A box stacked on top of another box of this level (the BoxArray is split in z, e.g.
    // amr.max_grid_size below the number of cells in z) must continue the integration of the
    // box below it.  Starting afresh from its own lowest cell would take that cell's density,
    // which has not been rebalanced, as given, and the base state would then depend on where
    // the boxes are split.  The boxes are therefore integrated in bands of equal lowest index,
    // bottom up.  The state each column reaches in a cell (pressure, theta, qv, qt and the
    // cell-centre height) is kept in "below", whose z ghost cells are filled from the bands
    // already integrated before the next band starts.  A column whose cell below the box is
    // not covered by this level (the bottom of the domain, or of a refined patch) starts from
    // its own lowest cell.  With a single band nothing is kept and nothing is filled.
    //
    std::set<int> band_klo;
    for (int ib = 0; ib < static_cast<int>(ba.size()); ++ib) {
        band_klo.insert(ba[ib].smallEnd(2));
    }
    const bool multi_band = (band_klo.size() > 1);

    enum { P_below = 0, Th_below, qv_below, qt_below, z_below, done_below, n_below };
    MultiFab below;
    if (multi_band) {
        below.define(ba, rho.DistributionMap(), n_below, IntVect(0,0,1));
        below.setVal(zero);
    }

    for (const int klo_band : band_klo) {

        if (klo_band != *band_klo.begin()) { below.FillBoundary(); }

        for (MFIter mfi(rho,TileNoZ()); mfi.isValid(); ++mfi) {
            Box bx  = mfi.tilebox();
            int klo = bx.smallEnd(2);
            int khi = bx.bigEnd(2);

            // NOTE: TileNoZ() above guarantees that klo/khi are the *box's* z extent rather
            //       than a tile's, so each column is integrated exactly once.
            if (klo != klo_band) { continue; }

            //
            // The use_sfc seeding marches up from p_0 at z = 0, so it is only meaningful for a
            // column that reaches the ground.  A box whose klo is in the interior can use it
            // only if every cell below it belongs to a box of this level to continue from.
            //
            if (use_sfc && klo > k_dom_lo) {
                Box slab_below = mfi.validbox();
                slab_below.setRange(2, klo-1);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ba.contains(slab_below),
                                                 "rebalance_columns with use_sfc requires every column "
                                                 "to reach the bottom of the domain: the integration is "
                                                 "seeded from p_0 at the surface.");
            }
            bx.makeSlab(2,klo);

            const Array4<      Real>& rho_arr = rho.array(mfi);
            const Array4<      Real>&  th_arr = theta.array(mfi);
            const Array4<const Real>&  qv_arr =    qv.const_array(mfi);
            const Array4<const Real>&  qt_arr =    qt.const_array(mfi);

            const Array4<const Real>&   z_arr = z_phys->const_array(mfi);

            const Array4<Real> below_arr = (multi_band) ? below.array(mfi) : Array4<Real>{};
            const bool can_continue = multi_band && (klo > k_dom_lo);

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
                int  k_start;

                // Continue the column from the cell below this box, in a band already done
                if (can_continue && below_arr(i,j,klo-1,done_below) > zero) {
                    P_lo  = below_arr(i,j,klo-1,P_below);
                    Th_lo = below_arr(i,j,klo-1,Th_below);
                    qv_lo = below_arr(i,j,klo-1,qv_below);
                    qt_lo = below_arr(i,j,klo-1,qt_below);
                    z_lo  = below_arr(i,j,klo-1,z_below);
                    P_hi  = P_lo;
                    k_start = klo;

                } else {
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

                    Th_lo = th_arr(i,j,klo);
                    qv_lo = qv_arr(i,j,klo);
                    qt_lo = qt_arr(i,j,klo);
                    k_start = klo+1;
                }

                for (int k(klo); k<=khi; ++k)
                {
                  if (k >= k_start)
                  {
                    z_hi = Real(0.125) * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                         +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1));
                    dz   = z_hi - z_lo;

                    // Establish known constant (the state in cell k-1)
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
                    P_lo  = P_hi;
                    z_lo  = z_hi;
                    Th_lo = th_arr(i,j,k);
                    qv_lo = qv_hi;
                    qt_lo = qt_hi;
                  }

                  // Keep the state in cell k for a box stacked on this one
                  if (multi_band) {
                      below_arr(i,j,k,P_below)    = P_lo;
                      below_arr(i,j,k,Th_below)   = Th_lo;
                      below_arr(i,j,k,qv_below)   = qv_lo;
                      below_arr(i,j,k,qt_below)   = qt_lo;
                      below_arr(i,j,k,z_below)    = z_lo;
                      below_arr(i,j,k,done_below) = one;
                  }
                }
            });
        } // mfi
    } // band
} // rebalance_columns
