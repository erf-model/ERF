#include "ERF_HSEUtils.H"
#include "ERF_Utils.H"

using namespace amrex;

void
rebalance_columns(MultiFab& cons, const MultiFab* z_phys, const MultiFab& qt, const Geometry& geom)
{

#ifdef AMREX_USE_FLOAT
    Real tol  = Real(1.0e-6);
#else
    Real tol  = Real(1.0e-10);
#endif
    Real grav = CONST_GRAV;

    int ncomp    = cons.nComp();
    int k_dom_lo = geom.Domain().smallEnd(2);
    int k_dom_hi = geom.Domain().bigEnd(2);

    for ( MFIter mfi(cons,TileNoZ()); mfi.isValid(); ++mfi ) {
        Box bx  = mfi.tilebox();
        int klo = bx.smallEnd(2);
        int khi = bx.bigEnd(2);
        AMREX_ALWAYS_ASSERT((klo == k_dom_lo) && (khi == k_dom_hi));
        bx.makeSlab(2,klo);

        const Array4<      Real>& con_arr = cons.array(mfi);
        const Array4<const Real>&  qt_arr = qt.const_array(mfi);
        const Array4<const Real>&   z_arr = z_phys->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            // integrate from surface to domain top
            Real Factor;
            Real dz, F, C;
            Real rho_tot_hi, rho_tot_lo;
            Real z_lo, z_hi;
            Real R_lo, R_hi;
            Real qv_lo, qv_hi;
            Real qt_lo, qt_hi;
            Real Th_lo, Th_hi;
            Real P_lo, P_hi;

            // First integrate from sea level to the height at klo
            {
                // Vertical grid spacing
                z_lo = zero; // corresponding to p_0
                z_hi = Real(0.125) * (z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                     +z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(i,j+1,klo+1) + z_arr(i+1,j+1,klo+1));
                dz = z_hi - z_lo;

                // Establish known constant
                qt_lo =  qt_arr(i,j,klo);
                qv_lo = con_arr(i,j,klo,RhoQ1_comp)    / con_arr(i,j,klo,Rho_comp);
                Th_lo = con_arr(i,j,klo,RhoTheta_comp) / con_arr(i,j,klo,Rho_comp);
                P_lo  = p_0;
                R_lo  = getRhogivenThetaPress(Th_lo, P_lo, R_d/Cp_d, qv_lo);
                rho_tot_lo = R_lo * (one + qt_lo);
                C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

                // Initial guess and residual
                qt_hi = qt_arr(i,j,klo);
                qv_hi = con_arr(i,j,klo,RhoQ1_comp)    / con_arr(i,j,klo,Rho_comp);
                Th_hi = con_arr(i,j,klo,RhoTheta_comp) / con_arr(i,j,klo,Rho_comp);
                P_hi  = p_0;
                R_hi  = getRhogivenThetaPress(Th_hi, P_hi, R_d/Cp_d, qv_hi);
                rho_tot_hi = R_hi * (one + qt_hi);
                F = P_hi + myhalf*rho_tot_hi*grav*dz + C;

                // Do iterations
                HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                             grav, C, Th_hi,
                                             qt_hi, qv_hi,
                                             P_hi, R_hi, F);

                // Assign data
                Factor = R_hi / con_arr(i,j,klo,Rho_comp);
                con_arr(i,j,klo,Rho_comp) = R_hi;
                for (int n(1); n<ncomp; ++n) { con_arr(i,j,klo,n) *= Factor; }
                P_lo = P_hi;
                z_lo = z_hi;
            }

            for (int k(klo+1); k<=khi; ++k) {
                // Vertical grid spacing
              z_hi = Real(0.125) * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                   +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1));
              dz   = z_hi - z_lo;

              // Establish known constant
              qt_lo = qt_arr(i,j,k-1);
              qv_lo = con_arr(i,j,k-1,RhoQ1_comp)    / con_arr(i,j,k-1,Rho_comp);
              Th_lo = con_arr(i,j,k-1,RhoTheta_comp) / con_arr(i,j,k-1,Rho_comp);
              R_lo  = getRhogivenThetaPress(Th_lo, P_lo, R_d/Cp_d, qv_lo);
              rho_tot_lo = R_lo * (one + qt_lo);
              C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

              // Initial guess and residual
              qt_hi = qt_arr(i,j,k);
              qv_hi = con_arr(i,j,k,RhoQ1_comp)    / con_arr(i,j,k,Rho_comp);
              Th_hi = con_arr(i,j,k,RhoTheta_comp) / con_arr(i,j,k,Rho_comp);
              R_hi  = getRhogivenThetaPress(Th_hi, P_hi, R_d/Cp_d, qv_hi);
              rho_tot_hi = R_hi * (one + qt_hi);
              F = P_hi + myhalf*rho_tot_hi*grav*dz + C;

              // Do iterations
              HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                           grav, C, Th_hi,
                                           qt_hi, qv_hi,
                                           P_hi, R_hi, F);

              // Assign data
              Factor = R_hi / con_arr(i,j,k,Rho_comp);
              con_arr(i,j,k,Rho_comp) = R_hi;
              for (int n(1); n<ncomp; ++n) { con_arr(i,j,k,n) *= Factor; }
              P_lo = P_hi;
              z_lo = z_hi;
            }
        });
    } // mfi
} // rebalance_columns
