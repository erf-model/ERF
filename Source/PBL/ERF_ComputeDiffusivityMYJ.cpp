#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"

using namespace amrex;

void
ComputeDiffusivityMYJ (const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& cons_in,
                       MultiFab& eddyViscosity,
                       const Geometry& geom,
                       const TurbChoice& /*turbChoice*/,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool /*use_moisture*/,
                       int level,
                       const BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<MultiFab>& z_phys_nd,
                       const MoistureComponentIndices& moisture_indices)
{
    // Dirichlet flags to switch derivative stencil
    bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
    bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir) );
    bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
    bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir) );
    bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
    bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir) );

    // Epsilon
    Real eps = std::numeric_limits<Real>::epsilon();

    // Expose constants
    Real d_kappa   = KAPPA;
    Real d_gravity = CONST_GRAV;

    // Closure coefficients (from Janjic (2002), NCEP Office Note 437)
    Real EPSRS   = 1.0e-7;
    Real EPSRU   = 1.0e-7;
    Real EPSL    = 0.32;
    Real EPSQ2   = 0.2;
    Real FH      = 1.01;

    //Real VKARMAN = KAPPA;
    Real G       = CONST_GRAV;
    //Real ALPH    = 0.30;
    Real BETA    = 1./273.;
    //Real EL0MAX  = 1000.;
    //Real EL0MIN  = 1.;
    Real ELFC    = 0.23*0.5;
    //Real GAM1    = 0.2222222222222222222;
    //Real PRT     = 1.;

    Real A1  =  0.659888514560862645;
    Real A2x =  0.6574209922667784586;
    Real B1  = 11.87799326209552761;
    Real B2  =  7.226971804046074028;
    Real C1  =  0.000830955950095854396;

    //Real A2S =  17.2693882;
    //Real A3S = 273.16;
    //Real A4S =  35.86;

    //Real ELZ0   = 0.;
    //Real ESQ    = 5.0;
    //Real EXCM   = 0.001;
    //Real FHNEU  = 0.8;
    //Real GLKBR  = 10.;
    //Real GLKBS  = 30.;
    //Real QVISC  = 2.1E-5;
    //Real RFC    = 0.191;
    //Real RIC    = 0.505;
    //Real SMALL  = 0.35;
    //Real SQPR   = 0.84;
    //Real SQSC   = 0.84;
    //Real SQVISC = 258.2;
    //Real TVISC  = 2.1E-5;
    //Real USTC   = 0.7;
    //Real USTR   = 0.225;
    //Real VISC   = 1.5E-5;
    //Real WOLD   = 0.15;
    //Real WWST   = 1.2;
    //Real ZTMAX  = 1.;
    //Real ZTFC   = 1.;
    //Real ZTMIN  =-5.;

    Real BTG    = BETA*G;
    //Real CZIV   = SMALL*GLKBS;
    //Real EP_1   = R_v/R_d - 1.0;
    //Real ESQHF  = 0.5*5.0;
    //Real GRRS   = GLKBR/GLKBS;
    //Real RB1    = 1./B1;
    //Real RTVISC = 1./TVISC;
    //Real RVISC  = 1./VISC;
    //Real ZQRZT  = SQSC/SQPR;

    Real ADNH  = 9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG;
    Real ADNM  = 18.*A1*A1*A2x*(B2-3.*A2x)*BTG;
    //Real ANMH  = -9.*A1*A2x*A2x*BTG*BTG;
    //Real ANMM  = -3.*A1*A2x*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)*BTG;
    Real BDNH  = 3.*A2x*(7.*A1+B2)*BTG;
    Real BDNM  = 6.*A1*A1;
    //Real BEQH  = A2x*B1*BTG+3.*A2x*(7.*A1+B2)*BTG;
    //Real BEQM  = -A1*B1*(1.-3.*C1)+6.*A1*A1;
    //Real BNMH  = -A2x*BTG;
    //Real BNMM  = A1*(1.-3.*C1);
    Real BSHH  = 9.*A1*A2x*A2x*BTG;
    Real BSHM  = 18.*A1*A1*A2x*C1;
    Real BSMH  = -3.*A1*A2x*(3.*A2x+3.*B2*C1+12.*A1*C1-B2)*BTG;
    Real CESH  = A2x;
    Real CESM  = A1*(1.-3.*C1);
    //Real CNV   = EP_1*G/BTG;
    //Real ELFCS = VKARMAN*BTG;
    //Real FZQ1  = RTVISC*QVISC*ZQRZT;
    //Real FZQ2  = RTVISC*QVISC*ZQRZT;
    //Real FZT1  = RVISC*TVISC*SQPR;
    //Real FZT2  = CZIV*GRRS*TVISC*SQPR;
    //Real FZU1  = CZIV*VISC;
    //Real PIHF  = 0.5*PI;
    //Real RFAC  = RIC/(FHNEU*RFC*RFC);
    //Real RQVISC = 1./QVISC;
    //Real RRIC   = 1./RIC;
    //Real USTFC  = 0.018/G;
    //Real WNEW   = 1.-WOLD;
    //Real WWST2  = WWST*WWST;

    Real AEQH = 9.*A1*A2x*A2x*B1*BTG*BTG
              + 9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG;
    Real AEQM =  3.*A1*A2x*B1*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)*BTG
              + 18.*A1*A1*A2x*(B2-3.*A2x)*BTG;

    Real REQU  = -AEQH/AEQM;
    Real EPSGH = 1.E-9;
    //Real EPSGM = REQU*EPSGH;

    Real UBRYL = (18.*REQU*A1*A1*A2x*B2*C1*BTG + 9.*A1*A2x*A2x*B2*BTG*BTG)
               / (REQU*ADNM+ADNH);
    Real UBRY  = (1.+EPSRS)*UBRYL;
    Real UBRY3 = 3.*UBRY;

    Real AUBH  = 27.*A1*A2x*A2x*B2*BTG*BTG-ADNH*UBRY3;
    Real AUBM  = 54.*A1*A1*A2x*B2*C1*BTG-ADNM*UBRY3;
    Real BUBH  = (9.*A1*A2x+3.*A2x*B2)*BTG-BDNH*UBRY3;
    Real BUBM  = 18.*A1*A1*C1-BDNM*UBRY3;
    Real CUBR  = 1.-UBRY3;
    Real RCUBR = 1./CUBR;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(eddyViscosity,false); mfi.isValid(); ++mfi) {

        const Box&  bx = mfi.growntilebox(IntVect(1,1,0));
        const Array4<Real const>& cell_data = cons_in.array(mfi);
        const Array4<Real      >& K_turb    = eddyViscosity.array(mfi);
        const Array4<Real const>& uvel      = xvel.array(mfi);
        const Array4<Real const>& vvel      = yvel.array(mfi);

        // Ensure the box spans the vertical domain
        const Box& dbx = geom.Domain();
        AMREX_ALWAYS_ASSERT(bx.smallEnd(2) == dbx.smallEnd(2) && bx.bigEnd(2) == dbx.bigEnd(2));

        // Create a plane box
        int klo = bx.smallEnd(2);
        int khi = bx.bigEnd(2);
        Box planexy = makeSlab(bx,2,klo);

        // Expose for GPU capture
        const GeometryData gdata = geom.data();

        // Allocate space for integrals
        const Box xybx = PerpendicularBox<ZDir>(bx, IntVect{0,0,0});
        FArrayBox qintegral(xybx,2);
        qintegral.setVal<RunOn::Device>(0.0);
        FArrayBox qturb(bx,1); FArrayBox qturb_old(bx,1);
        const Array4<Real> qint = qintegral.array();
        const Array4<Real> qvel = qturb.array();

        // Terrain and gradient calcs
        const Array4<Real const> &z_nd_arr = z_phys_nd->array(mfi);
        const auto& dxInv = geom.InvCellSizeArray();
        int izmin = geom.Domain().smallEnd(2);
        int izmax = geom.Domain().bigEnd(2);

        // Theta surface
        const auto& t_mean_mf = SurfLayer->get_mac_avg(level,4); // theta_v
        const auto& tm_arr    = t_mean_mf->const_array(mfi);

        // Vertical integrals to compute l0
        if (use_terrain_fitted_coords) {
            ParallelFor(planexy, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
            {
                // Locate PBL k index
                int kpbl = klo;
                for (int k(klo); k<=khi; ++k) {
                    if (qvel(i,j,k)*qvel(i,j,k)<=EPSQ2*FH) {
                        kpbl = k;
                        continue;
                    }
                }

                for (int k(klo); k<=kpbl; ++k) {
                    // q^2 / 2 is the TKE
                    qvel(i,j,k) = std::sqrt(2.0 * cell_data(i,j,k,RhoKE_comp) / cell_data(i,j,k,Rho_comp));
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "KE must have a positive value");
                    const Real Zval = Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr);
                    const Real dz   = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd_arr);
                    Gpu::Atomic::Add(&qint(i,j,0,0), Zval*qvel(i,j,k)*dz);
                    Gpu::Atomic::Add(&qint(i,j,0,1),      qvel(i,j,k)*dz);
                }
            });
        } else {
            ParallelFor(planexy, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
            {
                // Locate PBL k index
                int kpbl = klo;
                for (int k(klo); k<=khi; ++k) {
                    if (qvel(i,j,k)*qvel(i,j,k)<=EPSQ2*FH) {
                        kpbl = k;
                        continue;
                    }
                }

                for (int k(klo); k<=kpbl; ++k) {
                    // q^2 / 2 is the TKE
                    qvel(i,j,k) = std::sqrt(2.0 * cell_data(i,j,k,RhoKE_comp) / cell_data(i,j,k,Rho_comp));
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "KE must have a positive value");
                    // Not multiplying by dz: it's constant and would fall out when we divide qint0/qint1 anyway
                    const Real Zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                    Gpu::Atomic::Add(&qint(i,j,0,0), Zval*qvel(i,j,k));
                    Gpu::Atomic::Add(&qint(i,j,0,1),      qvel(i,j,k));
                }
            });
        }

        // Main work to fill diffusivities
        ParallelFor(planexy, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
        {
            // Locate PBL k index
            int kpbl = klo;
            for (int k(klo); k<=khi; ++k) {
                if (qvel(i,j,k)*qvel(i,j,k)<=EPSQ2*FH) {
                    kpbl = k;
                    continue;
                }
            }

            // Compute diffusivities in each column
            for (int k(klo); k<=khi; ++k) {
                // Gradients for shear and buoy production
                const Real met_h_zeta = use_terrain_fitted_coords ? Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd_arr) : 1.0;
                Real dthetavdz, dudz, dvdz;
                ComputeVerticalDerivativesPBL(i, j, k,
                                              uvel, vvel, cell_data, izmin, izmax, dxInv[2]/met_h_zeta,
                                              c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                              u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                              v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                              dthetavdz, dudz, dvdz,
                                              moisture_indices);

                // Calculate dimensional production terms
                Real theta0 = tm_arr(i,j,0);
                Real GML    = dudz*dudz + dvdz*dvdz;
                Real GHL    = -(d_gravity/theta0) * dthetavdz;

                // Find the maximum mixing length
                Real ELM;
                if (GHL >= EPSGH) {
                    if (GML/GHL <= REQU) {
                        ELM = EPSL;
                    } else {
                      Real AUBR = (AUBM*GML+AUBH*GHL)*GHL;
                      Real BUBR = BUBM*GML+BUBH*GHL;
                      Real QOL2ST = (-0.5*BUBR+std::sqrt(BUBR*BUBR*0.25-AUBR*CUBR))*RCUBR;
                      Real ELOQ2X = 1./QOL2ST;
                      ELM = std::max(std::sqrt(ELOQ2X*qvel(i,j,k)*qvel(i,j,k)),EPSL);
                    }
                } else {
                    Real ADEN = (ADNM*GML+ADNH*GHL)*GHL;
                    Real BDEN = BDNM*GML+BDNH*GHL;
                    Real QOL2UN = -0.5*BDEN+std::sqrt(BDEN*BDEN*0.25-ADEN);
                    Real ELOQ2X = 1./(QOL2UN+EPSRU);
                    ELM = std::max(std::sqrt(ELOQ2X*qvel(i,j,k)*qvel(i,j,k)),EPSL);
                }

                // Compute master length scale
                Real L;
                if (k>=kpbl) {
                    L = std::min((met_h_zeta/dxInv[2])*ELFC, ELM);
                } else {
                    Real l0 = 0.23*qint(i,j,0,0)/qint(i,j,0,1);
                    const Real zval = use_terrain_fitted_coords ? Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr)
                                                                : gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                    L = std::min(l0*d_kappa*zval / (d_kappa*zval + l0), ELM);
                }

                // L^n/Q^n
                Real ELOQ2 = L*L/(qvel(i,j,k)*qvel(i,j,k) + eps);
                Real ELOQ4 = ELOQ2 * ELOQ2;

                // COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
                Real ADEN=(ADNM*GML+ADNH*GHL)*GHL;
                Real BDEN= BDNM*GML+BDNH*GHL;
                Real CDEN= 1.;

                // COEFFICIENTS FOR THE SM DETERMINANT
                Real BESM=BSMH*GHL;

                // COEFFICIENTS FOR THE SH DETERMINANT
                Real BESH=BSHM*GML+BSHH*GHL;

                // 1./DENOMINATOR
                Real RDEN=1./(ADEN*ELOQ4+BDEN*ELOQ2+CDEN);

                // SM, SH, SQ
                Real SM=(BESM*ELOQ2+CESM)*RDEN;
                Real SH=(BESH*ELOQ2+CESH)*RDEN;
                Real SQ=3.0*SH;

                // Finally, compute the eddy viscosity/diffusivities
                const Real rho = cell_data(i,j,k,Rho_comp);
                K_turb(i,j,k,EddyDiff::Mom_v  ) = rho * L * qvel(i,j,k) * SM;
                K_turb(i,j,k,EddyDiff::Theta_v) = rho * L * qvel(i,j,k) * SH;
                K_turb(i,j,k,EddyDiff::KE_v   ) = rho * L * qvel(i,j,k) * SQ;
                K_turb(i,j,k,EddyDiff::Q_v    ) = rho * L * qvel(i,j,k) * SH;
                K_turb(i,j,k,EddyDiff::Turb_lengthscale) = L;

                // FOEXTRAP the top and bottom ghost cells
                if (k==klo) {
                    K_turb(i,j,k-1,EddyDiff::Mom_v  ) = K_turb(i,j,k,EddyDiff::Mom_v  );
                    K_turb(i,j,k-1,EddyDiff::Theta_v) = K_turb(i,j,k,EddyDiff::Theta_v);
                    K_turb(i,j,k-1,EddyDiff::KE_v   ) = K_turb(i,j,k,EddyDiff::KE_v   );
                    K_turb(i,j,k-1,EddyDiff::Q_v    ) = K_turb(i,j,k,EddyDiff::Q_v);
                    K_turb(i,j,k-1,EddyDiff::Turb_lengthscale) = K_turb(i,j,k,EddyDiff::Turb_lengthscale);
                } else if (k==khi) {
                    K_turb(i,j,k+1,EddyDiff::Mom_v  ) = K_turb(i,j,k,EddyDiff::Mom_v  );
                    K_turb(i,j,k+1,EddyDiff::Theta_v) = K_turb(i,j,k,EddyDiff::Theta_v);
                    K_turb(i,j,k+1,EddyDiff::KE_v   ) = K_turb(i,j,k,EddyDiff::KE_v   );
                    K_turb(i,j,k+1,EddyDiff::Q_v    ) = K_turb(i,j,k,EddyDiff::Q_v);
                    K_turb(i,j,k+1,EddyDiff::Turb_lengthscale) = K_turb(i,j,k,EddyDiff::Turb_lengthscale);
                }
            } // for k
        });
    }
}
