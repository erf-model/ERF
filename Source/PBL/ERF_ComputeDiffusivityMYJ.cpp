#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"

#include <math.h>

using namespace amrex;

void
ComputeDiffusivityMYJ (Real dt,
                       const MultiFab& xvel,
                       const MultiFab& yvel,
                       MultiFab& cons_in,
                       MultiFab& eddyViscosity,
                       const Geometry& geom,
                       const TurbChoice& /*turbChoice*/,
                       std::unique_ptr<SurfaceLayer>& /*SurfLayer*/,
                       bool use_terrain_fitted_coords,
                       bool /*use_moisture*/,
                       int /*level*/,
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

    // Expose constants
    Real d_kappa = KAPPA;

    // Closure coefficients (from Janjic (2002), NCEP Office Note 437)
    Real EPS1    = 1.0e-12;
    Real EPS2    = 0.;
    Real EPSRS   = 1.0e-7;
    Real EPSRU   = 1.0e-7;
    Real EPSTRB  = 1.0e-24;
    Real EPSL    = 0.32;
    Real EPSQ2   = 0.2;
    Real EPSQ1   = std::sqrt(EPSQ2);
    Real ESQHF   = 2.5;
    Real FH      = 1.01;

    Real G       = CONST_GRAV;
    Real ALPHA   = 0.3;
    Real BETA    = 1./273.;
    Real EL0MAX  = 1000.;
    Real EL0MIN  = 1.;
    Real ELFC    = 0.23*0.5;

    Real A1  =  0.659888514560862645;
    Real A2x =  0.6574209922667784586;
    Real B1  = 11.87799326209552761;
    Real B2  =  7.226971804046074028;
    Real C1  =  0.000830955950095854396;

    Real BTG    = BETA*G;
    Real RB1    = 1./B1;

    Real ADNH  = 9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG;
    Real ADNM  = 18.*A1*A1*A2x*(B2-3.*A2x)*BTG;
    Real ANMH  = -9.*A1*A2x*A2x*BTG*BTG;
    Real ANMM  = -3.*A1*A2x*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)*BTG;
    Real BDNH  = 3.*A2x*(7.*A1+B2)*BTG;
    Real BDNM  = 6.*A1*A1;
    Real BEQH  = A2x*B1*BTG+3.*A2x*(7.*A1+B2)*BTG;
    Real BEQM  = -A1*B1*(1.-3.*C1)+6.*A1*A1;
    Real BNMH  = -A2x*BTG;
    Real BNMM  = A1*(1.-3.*C1);
    Real BSHH  = 9.*A1*A2x*A2x*BTG;
    Real BSHM  = 18.*A1*A1*A2x*C1;
    Real BSMH  = -3.*A1*A2x*(3.*A2x+3.*B2*C1+12.*A1*C1-B2)*BTG;
    Real CESH  = A2x;
    Real CESM  = A1*(1.-3.*C1);

    Real AEQH = 9.*A1*A2x*A2x*B1*BTG*BTG
              + 9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG;
    Real AEQM =  3.*A1*A2x*B1*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)*BTG
              + 18.*A1*A1*A2x*(B2-3.*A2x)*BTG;

    Real REQU  = -AEQH/AEQM;
    Real EPSGH = 1.E-9;
    Real EPSGM = REQU*EPSGH;

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
    for (MFIter mfi(eddyViscosity,false); mfi.isValid(); ++mfi) {

        const Box&  bx = mfi.validbox();
        const Array4<Real      >& cell_data = cons_in.array(mfi);
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
        FArrayBox qturb(bx,1);
        FArrayBox qintegral(xybx,2);
        IArrayBox pbl_k(xybx,1);
        qintegral.setVal<RunOn::Device>(0.0);
        pbl_k.setVal<RunOn::Device>(khi);
        const Array4<Real> qint  = qintegral.array();
        const Array4<Real> qvel  = qturb.array();
        const Array4<int>  k_arr = pbl_k.array();

        // Terrain and gradient calcs
        const Array4<Real const> &z_nd_arr = z_phys_nd->array(mfi);
        const auto& dxInv = geom.InvCellSizeArray();
        int izmin = geom.Domain().smallEnd(2);
        int izmax = geom.Domain().bigEnd(2);

        // Ustar for BC
        //MultiFab* ustar = SurfLayer->get_u_star(level);
        //const Array4<Real const>& ustar_arr = ustar->array(mfi);

        // Vertical integrals to compute l0
        if (use_terrain_fitted_coords) {
            ParallelFor(planexy, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
            {
                // Locate PBL k index and set qvel
                for (int k(klo); k<=khi; ++k) {
                    Real q2 = 2.0 * cell_data(i,j,k,RhoKE_comp) / cell_data(i,j,k,Rho_comp);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(q2 > 0.0, "KE must have a positive value");
                    qvel(i,j,k) = std::sqrt(q2);
                    if (q2<=EPSQ2*FH) {
                        k_arr(i,j,0) = std::min(k,k_arr(i,j,0));
                    }
                }

                // Perform integral over PBL height
                for (int k(klo); k<=k_arr(i,j,0); ++k) {
                    const Real dz   = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd_arr);
                    const Real Zval = Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr);
                    Gpu::Atomic::Add(&qint(i,j,0,0), Zval*qvel(i,j,k)*dz);
                    Gpu::Atomic::Add(&qint(i,j,0,1),      qvel(i,j,k)*dz);
                }
            });
        } else {
            ParallelFor(planexy, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept
            {
                // Locate PBL k index and set qvel
                for (int k(klo); k<=khi; ++k) {
                    Real q2 = 2.0 * cell_data(i,j,k,RhoKE_comp) / cell_data(i,j,k,Rho_comp);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(q2 > 0.0, "KE must have a positive value");
                    qvel(i,j,k) = std::sqrt(q2);
                    if (q2<=EPSQ2*FH) {
                        k_arr(i,j,0) = std::min(k,k_arr(i,j,0));
                    }
                }

                // Perform integral over PBL height
                for (int k(klo); k<=k_arr(i,j,0); ++k) {
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
            // Get the PBL k index
            int kpbl = k_arr(i,j,0);

            // Compute the integral length scale
            Real l0 = std::max(std::min(ALPHA*qint(i,j,0,0)/qint(i,j,0,1),EL0MAX),EL0MIN);

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
                Real GML = std::max(dudz*dudz + dvdz*dvdz, EPSGM);
                // NOTE: model uses BTG = beta*g in coeffs above
                // NOTE: sign convention follows code but theory differs
                Real GHL = dthetavdz;
                if (std::fabs(GHL)<=EPSGH) { GHL=EPSGH; }

                // Find the maximum mixing length
                Real ELM;
                if (GHL >= EPSGH) {
                    if (GML/GHL <= REQU) {
                        ELM = EPSL;
                    } else {
                      Real AUBR   = (AUBM*GML+AUBH*GHL)*GHL;
                      Real BUBR   = BUBM*GML+BUBH*GHL;
                      Real QOL2ST = (-0.5*BUBR+std::sqrt(BUBR*BUBR*0.25-AUBR*CUBR))*RCUBR;
                      Real ELOQ2X = 1./QOL2ST;
                      ELM = std::max(std::sqrt(ELOQ2X*qvel(i,j,k)*qvel(i,j,k)),EPSL);
                    }
                } else {
                    Real ADEN   = (ADNM*GML+ADNH*GHL)*GHL;
                    Real BDEN   = BDNM*GML+BDNH*GHL;
                    Real QOL2UN = -0.5*BDEN+std::sqrt(BDEN*BDEN*0.25-ADEN);
                    Real ELOQ2X = 1./(QOL2UN+EPSRU);
                    ELM = std::max(std::sqrt(ELOQ2X*qvel(i,j,k)*qvel(i,j,k)),EPSL);
                }

                // Compute master length scale
                Real L;
                if (k>kpbl) {
                    L = std::min((met_h_zeta/dxInv[2])*ELFC, ELM);
                } else {
                    const Real zval = use_terrain_fitted_coords ? Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr)
                                                                : gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                    L = std::min(l0*d_kappa*zval / (d_kappa*zval + l0), ELM);
                }

                // Update qvel from production and dissipation
                Real AEQU  = (AEQM*GML+AEQH*GHL)*GHL;
                Real BEQU  = BEQM*GML+BEQH*GHL;

                Real EQOL2 = -0.5*BEQU+std::sqrt(BEQU*BEQU*0.25-AEQU);

                if ( ((GML+GHL*GHL)<=EPSTRB) ||
                     ((GHL>=EPSGH) && ((GML/GHL)<=REQU)) ||
                     (EQOL2<=EPS2) ) {
                         L = EPSL;
                         qvel(i,j,k) = EPSQ1;
                } else {
                    Real ANUM=(ANMM*GML+ANMH*GHL)*GHL;
                    Real BNUM= BNMM*GML+BNMH*GHL;

                    Real ADEN=(ADNM*GML+ADNH*GHL)*GHL;
                    Real BDEN= BDNM*GML+BDNH*GHL;
                    Real CDEN= 1.;

                    Real ARHS=-(ANUM*BDEN-BNUM*ADEN)*2.;
                    Real BRHS=- ANUM*4.;
                    Real CRHS=- BNUM*2.;

                    Real DLOQ1=L/qvel(i,j,k);

                    Real ELOQ21=1./EQOL2;
                    Real ELOQ11=std::sqrt(ELOQ21);
                    Real ELOQ31=ELOQ21*ELOQ11;
                    Real ELOQ41=ELOQ21*ELOQ21;
                    Real ELOQ51=ELOQ21*ELOQ31;

                    Real RDEN1=1./(ADEN*ELOQ41+BDEN*ELOQ21+CDEN);

                    Real RHSP1=(ARHS*ELOQ51+BRHS*ELOQ31+CRHS*ELOQ11)*RDEN1*RDEN1;

                    Real DTTURBL=dt;
                    Real ELOQ12=std::max(ELOQ11+(DLOQ1-ELOQ11)*exp(RHSP1*DTTURBL),EPS1);

                    Real ELOQ22=ELOQ12*ELOQ12;
                    Real ELOQ32=ELOQ22*ELOQ12;
                    Real ELOQ42=ELOQ22*ELOQ22;
                    Real ELOQ52=ELOQ22*ELOQ32;

                    Real RDEN2=1./(ADEN*ELOQ42+BDEN*ELOQ22+CDEN);
                    Real RHS2 =-(ANUM*ELOQ42+BNUM*ELOQ22)*RDEN2+RB1;
                    Real RHSP2= (ARHS*ELOQ52+BRHS*ELOQ32+CRHS*ELOQ12)*RDEN2*RDEN2;
                    Real RHST2=RHS2/RHSP2;

                    Real ELOQ13=std::max(ELOQ12-RHST2+(RHST2+DLOQ1-ELOQ12)*exp(RHSP2*DTTURBL),EPS1);

                    Real ELOQN=ELOQ13;
                    if (ELOQN>EPS1) {
                        qvel(i,j,k) = std::max(L/ELOQN,EPSQ1);
                        if (qvel(i,j,k)==EPSQ1) {L = EPSL; }
                    } else {
                        L = EPSL;
                        qvel(i,j,k) = EPSQ1;
                    }
                }
                /*
                // Boundary condition
                if (k==klo) {
                    Real q2 = std::pow(B1,(2./3.))*ustar_arr(i,j,k)*ustar_arr(i,j,k);
                    Real q  = std::max(std::sqrt(q2),EPSQ1);
                    qvel(i,j,k) = 0.5 * (q + qvel(i,j,k));
                }
                */
                cell_data(i,j,k,RhoKE_comp) = 0.5*cell_data(i,j,k,Rho_comp)*qvel(i,j,k)*qvel(i,j,k);

                // L^n/Q^n
                Real ELOQ2 = L*L/(qvel(i,j,k)*qvel(i,j,k));
                Real ELOQ4 = ELOQ2*ELOQ2;

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
                Real SQ=ESQHF*SH;

                // Finally, compute the eddy viscosity/diffusivities
                const Real rho = cell_data(i,j,k,Rho_comp);
                K_turb(i,j,k,EddyDiff::Mom_v  ) = rho * L * qvel(i,j,k) * SM;
                K_turb(i,j,k,EddyDiff::Theta_v) = rho * L * qvel(i,j,k) * SH;
                K_turb(i,j,k,EddyDiff::KE_v   ) = rho * L * qvel(i,j,k) * SQ;
                K_turb(i,j,k,EddyDiff::Q_v    ) = rho * L * qvel(i,j,k) * SH;
                K_turb(i,j,k,EddyDiff::Turb_lengthscale) = L;

                // NOTE: Ghost cells are handled at end of ERF_ComputeTurbulentViscosity.cpp

            } // for k
        }); // ParFor
    } // mfi
}
