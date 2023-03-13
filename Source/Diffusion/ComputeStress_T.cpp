#include <Diffusion.H>
#include <TerrainMetrics.H>

using namespace amrex;

void
ComputeStressConsVisc_T(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau13,
                        Array4<Real>& tau21, Array4<Real>& tau23,
                        Array4<Real>& tau31, Array4<Real>& tau32,
                        const Array4<const Real>& er_arr,
                        const Array4<const Real>& z_nd  ,
                        const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D)
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;

    // We don't need x/y ghost cells for tau33 (avoids linear comb issues)
    bxcc2.grow(IntVect(-1,-1,0));

    // First block: compute S-D
    //***********************************************************************************
    Real OneThird   = (1./3.);
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau11(i,j,k) -= OneThird*er_arr(i,j,k);
        tau22(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau33 first (no linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    amrex::ParallelFor(bxcc2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);

        Real tau31bar = 0.25 * ( tau31(i  , j  , k  ) + tau31(i+1, j  , k  )
                               + tau31(i  , j  , k+1) + tau31(i+1, j  , k+1) );
        Real tau32bar = 0.25 * ( tau32(i  , j  , k  ) + tau32(i  , j+1, k  )
                               + tau32(i  , j  , k+1) + tau32(i  , j+1, k+1) );
        Real mu_tot   = mu_eff;

        tau33(i,j,k) -= met_h_xi*tau31bar + met_h_eta*tau32bar;
        tau33(i,j,k) *= mu_tot;
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  ) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k+1) + tau11(i-1, j  , k+1) );
            Real tau11bar = 1.5*tau11lo - 0.5*tau11hi;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  ) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k+1) + tau12(i  , j+1, k+1) );
            Real tau12bar = 1.5*tau12lo - 0.5*tau12hi;

            Real mu_tot = mu_eff;

            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
            tau13(i,j,k) *= mu_tot;

            tau31(i,j,k) *= mu_tot*met_h_zeta;
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  ) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k+1) + tau21(i+1, j  , k+1) );
            Real tau21bar = 1.5*tau21lo - 0.5*tau21hi;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  ) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k+1) + tau22(i  , j-1, k+1) );
            Real tau22bar = 1.5*tau22lo - 0.5*tau22hi;

            Real mu_tot   = mu_eff;

            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
            tau23(i,j,k) *= mu_tot;

            tau32(i,j,k) *= mu_tot*met_h_zeta;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k-2) + tau11(i-1, j  , k-2) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
            Real tau11bar = 1.5*tau11hi - 0.5*tau11lo;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k-2) + tau12(i  , j+1, k-2) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
            Real tau12bar = 1.5*tau12hi - 0.5*tau12lo;

            Real mu_tot   = mu_eff;

            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
            tau13(i,j,k) *= mu_tot;

            tau31(i,j,k) *= mu_tot*met_h_zeta;
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k-2) + tau21(i+1, j  , k-2) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
            Real tau21bar = 1.5*tau21hi - 0.5*tau21lo;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k-2) + tau22(i  , j-1, k-2) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
            Real tau22bar = 1.5*tau22hi - 0.5*tau22lo;

            Real mu_tot   = mu_eff;

            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
            tau23(i,j,k) *= mu_tot;

            tau32(i,j,k) *= mu_tot*met_h_zeta;
        });
    }

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (valid averaging region)
    //-----------------------------------------------------------------------------------
    amrex::ParallelFor(tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        Real tau11bar = 0.25 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  )
                               + tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
        Real tau12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                               + tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
        Real mu_tot   = mu_eff;

        tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
        tau13(i,j,k) *= mu_tot;

        tau31(i,j,k) *= mu_tot*met_h_zeta;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = 0.25 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  )
                               + tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
        Real tau22bar = 0.25 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  )
                               + tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
        Real mu_tot   = mu_eff;

        tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
        tau23(i,j,k) *= mu_tot;

        tau32(i,j,k) *= mu_tot*met_h_zeta;
    });

    // Fill the remaining components: tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    amrex::ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);
        Real mu_tot = mu_eff;

        tau11(i,j,k) *= mu_tot*met_h_zeta;
        tau22(i,j,k) *= mu_tot*met_h_zeta;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        Real mu_tot = mu_eff;

        tau12(i,j,k) *= mu_tot*met_h_zeta;
        tau21(i,j,k) *= mu_tot*met_h_zeta;
    });
}


void
ComputeStressVarVisc_T(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                       const Array4<const Real>& mu_turb,
                       Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                       Array4<Real>& tau12, Array4<Real>& tau13,
                       Array4<Real>& tau21, Array4<Real>& tau23,
                       Array4<Real>& tau31, Array4<Real>& tau32,
                       const Array4<const Real>& er_arr,
                       const Array4<const Real>& z_nd  ,
                       const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D)
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;

    // First block: compute S-D
    //***********************************************************************************
    Real OneThird   = (1./3.);
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau11(i,j,k) -= OneThird*er_arr(i,j,k);
        tau22(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau33 first (no linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    amrex::ParallelFor(bxcc2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);

        Real tau31bar = 0.25 * ( tau31(i  , j  , k  ) + tau31(i+1, j  , k  )
                               + tau31(i  , j  , k+1) + tau31(i+1, j  , k+1) );
        Real tau32bar = 0.25 * ( tau32(i  , j  , k  ) + tau32(i  , j+1, k  )
                               + tau32(i  , j  , k+1) + tau32(i  , j+1, k+1) );

        Real mu_tot   = mu_eff + 2.0*mu_turb(i, j, k, EddyDiff::Mom_v);

        tau33(i,j,k) -= met_h_xi*tau31bar + met_h_eta*tau32bar;
        tau33(i,j,k) *= mu_tot;
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  ) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k+1) + tau11(i-1, j  , k+1) );
            Real tau11bar = 1.5*tau11lo - 0.5*tau11hi;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  ) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k+1) + tau12(i  , j+1, k+1) );
            Real tau12bar = 1.5*tau12lo - 0.5*tau12hi;

            Real mu_bar = 0.25*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real mu_tot = mu_eff + 2.0*mu_bar;

            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
            tau13(i,j,k) *= mu_tot;

            tau31(i,j,k) *= mu_tot*met_h_zeta;
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  ) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k+1) + tau21(i+1, j  , k+1) );
            Real tau21bar = 1.5*tau21lo - 0.5*tau21hi;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  ) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k+1) + tau22(i  , j-1, k+1) );
            Real tau22bar = 1.5*tau22lo - 0.5*tau22hi;

            Real mu_bar = 0.25*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real mu_tot = mu_eff + 2.0*mu_bar;

            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
            tau23(i,j,k) *= mu_tot;

            tau32(i,j,k) *= mu_tot*met_h_zeta;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k-2) + tau11(i-1, j  , k-2) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
            Real tau11bar = 1.5*tau11hi - 0.5*tau11lo;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k-2) + tau12(i  , j+1, k-2) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
            Real tau12bar = 1.5*tau12hi - 0.5*tau12lo;

            Real mu_bar = 0.25*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real mu_tot = mu_eff + 2.0*mu_bar;

            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
            tau13(i,j,k) *= mu_tot;

            tau31(i,j,k) *= mu_tot*met_h_zeta;
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k-2) + tau21(i+1, j  , k-2) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
            Real tau21bar = 1.5*tau21hi - 0.5*tau21lo;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k-2) + tau22(i  , j-1, k-2) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
            Real tau22bar = 1.5*tau22hi - 0.5*tau22lo;

            Real mu_bar = 0.25*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real mu_tot = mu_eff + 2.0*mu_bar;

            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
            tau23(i,j,k) *= mu_tot;

            tau32(i,j,k) *= mu_tot*met_h_zeta;
        });
    }

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (valid averaging region)
    //-----------------------------------------------------------------------------------
    amrex::ParallelFor(tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        Real tau11bar = 0.25 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  )
                               + tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
        Real tau12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                               + tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );

        Real mu_bar = 0.25 * ( mu_turb(i-1, j  , k  , EddyDiff::Mom_v) + mu_turb(i  , j  , k  , EddyDiff::Mom_v)
                              + mu_turb(i-1, j  , k-1, EddyDiff::Mom_v) + mu_turb(i  , j  , k-1, EddyDiff::Mom_v) );
        Real mu_tot = mu_eff + 2.0*mu_bar;

        tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
        tau13(i,j,k) *= mu_tot;

        tau31(i,j,k) *= mu_tot*met_h_zeta;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = 0.25 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  )
                               + tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
        Real tau22bar = 0.25 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  )
                               + tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );

        Real mu_bar = 0.25 * ( mu_turb(i  , j-1, k  , EddyDiff::Mom_v) + mu_turb(i  , j  , k  , EddyDiff::Mom_v)
                             + mu_turb(i  , j-1, k-1, EddyDiff::Mom_v) + mu_turb(i  , j  , k-1, EddyDiff::Mom_v) );
        Real mu_tot = mu_eff + 2.0*mu_bar;

        tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
        tau23(i,j,k) *= mu_tot;

        tau32(i,j,k) *= mu_tot*met_h_zeta;
    });

    // Fill the remaining components: tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    amrex::ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);

        Real mu_tot = mu_eff + 2.0*mu_turb(i, j, k, EddyDiff::Mom_h);

        tau11(i,j,k) *= mu_tot*met_h_zeta;
        tau22(i,j,k) *= mu_tot*met_h_zeta;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        Real mu_bar = 0.25*( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                           + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real mu_tot = mu_eff + 2.0*mu_bar;

        tau12(i,j,k) *= mu_tot*met_h_zeta;
        tau21(i,j,k) *= mu_tot*met_h_zeta;
    });
}
