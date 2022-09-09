#include <SpatialStencils.H>
#include <StressTerm.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
DiffusionSrcForXMomWithTerrain (int i, int j, int k,
                                const amrex::Array4<const amrex::Real>& u,
                                const amrex::Array4<const amrex::Real>& v,
                                const amrex::Array4<const amrex::Real>& w,
                                const amrex::Array4<const amrex::Real>& cons,
                                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                const amrex::Array4<const amrex::Real>& K_turb,
                                const SolverChoice &solverChoice,
                                const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                                const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                                const amrex::Array4<const amrex::Real>& er_arr)
{
    amrex::Real diffContrib = 0.0;

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

        const auto& dom_lo = amrex::lbound(domain);

        amrex::Real met_h_xi,met_h_eta,met_h_zeta;

        amrex::Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
        amrex::Real tau11BarN, tau11BarP, tau12BarN, tau12BarP;
        amrex::Real Tmp11, Tmp12;

        // 11 Next
        // Metric at cell center
        ComputeMetricAtCellCenter(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                  cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau11Next = ComputeStressTermWithTerrain(i+1,j  ,k  , u, v, w, MomentumEqn::x,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, er_arr(i,j,k));
        // Save for average
        Tmp11 = tau11Next;
        // Scale by metric
        tau11Next *= met_h_zeta;
        //-----------------------------------------------------------------------------------
        // 11 Prev
        // Metric at cell center
        ComputeMetricAtCellCenter(i-1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                  cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau11Prev = ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, er_arr(i-1,j,k));
        // Accumulate for average
        Tmp11 += tau11Prev;
        // Scale by metric
        tau11Prev *= met_h_zeta;
        //************************************************************************************
        // 12 Next
        // Metric at EdgeCenterK
        ComputeMetricAtEdgeCenterK(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                     cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau12Next = ComputeStressTermWithTerrain(i  ,j+1,k  , u, v, w, MomentumEqn::x,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        // Save for average
        Tmp12 = tau12Next;
        // Scale by metric
        tau12Next *= met_h_zeta;
        //-----------------------------------------------------------------------------------
        // 12 Prev
        // Metric at EdgeCenterK
        ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau12Prev = ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        // Accumulate for average
        Tmp12 += tau12Prev;
        // Scale by metric
        tau12Prev *= met_h_zeta;
        //************************************************************************************
        // 13 Next
        // Accumulate averages
        tau11BarN  = Tmp11;
        tau12BarN  = Tmp12;
        tau11BarN += ComputeStressTermWithTerrain(i+1,j  ,k+1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, detJ, domain, bc_ptr, er_arr(i  ,j,k+1));
        tau11BarN += ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                  z_nd, detJ, domain, bc_ptr, er_arr(i-1,j,k+1));
        tau11BarN *= 0.25;
        tau12BarN += ComputeStressTermWithTerrain(i  ,j+1,k+1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, detJ, domain, bc_ptr, 0.);
        tau12BarN += ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::x,
                                  DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                  z_nd, detJ, domain, bc_ptr, 0.);
        tau12BarN *= 0.25;

        // Metric at EdgeCenterJ
        ComputeMetricAtEdgeCenterJ(i  ,j  ,k+1,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_xi_eta);
        tau13Next = -met_h_xi * tau11BarN - met_h_eta * tau12BarN
                   + ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::x,
                                       DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, 0.);
        //-----------------------------------------------------------------------------------
        // 13 Prev
        // Accumulate averages
        tau11BarP = Tmp11;
        tau12BarP = Tmp12;
        amrex::Real e1 = (k == dom_lo.z) ? er_arr(i,j,k) : er_arr(i,j,k-1);
        tau11BarP += ComputeStressTermWithTerrain(i+1,j  ,k-1, u, v, w, MomentumEqn::x,
                                   DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                   z_nd, detJ, domain, bc_ptr, e1);

        amrex::Real e2 = (k == dom_lo.z) ? er_arr(i-1,j,k) : er_arr(i-1,j,k-1);
        tau11BarP += ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::x,
                                   DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                   z_nd, detJ, domain, bc_ptr, e2);
        tau11BarP *= 0.25;
        tau12BarP += ComputeStressTermWithTerrain(i  ,j+1,k-1, u, v, w, MomentumEqn::x,
                                   DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                   z_nd, detJ, domain, bc_ptr, 0.);
        tau12BarP += ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::x,
                                   DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                   z_nd, detJ, domain, bc_ptr, 0.);
        tau12BarP *= 0.25;

        // Metric at EdgeCenterJ
        ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_xi_eta);
        tau13Prev = -met_h_xi * tau11BarP - met_h_eta * tau12BarP
                   + ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::x,
                                       DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, 0.);
        //************************************************************************************
        diffContrib = (tau11Next - tau11Prev) * dxInv  // Contribution to x-mom eqn from diffusive flux in x-dir
                    + (tau12Next - tau12Prev) * dyInv  // Contribution to x-mom eqn from diffusive flux in y-dir
                    + (tau13Next - tau13Prev) * dzInv; // Contribution to x-mom eqn from diffusive flux in z-dir

        diffContrib /= 0.5*(detJ(i,j,k) + detJ(i-1,j,k)); // Terrain grid stretching

        if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
        {
          diffContrib *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i-1,j,k,Rho_comp)) / solverChoice.rho0_trans;
        }
    }
    return diffContrib;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
DiffusionSrcForYMomWithTerrain (int i, int j, int k,
                                const amrex::Array4<const amrex::Real>& u,
                                const amrex::Array4<const amrex::Real>& v,
                                const amrex::Array4<const amrex::Real>& w,
                                const amrex::Array4<const amrex::Real>& cons,
                                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                const amrex::Array4<const amrex::Real>& K_turb,
                                const SolverChoice &solverChoice,
                                const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                                const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                                const amrex::Array4<const amrex::Real>& er_arr)
{
    amrex::Real diffContrib = 0.0;

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

        const auto& dom_lo = amrex::lbound(domain);

        amrex::Real met_h_xi,met_h_eta,met_h_zeta;

        amrex::Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
        amrex::Real tau21BarN, tau22BarN, tau21BarP, tau22BarP;
        amrex::Real Tmp21, Tmp22;

        // 21 Next
        ComputeMetricAtEdgeCenterK(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau21Next = ComputeStressTermWithTerrain(i+1,j  ,k  , u, v, w, MomentumEqn::y,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        Tmp21 = tau21Next;
        tau21Next *= met_h_zeta;
        //-----------------------------------------------------------------------------------
        // 21 Prev
        ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau21Prev = ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        Tmp21 += tau21Prev;
        tau21Prev *= met_h_zeta;
        //************************************************************************************
        // 22 Next
        ComputeMetricAtCellCenter(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                  cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau22Next = ComputeStressTermWithTerrain(i  ,j+1,k  , u, v, w, MomentumEqn::y,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, er_arr(i,j,k));
        Tmp22 = tau22Next;
        tau22Next *= met_h_zeta;
        //-----------------------------------------------------------------------------------
        // 22 Prev
        ComputeMetricAtCellCenter(i  ,j-1,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                  cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau22Prev = ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, er_arr(i,j-1,k));
        Tmp22 += tau22Prev;
        tau22Prev *= met_h_zeta;
        //************************************************************************************
        // 23 Next
        // Accumulate averages
        tau21BarN  = Tmp21;
        tau22BarN  = Tmp22;
        tau21BarN += ComputeStressTermWithTerrain(i+1,j  ,k+1, u, v, w, MomentumEqn::y,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        tau21BarN += ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::y,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        tau21BarN *= 0.25;
        tau22BarN += ComputeStressTermWithTerrain(i  ,j+1,k+1, u, v, w, MomentumEqn::y,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, er_arr(i,j  ,k+1));
        tau22BarN += ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::y,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, er_arr(i,j-1,k+1));
        tau22BarN *= 0.25;

        // Metric at EdgeCenterI
        ComputeMetricAtEdgeCenterI(i  ,j  ,k+1,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_xi_eta);
        tau23Next = -met_h_xi * tau21BarN - met_h_eta * tau22BarN
                   + ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::y,
                                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                                  z_nd, detJ, domain, bc_ptr, 0.);
        //-----------------------------------------------------------------------------------
        // 23 Prev
        // Accumulate averages
        tau21BarP = Tmp21;
        tau22BarP = Tmp22;
        tau21BarP += ComputeStressTermWithTerrain(i+1,j  ,k-1, u, v, w, MomentumEqn::y,
                                       DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, 0.);
        tau21BarP += ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::y,
                                       DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, 0.);
        tau21BarP *= 0.25;
        amrex::Real e1 = (k == dom_lo.z) ? er_arr(i,j,k) : er_arr(i,j,k-1);
        tau22BarP += ComputeStressTermWithTerrain(i  ,j+1,k-1, u, v, w, MomentumEqn::y,
                                       DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, e1);

        amrex::Real e2 = (k == dom_lo.z) ? er_arr(i,j-1,k) : er_arr(i,j-1,k-1);
        tau22BarP += ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::y,
                                       DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, e2);
        tau22BarP *= 0.25;

        // Metric at EdgeCenterI
        ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_xi_eta);
        tau23Prev = -met_h_xi * tau21BarP - met_h_eta * tau22BarP
                   + ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::y,
                                       DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, 0.);
        //************************************************************************************
        diffContrib = (tau21Next - tau21Prev) * dxInv  // Contribution to y-mom eqn from diffusive flux in x-dir
                    + (tau22Next - tau22Prev) * dyInv  // Contribution to y-mom eqn from diffusive flux in y-dir
                    + (tau23Next - tau23Prev) * dzInv; // Contribution to y-mom eqn from diffusive flux in z-dir

        diffContrib /= 0.5*(detJ(i,j,k) + detJ(i,j-1,k)); // Terrain grid stretching

        if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
        {
          diffContrib *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j-1,k,Rho_comp)) / solverChoice.rho0_trans;
        }
    }
    return diffContrib;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
DiffusionSrcForZMomWithTerrain (int i, int j, int k,
                                const amrex::Array4<const amrex::Real>& u,
                                const amrex::Array4<const amrex::Real>& v,
                                const amrex::Array4<const amrex::Real>& w,
                                const amrex::Array4<const amrex::Real>& cons,
                                const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                const amrex::Array4<const amrex::Real>& K_turb,
                                const SolverChoice &solverChoice,
                                const amrex::Array4<const amrex::Real>& z_nd, const amrex::Array4<const amrex::Real>& detJ,
                                const amrex::Box& domain, const amrex::BCRec* bc_ptr,
                                const amrex::Array4<const amrex::Real>& er_arr)
{
    amrex::Real diffContrib = 0.0;

    if ( (solverChoice.molec_diff_type != MolecDiffType::None) ||
         (solverChoice.les_type        !=       LESType::None) ||
         (solverChoice.pbl_type        !=       PBLType::None) )
    {
        auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

        amrex::Real met_h_xi,met_h_eta,met_h_zeta;

        // Nodal in k for w-momentum
        int k_extrap_lb = domain.smallEnd(2);
        int k_extrap_ub = domain.bigEnd(2) + 1;

        amrex::Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev, normv;

        amrex::Real tau31BarN, tau32BarN, tau31BarP, tau32BarP;
        amrex::Real Tmp31, Tmp32;

        // 31 Next
        ComputeMetricAtEdgeCenterJ(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau31Next = ComputeStressTermWithTerrain(i+1,j  ,k  , u, v, w, MomentumEqn::z,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        Tmp31 = tau31Next;
        tau31Next *= met_h_zeta;
        //-----------------------------------------------------------------------------------
        //31 Prev
        ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau31Prev = ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::z,
                                      DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        Tmp31 += tau31Prev;
        tau31Prev *= met_h_zeta;
        //************************************************************************************
        // 32 Next
        ComputeMetricAtEdgeCenterI(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau32Next = ComputeStressTermWithTerrain(i  ,j+1,k  , u, v, w, MomentumEqn::z,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        Tmp32 = tau32Next;
        tau32Next *= met_h_zeta;
        //-----------------------------------------------------------------------------------
        // 32 Prev
        ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                   cellSizeInv,z_nd,TerrainMet::h_zeta);
        tau32Prev = ComputeStressTermWithTerrain(i  ,j  ,k  , u, v, w, MomentumEqn::z,
                                      DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                      z_nd, detJ, domain, bc_ptr, 0.);
        Tmp32 += tau32Prev;
        tau32Prev *= met_h_zeta;
        //************************************************************************************
        // 33 Next
        tau31BarN  = Tmp31;
        tau32BarN  = Tmp32;
        if (k==k_extrap_ub) {
          // Extrapolate to upper edge center
          tau31BarN *= 0.75;
          tau31BarN -= 0.25 * ComputeStressTermWithTerrain(i+1,j  ,k-1, u, v, w, MomentumEqn::z,
                                            DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                            z_nd, detJ, domain, bc_ptr, 0.);
          tau31BarN -= 0.25 * ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::z,
                                            DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                            z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarN *= 0.75;
          tau32BarN -= 0.25 * ComputeStressTermWithTerrain(i  ,j+1,k-1, u, v, w, MomentumEqn::z,
                                            DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                            z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarN -= 0.25 * ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::z,
                                            DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                            z_nd, detJ, domain, bc_ptr, 0.);
        } else {
          // Accumulate averages
          tau31BarN += ComputeStressTermWithTerrain(i+1,j  ,k+1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau31BarN += ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau31BarN *= 0.25;
          tau32BarN += ComputeStressTermWithTerrain(i  ,j+1,k+1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarN += ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarN *= 0.25;
        }

        // Metric at cell center
        ComputeMetricAtCellCenter(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,
                                  cellSizeInv,z_nd,TerrainMet::h_xi_eta);

        tau33Next = -met_h_xi * tau31BarN - met_h_eta * tau32BarN
                   + ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::z,
                                       DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                       z_nd, detJ, domain, bc_ptr, er_arr(i,j,k));
        //-----------------------------------------------------------------------------------
        // 33 Prev
        tau31BarP  = Tmp31;
        tau32BarP  = Tmp32;
        if (k==k_extrap_lb) {
          // Extrapolate to lower edge center
          tau31BarP *= 0.75;
          tau31BarP -= 0.25 * ComputeStressTermWithTerrain(i+1,j  ,k+1, u, v, w, MomentumEqn::z,
                                                DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                                z_nd, detJ, domain, bc_ptr, 0.);
          tau31BarP -= 0.25 * ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::z,
                                                DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                                z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarP *= 0.75;
          tau32BarP -= 0.25 * ComputeStressTermWithTerrain(i  ,j+1,k+1, u, v, w, MomentumEqn::z,
                                                DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                                z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarP -= 0.25 * ComputeStressTermWithTerrain(i  ,j  ,k+1, u, v, w, MomentumEqn::z,
                                                DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                                z_nd, detJ, domain, bc_ptr, 0.);
        } else {
          // Accumulate average to tau31BarN
          tau31BarP += ComputeStressTermWithTerrain(i+1,j  ,k-1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau31BarP += ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::x, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau31BarP *= 0.25;
          tau32BarP += ComputeStressTermWithTerrain(i  ,j+1,k-1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarP += ComputeStressTermWithTerrain(i  ,j  ,k-1, u, v, w, MomentumEqn::z,
                                         DiffusionDir::y, cellSizeInv, K_turb, solverChoice,
                                         z_nd, detJ, domain, bc_ptr, 0.);
          tau32BarP *= 0.25;
        }

        // Metrics at cell center
        ComputeMetricAtCellCenter(i  ,j  ,k-1,met_h_xi,met_h_eta,met_h_zeta,
                                  cellSizeInv,z_nd,TerrainMet::h_xi_eta);

        tau33Prev = -met_h_xi * tau31BarP - met_h_eta * tau32BarP
                   + ComputeStressTermWithTerrain(i, j, k, u, v, w, MomentumEqn::z,
                                                  DiffusionDir::z, cellSizeInv, K_turb, solverChoice,
                                                  z_nd, detJ, domain, bc_ptr, er_arr(i,j,k-1));

        //************************************************************************************

        diffContrib = (tau31Next - tau31Prev) * dxInv  // Contribution to z-mom eqn from diffusive flux in x-dir
                    + (tau32Next - tau32Prev) * dyInv  // Contribution to z-mom eqn from diffusive flux in y-dir
                    + (tau33Next - tau33Prev) * dzInv; // Contribution to z-mom eqn from diffusive flux in z-dir

        normv = (k == 0) ? detJ(i,j,k) : 0.5*( detJ(i,j,k) + detJ(i,j,k-1) ); // Terrain grid stretching
        diffContrib /= normv;

        if (solverChoice.molec_diff_type == MolecDiffType::ConstantAlpha)
        {
          diffContrib *= 0.5 * (cons(i,j,k,Rho_comp) + cons(i,j,k-1)) / solverChoice.rho0_trans;
        }
    }
    return diffContrib;
}

