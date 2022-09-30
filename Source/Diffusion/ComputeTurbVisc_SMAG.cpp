#include <Diffusion.H>

using namespace amrex;

void ComputeTurbVisc_SMAG(Box& bxcc,
                          const Array4<Real>& K_turb,
                          const Array4<const Real>& cell_data,
                          const Array4<const Real>& tau11,
                          const Array4<const Real>& tau22,
                          const Array4<const Real>& tau33,
                          const Array4<const Real>& tau12,
                          const Array4<const Real>& tau13,
                          const Array4<const Real>& tau23,
                          const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                          const SolverChoice& solverChoice)
{
  const Real cellVol = 1.0 / (dxInv[0] * dxInv[1] * dxInv[2]);
  const Real Delta = std::pow(cellVol,1.0/3.0);
  Real Cs = solverChoice.Cs;
  Real CsDeltaSqr = Cs*Cs*Delta*Delta;
  ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real s11bar = tau11(i,j,k);
    Real s22bar = tau22(i,j,k);
    Real s33bar = tau33(i,j,k);
    Real s12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                         + tau12(i+1, j  , k  ) + tau12(i+1, j+1, k  ) );
    Real s13bar = 0.25 * ( tau13(i  , j  , k  ) + tau13(i  , j  , k+1)
                         + tau13(i+1, j  , k  ) + tau13(i+1, j  , k+1) );
    Real s23bar = 0.25 * ( tau23(i  , j  , k  ) + tau23(i  , j  , k+1)
                         + tau23(i  , j+1, k  ) + tau23(i  , j+1, k+1) );
    Real SmnSmn = s11bar*s11bar + s22bar*s22bar + s33bar*s33bar
                + 2.0*s12bar*s12bar + 2.0*s13bar*s13bar + 2.0*s23bar*s23bar;
    
    K_turb(i, j, k, EddyDiff::Mom_h) = 2.0 * CsDeltaSqr * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
    K_turb(i, j, k, EddyDiff::Mom_v) = K_turb(i, j, k, EddyDiff::Mom_h);
  });

  Real inv_Pr_t    = solverChoice.Pr_t_inv;
  Real inv_Sc_t    = solverChoice.Sc_t_inv;
  Real inv_sigma_k = 1.0 / solverChoice.sigma_k;
  Vector<Real> Factors = {inv_Pr_t, inv_Sc_t, inv_sigma_k, inv_sigma_k};
  Real* fac_ptr = Factors.data();
  
  bool use_KE  = (solverChoice.les_type == LESType::Deardorff);
  bool use_QKE = (solverChoice.use_QKE && solverChoice.diffuse_QKE_3D);

  int offset = EddyDiff::Theta_h;
  int ntot   = 4;
  
  if(use_QKE) {
    int ncomp  = 4;    
    ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      int indx = n + offset;
      int indx_v = indx + ntot;
      K_turb(i,j,k,indx)   = 0.5 * K_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
      K_turb(i,j,k,indx_v) = K_turb(i,j,k,indx);
    });
  } else if (use_KE) {
    int ncomp  = 3;
    ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      int indx   = n + offset;
      int indx_v = indx + ntot;
      K_turb(i,j,k,indx)   = 0.5 * K_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
      K_turb(i,j,k,indx_v) = K_turb(i,j,k,indx);
    });
  } else {
    int ncomp  = 2;
    ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      int indx   = n + offset;
      int indx_v = indx + ntot;
      K_turb(i,j,k,indx)   = 0.5 * K_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
      K_turb(i,j,k,indx_v) = K_turb(i,j,k,indx);
    });
  }
}
