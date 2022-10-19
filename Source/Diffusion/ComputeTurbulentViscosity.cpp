/** \file ComputeTurbulentViscosity.cpp */

#include <ABLMost.H>
#include <EddyViscosity.H>
#include <Diffusion.H>

using namespace amrex;

void
ComputeTurbulentViscosityPBL (const amrex::MultiFab& xvel,
                              const amrex::MultiFab& yvel,
                              const amrex::MultiFab& cons_in,
                              amrex::MultiFab& eddyViscosity,
                              const amrex::Geometry& geom,
                              const SolverChoice& solverChoice,
                              std::unique_ptr<ABLMost>& most,
                              bool /*vert_only*/);

/** Compute Eddy Viscosity */
void ComputeTurbulentViscosityLES (const amrex::MultiFab& xvel, const amrex::MultiFab& yvel, const amrex::MultiFab& zvel,
                                   const amrex::MultiFab& cons_in, amrex::MultiFab& eddyViscosity,
                                   const amrex::Geometry& geom,
                                   const SolverChoice& solverChoice,
                                   const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                                   bool vert_only)
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> cellSizeInv = geom.InvCellSizeArray();

    const amrex::Real cellVol = 1.0 / (cellSizeInv[0] * cellSizeInv[1] * cellSizeInv[2]);
    const amrex::Real Delta = std::pow(cellVol,1.0/3.0);

    const auto& domain = geom.Domain();
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    const int klo = dom_lo.z;

    bool l_vert_only = vert_only;
    bool l_use_terrain = solverChoice.use_terrain;

    // SMAGORINSKY: Fill Kturb for momentum in horizontal and vertical
    //***********************************************************************************
    if (solverChoice.les_type == LESType::Smagorinsky)
    {
      const Real cellVol = 1.0 / (dxInv[0] * dxInv[1] * dxInv[2]);
      const Real Delta = std::pow(cellVol,1.0/3.0);
      Real Cs = solverChoice.Cs;
      Real CsDeltaSqr = Cs*Cs*Delta*Delta;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(cons_in,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box bxcc  = mfi.growntilebox(IntVect(1,1,0));

        const Array4<Real>& K_turb = eddyDiffs.array(mfi);
        const amrex::Array4<amrex::Real const > &cell_data = cons_in.array(mfi);
        
        Array4<Real> tau11 = Tau11->array(mfi);
        Array4<Real> tau22 = Tau22->array(mfi);
        Array4<Real> tau33 = Tau33->array(mfi);
        Array4<Real> tau12 = Tau12->array(mfi);
        Array4<Real> tau13 = Tau13->array(mfi);
        Array4<Real> tau23 = Tau23->array(mfi);
        
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
      }     
    }
    // DEARDORFF: Fill Kturb for momentum in horizontal and vertical
    //***********************************************************************************
    else if (solverChoice.les_type == LESType::Deardorff)
    {
      amrex::Real l_C_k = solverChoice.Ck;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for ( amrex::MFIter mfi(cons_in,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box bxcc  = mfi.growntilebox(IntVect(1,1,0));

        const Array4<Real>& K_turb = eddyDiffs.array(mfi);
        const amrex::Array4<amrex::Real const > &cell_data = cons_in.array(mfi);
                
        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          // K = rho * C_k * Delta * KE^(1/2) = C_k * Delta * (rho * RhoKE)^1/2
          K_turb(i,j,k,EddyDiff::Mom_h) = l_C_k * Delta *
            std::sqrt(cell_data(i,j,k,RhoKE_comp) * cell_data(i,j,k,Rho_comp));
          K_turb(i, j, k, EddyDiff::Mom_v) = K_turb(i, j, k, EddyDiff::Mom_h);
        });
      }
    }

    
    // Extrapolate Kturb at top and bottom then fill remaining elements
    //***********************************************************************************
    {
      for ( amrex::MFIter mfi(cons_in,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
        Box bxcc  = mfi.growntilebox(IntVect(1,1,0));
        Box planecc = bxcc; planecc.setBig(2, planecc.smallEnd(2) );
        int k_lo = bxcc.smallEnd(2); int k_hi = bxcc.bigEnd(2);

        // Grow box top and bottom
        bxcc.growLo(2,1); bxcc.growHi(2,1);
      
        amrex::ParallelFor(planecc, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
          K_turb(i, j, k_lo-1, EddyDiff::Mom_h) = K_turb(i, j, k_lo, EddyDiff::Mom_h);
          K_turb(i, j, k_lo-1, EddyDiff::Mom_v) = K_turb(i, j, k_lo, EddyDiff::Mom_v);
          
          K_turb(i, j, k_hi+1, EddyDiff::Mom_h) = K_turb(i, j, k_hi, EddyDiff::Mom_h);
          K_turb(i, j, k_hi+1, EddyDiff::Mom_v) = K_turb(i, j, k_hi, EddyDiff::Mom_v);
        });


        Real inv_Pr_t    = solverChoice.Pr_t_inv;
        Real inv_Sc_t    = solverChoice.Sc_t_inv;
        Real inv_sigma_k = 1.0 / solverChoice.sigma_k;
        Vector<Real> Factors = {inv_Pr_t, inv_Sc_t, inv_sigma_k, inv_sigma_k};
        Gpu::AsyncVector<Real> d_Factors; d_Factors.resize(Factors.size());
        Gpu::copy(Gpu::hostToDevice, Factors.begin(), Factors.end(), d_Factors.begin());
        Real* fac_ptr = d_Factors.data();

        bool use_KE  = (solverChoice.les_type == LESType::Deardorff);
        bool use_QKE = (solverChoice.use_QKE && solverChoice.diffuse_QKE_3D);

        int offset = EddyDiff::Theta_h;
        int ntot   = 5;

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
    }
    

    // Fill in the rest of the turbulent diffusivities
    if (solverChoice.les_type == LESType::Smagorinsky ||
        solverChoice.les_type == LESType::Deardorff)
    {
        amrex::Real inv_Pr_t = solverChoice.Pr_t_inv;
        amrex::Real inv_Sc_t = solverChoice.Sc_t_inv;
        amrex::Real inv_sigma_k = 1.0 / solverChoice.sigma_k;
        bool use_KE = (solverChoice.les_type == LESType::Deardorff);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(1);
            if (l_vert_only)
                bx.setRange(2,klo,1);

            const amrex::Array4<amrex::Real> &K = eddyViscosity.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Get eddy diffusivities from the eddy viscosity
                // Additional factor of 0.5 because K = 2 mu_t
                K(i,j,k,EddyDiff::Theta_h) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_Pr_t;
                K(i,j,k,EddyDiff::Scalar_h) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_Sc_t;
                if (use_KE) {
                    K(i,j,k,EddyDiff::KE_h) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_sigma_k;
                }
                if (solverChoice.use_QKE && solverChoice.diffuse_QKE_3D) {
                    K(i,j,k,EddyDiff::QKE_h) = 0.5 * K(i,j,k,EddyDiff::Mom_h) * inv_sigma_k;
                }

                // For LES: vertical and horizontal components are the same
                K(i,j,k,EddyDiff::Mom_v) = K(i,j,k,EddyDiff::Mom_h);
                K(i,j,k,EddyDiff::Theta_v) = K(i,j,k,EddyDiff::Theta_h);
                K(i,j,k,EddyDiff::Scalar_v) = K(i,j,k,EddyDiff::Scalar_h);
                if (use_KE) {
                    K(i,j,k,EddyDiff::KE_v) = K(i,j,k,EddyDiff::KE_h);
                }
                // QKE: vertical diffusion from PBL model
            });
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const amrex::Box &bx = mfi.growntilebox(1);

        if (!(domain.contains(bx))) {
            // Fill values outside the domain by straightforward extrapolation.  Note this must be
            // done separately from the loop above so all the interior values are filled.  We also
            // do this in three separate loops so that we don't have any race conditions.
            const amrex::Array4<amrex::Real> &K = eddyViscosity.array(mfi);
            amrex::ParallelFor(bx, (int) EddyDiff::NumDiffs, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                if (i < dom_lo.x) {
                    if (j < dom_lo.y)
                        K(i,j,k,n) = K(dom_lo.x,dom_lo.y,k,n);
                    else if (j > dom_hi.y)
                        K(i,j,k,n) = K(dom_lo.x,dom_hi.y,k,n);
                    else
                       K(i,j,k,n) = K(dom_lo.x,j,k,n);
                } else if (i > dom_hi.x) {
                    if (j < dom_lo.y)
                        K(i,j,k,n) = K(dom_hi.x,dom_lo.y,k,n);
                    else if (j > dom_hi.y)
                        K(i,j,k,n) = K(dom_hi.x,dom_hi.y,k,n);
                    else
                        K(i,j,k,n) = K(dom_hi.x,j,k,n);
                } else if (j < dom_lo.y) {
                        K(i,j,k,n) = K(i,dom_lo.y,k,n);
                } else if (j > dom_hi.y) {
                        K(i,j,k,n) = K(i,dom_hi.y,k,n);
                }
           });
        }
    } //mfi

    // Fill interior ghost cells and any ghost cells outside a periodic domain
    eddyViscosity.FillBoundary(geom.periodicity());

    // Now extend to low and high in vertical
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const amrex::Box &bx = mfi.growntilebox(1);

        if (!(domain.contains(bx))) {
            const amrex::Array4<amrex::Real> &K = eddyViscosity.array(mfi);
            amrex::ParallelFor(bx, (int) EddyDiff::NumDiffs, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k < dom_lo.z) K(i,j,k,n) = K(i,j,dom_lo.z,n);
                if (k > dom_hi.z) K(i,j,k,n) = K(i,j,dom_hi.z,n);
            });
        }
    } //mfi

    eddyViscosity.FillBoundary(geom.periodicity());

} // ComputeTurbulentViscosityLES

void ComputeTurbulentViscosity (const amrex::MultiFab& xvel, const amrex::MultiFab& yvel, const amrex::MultiFab& zvel,
                                const amrex::MultiFab& cons_in, amrex::MultiFab& eddyViscosity,
                                const amrex::Geometry& geom,
                                const SolverChoice& solverChoice, std::unique_ptr<ABLMost>& most,
                                const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                                bool vert_only)
{
    BL_PROFILE_VAR("ComputeTurbulentViscosity()",ComputeTurbulentViscosity);
    //
    // In LES mode, the turbulent viscosity is isotropic, so the LES model sets both horizontal and vertical viscosities
    // In PBL mode, the primary purpose of the PBL model is to control vertical transport, so the PBL model sets the vertical viscosity.
    // Optionally, the PBL model can be run in conjunction with an LES model that sets the horizontal viscosity
    // (this isn’t truly LES, but the model form is the same as Smagorinsky).
    //
    // ComputeTurbulentViscosityLES populates the LES viscosity for both horizontal and vertical components.
    // ComputeTurbulentViscosityPBL computes the PBL viscosity just for the vertical component.
    //

    // We must initialize all the components to 0 because we may not set all of them below
    //    (which ones depends on which LES / PBL scheme we are using)
    eddyViscosity.setVal(0.0);

    if (most) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(solverChoice.les_type == LESType::Smagorinsky ||
                                         solverChoice.les_type == LESType::Deardorff   ||
                                         solverChoice.pbl_type == PBLType::MYNN25,
                                         "Must use an LES or PBL model to compute turbulent viscosity for MOST boundaries");
    } else {
        AMREX_ALWAYS_ASSERT(!vert_only);
    }

    if ( (  vert_only  &&
           (solverChoice.pbl_type == PBLType::None) ) ||
         (  !vert_only &&
          ((solverChoice.les_type == LESType::Smagorinsky) ||
           (solverChoice.les_type == LESType::Deardorff  )) ) ) {

            ComputeTurbulentViscosityLES(xvel, yvel, zvel, cons_in, eddyViscosity,
                                         geom, solverChoice,
                                         domain_bcs_type_d, vert_only);
    }

    if (solverChoice.pbl_type != PBLType::None) {
         ComputeTurbulentViscosityPBL(xvel, yvel, cons_in, eddyViscosity, geom, solverChoice, most, vert_only);
    }
}
