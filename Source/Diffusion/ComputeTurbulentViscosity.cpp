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
void ComputeTurbulentViscosityLES (const amrex::MultiFab& Tau11, const amrex::MultiFab& Tau22, const amrex::MultiFab& Tau33,
                                   const amrex::MultiFab& Tau12, const amrex::MultiFab& Tau13, const amrex::MultiFab& Tau23,
                                   const amrex::MultiFab& cons_in, amrex::MultiFab& eddyViscosity,
                                   const amrex::Geometry& geom,
                                   const SolverChoice& solverChoice)
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const Real cellVol = 1.0 / (dxInv[0] * dxInv[1] * dxInv[2]);
    const Real Delta   = std::pow(cellVol,1.0/3.0);

    // SMAGORINSKY: Fill Kturb for momentum in horizontal and vertical
    //***********************************************************************************
    if (solverChoice.les_type == LESType::Smagorinsky)
    {
      Real Cs = solverChoice.Cs;
      Real CsDeltaSqr = Cs*Cs*Delta*Delta;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          Box bxcc  = mfi.tilebox();

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
        const amrex::Array4<amrex::Real const > &cell_data = cons_in.array(mfi);

        Array4<Real const> tau11 = Tau11.array(mfi);
        Array4<Real const> tau22 = Tau22.array(mfi);
        Array4<Real const> tau33 = Tau33.array(mfi);
        Array4<Real const> tau12 = Tau12.array(mfi);
        Array4<Real const> tau13 = Tau13.array(mfi);
        Array4<Real const> tau23 = Tau23.array(mfi);

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

          mu_turb(i, j, k, EddyDiff::Mom_h) = CsDeltaSqr * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
          mu_turb(i, j, k, EddyDiff::Mom_v) = mu_turb(i, j, k, EddyDiff::Mom_h);
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
      for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          Box bxcc  = mfi.tilebox();

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
        const amrex::Array4<amrex::Real const > &cell_data = cons_in.array(mfi);

        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          // K = rho * C_k * Delta * KE^(1/2) = C_k * Delta * (rho * RhoKE)^1/2
          mu_turb(i,j,k,EddyDiff::Mom_h) = l_C_k * Delta *
            std::sqrt(cell_data(i,j,k,RhoKE_comp) * cell_data(i,j,k,Rho_comp));
          mu_turb(i, j, k, EddyDiff::Mom_v) = mu_turb(i, j, k, EddyDiff::Mom_h);
        });
      }
    }

    // Extrapolate Kturb in extrap x/y, fill remaining elements
    //***********************************************************************************
    int ngc(1);
    Real inv_Pr_t    = solverChoice.Pr_t_inv;
    Real inv_Sc_t    = solverChoice.Sc_t_inv;
    Real inv_sigma_k = 1.0 / solverChoice.sigma_k;
    Vector<Real> Factors = {inv_Pr_t, inv_Sc_t, inv_sigma_k, inv_sigma_k}; // alpha = mu/Pr
    Gpu::AsyncVector<Real> d_Factors; d_Factors.resize(Factors.size());
    Gpu::copy(Gpu::hostToDevice, Factors.begin(), Factors.end(), d_Factors.begin());
    Real* fac_ptr = d_Factors.data();

    bool use_KE  = (solverChoice.les_type == LESType::Deardorff);
    bool use_QKE = (solverChoice.use_QKE && solverChoice.diffuse_QKE_3D);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bxcc   = mfi.tilebox();
        Box planex = bxcc; planex.setSmall(0, 1); planex.setBig(0, ngc);
        Box planey = bxcc; planey.setSmall(1, 1); planey.setBig(1, ngc);
        int i_lo   = bxcc.smallEnd(0); int i_hi = bxcc.bigEnd(0);
        int j_lo   = bxcc.smallEnd(1); int j_hi = bxcc.bigEnd(1);
        bxcc.growLo(0,ngc); bxcc.growHi(0,ngc);
        bxcc.growLo(1,ngc); bxcc.growHi(1,ngc);

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);

        // Extrapolate x & y (over-written on interior ghost cells by FillBoundary)
        amrex::ParallelFor(planex, planey,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            mu_turb(i_lo-i, j, k, EddyDiff::Mom_h) = mu_turb(i_lo, j, k, EddyDiff::Mom_h);
            mu_turb(i_lo-i, j, k, EddyDiff::Mom_v) = mu_turb(i_lo, j, k, EddyDiff::Mom_v);

            mu_turb(i_hi+i, j, k, EddyDiff::Mom_h) = mu_turb(i_hi, j, k, EddyDiff::Mom_h);
            mu_turb(i_hi+i, j, k, EddyDiff::Mom_v) = mu_turb(i_hi, j, k, EddyDiff::Mom_v);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            mu_turb(i, j_lo-j, k, EddyDiff::Mom_h) = mu_turb(i, j_lo, k, EddyDiff::Mom_h);
            mu_turb(i, j_lo-j, k, EddyDiff::Mom_v) = mu_turb(i, j_lo, k, EddyDiff::Mom_v);

            mu_turb(i, j_hi+j, k, EddyDiff::Mom_h) = mu_turb(i, j_hi, k, EddyDiff::Mom_h);
            mu_turb(i, j_hi+j, k, EddyDiff::Mom_v) = mu_turb(i, j_hi, k, EddyDiff::Mom_v);
        });


        int ntot   = 7;
        int offset = EddyDiff::Theta_h;
        // Populate element other than mom_h/v on the whole grid
        if(use_QKE) {
          int ncomp  = 4;
          ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            int indx   = n + offset;
            int indx_v = indx + ntot;
            mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
            mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
          });
        } else if (use_KE) {
          int ncomp  = 3;
          ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            int indx   = n + offset;
            int indx_v = indx + ntot;
            mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
            mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
          });
        } else {
          int ncomp  = 2;
          ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            int indx   = n + offset;
            int indx_v = indx + ntot;
            mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
            mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
          });
#ifdef ERF_USE_MOISTURE
      { // Qt
            int ncomp  = 5;
            ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              int indx   = n + offset;
              int indx_v = indx + ntot;
              mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
              mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
            });
      }

      {// Qp
            int ncomp  = 6;
            ParallelFor(bxcc,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              int indx   = n + offset;
              int indx_v = indx + ntot;
              mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx];
              mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
            });
      }
#endif
        }
    }

    // Fill interior ghost cells and any ghost cells outside a periodic domain
    //***********************************************************************************
    eddyViscosity.FillBoundary(geom.periodicity());


    // Extrapolate top & bottom
    //***********************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( amrex::MFIter mfi(eddyViscosity,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bxcc   = mfi.tilebox();
        Box planez = bxcc; planez.setSmall(2, 1); planez.setBig(2, ngc);
        int k_lo   = bxcc.smallEnd(2); int k_hi = bxcc.bigEnd(2);
        planez.growLo(0,ngc); planez.growHi(0,ngc);
        planez.growLo(1,ngc); planez.growHi(1,ngc);

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);

        int ntot   = 7;
        int offset = EddyDiff::Mom_h;
        // Extrap all components at top & bottom
        if(use_QKE) {
          int ncomp  = 4;
          ParallelFor(planez,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            int indx = n + offset;
            int indx_v = indx + ntot;
            mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
            mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
            mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
            mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
          });
        } else if (use_KE) {
          int ncomp  = 3;
          ParallelFor(planez,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            int indx   = n + offset;
            int indx_v = indx + ntot;
            mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
            mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
            mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
            mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
          });
        } else {
          int ncomp  = 2;
          ParallelFor(planez,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            int indx   = n + offset;
            int indx_v = indx + ntot;
            mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
            mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
            mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
            mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
          });
#ifdef ERF_USE_MOISTURE
      { // Qt
            int ncomp  = 5;
            ParallelFor(planez,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              int indx   = n + offset;
              int indx_v = indx + ntot;
              mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
              mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
              mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
              mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
            });
      }

      { // Qp
            int ncomp  = 6;
            ParallelFor(planez,ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
              int indx   = n + offset;
              int indx_v = indx + ntot;
              mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
              mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
              mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
              mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
            });
      }
#endif
        }
    }
}

void ComputeTurbulentViscosity (const amrex::MultiFab& xvel , const amrex::MultiFab& yvel ,
                                const amrex::MultiFab& Tau11, const amrex::MultiFab& Tau22, const amrex::MultiFab& Tau33,
                                const amrex::MultiFab& Tau12, const amrex::MultiFab& Tau13, const amrex::MultiFab& Tau23,
                                const amrex::MultiFab& cons_in,
                                amrex::MultiFab& eddyViscosity,
                                const amrex::Geometry& geom,
                                const SolverChoice& solverChoice,
                                std::unique_ptr<ABLMost>& most,
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
        bool l_use_turb = ( solverChoice.les_type == LESType::Smagorinsky ||
                            solverChoice.les_type == LESType::Deardorff   ||
                            solverChoice.pbl_type == PBLType::MYNN25 );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l_use_turb,
          "An LES or PBL model must be utilized with MOST boundaries to compute the turbulent viscosity");
    } else {
        AMREX_ALWAYS_ASSERT(!vert_only);
    }

    if (solverChoice.les_type != LESType::None) {
        ComputeTurbulentViscosityLES(Tau11, Tau22, Tau33,
                                     Tau12, Tau13, Tau23,
                                     cons_in, eddyViscosity,
                                     geom, solverChoice);
    }

    if (solverChoice.pbl_type != PBLType::None) {
        ComputeTurbulentViscosityPBL(xvel, yvel, cons_in, eddyViscosity,
                                     geom, solverChoice, most, vert_only);
    }
}
