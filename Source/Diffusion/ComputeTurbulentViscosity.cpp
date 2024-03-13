/** \file ComputeTurbulentViscosity.cpp */

#include <ABLMost.H>
#include <EddyViscosity.H>
#include <Diffusion.H>
#include <TileNoZ.H>
#include <TerrainMetrics.H>

using namespace amrex;

void
ComputeTurbulentViscosityPBL (const MultiFab& xvel,
                              const MultiFab& yvel,
                              const MultiFab& cons_in,
                              MultiFab& eddyViscosity,
                              const Geometry& geom,
                              const TurbChoice& turbChoice,
                              std::unique_ptr<ABLMost>& most,
                              const BCRec* bc_ptr,
                              bool /*vert_only*/,
                              const std::unique_ptr<MultiFab>& z_phys_nd);

/**
 * Function for computing the turbulent viscosity with LES.
 *
 * @param[in]  Tau11 11 strain
 * @param[in]  Tau22 22 strain
 * @param[in]  Tau33 33 strain
 * @param[in]  Tau12 12 strain
 * @param[in]  Tau13 13 strain
 * @param[in]  Tau23 23 strain
 * @param[in]  cons_in cell center conserved quantities
 * @param[out] eddyViscosity turbulent viscosity
 * @param[in]  Hfx1 heat flux in x-dir
 * @param[in]  Hfx2 heat flux in y-dir
 * @param[in]  Hfx3 heat flux in z-dir
 * @param[in]  Diss dissipation of turbulent kinetic energy
 * @param[in]  geom problem geometry
 * @param[in]  mapfac_u map factor at x-face
 * @param[in]  mapfac_v map factor at y-face
 * @param[in]  turbChoice container with turbulence parameters
 */
void ComputeTurbulentViscosityLES (const MultiFab& Tau11, const MultiFab& Tau22, const MultiFab& Tau33,
                                   const MultiFab& Tau12, const MultiFab& Tau13, const MultiFab& Tau23,
                                   const MultiFab& cons_in, MultiFab& eddyViscosity,
                                   MultiFab& Hfx1, MultiFab& Hfx2, MultiFab& Hfx3, MultiFab& Diss,
                                   const Geometry& geom,
                                   const MultiFab& mapfac_u, const MultiFab& mapfac_v,
                                   const std::unique_ptr<MultiFab>& z_phys_nd,
                                   const TurbChoice& turbChoice, const Real const_grav, std::unique_ptr<ABLMost>& most)
{
    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv = geom.InvCellSizeArray();
    const Box& domain = geom.Domain();
    const int& klo    = domain.smallEnd(2);
    const bool use_most = (most != nullptr);
    const bool use_terrain = (z_phys_nd != nullptr);

    // SMAGORINSKY: Fill Kturb for momentum in horizontal and vertical
    //***********************************************************************************
    if (turbChoice.les_type == LESType::Smagorinsky)
    {
      Real Cs = turbChoice.Cs;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          // NOTE: This gets us the lateral ghost cells for lev>0; which
          //       have been filled from FP Two Levels.
          Box bxcc  = mfi.growntilebox() & domain;

          const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
          const Array4<Real const > &cell_data = cons_in.array(mfi);

          Array4<Real const> tau11 = Tau11.array(mfi);
          Array4<Real const> tau22 = Tau22.array(mfi);
          Array4<Real const> tau33 = Tau33.array(mfi);
          Array4<Real const> tau12 = Tau12.array(mfi);
          Array4<Real const> tau13 = Tau13.array(mfi);
          Array4<Real const> tau23 = Tau23.array(mfi);

          Array4<Real const> mf_u = mapfac_u.array(mfi);
          Array4<Real const> mf_v = mapfac_v.array(mfi);

          Array4<Real const> z_nd_arr = (use_terrain) ? z_phys_nd->const_array(mfi) : Array4<Real const>{};

          ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              Real SmnSmn = ComputeSmnSmn(i,j,k,tau11,tau22,tau33,tau12,tau13,tau23,klo,use_most);
              Real dxInv = cellSizeInv[0];
              Real dyInv = cellSizeInv[1];
              Real dzInv = cellSizeInv[2];
              if (use_terrain) {
                  // the terrain grid is only deformed in z for now
                  dzInv /= Compute_h_zeta_AtCellCenter(i,j,k, cellSizeInv, z_nd_arr);
              }
              Real cellVolMsf = 1.0 / (dxInv * mf_u(i,j,0) * dyInv * mf_v(i,j,0) * dzInv);
              Real DeltaMsf   = std::pow(cellVolMsf,1.0/3.0);
              Real CsDeltaSqrMsf = Cs*Cs*DeltaMsf*DeltaMsf;

              mu_turb(i, j, k, EddyDiff::Mom_h) = CsDeltaSqrMsf * cell_data(i, j, k, Rho_comp) * std::sqrt(2.0*SmnSmn);
              mu_turb(i, j, k, EddyDiff::Mom_v) = mu_turb(i, j, k, EddyDiff::Mom_h);
          });
      }
    }
    // DEARDORFF: Fill Kturb for momentum in horizontal and vertical
    //***********************************************************************************
    else if (turbChoice.les_type == LESType::Deardorff)
    {
      const Real l_C_k        = turbChoice.Ck;
      const Real l_C_e        = turbChoice.Ce;
      const Real l_C_e_wall   = turbChoice.Ce_wall;
      const Real Ce_lcoeff    = amrex::max(0.0, l_C_e - 1.9*l_C_k);
      const Real l_abs_g      = const_grav;
      const Real l_inv_theta0 = 1.0 / turbChoice.theta_ref;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          Box bxcc  = mfi.tilebox();

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
        const Array4<Real>& hfx_x   = Hfx1.array(mfi);
        const Array4<Real>& hfx_y   = Hfx2.array(mfi);
        const Array4<Real>& hfx_z   = Hfx3.array(mfi);
        const Array4<Real>& diss    = Diss.array(mfi);

        const Array4<Real const > &cell_data = cons_in.array(mfi);

        Array4<Real const> mf_u = mapfac_u.array(mfi);
        Array4<Real const> mf_v = mapfac_v.array(mfi);

        Array4<Real const> z_nd_arr = (use_terrain) ? z_phys_nd->const_array(mfi) : Array4<Real const>{};

        ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          Real dxInv = cellSizeInv[0];
          Real dyInv = cellSizeInv[1];
          Real dzInv = cellSizeInv[2];
          if (use_terrain) {
              // the terrain grid is only deformed in z for now
              dzInv /= Compute_h_zeta_AtCellCenter(i,j,k, cellSizeInv, z_nd_arr);
          }
          Real cellVolMsf = 1.0 / (dxInv * mf_u(i,j,0) * dyInv * mf_v(i,j,0) * dzInv);
          Real DeltaMsf   = std::pow(cellVolMsf,1.0/3.0);

          // Calculate stratification-dependent mixing length (Deardorff 1980)
          Real eps       = std::numeric_limits<Real>::epsilon();
          Real dtheta_dz;
          if (use_most && k==klo) {
#ifdef ERF_EXPLICIT_MOST_STRESS
              dtheta_dz = ( cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp)
                          - cell_data(i,j,k  ,RhoTheta_comp)/cell_data(i,j,k  ,Rho_comp) )*dzInv;
#else
              dtheta_dz = 0.5 * (-3 * cell_data(i,j,k  ,RhoTheta_comp)
                                    / cell_data(i,j,k  ,Rho_comp)
                                + 4 * cell_data(i,j,k+1,RhoTheta_comp)
                                    / cell_data(i,j,k+1,Rho_comp)
                                -     cell_data(i,j,k+2,RhoTheta_comp)
                                    / cell_data(i,j,k+2,Rho_comp) ) * dzInv;
#endif
          } else {
              dtheta_dz = 0.5 * ( cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp)
                                - cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp) )*dzInv;
          }
          Real E         = cell_data(i,j,k,RhoKE_comp) / cell_data(i,j,k,Rho_comp);
          Real strat     = l_abs_g * dtheta_dz * l_inv_theta0; // stratification
          Real length;
          if (strat <= eps) {
              length = DeltaMsf;
          } else {
              length = 0.76 * std::sqrt(E / strat);
              // mixing length should be _reduced_ for stable stratification
              length = amrex::min(length, DeltaMsf);
              // following WRF, make sure the mixing length isn't too small
              length = amrex::max(length, 0.001 * DeltaMsf);
          }

          // Calculate eddy diffusivities
          // K = rho * C_k * l * KE^(1/2)
          mu_turb(i,j,k,EddyDiff::Mom_h) = cell_data(i,j,k,Rho_comp) * l_C_k * length * std::sqrt(E);
          mu_turb(i,j,k,EddyDiff::Mom_v) = mu_turb(i,j,k,EddyDiff::Mom_h);
          // KH = (1 + 2*l/delta) * mu_turb
          mu_turb(i,j,k,EddyDiff::Theta_v) = (1.+2.*length/DeltaMsf) * mu_turb(i,j,k,EddyDiff::Mom_v);

          // Calculate SFS quantities
          // - dissipation
          Real Ce;
          if ((l_C_e_wall > 0) && (k==0))
              Ce = l_C_e_wall;
          else
              Ce = 1.9*l_C_k + Ce_lcoeff*length / DeltaMsf;
          diss(i,j,k) = cell_data(i,j,k,Rho_comp) * Ce * std::pow(E,1.5) / length;
          // - heat flux
          //   (Note: If using ERF_EXPLICIT_MOST_STRESS, the value at k=0 will
          //    be overwritten when BCs are applied)
          hfx_x(i,j,k) = 0.0;
          hfx_y(i,j,k) = 0.0;
          hfx_z(i,j,k) = -mu_turb(i,j,k,EddyDiff::Theta_v) * dtheta_dz; // (rho*w)' theta' [kg m^-2 s^-1 K]
        });
      }
    }

    // Extrapolate Kturb in x/y, fill remaining elements (relevent to lev==0)
    //***********************************************************************************
    int ngc(1);
    Real inv_Pr_t    = turbChoice.Pr_t_inv;
    Real inv_Sc_t    = turbChoice.Sc_t_inv;
    Real inv_sigma_k = 1.0 / turbChoice.sigma_k;
    // EddyDiff mapping :   Theta_h     KE_h         QKE_h      Scalar_h    Q_h
    Vector<Real> Factors = {inv_Pr_t, inv_sigma_k, inv_sigma_k, inv_Sc_t, inv_Sc_t}; // alpha = mu/Pr
    Gpu::AsyncVector<Real> d_Factors; d_Factors.resize(Factors.size());
    Gpu::copy(Gpu::hostToDevice, Factors.begin(), Factors.end(), d_Factors.begin());
    Real* fac_ptr = d_Factors.data();

    bool use_KE  = (turbChoice.les_type == LESType::Deardorff);
    bool use_QKE = (turbChoice.use_QKE && turbChoice.diffuse_QKE_3D);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bxcc   = mfi.tilebox();
        Box planex = bxcc; planex.setSmall(0, 1); planex.setBig(0, ngc); planex.grow(1,1);
        Box planey = bxcc; planey.setSmall(1, 1); planey.setBig(1, ngc); planey.grow(0,1);
        int i_lo   = bxcc.smallEnd(0); int i_hi = bxcc.bigEnd(0);
        int j_lo   = bxcc.smallEnd(1); int j_hi = bxcc.bigEnd(1);
        bxcc.growLo(0,ngc); bxcc.growHi(0,ngc);
        bxcc.growLo(1,ngc); bxcc.growHi(1,ngc);

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);

        // Extrapolate outside the domain in lateral directions (planex owns corner cells)
        if (i_lo == domain.smallEnd(0)) {
            ParallelFor(planex, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int lj = amrex::min(amrex::max(j, domain.smallEnd(1)), domain.bigEnd(1));
                mu_turb(i_lo-i, j, k, EddyDiff::Mom_h) = mu_turb(i_lo, lj, k, EddyDiff::Mom_h);
                mu_turb(i_lo-i, j, k, EddyDiff::Mom_v) = mu_turb(i_lo, lj, k, EddyDiff::Mom_v);
            });
        }
        if (i_hi == domain.bigEnd(0)) {
            ParallelFor(planex, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int lj = amrex::min(amrex::max(j, domain.smallEnd(1)), domain.bigEnd(1));
                mu_turb(i_hi+i, j, k, EddyDiff::Mom_h) = mu_turb(i_hi, lj, k, EddyDiff::Mom_h);
                mu_turb(i_hi+i, j, k, EddyDiff::Mom_v) = mu_turb(i_hi, lj, k, EddyDiff::Mom_v);
            });
        }
        if (j_lo == domain.smallEnd(1)) {
            ParallelFor(planey, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int li = amrex::min(amrex::max(i, domain.smallEnd(0)), domain.bigEnd(0));
                mu_turb(i, j_lo-j, k, EddyDiff::Mom_h) = mu_turb(li, j_lo, k, EddyDiff::Mom_h);
                mu_turb(i, j_lo-j, k, EddyDiff::Mom_v) = mu_turb(li, j_lo, k, EddyDiff::Mom_v);
            });
        }
        if (j_hi == domain.bigEnd(1)) {
            ParallelFor(planey, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int li = amrex::min(amrex::max(i, domain.smallEnd(0)), domain.bigEnd(0));
                mu_turb(i, j_hi+j, k, EddyDiff::Mom_h) = mu_turb(li, j_hi, k, EddyDiff::Mom_h);
                mu_turb(i, j_hi+j, k, EddyDiff::Mom_v) = mu_turb(li, j_hi, k, EddyDiff::Mom_v);
            });
        }

        // refactor the code to eliminate the need for ifdef's
        for (auto n = 1; n < (EddyDiff::NumDiffs-1)/2; ++n) {
            int offset = (EddyDiff::NumDiffs-1)/2;
            switch (n)
            {
              case EddyDiff::QKE_h:
                 // Populate element other than mom_h/v on the whole grid
                 if(use_QKE) {
                   ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                   {
                     int indx   = n;
                     int indx_v = indx + offset;
                     mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx-1];
                     mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
                  });
                 }
                 break;
             case EddyDiff::KE_h:
                if (use_KE) {
                   ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                   {
                     int indx   = n;
                     int indx_v = indx + offset;
                     mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx-1];
                     mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
                   });
                }
                break;
            default:
                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                  int indx   = n;
                  int indx_v = indx + offset;
                  mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx-1];
                  mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
                });
                break;
          }
       }
    }

    // Fill interior ghost cells and any ghost cells outside a periodic domain
    //***********************************************************************************
    eddyViscosity.FillBoundary(geom.periodicity());


    // Extrapolate top & bottom
    //***********************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(eddyViscosity,TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bxcc   = mfi.tilebox();
        Box planez = bxcc; planez.setSmall(2, 1); planez.setBig(2, ngc);
        int k_lo   = bxcc.smallEnd(2); int k_hi = bxcc.bigEnd(2);
        planez.growLo(0,ngc); planez.growHi(0,ngc);
        planez.growLo(1,ngc); planez.growHi(1,ngc);

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);

        // refactor the code to eliminate the need for ifdef's
        for (auto n = 0; n < (EddyDiff::NumDiffs-1)/2; ++n) {
            int offset = (EddyDiff::NumDiffs-1)/2;
            switch (n)
            {
              case EddyDiff::QKE_h:
                 // Extrap all components at top & bottom
                 if(use_QKE) {
                    ParallelFor(planez, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                      int indx = n;
                      int indx_v = indx + offset;
                      mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
                      mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
                      mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
                      mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
                    });
                 }
                 break;
              case EddyDiff::KE_h:
                 if (use_KE) {
                    ParallelFor(planez, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                      int indx   = n;
                      int indx_v = indx + offset;
                      mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
                      mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
                      mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
                      mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
                    });
                 }
                 break;
              default:
                 ParallelFor(planez, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                   int indx   = n ;
                   int indx_v = indx + offset;
                   mu_turb(i, j, k_lo-k, indx  ) = mu_turb(i, j, k_lo, indx  );
                   mu_turb(i, j, k_hi+k, indx  ) = mu_turb(i, j, k_hi, indx  );
                   mu_turb(i, j, k_lo-k, indx_v) = mu_turb(i, j, k_lo, indx_v);
                   mu_turb(i, j, k_hi+k, indx_v) = mu_turb(i, j, k_hi, indx_v);
                 });
                 break;
            }
        }
    }
}

/**
 * Wrapper to compute turbulent viscosity with LES or PBL.
 *
 * @param[in]  xvel velocity in x-dir
 * @param[in]  yvel velocity in y-dir
 * @param[in]  Tau11 11 strain
 * @param[in]  Tau22 22 strain
 * @param[in]  Tau33 33 strain
 * @param[in]  Tau12 12 strain
 * @param[in]  Tau13 13 strain
 * @param[in]  Tau23 23 strain
 * @param[in]  cons_in cell center conserved quantities
 * @param[out] eddyViscosity turbulent viscosity
 * @param[in]  Hfx1 heat flux in x-dir
 * @param[in]  Hfx2 heat flux in y-dir
 * @param[in]  Hfx3 heat flux in z-dir
 * @param[in]  Diss dissipation of turbulent kinetic energy
 * @param[in]  geom problem geometry
 * @param[in]  mapfac_u map factor at x-face
 * @param[in]  mapfac_v map factor at y-face
 * @param[in]  turbChoice container with turbulence parameters
 * @param[in]  most pointer to Monin-Obukhov class if instantiated
 * @param[in]  vert_only flag for vertical components of eddyViscosity
 */
void ComputeTurbulentViscosity (const MultiFab& xvel , const MultiFab& yvel ,
                                const MultiFab& Tau11, const MultiFab& Tau22, const MultiFab& Tau33,
                                const MultiFab& Tau12, const MultiFab& Tau13, const MultiFab& Tau23,
                                const MultiFab& cons_in,
                                MultiFab& eddyViscosity,
                                MultiFab& Hfx1, MultiFab& Hfx2, MultiFab& Hfx3, MultiFab& Diss,
                                const Geometry& geom,
                                const MultiFab& mapfac_u, const MultiFab& mapfac_v,
                                const std::unique_ptr<MultiFab>& z_phys_nd,
                                const TurbChoice& turbChoice, const Real const_grav,
                                std::unique_ptr<ABLMost>& most,
                                const BCRec* bc_ptr,
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

    if (most) {
        bool l_use_turb = ( turbChoice.les_type == LESType::Smagorinsky ||
                            turbChoice.les_type == LESType::Deardorff   ||
                            turbChoice.pbl_type == PBLType::MYNN25      ||
                            turbChoice.pbl_type == PBLType::YSU );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(l_use_turb,
          "An LES or PBL model must be utilized with MOST boundaries to compute the turbulent viscosity");
    } else {
        AMREX_ALWAYS_ASSERT(!vert_only);
    }

    if (turbChoice.les_type != LESType::None) {
        ComputeTurbulentViscosityLES(Tau11, Tau22, Tau33,
                                     Tau12, Tau13, Tau23,
                                     cons_in, eddyViscosity,
                                     Hfx1, Hfx2, Hfx3, Diss,
                                     geom, mapfac_u, mapfac_v,
                                     z_phys_nd, turbChoice, const_grav, most);
    }

    if (turbChoice.pbl_type != PBLType::None) {
        // NOTE: state_new is passed in for Cons_old (due to ptr swap in advance)
        ComputeTurbulentViscosityPBL(xvel, yvel, cons_in, eddyViscosity,
                                     geom, turbChoice, most, bc_ptr, vert_only, z_phys_nd);
    }
}
