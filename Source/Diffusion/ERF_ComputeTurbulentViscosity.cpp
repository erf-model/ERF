/** \file ERF_ComputeTurbulentViscosity.cpp */

#include "ERF_SurfaceLayer.H"
#include "ERF_EddyViscosity.H"
#include "ERF_Diffusion.H"
#include "ERF_PBLModels.H"
#include "ERF_TileNoZ.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_MoistUtils.H"
#include "ERF_RichardsonNumber.H"

using namespace amrex;

/**
 * Function for computing the turbulent viscosity with LES.
 *
 * @param[in]  Tau_lev strain at this level
 * @param[in]  cons_in cell center conserved quantities
 * @param[out] eddyViscosity turbulent viscosity
 * @param[in]  Hfx1 heat flux in x-dir
 * @param[in]  Hfx2 heat flux in y-dir
 * @param[in]  Hfx3 heat flux in z-dir
 * @param[in]  Diss dissipation of turbulent kinetic energy
 * @param[in]  geom problem geometry
 * @param[in]  mapfac map factors
 * @param[in]  turbChoice container with turbulence parameters
 * @param[in]  xvel x-direction velocity (for moist Ri correction)
 * @param[in]  yvel y-direction velocity (for moist Ri correction)
 */
void ComputeTurbulentViscosityLES (Vector<std::unique_ptr<MultiFab>>& Tau_lev,
                                   const MultiFab& cons_in, MultiFab& eddyViscosity,
                                   MultiFab& Hfx1, MultiFab& Hfx2, MultiFab& Hfx3, MultiFab& Diss,
                                  const Geometry& geom, bool use_terrain_fitted_coords,
                                  Vector<std::unique_ptr<MultiFab>>& mapfac,
                                  const std::unique_ptr<MultiFab>& z_phys_nd,
                                  const TurbChoice& turbChoice, const Real const_grav,
                                  std::unique_ptr<SurfaceLayer>& /*SurfLayer*/,
                                  const MoistureComponentIndices& moisture_indices,
                                  const MultiFab* xvel,
                                  const MultiFab* yvel)
{
    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv = geom.InvCellSizeArray();
    const Box& domain = geom.Domain();

    Real inv_Pr_t    = turbChoice.Pr_t_inv;
    Real inv_Sc_t    = turbChoice.Sc_t_inv;
    Real inv_sigma_k = 1.0 / turbChoice.sigma_k;

    bool use_thetav_grad = (turbChoice.strat_type == StratType::thetav);
    bool use_thetal_grad = (turbChoice.strat_type == StratType::thetal);

    bool isotropic = turbChoice.mix_isotropic;

    // SMAGORINSKY: Fill Kturb for momentum in horizontal and vertical
    //***********************************************************************************
    if (turbChoice.les_type == LESType::Smagorinsky)
    {
        Real Cs = turbChoice.Cs;
        bool smag2d = turbChoice.smag2d;

        // Define variables required inside device lambdas (scalars only)
        Real l_abs_g = const_grav;
        Real l_Ri_crit = turbChoice.Ri_crit;
        bool l_use_Ri_corr = turbChoice.use_Ri_correction;
        bool l_has_xvel = (xvel != nullptr);
        bool l_has_yvel = (yvel != nullptr);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bxcc  = mfi.growntilebox(1) & domain;
            const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
            const Array4<Real>& hfx_x   = Hfx1.array(mfi);
            const Array4<Real>& hfx_y   = Hfx2.array(mfi);
            const Array4<Real>& hfx_z   = Hfx3.array(mfi);
            const Array4<Real const >& cell_data = cons_in.array(mfi);
            Array4<Real const> tau11 = Tau_lev[TauType::tau11]->array(mfi);
            Array4<Real const> tau22 = Tau_lev[TauType::tau22]->array(mfi);
            Array4<Real const> tau33 = Tau_lev[TauType::tau33]->array(mfi);
            Array4<Real const> tau12 = Tau_lev[TauType::tau12]->array(mfi);
            Array4<Real const> tau13 = Tau_lev[TauType::tau13]->array(mfi);
            Array4<Real const> tau23 = Tau_lev[TauType::tau23]->array(mfi);
            Array4<Real const> mf_u = mapfac[MapFacType::u_x]->const_array(mfi);
            Array4<Real const> mf_v = mapfac[MapFacType::v_y]->const_array(mfi);
            Array4<Real const> z_nd_arr = z_phys_nd->const_array(mfi);

            Array4<Real const> u_arr = (l_has_xvel) ? xvel->const_array(mfi) : Array4<Real const>{};
            Array4<Real const> v_arr = (l_has_yvel) ? yvel->const_array(mfi) : Array4<Real const>{};

            ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // =====================================================================
                // 1. STRAIN RATE MAGNITUDE CALCULATION
                // =====================================================================
                Real SmnSmn;
                if (smag2d) {
                    SmnSmn = ComputeSmnSmn2D(i,j,k,tau11,tau22,tau12);
                } else {
                    SmnSmn = ComputeSmnSmn(i,j,k,tau11,tau22,tau33,tau12,tau13,tau23);
                }
                Real strain_rate_magnitude = std::sqrt(2.0 * SmnSmn);

                // =====================================================================
                // 2. GRID SCALE CALCULATION (filter width Δ)
                // =====================================================================
                Real dxInv = cellSizeInv[0];
                Real dyInv = cellSizeInv[1];
                Real dzInv = cellSizeInv[2];
                if (use_terrain_fitted_coords) {
                    dzInv /= Compute_h_zeta_AtCellCenter(i,j,k, cellSizeInv, z_nd_arr);
                }

                Real Delta;
                Real DeltaH;
                if (isotropic) {
                    Real cellVolMsf = 1.0 / (dxInv * mf_u(i,j,0) * dyInv * mf_v(i,j,0) * dzInv);
                    Delta = std::cbrt(cellVolMsf);
                    DeltaH = Delta;
                } else {
                    Delta = 1.0 / dzInv;
                    DeltaH = std::sqrt(1.0 / (dxInv * mf_u(i,j,0) * dyInv * mf_v(i,j,0)));
                }

                Real rho = cell_data(i, j, k, Rho_comp);
                Real CsDeltaSqr_h = Cs * Cs * DeltaH * DeltaH;
                Real CsDeltaSqr_v = Cs * Cs * Delta * Delta;

                Real nu_turb_base_h = CsDeltaSqr_h * strain_rate_magnitude;
                Real nu_turb_base_v = CsDeltaSqr_v * strain_rate_magnitude;

                Real stability_factor = 1.0;

                if (l_use_Ri_corr && l_has_xvel && l_has_yvel) {
                    Real N2 = ComputeN2(i, j, k, dzInv, l_abs_g, cell_data, moisture_indices);
                    Real S2_vert = ComputeVerticalShear2(i, j, k, dzInv, u_arr, v_arr);
                    Real Ri = ComputeRichardson(N2, S2_vert);
                    stability_factor = StabilityFunction(Ri, l_Ri_crit);
                }

                if (isotropic) {
                    mu_turb(i, j, k, EddyDiff::Mom_h) = rho * nu_turb_base_h * stability_factor;
                    mu_turb(i, j, k, EddyDiff::Mom_v) = rho * nu_turb_base_v * stability_factor;
                } else {
                    mu_turb(i, j, k, EddyDiff::Mom_h) = rho * nu_turb_base_h;
                    mu_turb(i, j, k, EddyDiff::Mom_v) = rho * nu_turb_base_v * stability_factor;
                }

                Real dtheta_dz = 0.5 * ( cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp)
                                       - cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp) )*dzInv;
                hfx_x(i,j,k) = 0.0;
                hfx_y(i,j,k) = 0.0;
                hfx_z(i,j,k) = -inv_Pr_t * mu_turb(i,j,k,EddyDiff::Mom_v) * dtheta_dz;
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

        const bool use_ref_theta = (turbChoice.theta_ref > 0);
        const Real l_inv_theta0  = (use_ref_theta) ? 1.0 / turbChoice.theta_ref : 1.0;

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

            Array4<Real const> mf_u = mapfac[MapFacType::u_x]->const_array(mfi);
            Array4<Real const> mf_v = mapfac[MapFacType::v_y]->const_array(mfi);

            Array4<Real const> z_nd_arr = z_phys_nd->const_array(mfi);

            ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real dxInv = cellSizeInv[0];
                Real dyInv = cellSizeInv[1];
                Real dzInv = cellSizeInv[2];
                if (use_terrain_fitted_coords) {
                    // the terrain grid is only deformed in z for now
                    dzInv /= Compute_h_zeta_AtCellCenter(i,j,k, cellSizeInv, z_nd_arr);
                }
                Real Delta;
                if (isotropic) {
                    Real cellVolMsf = 1.0 / (dxInv * mf_u(i,j,0) * dyInv * mf_v(i,j,0) * dzInv);
                    Delta = std::cbrt(cellVolMsf);
                } else {
                    Delta = 1.0 / dzInv;
                }

                Real dtheta_dz;
                if (use_thetav_grad) {
                    dtheta_dz = 0.5 * ( GetThetav(i, j, k+1, cell_data, moisture_indices)
                                      - GetThetav(i, j, k-1, cell_data, moisture_indices) )*dzInv;
                } else if (use_thetal_grad) {
                    dtheta_dz = 0.5 * ( GetThetal(i, j, k+1, cell_data, moisture_indices)
                                      - GetThetal(i, j, k-1, cell_data, moisture_indices) )*dzInv;
                } else {
                    dtheta_dz = 0.5 * ( cell_data(i, j, k+1, RhoTheta_comp) / cell_data(i, j, k+1, Rho_comp)
                                      - cell_data(i, j, k-1, RhoTheta_comp) / cell_data(i, j, k-1, Rho_comp) )*dzInv;
                }

                // Calculate stratification-dependent mixing length (Deardorff 1980, Eqn. 10a)
                Real E              = amrex::max(cell_data(i,j,k,RhoKE_comp)/cell_data(i,j,k,Rho_comp),Real(0.0));
                Real stratification = l_abs_g * dtheta_dz * l_inv_theta0;
                if (!use_ref_theta) {
                    // l_inv_theta0 == 1, divide by actual theta
                    stratification *= cell_data(i,j,k,Rho_comp) /
                                      cell_data(i,j,k,RhoTheta_comp);
                }

                // Following WRF, the stratification effects are applied to the vertical length scales
                // in the case of anistropic mixing
                Real length;
                Real eps       = std::numeric_limits<Real>::epsilon();
                if (stratification <= eps) {
                    length = Delta;  // cbrt(dx*dy*dz) -or- dz
                } else {
                    length = 0.76 * std::sqrt(E / amrex::max(stratification,eps));
                    // mixing length should be _reduced_ for stable stratification
                    length = amrex::min(length, Delta);
                    // following WRF, make sure the mixing length isn't too small
                    length = amrex::max(length, 0.001 * Delta);
                }

                Real DeltaH = (isotropic) ? length : std::sqrt(1.0 / (dxInv * mf_u(i,j,0) * dyInv * mf_v(i,j,0)));

                Real Pr_inv_v = (1. + 2.*length/Delta);
                Real Pr_inv_h  = (isotropic) ? Pr_inv_v : inv_Pr_t;

                // Calculate eddy diffusivities
                // K = rho * C_k * l * KE^(1/2)
                mu_turb(i,j,k,EddyDiff::Mom_h) = cell_data(i,j,k,Rho_comp) * l_C_k * DeltaH * std::sqrt(E);
                mu_turb(i,j,k,EddyDiff::Mom_v) = cell_data(i,j,k,Rho_comp) * l_C_k * length * std::sqrt(E);
                // KH = (1 + 2*l/delta) * mu_turb
                mu_turb(i,j,k,EddyDiff::Theta_h) = Pr_inv_h * mu_turb(i,j,k,EddyDiff::Mom_h);
                mu_turb(i,j,k,EddyDiff::Theta_v) = Pr_inv_v * mu_turb(i,j,k,EddyDiff::Mom_v);
                // Store lengthscale for TKE source terms
                mu_turb(i,j,k,EddyDiff::Turb_lengthscale) = length;

                // Calculate SFS quantities
                // - dissipation
                Real Ce;
                if ((l_C_e_wall > 0) && (k==0)) {
                    Ce = l_C_e_wall;
                } else {
                    Ce = 1.9*l_C_k + Ce_lcoeff*length / Delta;
                }
                diss(i,j,k) = cell_data(i,j,k,Rho_comp) * Ce * std::pow(E,1.5) / length;

                // - heat flux
                //   (Note: If using SurfaceLayer, the value at k=0 will
                //    be overwritten)
                hfx_x(i,j,k) = 0.0;
                hfx_y(i,j,k) = 0.0;
                hfx_z(i,j,k) = -mu_turb(i,j,k,EddyDiff::Theta_v) * dtheta_dz; // (rho*w)' theta' [kg m^-2 s^-1 K]
            });
        }
    }

    // Extrapolate Kturb in x/y, fill remaining elements (relevant to lev==0)
    //***********************************************************************************
    int ngc(1);
    // EddyDiff mapping :   Theta_h     KE_h       Scalar_h    Q_h
    Vector<Real> Factors = {inv_Pr_t, inv_sigma_k, inv_Sc_t, inv_Sc_t}; // alpha = mu/Pr
    Gpu::AsyncVector<Real> d_Factors; d_Factors.resize(Factors.size());
    Gpu::copy(Gpu::hostToDevice, Factors.begin(), Factors.end(), d_Factors.begin());
    Real* fac_ptr = d_Factors.data();

    const bool use_KE = ( turbChoice.les_type == LESType::Deardorff );

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bxcc   = mfi.tilebox();
        Box planex = bxcc; planex.setSmall(0, 1); planex.setBig(0, ngc); planex.grow(1,1);
        Box planey = bxcc; planey.setSmall(1, 1); planey.setBig(1, ngc); planey.grow(0,1);
        bxcc.growLo(0,ngc); bxcc.growHi(0,ngc);
        bxcc.growLo(1,ngc); bxcc.growHi(1,ngc);

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);

        for (auto n = 1; n < (EddyDiff::NumDiffs-1)/2; ++n) {
            int offset = (EddyDiff::NumDiffs-1)/2;
            switch (n)
            {
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

                    // NOTE: Theta_h, Theta_v have already been set for Deardorff
                    if (!(indx_v == EddyDiff::Theta_v && use_KE)) {
                        mu_turb(i,j,k,indx)   = mu_turb(i,j,k,EddyDiff::Mom_h) * fac_ptr[indx-1];
                        mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
                    }
                });
                break;
          }
       }
    }
}

/**
 * Function for computing the eddy viscosity with RANS.
 *
 * @param[in]  Tau_lev strain at this level
 * @param[in]  cons_in cell center conserved quantities
 * @param[out] eddyViscosity turbulent viscosity
 * @param[in]  Hfx1 heat flux in x-dir
 * @param[in]  Hfx2 heat flux in y-dir
 * @param[in]  Hfx3 heat flux in z-dir
 * @param[in]  Diss dissipation of turbulent kinetic energy
 * @param[in]  geom problem geometry
 * @param[in]  mapfac map factor
 * @param[in]  turbChoice container with turbulence parameters
 */
void ComputeTurbulentViscosityRANS (Vector<std::unique_ptr<MultiFab>>& /*Tau_lev*/,
                                    const MultiFab& cons_in,
                                    const MultiFab& wdist,
                                    MultiFab& eddyViscosity,
                                    MultiFab& Hfx1,
                                    MultiFab& Hfx2,
                                    MultiFab& Hfx3,
                                    MultiFab& Diss,
                                    const Geometry& geom,
                                    bool use_terrain_fitted_coords,
                                    Vector<std::unique_ptr<MultiFab>>& /*mapfac*/,
                                    const std::unique_ptr<MultiFab>& z_phys_nd,
                                    const TurbChoice& turbChoice,
                                    const Real const_grav,
                                    std::unique_ptr<SurfaceLayer>& SurfLayer,
                                    const MultiFab* z_0)
{
    const GpuArray<Real, AMREX_SPACEDIM> cellSizeInv = geom.InvCellSizeArray();
    const bool use_SurfLayer = (SurfLayer != nullptr);

    Real inv_Pr_t    = turbChoice.Pr_t_inv;
    Real inv_Sc_t    = turbChoice.Sc_t_inv;
    Real inv_sigma_k = 1.0 / turbChoice.sigma_k;

    // One-Equation k model (Axell & Liungman 2001, Environ Fluid Mech)
    //***********************************************************************************
    if (turbChoice.rans_type == RANSType::kEqn)
    {
        const Real Cmu0       = turbChoice.Cmu0;
        const Real Cmu0_pow3  = Cmu0 * Cmu0 * Cmu0;
        const Real inv_Cb_sq  = 1.0 / (turbChoice.Cb * turbChoice.Cb);
        const Real Rt_crit    = turbChoice.Rt_crit;
        const Real Rt_min     = turbChoice.Rt_min;
        const Real l_g_max    = turbChoice.l_g_max;
        const Real abs_g      = const_grav;

        const bool use_ref_theta = (turbChoice.theta_ref > 0);
        const Real inv_theta0  = (use_ref_theta) ? 1.0 / turbChoice.theta_ref : 1.0;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bxcc = mfi.tilebox();

            const Array4<Real const>& d_arr  = wdist.const_array(mfi);
            const Array4<Real const>& z0_arr = (use_SurfLayer) ? z_0->const_array(mfi) : Array4<Real const>{};

            const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
            const Array4<Real>& hfx_x   = Hfx1.array(mfi);
            const Array4<Real>& hfx_y   = Hfx2.array(mfi);
            const Array4<Real>& hfx_z   = Hfx3.array(mfi);
            const Array4<Real>& diss    = Diss.array(mfi);

            const Array4<Real const>& cell_data = cons_in.array(mfi);

            const Array4<Real const>& z_nd_arr = z_phys_nd->const_array(mfi);

            ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real eps = std::numeric_limits<Real>::epsilon();
                Real tke = amrex::max(cell_data(i,j,k,RhoKE_comp)/cell_data(i,j,k,Rho_comp), eps);

                // Estimate stratification
                Real dzInv = cellSizeInv[2];
                if (use_terrain_fitted_coords) {
                    // the terrain grid is only deformed in z for now
                    dzInv /= Compute_h_zeta_AtCellCenter(i,j,k, cellSizeInv, z_nd_arr);
                }
                Real dtheta_dz = 0.5 * ( cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp)
                                       - cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp) )*dzInv;
                Real N2 = abs_g * inv_theta0 * dtheta_dz; // Brunt–Väisälä frequency squared
                if (!use_ref_theta) {
                    // inv_theta0 == 1, divide by actual theta
                    N2 *= cell_data(i,j,k,Rho_comp) /
                          cell_data(i,j,k,RhoTheta_comp);
                }

                // Geometric length scale (AL01, Eqn. 22)
                Real l_g = (z0_arr) ? KAPPA * (d_arr(i, j, k) + z0_arr(i, j, 0))
                                    : KAPPA * d_arr(i, j, k);

                // Enforce a maximum value
                l_g = l_g_max * l_g / (l_g_max + l_g);

                // Turbulent Richardson number (AL01, Eqn. 29)
                // using the old dissipation value
                Real diss0 = max(diss(i, j, k) / cell_data(i, j, k, Rho_comp),
                                 eps);
                Real Rt = tke*tke * N2 / diss0;

                // Turbulent length scale
                Real length;
                if (std::abs(N2) <= eps) {
                    length = l_g;
                } else if (N2 > eps) {
                    // Stable (AL01, Eqn. 26)
                    length = std::sqrt(1.0 /
                            (1.0 / (l_g * l_g) + inv_Cb_sq * N2 / tke));
                } else {
                    // Unstable (AL01, Eqn. 28)
                    length = l_g * std::sqrt(1.0 - Cmu0_pow3*Cmu0_pow3 * inv_Cb_sq * Rt);
                }
                mu_turb(i, j, k, EddyDiff::Turb_lengthscale) = length;

                // Burchard & Petersen smoothing function
                Rt = (Rt >= Rt_crit) ? Rt
                                     : std::max(Rt,
                                                Rt - std::pow(Rt - Rt_crit, 2) /
                                                     (Rt + Rt_min - 2*Rt_crit));

                // Stability functions
                // Note: These use the smoothed turbulent Richardson number
                Real cmu = (Cmu0 + 0.108*Rt)
                         / (1.0 + 0.308*Rt + 0.00837*Rt*Rt); // (AL01, Eqn. 31)
                Real cmu_prime = Cmu0 / (1 + 0.277*Rt); // (AL01, Eqn. 32)

                // Calculate eddy diffusivities
                // K = rho * nu_t = rho * c_mu * tke^(1/2) * length
                Real nut = cmu * std::sqrt(tke) * length; // eddy viscosity
                Real nut_prime = cmu_prime / cmu * nut;   // eddy diffusivity
                mu_turb(i, j, k, EddyDiff::Mom_h) = cell_data(i, j, k, Rho_comp) * nut;
                mu_turb(i, j, k, EddyDiff::Mom_v) = mu_turb(i, j, k, EddyDiff::Mom_h);
                mu_turb(i, j, k, EddyDiff::Theta_v) = cell_data(i, j, k, Rho_comp) * nut_prime;

                // Calculate SFS quantities
                // - dissipation (AL01 Eqn. 19)
                diss(i, j, k) = cell_data(i, j, k, Rho_comp) * Cmu0_pow3 * std::pow(tke,1.5) / length;

                // - heat flux
                //   (Note: If using SurfaceLayer, the value at k=0 will
                //    be overwritten)
                hfx_x(i, j, k) = 0.0;
                hfx_y(i, j, k) = 0.0;
                // Note: buoyant production = g/theta0 * hfx == -nut_prime * N^2 (c.f. AL01 Eqn. 15)
                //                          = nut_prime * g/theta0 * dtheta/dz
                //                  ==> hfx = nut_prime * dtheta/dz
                //   Our convention is such that dtheta/dz < 0 gives a positive
                //   (upward) heat flux.
                hfx_z(i, j, k) = -mu_turb(i, j, k, EddyDiff::Theta_v) * dtheta_dz; // (rho*w)' theta' [kg m^-2 s^-1 K]
            });
        }
    }

    // Extrapolate Kturb in x/y, fill remaining elements (relevant to lev==0)
    //***********************************************************************************
    int ngc(1);
    // EddyDiff mapping :   Theta_h     KE_h       Scalar_h    Q_h
    Vector<Real> Factors = {inv_Pr_t, inv_sigma_k, inv_Sc_t, inv_Sc_t}; // alpha = mu/Pr
    Gpu::AsyncVector<Real> d_Factors; d_Factors.resize(Factors.size());
    Gpu::copy(Gpu::hostToDevice, Factors.begin(), Factors.end(), d_Factors.begin());
    Real* fac_ptr = d_Factors.data();

    const bool use_KE = ( turbChoice.rans_type == RANSType::kEqn );

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bxcc   = mfi.tilebox();
        Box planex = bxcc; planex.setSmall(0, 1); planex.setBig(0, ngc); planex.grow(1,1);
        Box planey = bxcc; planey.setSmall(1, 1); planey.setBig(1, ngc); planey.grow(0,1);
        bxcc.growLo(0,ngc); bxcc.growHi(0,ngc);
        bxcc.growLo(1,ngc); bxcc.growHi(1,ngc);

        const Array4<Real>& mu_turb = eddyViscosity.array(mfi);

        for (auto n = 1; n < (EddyDiff::NumDiffs-1)/2; ++n) {
            int offset = (EddyDiff::NumDiffs-1)/2;
            switch (n)
            {
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

                    // NOTE: Theta_v has already been set for Deardorff
                    if (!(indx_v == EddyDiff::Theta_v && use_KE)) {
                        mu_turb(i,j,k,indx_v) = mu_turb(i,j,k,indx);
                    }
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
 * @param[in]  Tau_lev strain at this level
 * @param[in]  cons_in cell center conserved quantities
 * @param[out] eddyViscosity turbulent viscosity
 * @param[in]  Hfx1 heat flux in x-dir
 * @param[in]  Hfx2 heat flux in y-dir
 * @param[in]  Hfx3 heat flux in z-dir
 * @param[in]  Diss dissipation of turbulent kinetic energy
 * @param[in]  geom problem geometry
 * @param[in]  mapfac map factors
 * @param[in]  turbChoice container with turbulence parameters
 * @param[in]  most pointer to Monin-Obukhov class if instantiated
 * @param[in]  vert_only flag for vertical components of eddyViscosity
 */
void ComputeTurbulentViscosity (Real dt,
                                const MultiFab& xvel, const MultiFab& yvel,
                                Vector<std::unique_ptr<MultiFab>>& Tau_lev,
                                MultiFab& cons_in,
                                const MultiFab& wdist,
                                MultiFab& eddyViscosity,
                                MultiFab& Hfx1, MultiFab& Hfx2, MultiFab& Hfx3, MultiFab& Diss,
                                const Geometry& geom,
                                Vector<std::unique_ptr<MultiFab>>& mapfac,
                                const std::unique_ptr<MultiFab>& z_phys_nd,
                                const SolverChoice& solverChoice,
                                std::unique_ptr<SurfaceLayer>& SurfLayer,
                                const MultiFab* z_0,
                                const bool& use_terrain_fitted_coords,
                                const bool& use_moisture,
                                int level,
                                const BCRec* bc_ptr,
                                bool vert_only)
{
    BL_PROFILE_VAR("ComputeTurbulentViscosity()",ComputeTurbulentViscosity);
    //
    // In LES mode, the turbulent viscosity is isotropic (unless mix_isotropic is set to false), so
    // the LES model sets both horizontal and vertical viscosities
    //
    // In PBL mode, the primary purpose of the PBL model is to control vertical transport, so the PBL model sets the vertical viscosity.
    // Optionally, the PBL model can be run in conjunction with an LES model that sets the horizontal viscosity
    // (this isn’t truly LES, but the model form is the same as Smagorinsky).
    //
    // ComputeTurbulentViscosityLES populates the LES viscosity for both horizontal and vertical components.
    // ComputeTurbulentViscosityPBL computes the PBL viscosity just for the vertical component.
    //

    TurbChoice turbChoice = solverChoice.turbChoice[level];
    const Real const_grav = solverChoice.gravity;

    if (!SurfLayer) {
        AMREX_ALWAYS_ASSERT(!vert_only);
    }

    bool impose_phys_bcs = false;

    if (turbChoice.les_type != LESType::None) {
        impose_phys_bcs = true;
        ComputeTurbulentViscosityLES(Tau_lev,
                                     cons_in, eddyViscosity,
                                     Hfx1, Hfx2, Hfx3, Diss,
                                     geom, use_terrain_fitted_coords,
                                     mapfac, z_phys_nd, turbChoice, const_grav,
                                     SurfLayer, solverChoice.moisture_indices,
                                     &xvel, &yvel);
    }

    if (turbChoice.rans_type != RANSType::None) {
        impose_phys_bcs = true;
        ComputeTurbulentViscosityRANS(Tau_lev,
                                      cons_in, wdist,
                                      eddyViscosity,
                                      Hfx1, Hfx2, Hfx3, Diss,
                                      geom, use_terrain_fitted_coords,
                                      mapfac, z_phys_nd, turbChoice, const_grav,
                                      SurfLayer, z_0);
    }

    if (turbChoice.pbl_type == PBLType::MYJ) {
        ComputeDiffusivityMYJ(dt, xvel, yvel, cons_in, eddyViscosity,
                              geom, turbChoice, SurfLayer,
                              use_terrain_fitted_coords, use_moisture,
                              level, bc_ptr, vert_only, z_phys_nd,
                              solverChoice.moisture_indices);
    } else if (turbChoice.pbl_type == PBLType::MYNN25) {
        ComputeDiffusivityMYNN25(xvel, yvel, cons_in, eddyViscosity,
                                 geom, turbChoice, SurfLayer,
                                 use_terrain_fitted_coords, use_moisture,
                                 level, bc_ptr, vert_only, z_phys_nd,
                                 solverChoice.moisture_indices);
    } else if (turbChoice.pbl_type == PBLType::MYNNEDMF) {
        ComputeDiffusivityMYNNEDMF(xvel, yvel, cons_in, eddyViscosity,
                                   geom, turbChoice, SurfLayer,
                                   use_terrain_fitted_coords, use_moisture,
                                   level, bc_ptr, vert_only, z_phys_nd,
                                   solverChoice.moisture_indices);
    } else if (turbChoice.pbl_type == PBLType::YSU) {
        ComputeDiffusivityYSU(xvel, yvel, cons_in, eddyViscosity,
                              geom, turbChoice, SurfLayer,
                              use_terrain_fitted_coords, use_moisture,
                              level, bc_ptr, vert_only, z_phys_nd,
                              solverChoice.moisture_indices);
    } else if (turbChoice.pbl_type == PBLType::MRF) {
        ComputeDiffusivityMRF(xvel, yvel, cons_in, eddyViscosity,
                              geom, turbChoice, SurfLayer,
                              use_terrain_fitted_coords, use_moisture,
                              level, bc_ptr, vert_only, z_phys_nd,
                              solverChoice.moisture_indices);
    } else if (turbChoice.pbl_type == PBLType::SHOC) {
        // NOTE: Nothing to do here. The SHOC class handles setting the vertical
        //       components of eddyDiffs in slow RHS pre.
    }

    //
    // At all levels we need to fill values outside the physical boundary for the LES coeffs.
    // In addition, for all cases, if at level > 0, we want to fill fine ghost cell values that
    // overlie coarse grid cells (and that are not in another fine valid region) with
    // extrapolated values from the interior, rather than interpolating from the coarser level,
    // since we may be using a different turbulence model there.
    //
    // Note: here "covered" refers to "covered by valid region of another grid at this level"
    // Note: here "physbnd" refers to "cells outside the domain if not periodic"
    // Note: here "interior" refers to "valid cells, i.e. inside 'my' grid"
    //
    {
        int is_covered    = 0;
        int is_notcovered = 1;
        int is_physbnd    = 2;
        int is_interior   = 3;
        iMultiFab cc_mask(eddyViscosity.boxArray(),eddyViscosity.DistributionMap(),1,1);
        cc_mask.BuildMask(geom.Domain(), geom.periodicity(), is_covered, is_notcovered, is_physbnd, is_interior);

        Box domain = geom.Domain();
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (geom.isPeriodic(i)) {
                domain.grow(i,1);
            }
        }

        eddyViscosity.FillBoundary(geom.periodicity());

        int ncomp  = eddyViscosity.nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(eddyViscosity); mfi.isValid(); ++mfi)
        {
            Box vbx = mfi.validbox();

            Box planex_lo = mfi.growntilebox(1); planex_lo.setBig(0, vbx.smallEnd(0)-1);
            Box planey_lo = mfi.growntilebox(1); planey_lo.setBig(1, vbx.smallEnd(1)-1);
            Box planez_lo = mfi.growntilebox(1); planez_lo.setBig(2, vbx.smallEnd(2)-1);

            Box planex_hi = mfi.growntilebox(1); planex_hi.setSmall(0, vbx.bigEnd(0)+1);
            Box planey_hi = mfi.growntilebox(1); planey_hi.setSmall(1, vbx.bigEnd(1)+1);
            Box planez_hi = mfi.growntilebox(1); planez_hi.setSmall(2, vbx.bigEnd(2)+1);

            int i_lo   = vbx.smallEnd(0); int i_hi = vbx.bigEnd(0);
            int j_lo   = vbx.smallEnd(1); int j_hi = vbx.bigEnd(1);
            int k_lo   = vbx.smallEnd(2); int k_hi = vbx.bigEnd(2);

            const Array4<Real>& mu_turb = eddyViscosity.array(mfi);
            const Array4<int>& mask_arr = cc_mask.array(mfi);

            auto domlo = lbound(domain);
            auto domhi = ubound(domain);

            ParallelFor(planex_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int lj = amrex::min(amrex::max(j, domlo.y), domhi.y);
                int lk = amrex::min(amrex::max(k, domlo.z), domhi.z);
                if ((mask_arr(i,j,k) == is_notcovered && mask_arr(i_lo,lj,lk) != is_notcovered) ||
                    (mask_arr(i,j,k) == is_physbnd     && i < domlo.x && impose_phys_bcs)) {
                    for (int n = 0; n < ncomp; n++) {
                        mu_turb(i,j,k,n) = mu_turb(i_lo,lj,lk,n);
                    }
                }
            });
            ParallelFor(planex_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int lj = amrex::min(amrex::max(j, domlo.y), domhi.y);
                int lk = amrex::min(amrex::max(k, domlo.z), domhi.z);
                if ((mask_arr(i,j,k) == is_notcovered && mask_arr(i_hi,lj,lk) != is_notcovered) ||
                    (mask_arr(i,j,k) == is_physbnd     && i > domhi.x && impose_phys_bcs)) {
                    for (int n = 0; n < ncomp; n++) {
                        mu_turb(i,j,k,n) = mu_turb(i_hi,lj,lk,n);
                    }
                }
            });
            ParallelFor(planey_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int lk = amrex::min(amrex::max(k, domlo.z), domhi.z);
                if ((mask_arr(i,j,k) == is_notcovered && mask_arr(i,j_lo,lk) != is_notcovered) ||
                    (mask_arr(i,j,k) == is_physbnd     && j < domlo.y && impose_phys_bcs)) {
                    for (int n = 0; n < ncomp; n++) {
                        mu_turb(i,j,k,n) = mu_turb(i,j_lo,lk,n);
                    }
                }
            });
            ParallelFor(planey_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int lk = amrex::min(amrex::max(k, domlo.z), domhi.z);
                if ((mask_arr(i,j,k) == is_notcovered && mask_arr(i,j_hi,lk) != is_notcovered)||
                    (mask_arr(i,j,k) == is_physbnd     && j > domhi.y && impose_phys_bcs)) {
                    for (int n = 0; n < ncomp; n++) {
                        mu_turb(i,j,k,n) = mu_turb(i,j_hi,lk,n);
                    }
                }
            });
            ParallelFor(planez_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if ((mask_arr(i,j,k) == is_notcovered && mask_arr(i,j,k_lo) != is_notcovered) ||
                    (mask_arr(i,j,k) == is_physbnd     && k < domlo.z && impose_phys_bcs)) {
                    if (mask_arr(i,j,k) == is_physbnd) {
                    }
                    for (int n = 0; n < ncomp; n++) {
                        mu_turb(i,j,k,n) = mu_turb(i,j,k_lo,n);
                    }
                }
            });
            ParallelFor(planez_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if ((mask_arr(i,j,k) == is_notcovered && mask_arr(i,j,k_hi) != is_notcovered) ||
                    (mask_arr(i,j,k) == is_physbnd     && k > domhi.z && impose_phys_bcs)) {
                    for (int n = 0; n < ncomp; n++) {
                        mu_turb(i,j,k,n) = mu_turb(i,j,k_hi,n);
                    }
                }
            });
        } // mfi
        eddyViscosity.FillBoundary(geom.periodicity());
    }
}
