#include "ERF_SurfaceLayer.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"
#include "ERF_TileNoZ.H"

using namespace amrex;

void
ComputeDiffusivityMRF (const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& cons_in,
                       MultiFab& eddyViscosity,
                       const Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       bool use_terrain_fitted_coords,
                       bool /*use_moisture*/,
                       int level,
                       const BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<MultiFab>& z_phys_nd,
                       const MoistureComponentIndices& moisture_indices)
{
    /*
    Implementation of the older MRF Scheme based on Hong and Pan (1996)
    " Nonlocal Boundary Layer Vertical Diffusion in a Medium-Range Forecast
    Model"
    */

    // Domain extent in z-dir
    int klo = geom.Domain().smallEnd(2);
    int khi = geom.Domain().bigEnd(2);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(eddyViscosity, TileNoZ()); mfi.isValid(); ++mfi) {

        // Box operated on must span fill domain in z-dir
        const Box& gbx = mfi.growntilebox(IntVect(1,1,0));
        AMREX_ALWAYS_ASSERT( gbx.smallEnd(2) == klo &&
                             gbx.bigEnd(2)   == khi );

        // Step 1: Compute the height of the PBL without thermal excess
        // h = Rib_cf * theta_va * | U(h) |^2 / (g * (theta_va - theta_surf))
        // create flattened boxes to store PBL height
        const GeometryData gdata = geom.data();
        const Box xybx = PerpendicularBox<ZDir>(gbx, IntVect{0, 0, 0});
        FArrayBox pbl_height_predictor(xybx, 1, The_Async_Arena());
        FArrayBox pbl_height_corrector(xybx, 1, The_Async_Arena());
        IArrayBox pbl_index(xybx, 1, The_Async_Arena());
        const auto& pblh_arr      = pbl_height_predictor.array();
        const auto& pblh_corr_arr = pbl_height_corrector.array();
        const auto& pbli_arr      = pbl_index.array();

        // Get some data in arrays
        const auto& cell_data = cons_in.const_array(mfi);
        const auto& uvel = xvel.const_array(mfi);
        const auto& vvel = yvel.const_array(mfi);

        const Real Ribcr = turbChoice.pbl_mrf_Ribcr;
        //const Real f0 = turbChoice.pbl_mrf_coriolis_freq;
        const auto& u_star_arr = SurfLayer->get_u_star(level)->const_array(mfi);
        const auto& t_star_arr = SurfLayer->get_t_star(level)->const_array(mfi);
        const auto& l_obuk_arr = SurfLayer->get_olen(level)->const_array(mfi);
        const auto& t10av_arr  = SurfLayer->get_mac_avg(level, 2)->const_array(mfi);
        const auto& t_surf_arr = SurfLayer->get_t_surf(level)->const_array(mfi);
        const Array4<Real const> z_nd_arr = z_phys_nd->array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_surf  = t_surf_arr(i, j, 0);
            const Real t_layer = t10av_arr(i, j, 0);

            int kpbl  = klo;
            Real zval = 10;
            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : gdata.ProbLo(2) + (kpbl + 0.5) * gdata.CellSize(2);
                kpbl += 1;

                const Real theta = cell_data(i, j, kpbl, RhoTheta_comp) /
                                   cell_data(i, j, kpbl, Rho_comp);
                const Real ws2 = 0.25 * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real Rib = CONST_GRAV * zval * (theta - t_surf) / (ws2 * t_layer);
                above_critical = (Rib >= Ribcr);
            }
            // Initial PBL Height
            // Avoiding detailed interpolation here and just using map-nearest
            // neighbor Empirical expression for PBLH is given by h = c u* / f
            // Garratt (1994) and Tennekes (1982)
            //const Real c_pblh = (l_obuk_arr(i, j, 0) > 0) ? 0.16 : 0.60;
            //const Real pblh_emp = c_pblh * u_star_arr(i, j, 0) / f0;
            const Real pblh_emp = gdata.ProbLo(2) + 0.5 * gdata.CellSize(2);
            pblh_arr(i, j, 0) = (above_critical) ? zval : pblh_emp;
            pbli_arr(i, j, 0) = (above_critical) ? kpbl : 1;
        });

        // Corrector PBL height for thermal excess
        const Real const_b = turbChoice.pbl_mrf_const_b;
        const Real sf = turbChoice.pbl_mrf_sf;
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
        {
            const Real t_layer  = t10av_arr(i, j, 0);
            const Real phiM     = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -1.0 / 3.0);
            const Real wstar    = u_star_arr(i, j, 0) / phiM;
            const Real t_excess = -const_b * u_star_arr(i, j, 0) * t_star_arr(i, j, 0) / wstar;
            const Real t_surf   = t_layer + std::max(std::min(t_excess, 3.0), 0.0);

            int kpbl  = klo;
            Real zval = 10;
            bool above_critical = false;
            while (!above_critical && ((kpbl + 1) <= khi)) {
                zval = (use_terrain_fitted_coords)
                     ? Compute_Zrel_AtCellCenter(i, j, kpbl, z_nd_arr)
                     : gdata.ProbLo(2) + (kpbl + 0.5) * gdata.CellSize(2);
                kpbl += 1;
                const Real theta = cell_data(i, j, kpbl, RhoTheta_comp) /
                                   cell_data(i, j, kpbl, Rho_comp);
                const Real ws2 = 0.25 * ( (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) *
                                          (uvel(i, j, kpbl) + uvel(i + 1, j, kpbl)) +
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) *
                                          (vvel(i, j, kpbl) + vvel(i, j + 1, kpbl)) );
                const Real Rib = CONST_GRAV * zval * (theta - t_surf) / (ws2 * t_layer);
                above_critical = (Rib >= Ribcr);
            }
            //const Real c_pblh = (l_obuk_arr(i, j, 0) > 0) ? 0.16 : 0.60;
            //const Real pblh_emp = c_pblh * u_star_arr(i, j, 0) / f0;
            const Real pblh_emp = gdata.ProbLo(2) + 0.5 * gdata.CellSize(2);
            pblh_corr_arr(i, j, 0) = (above_critical) ? zval : pblh_emp;
                 pbli_arr(i, j, 0) = (above_critical) ? kpbl : 1;
        });
        /*
          amrex::Print() << "PBL height computed for MRF scheme at level "
          << pblh_arr(2, 2, 0) << "  " << pblh_corr_arr(2, 2, 0)
          << std::endl;
          amrex::Print() << "PBL Temp:" << t_surf_arr(2, 2, 0) << "  "
          << t10av_arr(2, 2, 0) << std::endl;
        */

        // -- Compute diffusion coefficients --

        const Array4<Real>& K_turb = eddyViscosity.array(mfi);

        // Dirichlet flags to switch derivative stencil
        bool c_ext_dir_on_zlo = ((bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir));
        bool c_ext_dir_on_zhi = ((bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir));
        bool u_ext_dir_on_zlo = ((bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir));
        bool u_ext_dir_on_zhi = ((bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir));
        bool v_ext_dir_on_zlo = ((bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir));
        bool v_ext_dir_on_zhi = ((bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir));

        const auto& dxInv = geom.InvCellSizeArray();
        const Real dz_inv = geom.InvCellSize(2);
        const int izmin   = geom.Domain().smallEnd(2);
        const int izmax   = geom.Domain().bigEnd(2);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real zval = (use_terrain_fitted_coords)
                            ? Compute_Zrel_AtCellCenter(i, j, k, z_nd_arr)
                            : gdata.ProbLo(2) + (k + 0.5) * gdata.CellSize(2);
            const Real rho = cell_data(i, j, k, Rho_comp);
            const Real met_h_zeta = (use_terrain_fitted_coords)
                                  ? Compute_h_zeta_AtCellCenter(i, j, k, dxInv, z_nd_arr) : 1.0;
            const Real dz_terrain = met_h_zeta / dz_inv;
            if (k < pbli_arr(i, j, 0)) {
                const Real phiM = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 8 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -1.0 / 3.0);
                const Real phit = (l_obuk_arr(i, j, 0) > 0)
                                ? (1 + 5 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0))
                                : std::pow(
                                           (1 - 16 * sf * pblh_arr(i, j, 0) / l_obuk_arr(i, j, 0)),
                                           -1.0 / 2.0);
                const Real Prt = phit / phiM + const_b * KAPPA * sf;
                const Real wstar = u_star_arr(i, j, 0) / phiM;
                K_turb(i, j, k, EddyDiff::Mom_v) = rho * wstar * KAPPA * zval
                                                 * (1 - zval / pblh_corr_arr(i, j, 0))
                                                 * (1 - zval / pblh_corr_arr(i, j, 0));
                K_turb(i, j, k, EddyDiff::Theta_v) = K_turb(i, j, k, EddyDiff::Mom_v) / Prt;
            } else {
                const Real lambda = 150.0;
                const Real lscale = (KAPPA * zval * lambda) / (KAPPA * zval + lambda);
                Real dthetadz, dudz, dvdz;
                ComputeVerticalDerivativesPBL(i, j, k, uvel, vvel, cell_data, izmin, izmax, 1.0 / dz_terrain,
                                              c_ext_dir_on_zlo, c_ext_dir_on_zhi, u_ext_dir_on_zlo,
                                              u_ext_dir_on_zhi, v_ext_dir_on_zlo, v_ext_dir_on_zhi, dthetadz,
                                              dudz, dvdz, moisture_indices);
                const Real wind_shear = dudz * dudz + dvdz * dvdz + 1.0e-9;
                const Real theta   = cell_data(i, j, k, RhoTheta_comp) / cell_data(i, j, k, Rho_comp);
                const Real grad_Ri = std::max(CONST_GRAV / theta * dthetadz / wind_shear, -100.0);
                /*
                  const Real Pr = 1.5 + 3.08 * grad_Ri;
                  const Real fm =
                  (grad_Ri > 0)
                  ? (std::exp(-8.5 * grad_Ri) + (0.15 / (grad_Ri + 3.0)) * Pr)
                  : std::pow((1 - 12 * grad_Ri), -1.0 / 3.0);
                  const Real ft =
                  (grad_Ri > 0)
                  ? (std::exp(-8.5 * grad_Ri) + (0.15 / (grad_Ri + 3.0)))
                  : std::pow((1 - 16 * grad_Ri), -1.0 / 2.0);
                */
                // Using YSU model instead of MRF model
                const Real Pr = 1.0 + 2.1 * grad_Ri;
                const Real fm = (grad_Ri > 0)
                              ? 1.0 / ((1.0 + 5.0 * grad_Ri) * (1.0 + 5.0 * grad_Ri))
                              : 1 - 8 * grad_Ri / (1 + 1.746 * std::sqrt(-grad_Ri));
                const Real ft = (grad_Ri > 0)
                              ? 1.0 / ((1.0 + 5.0 * grad_Ri) * (1.0 + 5.0 * grad_Ri))
                              : 1 - 8 * grad_Ri / (1 + 1.286 * std::sqrt(-grad_Ri));
                const Real rl2wsp = rho * lscale * lscale * std::sqrt(wind_shear);
                K_turb(i, j, k, EddyDiff::Mom_v)   = rl2wsp * fm * Pr;
                K_turb(i, j, k, EddyDiff::Theta_v) = rl2wsp * ft;
            }

            // limit both diffusion coefficients - from WRF, not documented in
            // papers
            constexpr Real ckz  = 0.001;
            constexpr Real Kmax = 1000.0;
            const Real rhoKmin  = ckz * dz_terrain * rho;
            const Real rhoKmax  = rho * Kmax;
            K_turb(i, j, k, EddyDiff::Mom_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Mom_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Theta_v) = std::max(
                std::min(K_turb(i, j, k, EddyDiff::Theta_v), rhoKmax), rhoKmin);
            K_turb(i, j, k, EddyDiff::Q_v) = K_turb(i, j, k, EddyDiff::Theta_v);
            K_turb(i, j, k, EddyDiff::Turb_lengthscale) = pblh_arr(i, j, 0);
        });

        // FOEXTRAP top and bottom ghost cells
        ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
        {
            K_turb(i, j, klo-1, EddyDiff::Mom_v  ) = K_turb(i, j, klo, EddyDiff::Mom_v  );
            K_turb(i, j, klo-1, EddyDiff::Theta_v) = K_turb(i, j, klo, EddyDiff::Theta_v);
            K_turb(i, j, klo-1, EddyDiff::Q_v    ) = K_turb(i, j, klo, EddyDiff::Q_v    );
            K_turb(i, j, khi+1, EddyDiff::Mom_v  ) = K_turb(i, j, khi, EddyDiff::Mom_v  );
            K_turb(i, j, khi+1, EddyDiff::Theta_v) = K_turb(i, j, khi, EddyDiff::Theta_v);
            K_turb(i, j, khi+1, EddyDiff::Q_v    ) = K_turb(i, j, khi, EddyDiff::Q_v    );
        });
    }// mfi
}
