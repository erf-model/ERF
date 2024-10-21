#include "ERF_ABLMost.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"

using namespace amrex;

void
ComputeDiffusivityYSU (const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& cons_in,
                       MultiFab& eddyViscosity,
                       const Geometry& geom,
                       const TurbChoice& turbChoice,
                       std::unique_ptr<ABLMost>& most,
                       bool /*use_moisture*/,
                       int level,
                       const BCRec* bc_ptr,
                       bool /*vert_only*/,
                       const std::unique_ptr<MultiFab>& z_phys_nd)
{
    const bool use_terrain = (z_phys_nd != nullptr);

    {
        /*
          YSU PBL initially introduced by S.-Y. Hong, Y. Noh, and J. Dudhia, MWR, 2006 [HND06]

          Further Modifications from S.-Y. Hong, Q. J. R. Meteorol. Soc., 2010 [H10]

          Implementation follows WRF as of early 2024 with some simplifications
        */

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            // Pull out the box we're working on, make sure it covers full domain in z-direction
            const Box &bx = mfi.growntilebox(1);
            const Box &dbx = geom.Domain();
            Box sbx(bx.smallEnd(), bx.bigEnd());
            sbx.grow(2,-1);
            AMREX_ALWAYS_ASSERT(sbx.smallEnd(2) == dbx.smallEnd(2) && sbx.bigEnd(2) == dbx.bigEnd(2));

            // Get some data in arrays
            const auto& cell_data = cons_in.const_array(mfi);
            const auto& uvel = xvel.const_array(mfi);
            const auto& vvel = yvel.const_array(mfi);

            const auto& z0_arr = most->get_z0(level)->const_array();
            const auto& ws10av_arr = most->get_mac_avg(level,5)->const_array(mfi);
            const auto& t10av_arr  = most->get_mac_avg(level,2)->const_array(mfi);
            const auto& t_surf_arr = most->get_t_surf(level)->const_array(mfi);
            const auto& over_land_arr = (most->get_lmask(level)) ? most->get_lmask(level)->const_array(mfi) : Array4<int> {};
            const Array4<Real const> z_nd_arr = use_terrain ? z_phys_nd->array(mfi) : Array4<Real>{};
            const Real most_zref = most->get_zref();

            // Require that MOST zref is 10 m so we get the wind speed at 10 m from most
            bool invalid_zref = false;
            if (use_terrain) {
                invalid_zref = most_zref != 10.0;
            } else {
                // zref gets reset to nearest cell center, so assert that zref is in the same cell as the 10m point
                Real dz = geom.CellSize(2);
                invalid_zref = int((most_zref - 0.5*dz)/dz) != int((10.0 - 0.5*dz)/dz);
            }
            if (invalid_zref) {
                Print() << "most_zref = " << most_zref << std::endl;
                Abort("MOST Zref must be 10m for YSU PBL scheme");
            }

            // create flattened boxes to store PBL height
            const GeometryData gdata = geom.data();
            const Box xybx = PerpendicularBox<ZDir>(bx, IntVect{0,0,0});
            FArrayBox pbl_height(xybx,1);
            IArrayBox pbl_index(xybx,1);
            const auto& pblh_arr = pbl_height.array();
            const auto& pbli_arr = pbl_index.array();

            // -- Diagnose PBL height - starting out assuming non-moist --
            // loop is only over i,j in order to find height at each x,y
            const Real f0 = turbChoice.pbl_ysu_coriolis_freq;
            const bool force_over_water = turbChoice.pbl_ysu_force_over_water;
            const Real land_Ribcr = turbChoice.pbl_ysu_land_Ribcr;
            const Real unst_Ribcr = turbChoice.pbl_ysu_unst_Ribcr;
            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                // Reconstruct a surface bulk Richardson number from the surface layer model
                // In WRF, this value is supplied to YSU by the MM5 surface layer model
                const Real t_surf = t_surf_arr(i,j,0);
                const Real t_layer = t10av_arr(i,j,0);
                const Real ws_layer = ws10av_arr(i,j,0);
                const Real Rib_layer = CONST_GRAV * most_zref / (ws_layer*ws_layer) * (t_layer - t_surf)/(t_layer);

                // For now, we only support stable boundary layers
                if (Rib_layer < unst_Ribcr) {
                    Abort("For now, YSU PBL only supports stable conditions");
                }

                // TODO: unstable BLs

                // PBL Height: Stable Conditions
                Real Rib_cr;
                bool over_land = (over_land_arr) ? over_land_arr(i,j,0) : 1;
                if (over_land && !force_over_water) {
                    Rib_cr = land_Ribcr;
                } else { // over water
                    // Velocity at z=10 m comes from MOST -> currently the average using whatever averaging MOST uses.
                    // TODO: Revisit this calculation with local ws10?
                    const Real z0 = z0_arr(i,j,0);
                    const Real Rossby = ws_layer/(f0*z0);
                    Rib_cr = min(0.16*std::pow(1.0e-7*Rossby,-0.18),0.3); // Note: upper bound in WRF code, but not H10 paper
                }

                bool above_critical = false;
                int kpbl = 0;
                Real Rib_up = Rib_layer, Rib_dn;
                const Real base_theta = cell_data(i,j,0,RhoTheta_comp) / cell_data(i,j,0,Rho_comp);
                while (!above_critical and bx.contains(i,j,kpbl+1)) {
                    kpbl += 1;
                    const Real zval = use_terrain ? Compute_Zrel_AtCellCenter(i,j,kpbl,z_nd_arr) : gdata.ProbLo(2) + (kpbl + 0.5)*gdata.CellSize(2);
                    const Real ws2_level = 0.25*( (uvel(i,j,kpbl)+uvel(i+1,j  ,kpbl))*(uvel(i,j,kpbl)+uvel(i+1,j  ,kpbl))
                                                + (vvel(i,j,kpbl)+vvel(i  ,j+1,kpbl))*(vvel(i,j,kpbl)+vvel(i  ,j+1,kpbl)) );
                    const Real theta = cell_data(i,j,kpbl,RhoTheta_comp) / cell_data(i,j,kpbl,Rho_comp);
                    Rib_dn = Rib_up;
                    Rib_up = (theta-base_theta)/base_theta * CONST_GRAV * zval / ws2_level;
                    above_critical = Rib_up >= Rib_cr;
                }

                Real interp_fact;
                if (Rib_dn >= Rib_cr) {
                    interp_fact = 0.0;
                } else if (Rib_up <= Rib_cr)
                    interp_fact = 1.0;
                else {
                    interp_fact = (Rib_cr - Rib_dn) / (Rib_up - Rib_dn);
                }

                const Real zval_up = use_terrain ? Compute_Zrel_AtCellCenter(i,j,kpbl,z_nd_arr) : gdata.ProbLo(2) + (kpbl + 0.5)*gdata.CellSize(2);
                const Real zval_dn = use_terrain ? Compute_Zrel_AtCellCenter(i,j,kpbl-1,z_nd_arr) : gdata.ProbLo(2) + (kpbl-1 + 0.5)*gdata.CellSize(2);
                pblh_arr(i,j,0) = zval_dn + interp_fact*(zval_up-zval_dn);

                const Real zval_0 = use_terrain ? Compute_Zrel_AtCellCenter(i,j,0,z_nd_arr) : gdata.ProbLo(2) + (0.5)*gdata.CellSize(2);
                const Real zval_1 = use_terrain ? Compute_Zrel_AtCellCenter(i,j,1,z_nd_arr) : gdata.ProbLo(2) + (1.5)*gdata.CellSize(2);
                if (pblh_arr(i,j,0) < 0.5*(zval_0+zval_1) ) {
                    kpbl = 0;
                }
                pbli_arr(i,j,0) = kpbl;
            });

            // -- Compute nonlocal/countergradient mixing parameters --
            // Not included for stable so nothing to do until unstable treatment is added

            // -- Compute entrainment parameters --
            // 0 for stable so nothing to do?

            // -- Compute diffusion coefficients --

            const auto& u_star_arr = most->get_u_star(level)->const_array(mfi);
            const auto& l_obuk_arr = most->get_olen(level)->const_array(mfi);
            const Array4<Real      > &K_turb = eddyViscosity.array(mfi);

            // Dirichlet flags to switch derivative stencil
            bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
            bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir) );

            const auto& dxInv = geom.InvCellSizeArray();
            const Real dz_inv = geom.InvCellSize(2);
            const int izmin = geom.Domain().smallEnd(2);
            const int izmax = geom.Domain().bigEnd(2);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const Real zval = use_terrain ? Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr) : gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                const Real rho = cell_data(i,j,k,Rho_comp);
                const Real met_h_zeta = use_terrain ? Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd_arr) : 1.0;
                const Real dz_terrain = met_h_zeta/dz_inv;
                if (k < pbli_arr(i,j,0)) {
                    // -- Compute diffusion coefficients within PBL
                    constexpr Real zfacmin = 1e-8; // value from WRF
                    constexpr Real phifac = 8.0; // value from H10 and WRF
                    constexpr Real wstar3 = 0.0; // only nonzero for unstable
                    constexpr Real pfac = 2.0; // profile exponent
                    const Real zfac = std::min(std::max(1 - zval / pblh_arr(i,j,0), zfacmin ), 1.0);
                    // Not including YSU top down PBL term (not in H10, added to WRF later)
                    const Real ust3 = u_star_arr(i,j,0) * u_star_arr(i,j,0) * u_star_arr(i,j,0);
                    Real wscalek = ust3 + phifac * KAPPA * wstar3 * (1.0 - zfac);
                    wscalek = std::pow(wscalek, 1.0/3.0);
                    // stable only
                    const Real phi_term = 1 + 5 * zval / l_obuk_arr(i,j,0); // phi_term appears in WRF but not papers
                    wscalek = std::max(u_star_arr(i,j,0) / phi_term, 0.001); // 0.001 limit appears in WRF but not papers
                    K_turb(i,j,k,EddyDiff::Mom_v) = rho * wscalek * KAPPA * zval * std::pow(zfac, pfac);
                    K_turb(i,j,k,EddyDiff::Theta_v) = K_turb(i,j,k,EddyDiff::Mom_v);
                } else {
                    // -- Compute coefficients in free stream above PBL
                    constexpr Real lam0 = 30.0;
                    constexpr Real min_richardson = -100.0;
                    constexpr Real prandtl_max = 4.0;
                    Real dthetadz, dudz, dvdz;
                    const int RhoQv_comp = -1;
                    const int RhoQr_comp = -1;
                    ComputeVerticalDerivativesPBL(i, j, k,
                                                  uvel, vvel, cell_data, izmin, izmax, 1.0/dz_terrain,
                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                                  u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                                  v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                                  dthetadz, dudz, dvdz, RhoQv_comp, RhoQr_comp);
                    const Real shear_squared = dudz*dudz + dvdz*dvdz + 1.0e-9; // 1.0e-9 from WRF to avoid divide by zero
                    const Real theta = cell_data(i,j,k,RhoTheta_comp) / cell_data(i,j,k,Rho_comp);
                    Real richardson = CONST_GRAV / theta * dthetadz / shear_squared;
                    const Real lambdadz = std::min(std::max(0.1*dz_terrain , lam0), 300.0); // in WRF, H10 paper just says use lam0
                    const Real lengthscale = lambdadz * KAPPA * zval / (lambdadz + KAPPA * zval);
                    const Real turbfact = lengthscale * lengthscale * std::sqrt(shear_squared);

                    if (richardson < 0) {
                        richardson = max(richardson, min_richardson);
                        Real sqrt_richardson = std::sqrt(-richardson);
                        K_turb(i,j,k,EddyDiff::Mom_v) = rho * turbfact * (1.0 - 8.0 * richardson / (1.0 + 1.746 * sqrt_richardson));
                        K_turb(i,j,k,EddyDiff::Theta_v) = rho * turbfact * (1.0 - 8.0 * richardson / (1.0 + 1.286 * sqrt_richardson));
                    } else {
                        const Real oneplus5ri = 1.0 + 5.0 * richardson;
                        K_turb(i,j,k,EddyDiff::Theta_v) = rho * turbfact / (oneplus5ri * oneplus5ri);
                        const Real prandtl = std::min(1.0+2.1*richardson, prandtl_max); // limit from WRF
                        K_turb(i,j,k,EddyDiff::Mom_v) = K_turb(i,j,k,EddyDiff::Theta_v) * prandtl;
                    }
                }

                // limit both diffusion coefficients - from WRF, not documented in papers
                constexpr Real ckz = 0.001;
                constexpr Real Kmax = 1000.0;
                const Real rhoKmin = ckz * dz_terrain * rho;
                const Real rhoKmax = rho * Kmax;
                K_turb(i,j,k,EddyDiff::Mom_v) = std::max(std::min(K_turb(i,j,k,EddyDiff::Mom_v) ,rhoKmax), rhoKmin);
                K_turb(i,j,k,EddyDiff::Theta_v) = std::max(std::min(K_turb(i,j,k,EddyDiff::Theta_v) ,rhoKmax), rhoKmin);
                K_turb(i,j,k,EddyDiff::PBL_lengthscale) = pblh_arr(i,j,0);

                // Entrainment factors: placeholder for now
                K_turb(i,j,k,EddyDiff::Mom_ent_YSU) = 0.0;
                K_turb(i,j,k,EddyDiff::Theta_ent_YSU) = 0.0;

            });

            // HACK set bottom ghost cell to 1st cell
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k==-1) {
                    K_turb(i,j,k,EddyDiff::Mom_v) = K_turb(i,j,0,EddyDiff::Mom_v);
                    K_turb(i,j,k,EddyDiff::Theta_v) = K_turb(i,j,0,EddyDiff::Theta_v);
                }
            });
        }
    }
}
