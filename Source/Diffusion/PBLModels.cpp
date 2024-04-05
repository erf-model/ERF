#include "ABLMost.H"
#include "DirectionSelector.H"
#include "Diffusion.H"
#include "ERF_Constants.H"
#include "TurbStruct.H"
#include "PBLModels.H"

using namespace amrex;

/**
 * Function to compute turbulent viscosity with PBL.
 *
 * @param[in] xvel velocity in x-dir
 * @param[in] yvel velocity in y-dir
 * @param[in] cons_in cell center conserved quantities
 * @param[out] eddyViscosity holds turbulent viscosity
 * @param[in] geom problem geometry
 * @param[in] turbChoice container with turbulence parameters
 * @param[in] most pointer to Monin-Obukhov class if instantiated
 */
void
ComputeTurbulentViscosityPBL (const MultiFab& xvel,
                              const MultiFab& yvel,
                              const MultiFab& cons_in,
                              MultiFab& eddyViscosity,
                              const Geometry& geom,
                              const TurbChoice& turbChoice,
                              std::unique_ptr<ABLMost>& most,
                              int level,
                              const BCRec* bc_ptr,
                              bool /*vert_only*/,
                              const std::unique_ptr<MultiFab>& z_phys_nd)
{
    const bool use_terrain = (z_phys_nd != nullptr);

    // MYNN Level 2.5 PBL Model
    if (turbChoice.pbl_type == PBLType::MYNN25) {

        const Real A1 = turbChoice.pbl_mynn_A1;
        const Real A2 = turbChoice.pbl_mynn_A2;
        const Real B1 = turbChoice.pbl_mynn_B1;
        const Real B2 = turbChoice.pbl_mynn_B2;
        const Real C1 = turbChoice.pbl_mynn_C1;
        const Real C2 = turbChoice.pbl_mynn_C2;
        const Real C3 = turbChoice.pbl_mynn_C3;
        //const Real C4 = turbChoice.pbl_mynn_C4;
        const Real C5 = turbChoice.pbl_mynn_C5;

        // Dirichlet flags to switch derivative stencil
        bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
        bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir) );
        bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
        bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir) );
        bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
        bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir) );

        // Epsilon
        Real eps = std::numeric_limits<Real>::epsilon();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(eddyViscosity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box &bx = mfi.growntilebox(1);
            const Array4<Real const> &cell_data     = cons_in.array(mfi);
            const Array4<Real      > &K_turb = eddyViscosity.array(mfi);
            const Array4<Real const> &uvel = xvel.array(mfi);
            const Array4<Real const> &vvel = yvel.array(mfi);

            // Compute some quantities that are constant in each column
            // Sbox is shrunk to only include the interior of the domain in the vertical direction to compute integrals
            // Box includes one ghost cell in each direction
            const Box &dbx = geom.Domain();
            Box sbx(bx.smallEnd(), bx.bigEnd());
            sbx.grow(2,-1);
            AMREX_ALWAYS_ASSERT(sbx.smallEnd(2) == dbx.smallEnd(2) && sbx.bigEnd(2) == dbx.bigEnd(2));

            const GeometryData gdata = geom.data();

            const Box xybx = PerpendicularBox<ZDir>(bx, IntVect{0,0,0});
            FArrayBox qintegral(xybx,2);
            qintegral.setVal<RunOn::Device>(0.0);
            FArrayBox qturb(bx,1); FArrayBox qturb_old(bx,1);
            const Array4<Real> qint = qintegral.array();
            const Array4<Real> qvel = qturb.array();
            const Array4<Real> qvel_old = qturb_old.array();

            // vertical integrals to compute lengthscale
            if (use_terrain) {
                const Array4<Real const> &z_nd_arr = z_phys_nd->array(mfi);
                const auto invCellSize = geom.InvCellSizeArray();
                ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
                {
                    qvel(i,j,k)     = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
                    qvel_old(i,j,k) = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp) + eps);
                    AMREX_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "QKE must have a positive value");
                    AMREX_ASSERT_WITH_MESSAGE(qvel_old(i,j,k) > 0.0, "Old QKE must have a positive value");

                    if (sbx.contains(i,j,k)) {
                        const Real Zval = Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr);
                        const Real dz = Compute_h_zeta_AtCellCenter(i,j,k,invCellSize,z_nd_arr);
                        Gpu::deviceReduceSum(&qint(i,j,0,0), Zval*qvel(i,j,k)*dz, handler);
                        Gpu::deviceReduceSum(&qint(i,j,0,1), qvel(i,j,k)*dz, handler);
                    }
                });
            } else {
                ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
                                   [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
                {
                    qvel(i,j,k)     = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
                    qvel_old(i,j,k) = std::sqrt(cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp) + eps);
                    AMREX_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "QKE must have a positive value");
                    AMREX_ASSERT_WITH_MESSAGE(qvel_old(i,j,k) > 0.0, "Old QKE must have a positive value");

                    // Not multiplying by dz: its constant and would fall out when we divide qint0/qint1 anyway
                    if (sbx.contains(i,j,k)) {
                        const Real Zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                        Gpu::deviceReduceSum(&qint(i,j,0,0), Zval*qvel(i,j,k), handler);
                        Gpu::deviceReduceSum(&qint(i,j,0,1), qvel(i,j,k), handler);
                    }
                });
            }

            Real dz_inv = geom.InvCellSize(2);
            const auto& dxInv = geom.InvCellSizeArray();
            int izmin = geom.Domain().smallEnd(2);
            int izmax = geom.Domain().bigEnd(2);

            // Spatially varying MOST
            Real d_kappa   = KAPPA;
            Real d_gravity = CONST_GRAV;

            const auto& t_mean_mf = most->get_mac_avg(level,2); // TODO: IS THIS ACTUALLY RHOTHETA
            const auto& u_star_mf = most->get_u_star(level);    // Use desired level
            const auto& t_star_mf = most->get_t_star(level);    // Use desired level

            const auto& tm_arr     = t_mean_mf->array(mfi); // TODO: IS THIS ACTUALLY RHOTHETA
            const auto& u_star_arr = u_star_mf->array(mfi);
            const auto& t_star_arr = t_star_mf->array(mfi);

            const Array4<Real const> z_nd_arr = use_terrain ? z_phys_nd->array(mfi) : Array4<Real>{};

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Compute some partial derivatives that we will need (second order)
                // U and V derivatives are interpolated to account for staggered grid
                const Real met_h_zeta = use_terrain ? Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd_arr) : 1.0;
                Real dthetadz, dudz, dvdz;
                ComputeVerticalDerivativesPBL(i, j, k,
                                              uvel, vvel, cell_data, izmin, izmax, dz_inv/met_h_zeta,
                                              c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                              u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                              v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                              dthetadz, dudz, dvdz);

                // Spatially varying MOST
                Real surface_heat_flux = -u_star_arr(i,j,0) * t_star_arr(i,j,0);
                Real theta0            = tm_arr(i,j,0); // TODO: IS THIS ACTUALLY RHOTHETA
                Real l_obukhov;
                if (std::abs(surface_heat_flux) > eps) {
                    l_obukhov = ( theta0 * u_star_arr(i,j,0) * u_star_arr(i,j,0) ) /
                        ( d_kappa * d_gravity * t_star_arr(i,j,0) );
                } else {
                    l_obukhov = std::numeric_limits<Real>::max();
                }

                // First Length Scale
                AMREX_ASSERT(l_obukhov != 0);
                int lk = amrex::max(k,0);
                const Real zval = use_terrain ? Compute_Zrel_AtCellCenter(i,j,lk,z_nd_arr) : gdata.ProbLo(2) + (lk + 0.5)*gdata.CellSize(2);
                const Real zeta = zval/l_obukhov;
                Real l_S;
                if (zeta >= 1.0) {
                    l_S = KAPPA*zval/3.7;
                } else if (zeta >= 0) {
                    l_S = KAPPA*zval/(1+2.7*zeta);
                } else {
                    l_S = KAPPA*zval*std::pow(1.0 - 100.0 * zeta, 0.2);
                }

                // Second Length Scale
                Real l_T;
                if (qint(i,j,0,1) > 0.0) {
                    l_T = 0.23*qint(i,j,0,0)/qint(i,j,0,1);
                } else {
                    l_T = std::numeric_limits<Real>::max();
                }

                // Third Length Scale
                Real l_B;
                if (dthetadz > 0) {
                    Real N_brunt_vaisala = std::sqrt(CONST_GRAV/theta0 * dthetadz);
                    if (zeta < 0) {
                        Real qc = CONST_GRAV/theta0 * surface_heat_flux * l_T;
                        qc = std::pow(qc,1.0/3.0);
                        l_B = (1.0 + 5.0*std::sqrt(qc/(N_brunt_vaisala * l_T))) * qvel(i,j,k)/N_brunt_vaisala;
                    } else {
                        l_B = qvel(i,j,k) / N_brunt_vaisala;
                    }
                } else {
                    l_B = std::numeric_limits<Real>::max();
                }

                // Overall Length Scale
                Real l_comb = 1.0 / (1.0/l_S + 1.0/l_T + 1.0/l_B);

                // NOTE: Level 2 limiting from balance of production and dissipation.
                //       K_turb has a setval of 0.0 when the MF is created (NOT EACH STEP).
                //       We do this inline to avoid storing qe^2 at each cell.
                Real l_comb_old   = K_turb(i,j,k,EddyDiff::PBL_lengthscale);
                Real shearProd    = dudz*dudz + dvdz*dvdz;
                Real buoyProd     = -(CONST_GRAV/theta0) * dthetadz;
                Real lSM          = K_turb(i,j,k,EddyDiff::Mom_v)   / (qvel_old(i,j,k) + eps);
                Real lSH          = K_turb(i,j,k,EddyDiff::Theta_v) / (qvel_old(i,j,k) + eps);
                Real qe2          = B1 * l_comb_old * ( lSM * shearProd + lSH * buoyProd );
                Real qe           = (qe2 < 0.0) ? 0.0 : std::sqrt(qe2);
                Real one_m_alpha  = (qvel(i,j,k) > qe) ? 1.0 : qvel(i,j,k) / (qe + eps);
                Real one_m_alpha2 = one_m_alpha * one_m_alpha;

                // Compute non-dimensional parameters
                Real l2_over_q2   = l_comb*l_comb/(qvel(i,j,k)*qvel(i,j,k));
                Real GM = l2_over_q2 * shearProd;
                Real GH = l2_over_q2 * buoyProd;
                Real E1 = 1.0 + one_m_alpha2 * ( 6.0*A1*A1*GM - 9.0*A1*A2*(1.0-C2)*GH );
                Real E2 = one_m_alpha2 * ( -3.0*A1*(4.0*A1 + 3.0*A2*(1.0-C5))*(1.0-C2)*GH );
                Real E3 = one_m_alpha2 * ( 6.0*A1*A2*GM );
                Real E4 = 1.0 + one_m_alpha2 * ( -12.0*A1*A2*(1.0-C2)*GH - 3.0*A2*B2*(1.0-C3)*GH );
                Real R1 = one_m_alpha * ( A1*(1.0-3.0*C1) );
                Real R2 = one_m_alpha * A2;

                Real SM = (R2*E2 - R1*E4)/(E2*E3 - E1*E4);
                Real SH = (R1*E3 - R2*E1)/(E2*E3 - E1*E4);
                Real SQ = 3.0 * SM; // Nakanishi & Niino 2009

                // Finally, compute the eddy viscosity/diffusivities
                const Real rho = cell_data(i,j,k,Rho_comp);
                K_turb(i,j,k,EddyDiff::Mom_v)   = rho * l_comb * qvel(i,j,k) * SM * 0.5; // 0.5 for mu_turb
                K_turb(i,j,k,EddyDiff::Theta_v) = rho * l_comb * qvel(i,j,k) * SH;
                K_turb(i,j,k,EddyDiff::QKE_v)   = rho * l_comb * qvel(i,j,k) * SQ;

                K_turb(i,j,k,EddyDiff::PBL_lengthscale) = l_comb;
                // TODO: How should this be done for other components (scalars, moisture)
            });
        }
    } else if (turbChoice.pbl_type == PBLType::YSU) {
        // FIXME
        // Error("YSU Model not implemented yet");

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
            //const auto& K_turb = eddyViscosity.array(mfi);
            const auto& uvel = xvel.const_array(mfi);
            const auto& vvel = yvel.const_array(mfi);

            const auto& z0_arr = most->get_z0(level)->const_array();
            const auto& ws10av_arr = most->get_mac_avg(level,4)->const_array(mfi);
            const auto& t10av_arr  = most->get_mac_avg(level,2)->const_array(mfi);
            //const auto& u_star_arr = most->get_u_star(level)->const_array(mfi);
            //const auto& t_star_arr = most->get_t_star(level)->const_array(mfi);
            const auto& t_surf_arr = most->get_t_surf(level)->const_array(mfi);
            //const auto& l_obuk_arr = most->get_olen(level)->const_array(mfi);
            const bool use_terrain = (z_phys_nd != nullptr);
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
            const auto& pblh_arr = pbl_height.array();

            // -- Diagnose PBL height - starting out assuming non-moist
            // loop is only over i,j in order to find height at each x,y
            const Real f0 = turbChoice.pbl_ysu_coriolis_freq;
            const bool over_land = turbChoice.pbl_ysu_over_land; // TODO: make this local and consistent
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
                if (over_land) {
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
                    const Real ws2_level = (uvel(i,j,kpbl)+uvel(i+1,j,kpbl))*(uvel(i,j,kpbl)+uvel(i+1,j,kpbl)) + (vvel(i,j,kpbl)+vvel(i,j+1,kpbl))*(uvel(i,j,kpbl)+uvel(i,j+1,kpbl));
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
                if (pblh_arr(i,j,0) < 0.5*(zval_up+zval_dn) ) {
                    kpbl = 0;
                }

                // FIXME: DEBUG
                std::cout << i << " " << j << " " << pblh_arr(i,j,0) << std::endl;
            });

            // -- Compute nonlocal/countergradient mixing parameters

            // -- Compute entrainment parameters

            // -- Compute diffusion coefficients within PBL

            // -- Compute coefficients in free stream above PBL

            // FIXME: DEBUG
            Error("YSU Model not implemented yet");
        }

    }
}
