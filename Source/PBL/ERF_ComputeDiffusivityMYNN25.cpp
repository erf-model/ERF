#include "ERF_ABLMost.H"
#include "ERF_DirectionSelector.H"
#include "ERF_Diffusion.H"
#include "ERF_Constants.H"
#include "ERF_TurbStruct.H"
#include "ERF_PBLModels.H"

using namespace amrex;

void
ComputeDiffusivityMYNN25 (const MultiFab& xvel,
                          const MultiFab& yvel,
                          const MultiFab& cons_in,
                          MultiFab& eddyViscosity,
                          const Geometry& geom,
                          const TurbChoice& turbChoice,
                          std::unique_ptr<ABLMost>& most,
                          bool use_moisture,
                          int level,
                          const BCRec* bc_ptr,
                          bool /*vert_only*/,
                          const std::unique_ptr<MultiFab>& z_phys_nd,
                          const int RhoQv_comp,
                          const int RhoQr_comp)
{
    const bool use_terrain = (z_phys_nd != nullptr);

    auto mynn     = turbChoice.pbl_mynn;
    auto level2   = turbChoice.pbl_mynn_level2;

    Real Lt_alpha = (mynn.config == MYNNConfigType::CHEN2021) ? 0.1 : 0.23;

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
        const Array4<Real const> &cell_data = cons_in.array(mfi);
        const Array4<Real      > &K_turb    = eddyViscosity.array(mfi);
        const Array4<Real const> &uvel      = xvel.array(mfi);
        const Array4<Real const> &vvel      = yvel.array(mfi);

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
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                qvel(i,j,k)     = std::sqrt(2.0*cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
                qvel_old(i,j,k) = std::sqrt(2.0*cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp) + eps);
                AMREX_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "QKE must have a positive value");
                AMREX_ASSERT_WITH_MESSAGE(qvel_old(i,j,k) > 0.0, "Old QKE must have a positive value");

                // NOTE: This factor is to avoid an if statement that will break
                //       the devicereducesum since all threads won't participate.
                //       This more performant than Gpu::Atomic::Add.
                Real fac = (sbx.contains(i,j,k)) ? 1.0 : 0.0;
                const Real Zval = Compute_Zrel_AtCellCenter(i,j,k,z_nd_arr);
                const Real dz = Compute_h_zeta_AtCellCenter(i,j,k,invCellSize,z_nd_arr);
                Gpu::Atomic::Add(&qint(i,j,0,0), Zval*qvel(i,j,k)*dz*fac);
                Gpu::Atomic::Add(&qint(i,j,0,1),      qvel(i,j,k)*dz*fac);

            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                qvel(i,j,k)     = std::sqrt(2.0*cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp));
                qvel_old(i,j,k) = std::sqrt(2.0*cell_data(i,j,k,RhoQKE_comp) / cell_data(i,j,k,Rho_comp) + eps);
                AMREX_ASSERT_WITH_MESSAGE(qvel(i,j,k) > 0.0, "QKE must have a positive value");
                AMREX_ASSERT_WITH_MESSAGE(qvel_old(i,j,k) > 0.0, "Old QKE must have a positive value");

                // Not multiplying by dz: its constant and would fall out when we divide qint0/qint1 anyway

                // NOTE: This factor is to avoid an if statement that will break
                //       the devicereducesum since all threads won't participate.
                //       This more performant than Gpu::Atomic::Add.
                Real fac = (sbx.contains(i,j,k)) ? 1.0 : 0.0;
                const Real Zval = gdata.ProbLo(2) + (k + 0.5)*gdata.CellSize(2);
                Gpu::Atomic::Add(&qint(i,j,0,0), Zval*qvel(i,j,k)*fac);
                Gpu::Atomic::Add(&qint(i,j,0,1),      qvel(i,j,k)*fac);
            });
        }

        Real dz_inv = geom.InvCellSize(2);
        const auto& dxInv = geom.InvCellSizeArray();
        int izmin = geom.Domain().smallEnd(2);
        int izmax = geom.Domain().bigEnd(2);

        // Spatially varying MOST
        Real d_kappa   = KAPPA;
        Real d_gravity = CONST_GRAV;

        const auto& t_mean_mf = most->get_mac_avg(level,4); // theta_v
        const auto& q_mean_mf = most->get_mac_avg(level,3); // q_v
        const auto& u_star_mf = most->get_u_star(level);
        const auto& t_star_mf = most->get_t_star(level);
        const auto& q_star_mf = most->get_q_star(level);

        const auto& tm_arr     = t_mean_mf->const_array(mfi);
        const auto& qm_arr     = q_mean_mf->const_array(mfi);
        const auto& u_star_arr = u_star_mf->const_array(mfi);
        const auto& t_star_arr = t_star_mf->const_array(mfi);
        const auto& q_star_arr = (use_moisture) ? q_star_mf->const_array(mfi) : Array4<Real>{};

        const Array4<Real const> z_nd_arr = use_terrain ? z_phys_nd->const_array(mfi) : Array4<Real>{};

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
                                          dthetadz, dudz, dvdz,
                                          RhoQv_comp, RhoQr_comp);

            // Spatially varying MOST
            Real theta0 = tm_arr(i,j,0);
            Real qv0    = qm_arr(i,j,0);
            Real surface_heat_flux = -u_star_arr(i,j,0) * t_star_arr(i,j,0);
            Real surface_latent_heat{0};
            if (use_moisture) {
                // Compute buoyancy flux (Stull Eqn. 4.4.5d)
                surface_latent_heat = -u_star_arr(i,j,0) * q_star_arr(i,j,0);
                surface_heat_flux *= (1.0 + 0.61*qv0);
                surface_heat_flux += 0.61 * theta0 * surface_latent_heat;
            }

            Real l_obukhov;
            if (std::abs(surface_heat_flux) > eps) {
                l_obukhov = -( theta0 * u_star_arr(i,j,0)*u_star_arr(i,j,0)*u_star_arr(i,j,0) )
                           / ( d_kappa * d_gravity * surface_heat_flux );
            } else {
                l_obukhov = std::numeric_limits<Real>::max();
            }

            // Surface-layer length scale (NN09, Eqn. 53)
            AMREX_ASSERT(l_obukhov != 0);
            int lk = amrex::max(k,0);
            const Real zval = use_terrain ? Compute_Zrel_AtCellCenter(i,j,lk,z_nd_arr)
                                          : gdata.ProbLo(2) + (lk + 0.5)*gdata.CellSize(2);
            const Real zeta = zval/l_obukhov;
            Real l_S;
            if (zeta >= 1.0) {
                l_S = KAPPA*zval/3.7;
            } else if (zeta >= 0) {
                l_S = KAPPA*zval/(1+2.7*zeta);
            } else {
                l_S = KAPPA*zval*std::pow(1.0 - 100.0 * zeta, 0.2);
            }

            // ABL-depth length scale (NN09, Eqn. 54)
            Real l_T;
            if (qint(i,j,0,1) > 0.0) {
                l_T = Lt_alpha*qint(i,j,0,0)/qint(i,j,0,1);
            } else {
                l_T = std::numeric_limits<Real>::max();
            }

            // Buoyancy length scale (NN09, Eqn. 55)
            Real l_B;
            if (dthetadz > 0) {
                Real N_brunt_vaisala = std::sqrt(CONST_GRAV/theta0 * dthetadz);
                if (zeta < 0) {
                    Real qc = CONST_GRAV/theta0 * surface_heat_flux * l_T; // velocity scale
                    qc = std::pow(qc,1.0/3.0);
                    l_B = (1.0 + 5.0*std::sqrt(qc/(N_brunt_vaisala * l_T))) * qvel(i,j,k)/N_brunt_vaisala;
                } else {
                    l_B = qvel(i,j,k) / N_brunt_vaisala;
                }
            } else {
                l_B = std::numeric_limits<Real>::max();
            }

            // Master length scale
            Real Lm;
            if (mynn.config == MYNNConfigType::CHEN2021) {
                Lm = std::pow(1.0/(l_S*l_S) + 1.0/(l_T*l_T) + 1.0/(l_B*l_B), -0.5);
            } else {
                // NN09, Eqn 52
                Lm = 1.0 / (1.0/l_S + 1.0/l_T + 1.0/l_B);
            }

            // Calculate nondimensional production terms
            Real shearProd  = dudz*dudz + dvdz*dvdz;
            Real buoyProd   = -(CONST_GRAV/theta0) * dthetadz;
            Real L2_over_q2 = Lm*Lm/(qvel(i,j,k)*qvel(i,j,k));
            Real GM         = L2_over_q2 * shearProd;
            Real GH         = L2_over_q2 * buoyProd;

            // Equilibrium (Level-2) q calculation follows NN09, Appendix 2
            Real Rf  = level2.calc_Rf(GM, GH);
            Real SM2 = level2.calc_SM(Rf);
            Real qe2 = mynn.B1*Lm*Lm*SM2*(1.0-Rf)*shearProd;
            Real qe  = (qe2 < 0.0) ? 0.0 : std::sqrt(qe2);

            // Level 2 limiting (Helfand and Labraga 1988)
            Real alphac  = (qvel(i,j,k) > qe) ? 1.0 : qvel(i,j,k) / (qe + eps);

            // Level 2.5 stability functions
            Real SM, SH, SQ;
            mynn.calc_stability_funcs(SM,SH,SQ,GM,GH,alphac);

            // Clip SM, SH following WRF
            SM = amrex::min(amrex::max(SM,mynn.SMmin), mynn.SMmax);
            SH = amrex::min(amrex::max(SH,mynn.SHmin), mynn.SHmax);

            // Finally, compute the eddy viscosity/diffusivities
            const Real rho = cell_data(i,j,k,Rho_comp);
            K_turb(i,j,k,EddyDiff::Mom_v)   = rho * Lm * qvel(i,j,k) * SM * 0.5;
            K_turb(i,j,k,EddyDiff::Theta_v) = rho * Lm * qvel(i,j,k) * SH;
            K_turb(i,j,k,EddyDiff::QKE_v)   = rho * Lm * qvel(i,j,k) * SQ;

            // TODO: implement partial-condensation scheme?
            // Currently, implementation matches NN09 without rain (i.e.,
            // the liquid water potential temperature is equal to the
            // potential temperature.

            // NN09 gives the total water content flux; this assumes that
            // all the species have the same eddy diffusivity
            if (mynn.diffuse_moistvars) {
                K_turb(i,j,k,EddyDiff::Q_v) = rho * Lm * qvel(i,j,k) * SH;
            }

            K_turb(i,j,k,EddyDiff::PBL_lengthscale) = Lm;
        });
    }
}
