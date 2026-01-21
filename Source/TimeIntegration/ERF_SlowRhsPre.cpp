
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_ArrayLim.H"
#include "AMReX_BCRec.H"
#include "AMReX_GpuContainers.H"
#include "AMReX_GpuPrint.H"

#include "ERF_TI_slow_headers.H"
#include "ERF_EOS.H"
#include "ERF_Utils.H"
#include "ERF_Diffusion.H"
#include "ERF_EBAdvection.H"
#include "ERF_EB.H"
#include "ERF_SurfaceLayer.H"

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in   ] level level of resolution
 * @param[in   ] finest_level finest level of resolution
 * @param[in   ] nrk   which RK stage
 * @param[in   ] dt    slow time step
 * @param[  out] S_rhs RHS computed here
 * @param[in   ] S_old  old-time solution -- used only for anelastic
 * @param[in   ] S_data current solution
 * @param[in   ] S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[inout] avg_xmom
 * @param[inout] avg_ymom
 * @param[inout] avg_zmom
 * @param[in   ] xvel x-component of velocity
 * @param[in   ] yvel y-component of velocity
 * @param[in   ] zvel z-component of velocity
 * @param[in   ] z_t_ mf rate of change of grid height -- only relevant for moving terrain
 * @param[in   ] cc_src source terms for conserved variables
 * @param[in   ] xmom_src source terms for x-momentum
 * @param[in   ] ymom_src source terms for y-momentum
 * @param[in   ] zmom_src source terms for z-momentum
 * @param[in   ] buoyancy buoyancy source term for z-momentum
 * @param[in   ] zmom_crse_rhs update term from coarser level for z-momentum; non-zero on c/f boundary only
 * @param[in   ] Tau_lev components of stress tensor
 * @param[in   ] SmnSmn strain rate magnitude
 * @param[in   ] eddyDiffs diffusion coefficients for LES turbulence models
 * @param[in   ] Hfx3 heat flux in z-dir
 * @param[in   ] Diss dissipation of turbulent kinetic energy
 * @param[in   ] geom   Container for geometric information
 * @param[in   ] solverChoice  Container for solver parameters
 * @param[in   ] SurfLayer  Pointer to SurfaceLayer class for Monin-Obukhov Similarity Theory boundary condition
 * @param[in   ] domain_bcs_type_d device vector for domain boundary conditions
 * @param[in   ] domain_bcs_type_h   host vector for domain boundary conditions
 * @param[in   ] z_phys_nd height coordinate at nodes
 * @param[in   ] ax area fractions on x-faces
 * @param[in   ] ay area fractions on y-faces
 * @param[in   ] az area fractions on z-faces
 * @param[in   ] detJ Jacobian of the metric transformation (= 1 if use_terrain_fitted_coords is false)
 * @param[in   ] gradp  pressure gradient
 * @param[in   ] mapfac map factors
 * @param[in   ] ebfact EB factories for cell- and face-centered variables
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 */

void erf_slow_rhs_pre (int level, int finest_level,
                       int nrk,
                       Real dt,
                       Vector<MultiFab>& S_rhs,
                       Vector<MultiFab>& S_old,
                       Vector<MultiFab>& S_data,
                       const MultiFab& S_prim,
                       const MultiFab& qt,
                       MultiFab& avg_xmom,
                       MultiFab& avg_ymom,
                       MultiFab& avg_zmom,
                       const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& zvel,
                       std::unique_ptr<MultiFab>& z_t_mf,
                       const MultiFab& cc_src,
                       const MultiFab& xmom_src,
                       const MultiFab& ymom_src,
                       const MultiFab& zmom_src,
                       const MultiFab& buoyancy,
                       const MultiFab* zmom_crse_rhs,
                       Vector<std::unique_ptr<MultiFab>>& Tau_lev,
                       Vector<std::unique_ptr<MultiFab>>& Tau_corr_lev,
                       MultiFab* SmnSmn,
                       MultiFab* eddyDiffs,
                       MultiFab* Hfx1, MultiFab* Hfx2, MultiFab* Hfx3,
                       MultiFab* Q1fx1, MultiFab* Q1fx2,
                       MultiFab* Q1fx3, MultiFab* Q2fx3,
                       MultiFab* Diss,
                       const Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<SurfaceLayer>& SurfLayer,
                       const Gpu::DeviceVector<BCRec>& domain_bcs_type_d,
                       const Vector<BCRec>& domain_bcs_type_h,
                       const MultiFab& z_phys_nd,
                       const MultiFab& z_phys_cc,
                       const MultiFab& ax, const MultiFab& ay, const MultiFab& az,
                       const MultiFab& detJ,
                       Gpu::DeviceVector<Real>& stretched_dz_d,
                       Vector<MultiFab>& gradp,
                       Vector<std::unique_ptr<MultiFab>>& mapfac,
                       const eb_& ebfact,
#ifdef ERF_USE_SHOC
                       std::unique_ptr<SHOCInterface>& shoc_lev,
#endif
                       YAFluxRegister* fr_as_crse,
                       YAFluxRegister* fr_as_fine)
{
    BL_PROFILE_REGION("erf_slow_rhs_pre()");

    const BCRec* bc_ptr_d = domain_bcs_type_d.data();
    const BCRec* bc_ptr_h = domain_bcs_type_h.data();

    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const MultiFab*  t_mean_mf = nullptr;
    if (SurfLayer) { t_mean_mf = SurfLayer->get_mac_avg(level,2); }

    const Box& domain = geom.Domain();
    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);

    const AdvType l_horiz_adv_type = solverChoice.advChoice.dycore_horiz_adv_type;
    const AdvType l_vert_adv_type  = solverChoice.advChoice.dycore_vert_adv_type;
    const Real    l_horiz_upw_frac = solverChoice.advChoice.dycore_horiz_upw_frac;
    const Real    l_vert_upw_frac  = solverChoice.advChoice.dycore_vert_upw_frac;
    const bool    l_use_stretched_dz             = (solverChoice.mesh_type == MeshType::StretchedDz);
    const bool    l_use_terrain_fitted_coords    = (solverChoice.mesh_type == MeshType::VariableDz);
    const bool    l_moving_terrain               = (solverChoice.terrain_type == TerrainType::MovingFittedMesh);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT (l_use_stretched_dz || l_use_terrain_fitted_coords);

    const bool l_use_diff       = ( (dc.molec_diff_type != MolecDiffType::None) ||
                                    (tc.les_type        !=       LESType::None) ||
                                    (tc.rans_type       !=      RANSType::None) ||
                                    (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_turb       = tc.use_kturb;
    const bool l_need_SmnSmn    = tc.use_keqn;

    const Real l_vert_implicit_fac = (solverChoice.vert_implicit_fac[nrk] > 0 &&
                                      solverChoice.implicit_thermal_diffusion);

    const bool l_use_moisture  = (solverChoice.moisture_type != MoistureType::None);
    const bool l_use_SurfLayer = (SurfLayer != nullptr);
    const bool l_rotate        = (solverChoice.use_rotate_surface_flux);

    const bool l_anelastic = (solverChoice.anelastic[level]     == 1);
    const bool l_fixed_rho = (solverChoice.fixed_density[level] == 1);

    const bool l_reflux = ( (solverChoice.coupling_type == CouplingType::TwoWay) && (finest_level > 0) &&
                            ( (l_anelastic && nrk == 1) || (!l_anelastic && nrk == 2) ) );

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const Real* dx = geom.CellSize();

    // *****************************************************************************
    // Combine external forcing terms
    // *****************************************************************************
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // **************************************************************************************
    // If doing advection with EB we need the extra values for tangential interpolation
    // **************************************************************************************
    if (solverChoice.terrain_type == TerrainType::EB) {
        S_data[IntVars::xmom].FillBoundary(geom.periodicity());
        S_data[IntVars::ymom].FillBoundary(geom.periodicity());
        S_data[IntVars::zmom].FillBoundary(geom.periodicity());
    }

    // *****************************************************************************
    // Pre-computed quantities
    // *****************************************************************************
    int nvars                     = S_data[IntVars::cons].nComp();
    const BoxArray& ba            = S_data[IntVars::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVars::cons].DistributionMap();

    int nGhost = (solverChoice.terrain_type == TerrainType::EB) ? 2 : 1;
    MultiFab Omega(convert(ba,IntVect(0,0,1)), dm, 1, nGhost);

    std::unique_ptr<MultiFab> expr;
    std::unique_ptr<MultiFab> dflux_x;
    std::unique_ptr<MultiFab> dflux_y;
    std::unique_ptr<MultiFab> dflux_z;

    if (l_use_diff) {
#ifdef ERF_USE_SHOC
        if (solverChoice.use_shoc) {
            // Populate vertical component of eddyDiffs
            shoc_lev->set_eddy_diffs();
        }
#endif

        erf_make_tau_terms(level,nrk,domain_bcs_type_h,z_phys_nd,
                           S_data,xvel,yvel,zvel,
                           Tau_lev,Tau_corr_lev,
                           SmnSmn,eddyDiffs,geom,solverChoice,SurfLayer,
                           stretched_dz_d, detJ,mapfac, ax, ay, az, ebfact);

        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, nvars, 0);

#ifdef ERF_USE_SHOC
        if (solverChoice.use_shoc) {
            // Zero out the surface stresses of tau13/tau23
            shoc_lev->set_diff_stresses();
        } else if (l_use_SurfLayer) {
            // Set surface shear stresses, update heat and moisture fluxes
            // (fluxes will be later applied in the diffusion source update)
            Vector<const MultiFab*> mfs = {&S_data[IntVars::cons], &xvel, &yvel, &zvel};
            SurfLayer->impose_SurfaceLayer_bcs(level, mfs, Tau_lev,
                                               Hfx1, Hfx2, Hfx3,
                                               Q1fx1, Q1fx2, Q1fx3,
                                               &z_phys_nd);
        }
#else
        // This is computed pre step in Advance if we use SHOC
        if (l_use_SurfLayer) {
            // Set surface shear stresses, update heat and moisture fluxes
            // (fluxes will be later applied in the diffusion source update)
            Vector<const MultiFab*> mfs = {&S_data[IntVars::cons], &xvel, &yvel, &zvel};
            SurfLayer->impose_SurfaceLayer_bcs(level, mfs, Tau_lev,
                                               Hfx1, Hfx2, Hfx3,
                                               Q1fx1, Q1fx2, Q1fx3,
                                               &z_phys_nd);

            //if (l_vert_implicit_fac > 0 && solverChoice.implicit_momentum_diffusion) {
            //    copy_surface_tau_for_implicit(Tau_lev, Tau_corr_lev);
            //}
        }
#endif
    } // l_use_diff

    // This is just cautionary to deal with grid boundaries that aren't domain boundaries
    S_rhs[IntVars::zmom].setVal(0.0);

    // *****************************************************************************
    // Define updates and fluxes in the current RK stage
    // *****************************************************************************
    // Cell-centered masks for EB (used for flux interpolation)
    bool already_on_centroids = false;
    Vector<iMultiFab> physbnd_mask;
    physbnd_mask.resize(IntVars::NumTypes);
    if (solverChoice.terrain_type == TerrainType::EB) {
        physbnd_mask[IntVars::cons].define(S_data[IntVars::cons].boxArray(), S_data[IntVars::cons].DistributionMap(), 1, 1);
        physbnd_mask[IntVars::cons].BuildMask(geom.Domain(), geom.periodicity(), 1, 1, 0, 1);
        // physbnd_mask[IntVars::cons].FillBoundary(geom.periodicity());
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            physbnd_mask[1+dir].define(S_data[1+dir].boxArray(), S_data[1+dir].DistributionMap(), 1, 1);
            physbnd_mask[1+dir].BuildMask(geom.Domain(), geom.periodicity(), 1, 1, 0, 1);
            // physbnd_mask[1+dir].FillBoundary(geom.periodicity());
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
    for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        // Boxes for momentum fluxes
        Vector<Box> tbx_grown(AMREX_SPACEDIM);
        Vector<Box> tby_grown(AMREX_SPACEDIM);
        Vector<Box> tbz_grown(AMREX_SPACEDIM);
        if (solverChoice.terrain_type == TerrainType::EB) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                tbx_grown[dir] = tbx;
                tby_grown[dir] = tby;
                tbz_grown[dir] = tbz;
                IntVect iv(1, 1, 1);
                iv[dir] = 0;
                tbx_grown[dir] = (tbx_grown[dir].growHi(dir,1)).grow(iv);
                tby_grown[dir] = (tby_grown[dir].growHi(dir,1)).grow(iv);
                tbz_grown[dir] = (tbz_grown[dir].growHi(dir,1)).grow(iv);
            }
        }

        // We don't compute a source term for z-momentum on the bottom or top domain boundary
        if (tbz.smallEnd(2) == domain.smallEnd(2)) {
            tbz.growLo(2,-1);
        }
        if (tbz.bigEnd(2) == domain.bigEnd(2)+1) {
            tbz.growHi(2,-1);
        }

        const Array4<const Real> & cell_data  = S_data[IntVars::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real>       & cell_rhs   = S_rhs[IntVars::cons].array(mfi);

        const Array4<const Real> & cell_old   = S_old[IntVars::cons].array(mfi);

        const Array4<Real const>& xmom_src_arr   = xmom_src.const_array(mfi);
        const Array4<Real const>& ymom_src_arr   = ymom_src.const_array(mfi);
        const Array4<Real const>& zmom_src_arr   = zmom_src.const_array(mfi);
        const Array4<Real const>& buoyancy_arr   = buoyancy.const_array(mfi);

        const Array4<Real const>& gpx_arr   = gradp[GpVars::gpx].const_array(mfi);
        const Array4<Real const>& gpy_arr   = gradp[GpVars::gpy].const_array(mfi);
        const Array4<Real const>& gpz_arr   = gradp[GpVars::gpz].const_array(mfi);

        const Array4<Real const>&  qt_arr   = qt.const_array(mfi);

        const Array4<Real>& rho_u_old = S_old[IntVars::xmom].array(mfi);
        const Array4<Real>& rho_v_old = S_old[IntVars::ymom].array(mfi);

        if (l_anelastic) {
            // When anelastic we must reset these to 0 each RK step
            avg_xmom[mfi].template setVal<RunOn::Device>(0.0,tbx);
            avg_ymom[mfi].template setVal<RunOn::Device>(0.0,tby);
            avg_zmom[mfi].template setVal<RunOn::Device>(0.0,tbz);
        }

        Array4<Real> avg_xmom_arr = avg_xmom.array(mfi);
        Array4<Real> avg_ymom_arr = avg_ymom.array(mfi);
        Array4<Real> avg_zmom_arr = avg_zmom.array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<const Real>& rho_u = S_data[IntVars::xmom].array(mfi);
        const Array4<const Real>& rho_v = S_data[IntVars::ymom].array(mfi);
        const Array4<const Real>& rho_w = S_data[IntVars::zmom].array(mfi);

        // Map factors
        const Array4<const Real>& mf_mx  = mapfac[MapFacType::m_x]->const_array(mfi);
        const Array4<const Real>& mf_ux  = mapfac[MapFacType::u_x]->const_array(mfi);
        const Array4<const Real>& mf_vx  = mapfac[MapFacType::v_x]->const_array(mfi);
        const Array4<const Real>& mf_my  = mapfac[MapFacType::m_y]->const_array(mfi);
        const Array4<const Real>& mf_uy  = mapfac[MapFacType::u_y]->const_array(mfi);
        const Array4<const Real>& mf_vy  = mapfac[MapFacType::v_y]->const_array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        Array4<const Real> z_t;
        if (z_t_mf) {
            z_t = z_t_mf->array(mfi);
        } else {
            z_t = Array4<const Real>{};
        }

        const Array4<Real>& rho_u_rhs = S_rhs[IntVars::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVars::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVars::zmom].array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        // Terrain metrics
        const Array4<const Real>& z_nd = z_phys_nd.const_array(mfi);
        const Array4<const Real>& z_cc = z_phys_cc.const_array(mfi);

        // *****************************************************************************
        // Define flux arrays for use in advection
        // *****************************************************************************
        std::array<FArrayBox,AMREX_SPACEDIM> flux;
        std::array<FArrayBox,AMREX_SPACEDIM> flux_u;
        std::array<FArrayBox,AMREX_SPACEDIM> flux_v;
        std::array<FArrayBox,AMREX_SPACEDIM> flux_w;

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (solverChoice.terrain_type != TerrainType::EB) {
                flux[dir].resize(surroundingNodes(bx,dir),2,The_Async_Arena());
            } else {
                flux[dir].resize(surroundingNodes(bx,dir).grow(1),2,The_Async_Arena());
            }
            flux[dir].setVal<RunOn::Device>(0.);
        }
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>
            flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};

        // Define flux arrays for momentum variables (used only for EB now)
        GpuArray<Array4<Real>, AMREX_SPACEDIM> flx_u_arr{};
        GpuArray<Array4<Real>, AMREX_SPACEDIM> flx_v_arr{};
        GpuArray<Array4<Real>, AMREX_SPACEDIM> flx_w_arr{};

        if (solverChoice.terrain_type == TerrainType::EB) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                flux_u[dir].resize(tbx_grown[dir],1,The_Async_Arena());
                flux_v[dir].resize(tby_grown[dir],1,The_Async_Arena());
                flux_w[dir].resize(tbz_grown[dir],1,The_Async_Arena());
                flux_u[dir].setVal<RunOn::Device>(0.);
                flux_v[dir].setVal<RunOn::Device>(0.);
                flux_w[dir].setVal<RunOn::Device>(0.);
                flx_u_arr[dir] = flux_u[dir].array();
                flx_v_arr[dir] = flux_v[dir].array();
                flx_w_arr[dir] = flux_w[dir].array();
            }
        }

        // *****************************************************************************
        // Contravariant flux field
        // *****************************************************************************
        {
        BL_PROFILE("slow_rhs_making_omega");
            IntVect nGrowVect = (solverChoice.terrain_type == TerrainType::EB)
                                ? IntVect(AMREX_D_DECL(2, 2, 2)) : IntVect(AMREX_D_DECL(1, 1, 1));
            Box gbxo = surroundingNodes(bx,2); gbxo.grow(nGrowVect);
            //
            // Now create Omega with momentum (not velocity) with z_t subtracted if moving terrain
            // ONLY if not doing anelastic + terrain -- in that case Omega will be defined coming
            // out of the projection
            //
            if (!l_use_terrain_fitted_coords) {
                ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = rho_w(i,j,k);
                });

            } else {

                Box gbxo_lo = gbxo; gbxo_lo.setBig(2,domain.smallEnd(2));
                int lo_z_face = domain.smallEnd(2);
                if (gbxo_lo.smallEnd(2) <= lo_z_face) {
                    ParallelFor(gbxo_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = 0.;
                    });
                }
                Box gbxo_hi = gbxo; gbxo_hi.setSmall(2,gbxo.bigEnd(2));
                int hi_z_face = domain.bigEnd(2)+1;
                if (gbxo_hi.bigEnd(2) >= hi_z_face) {
                    ParallelFor(gbxo_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = rho_w(i,j,k);
                    });
                }

                if (z_t) { // Note we never do anelastic with moving terrain
                    Box gbxo_mid = gbxo; gbxo_mid.setSmall(2,1); gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
                    ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        // We define rho on the z-face the same way as in MomentumToVelocity/VelocityToMomentum
                        Real rho_at_face = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp));
                        omega_arr(i,j,k) = OmegaFromW(i,j,k,rho_w(i,j,k),
                                                      rho_u,rho_v,mf_ux,mf_vy,z_nd,dxInv) -
                            rho_at_face * z_t(i,j,k);
                    });
                } else {
                    Box gbxo_mid = gbxo;
                    if (gbxo_mid.smallEnd(2) <= domain.smallEnd(2)) {
                        gbxo_mid.setSmall(2,1);
                    }
                    if (gbxo_mid.bigEnd(2) >= domain.bigEnd(2)+1) {
                        gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
                    }
                    ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = OmegaFromW(i,j,k,rho_w(i,j,k),
                                                      rho_u,rho_v,mf_ux,mf_vy,z_nd,dxInv);
                    });
                }
            }
        } // end profile

        // *****************************************************************************
        // Diffusive terms (pre-computed above)
        // *****************************************************************************
        // No terrain diffusion
        Array4<Real> tau11,tau22,tau33;
        Array4<Real> tau12,tau13,tau23;
        if (Tau_lev[TauType::tau11]) {
            tau11 = Tau_lev[TauType::tau11]->array(mfi); tau22 = Tau_lev[TauType::tau22]->array(mfi);
            tau33 = Tau_lev[TauType::tau33]->array(mfi); tau12 = Tau_lev[TauType::tau12]->array(mfi);
            tau13 = Tau_lev[TauType::tau13]->array(mfi); tau23 = Tau_lev[TauType::tau23]->array(mfi);
        } else {
            tau11 = Array4<Real>{}; tau22 = Array4<Real>{}; tau33 = Array4<Real>{};
            tau12 = Array4<Real>{}; tau13 = Array4<Real>{}; tau23 = Array4<Real>{};
        }
        // Terrain diffusion
        Array4<Real> tau21,tau31,tau32;
        if (Tau_lev[TauType::tau21]) {
            tau21 = Tau_lev[TauType::tau21]->array(mfi);
            tau31 = Tau_lev[TauType::tau31]->array(mfi);
            tau32 = Tau_lev[TauType::tau32]->array(mfi);
        } else {
            tau21 = Array4<Real>{}; tau31 = Array4<Real>{}; tau32 = Array4<Real>{};
        }

        // Strain magnitude
        Array4<Real> SmnSmn_a;
        if (l_need_SmnSmn) {
            SmnSmn_a = SmnSmn->array(mfi);
        } else {
            SmnSmn_a = Array4<Real>{};
        }

        // *****************************************************************************
        // Define updates in the RHS of continuity and potential temperature equations
        // *****************************************************************************
        bool l_eb_terrain_cc = false; // EB terrain on cell-centered grid
        Array4<const int> mask_arr{};
        Array4<const EBCellFlag> cfg_arr{};
        Array4<const Real> ax_arr{};
        Array4<const Real> ay_arr{};
        Array4<const Real> az_arr{};
        Array4<const Real> fcx_arr{};
        Array4<const Real> fcy_arr{};
        Array4<const Real> fcz_arr{};
        Array4<const Real> detJ_arr{};

        if (solverChoice.terrain_type == TerrainType::EB)
        {
            EBCellFlagFab const& cfg = (ebfact.get_const_factory())->getMultiEBCellFlagFab()[mfi];
            cfg_arr  = cfg.const_array();
            if (cfg.getType(bx) == FabType::singlevalued) {
                l_eb_terrain_cc = true;
                ax_arr   = (ebfact.get_const_factory())->getAreaFrac()[0]->const_array(mfi);
                ay_arr   = (ebfact.get_const_factory())->getAreaFrac()[1]->const_array(mfi);
                az_arr   = (ebfact.get_const_factory())->getAreaFrac()[2]->const_array(mfi);
                fcx_arr  = (ebfact.get_const_factory())->getFaceCent()[0]->const_array(mfi);
                fcy_arr  = (ebfact.get_const_factory())->getFaceCent()[1]->const_array(mfi);
                fcz_arr  = (ebfact.get_const_factory())->getFaceCent()[2]->const_array(mfi);
                detJ_arr = (ebfact.get_const_factory())->getVolFrac().const_array(mfi);
                // if (!already_on_centroids) {mask_arr = physbnd_mask[IntVars::cons].const_array(mfi);}
                mask_arr = physbnd_mask[IntVars::cons].const_array(mfi);
            } else {
                ax_arr   = ax.const_array(mfi);
                ay_arr   = ay.const_array(mfi);
                az_arr   = az.const_array(mfi);
                detJ_arr = detJ.const_array(mfi);
            }
        } else {
            ax_arr   = ax.const_array(mfi);
            ay_arr   = ay.const_array(mfi);
            az_arr   = az.const_array(mfi);
            detJ_arr = detJ.const_array(mfi);
        }

        int icomp = RhoTheta_comp; int ncomp = 1;
        if (!l_eb_terrain_cc){
            AdvectionSrcForRho( bx, cell_rhs,
                                rho_u, rho_v, omega_arr,      // these are being used to build the fluxes
                                avg_xmom_arr, avg_ymom_arr, avg_zmom_arr, // these are being defined from the fluxes
                                ax_arr, ay_arr, az_arr, detJ_arr,
                                dxInv, mf_mx, mf_my, mf_uy, mf_vx,
                                flx_arr, l_fixed_rho);
            AdvectionSrcForScalars(bx, icomp, ncomp,
                                avg_xmom_arr, avg_ymom_arr, avg_zmom_arr,
                                cell_prim, cell_rhs,
                                detJ_arr, dxInv, mf_mx, mf_my,
                                l_horiz_adv_type, l_vert_adv_type,
                                l_horiz_upw_frac, l_vert_upw_frac,
                                flx_arr, domain, bc_ptr_h);
        } else {
            EBAdvectionSrcForRho(bx, cell_rhs,
                                rho_u, rho_v, omega_arr,
                                avg_xmom_arr, avg_ymom_arr, avg_zmom_arr,
                                mask_arr, cfg_arr,
                                ax_arr, ay_arr, az_arr,
                                fcx_arr, fcy_arr, fcz_arr, detJ_arr,
                                dxInv, mf_mx, mf_my, mf_uy, mf_vx,
                                flx_arr, l_fixed_rho,
                                already_on_centroids);
            EBAdvectionSrcForScalars(bx, icomp, ncomp,
                                avg_xmom_arr, avg_ymom_arr, avg_zmom_arr,
                                cell_prim, cell_rhs,
                                mask_arr, cfg_arr, ax_arr, ay_arr, az_arr,
                                fcx_arr, fcy_arr, fcz_arr,
                                detJ_arr, dxInv, mf_mx, mf_my,
                                l_horiz_adv_type, l_vert_adv_type,
                                l_horiz_upw_frac, l_vert_upw_frac,
                                flx_arr, domain, bc_ptr_h,
                                already_on_centroids);
        }

        if (l_use_diff) {
            Array4<Real> diffflux_x = dflux_x->array(mfi);
            Array4<Real> diffflux_y = dflux_y->array(mfi);
            Array4<Real> diffflux_z = dflux_z->array(mfi);

            Array4<Real> hfx_x = Hfx1->array(mfi);
            Array4<Real> hfx_y = Hfx2->array(mfi);
            Array4<Real> hfx_z = Hfx3->array(mfi);

            Array4<Real> q1fx_x = (Q1fx1) ? Q1fx1->array(mfi) : Array4<Real>{};
            Array4<Real> q1fx_y = (Q1fx2) ? Q1fx2->array(mfi) : Array4<Real>{};
            Array4<Real> q1fx_z = (Q1fx3) ? Q1fx3->array(mfi) : Array4<Real>{};

            Array4<Real> q2fx_z = (Q2fx3) ? Q2fx3->array(mfi) : Array4<Real>{};
            Array4<Real> diss  = Diss->array(mfi);

            const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

            // NOTE: No diffusion for continuity, so n starts at 1.
            int n_start = RhoTheta_comp;
            int n_comp  = 1;

            // For l_vert_implicit_fac > 0, we scale the rho*theta contribution
            // by (1 - implicit_fac) and add in the implicit contribution with
            // ERF_Implicit.H
            if (l_use_stretched_dz) {
                DiffusionSrcForState_S(bx, domain, n_start, n_comp, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       stretched_dz_d, dxInv, SmnSmn_a,
                                       mf_mx, mf_ux, mf_vx,
                                       mf_my, mf_uy, mf_vy,
                                       hfx_z, q1fx_z, q2fx_z, diss,
                                       mu_turb, solverChoice, level,
                                       tm_arr, grav_gpu, bc_ptr_d, l_use_SurfLayer, l_vert_implicit_fac);
            } else if (l_use_terrain_fitted_coords) {
                DiffusionSrcForState_T(bx, domain, n_start, n_comp, l_rotate, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       z_nd, z_cc, ax_arr, ay_arr, az_arr, detJ_arr,
                                       dxInv, SmnSmn_a,
                                       mf_mx, mf_ux, mf_vx,
                                       mf_my, mf_uy, mf_vy,
                                       hfx_x, hfx_y, hfx_z, q1fx_x, q1fx_y, q1fx_z, q2fx_z, diss,
                                       mu_turb, solverChoice, level,
                                       tm_arr, grav_gpu, bc_ptr_d, l_use_SurfLayer, l_vert_implicit_fac);
            } else {
                DiffusionSrcForState_N(bx, domain, n_start, n_comp, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       dxInv, SmnSmn_a,
                                       mf_mx, mf_ux, mf_vx,
                                       mf_my, mf_uy, mf_vy,
                                       hfx_z, q1fx_z, q2fx_z, diss,
                                       mu_turb, solverChoice, level,
                                       tm_arr, grav_gpu, bc_ptr_d, l_use_SurfLayer, l_vert_implicit_fac);
            }
        }

        const Array4<Real const>& source_arr   = cc_src.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cell_rhs(i,j,k,Rho_comp)      += source_arr(i,j,k,Rho_comp);
            cell_rhs(i,j,k,RhoTheta_comp) += source_arr(i,j,k,RhoTheta_comp);
        });

        // If anelastic and in second RK stage, take average of old-time and new-time source
        if ( l_anelastic && (nrk == 1) )
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_rhs(i,j,k,     Rho_comp) *= 0.5;
                cell_rhs(i,j,k,RhoTheta_comp) *= 0.5;

                cell_rhs(i,j,k,     Rho_comp) += 0.5 / dt * (cell_data(i,j,k,     Rho_comp) - cell_old(i,j,k,     Rho_comp));
                cell_rhs(i,j,k,RhoTheta_comp) += 0.5 / dt * (cell_data(i,j,k,RhoTheta_comp) - cell_old(i,j,k,RhoTheta_comp));
            });
        }

        // *****************************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *****************************************************************************
        int lo_z_face = domain.smallEnd(2);
        int hi_z_face = domain.bigEnd(2)+1;

        AdvectionSrcForMom(mfi, bx, tbx, tby, tbz, tbx_grown, tby_grown, tbz_grown,
                           rho_u_rhs, rho_v_rhs, rho_w_rhs,
                           cell_data, u, v, w,
                           rho_u, rho_v, omega_arr,
                           z_nd, ax_arr, ay_arr, az_arr,
                           detJ_arr, stretched_dz_d,
                           dxInv, mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy,
                           l_horiz_adv_type, l_vert_adv_type,
                           l_horiz_upw_frac, l_vert_upw_frac,
                           solverChoice.mesh_type, solverChoice.terrain_type,
                           ebfact, flx_u_arr, flx_v_arr, flx_w_arr,
                           physbnd_mask, already_on_centroids,
                           lo_z_face, hi_z_face, domain, bc_ptr_h);

        if (l_use_diff) {
            // Note: tau** were calculated with calls to
            // ComputeStress[Cons|Var]Visc_[N|S|T] in which ConsVisc ("constant
            // viscosity") means that there is no contribution from a
            // turbulence model. However, whether this field truly is constant
            // depends on whether MolecDiffType is Constant or ConstantAlpha.
            if (solverChoice.terrain_type != TerrainType::EB) {
                DiffusionSrcForMom(tbx, tby, tbz,
                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                    tau11, tau22, tau33,
                    tau12, tau21, tau13, tau31, tau23, tau32,
                    detJ_arr, stretched_dz_d, dxInv,
                    mf_mx, mf_ux, mf_vx,
                    mf_my, mf_uy, mf_vy,
                    l_use_stretched_dz,
                    l_use_terrain_fitted_coords);
            } else {
                DiffusionSrcForMom_EB(mfi, domain, tbx, tby, tbz,
                    rho_u_rhs, rho_v_rhs, rho_w_rhs,
                    u, v, w,
                    tau11, tau22, tau33,
                    tau12, tau13, tau23,
                    dx, dxInv,
                    mf_mx, mf_ux, mf_vx,
                    mf_my, mf_uy, mf_vy,
                    solverChoice, ebfact, bc_ptr_d);
            }
        }

        auto abl_pressure_grad    = solverChoice.abl_pressure_grad;

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // x-momentum equation

            // Note that gradp arrays now carry the map factor in them

            Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i-1,j,k)) : 0.0;

            rho_u_rhs(i, j, k) += (-gpx_arr(i,j,k) - abl_pressure_grad[0]) / (1.0 + q) + xmom_src_arr(i,j,k);

            if (l_moving_terrain) {
                Real h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);
                rho_u_rhs(i, j, k) *= h_zeta;
            }

            if ( l_anelastic && (nrk == 1) ) {
              rho_u_rhs(i,j,k) *= 0.5;
              rho_u_rhs(i,j,k) += 0.5 / dt * (rho_u(i,j,k) - rho_u_old(i,j,k));
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // y-momentum equation

            // Note that gradp arrays now carry the map factor in them

            Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i,j-1,k)) : 0.0;

            rho_v_rhs(i, j, k) += (-gpy_arr(i,j,k) - abl_pressure_grad[1]) / (1.0 + q) + ymom_src_arr(i,j,k);

            if (l_moving_terrain) {
                Real h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
                rho_v_rhs(i, j, k) *= h_zeta;
            }

            if ( l_anelastic && (nrk == 1) ) {
              rho_v_rhs(i,j,k) *= 0.5;
              rho_v_rhs(i,j,k) += 0.5 / dt * (rho_v(i,j,k) - rho_v_old(i,j,k));
            }
        });

        // *****************************************************************************
        // Zero out source terms for x- and y- momenta if at walls or inflow
        // We need to do this -- even though we call the boundary conditions later --
        // because the slow source is used to update the state in the fast interpolater.
        // *****************************************************************************
        if (bx.smallEnd(0) == domain.smallEnd(0)) {
            Box lo_x_dom_face(bx); lo_x_dom_face.setBig(0,bx.smallEnd(0));
            if (bc_ptr_h[BCVars::xvel_bc].lo(0) == ERFBCType::ext_dir) {
                ParallelFor(lo_x_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    rho_u_rhs(i,j,k) = 0.;
                });
            } else if (bc_ptr_h[BCVars::xvel_bc].lo(0) == ERFBCType::ext_dir_upwind) {
                ParallelFor(lo_x_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (u(i,j,k) >= 0.) {
                        rho_u_rhs(i,j,k) = 0.;
                    }
                });
            }
        }
        if (bx.bigEnd(0) == domain.bigEnd(0)) {
            Box hi_x_dom_face(bx); hi_x_dom_face.setSmall(0,bx.bigEnd(0)+1); hi_x_dom_face.setBig(0,bx.bigEnd(0)+1);
            if (bc_ptr_h[BCVars::xvel_bc].hi(0) == ERFBCType::ext_dir) {
                ParallelFor(hi_x_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    rho_u_rhs(i,j,k) = 0.;
                });
            } else if (bc_ptr_h[BCVars::xvel_bc].hi(0) == ERFBCType::ext_dir_upwind) {
                ParallelFor(hi_x_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (u(i,j,k) <= 0.) {
                        rho_u_rhs(i,j,k) = 0.;
                    }
                });
            }
        }
        if (bx.smallEnd(1) == domain.smallEnd(1)) {
            Box lo_y_dom_face(bx); lo_y_dom_face.setBig(1,bx.smallEnd(1));
            if (bc_ptr_h[BCVars::yvel_bc].lo(1) == ERFBCType::ext_dir) {
                ParallelFor(lo_y_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    rho_v_rhs(i,j,k) = 0.;
                });
            } else if (bc_ptr_h[BCVars::yvel_bc].lo(1) == ERFBCType::ext_dir_upwind) {
                ParallelFor(lo_y_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (v(i,j,k) >= 0.) {
                        rho_v_rhs(i,j,k) = 0.;
                    }
                });
            }
        }
        if (bx.bigEnd(1) == domain.bigEnd(1)) {
            Box hi_y_dom_face(bx); hi_y_dom_face.setSmall(1,bx.bigEnd(1)+1); hi_y_dom_face.setBig(1,bx.bigEnd(1)+1);
            if (bc_ptr_h[BCVars::yvel_bc].hi(1) == ERFBCType::ext_dir) {
                ParallelFor(hi_y_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    rho_v_rhs(i,j,k) = 0.;
                });
            } else if (bc_ptr_h[BCVars::yvel_bc].hi(1) == ERFBCType::ext_dir_upwind) {
                ParallelFor(hi_y_dom_face, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (v(i,j,k) <= 0.) {
                        rho_v_rhs(i,j,k) = 0.;
                    }
                });
            }
        }

        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        { // z-momentum equation

            Real gpz = gpz_arr(i,j,k);

            Real q = (l_use_moisture) ? 0.5 * (qt_arr(i,j,k) + qt_arr(i,j,k-1)) : 0.0;

            rho_w_rhs(i, j, k) += (-gpz - abl_pressure_grad[2] + buoyancy_arr(i,j,k)) / (1.0 + q) + zmom_src_arr(i,j,k);

            if (l_moving_terrain) {
                 rho_w_rhs(i, j, k) *= 0.5 * (detJ_arr(i,j,k) + detJ_arr(i,j,k-1));
            }
        });

        auto const lo = lbound(bx);
        auto const hi = ubound(bx);

        // Note: the logic below assumes no tiling in z!
        if (level > 0) {

            const Array4<const Real>& rho_w_rhs_crse = zmom_crse_rhs->const_array(mfi);

            Box b2d = bx; b2d.setRange(2,0);

            if (lo.z > klo) {
                ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int ) // bottom of box but not of domain
                {
                    rho_w_rhs(i,j,lo.z) = rho_w_rhs_crse(i,j,lo.z);
                });
            }

            if (hi.z < khi+1) {
                ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int ) // top of box but not of domain
                {
                    rho_w_rhs(i,j,hi.z+1) = rho_w_rhs_crse(i,j,hi.z+1);
                });
            }
        }

        {
        BL_PROFILE("slow_rhs_pre_fluxreg");
        // We only add to the flux registers in the final RK step
        // NOTE: for now we are only refluxing density not (rho theta) since the latter seems to introduce
        //       a problem at top and bottom boundaries
        if (l_reflux) {
            int strt_comp_reflux = (l_fixed_rho) ? 1 : 0;
            int  num_comp_reflux = 1;
            if (level < finest_level) {
                fr_as_crse->CrseAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dt, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, RunOn::Device);
            }
            if (level > 0) {
                fr_as_fine->FineAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dt, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, RunOn::Device);
            }

        } // two-way coupling
        } // end profile
    } // mfi
    } // OMP
}
