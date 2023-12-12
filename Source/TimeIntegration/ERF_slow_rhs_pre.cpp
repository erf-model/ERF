#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_GpuContainers.H>
#include <ERF_Constants.H>
#include <Advection.H>
#include <Diffusion.H>
#include <NumericalDiffusion.H>
#include <TI_headers.H>
#include <TileNoZ.H>
#include <EOS.H>
#include <ERF.H>

#include <TerrainMetrics.H>
#include <IndexDefines.H>
#include <PlaneAverage.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in]  level level of resolution
 * @param[in]  finest_level finest level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[out]  S_rhs RHS computed here
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in]  S_scratch scratch space
 * @param[in]  xvel x-component of velocity
 * @param[in]  yvel y-component of velocity
 * @param[in]  zvel z-component of velocity
 * @param[in]  qv   water vapor
 * @param[in]  z_t_ mf rate of change of grid height -- only relevant for moving terrain
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] source source terms for conserved variables
 * @param[in] buoyancy buoyancy source term
 * @param[in] Tau11 tau_11 component of stress tensor
 * @param[in] Tau22 tau_22 component of stress tensor
 * @param[in] Tau33 tau_33 component of stress tensor
 * @param[in] Tau12 tau_12 component of stress tensor
 * @param[in] Tau12 tau_13 component of stress tensor
 * @param[in] Tau21 tau_21 component of stress tensor
 * @param[in] Tau23 tau_23 component of stress tensor
 * @param[in] Tau31 tau_31 component of stress tensor
 * @param[in] Tau32 tau_32 component of stress tensor
 * @param[in] SmnSmn strain rate magnitude
 * @param[in] eddyDiffs diffusion coefficients for LES turbulence models
 * @param[in] Hfx3 heat flux in z-dir
 * @param[in] Diss dissipation of turbulent kinetic energy
 * @param[in]  geom   Container for geometric informaiton
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in]  most  Pointer to MOST class for Monin-Obukhov Similarity Theory boundary condition
 * @param[in]  domain_bcs_type_d device vector for domain boundary conditions
 * @param[in]  domain_bcs_type     host vector for domain boundary conditions
 * @param[in]  domain_bcs_type     host vector for domain boundary conditions
 * @param[in] z_phys_nd height coordinate at nodes
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain is false)
 * @param[in]  p0     Reference (hydrostatically stratified) pressure
 * @param[in] mapfac_m map factor at cell centers
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 * @param[inout] fr_as_crse YAFluxRegister at level l at level l   / l+1 interface
 * @param[inout] fr_as_fine YAFluxRegister at level l at level l-1 / l   interface
 * @param[in] dptr_rayleigh_tau  strength of Rayleigh damping
 * @param[in] dptr_rayleigh_ubar reference value for x-velocity used to define Rayleigh damping
 * @param[in] dptr_rayleigh_vbar reference value for y-velocity used to define Rayleigh damping
 * @param[in] dptr_rayleigh_wbar reference value for z-velocity used to define Rayleigh damping
 * @param[in] dptr_rayleigh_thetabar reference value for potential temperature used to define Rayleigh damping
 */

void erf_slow_rhs_pre (int level, int finest_level,
                       int nrk,
                       amrex::Real dt,
                       Vector<MultiFab>& S_rhs,
                       Vector<MultiFab>& S_data,
                       const MultiFab& S_prim,
                       Vector<MultiFab>& S_scratch,
                       const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& zvel,
                       const MultiFab* qv,
                       std::unique_ptr<MultiFab>& z_t_mf,
                       MultiFab& Omega,
                       const MultiFab& source,
                       const MultiFab& buoyancy,
                       MultiFab* Tau11, MultiFab* Tau22, MultiFab* Tau33,
                       MultiFab* Tau12, MultiFab* Tau13, MultiFab* Tau21,
                       MultiFab* Tau23, MultiFab* Tau31, MultiFab* Tau32,
                       MultiFab* SmnSmn,
                       MultiFab* eddyDiffs,
                       MultiFab* Hfx3, MultiFab* Diss,
                       const amrex::Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<ABLMost>& most,
                       const Gpu::DeviceVector<amrex::BCRec>& domain_bcs_type_d,
                       const Vector<amrex::BCRec>& domain_bcs_type,
                       std::unique_ptr<MultiFab>& z_phys_nd, std::unique_ptr<MultiFab>& detJ,
                       const MultiFab* p0,
                       std::unique_ptr<MultiFab>& mapfac_m,
                       std::unique_ptr<MultiFab>& mapfac_u,
                       std::unique_ptr<MultiFab>& mapfac_v,
                       YAFluxRegister* fr_as_crse,
                       YAFluxRegister* fr_as_fine,
                       const amrex::Real* dptr_rayleigh_tau, const amrex::Real* dptr_rayleigh_ubar,
                       const amrex::Real* dptr_rayleigh_vbar, const amrex::Real* dptr_rayleigh_wbar,
                       const amrex::Real* dptr_rayleigh_thetabar)
{
    BL_PROFILE_REGION("erf_slow_rhs_pre()");

    DiffChoice dc = solverChoice.diffChoice;
    TurbChoice tc = solverChoice.turbChoice[level];

    const MultiFab* t_mean_mf = nullptr;
    if (most) t_mean_mf = most->get_mac_avg(0,2);

    int start_comp = 0;
    int   num_comp = 2;
    int   end_comp = start_comp + num_comp - 1;

    const AdvType l_horiz_adv_type = solverChoice.advChoice.dycore_horiz_adv_type;
    const AdvType l_vert_adv_type  = solverChoice.advChoice.dycore_vert_adv_type;
    const bool    l_use_terrain    = solverChoice.use_terrain;
    const bool    l_moving_terrain = (solverChoice.terrain_type == TerrainType::Moving);
    if (l_moving_terrain) AMREX_ALWAYS_ASSERT (l_use_terrain);

    const bool l_reflux = (solverChoice.coupling_type == CouplingType::TwoWay);

    const bool l_use_ndiff      = solverChoice.use_NumDiff;
    const bool l_use_diff       = ( (dc.molec_diff_type != MolecDiffType::None) ||
                                    (tc.les_type        !=       LESType::None) ||
                                    (tc.pbl_type        !=       PBLType::None) );
    const bool l_use_turb       = ( tc.les_type == LESType::Smagorinsky ||
                                    tc.les_type == LESType::Deardorff   ||
                                    tc.pbl_type == PBLType::MYNN25      ||
                                    tc.pbl_type == PBLType::YSU );

    const bool use_moisture = (solverChoice.moisture_type != MoistureType::None);

    const amrex::BCRec* bc_ptr   = domain_bcs_type_d.data();
    const amrex::BCRec* bc_ptr_h = domain_bcs_type.data();

    const Box& domain = geom.Domain();
    const int domhi_z = domain.bigEnd()[2];

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const Real* dx = geom.CellSize();

    // *************************************************************************
    // Combine external forcing terms
    // *************************************************************************
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *************************************************************************
    // Pre-computed quantities
    // *************************************************************************
    int nvars                     = S_data[IntVar::cons].nComp();
    const BoxArray& ba            = S_data[IntVar::cons].boxArray();
    const DistributionMapping& dm = S_data[IntVar::cons].DistributionMap();

    std::unique_ptr<MultiFab> expr;
    std::unique_ptr<MultiFab> dflux_x;
    std::unique_ptr<MultiFab> dflux_y;
    std::unique_ptr<MultiFab> dflux_z;

    if (l_use_diff) {
        expr    = std::make_unique<MultiFab>(ba  , dm, 1, IntVect(1,1,0));
        dflux_x = std::make_unique<MultiFab>(convert(ba,IntVect(1,0,0)), dm, nvars, 0);
        dflux_y = std::make_unique<MultiFab>(convert(ba,IntVect(0,1,0)), dm, nvars, 0);
        dflux_z = std::make_unique<MultiFab>(convert(ba,IntVect(0,0,1)), dm, nvars, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(S_data[IntVar::cons],TileNoZ()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& valid_bx = mfi.validbox();

            // Velocities
            const Array4<const Real> & u = xvel.array(mfi);
            const Array4<const Real> & v = yvel.array(mfi);
            const Array4<const Real> & w = zvel.array(mfi);

            // Contravariant velocity
            const Array4<Real>& omega_arr = Omega.array(mfi);

            // Map factors
            const Array4<const Real>& mf_m   = mapfac_m->const_array(mfi);
            const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
            const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

            // Eddy viscosity
            const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

            // Terrain metrics
            const Array4<const Real>& z_nd     = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
            const Array4<const Real>& detJ_arr = l_use_terrain ?      detJ->const_array(mfi) : Array4<const Real>{};

            //-------------------------------------------------------------------------------
            // NOTE: Tile boxes with terrain are not intuitive. The linear combination of
            //       stress terms requires care. Create a tile box that intersects the
            //       valid box, then grow the box in x/y. Compute the strain on the local
            //       FAB over this grown tile box. Compute the stress over the tile box,
            //       except tau_ii which still needs the halo cells. Finally, write from
            //       the local FAB to the Tau MF but only on the tile box.
            //-------------------------------------------------------------------------------

            //-------------------------------------------------------------------------------
            // TODO: Avoid recomputing strain on the first RK stage. One could populate
            //       the FABs with tau_ij, compute stress, and then write to tau_ij. The
            //       problem with this approach is you will over-write the needed halo layer
            //       needed by subsequent tile boxes (particularly S_ii becomes Tau_ii).
            //-------------------------------------------------------------------------------

            // Strain/Stress tile boxes
            Box bxcc  = mfi.tilebox();
            Box tbxxy = mfi.tilebox(IntVect(1,1,0));
            Box tbxxz = mfi.tilebox(IntVect(1,0,1));
            Box tbxyz = mfi.tilebox(IntVect(0,1,1));
            // We need a halo cell for terrain
             bxcc.grow(IntVect(1,1,0));
            tbxxy.grow(IntVect(1,1,0));
            tbxxz.grow(IntVect(1,1,0));
            tbxyz.grow(IntVect(1,1,0));

            // Expansion rate
            Array4<Real> er_arr = expr->array(mfi);

            // Temporary storage for tiling/OMP
            FArrayBox S11,S22,S33;
            FArrayBox S12,S13,S23;
            S11.resize(bxcc,1);  S22.resize(bxcc,1);  S33.resize(bxcc,1);
            S12.resize(tbxxy,1); S13.resize(tbxxz,1); S23.resize(tbxyz,1);
            Elixir S11_eli   = S11.elixir(); Elixir S22_eli   = S22.elixir(); Elixir S33_eli   = S33.elixir();
            Elixir S12_eli   = S12.elixir(); Elixir S13_eli   = S13.elixir(); Elixir S23_eli   = S23.elixir();
            Array4<Real> s11 = S11.array();  Array4<Real> s22 = S22.array();  Array4<Real> s33 = S33.array();
            Array4<Real> s12 = S12.array();  Array4<Real> s13 = S13.array();  Array4<Real> s23 = S23.array();

            // Symmetric strain/stresses
            Array4<Real> tau11 = Tau11->array(mfi); Array4<Real> tau22 = Tau22->array(mfi); Array4<Real> tau33 = Tau33->array(mfi);
            Array4<Real> tau12 = Tau12->array(mfi); Array4<Real> tau13 = Tau13->array(mfi); Array4<Real> tau23 = Tau23->array(mfi);

            // Strain magnitude
            Array4<Real> SmnSmn_a;

            if (l_use_terrain) {
                // Terrain non-symmetric terms
                FArrayBox S21,S31,S32;
                S21.resize(tbxxy,1); S31.resize(tbxxz,1); S32.resize(tbxyz,1);
                Elixir S21_eli     = S21.elixir();      Elixir S31_eli     = S31.elixir();      Elixir S32_eli     = S32.elixir();
                Array4<Real> s21   = S21.array();       Array4<Real> s31   = S31.array();       Array4<Real> s32   = S32.array();
                Array4<Real> tau21 = Tau21->array(mfi); Array4<Real> tau31 = Tau31->array(mfi); Array4<Real> tau32 = Tau32->array(mfi);

                //-----------------------------------------
                // Expansion rate compute terrain
                //-----------------------------------------
                {
                BL_PROFILE("slow_rhs_making_er_T");
                // First create Omega using velocity (not momentum)
                Box gbxo = surroundingNodes(bxcc,2);
                ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    omega_arr(i,j,k) = (k == 0) ? 0. : OmegaFromW(i,j,k,w(i,j,k),u,v,z_nd,dxInv);
                });

                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {

                    Real met_u_h_zeta_hi = Compute_h_zeta_AtIface(i+1, j  , k, dxInv, z_nd);
                    Real met_u_h_zeta_lo = Compute_h_zeta_AtIface(i  , j  , k, dxInv, z_nd);

                    Real met_v_h_zeta_hi = Compute_h_zeta_AtJface(i  , j+1, k, dxInv, z_nd);
                    Real met_v_h_zeta_lo = Compute_h_zeta_AtJface(i  , j  , k, dxInv, z_nd);

                    Real Omega_hi = omega_arr(i,j,k+1);
                    Real Omega_lo = omega_arr(i,j,k  );

                    Real mfsq = mf_m(i,j,0)*mf_m(i,j,0);

                    Real expansionRate = (u(i+1,j  ,k)/mf_u(i+1,j,0)*met_u_h_zeta_hi - u(i,j,k)/mf_u(i,j,0)*met_u_h_zeta_lo)*dxInv[0]*mfsq +
                                         (v(i  ,j+1,k)/mf_v(i,j+1,0)*met_v_h_zeta_hi - v(i,j,k)/mf_v(i,j,0)*met_v_h_zeta_lo)*dxInv[1]*mfsq +
                                         (Omega_hi - Omega_lo)*dxInv[2];

                    er_arr(i,j,k) = expansionRate / detJ_arr(i,j,k);
                });
                } // end profile

                //-----------------------------------------
                // Strain tensor compute terrain
                //-----------------------------------------
                {
                BL_PROFILE("slow_rhs_making_strain_T");
                ComputeStrain_T(bxcc, tbxxy, tbxxz, tbxyz,
                                u, v, w,
                                s11, s22, s33,
                                s12, s13,
                                s21, s23,
                                s31, s32,
                                z_nd, bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
                } // profile

                // Populate SmnSmn if using Deardorff (used as diff src in post)
                // and in the first RK stage (TKE tendencies constant for nrk>0, following WRF)
                if ((nrk==0) && (tc.les_type == LESType::Deardorff)) {
                    SmnSmn_a = SmnSmn->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,s11,s22,s33,s12,s13,s23);
                    });
                }

                //-----------------------------------------
                // Stress tensor compute terrain
                //-----------------------------------------
                {
                BL_PROFILE("slow_rhs_making_stress_T");

                // Remove Halo cells just for tau_ij comps
                tbxxy.grow(IntVect(-1,-1,0));
                tbxxz.grow(IntVect(-1,-1,0));
                tbxyz.grow(IntVect(-1,-1,0));

                Real mu_eff = 2.0 * dc.dynamicViscosity; // Initialized to 0
                if (!l_use_turb) {
                    ComputeStressConsVisc_T(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                            s11, s22, s33,
                                            s12, s13,
                                            s21, s23,
                                            s31, s32,
                                            er_arr, z_nd, dxInv);
                } else {
                    ComputeStressVarVisc_T(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           s11, s22, s33,
                                           s12, s13,
                                           s21, s23,
                                           s31, s32,
                                           er_arr, z_nd, dxInv);
                }

                // Remove halo cells from tau_ii but extend across valid_box bdry
                bxcc.grow(IntVect(-1,-1,0));
                if (bxcc.smallEnd(0) == valid_bx.smallEnd(0)) bxcc.growLo(0, 1);
                if (bxcc.bigEnd(0)   == valid_bx.bigEnd(0))   bxcc.growHi(0, 1);
                if (bxcc.smallEnd(1) == valid_bx.smallEnd(1)) bxcc.growLo(1, 1);
                if (bxcc.bigEnd(1)   == valid_bx.bigEnd(1))   bxcc.growHi(1, 1);

                // Copy from temp FABs back to tau
                ParallelFor(bxcc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau11(i,j,k) = s11(i,j,k);
                    tau22(i,j,k) = s22(i,j,k);
                    tau33(i,j,k) = s33(i,j,k);
                });

                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                    tau21(i,j,k) = s21(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                    tau31(i,j,k) = s31(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                    tau32(i,j,k) = s32(i,j,k);
                });
                } // end profile

            } else {

                //-----------------------------------------
                // Expansion rate compute no terrain
                //-----------------------------------------
                {
                BL_PROFILE("slow_rhs_making_er_N");
                ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    Real mfsq = mf_m(i,j,0)*mf_m(i,j,0);
                    er_arr(i,j,k) = (u(i+1, j  , k  )/mf_u(i+1,j,0) - u(i, j, k)/mf_u(i,j,0))*dxInv[0]*mfsq +
                                    (v(i  , j+1, k  )/mf_v(i,j+1,0) - v(i, j, k)/mf_v(i,j,0))*dxInv[1]*mfsq +
                                    (w(i  , j  , k+1) - w(i, j, k))*dxInv[2];
                });
                } // end profile


                //-----------------------------------------
                // Strain tensor compute no terrain
                //-----------------------------------------
                {
                BL_PROFILE("slow_rhs_making_strain_N");
                ComputeStrain_N(bxcc, tbxxy, tbxxz, tbxyz,
                                u, v, w,
                                s11, s22, s33,
                                s12, s13, s23,
                                bc_ptr_h, dxInv,
                                mf_m, mf_u, mf_v);
                } // end profile

                // Populate SmnSmn if using Deardorff (used as diff src in post)
                // and in the first RK stage (TKE tendencies constant for nrk>0, following WRF)
                if ((nrk==0) && (tc.les_type == LESType::Deardorff)) {
                    SmnSmn_a = SmnSmn->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        SmnSmn_a(i,j,k) = ComputeSmnSmn(i,j,k,s11,s22,s33,s12,s13,s23);
                    });
                }

                //-----------------------------------------
                // Stress tensor compute no terrain
                //-----------------------------------------
                {
                BL_PROFILE("slow_rhs_making_stress_N");

                // Remove Halo cells just for tau_ij comps
                tbxxy.grow(IntVect(-1,-1,0));
                tbxxz.grow(IntVect(-1,-1,0));
                tbxyz.grow(IntVect(-1,-1,0));

                Real mu_eff = 2.0 * dc.dynamicViscosity; // Initialized to 0
                if (!l_use_turb) {
                    ComputeStressConsVisc_N(bxcc, tbxxy, tbxxz, tbxyz, mu_eff,
                                            s11, s22, s33,
                                            s12, s13, s23,
                                            er_arr);
                } else {
                    ComputeStressVarVisc_N(bxcc, tbxxy, tbxxz, tbxyz, mu_eff, mu_turb,
                                           s11, s22, s33,
                                           s12, s13, s23,
                                           er_arr);
                }

                // Remove halo cells from tau_ii but extend across valid_box bdry
                bxcc.grow(IntVect(-1,-1,0));
                if (bxcc.smallEnd(0) == valid_bx.smallEnd(0)) bxcc.growLo(0, 1);
                if (bxcc.bigEnd(0)   == valid_bx.bigEnd(0))   bxcc.growHi(0, 1);
                if (bxcc.smallEnd(1) == valid_bx.smallEnd(1)) bxcc.growLo(1, 1);
                if (bxcc.bigEnd(1)   == valid_bx.bigEnd(1))   bxcc.growHi(1, 1);

                // Copy from temp FABs back to tau
                ParallelFor(bxcc,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau11(i,j,k) = s11(i,j,k);
                    tau22(i,j,k) = s22(i,j,k);
                    tau33(i,j,k) = s33(i,j,k);
                });
                ParallelFor(tbxxy, tbxxz, tbxyz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau12(i,j,k) = s12(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau13(i,j,k) = s13(i,j,k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    tau23(i,j,k) = s23(i,j,k);
                });
                } // end profile
            } // l_use_terrain
        } // MFIter
    } // l_use_diff

    // *************************************************************************
    // Define updates and fluxes in the current RK stage
    // *************************************************************************
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
    std::array<FArrayBox,AMREX_SPACEDIM> flux;

    for ( MFIter mfi(S_data[IntVar::cons],TileNoZ()); mfi.isValid(); ++mfi)
    {
        Box bx  = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        // We don't compute a source term for z-momentum on the bottom or top boundary
        tbz.growLo(2,-1);
        tbz.growHi(2,-1);

        const Array4<const Real> & cell_data  = S_data[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_prim  = S_prim.array(mfi);
        const Array4<Real>       & cell_rhs   = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & buoyancy_fab = buoyancy.const_array(mfi);

        Array4<Real> avg_xmom = S_scratch[IntVar::xmom].array(mfi);
        Array4<Real> avg_ymom = S_scratch[IntVar::ymom].array(mfi);
        Array4<Real> avg_zmom = S_scratch[IntVar::zmom].array(mfi);

        const Array4<const Real> & u = xvel.array(mfi);
        const Array4<const Real> & v = yvel.array(mfi);
        const Array4<const Real> & w = zvel.array(mfi);

        const Array4<const Real>& rho_u = S_data[IntVar::xmom].array(mfi);
        const Array4<const Real>& rho_v = S_data[IntVar::ymom].array(mfi);
        const Array4<const Real>& rho_w = S_data[IntVar::zmom].array(mfi);

        // Map factors
        const Array4<const Real>& mf_m   = mapfac_m->const_array(mfi);
        const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
        const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

        const Array4<      Real>& omega_arr = Omega.array(mfi);

        Array4<const Real> z_t;
        if (z_t_mf)
            z_t = z_t_mf->array(mfi);
        else
            z_t = Array4<const Real>{};

        const Array4<Real>& rho_u_rhs = S_rhs[IntVar::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVar::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVar::zmom].array(mfi);

        const Array4<Real const>& mu_turb = l_use_turb ? eddyDiffs->const_array(mfi) : Array4<const Real>{};

        // Terrain metrics
        const Array4<const Real>& z_nd     = l_use_terrain ? z_phys_nd->const_array(mfi) : Array4<const Real>{};
        const Array4<const Real>& detJ_arr = l_use_terrain ?      detJ->const_array(mfi) : Array4<const Real>{};

        // Base state
        const Array4<const Real>& p0_arr = p0->const_array(mfi);

        // *************************************************************************
        // Define flux arrays for use in advection
        // *************************************************************************
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flux[dir].resize(amrex::surroundingNodes(bx,dir),2);
            flux[dir].setVal<RunOn::Device>(0.);
        }
        const GpuArray<const Array4<Real>, AMREX_SPACEDIM>
            flx_arr{{AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};

        //-----------------------------------------
        // Perturbational pressure field
        //-----------------------------------------
        Box gbx = mfi.tilebox(); gbx.grow(IntVect(1,1,0));
        FArrayBox pprime; pprime.resize(gbx,1);
        Elixir pp_eli = pprime.elixir();
        const Array4<Real      > & pp_arr = pprime.array();
        const Array4<Real const> & qv_arr = (use_moisture) ? qv->const_array(mfi) : Array4<Real> {};
        {
        BL_PROFILE("slow_rhs_pre_pprime");
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            //if (cell_data(i,j,k,RhoTheta_comp) < 0.) printf("BAD THETA AT %d %d %d %e %e \n",
            //    i,j,k,cell_data(i,j,k,RhoTheta_comp),cell_data(i,j,k+1,RhoTheta_comp));
            AMREX_ASSERT(cell_data(i,j,k,RhoTheta_comp) > 0.);
            Real qv_for_p = 0.0;
            if (use_moisture) {
                qv_for_p = qv_arr(i,j,k);
            }
            pp_arr(i,j,k) = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp),qv_for_p) - p0_arr(i,j,k);
        });
        } // end profile

        //-----------------------------------------
        // Contravariant flux field
        //-----------------------------------------
        {
        BL_PROFILE("slow_rhs_making_omega");
            Box gbxo = surroundingNodes(bx,2); gbxo.grow(IntVect(1,1,0));
            // Now create Omega with momentum (not velocity) with z_t subtracted if moving terrain
            if (l_use_terrain) {

                Box gbxo_lo = gbxo; gbxo_lo.setBig(2,0);
                ParallelFor(gbxo_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = 0.;
                });
                Box gbxo_hi = gbxo; gbxo_hi.setSmall(2,gbxo.bigEnd(2));
                ParallelFor(gbxo_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = rho_w(i,j,k);
                });

                if (z_t) {
                    Box gbxo_mid = gbxo; gbxo_mid.setSmall(2,1); gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
                    ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        // We define rho on the z-face the same way as in MomentumToVelocity/VelocityToMomentum
                        Real rho_at_face = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp));
                        omega_arr(i,j,k) = OmegaFromW(i,j,k,rho_w(i,j,k),rho_u,rho_v,z_nd,dxInv) -
                            rho_at_face * z_t(i,j,k);
                    });
                } else {
                    Box gbxo_mid = gbxo; gbxo_mid.setSmall(2,1); gbxo_mid.setBig(2,gbxo.bigEnd(2)-1);
                    ParallelFor(gbxo_mid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        omega_arr(i,j,k) = OmegaFromW(i,j,k,rho_w(i,j,k),rho_u,rho_v,z_nd,dxInv);
                    });
                }
            } else {
                ParallelFor(gbxo, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    omega_arr(i,j,k) = rho_w(i,j,k);
                });
            }
        } // end profile


        //-----------------------------------------
        // Diffusive terms (pre-computed above)
        //-----------------------------------------
        // No terrain diffusion
        Array4<Real> tau11,tau22,tau33;
        Array4<Real> tau12,tau13,tau23;
        if (Tau11) {
            tau11 = Tau11->array(mfi); tau22 = Tau22->array(mfi); tau33 = Tau33->array(mfi);
            tau12 = Tau12->array(mfi); tau13 = Tau13->array(mfi); tau23 = Tau23->array(mfi);
        } else {
            tau11 = Array4<Real>{}; tau22 = Array4<Real>{}; tau33 = Array4<Real>{};
            tau12 = Array4<Real>{}; tau13 = Array4<Real>{}; tau23 = Array4<Real>{};
        }
        // Terrain diffusion
        Array4<Real> tau21,tau31,tau32;
        if (Tau21) {
            tau21 = Tau21->array(mfi); tau31 = Tau31->array(mfi); tau32 = Tau32->array(mfi);
        } else {
            tau21 = Array4<Real>{}; tau31 = Array4<Real>{}; tau32 = Array4<Real>{};
        }

        // Strain magnitude
        Array4<Real> SmnSmn_a;
        if (tc.les_type == LESType::Deardorff) {
            SmnSmn_a = SmnSmn->array(mfi);
        } else {
            SmnSmn_a = Array4<Real>{};
        }

        // **************************************************************************
        // Define updates in the RHS of continuity and potential temperature equations
        // **************************************************************************
        AdvectionSrcForRho(bx, cell_rhs,
                           rho_u, rho_v, omega_arr,      // these are being used to build the fluxes
                           avg_xmom, avg_ymom, avg_zmom, // these are being defined from the fluxes
                           z_nd, detJ_arr, dxInv, mf_m, mf_u, mf_v,
                           l_use_terrain, flx_arr);

        int icomp = RhoTheta_comp; int ncomp = 1;
        AdvectionSrcForScalars(bx, icomp, ncomp,
                               avg_xmom, avg_ymom, avg_zmom,
                               cell_prim, cell_rhs, detJ_arr, dxInv, mf_m,
                               l_horiz_adv_type, l_vert_adv_type, l_use_terrain, flx_arr);

        if (l_use_diff) {
            Array4<Real> diffflux_x = dflux_x->array(mfi);
            Array4<Real> diffflux_y = dflux_y->array(mfi);
            Array4<Real> diffflux_z = dflux_z->array(mfi);

            Array4<Real> hfx_z = Hfx3->array(mfi);
            Array4<Real> diss  = Diss->array(mfi);

            const Array4<const Real> tm_arr = t_mean_mf ? t_mean_mf->const_array(mfi) : Array4<const Real>{};

            // NOTE: No diffusion for continuity, so n starts at 1.
            int n_start = amrex::max(start_comp,RhoTheta_comp);
            int n_comp  = end_comp - n_start + 1;

            if (l_use_terrain) {
                DiffusionSrcForState_T(bx, domain, n_start, n_comp, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z, z_nd, detJ_arr,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       hfx_z, diss, mu_turb, dc, tc,
                                       tm_arr, grav_gpu, bc_ptr);
            } else {
                DiffusionSrcForState_N(bx, domain, n_start, n_comp, u, v,
                                       cell_data, cell_prim, cell_rhs,
                                       diffflux_x, diffflux_y, diffflux_z,
                                       dxInv, SmnSmn_a, mf_m, mf_u, mf_v,
                                       hfx_z, diss,
                                       mu_turb, dc, tc,
                                       tm_arr, grav_gpu, bc_ptr);
            }
        }

        if (l_use_ndiff) {
            NumericalDiffusion(bx, start_comp, num_comp, dt, solverChoice.NumDiffCoeff,
                               cell_data, cell_rhs, mf_u, mf_v, false, false);
       }

        // Add source terms for (rho theta)
        {
            auto const& src_arr = source.const_array(mfi);
            if (l_use_terrain && l_moving_terrain) {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_rhs(i,j,k,RhoTheta_comp) += src_arr(i,j,k,RhoTheta_comp) / detJ_arr(i,j,k);
                });
            } else {
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    cell_rhs(i,j,k,RhoTheta_comp) += src_arr(i,j,k,RhoTheta_comp);
                });
            }
        }

        // Add Rayleigh damping
        if (solverChoice.use_rayleigh_damping && solverChoice.rayleigh_damp_T) {
            int n  = RhoTheta_comp;
            int nr = Rho_comp;
            int np = PrimTheta_comp;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real theta = cell_prim(i,j,k,np);
                cell_rhs(i, j, k, n) -= dptr_rayleigh_tau[k] * (theta - dptr_rayleigh_thetabar[k]) * cell_data(i,j,k,nr);
            });
        }

        // Multiply the slow RHS for rho and rhotheta by detJ here so we don't have to later
        if (l_use_terrain && l_moving_terrain) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cell_rhs(i,j,k,Rho_comp)      *= detJ_arr(i,j,k);
                cell_rhs(i,j,k,RhoTheta_comp) *= detJ_arr(i,j,k);
            });
        }

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        AdvectionSrcForMom(tbx, tby, tbz,
                           rho_u_rhs, rho_v_rhs, rho_w_rhs, u, v, w,
                           rho_u    , rho_v    , omega_arr,
                           z_nd, detJ_arr, dxInv, mf_m, mf_u, mf_v,
                           l_horiz_adv_type, l_vert_adv_type, l_use_terrain, domhi_z);

        if (l_use_diff) {
            if (l_use_terrain) {
                DiffusionSrcForMom_T(tbx, tby, tbz,
                                     rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                     tau11, tau22, tau33,
                                     tau12, tau13,
                                     tau21, tau23,
                                     tau31, tau32,
                                     cell_data, detJ_arr, dc, dxInv,
                                     mf_m, mf_u, mf_v);
            } else {
                DiffusionSrcForMom_N(tbx, tby, tbz,
                                     rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                     tau11, tau22, tau33,
                                     tau12, tau13, tau23,
                                     cell_data, dc, dxInv,
                                     mf_m, mf_u, mf_v);
            }
        }

        if (l_use_ndiff) {
            NumericalDiffusion(tbx, 0, 1, dt, solverChoice.NumDiffCoeff,
                               rho_u, rho_u_rhs, mf_m, mf_v, false, true);
            NumericalDiffusion(tby, 0, 1, dt, solverChoice.NumDiffCoeff,
                               rho_v, rho_v_rhs, mf_u, mf_m, true, false);
            NumericalDiffusion(tbz, 0, 1, dt, solverChoice.NumDiffCoeff,
                               rho_w, rho_w_rhs, mf_u, mf_v, false, false);
        }

        auto abl_pressure_grad    = solverChoice.abl_pressure_grad;
        auto abl_geo_forcing      = solverChoice.abl_geo_forcing;
        auto use_coriolis         = solverChoice.use_coriolis;
        auto coriolis_factor      = solverChoice.coriolis_factor;
        auto cosphi               = solverChoice.cosphi;
        auto sinphi               = solverChoice.sinphi;
        auto use_rayleigh_damping = solverChoice.use_rayleigh_damping;

        {
        BL_PROFILE("slow_rhs_pre_xmom");
        auto rayleigh_damp_U      = solverChoice.rayleigh_damp_U;
        // ******************************************************************
        // TERRAIN VERSION
        // ******************************************************************
        if (l_use_terrain) {
          ParallelFor(tbx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          { // x-momentum equation

            Real rho_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
            Real met_h_xi   = Compute_h_xi_AtIface  (i, j, k, dxInv, z_nd);
            Real met_h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);

            //Note : mx/my == 1, so no map factor needed here
            Real gp_xi = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));
            Real gp_zeta_on_iface;
            if (k==0) {
                gp_zeta_on_iface = 0.5 * dxInv[2] * (
                                                    pp_arr(i-1,j,k+1) + pp_arr(i,j,k+1)
                                                  - pp_arr(i-1,j,k  ) - pp_arr(i,j,k  ) );
            } else if (k==domhi_z) {
                gp_zeta_on_iface = 0.5 * dxInv[2] * (
                                                    pp_arr(i-1,j,k  ) + pp_arr(i,j,k  )
                                                  - pp_arr(i-1,j,k-1) - pp_arr(i,j,k-1) );
            } else {
                gp_zeta_on_iface = 0.25 * dxInv[2] * (
                                                     pp_arr(i-1,j,k+1) + pp_arr(i,j,k+1)
                                                   - pp_arr(i-1,j,k-1) - pp_arr(i,j,k-1) );
            }
            Real gpx = gp_xi - (met_h_xi/ met_h_zeta) * gp_zeta_on_iface;
            gpx *= mf_u(i,j,0);

            Real q = 0.0;
            if (use_moisture) {
                q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i-1,j,k,PrimQ1_comp)
                           +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i-1,j,k,PrimQ2_comp) );
            }

            rho_u_rhs(i, j, k) += (-gpx - abl_pressure_grad[0]) / (1.0 + q)
                                  + rho_u_face * abl_geo_forcing[0];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (use_coriolis)
            {
                Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                rho_u_rhs(i, j, k) += coriolis_factor * (rho_v_loc * sinphi - rho_w_loc * cosphi);
            }

            // Add Rayleigh damping
            if (use_rayleigh_damping && rayleigh_damp_U)
            {
                Real uu = rho_u(i,j,k) / rho_u_face;
                rho_u_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (uu - dptr_rayleigh_ubar[k]) * rho_u_face;
            }

            if (l_moving_terrain) {
                Real h_zeta = Compute_h_zeta_AtIface(i, j, k, dxInv, z_nd);
                rho_u_rhs(i, j, k) *= h_zeta;
            }
        });

        } else {
        // ******************************************************************
        // NON-TERRAIN VERSION
        // ******************************************************************
          ParallelFor(tbx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          { // x-momentum equation

              Real rho_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
              Real gpx = dxInv[0] * (pp_arr(i,j,k) - pp_arr(i-1,j,k));
              gpx *= mf_u(i,j,0);

              Real q = 0.0;
              if (use_moisture) {
                  q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i-1,j,k,PrimQ1_comp)
                             +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i-1,j,k,PrimQ2_comp) );
              }

              rho_u_rhs(i, j, k) += (-gpx - abl_pressure_grad[0]) / (1.0 + q)
                                    + rho_u_face * abl_geo_forcing[0];

              // Add Coriolis forcing (that assumes east is +x, north is +y)
              if (use_coriolis)
              {
                  Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                  Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                  rho_u_rhs(i, j, k) += coriolis_factor * (rho_v_loc * sinphi - rho_w_loc * cosphi);
              }

              // Add Rayleigh damping
              if (use_rayleigh_damping && rayleigh_damp_U)
              {
                  Real uu = rho_u(i,j,k) / cell_data(i,j,k,Rho_comp);
                  rho_u_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (uu - dptr_rayleigh_ubar[k]) * cell_data(i,j,k,Rho_comp);
              }
          });
        } // no terrain
        } // end profile

        {
        BL_PROFILE("slow_rhs_pre_ymom");
        auto rayleigh_damp_V      = solverChoice.rayleigh_damp_V;
        // ******************************************************************
        // TERRAIN VERSION
        // ******************************************************************
        if (l_use_terrain) {
          ParallelFor(tby,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          { // y-momentum equation

              Real rho_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
              Real met_h_eta  = Compute_h_eta_AtJface (i, j, k, dxInv, z_nd);
              Real met_h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);

              //Note : mx/my == 1, so no map factor needed here
              Real gp_eta = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));
              Real gp_zeta_on_jface;
              if(k==0) {
                  gp_zeta_on_jface = 0.5 * dxInv[2] * (
                                                      pp_arr(i,j,k+1) + pp_arr(i,j-1,k+1)
                                                    - pp_arr(i,j,k  ) - pp_arr(i,j-1,k  ) );
              } else if (k==domhi_z) {
                  gp_zeta_on_jface = 0.5 * dxInv[2] * (
                                                      pp_arr(i,j,k  ) + pp_arr(i,j-1,k  )
                                                    - pp_arr(i,j,k-1) - pp_arr(i,j-1,k-1) );
              } else {
                  gp_zeta_on_jface = 0.25 * dxInv[2] * (
                                                       pp_arr(i,j,k+1) + pp_arr(i,j-1,k+1)
                                                     - pp_arr(i,j,k-1) - pp_arr(i,j-1,k-1) );
              }

              Real gpy = gp_eta - (met_h_eta / met_h_zeta) * gp_zeta_on_jface;
              gpy *= mf_v(i,j,0);

              Real q = 0.0;
              if (use_moisture) {
                  q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i,j-1,k,PrimQ1_comp)
                             +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i,j-1,k,PrimQ2_comp) );
              }
              rho_v_rhs(i, j, k) += (-gpy - abl_pressure_grad[1]) / (1.0_rt + q)
                                    + rho_v_face * abl_geo_forcing[1];

              // Add Coriolis forcing (that assumes east is +x, north is +y) if (use_coriolis)
              {
                  Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                  rho_v_rhs(i, j, k) += -coriolis_factor * rho_u_loc * sinphi;
              }

              // Add Rayleigh damping
              if (use_rayleigh_damping && rayleigh_damp_V)
              {
                  Real vv = rho_v(i,j,k) / rho_v_face;
                  rho_v_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (vv - dptr_rayleigh_vbar[k]) * rho_v_face;
              }

              if (l_moving_terrain) {
                  Real h_zeta = Compute_h_zeta_AtJface(i, j, k, dxInv, z_nd);
                  rho_v_rhs(i, j, k) *= h_zeta;
              }
          });

        // ******************************************************************
        // NON-TERRAIN VERSION
        // ******************************************************************
        } else {
          ParallelFor(tby,
          [=] AMREX_GPU_DEVICE (int i, int j, int k)
          { // y-momentum equation

              Real rho_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
              Real gpy = dxInv[1] * (pp_arr(i,j,k) - pp_arr(i,j-1,k));
              gpy *= mf_v(i,j,0);

              Real q = 0.0;
              if (use_moisture) {
                  q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i,j-1,k,PrimQ1_comp)
                             +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i,j-1,k,PrimQ2_comp) );
              }

              rho_v_rhs(i, j, k) += (-gpy - abl_pressure_grad[1]) / (1.0_rt + q)
                                    + rho_v_face * abl_geo_forcing[1];

              // Add Coriolis forcing (that assumes east is +x, north is +y)
              if (use_coriolis)
              {
                  Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                  rho_v_rhs(i, j, k) += -coriolis_factor * rho_u_loc * sinphi;
              }

              // Add Rayleigh damping
              if (use_rayleigh_damping && rayleigh_damp_V)
              {
                  Real vv = rho_v(i,j,k) / rho_v_face;
                  rho_v_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (vv - dptr_rayleigh_vbar[k]) * rho_v_face;
              }
          });
        } // no terrain
        } // end profile

        {
        BL_PROFILE("slow_rhs_pre_zmom_2d");
        amrex::Box b2d = tbz;
        b2d.setSmall(2,0);
        b2d.setBig(2,0);
        // Enforce no forcing term at top and bottom boundaries
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) {
            rho_w_rhs(i,j,        0) = 0.;
            rho_w_rhs(i,j,domhi_z+1) = 0.; // TODO: generalize this
        });
        } // end profile

        {
        BL_PROFILE("slow_rhs_pre_zmom");
        auto rayleigh_damp_W      = solverChoice.rayleigh_damp_W;

        // ******************************************************************
        // TERRAIN VERSION
        // ******************************************************************
        if (l_use_terrain) {
          ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation

                Real rho_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                Real met_h_zeta = Compute_h_zeta_AtKface(i, j, k, dxInv, z_nd);
                Real gpz = dxInv[2] * ( pp_arr(i,j,k)-pp_arr(i,j,k-1) )  / met_h_zeta;

                Real q = 0.0;
                if (use_moisture) {
                    q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i,j,k-1,PrimQ1_comp)
                               +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i,j,k-1,PrimQ2_comp) );
                }
                rho_w_rhs(i, j, k) += (buoyancy_fab(i,j,k) - gpz - abl_pressure_grad[2]) / (1.0_rt + q)
                                     + rho_w_face * abl_geo_forcing[2];

                // Add Coriolis forcing (that assumes east is +x, north is +y)
                if (use_coriolis)
                {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                    rho_w_rhs(i, j, k) += coriolis_factor * rho_u_loc * cosphi;
                }

                // Add Rayleigh damping
                if (use_rayleigh_damping && rayleigh_damp_W)
                {
                    Real ww = rho_w(i,j,k) / rho_w_face;
                    rho_w_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (ww - dptr_rayleigh_wbar[k]) * rho_w_face;
                }

                if (l_use_terrain && l_moving_terrain) {
                     rho_w_rhs(i, j, k) *= 0.5 * (detJ_arr(i,j,k) + detJ_arr(i,j,k-1));
                }
          });

        // ******************************************************************
        // NON-TERRAIN VERSION
        // ******************************************************************
        } else {
          ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
          { // z-momentum equation

                Real rho_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                Real gpz = dxInv[2] * ( pp_arr(i,j,k)-pp_arr(i,j,k-1) );

                Real q = 0.0;
                if (use_moisture) {
                    q = 0.5 * ( cell_prim(i,j,k,PrimQ1_comp) + cell_prim(i,j,k-1,PrimQ1_comp)
                               +cell_prim(i,j,k,PrimQ2_comp) + cell_prim(i,j,k-1,PrimQ2_comp) );
                }
                rho_w_rhs(i, j, k) += (buoyancy_fab(i,j,k) - gpz - abl_pressure_grad[2]) / (1.0_rt + q)
                                     + rho_w_face * abl_geo_forcing[2];

                // Add Coriolis forcing (that assumes east is +x, north is +y)
                if (use_coriolis)
                {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                    rho_w_rhs(i, j, k) += coriolis_factor * rho_u_loc * cosphi;
                }

                // Add Rayleigh damping
                if (use_rayleigh_damping && rayleigh_damp_W)
                {
                    Real ww = rho_w(i,j,k) / rho_w_face;
                    rho_w_rhs(i, j, k) -= dptr_rayleigh_tau[k] * (ww - dptr_rayleigh_wbar[k]) * rho_w_face;
                }
        });
        } // no terrain
            ApplySpongeZoneBCs(solverChoice.spongeChoice, geom, tbx, tby, tbz, rho_u_rhs, rho_v_rhs, rho_w_rhs, rho_u, rho_v,
                               rho_w, bx, cell_rhs, cell_data);
        } // end profile

        {
        BL_PROFILE("slow_rhs_pre_fluxreg");
        // We only add to the flux registers in the final RK step
        if (l_reflux && nrk == 2) {
            int strt_comp_reflux = 0;
            int  num_comp_reflux = 2;
            if (level < finest_level) {
                fr_as_crse->CrseAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dt, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, amrex::RunOn::Device);
            }
            if (level > 0) {
                fr_as_fine->FineAdd(mfi,
                    {{AMREX_D_DECL(&(flux[0]), &(flux[1]), &(flux[2]))}},
                    dx, dt, strt_comp_reflux, strt_comp_reflux, num_comp_reflux, amrex::RunOn::Device);
            }
        } // two-way coupling
        } // end profile
    } // mfi
    } // OMP
}
