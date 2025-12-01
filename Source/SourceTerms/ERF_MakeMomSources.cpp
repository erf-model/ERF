#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_TableData.H>
#include <AMReX_GpuContainers.H>

#include "ERF_NumericalDiffusion.H"
#include "ERF_PlaneAverage.H"
#include "ERF_TI_slow_headers.H"
#include "ERF_SrcHeaders.H"
#include "ERF_Utils.H"

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in] time current time
 * @param[in] dt current slow or fast timestep size
 * @param[in]  S_data current solution
 * @param[in]  xvel x-component of velocity
 * @param[in]  yvel y-component of velocity
 * @param[in] xmom_src source terms for x-momentum
 * @param[in] ymom_src source terms for y-momentum
 * @param[in] zmom_src source terms for z-momentum
 * @param[in]  geom   Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in] mapfac map factors
 * @param[in] dptr_u_geos  custom geostrophic wind profile
 * @param[in] dptr_v_geos  custom geostrophic wind profile
 * @param[in] dptr_wbar_sub  subsidence source term
 * @param[in] d_rayleigh_ptrs_at_lev  Vector of {strength of Rayleigh damping, reference value for xvel/yvel/zvel/theta} used to define Rayleigh damping
 * @param[in] d_sinesq_at_lev  sin( (pi/2) (z-z_t)/(damping depth)) at cell centers
 * @param[in] d_sinesq_stag_at_lev  sin( (pi/2) (z-z_t)/(damping depth)) at z-faces
 */

void make_mom_sources (Real time,
                       Real /*dt*/,
                       const Vector<MultiFab>& S_data,
                       const MultiFab* z_phys_nd,
                       const MultiFab* z_phys_cc,
                             Vector<Real>& stretched_dz_h,
                       const MultiFab& xvel,
                       const MultiFab& yvel,
                       const MultiFab& wvel,
                             MultiFab& xmom_src,
                             MultiFab& ymom_src,
                             MultiFab& zmom_src,
                       const MultiFab& base_state,
                             MultiFab* forest_drag,
                             MultiFab* terrain_blank,
                             MultiFab* cosPhi_mf,
                             MultiFab* sinPhi_mf,
                       const Geometry geom,
                       const SolverChoice& solverChoice,
                             Vector<std::unique_ptr<MultiFab>>& /*mapfac*/,
                       const Real* dptr_u_geos,
                       const Real* dptr_v_geos,
                       const Real* dptr_wbar_sub,
                       const Vector<Real*> d_rayleigh_ptrs_at_lev,
                       const amrex::Real* d_sinesq_at_lev,
                       const amrex::Real* d_sinesq_stag_at_lev,
                       const Vector<Real*> d_sponge_ptrs_at_lev,
                       const Vector<MultiFab>* forecast_state_at_lev,
                             InputSoundingData& input_sounding_data,
                             bool is_slow_step)
{
    BL_PROFILE_REGION("erf_make_mom_sources()");

    Box domain(geom.Domain());
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // Initialize sources to zero each time we may use them
    xmom_src.setVal(0.0);
    ymom_src.setVal(0.0);
    zmom_src.setVal(0.0);

    MultiFab r_hse (base_state, make_alias, BaseState::r0_comp , 1);

    // flags to apply certain source terms in substep call only
    bool use_Rayleigh_fast_uv = ( (solverChoice.dampingChoice.rayleigh_damping_type == RayleighDampingType::FastExplicit) ||
                                  (solverChoice.dampingChoice.rayleigh_damping_type == RayleighDampingType::FastImplicit) );
    bool use_Rayleigh_fast_w  = (solverChoice.dampingChoice.rayleigh_damping_type == RayleighDampingType::FastExplicit);
    bool use_canopy_fast = solverChoice.forest_substep;
    bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;

    // *****************************************************************************
    // Define source term for all three components of momenta from
    //    1. Coriolis forcing for (xmom,ymom,zmom)
    //    2. Rayleigh damping for (xmom,ymom,zmom)
    //    3. Constant / height-dependent geostrophic forcing
    //    4. Subsidence
    //    5. Nudging towards input sounding data
    //    6. Numerical diffusion for (xmom,ymom,zmom)
    //    7. Sponge
    //    8. Forest canopy
    //    9a. Immersed forcing for terrain
    //    9b. Immersed forcing for buildings
    //   10. Constant mass flux
    // *****************************************************************************
    // NOTE: buoyancy is now computed in a separate routine - it should not appear here
    // *****************************************************************************
    //const bool l_use_ndiff       = solverChoice.use_num_diff;

    if (solverChoice.terrain_type == TerrainType::ImmersedForcing) {
        if (solverChoice.do_forest_drag) {
            amrex::Error(" Currently forest canopy cannot be used with immersed forcing");
        }
    }


    // *****************************************************************************
    // Data for Coriolis forcing
    // *****************************************************************************
    auto use_coriolis         = solverChoice.use_coriolis;
    auto coriolis_factor      = solverChoice.coriolis_factor;
    auto cosphi               = solverChoice.cosphi;
    auto sinphi               = solverChoice.sinphi;
    auto var_coriolis         = solverChoice.variable_coriolis;
    auto has_lat_lon          = solverChoice.has_lat_lon;

    // *****************************************************************************
    // Flag for Geostrophic forcing
    // *****************************************************************************
    auto abl_geo_forcing  = solverChoice.abl_geo_forcing;
    auto geo_wind_profile = solverChoice.have_geo_wind_profile;

    // *****************************************************************************
    // Data for Rayleigh damping
    // *****************************************************************************
    auto rayleigh_damp_U  = solverChoice.dampingChoice.rayleigh_damp_U;
    auto rayleigh_damp_V  = solverChoice.dampingChoice.rayleigh_damp_V;
    auto rayleigh_damp_W  = solverChoice.dampingChoice.rayleigh_damp_W;

    Real*     ubar = d_rayleigh_ptrs_at_lev[Rayleigh::ubar];
    Real*     vbar = d_rayleigh_ptrs_at_lev[Rayleigh::vbar];
    Real*     wbar = d_rayleigh_ptrs_at_lev[Rayleigh::wbar];

    // *****************************************************************************
    // Data for constant mass flux
    // *****************************************************************************
    bool enforce_massflux_x = (solverChoice.const_massflux_u != 0);
    bool enforce_massflux_y = (solverChoice.const_massflux_v != 0);
    Real U_target = solverChoice.const_massflux_u;
    Real V_target = solverChoice.const_massflux_v;
    int massflux_klo = solverChoice.massflux_klo;
    int massflux_khi = solverChoice.massflux_khi;

    // These will be updated by integrating through the planar average profiles
    Real rhoUA_target{0};
    Real rhoVA_target{0};
    Real rhoUA{0};
    Real rhoVA{0};

    // *****************************************************************************
    // Planar averages for subsidence, nudging, or constant mass flux
    // *****************************************************************************
    Table1D<Real>     dptr_r_plane, dptr_u_plane, dptr_v_plane;
    TableData<Real, 1> r_plane_tab,  u_plane_tab,  v_plane_tab;

    if (is_slow_step && (dptr_wbar_sub || solverChoice.nudging_from_input_sounding ||
                         enforce_massflux_x || enforce_massflux_y))
    {
        const int offset = 1;
        const int u_offset = 1;
        const int v_offset = 1;

        //
        // We use the alias here to control ncomp inside the PlaneAverage
        //
        MultiFab cons(S_data[IntVars::cons], make_alias, 0, 1);

        IntVect ng_c = S_data[IntVars::cons].nGrowVect(); ng_c[2] = offset;
        PlaneAverage r_ave(&cons, geom, solverChoice.ave_plane, ng_c);
        r_ave.compute_averages(ZDir(), r_ave.field());

        int ncell = r_ave.ncell_line();
        Gpu::HostVector<    Real> r_plane_h(ncell);
        Gpu::DeviceVector<  Real> r_plane_d(ncell);

        r_ave.line_average(Rho_comp, r_plane_h);

        Gpu::copyAsync(Gpu::hostToDevice, r_plane_h.begin(), r_plane_h.end(), r_plane_d.begin());

        Real* dptr_r = r_plane_d.data();

        Box tdomain  = domain; tdomain.grow(2,ng_c[2]);
        r_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});

        dptr_r_plane = r_plane_tab.table();
        ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_r_plane(k-offset) = dptr_r[k];
        });

        // U and V momentum
        IntVect ng_u = S_data[IntVars::xmom].nGrowVect(); ng_u[2] = u_offset;
        PlaneAverage u_ave(&(S_data[IntVars::xmom]), geom, solverChoice.ave_plane, ng_u);

        IntVect ng_v = S_data[IntVars::ymom].nGrowVect(); ng_v[2] = v_offset;
        PlaneAverage v_ave(&(S_data[IntVars::ymom]), geom, solverChoice.ave_plane, ng_v);

        u_ave.compute_averages(ZDir(), u_ave.field());
        v_ave.compute_averages(ZDir(), v_ave.field());

        int u_ncell = u_ave.ncell_line();
        int v_ncell = v_ave.ncell_line();
        Gpu::HostVector<    Real> u_plane_h(u_ncell), v_plane_h(v_ncell);
        Gpu::DeviceVector<  Real> u_plane_d(u_ncell), v_plane_d(v_ncell);

        u_ave.line_average(0, u_plane_h);
        v_ave.line_average(0, v_plane_h);

        Gpu::copyAsync(Gpu::hostToDevice, u_plane_h.begin(), u_plane_h.end(), u_plane_d.begin());
        Gpu::copyAsync(Gpu::hostToDevice, v_plane_h.begin(), v_plane_h.end(), v_plane_d.begin());

        Real* dptr_u = u_plane_d.data();
        Real* dptr_v = v_plane_d.data();

        Box udomain = domain; udomain.grow(2,ng_u[2]);
        Box vdomain = domain; vdomain.grow(2,ng_v[2]);
        u_plane_tab.resize({udomain.smallEnd(2)}, {udomain.bigEnd(2)});
        v_plane_tab.resize({vdomain.smallEnd(2)}, {vdomain.bigEnd(2)});

        dptr_u_plane = u_plane_tab.table();
        ParallelFor(u_ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_u_plane(k-u_offset) = dptr_u[k];
        });

        dptr_v_plane = v_plane_tab.table();
        ParallelFor(v_ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_v_plane(k-v_offset) = dptr_v[k];
        });

        // sum in z for massflux adjustment
        if (enforce_massflux_x || enforce_massflux_y) {
            Real Lx = geom.ProbHi(0) - geom.ProbLo(0);
            Real Ly = geom.ProbHi(1) - geom.ProbLo(1);

            if (solverChoice.mesh_type == MeshType::ConstantDz) {
                // note: massflux_khi corresponds to unstaggered indices in this case
                rhoUA        = std::accumulate(u_plane_h.begin() + u_offset + massflux_klo,
                                               u_plane_h.begin() + u_offset + massflux_khi+1, 0.0);
                rhoVA        = std::accumulate(v_plane_h.begin() + v_offset + massflux_klo,
                                               v_plane_h.begin() + v_offset + massflux_khi+1, 0.0);
                rhoUA_target = std::accumulate(r_plane_h.begin() +   offset + massflux_klo,
                                               r_plane_h.begin() +   offset + massflux_khi+1, 0.0);
                rhoVA_target = rhoUA_target;

                rhoUA        *= geom.CellSize(2) * Ly;
                rhoVA        *= geom.CellSize(2) * Lx;
                rhoUA_target *= geom.CellSize(2) * Ly;
                rhoVA_target *= geom.CellSize(2) * Lx;

            } else if (solverChoice.mesh_type == MeshType::StretchedDz) {
                // note: massflux_khi corresponds to staggered indices in this case
                for (int k=massflux_klo; k < massflux_khi; ++k) {
                    rhoUA        += u_plane_h[k + u_offset] * stretched_dz_h[k];
                    rhoVA        += v_plane_h[k + v_offset] * stretched_dz_h[k];
                    rhoUA_target += r_plane_h[k +   offset] * stretched_dz_h[k];
                }
                rhoVA_target = rhoUA_target;

                rhoUA        *= Ly;
                rhoVA        *= Lx;
                rhoUA_target *= Ly;
                rhoVA_target *= Lx;
            }

            // at this point, this is integrated rho*dA
            rhoUA_target *= U_target;
            rhoVA_target *= V_target;

            Print() << "Integrated mass flux : " << rhoUA << " " << rhoVA
                << " (target: " << rhoUA_target << " " << rhoVA_target << ")"
                << std::endl;
        }
    }

    // *****************************************************************************
    // Add all the other forcings
    // *****************************************************************************
    for ( MFIter mfi(S_data[IntVars::cons]); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);
        if (tbz.bigEnd(2) == domain.bigEnd(2)+1) tbz.growHi(2,-1);

        const Array4<const Real>& cell_data = S_data[IntVars::cons].array(mfi);
        const Array4<const Real>&     rho_u = S_data[IntVars::xmom].array(mfi);
        const Array4<const Real>&     rho_v = S_data[IntVars::ymom].array(mfi);
        const Array4<const Real>&     rho_w = S_data[IntVars::zmom].array(mfi);

        const Array4<const Real>& u = xvel.array(mfi);
        const Array4<const Real>& v = yvel.array(mfi);
        const Array4<const Real>& w = wvel.array(mfi);

        const Array4<      Real>& xmom_src_arr = xmom_src.array(mfi);
        const Array4<      Real>& ymom_src_arr = ymom_src.array(mfi);
        const Array4<      Real>& zmom_src_arr = zmom_src.array(mfi);

        const Array4<const Real>& r0 = r_hse.const_array(mfi);

        const Array4<const Real>& f_drag_arr = (forest_drag) ? forest_drag->const_array(mfi) :
                                                               Array4<const Real>{};
        const Array4<const Real>& t_blank_arr = (terrain_blank) ? terrain_blank->const_array(mfi) :
                                                               Array4<const Real>{};

        const Array4<const Real>& cphi_arr = (cosPhi_mf) ? cosPhi_mf->const_array(mfi) :
                                                           Array4<const Real>{};
        const Array4<const Real>& sphi_arr = (sinPhi_mf) ? sinPhi_mf->const_array(mfi) :
                                                           Array4<const Real>{};

        const Array4<const Real>& z_nd_arr =  z_phys_nd->const_array(mfi);
        const Array4<const Real>& z_cc_arr =  z_phys_cc->const_array(mfi);


        // *****************************************************************************
        // 1. Add CORIOLIS forcing (this assumes east is +x, north is +y)
        // *****************************************************************************
        if (use_coriolis && is_slow_step) {
            if(solverChoice.init_type == InitType::HindCast) {
                const Array4<const Real>& latlon_arr = (*forecast_state_at_lev)[4].array(mfi);
                ParallelFor(tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                    Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                    Real latitude = latlon_arr(i,j,k,0);
                    Real sphi_loc = std::sin(latitude*PI/180.0);
                    Real cphi_loc = std::cos(latitude*PI/180.0);
                    xmom_src_arr(i, j, k) += coriolis_factor * (rho_v_loc * sphi_loc - rho_w_loc * cphi_loc);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                    Real latitude = latlon_arr(i,j,k,0);
                    Real sphi_loc = std::sin(latitude*PI/180.0);
                    ymom_src_arr(i, j, k) += -coriolis_factor * rho_u_loc * sphi_loc;
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                    Real latitude = latlon_arr(i,j,k,0);
                    Real cphi_loc = std::cos(latitude*PI/180.0);
                    zmom_src_arr(i, j, k) += coriolis_factor * rho_u_loc * cphi_loc;
                });
            }
            else if (var_coriolis && has_lat_lon) {
                ParallelFor(tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                    Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                    Real sphi_loc  = 0.5  * (sphi_arr(i,j,0) + sphi_arr(i-1,j,0));
                    Real cphi_loc  = 0.5  * (cphi_arr(i,j,0) + cphi_arr(i-1,j,0));
                    xmom_src_arr(i, j, k) += coriolis_factor * (rho_v_loc * sphi_loc - rho_w_loc * cphi_loc);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                    Real sphi_loc  = 0.5  * (sphi_arr(i,j,0) + sphi_arr(i,j-1,0));
                    ymom_src_arr(i, j, k) += -coriolis_factor * rho_u_loc * sphi_loc;
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                    zmom_src_arr(i, j, k) += coriolis_factor * rho_u_loc * cphi_arr(i,j,0);
                });
            } else {
                ParallelFor(tbx, tby, tbz,
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rho_v_loc = 0.25 * (rho_v(i,j+1,k) + rho_v(i,j,k) + rho_v(i-1,j+1,k) + rho_v(i-1,j,k));
                    Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                    xmom_src_arr(i, j, k) += coriolis_factor * (rho_v_loc * sinphi - rho_w_loc * cosphi);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                    ymom_src_arr(i, j, k) += -coriolis_factor * rho_u_loc * sinphi;
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                    zmom_src_arr(i, j, k) += coriolis_factor * rho_u_loc * cosphi;
                });
            } // var_coriolis
        } // use_coriolis

        // *****************************************************************************
        // 2. Add RAYLEIGH damping
        // *****************************************************************************
        Real dampcoef = solverChoice.dampingChoice.rayleigh_dampcoef;

        if ( (is_slow_step && !use_Rayleigh_fast_uv) || (!is_slow_step && use_Rayleigh_fast_uv)) {
            if (rayleigh_damp_U) {
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                    Real uu = rho_u(i,j,k) / rho_on_u_face;
                    Real sinesq = d_sinesq_at_lev[k];
                    xmom_src_arr(i, j, k) -= dampcoef*sinesq * (uu - ubar[k]) * rho_on_u_face;
                });
            }

            if (rayleigh_damp_V) {
                ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                    Real vv = rho_v(i,j,k) / rho_on_v_face;
                    Real sinesq = d_sinesq_at_lev[k];
                    ymom_src_arr(i, j, k) -= dampcoef*sinesq * (vv - vbar[k]) * rho_on_v_face;
                });
            }
        } // fast or slow step

        if ( (is_slow_step && !use_Rayleigh_fast_w) || (!is_slow_step && use_Rayleigh_fast_w)) {
            if (rayleigh_damp_W) {
                    ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real rho_on_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                    Real ww = rho_w(i,j,k) / rho_on_w_face;
                    Real sinesq = d_sinesq_stag_at_lev[k];
                    zmom_src_arr(i, j, k) -= dampcoef*sinesq * (ww - wbar[k]) * rho_on_w_face;
                });
            }
        } // fast or slow step

        // *****************************************************************************
        // 3a. Add constant GEOSTROPHIC forcing
        // *****************************************************************************
        if (is_slow_step) {
            ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                xmom_src_arr(i, j, k) += rho_on_u_face * abl_geo_forcing[0];
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                ymom_src_arr(i, j, k) += rho_on_v_face * abl_geo_forcing[1];
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                zmom_src_arr(i, j, k) += rho_on_w_face * abl_geo_forcing[2];
            });
        }

        // *****************************************************************************
        // 3b. Add height-dependent GEOSTROPHIC forcing
        // *****************************************************************************
        if (geo_wind_profile && is_slow_step) {
            ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                xmom_src_arr(i, j, k) -= coriolis_factor * rho_on_u_face * dptr_v_geos[k] * sinphi;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                ymom_src_arr(i, j, k) += coriolis_factor * rho_on_v_face * dptr_u_geos[k] * sinphi;
            });
        } // geo_wind_profile

        // *****************************************************************************
        // 4. Add custom SUBSIDENCE terms
        // *****************************************************************************
        if (solverChoice.custom_w_subsidence && is_slow_step) {
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(tbx, tby,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = 0.5*dxInv[2];
                    if (z_nd_arr) {
                        Real z_xf_lo = 0.25 * ( z_nd_arr(i,j,k  ) + z_nd_arr(i,j+1,k  )
                                              + z_nd_arr(i,j,k-1) + z_nd_arr(i,j+1,k-1) );
                        Real z_xf_hi = 0.25 * ( z_nd_arr(i,j,k+1) + z_nd_arr(i,j+1,k+1)
                                              + z_nd_arr(i,j,k+2) + z_nd_arr(i,j+1,k+2) );
                        dzInv = 1.0 / (z_xf_hi - z_xf_lo);
                    }
                    Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,nr) + cell_data(i-1,j,k,nr) );
                    Real U_hi = dptr_u_plane(k+1) / dptr_r_plane(k+1);
                    Real U_lo = dptr_u_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_xf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    xmom_src_arr(i, j, k) -= rho_on_u_face * wbar_xf * (U_hi - U_lo) * dzInv;
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = 0.5*dxInv[2];
                    if (z_nd_arr) {
                        Real z_yf_lo = 0.25 * ( z_nd_arr(i,j,k  ) + z_nd_arr(i+1,j,k  )
                                              + z_nd_arr(i,j,k-1) + z_nd_arr(i+1,j,k-1) );
                        Real z_yf_hi = 0.25 * ( z_nd_arr(i,j,k+1) + z_nd_arr(i+1,j,k+1)
                                              + z_nd_arr(i,j,k+2) + z_nd_arr(i+1,j,k+2) );
                        dzInv = 1.0 / (z_yf_hi - z_yf_lo);
                    }
                    Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,nr) + cell_data(i,j-1,k,nr) );
                    Real V_hi = dptr_v_plane(k+1) / dptr_r_plane(k+1);
                    Real V_lo = dptr_v_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_yf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    ymom_src_arr(i, j, k) -= rho_on_v_face * wbar_yf * (V_hi - V_lo) * dzInv;
                });
            } else {
                ParallelFor(tbx, tby,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = 0.5*dxInv[2];
                    if (z_nd_arr) {
                        Real z_xf_lo = 0.25 * ( z_nd_arr(i,j,k  ) + z_nd_arr(i,j+1,k  )
                                              + z_nd_arr(i,j,k-1) + z_nd_arr(i,j+1,k-1) );
                        Real z_xf_hi = 0.25 * ( z_nd_arr(i,j,k+1) + z_nd_arr(i,j+1,k+1)
                                              + z_nd_arr(i,j,k+2) + z_nd_arr(i,j+1,k+2) );
                        dzInv = 1.0 / (z_xf_hi - z_xf_lo);
                    }
                    Real U_hi = dptr_u_plane(k+1) / dptr_r_plane(k+1);
                    Real U_lo = dptr_u_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_xf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    xmom_src_arr(i, j, k) -= wbar_xf * (U_hi - U_lo) * dzInv;
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real dzInv = 0.5*dxInv[2];
                    if (z_nd_arr) {
                        Real z_yf_lo = 0.25 * ( z_nd_arr(i,j,k  ) + z_nd_arr(i+1,j,k  )
                                              + z_nd_arr(i,j,k-1) + z_nd_arr(i+1,j,k-1) );
                        Real z_yf_hi = 0.25 * ( z_nd_arr(i,j,k+1) + z_nd_arr(i+1,j,k+1)
                                              + z_nd_arr(i,j,k+2) + z_nd_arr(i+1,j,k+2) );
                        dzInv = 1.0 / (z_yf_hi - z_yf_lo);
                    }
                    Real V_hi = dptr_v_plane(k+1) / dptr_r_plane(k+1);
                    Real V_lo = dptr_v_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_yf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    ymom_src_arr(i, j, k) -= wbar_yf * (V_hi - V_lo) * dzInv;
                });
            }
        }

        // *************************************************************************************
        // 5. Add nudging towards value specified in input sounding
        // *************************************************************************************
        if (solverChoice.nudging_from_input_sounding && is_slow_step)
        {
            int itime_n    = 0;
            int itime_np1  = 0;
            Real coeff_n   = Real(1.0);
            Real coeff_np1 = Real(0.0);

            Real tau_inv = Real(1.0) / input_sounding_data.tau_nudging;

            int n_sounding_times = input_sounding_data.input_sounding_time.size();

            for (int nt = 1; nt < n_sounding_times; nt++) {
                if (time > input_sounding_data.input_sounding_time[nt]) itime_n = nt;
            }
            if (itime_n == n_sounding_times-1) {
                itime_np1 = itime_n;
            } else {
                itime_np1 = itime_n+1;
                coeff_np1 = (time                                               - input_sounding_data.input_sounding_time[itime_n]) /
                            (input_sounding_data.input_sounding_time[itime_np1] - input_sounding_data.input_sounding_time[itime_n]);
                coeff_n   = Real(1.0) - coeff_np1;
            }

            int nr = Rho_comp;

            const Real* u_inp_sound_n   = input_sounding_data.U_inp_sound_d[itime_n].dataPtr();
            const Real* u_inp_sound_np1 = input_sounding_data.U_inp_sound_d[itime_np1].dataPtr();
            const Real* v_inp_sound_n   = input_sounding_data.V_inp_sound_d[itime_n].dataPtr();
            const Real* v_inp_sound_np1 = input_sounding_data.V_inp_sound_d[itime_np1].dataPtr();
            ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real nudge_u = (coeff_n*u_inp_sound_n[k] + coeff_np1*u_inp_sound_np1[k]) - (dptr_u_plane(k)/dptr_r_plane(k));
                nudge_u *= tau_inv;
                xmom_src_arr(i, j, k) += cell_data(i, j, k, nr) * nudge_u;
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real nudge_v = (coeff_n*v_inp_sound_n[k] + coeff_np1*v_inp_sound_np1[k]) - (dptr_v_plane(k)/dptr_r_plane(k));
                nudge_v *= tau_inv;
                ymom_src_arr(i, j, k) += cell_data(i, j, k, nr) * nudge_v;
            });
        }

        // *****************************************************************************
        // 6. Add NUMERICAL DIFFUSION terms
        // *****************************************************************************
#if 0
        if (l_use_ndiff) {
            const Array4<const Real>& mf_ux   = mapfac[MapFac::ux]->const_array(mfi);
            const Array4<const Real>& mf_uy   = mapfac[MapFac::uy]->const_array(mfi);
            const Array4<const Real>& mf_vx   = mapfac[MapFac::vx]->const_array(mfi);
            const Array4<const Real>& mf_vy   = mapfac[MapFac::vy]->const_array(mfi);
            NumericalDiffusion_Xmom(tbx, dt, solverChoice.num_diff_coeff,
                                    u, cell_data, xmom_src_arr, mf_ux, mf_uy);
            NumericalDiffusion_Ymom(tby, dt, solverChoice.num_diff_coeff,
                                    v, cell_data, ymom_src_arr, mf_vx, mf_vy);
        }
#endif

        // *****************************************************************************
        // 7. Add SPONGING
        // *****************************************************************************
        if (is_slow_step) {
            if (solverChoice.spongeChoice.sponge_type == "input_sponge")
            {
                ApplySpongeZoneBCsForMom_ReadFromFile(solverChoice.spongeChoice, geom, tbx, tby, cell_data,
                                                    z_cc_arr, xmom_src_arr, ymom_src_arr,
                                                    rho_u, rho_v, d_sponge_ptrs_at_lev);
            }
            else
            {
                ApplySpongeZoneBCsForMom(solverChoice.spongeChoice, geom, tbx, tby, tbz,
                                        xmom_src_arr, ymom_src_arr, zmom_src_arr, rho_u, rho_v, rho_w,
                                        r0, z_nd_arr, z_cc_arr);
            }

            if(solverChoice.init_type == InitType::HindCast and solverChoice.hindcast_lateral_forcing){

                const Array4<const Real>& rho_u_forecast_state  = (*forecast_state_at_lev)[IntVars::xmom].array(mfi);
                const Array4<const Real>& rho_v_forecast_state  = (*forecast_state_at_lev)[IntVars::ymom].array(mfi);
                const Array4<const Real>& rho_w_forecast_state  = (*forecast_state_at_lev)[IntVars::zmom].array(mfi);
                const Array4<const Real>& cons_forecast_state   = (*forecast_state_at_lev)[IntVars::cons].array(mfi);
                ApplyBndryForcing_Forecast(solverChoice, geom, tbx, tby, tbz, z_nd_arr,
                                           xmom_src_arr, ymom_src_arr, zmom_src_arr,
                                           rho_u, rho_v, rho_w,
                                           rho_u_forecast_state, rho_v_forecast_state, rho_w_forecast_state,
                                           cons_forecast_state);
            }
        }

        // *****************************************************************************
        // 8. Add CANOPY source terms
        // *****************************************************************************
        if (solverChoice.do_forest_drag &&
           ((is_slow_step && !use_canopy_fast) || (!is_slow_step && use_canopy_fast))) {
            ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = u(i, j, k);
                const Real uy = 0.25 * ( v(i, j  , k  ) + v(i-1, j  , k  )
                                       + v(i, j+1, k  ) + v(i-1, j+1, k  ) );
                const Real uz = 0.25 * ( w(i, j  , k  ) + w(i-1, j  , k  )
                                       + w(i, j  , k+1) + w(i-1, j  , k+1) );
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
                const Real f_drag = 0.5 * (f_drag_arr(i, j, k) + f_drag_arr(i-1, j, k));
                xmom_src_arr(i, j, k) -= f_drag * ux * windspeed;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = 0.25 * ( u(i  , j  , k  ) + u(i  , j-1, k  )
                                       + u(i+1, j  , k  ) + u(i+1, j-1, k  ) );
                const Real uy = v(i, j, k);
                const Real uz = 0.25 * ( w(i  , j  , k  ) + w(i  , j-1, k  )
                                       + w(i  , j  , k+1) + w(i  , j-1, k+1) );
                const amrex::Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
                const Real f_drag = 0.5 * (f_drag_arr(i, j, k) + f_drag_arr(i, j-1, k));
                ymom_src_arr(i, j, k) -= f_drag * uy * windspeed;
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const amrex::Real ux = 0.25 * ( u(i  , j  , k  ) + u(i+1, j  , k  )
                                              + u(i  , j  , k-1) + u(i+1, j  , k-1) );
                const amrex::Real uy = 0.25 * ( v(i  , j  , k  ) + v(i  , j+1, k  )
                                              + v(i  , j  , k-1) + v(i  , j+1, k-1) );
                const amrex::Real uz = w(i, j, k);
                const amrex::Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
                const Real f_drag = 0.5 * (f_drag_arr(i, j, k) + f_drag_arr(i, j, k-1));
                zmom_src_arr(i, j, k) -= f_drag * uz * windspeed;
            });
        }
        // *****************************************************************************
        // 9a. Add immersed source terms for terrain
        // *****************************************************************************
        if (solverChoice.terrain_type == TerrainType::ImmersedForcing &&
           ((is_slow_step && !use_ImmersedForcing_fast) || (!is_slow_step && use_ImmersedForcing_fast))) {
            // geometric properties
            const Real* dx_arr = geom.CellSize();
            const Real dx_x = dx_arr[0];
            const Real dx_y = dx_arr[1];
            const Real dx_z = dx_arr[2];

            const Real alpha_m = solverChoice.if_Cd_momentum;
            const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, 1./3.);
            const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
            const Real U_s = 1.0; // unit velocity scale

            // MOST parameters
            similarity_funs sfuns;
            const Real ggg        = CONST_GRAV;
            const Real kappa      = KAPPA;
            const Real z0                 = solverChoice.if_z0;
            const Real tflux_in           = solverChoice.if_surf_temp_flux;
            const Real Olen_in            = solverChoice.if_Olen_in;
            const bool l_use_most         = solverChoice.if_use_most;

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = u(i, j, k);
                const Real uy = 0.25 * ( v(i, j  , k  ) + v(i-1, j  , k  )
                                       + v(i, j+1, k  ) + v(i-1, j+1, k  ) );
                const Real uz = 0.25 * ( w(i, j  , k  ) + w(i-1, j  , k  )
                                       + w(i, j  , k+1) + w(i-1, j  , k+1) );
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
                const Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i-1, j, k));
                const Real t_blank_above = 0.5 * (t_blank_arr(i, j, k+1) + t_blank_arr(i-1, j, k+1));
                const Real CdM = std::min(drag_coefficient / (windspeed + tiny), 1000.0);
                const Real rho_xface = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );

                if ((t_blank > 0 && (t_blank_above == 0.0)) && l_use_most) { // force to MOST value
                    // calculate tangential velocity one cell above
                    const Real ux2r = u(i, j, k+1) ;
                    const Real uy2r = 0.25 * ( v(i, j  , k+1) + v(i-1, j  , k+1)
                                       + v(i, j+1, k+1) + v(i-1, j+1, k+1) ) ;
                    const Real h_windspeed2r = std::sqrt(ux2r * ux2r + uy2r * uy2r);

                    // MOST
                    const Real theta_xface = (0.5 * (cell_data(i,j,k  ,RhoTheta_comp) + cell_data(i-1,j,k, RhoTheta_comp))) / rho_xface;
                    const Real rho_xface_below    = 0.5 * ( cell_data(i,j,k-1,Rho_comp) + cell_data(i-1,j,k-1,Rho_comp) );
                    const Real theta_xface_below  = (0.5 * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i-1,j,k-1, RhoTheta_comp))) / rho_xface_below;
                    const Real theta_surf         = theta_xface_below;

                    Real psi_m = 0.0;
                    Real psi_h = 0.0;
                    Real ustar = h_windspeed2r * kappa / (std::log(1.5 * dx_z / z0) - psi_m); // calculated from bottom of cell. Maintains flexibility for different Vf values
                    Real tflux = (tflux_in != 1e-8) ? tflux_in : -(theta_xface - theta_surf) * ustar * kappa / (std::log(1.5 * dx_z / z0) - psi_h);
                    Real Olen = (Olen_in != 1e-8)   ? Olen_in  : -ustar * ustar * ustar * theta_xface / (kappa * ggg * tflux + tiny);
                    Real zeta  = 1.5 * dx_z / Olen;

                    // similarity functions
                    psi_m          = sfuns.calc_psi_m(zeta);
                    psi_h          = sfuns.calc_psi_h(zeta);
                    ustar = h_windspeed2r * kappa / (std::log(1.5 * dx_z / z0) - psi_m);

                    // prevent some unphysical math
                    if (!(ustar > 0.0 && !std::isnan(ustar))) { ustar = 0.0; }
                    if (!(ustar < 2.0 && !std::isnan(ustar))) { ustar = 2.0; }
                    if (psi_m > std::log(0.5 * dx_z / z0)) { psi_m = std::log(0.5 * dx_z / z0); }

                    // determine target velocity
                    const Real uTarget  = ustar / kappa * (std::log(0.5 * dx_z / z0) - psi_m);
                    Real uxTarget = uTarget * ux2r / (tiny + h_windspeed2r);
                    const Real bc_forcing_x = -(uxTarget - ux); // BC forcing pushes nonrelative velocity toward target velocity
                    xmom_src_arr(i, j, k) -= (1-t_blank) * rho_xface * CdM * U_s * bc_forcing_x; // if Vf low, force more strongly to MOST. If high, less forcing.
                } else {
                    xmom_src_arr(i, j, k) -= t_blank * rho_xface * CdM * ux * windspeed;
                }
            });
            ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = 0.25 * ( u(i  , j  , k  ) + u(i  , j-1, k  )
                                       + u(i+1, j  , k  ) + u(i+1, j-1, k  ) );
                const Real uy = v(i, j, k);
                const Real uz = 0.25 * ( w(i  , j  , k  ) + w(i  , j-1, k  )
                                       + w(i  , j  , k+1) + w(i  , j-1, k+1) );
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
                const Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i, j-1, k));
                const Real t_blank_above = 0.5 * (t_blank_arr(i, j, k+1) + t_blank_arr(i, j-1, k+1));
                const Real CdM = std::min(drag_coefficient / (windspeed + tiny), 1000.0);
                const Real rho_yface =  0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );

                if ((t_blank > 0 && (t_blank_above == 0.0)) && l_use_most) { // force to MOST value
                    // calculate tangential velocity one cell above
                    const Real ux2r = 0.25 * ( u(i  , j  , k+1) + u(i  , j-1, k+1)
                                       + u(i+1, j  , k+1) + u(i+1, j-1, k+1) );
                    const Real uy2r = v(i, j, k+1) ;
                    const Real h_windspeed2r = std::sqrt(ux2r * ux2r + uy2r * uy2r);

                    // MOST
                    const Real theta_yface = (0.5 * (cell_data(i,j,k  ,RhoTheta_comp) + cell_data(i,j-1,k, RhoTheta_comp))) / rho_yface;
                    const Real rho_yface_below    =  0.5 * ( cell_data(i,j,k-1,Rho_comp) + cell_data(i,j-1,k-1,Rho_comp) );
                    const Real theta_yface_below  = (0.5 * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i,j-1,k-1, RhoTheta_comp))) / rho_yface_below;
                    const Real theta_surf         = theta_yface_below;

                    Real psi_m = 0.0;
                    Real psi_h = 0.0;
                    Real ustar = h_windspeed2r * kappa / (std::log(1.5 * dx_z / z0) - psi_m); // calculated from bottom of cell. Maintains flexibility for different Vf values
                    Real tflux = (tflux_in != 1e-8) ? tflux_in : -(theta_yface - theta_surf) * ustar * kappa / (std::log(1.5 * dx_z / z0) - psi_h);
                    Real Olen = (Olen_in != 1e-8)   ? Olen_in  : -ustar * ustar * ustar * theta_yface / (kappa * ggg * tflux + tiny);
                    Real zeta  = 1.5 * dx_z / Olen;

                    // similarity functions
                    psi_m          = sfuns.calc_psi_m(zeta);
                    psi_h          = sfuns.calc_psi_h(zeta);
                    ustar = h_windspeed2r * kappa / (std::log(1.5 * dx_z / z0) - psi_m);

                    // prevent some unphysical math
                    if (!(ustar > 0.0 && !std::isnan(ustar))) { ustar = 0.0; }
                    if (!(ustar < 2.0 && !std::isnan(ustar))) { ustar = 2.0; }
                    if (psi_m > std::log(0.5 * dx_z / z0)) { psi_m = std::log(0.5 * dx_z / z0); }

                    // determine target velocity
                    const Real uTarget  = ustar / kappa * (std::log(0.5 * dx_z / z0) - psi_m);
                    Real uyTarget = uTarget * uy2r / (tiny + h_windspeed2r);
                    const Real bc_forcing_y = -(uyTarget - uy);  // BC forcing pushes nonrelative velocity toward target velocity
                    ymom_src_arr(i, j, k) -= (1 - t_blank) * rho_yface * CdM * U_s * bc_forcing_y; // if Vf low, force more strongly to MOST. If high, less forcing.
                } else {
                    ymom_src_arr(i, j, k) -= t_blank * rho_yface * CdM * uy * windspeed;
                }
            });
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = 0.25 * ( u(i  , j  , k  ) + u(i+1, j  , k  )
                                       + u(i  , j  , k-1) + u(i+1, j  , k-1) );
                const Real uy = 0.25 * ( v(i  , j  , k  ) + v(i  , j+1, k  )
                                       + v(i  , j  , k-1) + v(i  , j+1, k-1) );
                const Real uz = w(i, j, k);
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
                const Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i, j, k-1));
                const Real CdM = std::min(drag_coefficient / (windspeed + tiny), 1000.0);
                const Real rho_zface =  0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                zmom_src_arr(i, j, k) -= t_blank * rho_zface * CdM * uz * windspeed;
            });
        }

        // *****************************************************************************
        // 9b. Add immersed source terms for buildings
        // *****************************************************************************
        if ((solverChoice.buildings_type == BuildingsType::ImmersedForcing ) &&
           ((is_slow_step && !use_ImmersedForcing_fast) || (!is_slow_step && use_ImmersedForcing_fast)))
        {
            // geometric properties
            const Real* dx_arr = geom.CellSize();
            const Real dx_x = dx_arr[0];
            const Real dx_y = dx_arr[1];

            const Real alpha_m          = solverChoice.if_Cd_momentum;
            const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
            const Real min_t_blank      = 0.005; // threshold for where immersed forcing acts

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = u(i, j, k);
                const Real uy = 0.25 * ( v(i, j  , k  ) + v(i-1, j  , k  )
                                       + v(i, j+1, k  ) + v(i-1, j+1, k  ) );
                const Real uz = 0.25 * ( w(i, j  , k  ) + w(i-1, j  , k  )
                                       + w(i, j  , k+1) + w(i-1, j  , k+1) );
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);

                Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i-1, j, k));
                if (t_blank < min_t_blank) { t_blank = 0.0; }
                const Real dx_z    = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
                const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, 1./3.);
                const Real CdM = std::min(drag_coefficient / (windspeed + tiny), 1000.0);
                const Real rho_xface = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                xmom_src_arr(i, j, k) -= t_blank * rho_xface * CdM * ux * windspeed;
            });
            ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = 0.25 * ( u(i  , j  , k  ) + u(i  , j-1, k  )
                                       + u(i+1, j  , k  ) + u(i+1, j-1, k  ) );
                const Real uy = v(i, j, k);
                const Real uz = 0.25 * ( w(i  , j  , k  ) + w(i  , j-1, k  )
                                       + w(i  , j  , k+1) + w(i  , j-1, k+1) );
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);

                Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i, j-1, k));
                if (t_blank < min_t_blank) { t_blank = 0.0; }
                const Real dx_z    = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
                const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, 1./3.);
                const Real CdM = std::min(drag_coefficient / (windspeed + tiny), 1000.0);
                const Real rho_yface =  0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                ymom_src_arr(i, j, k) -= t_blank * rho_yface * CdM * uy * windspeed;
            });
            ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                const Real ux = 0.25 * ( u(i  , j  , k  ) + u(i+1, j  , k  )
                                       + u(i  , j  , k-1) + u(i+1, j  , k-1) );
                const Real uy = 0.25 * ( v(i  , j  , k  ) + v(i  , j+1, k  )
                                       + v(i  , j  , k-1) + v(i  , j+1, k-1) );
                const Real uz = w(i, j, k);
                const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);

                Real t_blank = 0.5 * (t_blank_arr(i, j, k) + t_blank_arr(i, j, k-1));
                if (t_blank < min_t_blank) { t_blank = 0.0; }
                const Real dx_z    = (z_nd_arr) ? (z_nd_arr(i,j,k) - z_nd_arr(i,j,k-1)) : dx_arr[2]; // ASW double check
                const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, 1./3.);
                const Real CdM = std::min(drag_coefficient / (windspeed + tiny), 1000.0);
                const Real rho_zface =  0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                zmom_src_arr(i, j, k) -= t_blank * rho_zface * CdM * uz * windspeed;
            });
        }

        // *****************************************************************************
        // 10. Enforce constant mass flux
        // *****************************************************************************
        if (is_slow_step && (enforce_massflux_x || enforce_massflux_y)) {
            Real tau_inv = Real(1.0) / solverChoice.const_massflux_tau;

            ParallelFor(tbx, tby,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                xmom_src_arr(i, j, k) += tau_inv * (rhoUA_target - rhoUA);
            },
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                ymom_src_arr(i, j, k) += tau_inv * (rhoVA_target - rhoVA);
            });
        }

    } // mfi
}
