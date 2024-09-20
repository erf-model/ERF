#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_TableData.H>
#include <AMReX_GpuContainers.H>

#include <ERF_NumericalDiffusion.H>
#include <ERF_PlaneAverage.H>
#include <ERF_TI_slow_headers.H>
#include <ERF_Src_headers.H>
#include <ERF_Utils.H>

using namespace amrex;

/**
 * Function for computing the slow RHS for the evolution equations for the density, potential temperature and momentum.
 *
 * @param[in]  level level of resolution
 * @param[in]  nrk   which RK stage
 * @param[in]  dt    slow time step
 * @param[in]  S_data current solution
 * @param[in]  S_prim primitive variables (i.e. conserved variables divided by density)
 * @param[in]  xvel x-component of velocity
 * @param[in]  yvel y-component of velocity
 * @param[in] xmom_src source terms for x-momentum
 * @param[in] ymom_src source terms for y-momentum
 * @param[in] zmom_src source terms for z-momentum
 * @param[in]  geom   Container for geometric information
 * @param[in]  solverChoice  Container for solver parameters
 * @param[in] mapfac_m map factor at cell centers
 * @param[in] mapfac_u map factor at x-faces
 * @param[in] mapfac_v map factor at y-faces
 * @param[in] dptr_u_geos  custom geostrophic wind profile
 * @param[in] dptr_v_geos  custom geostrophic wind profile
 * @param[in] dptr_wbar_sub  subsidence source term
 * @param[in] d_rayleigh_ptrs_at_lev  Vector of {strength of Rayleigh damping, reference value for xvel/yvel/zvel/theta} used to define Rayleigh damping
 * @param[in] n_qstate number of moisture components
 */

void make_mom_sources (int /*level*/,
                       int /*nrk*/, Real dt, Real time,
                       Vector<MultiFab>& S_data,
                       const  MultiFab & S_prim,
                       std::unique_ptr<MultiFab>& z_phys_nd,
                       std::unique_ptr<MultiFab>& z_phys_cc,
                       const  MultiFab & /*xvel*/,
                       const  MultiFab & /*yvel*/,
                              MultiFab & xmom_src,
                              MultiFab & ymom_src,
                              MultiFab & zmom_src,
                       MultiFab* r0,
                       const Geometry geom,
                       const SolverChoice& solverChoice,
                       std::unique_ptr<MultiFab>& mapfac_m,
                       std::unique_ptr<MultiFab>& mapfac_u,
                       std::unique_ptr<MultiFab>& mapfac_v,
                       const Real* dptr_u_geos,
                       const Real* dptr_v_geos,
                       const Real* dptr_wbar_sub,
                       const Vector<Real*> d_rayleigh_ptrs_at_lev,
                       const Vector<Real*> d_sponge_ptrs_at_lev,
                       InputSoundingData& input_sounding_data,
                       int n_qstate)
{
    BL_PROFILE_REGION("erf_make_mom_sources()");

    Box domain(geom.Domain());
    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();

    // Initialize sources to zero since we re-compute them ever RK stage
    xmom_src.setVal(0.0);
    ymom_src.setVal(0.0);
    zmom_src.setVal(0.0);

    // *****************************************************************************
    // Define source term for all three components of momenta from
    //    1. buoyancy           for (zmom)
    //    2. Coriolis forcing   for (xmom,ymom,zmom)
    //    3. Rayleigh damping   for (xmom,ymom,zmom)
    //    4. Constant / height-dependent geostrophic forcing
    //    5. subsidence
    //    6. nudging towards input sounding data
    //    7. numerical diffusion for (xmom,ymom,zmom)
    //    8. sponge
    // *****************************************************************************
    const bool l_use_ndiff    = solverChoice.use_NumDiff;
    const bool use_terrain    = solverChoice.use_terrain;

    // *****************************************************************************
    // Data for Coriolis forcing
    // *****************************************************************************
    auto use_coriolis         = solverChoice.use_coriolis;
    auto coriolis_factor      = solverChoice.coriolis_factor;
    auto cosphi               = solverChoice.cosphi;
    auto sinphi               = solverChoice.sinphi;

    // *****************************************************************************
    // Flag for Geostrophic forcing
    // *****************************************************************************
    auto abl_geo_forcing  = solverChoice.abl_geo_forcing;
    auto geo_wind_profile = solverChoice.have_geo_wind_profile;

    // *****************************************************************************
    // Data for Rayleigh damping
    // *****************************************************************************
    auto rayleigh_damp_U  = solverChoice.rayleigh_damp_U;
    auto rayleigh_damp_V  = solverChoice.rayleigh_damp_V;
    auto rayleigh_damp_W  = solverChoice.rayleigh_damp_W;

    Real*     ubar = d_rayleigh_ptrs_at_lev[Rayleigh::ubar];
    Real*     vbar = d_rayleigh_ptrs_at_lev[Rayleigh::vbar];
    Real*     wbar = d_rayleigh_ptrs_at_lev[Rayleigh::wbar];

    // *****************************************************************************
    // Planar averages for subsidence terms
    // *****************************************************************************
    Table1D<Real>     dptr_r_plane, dptr_u_plane, dptr_v_plane;
    TableData<Real, 1> r_plane_tab,  u_plane_tab,  v_plane_tab;

    if (dptr_wbar_sub || solverChoice.nudging_from_input_sounding)
    {
        // Rho
        PlaneAverage r_ave(&(S_data[IntVars::cons]), geom, solverChoice.ave_plane, true);
        r_ave.compute_averages(ZDir(), r_ave.field());

        int ncell = r_ave.ncell_line();
        Gpu::HostVector<    Real> r_plane_h(ncell);
        Gpu::DeviceVector<  Real> r_plane_d(ncell);

        r_ave.line_average(Rho_comp, r_plane_h);

        Gpu::copyAsync(Gpu::hostToDevice, r_plane_h.begin(), r_plane_h.end(), r_plane_d.begin());

        Real* dptr_r = r_plane_d.data();

        IntVect ng_c = S_data[IntVars::cons].nGrowVect();
        Box tdomain  = domain; tdomain.grow(2,ng_c[2]);
        r_plane_tab.resize({tdomain.smallEnd(2)}, {tdomain.bigEnd(2)});

        int offset = ng_c[2];
        dptr_r_plane = r_plane_tab.table();
        ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_r_plane(k-offset) = dptr_r[k];
        });

        // U and V momentum
        PlaneAverage u_ave(&(S_data[IntVars::xmom]), geom, solverChoice.ave_plane, true);
        PlaneAverage v_ave(&(S_data[IntVars::ymom]), geom, solverChoice.ave_plane, true);

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

        IntVect ng_u = S_data[IntVars::xmom].nGrowVect();
        IntVect ng_v = S_data[IntVars::ymom].nGrowVect();
        Box udomain = domain; udomain.grow(2,ng_u[2]);
        Box vdomain = domain; vdomain.grow(2,ng_v[2]);
        u_plane_tab.resize({udomain.smallEnd(2)}, {udomain.bigEnd(2)});
        v_plane_tab.resize({vdomain.smallEnd(2)}, {vdomain.bigEnd(2)});

        int u_offset = ng_u[2];
        dptr_u_plane = u_plane_tab.table();
        ParallelFor(u_ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_u_plane(k-u_offset) = dptr_u[k];
        });

        int v_offset = ng_v[2];
        dptr_v_plane = v_plane_tab.table();
        ParallelFor(v_ncell, [=] AMREX_GPU_DEVICE (int k) noexcept
        {
            dptr_v_plane(k-v_offset) = dptr_v[k];
        });
    }

    // *****************************************************************************
    // 1. Create the BUOYANCY forcing term in the z-direction
    // *****************************************************************************
    make_buoyancy(S_data, S_prim, zmom_src, geom, solverChoice, r0, n_qstate);

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

        const Array4<      Real>& xmom_src_arr = xmom_src.array(mfi);
        const Array4<      Real>& ymom_src_arr = ymom_src.array(mfi);
        const Array4<      Real>& zmom_src_arr = zmom_src.array(mfi);

        // Map factors
        const Array4<const Real>& mf_m   = mapfac_m->const_array(mfi);
        const Array4<const Real>& mf_u   = mapfac_u->const_array(mfi);
        const Array4<const Real>& mf_v   = mapfac_v->const_array(mfi);

        const Array4<const Real>& z_nd_arr = (use_terrain) ? z_phys_nd->const_array(mfi) : Array4<Real>{};
        const Array4<const Real>& z_cc_arr = (use_terrain) ? z_phys_cc->const_array(mfi) : Array4<Real>{};

        // *****************************************************************************
        // 2. Add CORIOLIS forcing (this assumes east is +x, north is +y)
        // *****************************************************************************
        if (use_coriolis) {
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
        } // use_coriolis

        // *****************************************************************************
        // 3. Add RAYLEIGH damping
        // *****************************************************************************
        Real zlo      = geom.ProbLo(2);
        Real dz       = geom.CellSize(2);
        Real ztop     = solverChoice.rayleigh_ztop;
        Real zdamp    = solverChoice.rayleigh_zdamp;
        Real dampcoef = solverChoice.rayleigh_dampcoef;

        if (rayleigh_damp_U) {
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real zcc = (z_cc_arr) ? z_cc_arr(i,j,k) : zlo + (k+0.5)*dz;
                Real zfrac = 1 - (ztop - zcc) / zdamp;
                if (zfrac > 0) {
                    Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                    Real uu = rho_u(i,j,k) / rho_on_u_face;
                    Real sinefac = std::sin(PIoTwo*zfrac);
                    xmom_src_arr(i, j, k) -= dampcoef*sinefac*sinefac * (uu - ubar[k]) * rho_on_u_face;
                }
            });
        }

        if (rayleigh_damp_V) {
            ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real zcc = (z_cc_arr) ? z_cc_arr(i,j,k) : zlo + (k+0.5)*dz;
                Real zfrac = 1 - (ztop - zcc) / zdamp;
                if (zfrac > 0) {
                    Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                    Real vv = rho_v(i,j,k) / rho_on_v_face;
                    Real sinefac = std::sin(PIoTwo*zfrac);
                    ymom_src_arr(i, j, k) -= dampcoef*sinefac*sinefac * (vv - vbar[k]) * rho_on_v_face;
                }
            });
        }

        if (rayleigh_damp_W) {
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real zstag = (z_nd_arr) ? z_nd_arr(i,j,k) : zlo + k*dz;
                Real zfrac = 1 - (ztop - zstag) / zdamp;
                if (zfrac > 0) {
                    Real rho_on_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                    Real ww = rho_w(i,j,k) / rho_on_w_face;
                    Real sinefac = std::sin(PIoTwo*zfrac);
                    zmom_src_arr(i, j, k) -= dampcoef*sinefac*sinefac * (ww - wbar[k]) * rho_on_w_face;
                }
            });
        }

        // *****************************************************************************
        // 4. Add constant GEOSTROPHIC forcing
        // *****************************************************************************
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
            xmom_src_arr(i, j, k) += rho_on_u_face * abl_geo_forcing[0];
        });
        ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
            ymom_src_arr(i, j, k) += rho_on_v_face * abl_geo_forcing[1];
        });
        ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real rho_on_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
            zmom_src_arr(i, j, k) += rho_on_w_face * abl_geo_forcing[2];
        });

        // *****************************************************************************
        // 4. Add height-dependent GEOSTROPHIC forcing
        // *****************************************************************************
        if (geo_wind_profile) {
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                xmom_src_arr(i, j, k) -= coriolis_factor * rho_on_u_face * dptr_v_geos[k] * sinphi;
            });
            ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                ymom_src_arr(i, j, k) += coriolis_factor * rho_on_v_face * dptr_u_geos[k] * sinphi;
            });
        } // geo_wind_profile

        // *****************************************************************************
        // 5. Add custom SUBSIDENCE terms
        // *****************************************************************************
        if (solverChoice.custom_w_subsidence) {
            if (solverChoice.custom_forcing_prim_vars) {
                const int nr = Rho_comp;
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,nr) + cell_data(i-1,j,k,nr) );
                    Real U_hi = dptr_u_plane(k+1) / dptr_r_plane(k+1);
                    Real U_lo = dptr_u_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_xf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    xmom_src_arr(i, j, k) -= rho_on_u_face * wbar_xf *
                                             0.5 * (U_hi - U_lo) * dxInv[2];
                });
                ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,nr) + cell_data(i,j-1,k,nr) );
                    Real V_hi = dptr_v_plane(k+1) / dptr_r_plane(k+1);
                    Real V_lo = dptr_v_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_yf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    ymom_src_arr(i, j, k) -= rho_on_v_face * wbar_yf *
                                             0.5 * (V_hi - V_lo) * dxInv[2];
                });
            } else {
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real U_hi = dptr_u_plane(k+1) / dptr_r_plane(k+1);
                    Real U_lo = dptr_u_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_xf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    xmom_src_arr(i, j, k) -= wbar_xf *
                                             0.5 * (U_hi - U_lo) * dxInv[2];
                });
                ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real V_hi = dptr_v_plane(k+1) / dptr_r_plane(k+1);
                    Real V_lo = dptr_v_plane(k-1) / dptr_r_plane(k-1);
                    Real wbar_yf = 0.5 * (dptr_wbar_sub[k] + dptr_wbar_sub[k+1]);
                    ymom_src_arr(i, j, k) -= wbar_yf *
                                             0.5 * (V_hi - V_lo) * dxInv[2];
                });
            }
        }

        // *************************************************************************************
        // 6. Add nudging towards value specified in input sounding
        // *************************************************************************************
        if (solverChoice.nudging_from_input_sounding)
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
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real nudge_u = (coeff_n*u_inp_sound_n[k] + coeff_np1*u_inp_sound_np1[k]) - (dptr_u_plane(k)/dptr_r_plane(k));
                nudge_u *= tau_inv;
                xmom_src_arr(i, j, k) += cell_data(i, j, k, nr) * nudge_u;
            });

            const Real* v_inp_sound_n   = input_sounding_data.V_inp_sound_d[itime_n].dataPtr();
            const Real* v_inp_sound_np1 = input_sounding_data.V_inp_sound_d[itime_np1].dataPtr();
            ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real nudge_v = (coeff_n*v_inp_sound_n[k] + coeff_np1*v_inp_sound_np1[k]) - (dptr_v_plane(k)/dptr_r_plane(k));
                nudge_v *= tau_inv;
                ymom_src_arr(i, j, k) += cell_data(i, j, k, nr) * nudge_v;
            });
        }

        // *****************************************************************************
        // 7. Add NUMERICAL DIFFUSION terms
        // *****************************************************************************
        if (l_use_ndiff) {
            NumericalDiffusion(tbx, 0, 1, dt, solverChoice.NumDiffCoeff,
                               rho_u, xmom_src_arr, mf_m, mf_v, false, true);
            NumericalDiffusion(tby, 0, 1, dt, solverChoice.NumDiffCoeff,
                               rho_v, ymom_src_arr, mf_u, mf_m, true, false);
            NumericalDiffusion(tbz, 0, 1, dt, solverChoice.NumDiffCoeff,
                               rho_w, zmom_src_arr, mf_u, mf_v, false, false);
        }

        // *****************************************************************************
        // 8. Add SPONGING
        // *****************************************************************************
        if(solverChoice.spongeChoice.sponge_type == "input_sponge")
        {
             ApplySpongeZoneBCsForMom_ReadFromFile(solverChoice.spongeChoice, geom, tbx, tby, cell_data,
                                 xmom_src_arr, ymom_src_arr, rho_u, rho_v, d_sponge_ptrs_at_lev);
        }
        else
        {
            ApplySpongeZoneBCsForMom(solverChoice.spongeChoice, geom, tbx, tby, tbz,
                                 xmom_src_arr, ymom_src_arr, zmom_src_arr, rho_u, rho_v, rho_w);
        }

    } // mfi
}
