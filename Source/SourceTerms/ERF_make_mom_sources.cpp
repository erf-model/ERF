#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BCRec.H>
#include <AMReX_TableData.H>
#include <AMReX_GpuContainers.H>

#include <NumericalDiffusion.H>
#include <PlaneAverage.H>
#include <TI_slow_headers.H>
#include <Src_headers.H>
#include <Utils.H>

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
 * @param[in] xflux_imask_lev thin-body mask on x-faces
 * @param[in] yflux_imask_lev thin-body mask on y-faces
 * @param[in] zflux_imask_lev thin-body mask on z-faces
 * @param[in] thin_xforce_lev x-component of forces on thin immersed bodies
 * @param[in] thin_yforce_lev y-component of forces on thin immersed bodies
 * @param[in] thin_zforce_lev z-component of forces on thin immersed bodies
 */

void make_mom_sources (int /*level*/,
                       int /*nrk*/, Real dt,
                       Vector<MultiFab>& S_data,
                       const  MultiFab & S_prim,
                       const  MultiFab & xvel,
                       const  MultiFab & yvel,
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
                       int n_qstate,
                       std::unique_ptr<iMultiFab>& xflux_imask_lev,
                       std::unique_ptr<iMultiFab>& yflux_imask_lev,
                       std::unique_ptr<iMultiFab>& zflux_imask_lev,
                       std::unique_ptr<MultiFab>& thin_xforce_lev,
                       std::unique_ptr<MultiFab>& thin_yforce_lev,
                       std::unique_ptr<MultiFab>& thin_zforce_lev)
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
    //    6. numerical diffusion for (xmom,ymom,zmom)
    //    7. sponge
    //    8. thin body immersed boundary forcing
    // *****************************************************************************
    const bool l_use_ndiff      = solverChoice.use_NumDiff;

    const bool l_have_thin_xforce = (thin_xforce_lev != nullptr);
    const bool l_have_thin_yforce = (thin_yforce_lev != nullptr);
    const bool l_have_thin_zforce = (thin_zforce_lev != nullptr);

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
    auto geostrophic_wind = solverChoice.custom_geostrophic_profile;

    // *****************************************************************************
    // Data for Rayleigh damping
    // *****************************************************************************
    auto rayleigh_damp_U  = solverChoice.rayleigh_damp_U;
    auto rayleigh_damp_V  = solverChoice.rayleigh_damp_V;
    auto rayleigh_damp_W  = solverChoice.rayleigh_damp_W;

    Real*      tau = d_rayleigh_ptrs_at_lev[Rayleigh::tau];
    Real*     ubar = d_rayleigh_ptrs_at_lev[Rayleigh::ubar];
    Real*     vbar = d_rayleigh_ptrs_at_lev[Rayleigh::vbar];
    Real*     wbar = d_rayleigh_ptrs_at_lev[Rayleigh::wbar];

    // *****************************************************************************
    // Planar averages for subsidence terms
    // *****************************************************************************
    Table1D<Real> dptr_u_plane, dptr_v_plane;
    TableData<Real, 1> u_plane_tab, v_plane_tab;

    if (dptr_wbar_sub) {
        PlaneAverage u_ave(&xvel  , geom, solverChoice.ave_plane, true);
        PlaneAverage v_ave(&yvel  , geom, solverChoice.ave_plane, true);

        u_ave.compute_averages(ZDir(), u_ave.field());
        v_ave.compute_averages(ZDir(), v_ave.field());

        int u_ncell = u_ave.ncell_line();
        int v_ncell = v_ave.ncell_line();
        Gpu::HostVector<    Real> u_plane_h(u_ncell), v_plane_h(v_ncell);
        Gpu::DeviceVector<  Real> u_plane_d(u_ncell), v_plane_d(v_ncell);

        u_ave.line_average(0             , u_plane_h);
        v_ave.line_average(0             , v_plane_h);

        Gpu::copyAsync(Gpu::hostToDevice, u_plane_h.begin(), u_plane_h.end(), u_plane_d.begin());
        Gpu::copyAsync(Gpu::hostToDevice, v_plane_h.begin(), v_plane_h.end(), v_plane_d.begin());

        Real* dptr_u = u_plane_d.data();
        Real* dptr_v = v_plane_d.data();

        IntVect ng_u = xvel.nGrowVect();
        IntVect ng_v = yvel.nGrowVect();
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
    // Create the BUOYANCY forcing term in the z-direction
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

        // *****************************************************************************
        // Add CORIOLIS forcing (this assumes east is +x, north is +y)
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
        // Add RAYLEIGH damping
        // *****************************************************************************
        if (rayleigh_damp_U) {
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                Real uu = rho_u(i,j,k) / rho_on_u_face;
                xmom_src_arr(i, j, k) -= tau[k] * (uu - ubar[k]) * rho_on_u_face;
            });
        }

        if (rayleigh_damp_V) {
            ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                Real vv = rho_v(i,j,k) / rho_on_v_face;
                ymom_src_arr(i, j, k) -= tau[k] * (vv - vbar[k]) * rho_on_v_face;
            });
        }

        if (rayleigh_damp_W) {
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_w_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
                Real ww = rho_w(i,j,k) / rho_on_w_face;
                zmom_src_arr(i, j, k) -= tau[k] * (ww - wbar[k]) * rho_on_w_face;
            });
        }

        // *****************************************************************************
        // Add constant GEOSTROPHIC forcing
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
        // Add height-dependent GEOSTROPHIC forcing
        // *****************************************************************************
        if (geostrophic_wind) {
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_u_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
                xmom_src_arr(i, j, k) += rho_on_u_face * dptr_u_geos[k];
            });
            ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real rho_on_v_face = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
                ymom_src_arr(i, j, k) += rho_on_v_face * dptr_v_geos[k];
            });
        } // geostrophic_wind

        // *****************************************************************************
        // Add custom SUBSIDENCE terms
        // *****************************************************************************
        if (solverChoice.custom_w_subsidence) {
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                xmom_src_arr(i, j, k) += dptr_wbar_sub[k] *
                    0.5 * (dptr_u_plane(k+1) - dptr_u_plane(k-1)) * dxInv[2];
            });
            ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                ymom_src_arr(i, j, k) += dptr_wbar_sub[k] *
                    0.5 * (dptr_v_plane(k+1) - dptr_v_plane(k-1)) * dxInv[2];
            });
        }

        // *****************************************************************************
        // Add NUMERICAL DIFFUSION terms
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
        // Add SPONGING
        // *****************************************************************************
        ApplySpongeZoneBCsForMom(solverChoice.spongeChoice, geom, tbx, tby, tbz,
                                 xmom_src_arr, ymom_src_arr, zmom_src_arr, rho_u, rho_v, rho_w);

    } // mfi

    // *****************************************************************************
    // If a thin immersed body is present, add forcing terms
    // *****************************************************************************
    if (l_have_thin_xforce) {
        MultiFab::Copy(*thin_xforce_lev, xmom_src, 0, 0, 1, 0);
        thin_xforce_lev->mult(-1., 0, 1, 0);
        ApplyInvertedMask(*thin_xforce_lev, *xflux_imask_lev, 0);
        MultiFab::Add(xmom_src, *thin_xforce_lev, 0, 0, 1, 0);
    }

    if (l_have_thin_yforce) {
        MultiFab::Copy(*thin_yforce_lev, ymom_src, 0, 0, 1, 0);
        thin_yforce_lev->mult(-1., 0, 1, 0);
        ApplyInvertedMask(*thin_yforce_lev, *yflux_imask_lev, 0);
        MultiFab::Add(ymom_src, *thin_yforce_lev, 0, 0, 1, 0);
    }

    if (l_have_thin_zforce) {
        MultiFab::Copy(*thin_zforce_lev, zmom_src, 0, 0, 1, 0);
        thin_zforce_lev->mult(-1., 0, 1, 0);
        ApplyInvertedMask(*thin_zforce_lev, *zflux_imask_lev, 0);
        MultiFab::Add(zmom_src, *thin_zforce_lev, 0, 0, 1, 0);
    }

#if 0
#ifndef AMREX_USE_GPU
    if (l_have_thin_xforce) {
        // TODO: Implement particles to better track and output these data
        if (nrk==2) {
            for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.nodaltilebox(0);
                const Array4<const Real> & fx = thin_xforce_lev->const_array(mfi);
                const Array4<const int> & mask = xflux_imask_lev->const_array(mfi);
                ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (mask(i,j,k)==0) {
                        amrex::AllPrint() << "thin body fx"<<IntVect(i,j,k)<<" = " << fx(i,j,k) << std::endl;
                    }
                });
            }
        }
    }
#endif
#endif

#if 0
#ifndef AMREX_USE_GPU
    if (l_have_thin_yforce) {
        // TODO: Implement particles to better track and output these data
        if (nrk==2) {
            for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& tby = mfi.nodaltilebox(1);
                const Array4<const Real> & fy = thin_yforce_lev->const_array(mfi);
                const Array4<const int> & mask = yflux_imask_lev->const_array(mfi);
                ParallelFor(tby, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (mask(i,j,k)==0) {
                        amrex::AllPrint() << "thin body fy"<<IntVect(i,j,k)<<" = " << fy(i,j,k) << std::endl;
                    }
                });
            }
        }
    }
#endif
#endif

#if 0
#ifndef AMREX_USE_GPU
    if (l_have_thin_zforce) {
        // TODO: Implement particles to better track and output these data
        if (nrk==2) {
            for ( MFIter mfi(S_data[IntVars::cons],TileNoZ()); mfi.isValid(); ++mfi)
            {
                const Box& tbz = mfi.nodaltilebox(2);
                const Array4<const Real> & fz = thin_zforce_lev->const_array(mfi);
                const Array4<const int> & mask = zflux_imask_lev->const_array(mfi);
                ParallelFor(tbz, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (mask(i,j,k)==0) {
                        amrex::AllPrint() << "thin body fz"<<IntVect(i,j,k)<<" = " << fz(i,j,k) << std::endl;
                    }
                });
            }
        }
    }
#endif
#endif
}
