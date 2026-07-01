#include "ERF.H"
#include "ERF_Plotfile2DCatalog.H"
#include "ERF_NCPlotFile.H"
#include "ERF_Plotfile2DUtils.H"
#include "Diagnostics/ERF_SurfaceFluxDiagnostics.H"
#include "ERF_EpochTime.H"
#include "ERF_SrcHeaders.H"
#include "ERF_StormDiagnostics.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_Utils.H"

using namespace amrex;

namespace
{

// ERF's 2D plotfile path is intentionally an output assembly layer: it
// validates the requested names, assembles the slab geometry, fills existing
// diagnostics, and dispatches to AMReX or NetCDF output. Nontrivial science
// diagnostics should live in dedicated modules and be called from here.
// Keep the fill order below synchronized with plotfile2d::diagnostic_catalog()
// until the fill blocks move into dedicated diagnostic modules.

std::string make_2d_plotfile_name (int which,
                                   int plot_step,
                                   const std::string& plot2d_file_1,
                                   const std::string& plot2d_file_2,
                                   int file_name_digits)
{
    if (which == 1) {
        return Concatenate(plot2d_file_1, plot_step, file_name_digits);
    }

    if (which == 2) {
        return Concatenate(plot2d_file_2, plot_step, file_name_digits);
    }

    Abort(plotfile2d::format_invalid_2d_stream_error(which));
    return {};
}

Vector<Geometry> make_2d_plot_geometries (const Vector<Geometry>& geom,
                                          int finest_level)
{
    Vector<Geometry> my_geom(finest_level+1);

    Array<int,AMREX_SPACEDIM> is_per;
    is_per[0] = 0;
    is_per[1] = 0;
    is_per[2] = 0;
    if (geom[0].isPeriodic(0)) { is_per[0] = 1; }
    if (geom[0].isPeriodic(1)) { is_per[1] = 1; }

    int coord_sys = 0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        Box slab = makeSlab(geom[lev].Domain(), 2, 0);
        auto const slab_lo = lbound(slab);
        auto const slab_hi = ubound(slab);

        // The 2D slab uses the same horizontal geometry as the 3D domain but
        // collapses the vertical interval to the single slab that is written.
        Real dz = geom[lev].CellSize(2);
        RealBox rb = geom[lev].ProbDomain();
        rb.setLo(2,  slab_lo.z   * dz);
        rb.setHi(2, (slab_hi.z+1) * dz);
        my_geom[lev].define(slab, rb, coord_sys, is_per);
    }

    return my_geom;
}

void warn_for_unavailable_2d_plot_vars (const std::string& parameter_name,
                                        const Vector<std::string>& unavailable_names,
                                        const Vector<std::string>& available_names)
{
    if (unavailable_names.empty() || !ParallelDescriptor::IOProcessor()) {
        return;
    }

    for (const auto& name : unavailable_names) {
        Warning(plotfile2d::format_unavailable_2d_plot_var_warning(parameter_name, name, available_names));
    }
}

} // namespace
void
ERF::setPlotVariables2D (const std::string& pp_plot_var_names, Vector<std::string>& plot_var_names)
{
    ParmParse pp(pp_prefix);

    if (!pp.contains(pp_plot_var_names.c_str())) {
        //
        // The default is to add none of the variables to the list
        //
        plot_var_names.clear();
        return;
    }

    Vector<std::string> requested_plot_names;
    std::string nm;
    const int nPltVars = pp.countval(pp_plot_var_names.c_str());
    for (int i = 0; i < nPltVars; ++i) {
        pp.get(pp_plot_var_names.c_str(), nm, i);
        requested_plot_names.push_back(nm);
    }

    const auto available_names = plotfile2d::diagnostic_names();

    // Keep the canonical built-in 2D ordering so the plotfile component layout
    // stays stable even if the input request order changes.
    const auto selection = plotfile2d::select_requested_plot_variables(requested_plot_names,
                                                                        available_names);
    plot_var_names = selection.accepted;

    // Unknown 2D names are skipped rather than aborting because the 2D plot
    // list is intentionally user-configurable and may include names that are not
    // compiled into a given build. The warning is still explicit so the user
    // can correct the input deck.
    warn_for_unavailable_2d_plot_vars(plotfile2d::format_plot2d_parameter_name(pp_prefix, pp_plot_var_names),
                                      selection.unavailable, available_names);
}

void
ERF::Write2DPlotFile (int which, PlotFileType plotfile_type, Vector<std::string> plot_var_names)
{
    const int plot_step = istep[0];
    const std::string plotfilename = make_2d_plotfile_name(which, plot_step,
                                                           plot2d_file_1, plot2d_file_2,
                                                           file_name_digits);

    const Vector<std::string> varnames = PlotFileVarNames(plot_var_names);
    const int ncomp_mf = static_cast<int>(varnames.size());

    if (ncomp_mf == 0) return;

    // Vector of MultiFabs for cell-centered data
    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(ba2d[lev], dmap[lev], ncomp_mf, 0);
    }


    // **********************************************************************************************
    // (Effectively) 2D arrays
    // **********************************************************************************************
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Make sure getPgivenRTh and getTgivenRandRTh don't fail
        if (check_for_nans) {
            check_for_negative_theta(vars_new[lev][Vars::cons]);
        }

        int mf_comp = 0;

        // Set all components to zero in case they aren't defined below
        mf[lev].setVal(0.0);

        // Expose domain khi and klo at each level
        int klo = geom[lev].Domain().smallEnd(2);
        int khi = geom[lev].Domain().bigEnd(2);

        if (containerHasElement(plot_var_names, "z_surf")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<const Real>& z_phys_arr = z_phys_nd[lev]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   derdat(i, j, k, mf_comp) = Compute_Z_AtWFace(i, j, 0, z_phys_arr);
                });
            }
            mf_comp++;
        }

        if (containerHasElement(plot_var_names, "landmask")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<const int>& lmask_arr = lmask_lev[lev][0]->const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   derdat(i, j, k, mf_comp) = static_cast<Real>(lmask_arr(i, j, 0));
                });
            }
            mf_comp++;
        }

        if (containerHasElement(plot_var_names, "mapfac")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                const Array4<Real>& mf_m   = mapfac[lev][MapFacType::m_x]->array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                   derdat(i ,j ,k, mf_comp) = mf_m(i,j,0);
                });
            }
            mf_comp++;
        }

        if (containerHasElement(plot_var_names, "lat_m")) {
            if (lat_m[lev]) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = lat_m[lev]->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = data(i,j,0);
                    });
                }
            }
            mf_comp++;
        } // lat_m

        if (containerHasElement(plot_var_names, "lon_m")) {
            if (lon_m[lev]) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>& derdat = mf[lev].array(mfi);
                    const Array4<Real>& data   = lon_m[lev]->array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = data(i,j,0);
                    });
                }
            } else {
                mf[lev].setVal(0.0,mf_comp,1,0);
            }

            mf_comp++;

        } // lon_m

        ///////////////////////////////////////////////////////////////////////
        // These quantities are diagnosed by the surface layer
        if (containerHasElement(plot_var_names, "u_star")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_u_star(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // ustar

        if (containerHasElement(plot_var_names, "w_star")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_w_star(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // wstar

        if (containerHasElement(plot_var_names, "t_star")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_t_star(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // tstar

        if (containerHasElement(plot_var_names, "q_star")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_q_star(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // qstar

        if (containerHasElement(plot_var_names, "Olen")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_olen(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // Olen

        if (containerHasElement(plot_var_names, "pblh")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_pblh(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // pblh

        if (containerHasElement(plot_var_names, "t_surf")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& tsurf  = m_SurfaceLayer->get_t_surf(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = tsurf(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // tsurf

        if (containerHasElement(plot_var_names, "q_surf")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_q_surf(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // qsurf

        if (containerHasElement(plot_var_names, "z0")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& ustar  = m_SurfaceLayer->get_z0(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                       derdat(i, j, k, mf_comp) = ustar(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // z0

        if (containerHasElement(plot_var_names, "OLR")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (solverChoice.rad_type != RadiationType::None) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& olr    = rad_fluxes[lev]->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        derdat(i, j, k, mf_comp) = olr(i, j, khi, 2);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // OLR

        if (containerHasElement(plot_var_names, "sens_flux")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (SFS_hfx3_lev[lev]) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& hfx_arr = SFS_hfx3_lev[lev]->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        derdat(i, j, k, mf_comp) = hfx_arr(i, j, klo);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // sens_flux

        // Keep the legacy output name "laten_flux"; it maps to the vertical
        // water-vapor surface flux field.
        if (containerHasElement(plot_var_names, "laten_flux")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (SFS_q1fx3_lev[lev]) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& qfx_arr = SFS_q1fx3_lev[lev]->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        derdat(i, j, k, mf_comp) = qfx_arr(i, j, klo);
                    });
                }
            } else {
                mf[lev].setVal(-999,mf_comp,1,0);
            }
            mf_comp++;
        } // laten_flux

        if (containerHasElement(plot_var_names, "surf_pres")) {
            bool moist = (solverChoice.moisture_type != MoistureType::None);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const auto& derdat   = mf[lev].array(mfi);
                const auto& cons_arr = vars_new[lev][Vars::cons].const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    auto rt = cons_arr(i,j,klo,RhoTheta_comp);
                    auto qv = (moist) ? cons_arr(i,j,klo,RhoQ1_comp)/cons_arr(i,j,klo,Rho_comp)
                                      : zero;
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rt, qv);
                });
            }
            mf_comp++;
        } // surf_pres

        if (containerHasElement(plot_var_names, "integrated_qv")) {
            MultiFab mf_qv_int(mf[lev],make_alias,mf_comp,1);
            if (solverChoice.moisture_type != MoistureType::None) {
                volWgtColumnSum(lev, vars_new[lev][Vars::cons], RhoQ1_comp, mf_qv_int, *detJ_cc[lev]);
            } else {
                mf_qv_int.setVal(0.);
            }
            mf_comp++;
        }

        if (containerHasElement(plot_var_names, "surface_diagnostic_source")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (m_SurfaceLayer) {
                for (MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& source = m_SurfaceLayer->get_surface_diagnostic_source(lev)->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        derdat(i, j, k, mf_comp) = source(i, j, 0);
                    });
                }
            } else {
                mf[lev].setVal(-999, mf_comp, 1, 0);
            }
            mf_comp++;
        } // surface_diagnostic_source

        if (containerHasElement(plot_var_names, "sensible_heat_flux")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (SFS_hfx3_lev[lev]) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& hfx_arr = SFS_hfx3_lev[lev]->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Delegate unit semantics to the surface flux diagnostics helper.
                        derdat(i, j, k, mf_comp) =
                            surface_flux_diagnostics::sensible_heat_flux_wm2_from_rhotheta_flux(
                                hfx_arr(i, j, klo));
                    });
                }
            } else {
                mf[lev].setVal(-999, mf_comp, 1, 0);
            }
            mf_comp++;
        } // sensible_heat_flux

        if (containerHasElement(plot_var_names, "latent_heat_flux")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            if (SFS_q1fx3_lev[lev]) {
                for ( MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    const auto& derdat = mf[lev].array(mfi);
                    const auto& qfx_arr = SFS_q1fx3_lev[lev]->const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                        // Delegate unit semantics to the surface flux diagnostics helper.
                        derdat(i, j, k, mf_comp) =
                            surface_flux_diagnostics::latent_heat_flux_wm2_from_rhoqv_flux(
                                qfx_arr(i, j, klo));
                    });
                }
            } else {
                mf[lev].setVal(-999, mf_comp, 1, 0);
            }
            mf_comp++;
        } // latent_heat_flux

        if (mf_comp != ncomp_mf) {
            Abort(plotfile2d::format_2d_component_count_error(lev, mf_comp, ncomp_mf));
        }
    } // lev

    Vector<Geometry> my_geom = make_2d_plot_geometries(geom, finest_level);

    if (plotfile_type == PlotFileType::Amrex)
    {
        Print() << "Writing 2D native plotfile " << plotfilename << "\n";
        WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                GetVecOfConstPtrs(mf),
                                varnames, my_geom, t_new[0], istep, refRatio());
        writeJobInfo(plotfilename);

#ifdef ERF_USE_NETCDF
    } else if (plotfile_type == PlotFileType::Netcdf) {
         int lev   = 0;
         int l_which = 0;
         const Real* p_lo = my_geom[lev].ProbLo();
         const Real* p_hi = my_geom[lev].ProbHi();
         const auto dx    = my_geom[lev].CellSize();
         writeNCPlotFile(lev, l_which, plotfilename, GetVecOfConstPtrs(mf), varnames, istep,
                         {p_lo[0],p_lo[1],p_lo[2]},{p_hi[0],p_hi[1],dx[2]}, {dx[0],dx[1],dx[2]},
                         my_geom[lev].Domain(), t_new[0], start_bdy_time, solverChoice, zlevels_stag[lev]);
#endif
    } else {
        // Here we assume the plotfile_type is PlotFileType::None
        Print() << "Writing no 2D plotfile since plotfile_type is none" << std::endl;
    }
}
