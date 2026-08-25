/**
 * \file ERF_Plotfile2D.cpp
 */
#include "ERF.H"
#include "ERF_Plotfile2DCatalog.H"
#include "ERF_Plotfile2DFill.H"
#include "ERF_Plotfile2DMetadata.H"
#include "ERF_Plotfile2DInterpolator.H"
#include "ERF_Plotfile2DPrecip.H"
#include "ERF_Plotfile2DWaterPath.H"
#include "ERF_NCPlotFile.H"
#include "ERF_Plotfile2DUtils.H"
#include "Diagnostics/ERF_NearSurfaceDiagnostics.H"
#include "ERF_EpochTime.H"
#include "ERF_SrcHeaders.H"
#include "ERF_StormDiagnostics.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_Utils.H"
#include "Diagnostics/ERF_SeaLevelPressure.H"

using namespace amrex;

namespace
{

// ERF's 2D plotfile path is intentionally an output assembly layer: it
// validates the requested names, assembles the slab geometry, fills existing
// diagnostics, and dispatches to AMReX or NetCDF output. Nontrivial science
// diagnostics should live in dedicated modules and be called from here.
// The assembly path should stay explicit about three diagnostic shapes:
// surface/single-level fields, column reductions, and future interpolated
// horizontal surfaces.
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

bool is_unified_near_surface_diagnostic (plotfile2d::DiagnosticID id) noexcept
{
    return id == plotfile2d::DiagnosticID::Temperature2m ||
           id == plotfile2d::DiagnosticID::WaterVaporMixingRatio2m ||
           id == plotfile2d::DiagnosticID::NearSurfaceDiagnosticSource;
}

bool is_land_surface_diagnostic (const plotfile2d::DiagnosticDescriptor* descriptor) noexcept
{
    return descriptor && descriptor->category == plotfile2d::DiagnosticCategory::LandSurface;
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

    const bool has_surface_layer =
        phys_bc_type[Orientation(Direction::z, Orientation::low)] == ERF_BC::surface_layer;
    const auto active_lsm_names = lsm.Get_DataNames();
    const auto available_names = plotfile2d::available_diagnostic_names(solverChoice,
                                                                         has_surface_layer,
                                                                         active_lsm_names);

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

    const auto output_descriptors =
        plotfile2d::build_sampled_level_output_descriptors(pp_prefix, which,
                                                           plot_var_names, solverChoice);

    Vector<std::string> varnames;
    varnames.reserve(output_descriptors.size());
    for (const auto& descriptor : output_descriptors) {
        varnames.push_back(descriptor.name);
    }
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

        const MultiFab* pblh_source = nullptr;
        const MultiFab* sens_flux_source = SFS_hfx3_lev[lev].get();
        const MultiFab* laten_flux_source = SFS_q1fx3_lev[lev].get();
        const MultiFab* shoc_ustar_source = nullptr;
        const MultiFab* shoc_olen_source = nullptr;
        const MultiFab* shoc_wthv_source = nullptr;
        const ShocDriver* native_shoc = native_shoc_driver[lev].get();
        const bool native_shoc_owns_scalar_fluxes =
            native_shoc && native_shoc->owns_scalar_surface_fluxes();
        const bool native_shoc_has_consumed_flux_diagnostics =
            native_shoc && native_shoc->has_consumed_surface_flux_diagnostics();
        if (native_shoc && native_shoc->has_native_diagnostics()) {
            pblh_source = &native_shoc->pblh_diagnostics();
        }
        // Native SHOC state_update clears the host SFS arrays after consuming
        // them. Use SHOC's preserved snapshots only for flux components whose
        // corresponding host SFS field existed; otherwise keep the source null
        // so the 2D writer emits the documented -999 missing value.
        if (plotfile2d::use_native_shoc_consumed_flux_source(
                native_shoc_owns_scalar_fluxes,
                native_shoc_has_consumed_flux_diagnostics,
                SFS_hfx3_lev[lev] != nullptr)) {
            sens_flux_source = &native_shoc->consumed_sens_flux_diagnostics();
        }
        if (plotfile2d::use_native_shoc_consumed_flux_source(
                native_shoc_owns_scalar_fluxes,
                native_shoc_has_consumed_flux_diagnostics,
                SFS_q1fx3_lev[lev] != nullptr)) {
            laten_flux_source = &native_shoc->consumed_laten_flux_diagnostics();
        }
        if (native_shoc && native_shoc->has_native_diagnostics()) {
            shoc_ustar_source = &native_shoc->shoc_ustar_diagnostics();
            shoc_olen_source = &native_shoc->shoc_olen_diagnostics();
            shoc_wthv_source = &native_shoc->wthv_sec_diagnostics();
        }
        // pblh should follow the active PBL diagnostic provider. Native SHOC
        // diagnoses its own PBL height in state_update mode; SurfaceLayer
        // remains the fallback for non-SHOC configurations.
        //
        // NOTE: the SurfaceLayer pblh holds only the bogus value it was initialized
        //       with unless erf.most.pblh_calc was set (the default is "none"), so we
        //       leave the source null in that case and emit the missing value.
        if (!pblh_source && m_SurfaceLayer && m_SurfaceLayer->computes_pblh()) {
            pblh_source = m_SurfaceLayer->get_pblh(lev);
        }

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
        // Surface and single-level diagnostics use the fill helpers below.
        // Column reductions and future interpolated surfaces keep explicit
        // assembly paths until their contracts are isolated in dedicated
        // helpers.
        if (containerHasElement(plot_var_names, "u_star")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_u_star(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // u_star

        if (containerHasElement(plot_var_names, "w_star")) {
            // NOTE: w_star holds only the bogus value it was initialized with unless
            //       erf.most.include_wstar is on (it is off by default), so we must
            //       pass a null source in that case and emit the missing value
            const MultiFab* w_star_source = (m_SurfaceLayer && m_SurfaceLayer->computes_w_star())
                                          ? m_SurfaceLayer->get_w_star(lev) : nullptr;
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, w_star_source, 0, -999);
            mf_comp++;
        } // w_star

        if (containerHasElement(plot_var_names, "t_star")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_t_star(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // t_star

        if (containerHasElement(plot_var_names, "q_star")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_q_star(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // q_star

        if (containerHasElement(plot_var_names, "Olen")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_olen(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // Olen

        if (containerHasElement(plot_var_names, "pblh")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, pblh_source, 0, -999);
            mf_comp++;
        } // pblh

        if (containerHasElement(plot_var_names, "t_surf")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_t_surf(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // t_surf

        if (containerHasElement(plot_var_names, "q_surf")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_q_surf(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // q_surf

        if (containerHasElement(plot_var_names, "z0")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, m_SurfaceLayer ? m_SurfaceLayer->get_z0(lev) : nullptr,
                0, -999);
            mf_comp++;
        } // z0

        if (containerHasElement(plot_var_names, "OLR")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, rad_fluxes[lev].get(), khi, -999, 2);
            mf_comp++;
        } // OLR

        if (containerHasElement(plot_var_names, "sens_flux")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, sens_flux_source, klo, -999);
            mf_comp++;
        } // sens_flux

        // Keep the legacy output name "laten_flux"; it maps to the vertical
        // water-vapor surface flux field.
        if (containerHasElement(plot_var_names, "laten_flux")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, laten_flux_source, klo, -999);
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

        if (containerHasElement(plot_var_names, "sea_level_pressure")) {
            sea_level_pressure_diagnostics::fill_sea_level_pressure(
                mf[lev], mf_comp, vars_new[lev][Vars::cons], *z_phys_nd[lev],
                solverChoice.moisture_indices.qv, klo);
            mf_comp++;
        } // sea_level_pressure

        const auto precip_sources =
            micro ? micro->Get_Surface_Precip_Accumulation_Ptrs(lev)
                  : SurfacePrecipAccumulationSources{};
        const auto selected_precipitation =
            plotfile2d::selected_precipitation_accumulation_components(varnames,
                                                                       precip_sources);
        if (selected_precipitation.n > 0) {
            plotfile2d::fill_precipitation_accumulations(mf[lev],
                                                         precip_sources,
                                                         selected_precipitation,
                                                         klo);
        }

        if (containerHasElement(plot_var_names, "precip_total_accum"))   { mf_comp++; }
        if (containerHasElement(plot_var_names, "precip_rain_accum"))    { mf_comp++; }
        if (containerHasElement(plot_var_names, "precip_snow_accum"))    { mf_comp++; }
        if (containerHasElement(plot_var_names, "precip_graupel_accum")) { mf_comp++; }
        if (containerHasElement(plot_var_names, "precip_hail_accum"))    { mf_comp++; }
        if (containerHasElement(plot_var_names, "precip_frozen_accum"))  { mf_comp++; }

        if (containerHasElement(plot_var_names, "integrated_qv")) {
            // integrated_qv remains the legacy column-reduction example.
            // Scheme-aware condensed water paths use a dedicated helper below
            // rather than the single-k copy helpers.
            MultiFab mf_qv_int(mf[lev],make_alias,mf_comp,1);
            if (solverChoice.moisture_type != MoistureType::None) {
                volWgtColumnSum(lev, vars_new[lev][Vars::cons], RhoQ1_comp, mf_qv_int, *detJ_cc[lev]);
            } else {
                mf_qv_int.setVal(0.);
            }
            mf_comp++;
        }

        const auto selected_water_paths =
            plotfile2d::selected_condensed_water_path_components(varnames, solverChoice);
        if (selected_water_paths.n > 0) {
            // Condensed water paths are scheme-aware column reductions. Their
            // availability comes from solverChoice.moisture_indices, and the
            // reduction uses the same metric convention as integrated_qv.
            plotfile2d::fill_condensed_water_paths(mf[lev],
                                                   vars_new[lev][Vars::cons],
                                                   selected_water_paths,
                                                   geom[lev],
                                                   *detJ_cc[lev]);
        }

        if (containerHasElement(plot_var_names, "integrated_qc")) { mf_comp++; }
        if (containerHasElement(plot_var_names, "integrated_qi")) { mf_comp++; }
        if (containerHasElement(plot_var_names, "integrated_qr")) { mf_comp++; }
        if (containerHasElement(plot_var_names, "integrated_qs")) { mf_comp++; }
        if (containerHasElement(plot_var_names, "integrated_qg")) { mf_comp++; }

        if (containerHasElement(plot_var_names, "surface_diagnostic_source")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp,
                m_SurfaceLayer ? m_SurfaceLayer->get_surface_diagnostic_source(lev) : nullptr,
                0, 0);
            mf_comp++;
        } // surface_diagnostic_source

        if (containerHasElement(plot_var_names, "sensible_heat_flux")) {
            plotfile2d::fill_sensible_heat_flux_from_klevel_or_missing(
                mf[lev], mf_comp, sens_flux_source, klo, -999);
            mf_comp++;
        } // sensible_heat_flux

        if (containerHasElement(plot_var_names, "latent_heat_flux")) {
            plotfile2d::fill_latent_heat_flux_from_klevel_or_missing(
                mf[lev], mf_comp, laten_flux_source, klo, -999);
            mf_comp++;
        } // latent_heat_flux

        if (containerHasElement(plot_var_names, "shoc_u_star")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, shoc_ustar_source, klo, -999);
            mf_comp++;
        } // shoc_u_star

        if (containerHasElement(plot_var_names, "shoc_Olen")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, shoc_olen_source, klo, -999);
            mf_comp++;
        } // shoc_Olen

        if (containerHasElement(plot_var_names, "shoc_wthv_sfc")) {
            plotfile2d::fill_component_from_klevel_or_value(
                mf[lev], mf_comp, shoc_wthv_source, klo, -999);
            mf_comp++;
        } // shoc_wthv_sfc

        // Land-surface provider fields use the generic LandSurface name/index
        // interface. The catalog supplies metadata; this block only assembles
        // stored values and translates provider sentinels.
        int temperature_2m_comp = -1;
        int mixing_ratio_2m_comp = -1;
        int near_surface_source_comp = -1;
        for (const auto& name : plot_var_names) {
            const auto* descriptor = plotfile2d::find_diagnostic(name);
            if (!descriptor) {
                descriptor = plotfile2d::find_dynamic_soil_diagnostic(name);
            }
            if (!is_land_surface_diagnostic(descriptor)) {
                continue;
            }

            if (is_unified_near_surface_diagnostic(descriptor->id)) {
                if (descriptor->id == plotfile2d::DiagnosticID::Temperature2m) {
                    temperature_2m_comp = mf_comp;
                } else if (descriptor->id == plotfile2d::DiagnosticID::WaterVaporMixingRatio2m) {
                    mixing_ratio_2m_comp = mf_comp;
                } else {
                    near_surface_source_comp = mf_comp;
                }
                ++mf_comp;
                continue;
            }

            plotfile2d::fill_land_surface_component_from_klevel_or_missing(
                mf[lev], mf_comp, lsm.Get_Data_Ptr(lev, name), 0, -999);
            ++mf_comp;
        }

        if (temperature_2m_comp >= 0 || mixing_ratio_2m_comp >= 0 ||
            near_surface_source_comp >= 0) {
            near_surface_diagnostics::Sources near_surface_sources;
            near_surface_sources.native_temperature_vegetated =
                lsm.Get_Data_Ptr(lev, "noahmp_temperature_2m_vegetated");
            near_surface_sources.native_temperature_bare =
                lsm.Get_Data_Ptr(lev, "noahmp_temperature_2m_bare");
            near_surface_sources.native_mixing_ratio_vegetated =
                lsm.Get_Data_Ptr(lev, "noahmp_water_vapor_mixing_ratio_2m_vegetated");
            near_surface_sources.native_mixing_ratio_bare =
                lsm.Get_Data_Ptr(lev, "noahmp_water_vapor_mixing_ratio_2m_bare");
            near_surface_sources.native_vegetation_fraction =
                lsm.Get_Data_Ptr(lev, "noahmp_vegetation_fraction");
            near_surface_sources.theta_surface =
                m_SurfaceLayer ? m_SurfaceLayer->get_t_surf(lev) : nullptr;
            near_surface_sources.theta_star =
                m_SurfaceLayer ? m_SurfaceLayer->get_t_star(lev) : nullptr;
            near_surface_sources.mixing_ratio_surface =
                m_SurfaceLayer ? m_SurfaceLayer->get_q_surf(lev) : nullptr;
            near_surface_sources.mixing_ratio_star =
                m_SurfaceLayer ? m_SurfaceLayer->get_q_star(lev) : nullptr;
            near_surface_sources.roughness_height =
                m_SurfaceLayer ? m_SurfaceLayer->get_z0(lev) : nullptr;
            near_surface_sources.obukhov_length =
                m_SurfaceLayer ? m_SurfaceLayer->get_olen(lev) : nullptr;
            near_surface_sources.source_mask =
                m_SurfaceLayer ? m_SurfaceLayer->get_surface_diagnostic_source(lev) : nullptr;
            near_surface_sources.land_mask =
                (lmask_lev[lev].empty()) ? nullptr : lmask_lev[lev][0].get();
            near_surface_sources.cons = &vars_new[lev][Vars::cons];
            near_surface_sources.z_phys_nd = z_phys_nd[lev].get();
            near_surface_sources.dz = geom[lev].CellSize(2);
            near_surface_sources.klo = klo;
            near_surface_sources.moist = solverChoice.moisture_type != MoistureType::None;
            near_surface_sources.has_lsm =
                near_surface_sources.native_temperature_vegetated != nullptr &&
                near_surface_sources.native_temperature_bare != nullptr &&
                near_surface_sources.native_mixing_ratio_vegetated != nullptr &&
                near_surface_sources.native_mixing_ratio_bare != nullptr &&
                near_surface_sources.native_vegetation_fraction != nullptr;
            near_surface_diagnostics::fill(mf[lev], temperature_2m_comp,
                                           mixing_ratio_2m_comp,
                                           near_surface_source_comp,
                                           near_surface_sources);
        }

        const int static_output_count = static_cast<int>(plot_var_names.size());
        // If this level is anelastic then the base state pressure -- not the compressible
        //    EOS -- defines the pressure, both for the sampled "pressure" field and for
        //    the pressure vertical coordinate.  This matches the 3-D plotfile path.
        MultiFab p_hse(base_state[lev], make_alias, BaseState::p0_comp, 1);
        const MultiFab* p_hse_ptr = (solverChoice.anelastic[lev] == 1) ? &p_hse : nullptr;
        for (int out_idx = static_output_count; out_idx < static_cast<int>(output_descriptors.size()); ++out_idx) {
            const auto& descriptor = output_descriptors[out_idx];
            // Sampled-level descriptors are dynamic. The interpolator owns
            // their vertical sampling so the writer remains an output
            // assembly layer.
            plotfile2d::SampledWindSources wind_sources;
            wind_sources.xvel = &vars_new[lev][Vars::xvel];
            wind_sources.yvel = &vars_new[lev][Vars::yvel];
            wind_sources.zvel = &vars_new[lev][Vars::zvel];
            plotfile2d::fill_sampled_level_component(
                mf[lev], mf_comp, descriptor,
                vars_new[lev][Vars::cons],
                z_phys_cc[lev].get(),
                *z_phys_nd[lev],
                z_phys_cc[lev] != nullptr,
                solverChoice.moisture_indices,
                klo, khi,
                wind_sources,
                p_hse_ptr);
            mf_comp++;
        }

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
                                varnames, my_geom, static_cast<Real>(t_new[0]), istep, refRatio());
        // Native AMReX 2D plotfiles write a JSON sidecar with catalog
        // metadata for the selected output variables only.
        plotfile2d::write_2d_metadata_json(plotfilename, output_descriptors);
        writeJobInfo(plotfilename, erf_provenance::ArtifactType::Plotfile2D,
                     istep[0], t_new[0]);

#ifdef ERF_USE_NETCDF
    } else if (plotfile_type == PlotFileType::Netcdf) {
         int lev   = 0;
         int l_which = 0;
         const Real* p_lo = my_geom[lev].ProbLo();
         const Real* p_hi = my_geom[lev].ProbHi();
         const auto dx    = my_geom[lev].CellSize();
         writeNCPlotFile(lev, l_which, plotfilename, GetVecOfConstPtrs(mf), varnames, istep,
                         {p_lo[0],p_lo[1],p_lo[2]},{p_hi[0],p_hi[1],dx[2]}, {dx[0],dx[1],dx[2]},
                         my_geom[lev].Domain(), static_cast<Real>(t_new[0]),
                         static_cast<Real>(start_bdy_time), solverChoice, zlevels_stag[lev]);
#endif
    } else {
        // Here we assume the plotfile_type is PlotFileType::None
        Print() << "Writing no 2D plotfile since plotfile_type is none" << std::endl;
    }
}
