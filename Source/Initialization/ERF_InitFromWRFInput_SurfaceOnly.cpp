/**
 * \file ERF_InitFromWRFInput_SurfaceOnly.cpp
 */

#include "ERF.H"
#include "ERF_EOS.H"
#include "ERF_Constants.H"
#include "ERF_Utils.H"
#include "ERF_ProbCommon.H"
#include "ERF_DataStruct.H"
#include "ERF_TerrainMetrics.H"

#include "ERF_ReadFromWRFInput.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

// Forward declarations for helper functions defined in ERF_InitFromWRFInput.cpp
Box read_subdomain_from_wrfinput(int lev, const std::string& fname, int& ratio);
Real compute_terrain_top_and_bottom(const MultiFab& mf_PH, const MultiFab& mf_PHB, const Box& domain);
void init_terrain_from_wrfinput(int lev, Geometry& geom, const Real& z_top, const Box& subdomain,
                                MultiFab* z_phys_nd, const MultiFab& NC_PH_fab, const MultiFab& NC_PHB_fab,
                                Real& dz0_max, const bool& avg_grid_faces_to_nodes);
void refine_fab_in_z(const FArrayBox& crse_fab, FArrayBox& fine_fab, int rr_z);
void read_base_state_params_from_wrfinput(const std::string& fname, BaseStateParams& bsp);

/**
 * Initialize surface fields from a WRF dataset (for use with atmospheric interpolation from coarse).
 *
 * This function reads only surface fields from wrfinput (terrain, surface temps, land masks,
 * map factors, LSM variables) while leaving atmospheric state (U, V, W, THM, moisture, density)
 * to be filled by interpolation from the coarse level via FillCoarsePatch.
 *
 * @param lev Integer specifying the current level
 * @param mf_PSFC MultiFab storing surface pressure for this level
 */
void
ERF::init_from_wrfinput_surface_only (int lev, MultiFab& mf_PSFC_lev)
{
    if (nc_init_file.empty()) {
        amrex::Error("NetCDF initialization file name must be provided via input");
    }

    bool use_lsm = (solverChoice.lsm_type != LandSurfaceType::None);

    // List of surface-only variables to read
    Vector<std::string> NC_names;
    NC_names.push_back("PH");        // Terrain
    NC_names.push_back("PHB");       // Terrain base state
    NC_names.push_back("PSFC");      // Surface pressure
    NC_names.push_back("MUB");       // Base state column mass
    NC_names.push_back("MAPFAC_U");  // Map factor U
    NC_names.push_back("MAPFAC_V");  // Map factor V
    NC_names.push_back("MAPFAC_M");  // Map factor M
    NC_names.push_back("SST");       // Sea surface temperature
    NC_names.push_back("TSK");       // Skin temperature
    NC_names.push_back("LANDMASK");  // Land mask
    NC_names.push_back("C1H");       // WRF vertical coordinate
    NC_names.push_back("C2H");       // WRF vertical coordinate
    NC_names.push_back("RDNW");      // WRF vertical coordinate
    NC_names.push_back("XLAT_V");    // Latitude
    NC_names.push_back("XLONG_U");   // Longitude
    NC_names.push_back("IVGTYP");    // Vegetation type
    NC_names.push_back("ISLTYP");    // Soil type

    // Add LSM variables if using land surface model
    if (use_lsm) {
        NC_names.push_back("TSLB");   // Soil temperature
        NC_names.push_back("SMOIS");  // Soil moisture
        NC_names.push_back("SH2O");   // Liquid soil moisture
        NC_names.push_back("LAI");    // Leaf area index
        NC_names.push_back("ZS");     // Soil layer depths
        NC_names.push_back("DZS");    // Soil layer thickness
        NC_names.push_back("VEGFRA"); // Vegetation fraction
        NC_names.push_back("TMN");    // Deep soil temperature
        NC_names.push_back("SHDMIN"); // Minimum green vegetation fraction
        NC_names.push_back("SHDMAX"); // Maximum green vegetation fraction
    }

    int nvar = NC_names.size();
    Vector<Vector<FArrayBox>> NC_fab_var;
    NC_fab_var.resize(num_boxes_at_level[lev]);
    for (int idx(0); idx < num_boxes_at_level[lev]; ++idx) {
        NC_fab_var[idx].resize(nvar);
    }

    auto& lev_new = vars_new[lev];

    // Temporary MultiFabs
    MultiFab* mf_C1H;
    MultiFab* mf_C2H;
    MultiFab* mf_RDNW;
    MultiFab* mf_MUB;
    MultiFab* mf_PHB;
    MultiFab  C1H_tmp;
    MultiFab  C2H_tmp;
    MultiFab  RDNW_tmp;
    MultiFab  MUB_tmp;
    MultiFab  PHB_tmp;

    MultiFab mf_PH;

    // Read base state parameters from level 0 file (same logic as full init)
    if (lev == 0) {
        read_base_state_params_from_wrfinput(nc_init_file[lev][0], wrf_bsp);
        wrf_bsp.set_layer_interfaces();
    } else {
        BaseStateParams lev_bsp;
        read_base_state_params_from_wrfinput(nc_init_file[lev][0], lev_bsp);
        if (!lev_bsp.same_params_as(wrf_bsp)) {
            Print() << "WARNING: the base state parameters in " << nc_init_file[lev][0]
                    << " differ from those at level 0; using level 0 values.\n";
        }
    }
    AMREX_ALWAYS_ASSERT(wrf_bsp.is_set);

    // Check for vertical refinement
    int rr_z = 1;
    for (int l = 0; l < lev; ++l) { rr_z *= ref_ratio[l][2]; }
    if (rr_z > 1) {
        Print() << "Level " << lev << " is refined by " << rr_z << " in the vertical relative to "
                << "the wrfinput files.\n";
    }

    auto& ba    = lev_new[Vars::cons].boxArray();
    auto& dm    = lev_new[Vars::cons].DistributionMap();
    IntVect ng  = lev_new[Vars::cons].nGrowVect();
    IntVect ngz = (z_phys_nd[lev]) ? z_phys_nd[lev]->nGrowVect() : IntVect(0);
    ngz[0] += 1; ngz[1] += 1;
    IntVect ngv = ng; ngv[2] = 0;

    bool compute_terrain_here = true;
    const Real l_rdOcp = solverChoice.rdOcp;

    Print() << "Loading surface data from NetCDF file at level " << lev << "\n";
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
        Print() << "Reading from file " << nc_init_file[lev][idx] << "\n";

        int ratio_from_file;
        Box subdomain_to_read = read_subdomain_from_wrfinput(lev, nc_init_file[lev][idx], ratio_from_file);
        Print() << "Box in file " << subdomain_to_read << "\n";

        Box subdomain_to_fill = boxes_at_level[lev][idx];
        Print() << "Box to fill " << subdomain_to_fill << "\n";

        for (int ivar = 0; ivar < nvar; ++ivar) {
            auto var_name = NC_names[ivar];
            Print() << "Checking for " << var_name << " ...";

            int success, use_theta_m;
            read_from_wrfinput(lev, subdomain_to_read, nc_init_file[lev][idx],
                               NC_fab_var[idx][ivar], var_name, geom[lev],
                               use_theta_m, success);

            auto& var_fab_from_file = NC_fab_var[idx][ivar];

            if (!success) {
                Print() << " not found (skipping)\n";
                continue;
            }

            FArrayBox var_fab;
            FArrayBox var_fab_crse;

            // Handle shifting for level > 1
            if (success && lev > 1) {
                Box shift_by_box(subdomains[lev][0].minimalBox());
                IntVect shift_by(shift_by_box.smallEnd());
                for (int i = 0; i < AMREX_SPACEDIM; i++) {
                    shift_by[i] -= var_fab_from_file.box().smallEnd(i);
                }
                if (rr_z > 1) {
                    shift_by[2] = amrex::coarsen(shift_by_box.smallEnd(2), rr_z)
                                - var_fab_from_file.box().smallEnd(2);
                }
                var_fab_from_file.shift(shift_by);
            }

            // Handle dimension reduction for 1D/2D variables
            int nx = var_fab_from_file.box().length(0);
            int ny = var_fab_from_file.box().length(1);
            int nz = var_fab_from_file.box().length(2);
            Box subdomain_tmp(subdomain_to_fill);
            if (nx == 1 && ny == 1) {
                subdomain_tmp.setBig(0, subdomain_tmp.smallEnd(0));
                subdomain_tmp.setBig(1, subdomain_tmp.smallEnd(1));
            }
            if (nz == 1) {
                subdomain_tmp.setBig(2, subdomain_tmp.smallEnd(2));
            }

            Box subdomain_to_fill_typed(convert(subdomain_tmp, var_fab_from_file.box().ixType()));

            // Handle staggered XLONG_U and XLAT_V
            if (var_name == "XLONG_U" &&
                var_fab_from_file.box().bigEnd(0) > subdomain_to_fill_typed.bigEnd(0)) {
                subdomain_to_fill_typed.growHi(0, 1);
            }
            if (var_name == "XLAT_V" &&
                var_fab_from_file.box().bigEnd(1) > subdomain_to_fill_typed.bigEnd(1)) {
                subdomain_to_fill_typed.growHi(1, 1);
            }

            // Check if vertical refinement is needed
            const bool is_soil_var = (var_name == "TSLB" || var_name == "SMOIS" ||
                                      var_name == "SH2O" || var_name == "ZS" ||
                                      var_name == "DZS");
            const bool refine_in_z = (rr_z > 1) && (nz > 1) && !is_soil_var;

            Box subdomain_crse(subdomain_to_fill_typed);
            if (refine_in_z) {
                subdomain_crse.coarsen(IntVect(1, 1, rr_z));
                subdomain_crse.grow(2, 1);
                subdomain_crse.setSmall(2, amrex::max(subdomain_crse.smallEnd(2),
                                                      var_fab_from_file.box().smallEnd(2)));
                subdomain_crse.setBig(2, amrex::min(subdomain_crse.bigEnd(2),
                                                    var_fab_from_file.box().bigEnd(2)));
            }

#ifdef AMREX_USE_GPU
            var_fab.resize(subdomain_to_fill_typed, 1, amrex::The_Pinned_Arena());
            if (refine_in_z) { var_fab_crse.resize(subdomain_crse, 1, amrex::The_Pinned_Arena()); }
#else
            var_fab.resize(subdomain_to_fill_typed, 1);
            if (refine_in_z) { var_fab_crse.resize(subdomain_crse, 1); }
#endif

            FArrayBox& read_fab = (refine_in_z) ? var_fab_crse : var_fab;

            Box intersection = read_fab.box() & var_fab_from_file.box();
            if (intersection.ok()) {
                if (refine_in_z) {
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(intersection.contains(read_fab.box()),
                                                     "ERF::init_from_wrfinput_surface_only: file doesn't "
                                                     "cover region needed for vertical refinement");
                }
                read_fab.template copy<RunOn::Device>(var_fab_from_file, intersection, 0, intersection, 0, 1);
            } else if (nx == 1 && ny == 1) {
                // Handle column data
                Print() << " Copying 1D FAB from " << var_fab_from_file.box()
                        << " to " << read_fab.box() << std::endl;
                IntVect shift_ij(AMREX_D_DECL(read_fab.box().smallEnd(0) - var_fab_from_file.box().smallEnd(0),
                                              read_fab.box().smallEnd(1) - var_fab_from_file.box().smallEnd(1),
                                              0));
                var_fab_from_file.shift(shift_ij);
                Box isect_1d = read_fab.box() & var_fab_from_file.box();
                if (isect_1d.ok()) {
                    read_fab.template copy<RunOn::Device>(var_fab_from_file, isect_1d, 0, isect_1d, 0, 1);
                }
                var_fab_from_file.shift(-shift_ij);
            } else {
                amrex::Error("ERF::init_from_wrfinput_surface_only: Region we want not contained in region we have");
            }

            if (refine_in_z) {
                refine_fab_in_z(var_fab_crse, var_fab, rr_z);
            }

            // Now process each specific variable (keeping all the original processing logic)
            bool lat_periodic = (geom[lev].isPeriodic(0) && geom[lev].isPeriodic(1));
            int i_lo = boxes_at_level[lev][0].smallEnd(0);
            int i_hi = boxes_at_level[lev][0].bigEnd(0);
            int j_lo = boxes_at_level[lev][0].smallEnd(1);
            int j_hi = boxes_at_level[lev][0].bigEnd(1);

            // PH and PHB for terrain
            if (var_name == "PH") {
                auto& ba_w = lev_new[Vars::zvel].boxArray();
                mf_PH.define(ba_w, dm, 1, IntVect(ngz[0], ngz[1], 0));
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf_PH, false); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<Real>& dst_arr = mf_PH.array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, k) = src_arr(li, lj, k);
                    });
                }
                var_fab.clear();
            } else if (var_name == "PHB") {
                auto& ba_w = lev_new[Vars::zvel].boxArray();
                if (lev == 0) {
                    mf_PHB = wrf_PHB.get();
                } else {
                    PHB_tmp.define(ba_w, dm, 1, IntVect(ngz[0], ngz[1], 0));
                    mf_PHB = &PHB_tmp;
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mf_PHB, false); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<Real>& dst_arr = mf_PHB->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, k) = src_arr(li, lj, k);
                    });
                }
                var_fab.clear();
            } else if (var_name == "PSFC") {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf_PSFC_lev, false); mfi.isValid(); ++mfi) {
                    FArrayBox &cur_fab = mf_PSFC_lev[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                }
                Real pmax = mf_PSFC_lev.max(0);
                if (pmax == zero) {
                    amrex::Print() << " PSFC read in had max of 0; replacing it by 1e5 everywhere" << std::endl;
                    mf_PSFC_lev.setVal(p_0);
                }
                var_fab.clear();
            } else if (var_name == "MUB") {
                if (lev == 0) {
                    mf_MUB = wrf_MUB.get();
                } else {
                    MUB_tmp.define(ba2d[lev], dm, 1, IntVect(ngz[0], ngz[1], 0));
                    mf_MUB = &MUB_tmp;
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mf_MUB, false); mfi.isValid(); ++mfi) {
                    FArrayBox &cur_fab = (*mf_MUB)[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                }
                var_fab.clear();
            } else if (var_name == "C1H") {
                if (lev == 0) {
                    mf_C1H = wrf_C1H.get();
                } else {
                    C1H_tmp.define(ba1d[lev], dm, 1, IntVect(ngz[0], ngz[1], 0));
                    mf_C1H = &C1H_tmp;
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mf_C1H, false); mfi.isValid(); ++mfi) {
                    FArrayBox &cur_fab = (*mf_C1H)[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                }
                var_fab.clear();
            } else if (var_name == "C2H") {
                if (lev == 0) {
                    mf_C2H = wrf_C2H.get();
                } else {
                    C2H_tmp.define(ba1d[lev], dm, 1, IntVect(ngz[0], ngz[1], 0));
                    mf_C2H = &C2H_tmp;
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mf_C2H, false); mfi.isValid(); ++mfi) {
                    FArrayBox &cur_fab = (*mf_C2H)[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                }
                var_fab.clear();
            } else if (var_name == "RDNW") {
                if (lev == 0) {
                    mf_RDNW = wrf_RDNW.get();
                } else {
                    RDNW_tmp.define(ba1d[lev], dm, 1, IntVect(ngz[0], ngz[1], 0));
                    mf_RDNW = &RDNW_tmp;
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mf_RDNW, false); mfi.isValid(); ++mfi) {
                    FArrayBox &cur_fab = (*mf_RDNW)[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                }
                var_fab.clear();
            }

            // Latitude & Coriolis factors
            if (var_name == "XLAT_V") {
                int vf_j_hi = var_fab.box().bigEnd(1);
                lat_m[lev] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, ngv);
                sinPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, ngv);
                cosPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, ngv);
                for (MFIter mfi(*(lat_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<Real>& sin_arr = (sinPhi_m[lev])->array(mfi);
                    const Array4<Real>& cos_arr = (cosPhi_m[lev])->array(mfi);
                    const Array4<Real>& dst_arr = (lat_m[lev])->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        int sj = amrex::min(amrex::max(j, j_lo), vf_j_hi);
                        dst_arr(i, j, 0) = src_arr(li, sj, 0);
                        Real lat_rad = src_arr(li, lj, 0) * (PI / Real(180.));
                        sin_arr(i, j, 0) = std::sin(lat_rad);
                        cos_arr(i, j, 0) = std::cos(lat_rad);
                    });
                }
            }

            // Longitude
            if (var_name == "XLONG_U") {
                int vf_i_hi = var_fab.box().bigEnd(0);
                lon_m[lev] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, ngv);
                for (MFIter mfi(*(lon_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<Real>& dst_arr = (lon_m[lev])->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), vf_i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, 0) = src_arr(li, lj, 0);
                    });
                }
            }

            // SST
            if (var_name == "SST") {
                sst_lev[lev][0] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, ngv);
                for (MFIter mfi(*(sst_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<Real>& dst_arr = sst_lev[lev][0]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    const Array4<const Real>& psfc_arr = mf_PSFC_lev.const_array(mfi);
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, 0) = getThgivenTandP(src_arr(li, lj, 0), psfc_arr(li, lj, 0), l_rdOcp);
                    });
                }
                (sst_lev[lev][0])->FillBoundary(geom[lev].periodicity());
            }

            // TSK
            if (var_name == "TSK") {
                tsk_lev[lev][0] = std::make_unique<MultiFab>(ba2d[lev], dm, 1, ngv);
                for (MFIter mfi(*(tsk_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<Real>& dst_arr = tsk_lev[lev][0]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    const Array4<const Real>& psfc_arr = mf_PSFC_lev.const_array(mfi);
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, 0) = getThgivenTandP(src_arr(li, lj, 0), psfc_arr(li, lj, 0), l_rdOcp);
                    });
                }
                (tsk_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
            }

            // Landmask
            if (var_name == "LANDMASK") {
                for (MFIter mfi(*(lmask_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<int>& dst_arr = lmask_lev[lev][0]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, 0) = static_cast<int>(src_arr(li, lj, 0));
                    });
                }
                (lmask_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
            }

            // Landtype
            if (var_name == "IVGTYP") {
                for (MFIter mfi(*(land_type_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<int>& dst_arr = land_type_lev[lev][0]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, 0) = static_cast<int>(src_arr(li, lj, 0));
                    });
                }
                (land_type_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
            }

            // Soil type
            if (var_name == "ISLTYP") {
                for (MFIter mfi(*(soil_type_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    const Array4<int>& dst_arr = soil_type_lev[lev][0]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, i_lo), i_hi);
                        int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                        dst_arr(i, j, 0) = static_cast<int>(src_arr(li, lj, 0));
                    });
                }
                (soil_type_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
            }

            // LSM variables
            if (use_lsm) {
                auto &lsm_wrfmap = lsm.Get_WRFInputNames();
                for (auto &var : lsm_wrfmap) {
                    if (var_name == var.first) {
                        bool is_3d = var_fab.box().length(2) > 1;
                        amrex::Print() << "   Reading " << ((is_3d) ? "3D" : "2D") << " LSM variable '"
                                      << var.first << "' (" << var.second << ")" << std::endl;
                        int lsm_idx = lsm.Get_DataIdx(lev, var.second);
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lsm_idx != -1, "LSM variable mapping invalid!");
                        AMREX_ALWAYS_ASSERT(lsm_data[lev][lsm_idx]);

                        int lsm_nsoil = lsm.Get_Lsm_Geom(lev).Domain().length(2);
                        if (is_3d) {
                            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lsm_nsoil == var_fab.box().length(2),
                                                             "Number of soil layers must match!");
                        }

                        bool is_column = var_fab.box().length(0) == 1 && var_fab.box().length(1) == 1;

                        for (MFIter mfi(*lsm_data[lev][lsm_idx], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                            Box gtbx = mfi.tilebox();
                            int lsm_khi = gtbx.bigEnd(2);
                            gtbx.setRange(2, 0, var_fab.box().length(2));
                            const Array4<Real>& dst_arr = lsm_data[lev][lsm_idx]->array(mfi);
                            const Array4<const Real>& src_arr = var_fab.const_array();
                            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                                if (is_column) {
                                    li = 0;
                                    lj = 0;
                                }
                                const int lsm_k = lsm_khi - k;
                                dst_arr(i, j, lsm_k) = src_arr(li, lj, k);
                            });
                        }
                        (lsm_data[lev][lsm_idx])->FillBoundary(geom[lev].periodicity());
                    }
                }
            }

            // MapFac U
            if (var_name == "MAPFAC_U") {
                Real max_val = var_fab.template max<RunOn::Device>();
                if (std::fabs(max_val) < std::numeric_limits<Real>::epsilon()) {
                    Print() << "MAPFAC_U cannot be 0, resetting to 1!\n";
                    var_fab.template setVal<RunOn::Device>(1);
                }
                if (lat_periodic) {
                    Print() << "MAPFAC_U resetting to 1 with lateral periodic BCs!\n";
                    var_fab.template setVal<RunOn::Device>(1);
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mapfac[lev][MapFacType::u_x], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    Box vbx = mfi.validbox();
                    int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
                    int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
                    const Array4<Real>& dst_arr = mapfac[lev][MapFacType::u_x]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, ilo), ihi);
                        int lj = amrex::min(amrex::max(j, jlo), jhi);
                        dst_arr(i, j, 0) = src_arr(li, lj, 0);
                    });
                }
                mapfac[lev][MapFacType::u_x]->FillBoundary(geom[lev].periodicity());
            }

            // MapFac V
            if (var_name == "MAPFAC_V") {
                Real max_val = var_fab.template max<RunOn::Device>();
                if (std::fabs(max_val) < std::numeric_limits<Real>::epsilon()) {
                    Print() << "MAPFAC_V cannot be 0, resetting to 1!\n";
                    var_fab.template setVal<RunOn::Device>(1);
                }
                if (lat_periodic) {
                    Print() << "MAPFAC_V resetting to 1 with lateral periodic BCs!\n";
                    var_fab.template setVal<RunOn::Device>(1);
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mapfac[lev][MapFacType::v_x], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    Box vbx = mfi.validbox();
                    int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
                    int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
                    const Array4<Real>& dst_arr = mapfac[lev][MapFacType::v_x]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, ilo), ihi);
                        int lj = amrex::min(amrex::max(j, jlo), jhi);
                        dst_arr(i, j, 0) = src_arr(li, lj, 0);
                    });
                }
                mapfac[lev][MapFacType::v_x]->FillBoundary(geom[lev].periodicity());
            }

            // MapFac M
            if (var_name == "MAPFAC_M") {
                Real max_val = var_fab.template max<RunOn::Device>();
                if (std::fabs(max_val) < std::numeric_limits<Real>::epsilon()) {
                    Print() << "MAPFAC_M cannot be 0, resetting to 1!\n";
                    var_fab.template setVal<RunOn::Device>(1);
                }
                if (lat_periodic) {
                    Print() << "MAPFAC_M resetting to 1 with lateral periodic BCs!\n";
                    var_fab.template setVal<RunOn::Device>(1);
                }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*mapfac[lev][MapFacType::m_x], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box gtbx = mfi.growntilebox();
                    Box vbx = mfi.validbox();
                    int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
                    int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
                    const Array4<Real>& dst_arr = mapfac[lev][MapFacType::m_x]->array(mfi);
                    const Array4<const Real>& src_arr = var_fab.const_array();
                    ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                        int li = amrex::min(amrex::max(i, ilo), ihi);
                        int lj = amrex::min(amrex::max(j, jlo), jhi);
                        dst_arr(i, j, 0) = src_arr(li, lj, 0);
                    });
                }
                mapfac[lev][MapFacType::m_x]->FillBoundary(geom[lev].periodicity());
            }

            if (success) {
                var_fab.clear();
            }
        } // ivar
        Print() << "\n";
        have_read_nc_init_file[lev][idx] = 1;
    } // idx

    // Compute terrain from PH and PHB (same as full init)
    if (compute_terrain_here) {
        if (lev == 0) {
            AMREX_ALWAYS_ASSERT(solverChoice.terrain_type == TerrainType::StaticFittedMesh);
            z_top = compute_terrain_top_and_bottom(mf_PH, *mf_PHB, geom[lev].Domain());
        } else {
            amrex::Print() << "Using top of domain set at level 0 which is " << z_top << std::endl;
        }

        // Initialize terrain
        ParmParse pp("erf");
        int terrain_smoothing = 0;
        pp.query("terrain_smoothing", terrain_smoothing);

        FineTerrain fine_terrain = FineTerrain::None;
        MultiFab z_phys_interp;

        if (lev > 0 && terrain_smoothing != 0) {
            fine_terrain = which_fine_terrain();
            if (fine_terrain != FineTerrain::Transform) {
                Abort("terrain_smoothing = " + std::to_string(terrain_smoothing) +
                      " with wrfinput initialization on level > 0 requires "
                      "erf.amr_terrain_refinement = transform (not interpolate)");
            }

            InterpFromCoarseLevel(*z_phys_nd[lev], z_phys_nd[lev]->nGrowVect(),
                                  IntVect(0, 0, 0),
                                  *z_phys_nd[lev-1], 0, 0, 1,
                                  geom[lev-1], geom[lev],
                                  refRatio(lev-1), &node_bilinear_interp,
                                  domain_bcs_type, BCVars::cons_bc);

            z_phys_interp.define(z_phys_nd[lev]->boxArray(), z_phys_nd[lev]->DistributionMap(),
                                 1, z_phys_nd[lev]->nGrowVect());
            MultiFab::Copy(z_phys_interp, *z_phys_nd[lev], 0, 0, 1, 0);
        }

        Real dz0_max;
        init_terrain_from_wrfinput(lev, geom[lev], z_top, boxes_at_level[lev][0], z_phys_nd[lev].get(),
                                   mf_PH, *mf_PHB, dz0_max, solverChoice.avg_grid_faces_to_nodes);
        z_phys_nd[lev]->FillBoundary(geom[lev].periodicity());

        if (!solverChoice.avg_grid_faces_to_nodes) {
#ifdef AMREX_USE_FLOAT
            const Real tol = Real(1.e-4);
#else
            const Real tol = Real(1.e-8);
#endif
            Real SFact = Real(1.03);
            Real Nz = static_cast<Real>(zlevels_stag[lev].size() - 1);

            if (dz0_max >= z_top / Nz) {
                SFact = one;
                dz0_max = z_top / Nz;
            } else {
                int max_iter = 50;
                int iter = 0;
                Real F = dz0_max * ((std::pow(SFact, Nz) - one) / (SFact - one)) - z_top;
                while (std::fabs(F) > tol && iter < max_iter) {
                    Real dFdSF = dz0_max * (Nz * std::pow(SFact, Nz - one) * (SFact - one)
                                           - std::pow(SFact, Nz) + one) /
                                           std::pow(SFact - one, two);
                    SFact -= F / dFdSF;
                    SFact = std::max(one + tol, SFact);
                    F = dz0_max * ((std::pow(SFact, Nz) - one) / (SFact - one)) - z_top;
                    ++iter;
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::fabs(F) <= tol,
                                                 "Newton iterations to determine the grid stretching factor failed!\n");
            }

            Print() << "Building an ERF grid with dz0: " << dz0_max <<
                " and stretching factor: " << SFact << "\n";
            Real dz = dz0_max;
            zlevels_stag[lev][0] = zero;
            for (int k(1); k < zlevels_stag[lev].size(); ++k) {
                zlevels_stag[lev][k] = zlevels_stag[lev][k-1] + dz;
                dz *= SFact;
            }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::fabs(zlevels_stag[lev].back() - z_top) <= tol,
                "Top of zlevels_stag does not match z_top!\n");

            update_stretched_dz(lev, zlevels_stag, stretched_dz_h, stretched_dz_d);
            make_terrain_fitted_coords(lev, geom[lev], *z_phys_nd[lev], zlevels_stag[lev], phys_bc_type,
                                       fine_terrain,
                                       (fine_terrain == FineTerrain::Transform) ? &z_phys_interp : nullptr);
        }

        // Initialize metric quantities
        make_J(geom[lev], *z_phys_nd[lev], *detJ_cc[lev]);
        make_areas(geom[lev], *z_phys_nd[lev], *ax[lev], *ay[lev], *az[lev]);
        make_zcc(geom[lev], *z_phys_nd[lev], *z_phys_cc[lev]);
    }

    Print() << "Surface-only initialization from wrfinput complete at level " << lev << ".\n";
    Print() << "Atmospheric state will be interpolated from coarse level via FillCoarsePatch.\n";
}

#endif // ERF_USE_NETCDF
