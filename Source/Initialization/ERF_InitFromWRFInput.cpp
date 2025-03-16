/**
 * \file ERF_InitFromWRFInput.cpp
 */

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_ProbCommon.H>
#include <ERF_DataStruct.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

void
read_from_wrfinput (int lev,
                    const Box& domain,
                    const std::string& fname,
                    FArrayBox& NC_fab,
                    const std::string& NC_name,
                    Geometry& geom,
                    int& use_theta_m,
                    int& success);

Real
read_from_wrfbdy (const std::string& nc_bdy_file,
                  const Box& domain,
                  Vector<Vector<FArrayBox>>& bdy_data_xlo,
                  Vector<Vector<FArrayBox>>& bdy_data_xhi,
                  Vector<Vector<FArrayBox>>& bdy_data_ylo,
                  Vector<Vector<FArrayBox>>& bdy_data_yhi,
                  int& width,
                  Real& start_bdy_time);

void
convert_wrfbdy_data (const Box& domain,
                     Vector<Vector<FArrayBox>>& bdy_data,
                     const MultiFab& mf_MUB,
                     const MultiFab& mf_C1H,
                     const MultiFab& mf_C2H,
                     const MultiFab& xvel,
                     const MultiFab& yvel,
                     const MultiFab& cons,
                     const Geometry& geom,
                     const bool& use_moist);

void
compute_terrain_top_and_bottom (Real& terrain_bottom_min,
                                Real& terrain_bottom_max,
                                Real& terrain_top_min,
                                Real& terrain_top_max,
                                const MultiFab& mf_PH,
                                const MultiFab& mf_PHB,
                                const Box& domain);

void
init_terrain_from_wrfinput (int lev,
                            const Real& z_top,
                            const Box& domain,
                            MultiFab* z_phys,
                            const MultiFab& NC_PH_fab,
                            const MultiFab& NC_PHB_fab);

void
init_base_state_from_wrfinput (const Box& domain,
                               Real l_rdOcp,
                               MoistureType moisture_type,
                               const int& n_qstate,
                               MultiFab& cons_fab,
                               MultiFab& p_hse,
                               MultiFab& pi_hse,
                               MultiFab& th_hse,
                               MultiFab& qv_hse,
                               MultiFab& r_hse,
                               const MultiFab& mf_PB,
                               const MultiFab& mf_P,
                               const bool& use_P_eos);

/**
 * ERF function that initializes data from a WRF dataset
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_wrfinput (int lev)
{
    const Box& domain = geom[lev].Domain();

    if (nc_init_file.empty()) {
        amrex::Error("NetCDF initialization file name must be provided via input");
    }

    bool use_moist = (solverChoice.moisture_type != MoistureType::None);

    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<std::string> NC_names;
    NC_names.push_back("ALB");       // 0 DO RHO FIRST
    NC_names.push_back("AL");        // 1 DO RHO FIRST
    NC_names.push_back("U");         // 2
    NC_names.push_back("V");         // 3
    NC_names.push_back("W");         // 4
    NC_names.push_back("THM");       // 5
    NC_names.push_back("PH");        // 6
    NC_names.push_back("PHB");       // 7
    NC_names.push_back("PB");        // 8
    NC_names.push_back("P");         // 9
    NC_names.push_back("MUB");       // 10
    NC_names.push_back("MAPFAC_UY"); // 11
    NC_names.push_back("MAPFAC_VY"); // 12
    NC_names.push_back("MAPFAC_MY"); // 13
    NC_names.push_back("SST");       // 14
    NC_names.push_back("LANDMASK");  // 15
    NC_names.push_back("C1H");       // 16
    NC_names.push_back("C2H");       // 17
    NC_names.push_back("XLAT_V");    // 18
    NC_names.push_back("XLONG_U");   // 19
    if (use_moist) {
        NC_names.push_back("QVAPOR"); // 20
        NC_names.push_back("QCLOUD"); // 21
        NC_names.push_back("QRAIN");  // 22
    }
    int nvar = NC_names.size();
    Vector<Vector<FArrayBox>> NC_fab_var_file;
    NC_fab_var_file.resize(num_boxes_at_level[lev]);
    for (int idx(0); idx < num_boxes_at_level[lev]; ++ idx) { NC_fab_var_file[idx].resize(nvar); }

    auto& lev_new = vars_new[lev];

    // NOTE: Following MFs must have an underlying BA that follows
    //       the shapes in ERF_ReadFromWRFInput.cpp
    //       Most are 3D but MU/MUB are 2D and C1/2H are 1D

    MultiFab mf_PH , mf_PHB;         // For geopotential height
    MultiFab mf_ALB, mf_PB , mf_P  ; // For base state
    MultiFab mf_MUB, mf_C1H, mf_C2H; // For bdy convert

    // Temporary MFs for derived quantities
    auto& ba    = lev_new[Vars::cons].boxArray();
    auto& dm    = lev_new[Vars::cons].DistributionMap();
    IntVect ng  = lev_new[Vars::cons].nGrowVect();
    IntVect ngz = (z_phys_nd[lev]) ? z_phys_nd[lev]->nGrowVect() : IntVect(0); ngz[0] +=1; ngz[1] += 1;
    IntVect ngv = ng; ngv[2] = 0;

    // Build 2D BA
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    // Build 1D BA
    BoxList bl1d = ba.boxList();
    for (auto& b : bl1d) {
        b.setRange(0,0);
        b.setRange(1,0);
    }
    BoxArray ba1d(std::move(bl1d));

    bool compute_terrain_here = true;

    Print() << "Loading initial data from NetCDF file at level " << lev << "\n";
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
        Print() << "Reading from file " << nc_init_file[lev][idx] << "\n";
        for (int ivar = 0; ivar < nvar; ++ ivar) {
            Print() << "Reading variable " << NC_names[ivar] << " ...";

            int success, use_theta_m;
            read_from_wrfinput(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                               NC_fab_var_file[idx][ivar], NC_names[ivar], geom[lev],
                               use_theta_m, success);

            auto var_name = NC_names[ivar];
            auto& var_fab = NC_fab_var_file[idx][ivar];

            // Initialize rho =  1/(ALB + AL)
            if ( var_name == "ALB" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                {
                    FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
                    cons_fab.template copy<RunOn::Device>(var_fab, 0, Rho_comp, 1);
                }

            } if ( var_name == "AL" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                {
                    FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
                    Box vbx = cons_fab.box(); vbx.grow(-ng);
                    cons_fab.template   plus<RunOn::Device>(var_fab, 0, Rho_comp, 1);
                    cons_fab.template invert<RunOn::Device>(1.0, vbx, Rho_comp, 1);
                }
            }

            // Initialize other vars (NOT rho)
            if ( var_name == "U"      ||
                 var_name == "V"      ||
                 var_name == "W"      ||
                 var_name == "THM"    ||
                 var_name == "QVAPOR" ||
                 var_name == "QCLOUD" ||
                 var_name == "QRAIN" ) {

              int n_qstate = micro->Get_Qstate_Size();
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              // INITIAL DATA common for "ideal" as well as "real" simulation
              // Don't tile this since we are operating on full FABs in this routine
              for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
              {
                  // Define fabs for holding the initial data
                  int icomp = 0;
                  bool mult_rho = false;
                  FArrayBox& cons_fab = lev_new[Vars::cons][mfi];
                  FArrayBox* cur_fab;
                  if (var_name == "U") {
                    cur_fab  = &lev_new[Vars::xvel][mfi];
                  } else if (var_name == "V") {
                    cur_fab  = &lev_new[Vars::yvel][mfi];
                  } else if (var_name == "W") {
                    cur_fab  = &lev_new[Vars::zvel][mfi];
                  } else if (var_name == "THM") {
                    const Real theta_ref = 300.0;
                    var_fab.template plus<RunOn::Device>(theta_ref);
                    cur_fab  = &lev_new[Vars::cons][mfi];
                    mult_rho = true;
                    icomp    = RhoTheta_comp;
                  } else if (var_name == "QVAPOR") {
                    cur_fab  = &lev_new[Vars::cons][mfi];
                    mult_rho = true;
                    icomp    = RhoQ1_comp;
                  } else if (var_name == "QCLOUD") {
                    cur_fab  = &lev_new[Vars::cons][mfi];
                    mult_rho = true;
                    icomp    = RhoQ2_comp;
                  } else if (var_name == "QRAIN") {
                    cur_fab  = &lev_new[Vars::cons][mfi];
                    mult_rho = true;
                    icomp    = RhoQ3_comp;
                    if (n_qstate > 3) { icomp = RhoQ4_comp; }
                  }

                  if (success) {
                      cur_fab->template copy<RunOn::Device>(var_fab, 0, icomp, 1);
                      if (mult_rho) { cur_fab->template mult<RunOn::Device>(cons_fab, Rho_comp, icomp, 1); }
                      if (use_theta_m && (var_name == "QVAPOR")) {
                          // Now, we can calculate theta = thm / (1 + R_v/R_d * Qv)
                          var_fab.template mult<RunOn::Device>(R_v/R_d);
                          var_fab.template plus<RunOn::Device>(1.0);
                          var_fab.template invert<RunOn::Device>(1.0);
                          cur_fab->template mult<RunOn::Device>(var_fab, 0, RhoTheta_comp, 1);
                      }
                  } else {
                      if (icomp < cur_fab->nComp()) {
                          amrex::Print() << "Setting " << var_name << " to 0 since we couldn't read it in ... DONE" << std::endl;
                          cur_fab->template setVal<RunOn::Device>(0.0,cur_fab->box(),icomp,1);
                          if (mult_rho) { cur_fab->template mult<RunOn::Device>(cons_fab, Rho_comp, icomp, 1); }
                      } else {
                          amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                      }
                  }
              } // mfi
              var_fab.clear();
          } // valid var (not rho)

          if ( var_name == "PH" ) {
              if (success) {
                  auto& ba_w = lev_new[Vars::zvel].boxArray();
                  mf_PH.define(ba_w, dm, 1, ngz);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                  for ( MFIter mfi(mf_PH, false); mfi.isValid(); ++mfi )
                  {
                    FArrayBox &cur_fab = mf_PH[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                  }
                  var_fab.clear();
              } else {
                  amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                  compute_terrain_here = false;
              }
          } else if ( var_name == "PHB" ) {
              if (success) {
                  auto& ba_w = lev_new[Vars::zvel].boxArray();
                  mf_PHB.define(ba_w, dm, 1, ngz);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                  for ( MFIter mfi(mf_PHB, false); mfi.isValid(); ++mfi )
                  {
                    FArrayBox &cur_fab = mf_PHB[mfi];
                    cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                  }
                  var_fab.clear();
              } else {
                  amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                  compute_terrain_here = false;
              }
          } else if ( var_name == "ALB" ) {
              mf_ALB.define(ba, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_ALB, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_ALB[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "PB" ) {
              mf_PB.define(ba, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_PB, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_PB[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "P" ) {
              mf_P.define(ba, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_P, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_P[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "MUB" ) {
              mf_MUB.define(ba2d, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_MUB, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_MUB[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "C1H" ) {
              mf_C1H.define(ba1d, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_C1H, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_C1H[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "C2H" ) {
              mf_C2H.define(ba1d, dm, 1, ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_C2H, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_C2H[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          }


          int i_lo = geom[lev].Domain().smallEnd(0); int i_hi = geom[lev].Domain().bigEnd(0);
          int j_lo = geom[lev].Domain().smallEnd(1); int j_hi = geom[lev].Domain().bigEnd(1);

          // Initialize Latitude
          if ( var_name == "XLAT_V" ) {
            lat_m[lev] = std::make_unique<MultiFab>(ba2d,dm,1,ngv);
            for ( MFIter mfi(*(lat_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
              Box gtbx = mfi.growntilebox();
              const Array4<      Real>& dst_arr = (lat_m[lev])->array(mfi);
              const Array4<const Real>& src_arr = var_fab.const_array();
              ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
              {
                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                dst_arr(i,j,0) = src_arr(li,lj,0);
              });
            }
          }

          // Initialize Longitude
          if ( var_name == "XLONG_U" ) {
            lon_m[lev] = std::make_unique<MultiFab>(ba2d,dm,1,ngv);
            for ( MFIter mfi(*(lon_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
              Box gtbx = mfi.growntilebox();
              const Array4<      Real>& dst_arr = (lon_m[lev])->array(mfi);
              const Array4<const Real>& src_arr = var_fab.const_array();
              ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
              {
                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                dst_arr(i,j,0) = src_arr(li,lj,0);
              });
            }
          }

          // Initialize Landmask
          if ( var_name == "LANDMASK" ) {
            for ( MFIter mfi(*(lmask_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
              Box gtbx = mfi.growntilebox();
              const Array4<       int>& dst_arr = lmask_lev[lev][0]->array(mfi);
              const Array4<const Real>& src_arr = var_fab.const_array();
              ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
              {
                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                dst_arr(i,j,0) = static_cast<int>(src_arr(li,lj,0));
              });
            }
            (lmask_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
          }

          // Initialize MapFac U
          if ( var_name == "MAPFAC_UY" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(*mapfac_u[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
              // Define fabs for holding the initial data
              FArrayBox &msf_fab = (*mapfac_u[lev])[mfi];
              msf_fab.template copy<RunOn::Device>(var_fab);
            }
          }

          // Initialize MapFac V
          if ( var_name == "MAPFAC_VY" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(*mapfac_v[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
              // Define fabs for holding the initial data
              FArrayBox &msf_fab = (*mapfac_v[lev])[mfi];
              msf_fab.template copy<RunOn::Device>(var_fab);
            }
          }

          // Initialize MapFac M
          if ( var_name == "MAPFAC_MY" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(*mapfac_m[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
              // Define fabs for holding the initial data
              FArrayBox &msf_fab = (*mapfac_m[lev])[mfi];
              msf_fab.template copy<RunOn::Device>(var_fab);
            }
          }

          if (success) {
              var_fab.clear();
              Print() << " DONE\n";
          }
        } // ivar
      Print() << "\n";
    } // idx

    // Convert the velocities using the map factors
    for ( MFIter mfi(lev_new[Vars::xvel], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();
        const Array4<      Real>& dst_arr = lev_new[Vars::xvel].array(mfi);
        const Array4<const Real>& src_arr = mapfac_u[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst_arr(i,j,k) /= src_arr(i,j,0);
        });
    }
    for ( MFIter mfi(lev_new[Vars::yvel], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();
        const Array4<      Real>& dst_arr = lev_new[Vars::yvel].array(mfi);
        const Array4<const Real>& src_arr = mapfac_v[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst_arr(i,j,k) /= src_arr(i,j,0);
        });
    }
    for ( MFIter mfi(lev_new[Vars::xvel], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();
        const Array4<      Real>& dst_arr = lev_new[Vars::zvel].array(mfi);
        const Array4<const Real>& src_arr = mapfac_m[lev]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst_arr(i,j,k) /= src_arr(i,j,0);
        });
    }

    if (compute_terrain_here) {
        Real terrain_bottom_min, terrain_bottom_max;
        Real terrain_top_min, terrain_top_max;
        compute_terrain_top_and_bottom(terrain_bottom_min, terrain_bottom_max,
                                       terrain_top_min, terrain_top_max,
                                       mf_PH, mf_PHB, domain);

        // **************************************************************************
        // Initialize the terrain itself and the metric quantities
        // **************************************************************************
        AMREX_ALWAYS_ASSERT(solverChoice.terrain_type == TerrainType::StaticFittedMesh);

        // FillBoundary to populate the internal ghost cells (for averaging)
         mf_PH.FillBoundary(geom[lev].periodicity());
        mf_PHB.FillBoundary(geom[lev].periodicity());
        Real z_top = 0.5 * (terrain_top_min + terrain_top_max);
        amrex::Print() << "Warning: ProbHi(2) will be ignored; we are setting top of domain to " << z_top << std::endl;
        init_terrain_from_wrfinput(lev, z_top, domain, z_phys_nd[lev].get(), mf_PH, mf_PHB);

        make_J  (geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
        make_areas(geom[lev],*z_phys_nd[lev],*ax[lev],*ay[lev],*az[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
    }

    // **************************************************************************
    // Rebalance the base state if needed
    // **************************************************************************
    if (solverChoice.rebalance_wrfinput) {
        int ncomp = lev_new[Vars::cons].nComp();
        int k_dom_lo = geom[lev].Domain().smallEnd(2);
        int k_dom_hi = geom[lev].Domain().bigEnd(2);
        Real tol = 1.0e-10;
        Real grav = CONST_GRAV;
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            Box bx = mfi.tilebox();
            int klo = bx.smallEnd(2);
            int khi = bx.bigEnd(2);
            AMREX_ALWAYS_ASSERT((klo == k_dom_lo) && (khi == k_dom_hi));
            bx.makeSlab(2,klo);

            const Array4<      Real>& con_arr = lev_new[Vars::cons].array(mfi);
            const Array4<const Real>& z_arr = z_phys_nd[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
            {
                // integrate from surface to domain top
                Real Factor;
                Real dz, F, C;
                Real rho_tot_hi, rho_tot_lo;
                Real z_lo, z_hi;
                Real R_lo, R_hi;
                Real qv_lo, qv_hi;
                Real Th_lo, Th_hi;
                Real P_lo, P_hi;

                // First integrate from sea level to the height at klo
                {
                    // Vertical grid spacing
                    z_lo = 0.0;
                    z_hi = 0.125 * (z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(1,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                   +z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(1,j+1,klo+1) + z_arr(i+1,j+1,klo+1));
                    dz   = z_hi - z_lo;

                    // Establish known constant
                    qv_lo = con_arr(i,j,klo,RhoQ1_comp)    / con_arr(i,j,klo,Rho_comp);
                    Th_lo = con_arr(i,j,klo,RhoTheta_comp) / con_arr(i,j,klo,Rho_comp);
                    P_lo  = p_0;
                    R_lo  = getRhogivenThetaPress(Th_lo, P_lo, R_d/Cp_d, qv_lo);
                    rho_tot_lo = R_lo * (1. + qv_lo);
                    C  = -P_lo + 0.5*rho_tot_lo*grav*dz;

                    // Initial guess and residual
                    qv_hi = con_arr(i,j,klo,RhoQ1_comp)    / con_arr(i,j,klo,Rho_comp);
                    Th_hi = con_arr(i,j,klo,RhoTheta_comp) / con_arr(i,j,klo,Rho_comp);
                    P_hi  = p_0;
                    R_hi  = getRhogivenThetaPress(Th_hi, P_hi, R_d/Cp_d, qv_hi);
                    rho_tot_hi = R_hi * (1. + qv_hi);
                    F = P_hi + 0.5*rho_tot_hi*grav*dz + C;

                    // Do iterations
                    HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                                 grav, C, Th_hi,
                                                 qv_hi, qv_hi,
                                                 P_hi, R_hi, F);

                    // Assign data
                    Factor = R_hi / con_arr(i,j,klo,Rho_comp);
                    con_arr(i,j,klo,Rho_comp) = R_hi;
                    for (int n(1); n<ncomp; ++n) { con_arr(i,j,klo,n) *= Factor; }
                    P_lo = P_hi;
                    z_lo = z_hi;
                }

                for (int k(klo+1); k<=khi; ++k) {
                    // Vertical grid spacing
                  z_hi = 0.125 * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(1,j+1,k  ) + z_arr(i+1,j+1,k  )
                                 +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(1,j+1,k+1) + z_arr(i+1,j+1,k+1));
                  dz   = z_hi - z_lo;

                  // Establish known constant
                  qv_lo = con_arr(i,j,k,RhoQ1_comp)    / con_arr(i,j,k,Rho_comp);
                  Th_lo = con_arr(i,j,k,RhoTheta_comp) / con_arr(i,j,k,Rho_comp);
                  R_lo  = getRhogivenThetaPress(Th_lo, P_lo, R_d/Cp_d, qv_lo);
                  rho_tot_lo = R_lo * (1. + qv_lo);
                  C  = -P_lo + 0.5*rho_tot_lo*grav*dz;

                  // Initial guess and residual
                  qv_hi = con_arr(i,j,k,RhoQ1_comp)    / con_arr(i,j,k,Rho_comp);
                  Th_hi = con_arr(i,j,k,RhoTheta_comp) / con_arr(i,j,k,Rho_comp);
                  R_hi  = getRhogivenThetaPress(Th_hi, P_hi, R_d/Cp_d, qv_hi);
                  rho_tot_hi = R_hi * (1. + qv_hi);
                  F = P_hi + 0.5*rho_tot_hi*grav*dz + C;

                  // Do iterations
                  HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                               grav, C, Th_hi,
                                               qv_hi, qv_hi,
                                               P_hi, R_hi, F);

                  // Assign data
                  Factor = R_hi / con_arr(i,j,k,Rho_comp);
                  con_arr(i,j,k,Rho_comp) = R_hi;
                  for (int n(1); n<ncomp; ++n) { con_arr(i,j,k,n) *= Factor; }
                  P_lo = P_hi;
                  z_lo = z_hi;
                }
            });
        } // mfi
    } // rebalance_wrfinput

    // **************************************************************************
    // Initialize the base state
    // **************************************************************************
    const Real l_rdOcp = solverChoice.rdOcp;
    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    if (solverChoice.init_type == InitType::WRFInput) {

        int n_qstate = micro->Get_Qstate_Size();

        bool use_P_eos = (solverChoice.rebalance_wrfinput);

        init_base_state_from_wrfinput(domain, l_rdOcp, solverChoice.moisture_type, n_qstate,
                                      lev_new[Vars::cons], p_hse, pi_hse, th_hse, qv_hse, r_hse,
                                      mf_PB, mf_P, use_P_eos);

        // FillBoundary to populate the internal ghost cells (no averaging in above call)
         r_hse.FillBoundary(geom[lev].periodicity());
         p_hse.FillBoundary(geom[lev].periodicity());
        pi_hse.FillBoundary(geom[lev].periodicity());
        th_hse.FillBoundary(geom[lev].periodicity());
        qv_hse.FillBoundary(geom[lev].periodicity());
    }

    // Initialize the bdy data
    if (solverChoice.init_type == InitType::WRFInput && solverChoice.use_real_bcs && (lev == 0))
    {
        if (nc_bdy_file.empty()) {
            amrex::Error("NetCDF boundary file name must be provided via input");
        }

        // Three points are necessary if a relaxation zone is present.
        if (real_width > real_set_width) {
            AMREX_ALWAYS_ASSERT(real_width-real_set_width >= 3);
        }

        bdy_time_interval = read_from_wrfbdy(nc_bdy_file,geom[0].Domain(),
                                             bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                             real_width, start_bdy_time);

        Print() << "Read in boundary data with width "  << real_width << std::endl;
        Print() << "Running with specification width: " << real_set_width
                << " and relaxation width: " << real_width - real_set_width << std::endl;

        convert_wrfbdy_data(domain,bdy_data_xlo,
                            mf_MUB, mf_C1H, mf_C2H,
                            lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                            geom[lev], use_moist);
        convert_wrfbdy_data(domain,bdy_data_xhi,
                            mf_MUB, mf_C1H, mf_C2H,
                            lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                            geom[lev], use_moist);
        convert_wrfbdy_data(domain,bdy_data_ylo,
                            mf_MUB, mf_C1H, mf_C2H,
                            lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                            geom[lev], use_moist);
        convert_wrfbdy_data(domain,bdy_data_yhi,
                            mf_MUB, mf_C1H, mf_C2H,
                            lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                            geom[lev], use_moist);
    } // init_type == Real && lev == 0

    if (solverChoice.init_type == InitType::WRFInput && solverChoice.use_real_bcs)
    {
        //
        // Start at the earliest time (read_from_wrfbdy)
        // Note we only have start_bdy_time if at level 0 and init_type == InitType:WRFInput
        //
        if (lev == 0) {
            Print() << "Setting start_time to "
                    << std::setprecision(timeprecision) << start_bdy_time
                    << " from wrfbdy" << std::endl;
            t_new[lev] = start_bdy_time;
            t_old[lev] = start_bdy_time - 1.e200;
        } else {
            t_new[lev] = t_new[0];
            t_old[lev] = t_old[0];
        }
    }
}


/**
 * Helper function to initialize hydrostatic base state data from WRF dataset
 *
 * @param lev Integer specifying current level
 * @param valid_bx Box specifying the index space we are to initialize
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param p_hse FArrayBox specifying the hydrostatic base state pressure we initialize
 * @param pi_hse FArrayBox specifying the hydrostatic base state Exner pressure we initialize
 * @param th_hse FArrayBox specifying the hydrostatic base state potential temperature
 * @param qv_hse FArrayBox specifying the hydrostatic base state qv
 * @param r_hse FArrayBox specifying the hydrostatic base state density we initialize
 * @param NC_ALB_fab Vector of FArrayBox objects containing WRF data specifying 1/density
 * @param NC_PB_fab Vector of FArrayBox objects containing WRF data specifying pressure
 */
void
init_base_state_from_wrfinput (const Box& domain,
                               const Real l_rdOcp,
                               MoistureType moisture_type,
                               const int& n_qstate,
                               MultiFab& cons,
                               MultiFab& p_hse,
                               MultiFab& pi_hse,
                               MultiFab& th_hse,
                               MultiFab& qv_hse,
                               MultiFab& r_hse,
                               const MultiFab& mf_PB,
                               const MultiFab& mf_P,
                               const bool& use_P_eos)
{
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    for ( MFIter mfi(p_hse, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        Box tbx  = mfi.tilebox();
        Box gtbx = mfi.growntilebox();

        const Array4<Real      >&   cons_arr = cons.array(mfi);
        const Array4<Real      >&  p_hse_arr = p_hse.array(mfi);
        const Array4<Real      >& pi_hse_arr = pi_hse.array(mfi);
        const Array4<Real      >& th_hse_arr = th_hse.array(mfi);
        const Array4<Real      >& qv_hse_arr = qv_hse.array(mfi);
        const Array4<Real      >&  r_hse_arr = r_hse.array(mfi);
        const Array4<Real const>&  nc_pb_arr = mf_PB.const_array(mfi);
        const Array4<Real const>&   nc_p_arr = mf_P.const_array(mfi);

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Base state needs ghost cells filled, protect FAB access
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            int kk = std::max(k , dom_lo.z);
                kk = std::min(kk, dom_hi.z);

            // Base plus perturbational pressure
            Real Ptot = nc_pb_arr(ii,jj,kk) + nc_p_arr(ii,jj,kk);

            // Compute pressure from EOS
            Real Qv    = (moisture_type != MoistureType::None) ?
                         cons_arr(ii,jj,kk,RhoQ1_comp) / cons_arr(ii,jj,kk,Rho_comp) : 0.0;
            Real RT    = cons_arr(ii,jj,kk,RhoTheta_comp);
            Real P_eos = getPgivenRTh(RT, Qv);
            if (use_P_eos) { Ptot = P_eos; }
            Real DelP  = std::fabs(Ptot - P_eos);

            // NOTE: Ghost cells don't contain valid data
            //       We want domain GCs and FB picks up interior GCs
            if (tbx.contains(i,j,k)) {
                if ( (DelP > 1.0) || (DelP/Ptot > 1e-6) ) {
                    printf("p (%i, %i, %i): %e; p_eos: %e; (qv = %e, rho = %e, rT = %e) \n",
                           i, j, k, Ptot, P_eos, Qv, cons_arr(ii,jj,kk,Rho_comp), RT);
                    amrex::Abort("Initial state is inconsistent with EOS!?");
                }
            }

            // Compute rhse
            Real Rhse_Sum = cons_arr(ii,jj,kk,Rho_comp);
            for (int q_offset(0); q_offset<n_qstate; ++q_offset) {
                Rhse_Sum += cons_arr(ii,jj,kk,RhoQ1_comp+q_offset);
            }

            r_hse_arr(i,j,k)  = Rhse_Sum;
            p_hse_arr(i,j,k)  = Ptot;
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
            th_hse_arr(i,j,k) = getRhoThetagivenP(p_hse_arr(i,j,k), Qv) / cons_arr(ii,jj,kk,Rho_comp);
            qv_hse_arr(i,j,k) = Qv;
        });
    }
}

/**
 * Helper function for verifying the top boundary is valid and computing the bottom boundary.
 *
 * @param z_top      Real imposed top boundary
 * @param NC_PH_fab  Vector of FArrayBox objects storing WRF terrain coordinate data (PH)
 * @param NC_PHB_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PHB)
 */
void
compute_terrain_top_and_bottom (Real& terrain_bottom_min,
                                Real& terrain_bottom_max,
                                Real& terrain_top_min,
                                Real& terrain_top_max,
                                const MultiFab& mf_PH,
                                const MultiFab& mf_PHB,
                                const Box& domain)
{
    //
    // For the bottom/top boundary (in that order)
    //
    Gpu::HostVector  <Real> Max_h(3,-1.0e16);
    Gpu::DeviceVector<Real> Max_d(3);
    Gpu::copy(Gpu::hostToDevice, Max_h.begin(), Max_h.end(), Max_d.begin());

    Gpu::HostVector  <Real> Min_h(1, 1.0e16);
    Gpu::DeviceVector<Real> Min_d(1);
    Gpu::copy(Gpu::hostToDevice, Min_h.begin(), Min_h.end(), Min_d.begin());

    Real* min_d = Min_d.data();
    Real* max_d = Max_d.data();

    //
    // ********************************************************************************
    //

    // Index type of (0,0,1)
    int klo = domain.smallEnd()[2];
    int khi = domain.bigEnd()[2]+1;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_PH, TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        Box vbx = mfi.validbox();

        Box nodal_box = amrex::surroundingNodes(vbx);
        int ilo = nodal_box.smallEnd()[0];
        int ihi = nodal_box.bigEnd()[0];
        int jlo = nodal_box.smallEnd()[1];
        int jhi = nodal_box.bigEnd()[1];

        // For the top boundary
        Box Fab2dBox_hi, Fab2dBox_hi_m1;
        if (vbx.bigEnd(2) == khi) {
            Fab2dBox_hi    = makeSlab(vbx,2,khi  );
            Fab2dBox_hi_m1 = makeSlab(vbx,2,khi-1);
        }

        // For the bottom boundary
        Box Fab2dBox_lo;
        if (vbx.smallEnd(2) == klo) {
            Fab2dBox_lo = makeSlab(vbx,2,klo);
        }

        auto const& phb = mf_PHB.const_array(mfi);
        auto const& ph  = mf_PH.const_array(mfi);

        //
        // This loop computes the min and max values of the bottom surface
        //
        ParallelFor(Fab2dBox_lo, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            Real z_calc_lo = 0.25 * ( ph (ii,jj  ,klo) + ph (ii-1,jj  ,klo) +
                                      ph (ii,jj-1,klo) + ph (ii-1,jj-1,klo) +
                                      phb(ii,jj  ,klo) + phb(ii-1,jj  ,klo) +
                                      phb(ii,jj-1,klo) + phb(ii-1,jj-1,klo) ) / CONST_GRAV;
            amrex::Gpu::Atomic::Min(&(min_d[0]),z_calc_lo);
            amrex::Gpu::Atomic::Max(&(max_d[0]),z_calc_lo);
        });

        //
        // This loop computes the max value of the top surface
        //
        ParallelFor(Fab2dBox_hi, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            Real z_calc_hi = 0.25 * ( ph (ii,jj  ,khi) + ph (ii-1,jj  ,khi) +
                                      ph (ii,jj-1,khi) + ph (ii-1,jj-1,khi) +
                                      phb(ii,jj  ,khi) + phb(ii-1,jj  ,khi) +
                                      phb(ii,jj-1,khi) + phb(ii-1,jj-1,khi) ) / CONST_GRAV;
            amrex::Gpu::Atomic::Max(&(max_d[1]),z_calc_hi);
        });

        //
        // This loop computes the max value of the layer just below the top surface
        //
        ParallelFor(Fab2dBox_hi_m1, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            Real z_calc_hi = 0.25 * ( ph (ii,jj  ,khi-1) + ph (ii-1,jj  ,khi-1) +
                                      ph (ii,jj-1,khi-1) + ph (ii-1,jj-1,khi-1) +
                                      phb(ii,jj  ,khi-1) + phb(ii-1,jj  ,khi-1) +
                                      phb(ii,jj-1,khi-1) + phb(ii-1,jj-1,khi-1) ) / CONST_GRAV;
            amrex::Gpu::Atomic::Max(&(max_d[2]),z_calc_hi);
        });
    } // mfi

    Gpu::copy(Gpu::deviceToHost, Min_d.begin(), Min_d.end(), Min_h.begin());
    Gpu::copy(Gpu::deviceToHost, Max_d.begin(), Max_d.end(), Max_h.begin());

    ParallelDescriptor::ReduceRealMin(Min_h[0]);

    ParallelDescriptor::ReduceRealMax(Max_h[0]);
    ParallelDescriptor::ReduceRealMax(Max_h[1]);
    ParallelDescriptor::ReduceRealMax(Max_h[2]);

    terrain_bottom_max = Max_h[0];
    terrain_bottom_min = Min_h[0];
    terrain_top_max    = Max_h[1];
    terrain_top_min    = Max_h[2];

    Print() << "Terrain     has min value = " << terrain_bottom_min << " and max value = " << terrain_bottom_max << std::endl;
    Print() << "Top of mesh has min value = " << terrain_top_min    << " and max value = " << terrain_top_max << std::endl;
}

/**
 * Helper function for initializing terrain coordinates from a WRF dataset.
 *
 * @param lev Integer specifying the current level
 * @param z_phys FArrayBox specifying the node-centered z coordinates of the terrain
 * @param NC_PH_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PH)
 * @param NC_PHB_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PHB)
 */
void
init_terrain_from_wrfinput (int /*lev*/,
                            const Real& z_top,
                            const Box& domain,
                            MultiFab* z_phys,
                            const MultiFab& mf_PH,
                            const MultiFab& mf_PHB)
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*z_phys, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box gnbx = mfi.growntilebox();

        // This copies from NC_zphys on z-faces to z_phys_nd on nodes
        const Array4<Real      >&      z_arr = z_phys->array(mfi);
        const Array4<Real const>& nc_phb_arr = mf_PHB.const_array(mfi);
        const Array4<Real const>& nc_ph_arr  = mf_PH.const_array(mfi);

        // PHB and PH are on z-faces (0.5 dx/y ahead of zphys)
        Box z_face_box = convert(domain,IntVect(0,0,1));

        // Prevent averaging from going into domain ghost cells
        int ilo = z_face_box.smallEnd()[0] + 1;
        int ihi = z_face_box.bigEnd()[0];
        int jlo = z_face_box.smallEnd()[1] + 1;
        int jhi = z_face_box.bigEnd()[1];

        int klo = domain.smallEnd()[2];
        int khi = domain.bigEnd()[2]+1;

        ParallelFor(gnbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int ii = std::max(std::min(i,ihi),ilo);
            int jj = std::max(std::min(j,jhi),jlo);
            if (k < klo) {
                Real z_klo   =  0.25 * ( nc_ph_arr (ii,jj  ,klo  ) + nc_ph_arr (ii-1,jj  ,klo  ) +
                                         nc_ph_arr (ii,jj-1,klo  ) + nc_ph_arr (ii-1,jj-1,klo) +
                                         nc_phb_arr(ii,jj  ,klo  ) + nc_phb_arr(ii-1,jj  ,klo  ) +
                                         nc_phb_arr(ii,jj-1,klo  ) + nc_phb_arr(ii-1,jj-1,klo) ) / CONST_GRAV;
                Real z_klop1 =  0.25 * ( nc_ph_arr (ii,jj  ,klo+1) + nc_ph_arr (ii-1,jj  ,klo+1) +
                                         nc_ph_arr (ii,jj-1,klo+1) + nc_ph_arr (ii-1,jj-1,klo+1) +
                                         nc_phb_arr(ii,jj  ,klo+1) + nc_phb_arr(ii-1,jj  ,klo+1) +
                                         nc_phb_arr(ii,jj-1,klo+1) + nc_phb_arr(ii-1,jj-1,klo+1) ) / CONST_GRAV;
                z_arr(i, j, k) = 2.0 * z_klo - z_klop1;
            } else if (k > khi) {
                Real z_khim1 =  0.25 * ( nc_ph_arr (ii,jj  ,khi-1) + nc_ph_arr (ii-1,jj  ,khi-1) +
                                         nc_ph_arr (ii,jj-1,khi-1) + nc_ph_arr (ii-1,jj-1,khi-1) +
                                         nc_phb_arr(ii,jj  ,khi-1) + nc_phb_arr(ii-1,jj  ,khi-1) +
                                         nc_phb_arr(ii,jj-1,khi-1) + nc_phb_arr(ii-1,jj-1,khi-1) ) / CONST_GRAV;
                z_arr(i, j, k) = 2.0 * z_top - z_khim1;
            } else if (k == khi) {
                z_arr(i, j, k) = z_top;
            } else {
                z_arr(i, j, k) = 0.25 * ( nc_ph_arr (ii,jj  ,k) + nc_ph_arr (ii-1,jj  ,k) +
                                          nc_ph_arr (ii,jj-1,k) + nc_ph_arr (ii-1,jj-1,k) +
                                          nc_phb_arr(ii,jj  ,k) + nc_phb_arr(ii-1,jj  ,k) +
                                          nc_phb_arr(ii,jj-1,k) + nc_phb_arr(ii-1,jj-1,k) ) / CONST_GRAV;
            } // k
        });
    } // mfi
}
#endif // ERF_USE_NETCDF
