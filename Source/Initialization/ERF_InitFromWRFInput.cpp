/**
 * \file ERF_InitFromWRFInput.cpp
 */

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_ProbCommon.H>
#include <ERF_DataStruct.H>

#include <ERF_ReadFromWRFInput.H>
#include <ERF_ReadFromWRFBdy.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

#include "ERF_NCWpsFile.H"

Real
compute_terrain_top_and_bottom (const MultiFab& mf_PH,
                                const MultiFab& mf_PHB,
                                const Box& domain);

void
init_terrain_from_wrfinput (int lev,
                            const Real& z_top,
                            const Box& subdomain,
                            MultiFab* z_phys,
                            const MultiFab& NC_PH_fab,
                            const MultiFab& NC_PHB_fab);

void
init_base_state_from_wrfinput (const Box& subdomain,
                               Real l_rdOcp,
                               MoistureType moisture_type,
                               const int& n_qstate_moist,
                               MultiFab& cons_fab,
                               MultiFab& p_hse,
                               MultiFab& pi_hse,
                               MultiFab& th_hse,
                               MultiFab& qv_hse,
                               MultiFab& r_hse,
                               const MultiFab& mf_PB,
                               const MultiFab& mf_P,
                               const bool& use_P_eos);

Real
read_start_time_from_wrfinput(int lev, const std::string& fname)
{
    std::string NC_dateTime;
    Real        NC_epochTime;
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);

        NC_dateTime = ncf.get_attr("SIMULATION_START_DATE");

        const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";
        NC_epochTime = getEpochTime(NC_dateTime, dateTimeFormat);

        ncf.close();

        amrex::Print() << "Have read start_time string at level "<< lev << " is " << NC_dateTime << std::endl;
        amrex::Print() << "Have read start_time number at level "<< lev << " is " << NC_epochTime << std::endl;
    }

    amrex::ParallelDescriptor::Bcast(&NC_epochTime,1,amrex::ParallelDescriptor::IOProcessorNumber());

    return NC_epochTime;
}

/**
 * ERF function that initializes data from a WRF dataset
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_wrfinput (int lev,
                         MultiFab& mf_C1H_lev,
                         MultiFab& mf_C2H_lev,
                         MultiFab& mf_MUB_lev,
                         MultiFab& mf_PSFC_lev)
{
    if (nc_init_file.empty()) {
        amrex::Error("NetCDF initialization file name must be provided via input");
    }

    bool use_moist = (solverChoice.moisture_type != MoistureType::None);
    bool use_lsm = (solverChoice.lsm_type != LandSurfaceType::None);

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
    NC_names.push_back("PSFC");      // 10
    NC_names.push_back("MUB");       // 11
    NC_names.push_back("MAPFAC_U");  // 12
    NC_names.push_back("MAPFAC_V");  // 13
    NC_names.push_back("MAPFAC_M");  // 14
    NC_names.push_back("SST");       // 15
    NC_names.push_back("TSK");       // 16
    NC_names.push_back("LANDMASK");  // 17
    NC_names.push_back("C1H");       // 18
    NC_names.push_back("C2H");       // 19
    NC_names.push_back("XLAT_V");    // 20
    NC_names.push_back("XLONG_U");   // 21
    if (use_moist) {
        NC_names.push_back("QVAPOR"); // 22
        NC_names.push_back("QCLOUD"); // 23
        NC_names.push_back("QRAIN");  // 24
    }
    NC_names.push_back("IVGTYP");     // 25
    NC_names.push_back("ISLTYP");     // 26
    if (use_lsm) {
        NC_names.push_back("TSLB");   // 27
        NC_names.push_back("SMOIS");  // 28
        NC_names.push_back("SH2O");   // 29
        NC_names.push_back("LAI");    // 30
        NC_names.push_back("ZS");     // 31
        NC_names.push_back("DZS");    // 32
        NC_names.push_back("VEGFRA"); // 33
        NC_names.push_back("TMN");    // 34
        NC_names.push_back("SHDMIN"); // 35
        NC_names.push_back("SHDMAX"); // 36

        // --- debugging ---
        // print LSM varname->WRF input name map
        auto &lsm_wrfmap = lsm.Get_WRFInputNames();
        for (const auto &[wrfname, lsmname] : lsm_wrfmap) {
            amrex::Print() << " LSM input for WRF name '" << wrfname << "' -> '" << lsmname << "'" << std::endl;
        }
        // ---

        if (lsm_wrfmap.size() == 0) {
            amrex::Print() << "Warning: LSM model is being used, but no mapping is defined to fill its variables from WRFinput!" << std::endl;
        }
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

    // Temporary MFs for derived quantities
    auto& ba    = lev_new[Vars::cons].boxArray();
    auto& dm    = lev_new[Vars::cons].DistributionMap();
    IntVect ng  = lev_new[Vars::cons].nGrowVect();
    IntVect ngz = (z_phys_nd[lev]) ? z_phys_nd[lev]->nGrowVect() : IntVect(0); ngz[0] +=1; ngz[1] += 1;
    IntVect ngv = ng; ngv[2] = 0;

    bool compute_terrain_here = true;

    const Real l_rdOcp = solverChoice.rdOcp;

    Print() << "Loading initial data from NetCDF file at level " << lev << "\n";
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++) {
        Print() << "Reading from file " << nc_init_file[lev][idx] << "\n";
        for (int ivar = 0; ivar < nvar; ++ ivar) {
            Print() << "Checking for " << NC_names[ivar] << " ...";

            int success, use_theta_m;
            read_from_wrfinput(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                               NC_fab_var_file[idx][ivar], NC_names[ivar], geom[lev],
                               use_theta_m, success);

            auto var_name = NC_names[ivar];
            auto& var_fab = NC_fab_var_file[idx][ivar];

            if (lev > 1) {
                Box shift_by_box(subdomains[lev][0].minimalBox());
                IntVect shift_by(shift_by_box.smallEnd());
                for (int i = 0; i < AMREX_SPACEDIM; i++) {
                    shift_by[i] -= var_fab.box().smallEnd(i);
                }
                var_fab.shift(shift_by);
            }

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

            if ( var_name == "THM") {
                const Real wrf_theta_ref = 300.0;
                var_fab.template plus<RunOn::Device>(wrf_theta_ref);
            }

            // Initialize velocities
            if ( var_name == "U"      ||
                 var_name == "V"      ||
                 var_name == "W")
            {
                for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                {
                    FArrayBox* cur_fab;
                    if (var_name == "U") {
                      cur_fab  = &lev_new[Vars::xvel][mfi];
                    } else if (var_name == "V") {
                      cur_fab  = &lev_new[Vars::yvel][mfi];
                    } else if (var_name == "W") {
                      cur_fab  = &lev_new[Vars::zvel][mfi];
                    }

                    if (success) {
                        cur_fab->template copy<RunOn::Device>(var_fab, 0, 0, 1);
                    } else {
                        amrex::Print() << "Setting " << var_name << " to 0 since we couldn't read it in ... DONE" << std::endl;
                        cur_fab->template setVal<RunOn::Device>(0.0,cur_fab->box(),0,1);
                    }
                } // mfi
            }

            // Initialize cell-centered variables that need to be density-weighted
            if ( var_name == "THM"    ||
                 var_name == "QVAPOR" ||
                 var_name == "QCLOUD" ||
                 var_name == "QRAIN" )
            {
                int n_qstate_moist = micro->Get_Qstate_Moist_Size();
                AMREX_ALWAYS_ASSERT(micro->Get_Qstate_NonMoist_Size() == 0);

                int icomp = -1;
                if (var_name == "THM") {
                    icomp    = RhoTheta_comp;
                } else if (var_name == "QVAPOR") {
                    icomp    = RhoQ1_comp;
                } else if (var_name == "QCLOUD") {
                    icomp    = RhoQ2_comp;
                } else if (var_name == "QRAIN") {
                    icomp    = RhoQ3_comp;
                    if (n_qstate_moist > 3) { icomp = RhoQ4_comp; }
                    if (n_qstate_moist < 3) { success = 0; }
                }

                // INITIAL DATA common for "ideal" as well as "real" simulation
                // Don't tile this since we are operating on full FABs in this routine
                if (success)
                {
                    for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                    {
                        lev_new[Vars::cons][mfi].template copy<RunOn::Device>(var_fab, 0, icomp, 1);
                    }

                    // Multiply by density
                    MultiFab::Multiply(lev_new[Vars::cons], lev_new[Vars::cons], Rho_comp, icomp, 1, lev_new[Vars::cons].nGrowVect());

                    if (use_theta_m && (var_name == "QVAPOR")) {
                        // Now, we can calculate theta = thm / (1 + R_v/R_d * Qv)
                        var_fab.template mult<RunOn::Device>(R_v/R_d);
                        var_fab.template plus<RunOn::Device>(1.0);
                        var_fab.template invert<RunOn::Device>(1.0);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                        for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                        {
                            lev_new[Vars::cons][mfi].template mult<RunOn::Device>(var_fab, 0, RhoTheta_comp, 1);
                        }
                    } // use_theta_m

                } else {
                    if (icomp < lev_new[Vars::cons].nComp()) {
                        amrex::Print() << "Setting " << var_name << " to 0 since we couldn't read it in ... DONE" << std::endl;
                        lev_new[Vars::cons].setVal(0.0,icomp,1);
                    } else {
                        amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                    }
                }

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
          } else if ( var_name == "PSFC" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_PSFC_lev, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_PSFC_lev[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "MUB" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_MUB_lev, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_MUB_lev[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "C1H" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_C1H_lev, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_C1H_lev[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "C2H" ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(mf_C2H_lev, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = mf_C2H_lev[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          }

          bool lat_periodic = (geom[lev].isPeriodic(0) && geom[lev].isPeriodic(1));
          int i_lo = boxes_at_level[lev][0].smallEnd(0); int i_hi = boxes_at_level[lev][0].bigEnd(0);
          int j_lo = boxes_at_level[lev][0].smallEnd(1); int j_hi = boxes_at_level[lev][0].bigEnd(1);

          // Initialize Latitude & Coriolis factors
          if ( var_name == "XLAT_V" ) {
              solverChoice.has_lat_lon = true;
              lat_m[lev]    = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
              sinPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
              cosPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
              for ( MFIter mfi(*(lat_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                  Box gtbx = mfi.growntilebox();
                  const Array4<      Real>& sin_arr = (sinPhi_m[lev])->array(mfi);
                  const Array4<      Real>& cos_arr = (cosPhi_m[lev])->array(mfi);
                  const Array4<      Real>& dst_arr = (lat_m[lev])->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, i_lo), i_hi);
                      int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                      dst_arr(i,j,0) = src_arr(li,lj,0);

                      Real lat_rad = dst_arr(i,j,0) * (PI/180.);
                      sin_arr(i,j,0) = std::sin(lat_rad);
                      cos_arr(i,j,0) = std::cos(lat_rad);
                  });
              }
          }

          // Initialize Longitude
          if ( var_name == "XLONG_U" ) {
              lon_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
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

          // Initialize SST
          if ( var_name == "SST" ) {
              sst_lev[lev][0] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
              for ( MFIter mfi(*(sst_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                  Box gtbx = mfi.growntilebox();
                  const Array4<      Real>& dst_arr = sst_lev[lev][0]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  const Array4<const Real>& psfc_arr = mf_PSFC_lev.const_array(mfi);
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, i_lo), i_hi);
                      int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                      dst_arr(i,j,0) = getThgivenTandP(src_arr(li,lj,0), psfc_arr(li,lj,0), l_rdOcp);
                  });
              }
              (sst_lev[lev][0])->FillBoundary(geom[lev].periodicity());
          }

          // Initialize TSK
          if ( var_name == "TSK" ) {
              tsk_lev[lev][0] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
              for ( MFIter mfi(*(tsk_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                  Box gtbx = mfi.growntilebox();
                  const Array4<      Real>& dst_arr = tsk_lev[lev][0]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  const Array4<const Real>& psfc_arr = mf_PSFC_lev.const_array(mfi);
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, i_lo), i_hi);
                      int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                      dst_arr(i,j,0) = getThgivenTandP(src_arr(li,lj,0), psfc_arr(li,lj,0), l_rdOcp);
                  });
              }
              (tsk_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
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

          // Initialize Landtype
          if ( var_name == "IVGTYP" ) {
              for ( MFIter mfi(*(land_type_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                  Box gtbx = mfi.growntilebox();
                  const Array4<       int>& dst_arr = land_type_lev[lev][0]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, i_lo), i_hi);
                      int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                      dst_arr(i,j,0) = static_cast<int>(src_arr(li,lj,0));
                  });
              }
              (land_type_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
          }

          // Initialize Soil type
          if ( var_name == "ISLTYP" ) {
              for ( MFIter mfi(*(soil_type_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                  Box gtbx = mfi.growntilebox();
                  const Array4<       int>& dst_arr = soil_type_lev[lev][0]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, i_lo), i_hi);
                      int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                      dst_arr(i,j,0) = static_cast<int>(src_arr(li,lj,0));
                  });
              }
              (soil_type_lev[lev])[0]->FillBoundary(geom[lev].periodicity());
          }

          // Initialize any LSM variables
          if (use_lsm) {
              auto &lsm_wrfmap = lsm.Get_WRFInputNames();
              for (auto &var : lsm_wrfmap) {
                  if (var_name == var.first) {
                      bool is_3d = var_fab.box().length(2) > 1;
                      amrex::Print() << "   Reading " << ((is_3d) ? "3D" : "2D") << " LSM variable '" << var.first << "' (" << var.second << ")" << std::endl;
                      int lsm_idx = lsm.Get_DataIdx(lev, var.second);
                      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lsm_idx != -1, "LSM variable mapping invalid!");
                      AMREX_ALWAYS_ASSERT(lsm_data[lev][lsm_idx]);

                      int lsm_nsoil = lsm.Get_Lsm_Geom(lev).Domain().length(2);
                      amrex::Print() << " LSM NZ = " << lsm_nsoil << " WRFINPUT NZ = " << var_fab.box().length(2) << std::endl;
                      if (is_3d) {
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lsm_nsoil == var_fab.box().length(2), "Number of soil layers must match!");
                      }

                      // check for special case of single column data (such as soil thickness ZS, DZS)
                      //  the single column is duplicated across all grid points
                      bool is_column = var_fab.box().length(0) == 1 && var_fab.box().length(1) == 1;

                      for ( MFIter mfi(*lsm_data[lev][lsm_idx], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                          Box gtbx = mfi.tilebox();
                          int lsm_khi = gtbx.bigEnd(2);
                          gtbx.setRange(2, 0, var_fab.box().length(2));
                          const Array4<      Real>& dst_arr = lsm_data[lev][lsm_idx]->array(mfi);
                          const Array4<const Real>& src_arr = var_fab.const_array();
                          ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                          {
                              int li = amrex::min(amrex::max(i, i_lo), i_hi);
                              int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                              if (is_column) {
                                // single column src input is copied to all i,j in dst
                                li = 0;
                                lj = 0;
                              }
                              // Note: LSM z levels are at negative k below surface
                              //  map [0, nsoil-1] to [-1, -nsoil]
                              const int lsm_k = lsm_khi - k;
                              dst_arr(i,j,lsm_k) = src_arr(li,lj,k);
                          });
                      }
                      (lsm_data[lev][lsm_idx])->FillBoundary(geom[lev].periodicity());
                  }
              }
          }

          // Initialize MapFac U
          if ( var_name == "MAPFAC_U" ) {
              Real max_val = var_fab.template max<RunOn::Device>();
              if (std::fabs(max_val) < std::numeric_limits<Real>::epsilon()) {
                  Print() << "MAPFAC_U cannot be 0, resetting to 1!\n";
                  var_fab.template setVal<RunOn::Device>(1.0);
              }
              if (lat_periodic) {
                  Print() << "MAPFAC_U resetting to 1 with lateral periodic BCs!\n";
                  var_fab.template setVal<RunOn::Device>(1.0);
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mapfac[lev][MapFacType::u_x], TilingIfNotGPU()); mfi.isValid(); ++mfi )
              {
                  Box gtbx = mfi.growntilebox();
                  Box vbx  = mfi.validbox();
                  int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
                  int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
                  const Array4<      Real>& dst_arr = mapfac[lev][MapFacType::u_x]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, ilo), ihi);
                      int lj = amrex::min(amrex::max(j, jlo), jhi);
                      dst_arr(i,j,0) = src_arr(li,lj,0);
                  });
              }
              mapfac[lev][MapFacType::u_x]->FillBoundary(geom[lev].periodicity());
          }

          // Initialize MapFac V
          if ( var_name == "MAPFAC_V" ) {
              Real max_val = var_fab.template max<RunOn::Device>();
              if (std::fabs(max_val) < std::numeric_limits<Real>::epsilon()) {
                  Print() << "MAPFAC_V cannot be 0, resetting to 1!\n";
                  var_fab.template setVal<RunOn::Device>(1.0);
              }
              if (lat_periodic) {
                  Print() << "MAPFAC_V resetting to 1 with lateral periodic BCs!\n";
                  var_fab.template setVal<RunOn::Device>(1.0);
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mapfac[lev][MapFacType::v_x], TilingIfNotGPU()); mfi.isValid(); ++mfi )
              {
                  Box gtbx = mfi.growntilebox();
                  Box vbx  = mfi.validbox();
                  int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
                  int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
                  const Array4<      Real>& dst_arr = mapfac[lev][MapFacType::v_x]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, ilo), ihi);
                      int lj = amrex::min(amrex::max(j, jlo), jhi);
                      dst_arr(i,j,0) = src_arr(li,lj,0);
                  });
              }
              mapfac[lev][MapFacType::v_x]->FillBoundary(geom[lev].periodicity());
          }

          // Initialize MapFac M
          if ( var_name == "MAPFAC_M" ) {
              Real max_val = var_fab.template max<RunOn::Device>();
              if (std::fabs(max_val) < std::numeric_limits<Real>::epsilon()) {
                  Print() << "MAPFAC_M cannot be 0, resetting to 1!\n";
                  var_fab.template setVal<RunOn::Device>(1.0);
              }
              if (lat_periodic) {
                  Print() << "MAPFAC_M resetting to 1 with lateral periodic BCs!\n";
                  var_fab.template setVal<RunOn::Device>(1.0);
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mapfac[lev][MapFacType::m_x], TilingIfNotGPU()); mfi.isValid(); ++mfi )
              {
               Box gtbx = mfi.growntilebox();
                  Box vbx  = mfi.validbox();
                  int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
                  int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
                  const Array4<      Real>& dst_arr = mapfac[lev][MapFacType::m_x]->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, ilo), ihi);
                      int lj = amrex::min(amrex::max(j, jlo), jhi);
                      dst_arr(i,j,0) = src_arr(li,lj,0);
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

    // Convert the velocities using the map factors
    for ( MFIter mfi(lev_new[Vars::xvel], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();
        const Array4<      Real>& dst_arr = lev_new[Vars::xvel].array(mfi);
        const Array4<const Real>& src_arr = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst_arr(i,j,k) /= src_arr(i,j,0);
        });
    }
    for ( MFIter mfi(lev_new[Vars::yvel], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();
        const Array4<      Real>& dst_arr = lev_new[Vars::yvel].array(mfi);
        const Array4<const Real>& src_arr = mapfac[lev][MapFacType::v_x]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst_arr(i,j,k) /= src_arr(i,j,0);
        });
    }
    for ( MFIter mfi(lev_new[Vars::zvel], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box bx = mfi.tilebox();
        const Array4<      Real>& dst_arr = lev_new[Vars::zvel].array(mfi);
        const Array4<const Real>& src_arr = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            dst_arr(i,j,k) /= src_arr(i,j,0);
        });
    }

    // **************************************************************************
    // Compute min and max of terrain
    // **************************************************************************
    if (compute_terrain_here) {
        if (lev == 0) {
            AMREX_ALWAYS_ASSERT(solverChoice.terrain_type == TerrainType::StaticFittedMesh);
            z_top = compute_terrain_top_and_bottom(mf_PH, mf_PHB, geom[lev].Domain());
        } else {
            amrex::Print() << "Warning: using top of domain set at level 0 which is " << z_top << std::endl;
        }

        // **************************************************************************
        // FillBoundary to populate the internal ghost cells (for averaging)
        // **************************************************************************
         mf_PH.FillBoundary(geom[lev].periodicity());
        mf_PHB.FillBoundary(geom[lev].periodicity());

        // **************************************************************************
        // Initialize the terrain itself
        // **************************************************************************
        init_terrain_from_wrfinput(lev, z_top, boxes_at_level[lev][0], z_phys_nd[lev].get(), mf_PH, mf_PHB);

        // **************************************************************************
        // Initialize the metric quantities
        // **************************************************************************
        make_J  (geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
        make_areas(geom[lev],*z_phys_nd[lev],*ax[lev],*ay[lev],*az[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);
    }

    // **************************************************************************
    // Rebalance the base state if needed
    // **************************************************************************
    if (solverChoice.rebalance_wrfinput) {
        Print() << "Rebalancing the HSE state!\n";
        int ncomp = lev_new[Vars::cons].nComp();
        int k_dom_lo = geom[lev].Domain().smallEnd(2);
        int k_dom_hi = geom[lev].Domain().bigEnd(2);
        Real tol = 1.0e-10;
        Real grav = CONST_GRAV;
        for ( MFIter mfi(lev_new[Vars::cons],TileNoZ()); mfi.isValid(); ++mfi ) {
            Box bx  = mfi.tilebox();
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
                    z_lo = 0.0; // corresponding to p_0
                    z_hi = 0.125 * (z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                   +z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(i,j+1,klo+1) + z_arr(i+1,j+1,klo+1));
                    dz = z_hi - z_lo;

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
                  z_hi = 0.125 * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                 +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1));
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
    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);


    int n_qstate_moist = micro->Get_Qstate_Moist_Size();
    AMREX_ALWAYS_ASSERT(micro->Get_Qstate_NonMoist_Size() == 0);

    bool use_P_eos = (solverChoice.rebalance_wrfinput);

    init_base_state_from_wrfinput(boxes_at_level[lev][0], l_rdOcp, solverChoice.moisture_type, n_qstate_moist,
                                  lev_new[Vars::cons], p_hse, pi_hse, th_hse, qv_hse, r_hse,
                                  mf_PB, mf_P, use_P_eos);

    // FillBoundary to populate the internal ghost cells (no averaging in above call)
     r_hse.FillBoundary(geom[lev].periodicity());
     p_hse.FillBoundary(geom[lev].periodicity());
    pi_hse.FillBoundary(geom[lev].periodicity());
    th_hse.FillBoundary(geom[lev].periodicity());
    qv_hse.FillBoundary(geom[lev].periodicity());

    // *******************************************************************************************
    // Initialize the bdy data
    // *******************************************************************************************
    if (solverChoice.use_real_bcs && (lev == 0))
    {
        if (geom[0].isPeriodic(0) || geom[0].isPeriodic(1) ) {
             amrex::Error("Cannot set periodic lateral boundary conditions when reading in real boundary values");
        }

        if (nc_bdy_file.empty()) {
            amrex::Error("NetCDF boundary file name must be provided via input");
        }

        // Three points are necessary if a relaxation zone is present.
        if (real_width > real_set_width) {
            AMREX_ALWAYS_ASSERT(real_width-real_set_width >= 3);
        }

        bdy_time_interval = read_times_from_wrfbdy(nc_bdy_file,
                                                   bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                                   start_bdy_time);

        // *******************************************************************************************
        // We intentionally only read in the first three slices here ... we will read the rest in
        // as needed during the time stepping procedure
        // *******************************************************************************************
        int ntimes = bdy_data_xlo.size(); ntimes = amrex::min(ntimes, 3);
        for (int itime = 0; itime < ntimes; itime++)
        {
            read_from_wrfbdy(itime,nc_bdy_file,geom[0].Domain(),
                             bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                             real_width);

            if (itime == 0) {
                Print() << "Read in boundary data with width "  << real_width << std::endl;
                Print() << "Running with specification width: " << real_set_width
                        << " and relaxation width: " << real_width - real_set_width << std::endl;
            }

            convert_all_wrfbdy_data(itime, geom[lev].Domain(), bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                    mf_MUB_lev, mf_C1H_lev, mf_C2H_lev,
                                    lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                                    geom[lev], use_moist);
        } // itime

        //
        // Start at the earliest time (read_from_wrfbdy)
        // Note we only have start_bdy_time if at level 0 and init_type == InitType::WRFInput or InitType::Metgrid
        //
        // Note that t_new and t_old carry *elapsed* time, not total time
        //
        if (lev == 0) {
            Print() << "start_bdy_time is " << std::setprecision(timeprecision) << start_bdy_time
                    << " from wrfbdy but note that time variable in simulation is elapsed time" << std::endl;
            t_new[lev] = 0.;
            t_old[lev] = -1.e200;
        } else {
            t_new[lev] = t_new[0];
            t_old[lev] = t_old[0];
        }
    } // use_real_bcs && (lev == 0)

    // *******************************************************************************************
    // Initialize the low data if available
    // *******************************************************************************************
    if ((lev == 0) && !nc_low_file.empty())
    {
        low_time_interval = read_times_from_wrflow(nc_low_file,
                                                   low_data_zlo,
                                                   start_low_time);

        int ntimes = low_data_zlo.size();
        sst_lev[lev].resize(ntimes);
        tsk_lev[lev].resize(ntimes);

        // We can possibly run out of memory if we load all of wrfbdy and all of wrflow
        // Thus we only load the first two time slices here and load more only if needed
        ntimes = amrex::min(ntimes, 2);

        for (int itime(0); itime < ntimes; ++itime) {
            read_from_wrflow(itime, nc_low_file, geom[0].Domain(), low_data_zlo);

            update_sst_tsk(itime, geom[lev], ba2d[lev],
                           sst_lev[lev], tsk_lev[lev],
                           m_SurfaceLayer, low_data_zlo,
                           lev_new[Vars::cons], *mf_PSFC[lev],
                           l_rdOcp, lmask_lev[lev][0], use_moist);
        }
    } // lev == 0 && nc_low_file exists
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
init_base_state_from_wrfinput (const Box& subdomain,
                               const Real l_rdOcp,
                               MoistureType moisture_type,
                               const int& n_qstate_moist,
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
    const auto& dom_lo = lbound(subdomain);
    const auto& dom_hi = ubound(subdomain);

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
                    AMREX_DEVICE_PRINTF("p (%i, %i, %i): %e; p_eos: %e; (qv = %e, rho = %e, rT = %e) \n",
                           i, j, k, Ptot, P_eos, Qv, cons_arr(ii,jj,kk,Rho_comp), RT);
                    amrex::Abort("Initial state is inconsistent with EOS!?");
                }
            }

            // Compute rhse
            Real Rhse_Sum = cons_arr(ii,jj,kk,Rho_comp);
            for (int q_offset(0); q_offset<n_qstate_moist; ++q_offset) {
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
Real
compute_terrain_top_and_bottom (const MultiFab& mf_PH,
                                const MultiFab& mf_PHB,
                                const Box& domain)
{
    Real z_top;

    //
    // For the bottom/top boundary (in that order)
    //
    Gpu::HostVector  <Real> Max_h(3,-1.0e16);
    Gpu::DeviceVector<Real> Max_d(3);
    Gpu::copy(Gpu::hostToDevice, Max_h.begin(), Max_h.end(), Max_d.begin());

    Gpu::HostVector  <Real> Min_h(2, 1.0e16);
    Gpu::DeviceVector<Real> Min_d(2);
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
            amrex::Gpu::Atomic::Min(&(min_d[1]),z_calc_hi);
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
    ParallelDescriptor::ReduceRealMin(Min_h[1]);

    ParallelDescriptor::ReduceRealMax(Max_h[0]);
    ParallelDescriptor::ReduceRealMax(Max_h[1]);
    ParallelDescriptor::ReduceRealMax(Max_h[2]);

    Real terrain_bottom_max = Max_h[0];
    Real terrain_bottom_min = Min_h[0];
    Real terrain_top_max    = Max_h[1];
    Real terrain_top_min    = Min_h[1];
    Real terrain_km1_max    = Max_h[2];

    Print() << "Terrain     has min value    = " << terrain_bottom_min << " and max value = " << terrain_bottom_max << std::endl;
    Print() << "Top of mesh has min value    = " << terrain_top_min    << " and max value = " << terrain_top_max << std::endl;

    // Average the top nodes to define a flat surface at the top
    z_top = 0.5 * (terrain_top_min + terrain_top_max);

    // If this creates a case where z_k < z_{k-1} then we do what we used to do
    if (terrain_km1_max > z_top) {
        amrex::Print() << "Max of second-to-highest row = " << terrain_km1_max <<
                          " which is greater than average of top row so defaulting to alternate approach " << std::endl;
        z_top = 0.5 * (terrain_km1_max + terrain_top_max);
    }

    amrex::Print() << "Warning: ProbHi(2) will be ignored; we are setting top of domain to " << z_top << std::endl;

    return z_top;
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
                            const Box& subdomain,
                            MultiFab* z_phys,
                            const MultiFab& mf_PH,
                            const MultiFab& mf_PHB)
{
    for ( MFIter mfi(*z_phys, false); mfi.isValid(); ++mfi )
    {
        Box gnbx = mfi.growntilebox();

        // This copies from NC_zphys on z-faces to z_phys_nd on nodes
        const Array4<Real      >&      z_arr = z_phys->array(mfi);
        const Array4<Real const>& nc_phb_arr = mf_PHB.const_array(mfi);
        const Array4<Real const>& nc_ph_arr  = mf_PH.const_array(mfi);

        // PHB and PH are on z-faces (0.5 dx/y ahead of zphys)
        Box z_face_box = convert(subdomain,IntVect(0,0,1));

        // Prevent averaging from going into domain ghost cells
        int ilo = z_face_box.smallEnd()[0] + 1;
        int ihi = z_face_box.bigEnd()[0];
        int jlo = z_face_box.smallEnd()[1] + 1;
        int jhi = z_face_box.bigEnd()[1];
        int klo = z_face_box.smallEnd()[2];
        int khi = z_face_box.bigEnd()[2];

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
                z_arr(i, j, k) = 0.25 * ( nc_ph_arr (ii,jj  ,k) + nc_ph_arr (ii-1,jj  ,k) +
                                          nc_ph_arr (ii,jj-1,k) + nc_ph_arr (ii-1,jj-1,k) +
                                          nc_phb_arr(ii,jj  ,k) + nc_phb_arr(ii-1,jj  ,k) +
                                          nc_phb_arr(ii,jj-1,k) + nc_phb_arr(ii-1,jj-1,k) ) / CONST_GRAV;
                z_arr(i, j, k) = z_top;
            } else {
                // Note: wrfinput geopotentials ph, phb are only staggered in the vertical, i.e.,
                //       they have dims (bottom_top_stag, south_north, west_east). On k==klo, we
                //       will end up smoothing the terrain as we average from surface face centers
                //       to nodes.
                z_arr(i, j, k) = 0.25 * ( nc_ph_arr (ii,jj  ,k) + nc_ph_arr (ii-1,jj  ,k) +
                                          nc_ph_arr (ii,jj-1,k) + nc_ph_arr (ii-1,jj-1,k) +
                                          nc_phb_arr(ii,jj  ,k) + nc_phb_arr(ii-1,jj  ,k) +
                                          nc_phb_arr(ii,jj-1,k) + nc_phb_arr(ii-1,jj-1,k) ) / CONST_GRAV;
            }
        });

        // Sanity check
        Print() << "Verifying grid integrity" << std::endl;
        const Box& vbox = mfi.validbox();
        if (vbox.smallEnd(2) == klo) {
            Box z_surf_faces = makeSlab(vbox, 2, klo);
            ParallelFor(z_surf_faces, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                if (z_arr(i,j,k+1) < z_arr(i,j,k)) {
#ifdef AMREX_USE_GPU
                    AMREX_DEVICE_PRINTF("z values at (%d,%d,%d) and k+1 are %f, %f\n",
                           i,j,k, z_arr(i,j,k), z_arr(i,j,k+1));
#else
                    printf("z values at (%d,%d,%d) and k+1 are %f, %f\n",
                           i,j,k, z_arr(i,j,k), z_arr(i,j,k+1));
#endif
                    Error("Grid integrity issue detected");
                }
            });
        } // tile includes zlo
    } // mfi
}
#endif // ERF_USE_NETCDF
