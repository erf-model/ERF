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
#include <ERF_WriteERFBdy.H>
#include <ERF_ReadFromERFBdy.H>

#include "ERF_NodalReconstruction.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

#include "ERF_NCWpsFile.H"

/**
 * Determine which WRF density variables are available in a NetCDF file.
 *
 * @param fname Path to the WRF input file
 * @return 1 when ALT density should be used, 0 when ALB plus AL should be used
 */
bool CheckForDensity (const std::string& fname)
{
    int failed = false;
    int use_alt_density = 0;
    if (ParallelDescriptor::IOProcessor())
    {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);
        int success_al  = ncf.has_var("AL");
        int success_alb = ncf.has_var("ALB");
        int success_alt = ncf.has_var("ALT");

        if (success_al && success_alb && !success_alt) {
            Print() << " Will read density from ALB + AL variables" << std::endl;
        } else if (!success_al && !success_alb && success_alt) {
            Print() << " Will read density from ALT variable" << std::endl;
            use_alt_density = 1;
        } else if ((success_al != success_alb) && !success_alt) {
            Print() << " Density must be read through both ALB and AL when ALT is absent" << std::endl;
            failed = true;
        } else if (!success_al && !success_alb && !success_alt) {
            Print() << " Density is not defined in the netcdf file!" << std::endl;
            failed = true;
        } else if ( (success_al || success_alb) && success_alt) {
            Print() << " Density must be read either through ALB/AL or through ALT, not both" << std::endl;
            failed = true;
        }
    }
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
    ParallelDescriptor::Bcast(&failed, 1, ioproc);
    if (failed) amrex::Abort();

    ParallelDescriptor::Bcast(&use_alt_density, 1, ioproc);
    return use_alt_density;
}

/**
 * Read the subdomain index space from a WRF input file.
 *
 * @param lev Current level (unused).
 * @param fname Path to the WRF input file.
 * @param[out] ratio Grid ratio read from the file.
 * @return Box specifying the subdomain index space.
 */
Box
read_subdomain_from_wrfinput(int /*lev*/,
                             const std::string& fname,
                             int& ratio);

/**
 * Compute the top height of the domain from WRF geopotential data.
 *
 * @param[in] mf_PH MultiFab storing WRF perturbation geopotential data.
 * @param[in] mf_PHB MultiFab storing WRF base-state geopotential data.
 * @param[in] domain Box holding the index space of the computational domain.
 * @return Height assigned to the top of the ERF domain.
 */
Real
compute_terrain_top_and_bottom (const MultiFab& mf_PH,
                                const MultiFab& mf_PHB,
                                const Box& domain);

/**
 * Initialize nodal terrain coordinates from WRF input data.
 *
 * @param lev Current level.
 * @param[in,out] geom Geometry object defining the domain.
 * @param[in] z_top Height assigned to the top of the ERF domain.
 * @param[in] subdomain Box specifying the index space to initialize.
 * @param[out] z_phys MultiFab specifying the node-centered z coordinates.
 * @param[in] NC_PH_fab MultiFab storing WRF perturbation geopotential data.
 * @param[in] NC_PHB_fab MultiFab storing WRF base-state geopotential data.
 * @param[out] dz0_max Maximum first-layer thickness.
 * @param[in] use_wrf_height_grid Whether to use the WRF height grid directly.
 */
void
init_terrain_from_wrfinput (int lev,
                            Geometry& geom,
                            const Real& z_top,
                            const Box& subdomain,
                            MultiFab* z_phys,
                            const MultiFab& NC_PH_fab,
                            const MultiFab& NC_PHB_fab,
                            Real& dz0_max,
                            const bool& use_wrf_height_grid);

/**
 * Initialize hydrostatic base state data from a WRF dataset.
 *
 * @param[in] subdomain Box specifying the index space to initialize.
 * @param[in] l_rdOcp Constant $R_d/c_p$.
 * @param[out] p_hse MultiFab holding the hydrostatic base state pressure.
 * @param[out] pi_hse MultiFab holding the hydrostatic base state Exner pressure.
 * @param[out] th_hse MultiFab holding the hydrostatic base state potential temperature.
 * @param[out] qv_hse MultiFab holding the hydrostatic base state qv.
 * @param[out] r_hse MultiFab holding the hydrostatic base state density.
 * @param[in] mf_PB MultiFab holding WRF data specifying base state pressure.
 * @param[in] mf_ALB Optional MultiFab holding inverse density perturbation data.
 * @param[in] z_phys Optional terrain nodal z-coordinate MultiFab.
 * @param[in] T00 Sea-level base-state temperature.
 * @param[in] P00 Sea-level base-state pressure.
 * @param[in] TLP Base-state lapse rate.
 * @param[in] TISO Isothermal stratosphere temperature.
 * @param[in] TLP_STRAT Stratospheric lapse rate.
 * @param[in] P_STRAT Pressure at the stratosphere transition.
 */
void
init_base_state_from_wrfinput (const Box& subdomain,
                               const Real& l_rdOcp,
                               MultiFab& p_hse,
                               MultiFab& pi_hse,
                               MultiFab& th_hse,
                               MultiFab& qv_hse,
                               MultiFab& r_hse,
                               MultiFab& mf_PB,
                               MultiFab* mf_ALB,
                               MultiFab* z_phys,
                               const Real& T00,
                               const Real& P00,
                               const Real& TLP,
                               const Real& TISO,
                               const Real& TLP_STRAT,
                               const Real& P_STRAT);

/**
 * Read start_time from the first WRF input file.
 *
 * @param lev Integer specifying the current level
 * @param fname Path to the WRF input file
 * @return Epoch time read from the file
 */
double
read_start_time_from_wrfinput (int lev, const std::string& fname)
{
    double NC_epochTime = 0.0;
    const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor()) {
        // Read the time stamps
        using CharArray = NDArray<char>;
        Vector<CharArray> array_ts(1);
        Vector<int> success(1);
        ReadNetCDFFile(fname, {"Times"}, array_ts, success);

        auto dateStrLen = array_ts[0].get_vshape()[1];
        const char* time_stamp_data = array_ts[0].get_data();

        std::string date(time_stamp_data, time_stamp_data + dateStrLen);

        auto epochTime = getEpochTime(date, dateTimeFormat);
        Print() << "  wrfinput datetime 0 : " << date << " " << epochTime << std::endl;
        NC_epochTime = static_cast<double>(epochTime);

        Print() << "Have read start_time string at level "<< lev << " is " << date << std::endl;
        Print() << "Have read start_time number at level "<< lev << " is " << NC_epochTime << std::endl;
    }

    ParallelDescriptor::Bcast(&NC_epochTime,1,ParallelDescriptor::IOProcessorNumber());

    return NC_epochTime;
}

/**
 * Read WRF base-state thermodynamic parameters from a NetCDF file.
 *
 * @param fname Path to the WRF input file
 * @param T00 Sea-level base-state temperature
 * @param P00 Sea-level base-state pressure
 * @param TLP Base-state lapse rate
 * @param TISO Isothermal stratosphere temperature
 * @param TLP_STRAT Stratospheric lapse rate
 * @param P_STRAT Pressure at the stratosphere transition
 */
void
read_base_state_params_from_wrfinput (const std::string& fname,
                                      Real& T00,
                                      Real& P00,
                                      Real& TLP,
                                      Real& TISO,
                                      Real& TLP_STRAT,
                                      Real& P_STRAT)
{
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);

        // Remember what was passed in so we can fall back to it below
        const Real T00_def = T00;
        const Real P00_def = P00;
        const Real TLP_def = TLP;
        const Real TISO_def = TISO;
        const Real TLP_STRAT_def = TLP_STRAT;
        const Real P_STRAT_def   = P_STRAT;

        std::vector<size_t> shape;
        std::vector<size_t> start;
        auto read_scalar = [&] (const std::string& name, Real& val)
        {
            if (!ncf.has_var(name)) return;
            Print() << "Reading " << name << " from wrfinput\n";
            shape = ncf.var(name).shape();
            start.clear();
            start.resize(shape.size(), 0);
            ncf.var(name).get(&val, start, shape);
        };

        read_scalar("T00"      , T00);
        read_scalar("P00"      , P00);
        read_scalar("TLP"      , TLP);
        read_scalar("TISO"     , TISO);
        read_scalar("TLP_STRAT", TLP_STRAT);
        read_scalar("P_STRAT"  , P_STRAT);

        ncf.close();

        // Idealized WRF cases (and some hand-built wrfinput files) declare these
        // variables but leave them zero-filled.  T00, P00 and TISO are absolute
        // temperatures/pressures, so a non-positive value is never meaningful; if
        // any of them is bad we discard the whole group and keep the defaults.
        // (TLP, TLP_STRAT and P_STRAT may legitimately be zero, so they are only
        // reset alongside a bad T00/P00/TISO.)
        const bool params_ok = std::isfinite(T00)  && (T00  > Real(0)) &&
                               std::isfinite(P00)  && (P00  > Real(0)) &&
                               std::isfinite(TISO) && (TISO > Real(0)) &&
                               std::isfinite(TLP)  && std::isfinite(TLP_STRAT) &&
                               std::isfinite(P_STRAT);

        if (!params_ok) {
            Print() << "WARNING: WRF base state parameters read from " << fname
                    << " are invalid: (T00, P00, TLP, TISO, TLP_STRAT, P_STRAT) = ("
                    << T00 << ", " << P00 << ", " << TLP << ", " << TISO << ", "
                    << TLP_STRAT << ", " << P_STRAT << ")\n";
            Print() << "         T00, P00 and TISO must all be positive and finite; "
                       "reverting to ERF defaults.\n";
            T00 = T00_def;   P00 = P00_def;   TLP = TLP_def;
            TISO = TISO_def; TLP_STRAT = TLP_STRAT_def; P_STRAT = P_STRAT_def;
        }

        Print() << "WRF base state parameters (T00, P00, TLP, TISO, TLP_STRAT, P_STRAT) are: ("
                << T00 << ", " << P00 << ", " << TLP << ", " << TISO << ", " << TLP_STRAT << ", "
                << P_STRAT << ") \n";
    }

    ParallelDescriptor::Bcast(&T00,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&P00,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&TLP,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&TISO,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&TLP_STRAT,1,ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&P_STRAT,1,ParallelDescriptor::IOProcessorNumber());
}

/**
 * ERF function that initializes data from a WRF dataset
 *
 * @param lev Integer specifying the current level
 * @param mf_PSFC MultiFab storing surface pressure for this level
 */
void
ERF::init_from_wrfinput (int lev, MultiFab& mf_PSFC_lev)
{
    if (nc_init_file.empty()) {
        amrex::Error("NetCDF initialization file name must be provided via input");
    }

    bool use_moist = (solverChoice.moisture_type != MoistureType::None);
    bool use_lsm   = (solverChoice.lsm_type != LandSurfaceType::None);

    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<std::string> NC_names;
    NC_names.push_back("ALB");       // 0 DO RHO FIRST
    NC_names.push_back("AL");        // 1 DO RHO FIRST
    NC_names.push_back("ALT");       // 2 DO RHO FIRST
    NC_names.push_back("U");         // 3
    NC_names.push_back("V");         // 4
    NC_names.push_back("W");         // 5
    NC_names.push_back("THM");       // 6
    NC_names.push_back("PH");        // 7
    NC_names.push_back("PHB");       // 8
    NC_names.push_back("PB");        // 9
    NC_names.push_back("P");         // 10
    NC_names.push_back("PSFC");      // 11
    NC_names.push_back("MUB");       // 12
    NC_names.push_back("MAPFAC_U");  // 13
    NC_names.push_back("MAPFAC_V");  // 14
    NC_names.push_back("MAPFAC_M");  // 15
    NC_names.push_back("SST");       // 16
    NC_names.push_back("TSK");       // 17
    NC_names.push_back("LANDMASK");  // 18
    NC_names.push_back("C1H");       // 19
    NC_names.push_back("C2H");       // 20
    NC_names.push_back("RDNW");      // 21
    NC_names.push_back("XLAT_V");    // 22
    NC_names.push_back("XLONG_U");   // 23
    if (use_moist) {
        NC_names.push_back("QVAPOR"); // 24
        NC_names.push_back("QCLOUD"); // 25
        NC_names.push_back("QRAIN");  // 26
    }
    NC_names.push_back("IVGTYP");     // 27
    NC_names.push_back("ISLTYP");     // 28
    if (use_lsm) {
        NC_names.push_back("TSLB");   // 29
        NC_names.push_back("SMOIS");  // 30
        NC_names.push_back("SH2O");   // 31
        NC_names.push_back("LAI");    // 32
        NC_names.push_back("ZS");     // 33
        NC_names.push_back("DZS");    // 34
        NC_names.push_back("VEGFRA"); // 35
        NC_names.push_back("TMN");    // 36
        NC_names.push_back("SHDMIN"); // 37
        NC_names.push_back("SHDMAX"); // 38

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
    Vector<Vector<FArrayBox>> NC_fab_var;
    NC_fab_var.resize(num_boxes_at_level[lev]);
    for (int idx(0); idx < num_boxes_at_level[lev]; ++ idx) { NC_fab_var[idx].resize(nvar); }

    auto& lev_new = vars_new[lev];

    // NOTE: These temporaries keep us from overwriting the lev==0 wrf data that is
    //       stored for the BDY operations.
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

    MultiFab mf_PH, mf_PB, mf_P;      // For geopotential and base state
    std::unique_ptr<MultiFab> mf_ALB; // For density

    // Read base state params (used if ALB is not read)
    Real T00 = Real(290.0);
    Real P00 = p_0;
    Real TLP = Real(50.0);
    Real TISO = Real(200.0);
    Real TLP_STRAT = Real(-11.0);
    Real P_STRAT   = zero;
    read_base_state_params_from_wrfinput(nc_init_file[lev][0],
                                         T00, P00, TLP, TISO,
                                         TLP_STRAT, P_STRAT);

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

        // Check if the density variable is available and uniquely defined
        // for this specific input file.
        int use_alt_density = CheckForDensity(nc_init_file[lev][idx]);

        int ratio_from_file;
        Box subdomain_to_read = read_subdomain_from_wrfinput(lev, nc_init_file[lev][idx], ratio_from_file);
        Print() << "Box in file " << subdomain_to_read << "\n";

        Box subdomain_to_fill = boxes_at_level[lev][idx];
        Print() << "Box to fill " << subdomain_to_fill << "\n";

        for (int ivar = 0; ivar < nvar; ++ ivar) {
            auto var_name = NC_names[ivar];

            // Once the density mode is selected, skip the unused path entirely.
            // This avoids trying to read optional density variables that are absent.
            bool skip_density_var =
                (use_alt_density == 1 && (var_name == "ALB" || var_name == "AL")) ||
                (use_alt_density == 0 &&  var_name == "ALT");
            if (skip_density_var) { continue; }

            Print() << "Checking for " << var_name << " ...";

            int success, use_theta_m;
            read_from_wrfinput(lev, subdomain_to_read, nc_init_file[lev][idx],
                               NC_fab_var[idx][ivar], var_name, geom[lev],
                               use_theta_m, success);

            auto& var_fab_from_file = NC_fab_var[idx][ivar];
            bool has_fallback_behavior =
                (var_name == "U")      || (var_name == "V")      || (var_name == "W")      ||
                (var_name == "THM")    || (var_name == "QVAPOR") || (var_name == "QCLOUD") ||
                (var_name == "QRAIN")  || (var_name == "PH")     || (var_name == "PHB");
            if (!success && !has_fallback_behavior) {
                amrex::Abort(std::string("ERF::init_from_wrfinput: failed to read required variable " + var_name).c_str());
            }

            FArrayBox var_fab;
            FArrayBox var_fab_crse;

            if (success) {
                // This shift occurs only when the coarser grid is at least at level 1,
                // because the indices are always given relative to that coarser "domain"
                if (lev > 1) {
                    Box shift_by_box(subdomains[lev][0].minimalBox());
                    IntVect shift_by(shift_by_box.smallEnd());
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        shift_by[i] -= var_fab_from_file.box().smallEnd(i);
                    }
                    var_fab_from_file.shift(shift_by);
                }


                // In the case where the array is 1D in the z-direction, the destination box needs to also
                //    be 1D in the z-direction, but with (i,j) corresponding to the low corner of the box
                //    to be filled
                int nx = var_fab_from_file.box().length(0);
                int ny = var_fab_from_file.box().length(1);
                int nz = var_fab_from_file.box().length(2);
                Box subdomain_tmp(subdomain_to_fill);
                if (nx == 1 and ny == 1) {
                    subdomain_tmp.setBig(0,subdomain_tmp.smallEnd(0));
                    subdomain_tmp.setBig(1,subdomain_tmp.smallEnd(1));
                }
                // For 2D variables like MAPFAC_U, MAPFAC_V, PSFC, MUB, etc.
                if (nz == 1) {
                    subdomain_tmp.setBig(2,subdomain_tmp.smallEnd(2));
                }

                Box subdomain_to_fill_typed(convert(subdomain_tmp,var_fab_from_file.box().ixType()));
                Box subdomain_crse(subdomain_to_fill_typed);
                if (lev > 0) {
                    subdomain_crse.coarsen(IntVect(1,1,ref_ratio[lev-1][2]));
                    if (ref_ratio[lev-1][2] > 1) {
                        amrex::Abort("This pathway in init_from_wrfinput not ready yet");
                    }
                }
#ifdef AMREX_USE_GPU
                var_fab.resize(subdomain_to_fill_typed, 1, amrex::The_Pinned_Arena());
                var_fab_crse.resize(subdomain_crse, 1, amrex::The_Pinned_Arena());
#else
                var_fab.resize(subdomain_to_fill_typed, 1);
                var_fab_crse.resize(subdomain_crse, 1);
#endif
                Box intersection = var_fab.box() & var_fab_from_file.box();
                if (intersection.ok()) {
#if 0
                    var_fab.template copy<RunOn::Device>(var_fab_from_file,intersection,0,intersection,0,1);
#else
                    if (lev == 0 || ref_ratio[lev-1][2] == 1) {
                        var_fab.template copy<RunOn::Device>(var_fab_from_file,intersection,0,intersection,0,1);
                    } else {
                        var_fab_crse.template copy<RunOn::Device>(var_fab_from_file,intersection,0,intersection,0,1);
                    }
#endif
                } else if (nx == 1 and ny == 1) {
                    Print() << " Copying 1D FAB from " << var_fab_from_file.box() << " to " << var_fab.box() << std::endl;
#if 0
                    var_fab.template copy<RunOn::Device>(var_fab_from_file,var_fab_from_file.box(),0,var_fab.box(),0,1);
#else
                    if (lev == 0 || ref_ratio[lev-1][2] == 1) {
                        var_fab.template copy<RunOn::Device>(var_fab_from_file,var_fab_from_file.box(),0,var_fab.box(),0,1);
                    } else {
                        var_fab_crse.template copy<RunOn::Device>(var_fab_from_file,var_fab_from_file.box(),0,var_fab.box(),0,1);
                    }
#endif
                } else {
                    Print() <<"var_fab_crse.box()      " << subdomain_crse << std::endl;
                    Print() <<"var_fab_from_file.box() " << var_fab_from_file.box() << std::endl;
                    amrex::Error("ERF::init_from_wrfinput: Region we want not contained in region we have");
                }
            }

            // Initialize rho =  1/(ALB + AL)
            if ( (use_alt_density == 0) && (var_name == "ALB") ) {
                mf_ALB = std::make_unique<MultiFab>(ba,dm,1,ng);
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                {
                    FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
                    cons_fab.template copy<RunOn::Device>(var_fab, 0, Rho_comp, 1);
                    FArrayBox &ALB_fab = (*mf_ALB)[mfi];
                    ALB_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
                }

            } if ( (use_alt_density == 0) && (var_name == "AL") ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                {
                    FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
                    Box vbx = cons_fab.box(); vbx.grow(-ng);

                    // Add "AL" to "ALB" before inverting
                    cons_fab.template   plus<RunOn::Device>(var_fab, 0, Rho_comp, 1);
                    cons_fab.template invert<RunOn::Device>(one, vbx, Rho_comp, 1);
                }
            }

            // OR Initialize rho =  1/ALT
            if ( (use_alt_density == 1) && (var_name == "ALT") ) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
                {
                    FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
                    Box vbx = cons_fab.box(); vbx.grow(-ng);

                    // "ALT" holds the full 1/density so we can invert here
                    cons_fab.template copy<RunOn::Device>(var_fab, 0, Rho_comp, 1);
                    cons_fab.template invert<RunOn::Device>(one, vbx, Rho_comp, 1);
                }
            }

            if ( var_name == "THM") {
                const Real wrf_theta_ref = Real(300.0);
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
                        cur_fab->template setVal<RunOn::Device>(0,cur_fab->box(),0,1);
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
                        var_fab.template mult<RunOn::Device>(RvoRd);
                        var_fab.template plus<RunOn::Device>(one);
                        var_fab.template invert<RunOn::Device>(one);
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
                        lev_new[Vars::cons].setVal(0,icomp,1);
                    } else {
                        amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                    }
                }

                var_fab.clear();
            } // valid var (not rho)

          bool lat_periodic = (geom[lev].isPeriodic(0) && geom[lev].isPeriodic(1));
          int i_lo = boxes_at_level[lev][0].smallEnd(0); int i_hi = boxes_at_level[lev][0].bigEnd(0);
          int j_lo = boxes_at_level[lev][0].smallEnd(1); int j_hi = boxes_at_level[lev][0].bigEnd(1);

          if ( var_name == "PH" ) {
              if (success) {
                  // NOTE: We call FillBoundary on mf_PH below
                  auto& ba_w = lev_new[Vars::zvel].boxArray();
                  mf_PH.define(ba_w, dm, 1, IntVect(ngz[0],ngz[1],0));
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                  for ( MFIter mfi(mf_PH, false); mfi.isValid(); ++mfi )
                  {
                      Box gtbx = mfi.growntilebox();
                      const Array4<      Real>& dst_arr = mf_PH.array(mfi);
                      const Array4<const Real>& src_arr = var_fab.const_array();
                      ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                          int li = amrex::min(amrex::max(i, i_lo), i_hi);
                          int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                          dst_arr(i,j,k) = src_arr(li,lj,k);
                      });
                  }
                  var_fab.clear();
              } else {
                  amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                  compute_terrain_here = false;
              }
          } else if ( var_name == "PHB" ) {
              if (success) {
                  // NOTE: We call FillBoundary on PHB below
                  auto& ba_w = lev_new[Vars::zvel].boxArray();
                  if (lev == 0) {
                      mf_PHB = wrf_PHB.get();
                  } else {
                      PHB_tmp.define(ba_w, dm, 1, IntVect(ngz[0],ngz[1],0));
                      mf_PHB = &PHB_tmp;
                  }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                  for ( MFIter mfi(*mf_PHB, false); mfi.isValid(); ++mfi )
                  {
                      Box gtbx = mfi.growntilebox();
                      const Array4<      Real>& dst_arr = mf_PHB->array(mfi);
                      const Array4<const Real>& src_arr = var_fab.const_array();
                      ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                      {
                          int li = amrex::min(amrex::max(i, i_lo), i_hi);
                          int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                          dst_arr(i,j,k) = src_arr(li,lj,k);
                      });
                  }
                  var_fab.clear();
              } else {
                  amrex::Print() << "Ignoring " << var_name << " since we aren't using it ... DONE" << std::endl;
                  compute_terrain_here = false;
              }
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
              Real pmax = mf_PSFC_lev.max(0);
              if (pmax == zero) {
                  amrex::Print() << " PSFC read in had max of 0; replacing it by 1e5 everywhere" << std::endl;
                  mf_PSFC_lev.setVal(p_0);
              }
              var_fab.clear();
          } else if ( var_name == "MUB" ) {
              if (lev == 0) {
                  mf_MUB = wrf_MUB.get();
              } else {
                  MUB_tmp.define(ba2d[lev], dm, 1, IntVect(ngz[0],ngz[1],0));
                  mf_MUB = &MUB_tmp;
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mf_MUB, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = (*mf_MUB)[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "C1H" ) {
              if (lev == 0) {
                  mf_C1H = wrf_C1H.get();
              } else {
                  C1H_tmp.define(ba1d[lev], dm, 1, IntVect(ngz[0],ngz[1],0));
                  mf_C1H = &C1H_tmp;
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mf_C1H, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = (*mf_C1H)[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "C2H" ) {
              if (lev == 0) {
                  mf_C2H = wrf_C2H.get();
              } else {
                  C2H_tmp.define(ba1d[lev], dm, 1, IntVect(ngz[0],ngz[1],0));
                  mf_C2H = &C2H_tmp;
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mf_C2H, false); mfi.isValid(); ++mfi )
              {
                FArrayBox &cur_fab = (*mf_C2H)[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          } else if ( var_name == "RDNW" ) {
              if (lev == 0) {
                  mf_RDNW = wrf_RDNW.get();
              } else {
                  RDNW_tmp.define(ba1d[lev], dm, 1, IntVect(ngz[0],ngz[1],0));
                  mf_RDNW = &RDNW_tmp;
              }
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
              for ( MFIter mfi(*mf_RDNW, false); mfi.isValid(); ++mfi )
              {
                  FArrayBox &cur_fab = (*mf_RDNW)[mfi];
                cur_fab.template copy<RunOn::Device>(var_fab, 0, 0, 1);
              }
              var_fab.clear();
          }

          // Initialize Latitude & Coriolis factors
          if ( var_name == "XLAT_V" ) {
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

                      Real lat_rad = dst_arr(i,j,0) * (PI/Real(180.));
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
                  var_fab.template setVal<RunOn::Device>(1);
              }
              if (lat_periodic) {
                  Print() << "MAPFAC_U resetting to 1 with lateral periodic BCs!\n";
                  var_fab.template setVal<RunOn::Device>(1);
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
                  var_fab.template setVal<RunOn::Device>(1);
              }
              if (lat_periodic) {
                  Print() << "MAPFAC_V resetting to 1 with lateral periodic BCs!\n";
                  var_fab.template setVal<RunOn::Device>(1);
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
                  var_fab.template setVal<RunOn::Device>(1);
              }
              if (lat_periodic) {
                  Print() << "MAPFAC_M resetting to 1 with lateral periodic BCs!\n";
                  var_fab.template setVal<RunOn::Device>(1);
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
            z_top = compute_terrain_top_and_bottom(mf_PH, *mf_PHB, geom[lev].Domain());
        } else {
            amrex::Print() << "Warning: using top of domain set at level 0 which is " << z_top << std::endl;
        }

        // **************************************************************************
        // Initialize the terrain itself
        // **************************************************************************
        Real dz0_max;
        init_terrain_from_wrfinput(lev, geom[lev], z_top, boxes_at_level[lev][0], z_phys_nd[lev].get(),
                                   mf_PH, *mf_PHB, dz0_max, solverChoice.use_wrf_height_grid);
        z_phys_nd[lev]->FillBoundary(geom[lev].periodicity());
        if (!solverChoice.use_wrf_height_grid) {
#ifdef AMREX_USE_FLOAT
            const Real tol = Real(1.e-4);
#else
            const Real tol = Real(1.e-8);
#endif
            int max_iter = 20;

            int iter   = 0;
            Real Nz    = static_cast<Real>(zlevels_stag[lev].size() - 1);
            Real SFact = Real(1.03);
            Real F     = dz0_max * ( (std::pow(SFact,Nz) - one) / (SFact - one) ) - z_top;
            while (std::fabs(F)>tol && iter<max_iter) {
                Real dFdSF = dz0_max * ( Nz * std::pow(SFact,Nz-one) * (SFact - one) - std::pow(SFact,Nz) + one )
                           / std::pow(SFact-one,two);
                SFact     -= F/dFdSF;
                SFact      = std::max(one+tol,SFact);
                F          = dz0_max * ( (std::pow(SFact,Nz) - one) / (SFact - one) ) - z_top;
                ++iter;
            }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::fabs(F) <= tol,
                "Newton iterations to determine the grid stretching factor failed!\n");

            Print() << "Building an ERF grid with dz0: " << dz0_max <<
                " and stretching factor: " << SFact << "\n";
            Real dz = dz0_max;
            zlevels_stag[lev][0] = zero;
            for (int k(1); k<zlevels_stag[lev].size(); ++k) {
                zlevels_stag[lev][k] = zlevels_stag[lev][k-1] + dz;
                dz *= SFact;
            }
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::fabs(zlevels_stag[lev].back() - z_top) <= tol,
                "Top of zlevels_stag does not match z_top!\n");
            make_terrain_fitted_coords(lev, geom[lev], *z_phys_nd[lev], zlevels_stag[lev], phys_bc_type);
        }

        // **************************************************************************
        // Initialize the metric quantities
        // **************************************************************************
        make_J    (geom[lev],*z_phys_nd[lev],*detJ_cc[lev]);
        make_areas(geom[lev],*z_phys_nd[lev],*ax[lev],*ay[lev],*az[lev]);
        make_zcc  (geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);

        // **************************************************************************
        // Interpolate the data to the grid
        //
        // NOTE: When keeping the WRF grid, we really only need to interpolate
        //       the highest level values, due to enforcing z_top. When using
        //       the ERF grid, interpolation is required everywhere.
        //
        // NOTE: z_cc must be averaged to the destination location when interpolating.
        //        This is due to the fact that z_cc is conserved WRT WRF heights.
        // **************************************************************************
        int ncons = lev_new[Vars::cons].nComp();
        int imin  = mf_PH.boxArray().minimalBox().smallEnd(0);
        int imax  = mf_PH.boxArray().minimalBox().bigEnd(0);
        int jmin  = mf_PH.boxArray().minimalBox().smallEnd(1);
        int jmax  = mf_PH.boxArray().minimalBox().bigEnd(1);

        int klo   = geom[lev].Domain().smallEnd(2);
        int khi   = geom[lev].Domain().bigEnd(2);

        MultiFab cons_tmp(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), ncons, 0);
        MultiFab xvel_tmp(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1    , 0);
        MultiFab yvel_tmp(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1    , 0);

        MultiFab::Copy(cons_tmp, lev_new[Vars::cons], 0, 0, ncons, 0);
        MultiFab::Copy(xvel_tmp, lev_new[Vars::xvel], 0, 0, 1    , 0);
        MultiFab::Copy(yvel_tmp, lev_new[Vars::yvel], 0, 0, 1    , 0);

        for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx  = mfi.tilebox();
            const Box& bxx = mfi.tilebox(IntVect(1,0,0));
            const Box& bxy = mfi.tilebox(IntVect(0,1,0));

            const Array4<      Real>& cons_arr = lev_new[Vars::cons].array(mfi);
            const Array4<      Real>& xvel_arr = lev_new[Vars::xvel].array(mfi);
            const Array4<      Real>& yvel_arr = lev_new[Vars::yvel].array(mfi);

            const Array4<const Real>& cons_tmp_arr = cons_tmp.array(mfi);
            const Array4<const Real>& xvel_tmp_arr = xvel_tmp.array(mfi);
            const Array4<const Real>& yvel_tmp_arr = yvel_tmp.array(mfi);

            const Array4<const Real>&   ph_arr = mf_PH.const_array(mfi);
            const Array4<const Real>&  phb_arr = mf_PHB->const_array(mfi);

            const Array4<const Real>& z_cc_arr = z_phys_cc[lev]->const_array(mfi);
            const Array4<const Real>& z_nd_arr = z_phys_nd[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int kstart = klo;
                Real z_dst = z_cc_arr(i,j,k);
                Real z_lo_src = Real(0.5) *
                    ( ph_arr(i,j,kstart  ) + phb_arr(i,j,kstart  ) +
                      ph_arr(i,j,kstart+1) + phb_arr(i,j,kstart+1)) / CONST_GRAV;

                bool found = false;
                int kend   = kstart;
                for (int lk(kstart+1); lk<=khi; ++lk) {
                    Real z_hi_src = Real(0.5) *
                        ( ph_arr(i,j,lk  ) + phb_arr(i,j,lk  ) +
                          ph_arr(i,j,lk+1) + phb_arr(i,j,lk+1)) / CONST_GRAV;
                    if (z_dst >= z_lo_src && z_dst < z_hi_src) {
                        found = true;
                        kend  = lk;
                        break;
                    }
                    z_lo_src = z_hi_src;
                    kstart   = lk;
                }

                if (found) {
                    Real z_hi_src = Real(0.5) *
                        (ph_arr(i,j,kend  ) + phb_arr(i,j,kend  ) +
                         ph_arr(i,j,kend+1) + phb_arr(i,j,kend+1)) / CONST_GRAV;
                    Real dz_rat = (z_dst - z_lo_src) / (z_hi_src - z_lo_src);
                    for (int icons(Rho_comp); icons<ncons; ++icons) {
                        Real cons_hi = cons_tmp_arr(i,j,kend  ,icons);
                        Real cons_lo = cons_tmp_arr(i,j,kstart,icons);
                        cons_arr(i,j,k,icons) = ( cons_hi - cons_lo ) * dz_rat + cons_lo;
                    }
                }
            });


            ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Prevent averaging outside domain
                int ii  = std::max(std::min(i  ,imax),imin);
                int iim = std::max(std::min(i-1,imax),imin);

                int kstart = klo;
                Real z_dst = Real(0.25) * ( z_nd_arr(i,j,k  ) + z_nd_arr(i,j+1,k  )
                                          + z_nd_arr(i,j,k+1) + z_nd_arr(i,j+1,k+1) );
                Real z_lo_src = Real(0.25) *
                    ( ph_arr(ii ,j,kstart  ) + phb_arr(ii ,j,kstart  ) +
                      ph_arr(ii ,j,kstart+1) + phb_arr(ii ,j,kstart+1) +
                      ph_arr(iim,j,kstart  ) + phb_arr(iim,j,kstart  ) +
                      ph_arr(iim,j,kstart+1) + phb_arr(iim,j,kstart+1) ) / CONST_GRAV;

                bool found = false;
                int kend   = kstart;
                for (int lk(kstart+1); lk<=khi; ++lk) {
                    Real z_hi_src = Real(0.25) *
                    ( ph_arr(ii ,j,lk  ) + phb_arr(ii ,j,lk  ) +
                      ph_arr(ii ,j,lk+1) + phb_arr(ii ,j,lk+1) +
                      ph_arr(iim,j,lk  ) + phb_arr(iim,j,lk  ) +
                      ph_arr(iim,j,lk+1) + phb_arr(iim,j,lk+1) ) / CONST_GRAV;
                    if (z_dst >= z_lo_src && z_dst < z_hi_src) {
                        found = true;
                        kend  = lk;
                        break;
                    }
                    z_lo_src = z_hi_src;
                    kstart   = lk;
                }

                if (found) {
                    Real z_hi_src = Real(0.25) *
                    ( ph_arr(ii ,j,kend  ) + phb_arr(ii ,j,kend  ) +
                      ph_arr(ii ,j,kend+1) + phb_arr(ii ,j,kend+1) +
                      ph_arr(iim,j,kend  ) + phb_arr(iim,j,kend  ) +
                      ph_arr(iim,j,kend+1) + phb_arr(iim,j,kend+1) ) / CONST_GRAV;
                    Real dz_rat = (z_dst - z_lo_src) / (z_hi_src - z_lo_src);
                    Real xvel_hi = xvel_tmp_arr(i,j,kend  );
                    Real xvel_lo = xvel_tmp_arr(i,j,kstart);
                    xvel_arr(i,j,k) = ( xvel_hi - xvel_lo ) * dz_rat + xvel_lo;
                }
            });

            ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Prevent averaging outside domain
                int jj  = std::max(std::min(j  ,jmax),jmin);
                int jjm = std::max(std::min(j-1,jmax),jmin);

                int kstart = klo;
                Real z_dst = Real(0.25) * ( z_nd_arr(i,j,k  ) + z_nd_arr(i+1,j,k  )
                                          + z_nd_arr(i,j,k+1) + z_nd_arr(i+1,j,k+1) );
                Real z_lo_src = Real(0.25) *
                    ( ph_arr(i,jj ,kstart  ) + phb_arr(i,jj ,kstart  ) +
                      ph_arr(i,jj ,kstart+1) + phb_arr(i,jj ,kstart+1) +
                      ph_arr(i,jjm,kstart  ) + phb_arr(i,jjm,kstart  ) +
                      ph_arr(i,jjm,kstart+1) + phb_arr(i,jjm,kstart+1) ) / CONST_GRAV;

                bool found = false;
                int kend   = kstart;
                for (int lk(kstart+1); lk<=khi; ++lk) {
                    Real z_hi_src = Real(0.25) *
                    ( ph_arr(i,jj ,lk  ) + phb_arr(i,jj ,lk  ) +
                      ph_arr(i,jj ,lk+1) + phb_arr(i,jj ,lk+1) +
                      ph_arr(i,jjm,lk  ) + phb_arr(i,jjm,lk  ) +
                      ph_arr(i,jjm,lk+1) + phb_arr(i,jjm,lk+1) ) / CONST_GRAV;
                    if (z_dst >= z_lo_src && z_dst < z_hi_src) {
                        found = true;
                        kend  = lk;
                        break;
                    }
                    z_lo_src = z_hi_src;
                    kstart   = lk;
                }

                if (found) {
                    Real z_hi_src = Real(0.25) *
                    ( ph_arr(i,jj ,kend  ) + phb_arr(i,jj ,kend  ) +
                      ph_arr(i,jj ,kend+1) + phb_arr(i,jj ,kend+1) +
                      ph_arr(i,jjm,kend  ) + phb_arr(i,jjm,kend  ) +
                      ph_arr(i,jjm,kend+1) + phb_arr(i,jjm,kend+1) ) / CONST_GRAV;
                    Real dz_rat = (z_dst - z_lo_src) / (z_hi_src - z_lo_src);
                    Real yvel_hi = yvel_tmp_arr(i,j,kend  );
                    Real yvel_lo = yvel_tmp_arr(i,j,kstart);
                    yvel_arr(i,j,k) = ( yvel_hi - yvel_lo ) * dz_rat + yvel_lo;
                }
            });

        } // mfi
    }

    // **************************************************************************
    // Rebalance the WRF state if needed
    // **************************************************************************
    if (solverChoice.rebalance_wrf_input) {
        Print() << "The state read from WRF is being rebalanced!\n";

        MultiFab rho (lev_new[Vars::cons], make_alias, Rho_comp, 1);

        MultiFab theta(rho.boxArray(), rho.DistributionMap(), 1, 1);
        MultiFab::Copy(theta, lev_new[Vars::cons], RhoTheta_comp, 0, 1, 1);
        MultiFab::Divide(theta, vars_new[lev][Vars::cons], Rho_comp , 0, 1, 1);

        MultiFab qv(rho.boxArray(), rho.DistributionMap(), 1, 1);
        MultiFab::Copy(qv, lev_new[Vars::cons], RhoQ1_comp, 0, 1, 1);
        MultiFab::Divide(qv, vars_new[lev][Vars::cons], Rho_comp , 0, 1, 1);

        MultiFab qt(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, 0);
        int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
        make_qt(lev_new[Vars::cons], qt, n_qstate_into_total);

        bool maintain_Th = true;
        rebalance_columns(rho, theta, qv, qt, z_phys_nd[lev].get(), geom[lev], maintain_Th);

        // Update (rho qv) in the state
        MultiFab::Multiply(qv, rho, 0, 0, 1, 1);
        MultiFab::Copy(lev_new[Vars::cons], qv, 0, RhoQ1_comp, 1, 1);

        // Update (rho theta) in the state
        MultiFab::Multiply(theta, rho, 0, 0, 1, 1);
        MultiFab::Copy(lev_new[Vars::cons], theta, 0, RhoTheta_comp, 1, 1);
    }

    // **************************************************************************
    // Initialize the base state
    // **************************************************************************
    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    init_base_state_from_wrfinput(boxes_at_level[lev][0], l_rdOcp,
                                  p_hse, pi_hse, th_hse, qv_hse, r_hse,
                                  mf_PB, mf_ALB.get(), z_phys_nd[lev].get(),
                                  T00, P00, TLP, TISO, TLP_STRAT, P_STRAT);

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

        // Check for erfbdy file.
        std::string erfbdy_header = erfbdy_file + "/Header";
        use_erfbdy = FileSystem::Exists(erfbdy_header);
        if (use_erfbdy || write_erfbdy) nvars_erfbdy = WRFBdyVars::NumTypes;

        // Path 1: Load from existing erfbdy file.
        if (use_erfbdy) {
            Print() << "Loading boundary data from erfbdy file: " << erfbdy_file << std::endl;

            // Read metadata and times from erfbdy.
            int ntimes_erfbdy;
            Vector<double> bdy_times;
            bdy_time_interval = read_times_from_erfbdy(erfbdy_file,
                                                       ntimes_erfbdy, nvars_erfbdy, real_width,
                                                       bdy_times, start_bdy_time, final_bdy_time);

            Print() << "erfbdy file contains " << ntimes_erfbdy << " time slices" << std::endl;
            Print() << "start_bdy_time = " << start_bdy_time << std::endl;
            Print() << "final_bdy_time = " << final_bdy_time << std::endl;
            Print() << "bdy_time_interval = " << bdy_time_interval << std::endl;

            bdy_data_xlo.resize(ntimes_erfbdy);
            bdy_data_xhi.resize(ntimes_erfbdy);
            bdy_data_ylo.resize(ntimes_erfbdy);
            bdy_data_yhi.resize(ntimes_erfbdy);

            // Load the first 2 times for simulation initialization.
            for (int itime = 0; itime < std::min(2, ntimes_erfbdy); ++itime) {
                read_from_erfbdy(itime, erfbdy_file,
                                 bdy_data_xlo, bdy_data_xhi,
                                 bdy_data_ylo, bdy_data_yhi,
                                 nvars_erfbdy, real_width);
                Print() << "Loaded erfbdy time slice " << itime << std::endl;
            }

            Print() << "Read in boundary data with width "  << real_width << std::endl;
            Print() << "Running with relaxation width: " << real_width << std::endl;
        }
        // Path 2: Load from wrfbdy and optionally write to erfbdy.
        else {
            if (nc_bdy_file.empty()) {
                amrex::Error("NetCDF boundary file name must be provided via input");
            }

            bdy_time_interval = read_times_from_wrfbdy(nc_bdy_file,
                                                       bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                                       start_bdy_time, final_bdy_time);

            int ntimes_total = bdy_data_xlo.size();
            Vector<double> bdy_times(ntimes_total);

            // Initialize erfbdy file.
            if (write_erfbdy) {
                for (int itime = 0; itime < ntimes_total; ++itime) {
                    bdy_times[itime] = start_bdy_time + itime * bdy_time_interval;
                }

                InitERFBdyFile(erfbdy_file, ntimes_total, bdy_times,
                               geom[lev].Domain(), nvars_erfbdy, real_width);
                Print() << "Initialized erfbdy file: " << erfbdy_file << std::endl;
            }

            // *******************************************************************************************
            // We intentionally only read in the first three slices here ... we will read the rest in
            // as needed during the time stepping procedure
            // *******************************************************************************************
            int ntimes = bdy_data_xlo.size(); ntimes = amrex::min(ntimes, 3);
            Array<MultiFab*, AMREX_SPACEDIM> area_vec = {ax[lev].get(), ay[lev].get(), az[lev].get()};
            bool is_anelastic = (solverChoice.anelastic[0] == 1);
            for (int itime = 0; itime < ntimes; itime++)
            {
                read_and_convert_from_wrfbdy(itime, nc_bdy_file,
                                             bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                             wrf_MUB, wrf_C1H, wrf_C2H, wrf_RDNW, wrf_PHB, z_phys_cc[lev], z_phys_nd[lev],
                                             lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                                             r_hse, area_vec, geom[lev], use_moist, solverChoice.rebalance_wrf_input, domain_bcs_type,
                                             real_width, bdy_time_interval, is_anelastic);

                // Write this time to erfbdy.
                if (write_erfbdy) {
                    WriteERFBdyTimeSlice(erfbdy_file, itime,
                                         bdy_data_xlo[itime], bdy_data_xhi[itime],
                                         bdy_data_ylo[itime], bdy_data_yhi[itime],
                                         WRFBdyVars::NumTypes);
                    Print() << "Wrote erfbdy time index " << itime << " of " << ntimes_total-1 << std::endl;
                }
            } // itime

            // If writing erfbdy and we have more than 3 times, then process the remaining times.
            if (write_erfbdy && ntimes_total > 3) {
                Print() << "Processing remaining " << ntimes_total - 3 << " boundary times..." << std::endl;
                for (int itime = 3; itime < ntimes_total; ++itime) {
                    read_and_convert_from_wrfbdy(itime, nc_bdy_file,
                                                 bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                                 wrf_MUB, wrf_C1H, wrf_C2H, wrf_RDNW, wrf_PHB, z_phys_cc[lev], z_phys_nd[lev],
                                                 lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                                                 r_hse, area_vec, geom[lev], use_moist, solverChoice.rebalance_wrf_input, domain_bcs_type,
                                                 real_width, bdy_time_interval, is_anelastic);

                    WriteERFBdyTimeSlice(erfbdy_file, itime,
                                         bdy_data_xlo[itime], bdy_data_xhi[itime],
                                         bdy_data_ylo[itime], bdy_data_yhi[itime],
                                         WRFBdyVars::NumTypes);
                    Print() << "Wrote erfbdy time index " << itime << " of " << ntimes_total-1 << std::endl;

                    bdy_data_xlo[itime].clear();
                    bdy_data_xhi[itime].clear();
                    bdy_data_ylo[itime].clear();
                    bdy_data_yhi[itime].clear();
                }
                Print() << "Completed writing erfbdy times" << std::endl;
            } // itime
        } // use_erfbdy

        //
        // Start at the earliest time (read_from_wrfbdy)
        // Note we only have start_bdy_time if at level 0 and init_type == InitType::WRFInput or InitType::Metgrid
        //
        // Note that t_new and t_old carry *elapsed* time, not total time
        //
        if (lev == 0) {
            Print() << "start_bdy_time is " << std::setprecision(timeprecision) << start_bdy_time
                    << " from wrfbdy but note that time variable in simulation is elapsed time" << std::endl;
            t_new[lev] = zero;
            t_old[lev] = -bogus_large_value;
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
        low_time_interval = read_times_from_wrflow(nc_low_file, low_data_zlo, start_low_time, final_low_time);

        int ntimes = low_data_zlo.size();
        sst_lev[lev].resize(ntimes);
        tsk_lev[lev].resize(ntimes);

        // We can possibly run out of memory if we load all of wrfbdy and all of wrflow
        // Thus we only load the first two time slices here and load more only if needed
        ntimes = amrex::min(ntimes, 3);

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
 * @param subdomain        Box specifying the index space we are to initialize
 * @param l_rdOcp          Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param p_hse            MultiFab holding the hydrostatic base state pressure to be initialized
 * @param pi_hse           MultiFab holding the hydrostatic base state Exner pressure to be initialized
 * @param th_hse           MultiFab holding the hydrostatic base state potential temperature to be initialized
 * @param qv_hse           MultiFab holding the hydrostatic base state qv to be initialized
 * @param r_hse            MultiFab holding the hydrostatic base state density to be initialized
 * @param mf_PB            MultiFab holding WRF data specifying base state pressure
 * @param mf_ALB           Optional MultiFab holding inverse density perturbation data
 * @param z_phys_nd        Optional terrain nodal z-coordinate MultiFab
 * @param T00              Sea-level base-state temperature
 * @param P00              Sea-level base-state pressure
 * @param TLP              Base-state lapse rate
 * @param TISO             Isothermal stratosphere temperature
 * @param TLP_STRAT        Stratospheric lapse rate
 * @param P_STRAT          Pressure at the stratosphere transition
 */
void
init_base_state_from_wrfinput (const Box& subdomain,
                               const Real& l_rdOcp,
                               MultiFab& p_hse,
                               MultiFab& pi_hse,
                               MultiFab& th_hse,
                               MultiFab& qv_hse,
                               MultiFab& r_hse,
                               MultiFab& mf_PB,
                               MultiFab* mf_ALB,
                               MultiFab* z_phys_nd,
                               const Real& T00,
                               const Real& P00,
                               const Real& TLP,
                               const Real& TISO,
                               const Real& TLP_STRAT,
                               const Real& P_STRAT)
{
    const auto& dom_lo = lbound(subdomain);
    const auto& dom_hi = ubound(subdomain);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(p_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        Box gtbx = mfi.growntilebox();

        const Array4<Real      >&  p_hse_arr = p_hse.array(mfi);
        const Array4<Real      >& pi_hse_arr = pi_hse.array(mfi);
        const Array4<Real      >& th_hse_arr = th_hse.array(mfi);
        const Array4<Real      >& qv_hse_arr = qv_hse.array(mfi);
        const Array4<Real      >&  r_hse_arr = r_hse.array(mfi);

        const Array4<Real const>&      PB_arr = mf_PB.const_array(mfi);
        const Array4<Real const>&     ALB_arr = (mf_ALB) ? mf_ALB->const_array(mfi) :
                                                           Array4<const Real> {};

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Base state needs ghost cells filled, protect FAB access
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            int kk = std::max(k , dom_lo.z);
                kk = std::min(kk, dom_hi.z);

            Real Rd, Td, Thd;
            Real Pd = PB_arr(ii,jj,kk);
            // Have inverse base density
            if (ALB_arr) {
                Rd  = Real(1.0) / ALB_arr(ii,jj,kk);
                Td  = Pd / (R_d * Rd);
                Thd = getThgivenTandP(Td, Pd, l_rdOcp);
            } else {
                Td  = std::max(TISO, T00 + TLP * std::log(Pd/P00));
                if (P_STRAT > Real(0.) && Pd <= P_STRAT) {
                    Td = TISO + TLP_STRAT * std::log(Pd/P_STRAT);
                }
                Thd = getThgivenTandP(Td, Pd, l_rdOcp);
                Rd  = getRhogivenThetaPress (Thd, Pd, l_rdOcp);
            }

            // Fill HSE arrays (FOEXTRAP ghost cells)
             r_hse_arr(i,j,k) = Rd;
            th_hse_arr(i,j,k) = Thd;
            qv_hse_arr(i,j,k) = Real(0.);
             p_hse_arr(i,j,k) = Pd;
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
        });
    }

    // **************************************************************************
    // Rebalance the base state since state from WRFInput does not discretely
    // satisfy dp0/dz = -rho0 g
    // **************************************************************************
    int k_dom_lo = dom_lo.z;
    int k_dom_hi = dom_hi.z;

    // The vertical integration below is seeded with the surface values (P00,T00),
    // so these must be valid regardless of whether ALB was present in the file
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(P00) && (P00 > Real(0)) &&
                                     std::isfinite(T00) && (T00 > Real(0)),
                                     "Cannot rebalance the WRF base state: P00 and T00 must be positive");

#ifdef AMREX_USE_FLOAT
    Real tol  = Real(1.0e-6);
#else
    Real tol  = Real(1.0e-10);
#endif

    Real grav = CONST_GRAV;

    for (MFIter mfi(th_hse,TileNoZ()); mfi.isValid(); ++mfi)
    {
            Box bx  = mfi.tilebox();
            int klo = bx.smallEnd(2);
            int khi = bx.bigEnd(2);

            AMREX_ALWAYS_ASSERT((klo == k_dom_lo) && (khi == k_dom_hi));
            bx.makeSlab(2,klo);

            const Array4<Real>&  p_hse_arr = p_hse.array(mfi);
            const Array4<Real>& pi_hse_arr = pi_hse.array(mfi);
            const Array4<Real>& th_hse_arr = th_hse.array(mfi);
            const Array4<Real>&  r_hse_arr = r_hse.array(mfi);

            const Array4<const Real>& z_arr = z_phys_nd->const_array(mfi);

            ParallelFor(bx, [=,RdoCp_d=RdoCp]
                        AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
            {
                // integrate from surface to domain top
                Real dz, F, C;
                Real rho_tot_hi, rho_tot_lo;
                Real z_lo, z_hi;
                Real R_lo, R_hi;
                Real Th_lo, Th_hi;
                Real T_hi;
                Real P_lo, P_hi;

                Real qv_lo = zero;
                Real qv_hi = zero;

                // First integrate from surface to first CC at klo
                {
                    // Vertical grid spacing
                    z_lo = zero; // corresponding to p_0
                    z_hi = Real(0.125) * ( z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                         + z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(i,j+1,klo+1) + z_arr(i+1,j+1,klo+1) );

                    // dz == height of first cell center
                    dz = z_hi - z_lo;

                    // Known surface values
                    P_lo  = P00;
                    Th_lo = getThgivenTandP(T00, P00, RdoCp_d);
                    R_lo  = getRhogivenThetaPress(Th_lo, P_lo, RdoCp_d);
                    rho_tot_lo = R_lo;
                    C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

                    // Initial guess and residual
                    P_hi  = P_lo;
                    Th_hi = th_hse_arr(i,j,klo);
                    T_hi  = getTgivenPandTh(P_hi, Th_hi, RdoCp_d);
                    R_hi  = getRhogivenThetaPress(Th_hi, P_hi, RdoCp_d);
                    rho_tot_hi = R_hi;
                    F = P_hi + myhalf*rho_tot_hi*grav*dz + C;

                    // Do iterations
                    bool maintain_Th = true;
                    HSEutils::Newton_Raphson_hse(tol, RdoCp_d, dz,
                                                 grav, C, Th_hi, T_hi,
                                                 qv_hi, qv_hi,
                                                 P_hi, R_hi, F, maintain_Th);

                    // At first cell center
                     r_hse_arr(i,j,klo) = R_hi;
                     p_hse_arr(i,j,klo) = P_hi;
                    pi_hse_arr(i,j,klo) = getExnergivenP(p_hse_arr(i,j,klo), l_rdOcp);

                    P_lo = P_hi;
                    z_lo = z_hi;
                }

                for (int k(klo+1); k<=khi; ++k) {

                  z_hi = Real(0.125) * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                       +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1));
                  dz   = z_hi - z_lo;

                  // Establish known constant
                  Th_lo = th_hse_arr(i,j,k-1);
                  R_lo  = getRhogivenThetaPress(Th_lo, P_lo, RdoCp_d, qv_lo);
                  rho_tot_lo = R_lo;
                  C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

                  // Initial guess and residual
                  Th_hi = th_hse_arr(i,j,k);
                  T_hi  = getTgivenPandTh(P_hi, Th_hi, RdoCp_d);
                  R_hi  = getRhogivenThetaPress(Th_hi, P_hi, RdoCp_d, qv_hi);
                  rho_tot_hi = R_hi;
                  F = P_hi + myhalf*rho_tot_hi*grav*dz + C;

                  // Do iterations
                  bool maintain_Th = true;
                  HSEutils::Newton_Raphson_hse(tol, RdoCp_d, dz,
                                               grav, C, Th_hi, T_hi,
                                               qv_hi, qv_hi,
                                               P_hi, R_hi, F, maintain_Th);

                  // Assign data
                   r_hse_arr(i,j,k) = R_hi;
                   p_hse_arr(i,j,k) = P_hi;
                  pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);

                  P_lo = P_hi;
                  z_lo = z_hi;
                }
            });
    } // mfi
}

/**
 * Helper function for verifying the top boundary is valid and computing the bottom boundary.
 *
 * @param mf_PH  MultiFab storing WRF terrain coordinate data (PH)
 * @param mf_PHB MultiFab storing WRF terrain coordinate data (PHB)
 * @param domain Box holding index space of computational domain
 * @return Height assigned to the top of the ERF domain
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
    Gpu::HostVector  <Real> Max_h(3,-bogus_large_value);
    Gpu::DeviceVector<Real> Max_d(3);
    Gpu::copy(Gpu::hostToDevice, Max_h.begin(), Max_h.end(), Max_d.begin());

    Gpu::HostVector  <Real> Min_h(2, bogus_large_value);
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
            Real z_calc_lo = Real(0.25) * ( ph (ii,jj  ,klo) + ph (ii-1,jj  ,klo) +
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
            Real z_calc_hi = Real(0.25) * ( ph (ii,jj  ,khi) + ph (ii-1,jj  ,khi) +
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
            Real z_calc_hi = Real(0.25) * ( ph (ii,jj  ,khi-1) + ph (ii-1,jj  ,khi-1) +
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
    z_top = myhalf * (terrain_top_min + terrain_top_max);

    // If this creates a case where z_k < z_{k-1} then we do what we used to do
    if (terrain_km1_max > z_top) {
        amrex::Print() << "Max of second-to-highest row = " << terrain_km1_max <<
                          " which is greater than average of top row so defaulting to alternate approach " << std::endl;
        z_top = myhalf * (terrain_km1_max + terrain_top_max);
    }

    amrex::Print() << "Warning: ProbHi(2) will be ignored; we are setting top of domain to " << z_top << std::endl;

    return z_top;
}

/**
 * Helper function for initializing terrain coordinates from a WRF dataset.
 *
 * The unused level argument is retained for interface compatibility.
 *
 * @param z_top Height assigned to the top of the ERF domain
 * @param subdomain Box specifying the index space we are initializing
 * @param z_phys MultiFab specifying the node-centered z coordinates of the terrain
 * @param mf_PH MultiFab storing WRF perturbation geopotential data
 * @param mf_PHB MultiFab storing WRF base-state geopotential data
 */
void
init_terrain_from_wrfinput (int /*lev*/,
                            Geometry& geom,
                            const Real& z_top,
                            const Box& subdomain,
                            MultiFab* z_phys,
                            const MultiFab& mf_PH,
                            const MultiFab& mf_PHB,
                            Real& dz0_max,
                            const bool& use_wrf_height_grid)
{
    Print() << "Constructing nodal heights (z_phys_nd)" << std::endl;

    if (use_wrf_height_grid) {
        for ( MFIter mfi(*z_phys, false); mfi.isValid(); ++mfi )
        {
            Box gnbx = mfi.growntilebox();

            // This copies from NC_zphys on z-faces to z_phys_nd on nodes
            const Array4<Real      >&      z_arr = z_phys->array(mfi);
            const Array4<Real const>& nc_phb_arr = mf_PHB.const_array(mfi);
            const Array4<Real const>& nc_ph_arr  = mf_PH.const_array(mfi);

            // PHB and PH are on z-faces (half dx / half dy ahead of zphys)
            Box z_face_box = convert(subdomain,IntVect(0,0,1));

            // Prevent averaging from going into ghost cells
            int ilo = z_face_box.smallEnd()[0];
            int ihi = z_face_box.bigEnd()[0];
            int jlo = z_face_box.smallEnd()[1];
            int jhi = z_face_box.bigEnd()[1];
            int klo = z_face_box.smallEnd()[2];
            int khi = z_face_box.bigEnd()[2];

            ParallelFor(gnbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int ii = std::max(std::min(i,ihi),ilo);
                int jj = std::max(std::min(j,jhi),jlo);

                int im = std::max(std::min(i-1,ihi),ilo);
                int jm = std::max(std::min(j-1,jhi),jlo);

                if (k < klo) {
                    Real z_klo   = Real(0.25) * ( nc_ph_arr (ii,jj,klo  ) + nc_ph_arr (im,jj,klo  ) +
                                                  nc_ph_arr (ii,jm,klo  ) + nc_ph_arr (im,jm,klo) +
                                                  nc_phb_arr(ii,jj,klo  ) + nc_phb_arr(im,jj,klo  ) +
                                                  nc_phb_arr(ii,jm,klo  ) + nc_phb_arr(im,jm,klo) ) / CONST_GRAV;
                    Real z_klop1 = Real(0.25) * ( nc_ph_arr (ii,jj,klo+1) + nc_ph_arr (im,jj,klo+1) +
                                                  nc_ph_arr (ii,jm,klo+1) + nc_ph_arr (im,jm,klo+1) +
                                                  nc_phb_arr(ii,jj,klo+1) + nc_phb_arr(im,jj,klo+1) +
                                                  nc_phb_arr(ii,jm,klo+1) + nc_phb_arr(im,jm,klo+1) ) / CONST_GRAV;
                    z_arr(i, j, k) = two * z_klo - z_klop1;
                } else if (k > khi) {
                    Real z_khim1 = Real(0.25) * ( nc_ph_arr (ii,jj,khi-1) + nc_ph_arr (im,jj,khi-1) +
                                                  nc_ph_arr (ii,jm,khi-1) + nc_ph_arr (im,jm,khi-1) +
                                                  nc_phb_arr(ii,jj,khi-1) + nc_phb_arr(im,jj,khi-1) +
                                                  nc_phb_arr(ii,jm,khi-1) + nc_phb_arr(im,jm,khi-1) ) / CONST_GRAV;
                    z_arr(i, j, k) = two * z_top - z_khim1;
                } else if (k == khi) {
                    z_arr(i, j, k) = Real(0.25) * ( nc_ph_arr (ii,jj,k) + nc_ph_arr (im,jj,k) +
                                                    nc_ph_arr (ii,jm,k) + nc_ph_arr (im,jm,k) +
                                                    nc_phb_arr(ii,jj,k) + nc_phb_arr(im,jj,k) +
                                                    nc_phb_arr(ii,jm,k) + nc_phb_arr(im,jm,k) ) / CONST_GRAV;
                    z_arr(i, j, k) = z_top;
                } else {
                    // Note: wrfinput geopotentials ph, phb are only staggered in the vertical, i.e.,
                    //       they have dims (bottom_top_stag, south_north, west_east). On k==klo, we
                    //       will end up smoothing the terrain as we average from surface face centers
                    //       to nodes.
                    z_arr(i, j, k) = Real(0.25) * ( nc_ph_arr (ii,jj,k) + nc_ph_arr (im,jj,k) +
                                                    nc_ph_arr (ii,jm,k) + nc_ph_arr (im,jm,k) +
                                                    nc_phb_arr(ii,jj,k) + nc_phb_arr(im,jj,k) +
                                                    nc_phb_arr(ii,jm,k) + nc_phb_arr(im,jm,k) ) / CONST_GRAV;
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
    } else {

        // Lateral ghost cells
        IntVect ngz = z_phys->nGrowVect(); int kghost = ngz[2]; ngz[2] = 0;

        // PHB and PH are on z-faces
        Box z_face_dom = convert(subdomain,IntVect(0,0,1));

        // Z-face FAB for each z_wrf slice
        Box z_face_dom_slice = makeSlab(z_face_dom, 2, 0);
        FArrayBox z_slice_wrf(z_face_dom_slice, 1, The_Managed_Arena());
        FArrayBox z_slice_wrf_sfc(z_face_dom_slice, 1, The_Managed_Arena());

        // Z_phys is nodal
        Box node_dom = convert(subdomain, IntVect(1,1,1));
        int ilo = node_dom.smallEnd(0);
        int jlo = node_dom.smallEnd(1);
        int ihi = node_dom.bigEnd(0);
        int jhi = node_dom.bigEnd(1);
        int klo = node_dom.smallEnd(2);
        int khi = node_dom.bigEnd(2);

        // Process each slice
        int kstart = klo;
        int kend   = klo + 2;
        for (int k(kstart); k<kend; ++k) {
            z_slice_wrf.setVal<RunOn::Device>(zero);

            const Array4<Real>& z_slice_wrf_arr     = z_slice_wrf.array();
            const Array4<Real>& z_slice_wrf_sfc_arr = z_slice_wrf_sfc.array();

            // Fill the z-face fab with wrf heights
            for ( MFIter mfi(mf_PH); mfi.isValid(); ++mfi ) {

                Box vbx = mfi.validbox(); vbx.makeSlab(2,0);

                const Array4<const Real>& nc_phb_arr = mf_PHB.const_array(mfi);
                const Array4<const Real>& nc_ph_arr  = mf_PH.const_array(mfi);

                ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
                {
                    z_slice_wrf_arr(i,j,0) = ( nc_ph_arr(i,j,k) + nc_phb_arr(i,j,k) ) / CONST_GRAV;
                });
            }

            // Get global slice of WRF heights
            ParallelAllReduce::Sum(z_slice_wrf.dataPtr(),
                                   z_slice_wrf.size(),
                                   ParallelContext::CommunicatorAll());

            // Solve for node values that reproduce the WRF z-face values as
            // closely as a bounded, smooth nodal field can.  We deliberately do
            // *not* invert the four-node averaging operator exactly: its symbol
            // vanishes at the grid Nyquist mode, so exact de-averaging amplifies
            // grid-scale content of the WRF terrain without bound and returns
            // nodal heights that are kilometers away from the WRF terrain.
            // See ERF_NodalReconstruction.H for the regularized least-squares
            // formulation used instead.
            const Real tol = Real(1.e-10);
            NodalReconstruction NR_solver(z_face_dom_slice, geom);
            FArrayBox z_slice_ref = NR_solver.makeReference(z_slice_wrf);
            std::pair<amrex::FArrayBox,SolveInfo> result = NR_solver.solve(z_slice_wrf, z_slice_ref,
                                                                           VariationOperator::FirstDeriv, tol);
            const Array4<Real>& z_slice_erf_arr = result.first.array();
            const SolveInfo& info = result.second;
            if (!info.converged) {
                Print() << "WARNING: Nodal reconstruction did not converge at k = " << k
                        << "; residual is: " << info.final_residual
                        << " and requested tolerance was: " << tol << "\n";
            }

            // Range check on the reconstructed heights.  This is the check that
            // has teeth: comparing the four-node average against WRF (done at
            // the end of this routine) is close to satisfied by construction and
            // cannot detect a blown-up reconstruction.
            {
                Real wrf_min =  std::numeric_limits<Real>::max();
                Real wrf_max = -std::numeric_limits<Real>::max();
                LoopOnCpu(z_face_dom_slice, [=,&wrf_min,&wrf_max] (int i, int j, int /*k*/) noexcept
                {
                    wrf_min = amrex::min(wrf_min, z_slice_wrf_arr(i,j,0));
                    wrf_max = amrex::max(wrf_max, z_slice_wrf_arr(i,j,0));
                });

                Print() << "Nodal reconstruction at k = " << k << ": "
                        << info.iterations << " CG iterations, " << info.refinements
                        << " refinements, regularization " << info.regularization
                        << "\n    nodal heights in [" << info.min_value << ", " << info.max_value
                        << "] m vs WRF z-faces in [" << wrf_min << ", " << wrf_max << "] m"
                        << "\n    max |avg4(nodal) - WRF| = " << info.max_average_error
                        << " m (direct interpolation gives " << info.interp_average_error << " m)"
                        << "\n    max deviation from direct interpolation = " << info.deviation
                        << " m (cap " << info.deviation_cap << " m)" << std::endl;

                // Nodes may legitimately over/undershoot the cell values where
                // the terrain is under-resolved, but only by a fraction of the
                // relief of the layer itself.
                const Real relief = amrex::max(wrf_max - wrf_min, Real(1.0));
                const Real slack  = amrex::max(Real(0.5) * relief, Real(10.0));

                if ( !std::isfinite(info.min_value) || !std::isfinite(info.max_value) ||
                     (info.min_value < wrf_min - slack) || (info.max_value > wrf_max + slack) )
                {
                    Error("Nodal reconstruction produced heights far outside the range of the "
                          "WRF z-face heights; the reconstruction is not usable as terrain.");
                }
            }

            // Store the surface
            if (k==kstart) {
                LoopOnCpu(z_face_dom_slice, [=] (int i, int j, int /*k*/) noexcept
                {
                    z_slice_wrf_sfc_arr(i,j,0) = z_slice_wrf_arr(i,j,0);
                });
            }

            // Compute dz0_max from current layer and surface layer
            if (k==kstart+1) {
                dz0_max = std::numeric_limits<Real>::min();
                LoopOnCpu(z_face_dom_slice, [=,&dz0_max] (int i, int j, int /*k*/) noexcept
                {
                    Real dz0 = z_slice_wrf_arr(i,j,0) - z_slice_wrf_sfc_arr(i,j,0);
                    dz0_max = amrex::max(dz0, dz0_max);
                });
            }

            // Copy back to z_phys and handle all ghost cells
            for ( MFIter mfi(*z_phys); mfi.isValid(); ++mfi ) {
                Box gbx = mfi.growntilebox();
                Box sbx = makeSlab(gbx, 2, 0);
                const Array4<Real>& z_arr = z_phys->array(mfi);
                ParallelFor(sbx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
                {
                    int ii  = std::max(std::min(i,ihi),ilo);
                    int jj  = std::max(std::min(j,jhi),jlo);
                    z_arr(i,j,k) = z_slice_erf_arr(ii,jj,0);
                    if (k == klo + 1) {
                        Real dz = z_arr(i,j,k) - z_arr(i,j,k-1);
                        for (int lk(1); lk<(kghost+1); ++lk) {
                            z_arr(i,j,klo-lk) = z_arr(i,j,klo) - dz * static_cast<Real>(lk);
                        }
                    }
                });
            }
        } // k

        // Sanity check.
        //
        // The nodal heights are a regularized least-squares fit to the WRF
        // z-face heights, not an exact de-averaging of them.  The four-node average
        // therefore no longer reproduces WRF to round-off, and the size of the
        // mismatch is a genuine diagnostic of how much grid-scale terrain the
        // WRF field carries.  We report it, and we abort only on failures that
        // make the mesh unusable: non-finite heights, or a layer that is not
        // strictly increasing in z.
        Print() << "Verifying reconstructed nodal heights" << std::endl;
        {
            ReduceOps<ReduceOpMax, ReduceOpMin> reduce_op;
            ReduceData<Real, Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            for ( MFIter mfi(mf_PHB); mfi.isValid(); ++mfi ) {
                Box vbx = mfi.validbox();
                if (vbx.smallEnd(2) > klo+1) { continue; }
                vbx.setBig(2,klo+1);

                const Array4<const Real>& nc_phb_arr = mf_PHB.const_array(mfi);
                const Array4<const Real>& nc_ph_arr  = mf_PH.const_array(mfi);
                const Array4<const Real>& z_arr      = z_phys->const_array(mfi);

                reduce_op.eval(vbx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    Real z_face_wrf = ( nc_ph_arr(i,j,k) + nc_phb_arr(i,j,k) ) / CONST_GRAV;
                    Real z_face_erf = fourth * ( z_arr(i, j  , k) + z_arr(i+1, j  , k)
                                               + z_arr(i, j+1, k) + z_arr(i+1, j+1, k) );
                    Real avg_err = (k < khi) ? std::fabs(z_face_erf-z_face_wrf) : zero;

                    // Nodal layer thickness; only defined once we are at klo
                    Real dz = (k == klo) ? (z_arr(i,j,klo+1) - z_arr(i,j,klo))
                                         : std::numeric_limits<Real>::max();
                    return {avg_err, dz};
                });
            } // mfi

            ReduceTuple hv = reduce_data.value(reduce_op);
            Real max_avg_err = amrex::max(amrex::get<0>(hv), zero);
            Real min_dz      = amrex::get<1>(hv);
            ParallelAllReduce::Max(max_avg_err, ParallelContext::CommunicatorAll());
            ParallelAllReduce::Min(min_dz     , ParallelContext::CommunicatorAll());

            Print() << "Max |avg4(nodal z) - WRF z-face| over the reconstructed levels: "
                    << max_avg_err << " m" << std::endl;
            Print() << "Min nodal thickness of the first layer: " << min_dz << " m" << std::endl;

            if (!std::isfinite(max_avg_err) || !std::isfinite(min_dz)) {
                Error("Non-finite nodal heights produced by the terrain reconstruction");
            }
            if (min_dz <= zero) {
                Error("Reconstructed nodal terrain gives a non-positive first layer thickness; "
                      "the WRF terrain is too rough to represent on the ERF nodal mesh. "
                      "Consider running with erf.use_wrf_height_grid = true.");
            }
        }
    } // use wrf grid
}
#endif // ERF_USE_NETCDF
