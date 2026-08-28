/**
 * \file ERF_InitFromWRFInput.cpp
 */

#include "ERF.H"
#include "ERF_EOS.H"
#include "ERF_Constants.H"
#include "ERF_Utils.H"
#include "ERF_ProbCommon.H"
#include "ERF_DataStruct.H"

#include "ERF_ReadFromWRFInput.H"
#include "ERF_ReadFromWRFBdy.H"
#include "ERF_WriteERFBdy.H"
#include "ERF_ReadFromERFBdy.H"

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
 * @param[in] avg_grid_faces_to_nodes Whether to average the wrfinput heights onto the nodes
 *                                   rather than reconstructing nodal heights from them.
 */
void
init_terrain_from_wrfinput (int lev,
                            Geometry& geom,
                            const Real& z_top,
                            const Box& subdomain,
                            MultiFab* z_phys_nd,
                            const MultiFab& NC_PH_fab,
                            const MultiFab& NC_PHB_fab,
                            Real& dz0_max,
                            const bool& avg_grid_faces_to_nodes);

/**
 * Initialize hydrostatic base state data from a WRF dataset.
 *
 * The profile is built analytically from the six WRF reference-state parameters
 * and the ERF cell-centered heights; PB and ALB from the file are not used.
 *
 * @param[in] subdomain Box specifying the index space to initialize.
 * @param[in] l_rdOcp Constant $R_d/c_p$ (currently unused; the constexpr RdoCp is used instead).
 * @param[out] p_hse MultiFab holding the hydrostatic base state pressure.
 * @param[out] pi_hse MultiFab holding the hydrostatic base state Exner pressure.
 * @param[out] th_hse MultiFab holding the hydrostatic base state potential temperature.
 * @param[out] qv_hse MultiFab holding the hydrostatic base state qv.
 * @param[out] r_hse MultiFab holding the hydrostatic base state density.
 * @param[in] mf_PB MultiFab holding WRF data specifying base state pressure (currently unused).
 * @param[in] mf_ALB MultiFab holding inverse density perturbation data (currently unused).
 * @param[in] z_phys_cc Cell-centered z-coordinate MultiFab; required, must not be null.
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
                               MultiFab* z_phys_cc,
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
        // variables but leave them zero-filled.  T00 and P00 are an absolute
        // temperature and pressure, and TLP is the lapse rate that the closed-form
        // p(z) inversion in init_base_state_from_wrfinput divides by, so a
        // non-positive value of any of the three leaves us with no usable
        // reference profile: TLP == 0 makes the inversion NaN and TLP < 0 flips the
        // sign of the quadratic root it takes.  If any of them is bad we discard
        // the whole group and keep the defaults.  (TLP_STRAT and P_STRAT may
        // legitimately be zero -- that just disables the stratospheric layer -- so
        // they are only reset alongside a bad T00/P00/TLP.)
        const bool params_ok = std::isfinite(T00)  && (T00  > Real(0)) &&
                               std::isfinite(P00)  && (P00  > Real(0)) &&
                               std::isfinite(TLP)  && (TLP  > Real(0)) &&
                               std::isfinite(TISO) && std::isfinite(TLP_STRAT) &&
                               std::isfinite(P_STRAT);

        if (!params_ok) {
            Print() << "WARNING: WRF base state parameters read from " << fname
                    << " are invalid: (T00, P00, TLP, TISO, TLP_STRAT, P_STRAT) = ("
                    << T00 << ", " << P00 << ", " << TLP << ", " << TISO << ", "
                    << TLP_STRAT << ", " << P_STRAT << ")\n";
            Print() << "         T00, P00 and TLP must all be positive and finite; "
                       "reverting to ERF defaults.\n";
            Print() << "         NOTE: the base state is built entirely from these six "
                       "parameters -- PB and ALB in the file are not used -- so the "
                       "resulting base state is synthetic and may be inconsistent with "
                       "the state read from " << fname << ".\n";
            T00 = T00_def;   P00 = P00_def;   TLP = TLP_def;
            TISO = TISO_def; TLP_STRAT = TLP_STRAT_def; P_STRAT = P_STRAT_def;
        }

        // WRF evaluates the reference temperature as max(TISO, T00 + TLP*ln(p/P00)),
        // so iso_temp == 0 there simply means "no isothermal cap".  ERF instead
        // inverts each layer analytically and divides by TISO inside the isothermal
        // layer, so a non-positive TISO would be a division by zero (and, taken
        // literally, a profile running down to 0 K).  The tropospheric parameters
        // this file supplies are still good, so keep them and cap with the ERF
        // default rather than discarding the whole group.
        if (params_ok && (TISO <= Real(0))) {
            Print() << "WARNING: TISO read from " << fname << " is " << TISO
                    << "; using the ERF default isothermal cap of " << TISO_def
                    << " K and keeping T00, P00 and TLP from the file.\n";
            TISO = TISO_def;
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
        NC_names.push_back("QVAPOR");
        NC_names.push_back("QCLOUD");
        if (solverChoice.moisture_indices.qi >= 0) {
            NC_names.push_back("QICE");
        }
        if (solverChoice.moisture_indices.qr >= 0) {
            NC_names.push_back("QRAIN");   // 26
        }
        if (solverChoice.moisture_indices.qs >= 0) {
            NC_names.push_back("QSNOW");
        }
        if (solverChoice.moisture_indices.qg >= 0) {
            NC_names.push_back("QGRAUP");
        }
    }
    NC_names.push_back("IVGTYP");
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
                (var_name == "U")       || (var_name == "V")       || (var_name == "W")      ||
                (var_name == "THM")     || (var_name == "QVAPOR")  || (var_name == "QCLOUD") ||
                (var_name == "QICE")    || (var_name == "QRAIN")   || (var_name == "QSNOW")  ||
                (var_name == "QGRAUP")  || (var_name == "PH")      || (var_name == "PHB");
            const bool required_hydrometeor = solverChoice.use_wrf_bdy_qc_qi &&
                ((var_name == "QCLOUD" && solverChoice.moisture_indices.qc >= 0) ||
                 (var_name == "QICE"   && solverChoice.moisture_indices.qi >= 0) ||
                 (var_name == "QRAIN"  && solverChoice.moisture_indices.qr >= 0) ||
                 (var_name == "QSNOW"  && solverChoice.moisture_indices.qs >= 0) ||
                 (var_name == "QGRAUP" && solverChoice.moisture_indices.qg >= 0));
            if (!success && required_hydrometeor) {
                amrex::Abort(std::string("erf.use_wrf_bdy_qc_qi requires " + var_name +
                                         " in wrfinput for the active moisture component").c_str());
            }
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

                // XLONG_U and XLAT_V are edge-staggered in the file (west_east_stag /
                // south_north_stag, so nx+1 / ny+1 entries) but their fabs carry CELL
                // index type, so the typed subdomain -- and hence the intersection copy
                // below -- would drop the last staggered column/row. Keep it here so the
                // ghost fill further down can pick up the true east/north edge: a coupled
                // ocean model consumes lon_m/lat_m as a corner mesh through
                // ERF::GetOceanToAtmosCornerCoordinates, and duplicating the neighbour
                // instead collapses the outermost corner quads to zero area.
                // NOTE: var_fab keeps CELL index type; only its extent is widened.
                if (var_name == "XLONG_U" &&
                    var_fab_from_file.box().bigEnd(0) > subdomain_to_fill_typed.bigEnd(0)) {
                    subdomain_to_fill_typed.growHi(0,1);
                }
                if (var_name == "XLAT_V" &&
                    var_fab_from_file.box().bigEnd(1) > subdomain_to_fill_typed.bigEnd(1)) {
                    subdomain_to_fill_typed.growHi(1,1);
                }

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

                    // Invert exactly the region that received data from the file: the copy
                    // of "ALB" and the plus of "AL" both act on cons_fab.box() & var_fab.box(),
                    // which includes the ghost cells lying inside the region read from the file.
                    Box vbx = cons_fab.box() & var_fab.box();

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

                    // See the note in the ALB/AL branch above: invert exactly the region
                    // the copy below writes, so the in-domain ghost cells hold density
                    // rather than specific volume, and the cells still holding zero are skipped.
                    Box vbx = cons_fab.box() & var_fab.box();

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
            if ( var_name == "THM"     ||
                 var_name == "QVAPOR"  ||
                 var_name == "QCLOUD"  ||
                 var_name == "QICE"    ||
                 var_name == "QRAIN"   ||
                 var_name == "QSNOW"   ||
                 var_name == "QGRAUP" )
            {
                AMREX_ALWAYS_ASSERT(micro->Get_Qstate_NonMoist_Size() == 0);

                int icomp = -1;
                if (var_name == "THM") {
                    icomp    = RhoTheta_comp;
                } else if (var_name == "QVAPOR") {
                    icomp    = solverChoice.moisture_indices.qv;
                } else if (var_name == "QCLOUD") {
                    icomp    = solverChoice.moisture_indices.qc;
                } else if (var_name == "QICE") {
                    icomp    = solverChoice.moisture_indices.qi;
                } else if (var_name == "QRAIN") {
                    icomp    = solverChoice.moisture_indices.qr;
                } else if (var_name == "QSNOW") {
                    icomp    = solverChoice.moisture_indices.qs;
                } else if (var_name == "QGRAUP") {
                    icomp    = solverChoice.moisture_indices.qg;
                }
                // Note: RhoQ7-RhoQ9 (nc, nn, nr for WDM6) or RhoQ7-RhoQ11 (nc, ni, nr, ns, ng for Morrison)
                // start at zero and are diagnosed/initialized by the microphysics scheme

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
                    if (icomp >= 0 && icomp < lev_new[Vars::cons].nComp()) {
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
              // var_fab retains XLAT_V's staggered row at j = ny (see the growHi above),
              // so clamp lat_m against var_fab's own extent to give the j = ny ghost the
              // true north edge instead of a copy of row ny-1.
              // sinPhi_m/cosPhi_m deliberately stay on the cell-domain clamp: they are
              // cell-centred Coriolis factors whose ghosts are read at the hi domain
              // faces by ERF_MakeMomSources.cpp, and this fix is not meant to move the
              // Coriolis source. So sin_arr/cos_arr do not track lat_m in that one row.
              int vf_j_hi = var_fab.box().bigEnd(1);
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
                      int sj = amrex::min(amrex::max(j, j_lo), vf_j_hi);
                      dst_arr(i,j,0) = src_arr(li,sj,0);

                      Real lat_rad = src_arr(li,lj,0) * (PI/Real(180.));
                      sin_arr(i,j,0) = std::sin(lat_rad);
                      cos_arr(i,j,0) = std::cos(lat_rad);
                  });
              }
          }

          // Initialize Longitude
          if ( var_name == "XLONG_U" ) {
              // var_fab retains XLONG_U's staggered column at i = nx (see the growHi
              // above), so clamp lon_m against var_fab's own extent to give the i = nx
              // ghost the true east edge instead of a copy of column nx-1.
              int vf_i_hi = var_fab.box().bigEnd(0);
              lon_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
              for ( MFIter mfi(*(lon_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                  Box gtbx = mfi.growntilebox();
                  const Array4<      Real>& dst_arr = (lon_m[lev])->array(mfi);
                  const Array4<const Real>& src_arr = var_fab.const_array();
                  ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                  {
                      int li = amrex::min(amrex::max(i, i_lo), vf_i_hi);
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

            // Default to uniform grid or solve for a stretched grid
            if (dz0_max >= z_top/Nz) {
                SFact   = one;
                dz0_max = z_top/Nz;
            } else {
                int max_iter = 50;
                int iter     = 0;
                Real F       = dz0_max * ( (std::pow(SFact,Nz) - one) / (SFact - one) ) - z_top;
                while (std::fabs(F)>tol && iter<max_iter) {
                    Real dFdSF = dz0_max * ( Nz * std::pow(SFact,Nz-one) * (SFact - one)
                                           - std::pow(SFact,Nz) + one ) /
                                           std::pow(SFact-one,two);
                    SFact     -= F/dFdSF;
                    SFact      = std::max(one+tol,SFact);
                    F          = dz0_max * ( (std::pow(SFact,Nz) - one) / (SFact - one) ) - z_top;
                    ++iter;
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::fabs(F) <= tol,
                                                 "Newton iterations to determine the grid stretching factor failed!\n");
            }

            // Build the zlevels
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

            // Update stretched dz and build terrain fitted coords
            update_stretched_dz(lev, zlevels_stag, stretched_dz_h, stretched_dz_d);
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
        //
        // NOTE: W is intentionally *not* remapped below, unlike the cell-centered state
        //       and the horizontal velocities.  ERF initializes from wrfinput files
        //       produced by real.exe / ideal.exe (see Docs/sphinx_doc/Initialization.rst),
        //       and those preprocessors leave W identically zero -- they carry no vertical
        //       velocity from the driving analysis.  If a file
        //       carrying a non-zero W were ever fed in here (a WRF history or restart file
        //       rather than a preprocessor wrfinput), W would need a zvel_tmp remap onto
        //       the z-faces exactly as xvel and yvel get below, or it would be left at
        //       WRF's z-face heights while everything else moved to ERF's.
        // **************************************************************************
#ifdef AMREX_DEBUG
        // Hold the assumption stated in the note above to account.  A non-zero W here
        // means the input is not a real.exe/ideal.exe wrfinput, and silently skipping
        // its remap would leave the velocity field inconsistent at initialization.
        // W is copied verbatim from the file and only ever divided by the map factor,
        // so a genuinely zero W stays exactly zero and the exact test is safe.
        {
            const Real w_norm = lev_new[Vars::zvel].norm0();
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(w_norm == Real(0),
                "init_from_wrfinput: the input file carries a non-zero W, which is not "
                "remapped onto the ERF grid (see issue 3655).  A zvel_tmp remap onto the "
                "z-faces must be added before such a file can be used here.");
        }
#endif

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

        MultiFab theta(rho.boxArray(), rho.DistributionMap(), 1, 0);
        MultiFab::Copy(theta, lev_new[Vars::cons], RhoTheta_comp, 0, 1, 0);
        MultiFab::Divide(theta, lev_new[Vars::cons], Rho_comp , 0, 1, 0);

        MultiFab qv(rho.boxArray(), rho.DistributionMap(), 1, 0);
        MultiFab::Copy(qv, lev_new[Vars::cons], RhoQ1_comp, 0, 1, 0);
        MultiFab::Divide(qv, lev_new[Vars::cons], Rho_comp , 0, 1, 0);

        MultiFab qt(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, 0);
        int n_qstate_into_total = micro->Get_Qstate_Moist_Size() - micro->Get_Qstate_Moist_NumConc_Size();
        make_qt(lev_new[Vars::cons], qt, n_qstate_into_total);

        bool maintain_Th = true;
        rebalance_columns(rho, theta, qv, qt, z_phys_nd[lev].get(), geom[lev], maintain_Th);

        // Update (rho qv) in the state
        MultiFab::Multiply(qv, rho, 0, 0, 1, 0);
        MultiFab::Copy(lev_new[Vars::cons], qv, 0, RhoQ1_comp, 1, 0);

        // Update (rho theta) in the state
        MultiFab::Multiply(theta, rho, 0, 0, 1, 0);
        MultiFab::Copy(lev_new[Vars::cons], theta, 0, RhoTheta_comp, 1, 0);
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
                                  mf_PB, mf_ALB.get(), z_phys_cc[lev].get(),
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
        const bool separate_hydrometeors = solverChoice.use_wrf_bdy_qc_qi &&
            wrf_bdy_has_separate_hydrometeors(solverChoice.moisture_indices);
        if (write_erfbdy) {
            nvars_erfbdy = separate_hydrometeors
                         ? WRFBdyHydrometeorVars::NumTypes
                         : (solverChoice.use_wrf_bdy_qc_qi
                            ? WRFBdyVars::NumTypes : WRFBdyVars::LegacyNumTypes);
        }
        auto repack_runtime_bdy = [&] (const int itime) {
            repack_wrfbdy_to_realbdy(bdy_data_xlo[itime], solverChoice.use_wrf_bdy_qc_qi, separate_hydrometeors);
            repack_wrfbdy_to_realbdy(bdy_data_xhi[itime], solverChoice.use_wrf_bdy_qc_qi, separate_hydrometeors);
            repack_wrfbdy_to_realbdy(bdy_data_ylo[itime], solverChoice.use_wrf_bdy_qc_qi, separate_hydrometeors);
            repack_wrfbdy_to_realbdy(bdy_data_yhi[itime], solverChoice.use_wrf_bdy_qc_qi, separate_hydrometeors);
        };

        // Path 1: Load from existing erfbdy file.
        if (use_erfbdy) {
            Print() << "Loading boundary data from erfbdy file: " << erfbdy_file << std::endl;

            // Read metadata and times from erfbdy.
            int ntimes_erfbdy;
            Vector<double> bdy_times;
            bdy_time_interval = read_times_from_erfbdy(erfbdy_file,
                                                       ntimes_erfbdy, nvars_erfbdy, real_width,
                                                       bdy_times, start_bdy_time, final_bdy_time);

            if (nvars_erfbdy != WRFBdyVars::LegacyNumTypes &&
                nvars_erfbdy != WRFBdyVars::NumTypes &&
                nvars_erfbdy != WRFBdyHydrometeorVars::NumTypes) {
                amrex::Error("ERFBdy cache has an unsupported boundary-variable layout");
            }
            const int expected_bdy_nvars = separate_hydrometeors
                ? WRFBdyHydrometeorVars::NumTypes
                : (solverChoice.use_wrf_bdy_qc_qi ? WRFBdyVars::NumTypes
                                                  : WRFBdyVars::LegacyNumTypes);
            if (nvars_erfbdy != expected_bdy_nvars) {
                amrex::Error("ERFBdy cache layout does not match the active WRF hydrometeor boundary mode; regenerate it from wrfbdy");
            }

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
                repack_runtime_bdy(itime);
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
                                             r_hse, area_vec, geom[lev], use_moist,
                                             solverChoice.use_wrf_bdy_qc_qi,
                                             solverChoice.moisture_indices.qi >= 0,
                                             separate_hydrometeors,
                                             solverChoice.rebalance_wrf_input, domain_bcs_type,
                                             real_width, bdy_time_interval, is_anelastic);

                // Write this time to erfbdy.
                if (write_erfbdy) {
                    WriteERFBdyTimeSlice(erfbdy_file, itime,
                                         bdy_data_xlo[itime], bdy_data_xhi[itime],
                                         bdy_data_ylo[itime], bdy_data_yhi[itime],
                                         nvars_erfbdy);
                    Print() << "Wrote erfbdy time index " << itime << " of " << ntimes_total-1 << std::endl;
                }
                if (itime == ntimes_total-1 && itime > 0) {
                    repack_runtime_bdy(itime-1);
                }
                repack_runtime_bdy(itime);
            } // itime

            // If writing erfbdy and we have more than 3 times, then process the remaining times.
            if (write_erfbdy && ntimes_total > 3) {
                Print() << "Processing remaining " << ntimes_total - 3 << " boundary times..." << std::endl;
                for (int itime = 3; itime < ntimes_total; ++itime) {
                    read_and_convert_from_wrfbdy(itime, nc_bdy_file,
                                                 bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                                 wrf_MUB, wrf_C1H, wrf_C2H, wrf_RDNW, wrf_PHB, z_phys_cc[lev], z_phys_nd[lev],
                                                 lev_new[Vars::xvel], lev_new[Vars::yvel], lev_new[Vars::cons],
                                                 r_hse, area_vec, geom[lev], use_moist,
                                                 solverChoice.use_wrf_bdy_qc_qi,
                                                 solverChoice.moisture_indices.qi >= 0,
                                                 separate_hydrometeors,
                                                 solverChoice.rebalance_wrf_input, domain_bcs_type,
                                                 real_width, bdy_time_interval, is_anelastic);

                    WriteERFBdyTimeSlice(erfbdy_file, itime,
                                         bdy_data_xlo[itime], bdy_data_xhi[itime],
                                         bdy_data_ylo[itime], bdy_data_yhi[itime],
                                         nvars_erfbdy);
                    Print() << "Wrote erfbdy time index " << itime << " of " << ntimes_total-1 << std::endl;

                    if (itime == ntimes_total-1 && itime > 0) {
                        repack_runtime_bdy(itime-1);
                    }
                    repack_runtime_bdy(itime);

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
 * The profile is built analytically from the six WRF reference-state parameters
 * and the ERF cell-centered heights; PB and ALB from the file are not used.
 *
 * @param subdomain        Box specifying the index space we are to initialize
 * @param l_rdOcp          Real constant specifying Rydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$); currently unused, the constexpr RdoCp is used instead
 * @param p_hse            MultiFab holding the hydrostatic base state pressure to be initialized
 * @param pi_hse           MultiFab holding the hydrostatic base state Exner pressure to be initialized
 * @param th_hse           MultiFab holding the hydrostatic base state potential temperature to be initialized
 * @param qv_hse           MultiFab holding the hydrostatic base state qv to be initialized
 * @param r_hse            MultiFab holding the hydrostatic base state density to be initialized
 * @param mf_PB            MultiFab holding WRF data specifying base state pressure; currently unused
 * @param mf_ALB           MultiFab holding inverse density perturbation data; currently unused
 * @param z_phys_cc        Cell-centered z-coordinate MultiFab; required, must not be null
 * @param T00              Sea-level base-state temperature
 * @param P00              Sea-level base-state pressure
 * @param TLP              Base-state lapse rate
 * @param TISO             Isothermal stratosphere temperature
 * @param TLP_STRAT        Stratospheric lapse rate
 * @param P_STRAT          Pressure at the stratosphere transition
 */
void
init_base_state_from_wrfinput (const Box& subdomain,
                               const Real& /*l_rdOcp*/,
                               MultiFab& p_hse,
                               MultiFab& pi_hse,
                               MultiFab& th_hse,
                               MultiFab& qv_hse,
                               MultiFab& r_hse,
                               MultiFab& /*mf_PB*/,
                               MultiFab* /*mf_ALB*/,
                               MultiFab* z_phys_cc,
                               const Real& T00,
                               const Real& P00,
                               const Real& TLP,
                               const Real& TISO,
                               const Real& TLP_STRAT,
                               const Real& P_STRAT)
{
    const auto& dom_lo = lbound(subdomain);
    const auto& dom_hi = ubound(subdomain);

    // The analytic inversion below is a function of the true cell-centered
    // height, so z_phys_cc is required here -- it is dereferenced unconditionally
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(z_phys_cc != nullptr,
                                     "init_base_state_from_wrfinput requires z_phys_cc");

    // **************************************************************************
    // The WRF reference state is piecewise in log-pressure:
    //   (1) troposphere      T = T00  + TLP       * ln(p/P00)         p >  P_iso
    //   (2) isothermal layer T = TISO                       P_STRAT < p <= P_iso
    //   (3) stratosphere     T = TISO + TLP_STRAT * ln(p/P_STRAT)     p <= P_STRAT
    // Each piece inverts to p(z) in closed form under dp/dz = -rho g; here we
    // precompute the two interface heights so the inversion can branch on z.
    // **************************************************************************
    const Real x_iso = (TISO - T00) / TLP;
    const Real P_iso = P00 * std::exp(x_iso);
    const Real z_iso = -(R_d/CONST_GRAV) * (T00*x_iso + myhalf*TLP*x_iso*x_iso);

    // The upper stratospheric layer is optional (P_STRAT == 0 or TLP_STRAT == 0
    // disables it) and is only meaningful if it begins above the isothermal layer,
    // i.e. if P_STRAT is below the pressure at which the isothermal layer starts.
    const bool want_strat = (P_STRAT > zero) && (TLP_STRAT != zero);
    const bool use_strat  = want_strat && (P_STRAT < P_iso);
    const Real z_strat    = (use_strat) ? z_iso + (R_d*TISO/CONST_GRAV)*std::log(P_iso/P_STRAT)
                                        : z_iso;

    // A configured stratospheric layer that lies at or below the isothermal
    // transition cannot be represented, so say so rather than dropping it quietly.
    if (want_strat && !use_strat) {
        Print() << "WARNING: the WRF stratospheric layer is being ignored: P_STRAT = "
                << P_STRAT << " Pa is not below the pressure at the base of the "
                << "isothermal layer, P_iso = " << P_iso << " Pa.\n";
        Print() << "         TLP_STRAT = " << TLP_STRAT << " will have no effect and "
                << "the atmosphere above z_iso will be isothermal at TISO = "
                << TISO << " K.\n";
    }

    Print() << "WRF base state layer interfaces: z_iso = " << z_iso << " m";
    if (use_strat) Print() << ", z_strat = " << z_strat << " m";
    Print() << "\n";

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(p_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // The base state must be valid in its ghost cells as well: physbcs_base is
        // not applied until after init_from_wrfinput returns, but r_hse is consumed
        // before that (read_and_convert_from_wrfbdy -> scale_bdy_normal_by_rho0
        // reads rho0 one cell outside the domain).  Leaving the ghost cells at the
        // setVal(0) from ERF_MakeNewArrays would divide by zero there, so fill the
        // grown box here just as the pre-analytic version of this routine did.
        Box gtbx = mfi.growntilebox();

        // z_phys_cc carries fewer ghost cells than the base state, so clamp the
        // height lookup into the region where it is defined; that gives the outer
        // ghost cells a zeroth-order extrapolation of the profile.
        const Box  zbx  = amrex::grow(mfi.validbox(), z_phys_cc->nGrowVect());
        const auto z_lo = lbound(zbx);
        const auto z_hi = ubound(zbx);

        const Array4<Real      >&  p_hse_arr =  p_hse.array(mfi);
        const Array4<Real      >& th_hse_arr = th_hse.array(mfi);
        const Array4<Real      >& qv_hse_arr = qv_hse.array(mfi);
        const Array4<Real      >&  r_hse_arr =  r_hse.array(mfi);
        const Array4<Real      >& pi_hse_arr = pi_hse.array(mfi);

        const Array4<const Real>& z_cc_arr   = z_phys_cc->const_array(mfi);

        ParallelFor(gtbx, [=,zero_d=zero,RdoCp_d=RdoCp]
                    AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const int ii = amrex::min(amrex::max(i, z_lo.x), z_hi.x);
            const int jj = amrex::min(amrex::max(j, z_lo.y), z_hi.y);
            const int kk = amrex::min(amrex::max(k, z_lo.z), z_hi.z);

            // Analytical inversion with true CC heights, branching on the layer
            Real Pd, Td;
            const Real z = z_cc_arr(ii,jj,kk);
            if (z <= z_iso) {
                // Troposphere: z = -(R_d/g) * (T00*x + TLP*x^2/2), x = ln(p/P00)
                const Real ToA  = T00 / TLP;
                const Real disc = amrex::max(ToA*ToA - two*CONST_GRAV*z/(TLP*R_d), zero_d);
                Pd = P00 * std::exp(-ToA + std::sqrt(disc));
                Td = T00 + TLP * std::log(Pd/P00);
            }
            else if (!use_strat || z <= z_strat) {
                // Isothermal layer: exponential decay with scale height R_d*TISO/g
                Pd = P_iso * std::exp(-CONST_GRAV*(z - z_iso)/(R_d*TISO));
                Td = TISO;
            }
            else {
                // Upper stratosphere. Same quadratic as the troposphere with
                // (TISO, TLP_STRAT, P_STRAT, z_strat) in place of (T00, TLP, P00, 0).
                // NOTE: TLP_STRAT is negative, so the root must NOT be folded as
                // sqrt(X)/TLP_STRAT -> sqrt(X/TLP_STRAT^2); that drops the sign.
                const Real disc = amrex::max(TISO*TISO
                                  - two*TLP_STRAT*CONST_GRAV*(z - z_strat)/R_d, zero_d);
                Pd = P_STRAT * std::exp((-TISO + std::sqrt(disc))/TLP_STRAT);
                Td = TISO + TLP_STRAT * std::log(Pd/P_STRAT);
            }

            // Fill HSE arrays for balancing
             p_hse_arr(i,j,k) = Pd;
            th_hse_arr(i,j,k) = getThgivenTandP(Td, Pd, RdoCp_d);
            qv_hse_arr(i,j,k) = zero;
             r_hse_arr(i,j,k) = getRhogivenThetaPress(th_hse_arr(i,j,k), Pd, RdoCp_d);
            pi_hse_arr(i,j,k) = getExnergivenP(Pd, RdoCp_d);
        });
    }

    // **************************************************************************
    // Rebalance the base state since state from WRFInput since it does not
    // discretely satisfy dp0/dz = -rho0 g on the ERF grid
    // **************************************************************************
    int k_dom_lo = dom_lo.z;
    int k_dom_hi = dom_hi.z;

    // The vertical integration below is seeded with the analytic profile in the
    // lowest cell of each column, z_cc(i,j,klo) and p_hse(i,j,klo), not with
    // (P00,T00); read_base_state_params_from_wrfinput has already guaranteed that
    // the parameters those were built from are positive and finite.

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

            const Array4<Real>&  r_hse_arr = r_hse.array(mfi);
            const Array4<Real>&  p_hse_arr = p_hse.array(mfi);
            const Array4<Real>& pi_hse_arr = pi_hse.array(mfi);

            const Array4<const Real>& th_hse_arr = th_hse.const_array(mfi);
            const Array4<const Real>& z_cc_arr   = z_phys_cc->const_array(mfi);

            ParallelFor(bx, [=,RdoCp_d=RdoCp]
                        AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
            {
                // Integrate from surface to domain top
                Real T_hi;
                Real z_lo, z_hi;
                Real R_lo, R_hi;
                Real Th_lo, Th_hi;
                Real P_lo, P_hi;
                Real rho_tot_hi, rho_tot_lo;

                Real dz, F, C;

                Real qv_lo = zero;
                Real qv_hi = zero;

                z_lo =  z_cc_arr(i,j,klo);
                P_lo = p_hse_arr(i,j,klo);
                for (int k(klo+1); k<=khi; ++k) {
                    z_hi = z_cc_arr(i,j,k);
                    dz   = z_hi - z_lo;

                    // Establish known constant
                    Th_lo = th_hse_arr(i,j,k-1);
                    R_lo  = getRhogivenThetaPress(Th_lo, P_lo, RdoCp_d, qv_lo);
                    rho_tot_lo = R_lo;
                    C  = -P_lo + myhalf*rho_tot_lo*grav*dz;

                    // Initial guess and residual
                    P_hi  =  p_hse_arr(i,j,k);
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
                    pi_hse_arr(i,j,k) = getExnergivenP(P_hi, RdoCp_d);

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
    // Reductions for the bottom/top boundary (in that order)
    //
    // NOTE: These must use the ReduceOps machinery rather than hand-rolled
    //       Gpu::Atomic::Min/Max calls.  The Gpu::Atomic operations are *not*
    //       atomic when running on the host, so they lose updates as soon as
    //       this loop is threaded with OpenMP.  ReduceData holds one
    //       accumulator per thread as long as it is constructed (and read)
    //       outside of the OpenMP parallel region, as it is here.
    //
    ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op_bot;
    ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op_top;
    ReduceOps<ReduceOpMax>              reduce_op_km1;

    ReduceData<Real, Real> reduce_data_bot(reduce_op_bot);
    ReduceData<Real, Real> reduce_data_top(reduce_op_top);
    ReduceData<Real>       reduce_data_km1(reduce_op_km1);

    using ReduceTupleMinMax = typename decltype(reduce_data_bot)::Type;
    using ReduceTupleMax    = typename decltype(reduce_data_km1)::Type;

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
        //
        // NOTE: The slabs below must be built from the *tilebox* so that each
        //       (i,j) is visited exactly once.  Using the validbox here made
        //       every tile of a box revisit the entire slab.
        //
        // NOTE: The clamping bounds must still come from the validbox since
        //       they exist to keep the stencil inside this box's data.
        //
        Box vbx = mfi.validbox();
        Box tbx = mfi.tilebox();

        Box nodal_box = amrex::surroundingNodes(vbx);
        int ilo = nodal_box.smallEnd()[0];
        int ihi = nodal_box.bigEnd()[0];
        int jlo = nodal_box.smallEnd()[1];
        int jhi = nodal_box.bigEnd()[1];

        // For the top boundary
        Box Fab2dBox_hi, Fab2dBox_hi_m1;
        if (tbx.bigEnd(2) == khi) {
            Fab2dBox_hi    = makeSlab(tbx,2,khi  );
            Fab2dBox_hi_m1 = makeSlab(tbx,2,khi-1);
        }

        // For the bottom boundary
        Box Fab2dBox_lo;
        if (tbx.smallEnd(2) == klo) {
            Fab2dBox_lo = makeSlab(tbx,2,klo);
        }

        auto const& phb = mf_PHB.const_array(mfi);
        auto const& ph  = mf_PH.const_array(mfi);

        //
        // This loop computes the min and max values of the bottom surface
        //
        if (Fab2dBox_lo.ok()) {
            reduce_op_bot.eval(Fab2dBox_lo, reduce_data_bot,
            [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept -> ReduceTupleMinMax
            {
                int ii = std::max(std::min(i,ihi-1),ilo+1);
                int jj = std::max(std::min(j,jhi-1),jlo+1);
                Real z_calc_lo = Real(0.25) * ( ph (ii,jj  ,klo) + ph (ii-1,jj  ,klo) +
                                                ph (ii,jj-1,klo) + ph (ii-1,jj-1,klo) +
                                                phb(ii,jj  ,klo) + phb(ii-1,jj  ,klo) +
                                                phb(ii,jj-1,klo) + phb(ii-1,jj-1,klo) ) / CONST_GRAV;
                return {z_calc_lo, z_calc_lo};
            });
        }

        //
        // This loop computes the max value of the top surface
        //
        if (Fab2dBox_hi.ok()) {
            reduce_op_top.eval(Fab2dBox_hi, reduce_data_top,
            [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept -> ReduceTupleMinMax
            {
                int ii = std::max(std::min(i,ihi-1),ilo+1);
                int jj = std::max(std::min(j,jhi-1),jlo+1);
                Real z_calc_hi = Real(0.25) * ( ph (ii,jj  ,khi) + ph (ii-1,jj  ,khi) +
                                                ph (ii,jj-1,khi) + ph (ii-1,jj-1,khi) +
                                                phb(ii,jj  ,khi) + phb(ii-1,jj  ,khi) +
                                                phb(ii,jj-1,khi) + phb(ii-1,jj-1,khi) ) / CONST_GRAV;
                return {z_calc_hi, z_calc_hi};
            });
        }

        //
        // This loop computes the max value of the layer just below the top surface
        //
        if (Fab2dBox_hi_m1.ok()) {
            reduce_op_km1.eval(Fab2dBox_hi_m1, reduce_data_km1,
            [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept -> ReduceTupleMax
            {
                int ii = std::max(std::min(i,ihi-1),ilo+1);
                int jj = std::max(std::min(j,jhi-1),jlo+1);
                Real z_calc_hi = Real(0.25) * ( ph (ii,jj  ,khi-1) + ph (ii-1,jj  ,khi-1) +
                                                ph (ii,jj-1,khi-1) + ph (ii-1,jj-1,khi-1) +
                                                phb(ii,jj  ,khi-1) + phb(ii-1,jj  ,khi-1) +
                                                phb(ii,jj-1,khi-1) + phb(ii-1,jj-1,khi-1) ) / CONST_GRAV;
                return {z_calc_hi};
            });
        }
    } // mfi

    ReduceTupleMinMax hv_bot = reduce_data_bot.value(reduce_op_bot);
    ReduceTupleMinMax hv_top = reduce_data_top.value(reduce_op_top);
    ReduceTupleMax    hv_km1 = reduce_data_km1.value(reduce_op_km1);

    Real terrain_bottom_min = amrex::get<0>(hv_bot);
    Real terrain_bottom_max = amrex::get<1>(hv_bot);
    Real terrain_top_min    = amrex::get<0>(hv_top);
    Real terrain_top_max    = amrex::get<1>(hv_top);
    Real terrain_km1_max    = amrex::get<0>(hv_km1);

    ParallelDescriptor::ReduceRealMin(terrain_bottom_min);
    ParallelDescriptor::ReduceRealMin(terrain_top_min);

    ParallelDescriptor::ReduceRealMax(terrain_bottom_max);
    ParallelDescriptor::ReduceRealMax(terrain_top_max);
    ParallelDescriptor::ReduceRealMax(terrain_km1_max);

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
                            MultiFab* z_phys_nd,
                            const MultiFab& mf_PH,
                            const MultiFab& mf_PHB,
                            Real& dz0_max,
                            const bool& avg_grid_faces_to_nodes)
{
    Print() << "Constructing nodal heights (z_phys_nd)" << std::endl;

    if (avg_grid_faces_to_nodes) {
        for ( MFIter mfi(*z_phys_nd); mfi.isValid(); ++mfi )
        {
            Box gtbx = mfi.growntilebox();

            // This copies from NC_zphys on z-faces to z_phys_nd on nodes
            const Array4<Real      >&      z_arr = z_phys_nd->array(mfi);
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

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
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
                    // Extrapolate linearly below the surface -- note that z_phys_nd
                    // has more than one ghost node in the vertical, so this must
                    // depend on k rather than filling every ghost node with one value
                    z_arr(i, j, k) = z_klo - static_cast<Real>(klo-k) * (z_klop1 - z_klo);
                } else if (k > khi) {
                    Real z_khim1 = Real(0.25) * ( nc_ph_arr (ii,jj,khi-1) + nc_ph_arr (im,jj,khi-1) +
                                                  nc_ph_arr (ii,jm,khi-1) + nc_ph_arr (im,jm,khi-1) +
                                                  nc_phb_arr(ii,jj,khi-1) + nc_phb_arr(im,jj,khi-1) +
                                                  nc_phb_arr(ii,jm,khi-1) + nc_phb_arr(im,jm,khi-1) ) / CONST_GRAV;
                    // Extrapolate linearly above the top of the domain -- note that
                    // z_phys_nd has more than one ghost node in the vertical, so this
                    // must depend on k rather than filling every ghost node with one value
                    z_arr(i, j, k) = z_top + static_cast<Real>(k-khi) * (z_top - z_khim1);
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
        IntVect ngz = z_phys_nd->nGrowVect(); int kghost = ngz[2]; ngz[2] = 0;

        // PHB and PH are on z-faces
        Box z_face_dom = convert(subdomain,IntVect(0,0,1));

        // Z-face FAB for each z_wrf slice
        Box z_face_dom_slice = makeSlab(z_face_dom, 2, 0);
        FArrayBox z_slice_wrf(z_face_dom_slice, 1, The_Managed_Arena());
        FArrayBox z_slice_wrf_sfc(z_face_dom_slice, 1, The_Managed_Arena());

        // Z_phys is nodal
        Box node_dom = convert(subdomain, IntVect(1,1,1));
        int klo = node_dom.smallEnd(2);
        int khi = node_dom.bigEnd(2);

        // Process the surface and the first level above it
        int kstart = klo;
        int kend   = klo + 2;
        for (int k(kstart); k<kend; ++k) {
            z_slice_wrf.setVal<RunOn::Device>(zero);

            const Array4<Real>& z_slice_wrf_arr     = z_slice_wrf.array();
            const Array4<Real>& z_slice_wrf_sfc_arr = z_slice_wrf_sfc.array();

            // Fill the z-face fab with wrf heights -- each rank fills only the
            // boxes it owns, so the AllReduce below gathers the global slice
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
            Gpu::streamSynchronize();
            ParallelAllReduce::Sum(z_slice_wrf.dataPtr(),
                                   z_slice_wrf.size(),
                                   ParallelContext::CommunicatorAll());

            // Solve for the nodal heights of this level
            FArrayBox z_slice_erf = reconstruct_nodal_height_slice(z_face_dom_slice, geom,
                                                                   z_slice_wrf, k, "WRF z-faces");

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

            // Copy back to z_phys, filling the lateral ghost nodes
            fill_nodal_level_from_slice(*z_phys_nd, k, z_slice_erf);
        } // k

        // Extrapolate linearly into the ghost nodes below the surface
        for ( MFIter mfi(*z_phys_nd); mfi.isValid(); ++mfi ) {
            Box gbx = mfi.growntilebox();
            if (klo < gbx.smallEnd(2) || klo+1 > gbx.bigEnd(2)) { continue; }

            const Array4<Real>& z_arr = z_phys_nd->array(mfi);
            ParallelFor(makeSlab(gbx,2,klo), [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
            {
                Real dz = z_arr(i,j,klo+1) - z_arr(i,j,klo);
                for (int lk(1); lk<(kghost+1); ++lk) {
                    z_arr(i,j,klo-lk) = z_arr(i,j,klo) - dz * static_cast<Real>(lk);
                }
            });
        }

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
                const Array4<const Real>& z_arr      = z_phys_nd->const_array(mfi);

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
                      "Consider running with erf.avg_grid_faces_to_nodes = true.");
            }
        }
    } // avg_grid_faces_to_nodes
}
#endif // ERF_USE_NETCDF
