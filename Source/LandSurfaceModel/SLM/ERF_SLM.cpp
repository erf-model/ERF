#include <ERF_SLM.H>
#include "ERF_EOS.H"
#include "ERF_TileNoZ.H"
#include "ERF_MicrophysicsUtils.H"
#include <AMReX_PlotFileUtil.H>
#include "ERF.H"

using namespace amrex;

/* Initialize lsm data structures */
void
SLM::Init (const int& /*lev*/,
           const MultiFab& cons_in,
           const MultiFab& u_in,
           const MultiFab& v_in,
           const Geometry& geom,
           const Geometry& /*geom0*/,
           Vector<BCRec>& /*domain_bcs_type*/,
           IntVect& /*refRatio*/,
           const Real& dt,
           std::unique_ptr<amrex::MultiFab>& z_phys_nd_in,
           Vector<Vector<std::string>>& /*nc_init_file*/)
{
    m_dt = dt;
    m_geom = geom;
    z_phys_nd = z_phys_nd_in.get();

    ParmParse pp("slm");
    pp.query("nsoil", m_nz_lsm);
    pp.queryarr("soil_dz", m_dz_lsm);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      m_dz_lsm.size() == m_nz_lsm,
      "Provided soil thicknesses most match number of soil layers");
    AMREX_ALWAYS_ASSERT(m_dz_lsm.size() > 0);

    Box domain = geom.Domain();
    khi_lsm    = domain.smallEnd(2) - 1; // index of z_r
    klo_lsm    = khi_lsm - m_nz_lsm + 1;

    LsmDataMap.resize(m_lsm_data_size);
    LsmDataMap = {
      LsmVar_SLM::tsurf,         LsmVar_SLM::ustar,        LsmVar_SLM::tstar,       LsmVar_SLM::qstar,
      LsmVar_SLM::tv,            LsmVar_SLM::mv,           LsmVar_SLM::soilt,       LsmVar_SLM::soilw,
      LsmVar_SLM::sand,          LsmVar_SLM::clay,         LsmVar_SLM::s_depth,
      LsmVar_SLM::flbu,          LsmVar_SLM::flbv,         LsmVar_SLM::flbq,
      LsmVar_SLM::flbt,          LsmVar_SLM::prsfc,        LsmVar_SLM::precipref,
      LsmVar_SLM::swdsvisxyref,  LsmVar_SLM::swdsnirxyref, LsmVar_SLM::swdsvisdxyref,
      LsmVar_SLM::swdsnirdxyref, LsmVar_SLM::lwref,        LsmVar_SLM::coszrsxy, LsmVar_SLM::tref,
      LsmVar_SLM::uref,          LsmVar_SLM::vref,         LsmVar_SLM::dref,
      LsmVar_SLM::qref,          LsmVar_SLM::pref,         LsmVar_SLM::node_z,
      LsmVar_SLM::soilt_nudge,   LsmVar_SLM::soilw_nudge,  LsmVar_SLM::lai,
      LsmVar_SLM::vegtype,       LsmVar_SLM::soiltype,     LsmVar_SLM::veg_frac, LsmVar_SLM::veg_frac_min, LsmVar_SLM::veg_frac_max,
      LsmVar_SLM::emis_sfc,      LsmVar_SLM::alb_nir_sfc,  LsmVar_SLM::alb_vis_sfc,
      LsmVar_SLM::alb_nir_sfc_diff, LsmVar_SLM::alb_vis_sfc_diff,
      LsmVar_SLM::soil_transp_frac};

    LsmDataName.resize(m_lsm_data_size);
    LsmDataName = {"tsurf",         "ustar",          "tstar",
                  "qstar",         "tveg",           "mv",
                  "tsoil",         "wsoil",          "sand",
                  "clay",          "soil_thickness", "surface_u",
                  "surface_v",     "surface_vapor",  "surface_heat",
                  "precip_soil",   "ref_precip",     "SW_dw_dir_vis",
                  "SW_dw_dir_nir", "SW_dw_dif_vis",  "SW_dw_dif_nir",
                  "LW_dw",         "cos_zenith",     "ref_t",          "ref_u",
                  "ref_v",         "ref_d",          "ref_q",
                  "ref_p",         "node_z",         "soilt_nudge",
                  "soilw_nudge",   "lai",            "vegtype",
                  "soiltype",      "veg_frac", "veg_frac_min", "veg_frac_max", "emis_sfc",
                  "alb_nir_sfc",   "alb_vis_sfc", "alb_nir_sfc_diff", "alb_vis_sfc_diff",
                  "soil_transp_frac"};

    AMREX_ALWAYS_ASSERT(LsmDataMap.size() == LsmDataName.size());
    AMREX_ALWAYS_ASSERT(LsmDataMap.size() == m_lsm_data_size);

    LsmFluxMap.resize(m_lsm_flux_size);
    LsmFluxMap = {LsmFlux_SLM::t_flux, LsmFlux_SLM::q_flux, LsmFlux_SLM::tau13, LsmFlux_SLM::tau23, LsmFlux_SLM::olen};

    LsmFluxName.resize(m_lsm_flux_size);
    LsmFluxName = {"t_flux", "q_flux", "tau13", "tau23", "olen"};

    AMREX_ALWAYS_ASSERT(LsmFluxMap.size() == LsmFluxName.size());
    AMREX_ALWAYS_ASSERT(LsmFluxMap.size() == m_lsm_flux_size);

    // NOTE: All boxes in ba extend from zlo to zhi, so this transform is valid.
    //       If that were to change, the dm and new ba are no longer valid and
    //       direct copying between lsm data/flux vars cannot be done in a parfor.

    // Set box array for lsm data
    IntVect ng(1,1,1);
    BoxArray ba = cons_in.boxArray();
    DistributionMapping dm = cons_in.DistributionMap();
    BoxList bl_lsm = ba.boxList();
    for (auto& b : bl_lsm) {
        b.setBig(2, khi_lsm);                  // First point below the surface
        b.setSmall(2, klo_lsm);                // Last point below the surface
    }
    BoxArray ba_lsm(std::move(bl_lsm));

    // Set up lsm geometry
    const RealBox& dom_rb = m_geom.ProbDomain();
    const Real*    dom_dx = m_geom.CellSize();
    RealBox lsm_rb = dom_rb;
    Real lsm_z_hi = dom_rb.lo(2); // z_r
    Real lsm_z_lo = lsm_z_hi;
    for (int k = 0; k < m_nz_lsm; k++)
    {
        lsm_z_lo -= m_dz_lsm[k];
    }
    lsm_rb.setHi(2,lsm_z_hi); lsm_rb.setLo(2,lsm_z_lo);

    amrex::Box lsm_dom = m_geom.Domain();
    lsm_dom.setSmall(2, klo_lsm);
    lsm_dom.setBig(2, khi_lsm);
    m_lsm_geom.define(lsm_dom, lsm_rb, m_geom.Coord(), m_geom.isPeriodic());

    BoxList bl_lsm_2d = ba_lsm.boxList();
    for (auto& b : bl_lsm_2d) {
        b.setRange(2, 0, 1);
    }
    ba_lsm_2d = BoxArray(std::move(bl_lsm_2d));
    IntVect ng_2d(0, 0, 0);

    // Create the data
    for (auto ivar = 0; ivar < LsmVar_SLM::NumVars; ++ivar) {
        // State vars are CC
        lsm_fab_vars[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_vars[ivar]->setVal(0.0);
    }

    // Create the fluxes
    for (auto ivar = 0; ivar < LsmFlux_SLM::NumVars; ++ivar) {
        // NOTE: Fluxes are CC with ghost cells for averaging
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(ba_lsm_2d, dm, 1, IntVect(1,1,0));
        lsm_fab_flux[ivar]->setVal(0.0);
    }

    // build list for checkpointing extra variables not mapped to ERF in lsm_fab_vars
    for (int i = 0; i < LsmVar_SLM::NumVars; i++) {
        int found = 0;
        for (int j = 0; j < m_lsm_data_size; j++) {
            if (LsmDataMap[j] == i) {
                found = 1;
                break;
            }
        }
        if (!found) {
            unmapped_fields.push_back(i);
        }
    }

    // Initial olen to neutral condition
    lsm_fab_flux[LsmFlux_SLM::olen]->setVal(1.0E34);

    // packed temporary 1D arrays in soil water and soil temperature
    soilt_vars.define(ba_lsm, dm, SLM_DST::NumVars, 0);
    soilw_vars.define(ba_lsm, dm, SLM_DSW::NumVars, 0);
    soilt_vars.setVal(0.0);
    soilw_vars.setVal(0.0);

    // Create local 2D data
    // TODO: Placeholder landmask array - fix!
    landmask.define(ba_lsm_2d, dm, 1, ng_2d);
    landmask.setVal(1);

    landtype.define(ba_lsm_2d, dm, 1, ng_2d);
    nroot.define(ba_lsm_2d, dm, 1, ng_2d);
    LAI.define(ba_lsm_2d, dm, 1, ng_2d);
    SAI.define(ba_lsm_2d, dm, 1, ng_2d);
    sstxy.define(ba_lsm_2d, dm, 1, ng_2d);

    t_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    mw.define(ba_lsm_2d, dm, 1, ng_2d);
    mws.define(ba_lsm_2d, dm, 1, ng_2d);
    t_skin.define(ba_lsm_2d, dm, 1, ng_2d);
    t_ground_skin.define(ba_lsm_2d, dm, 1, ng_2d);
    t_cas.define(ba_lsm_2d, dm, 1, ng_2d);
    q_cas.define(ba_lsm_2d, dm, 1, ng_2d);
    
    t_sfc.define(ba_lsm_2d, dm, 1, ng_2d);
    q_sfc.define(ba_lsm_2d, dm, 1, ng_2d);
    
    q_gr.define(ba_lsm_2d, dm, 1, ng_2d);
    sdew.define(ba_lsm_2d, dm, 1, ng_2d);

    // NOAHMP radiation state variables
    albold_noahmp.define(ba_lsm_2d, dm, 1, ng_2d);
    tauss_noahmp.define(ba_lsm_2d, dm, 1, ng_2d);

    vegetype.define(ba_lsm_2d, dm, 1, ng_2d);
    vege_YES.define(ba_lsm_2d, dm, 1, ng_2d);
    cp_vege.define(ba_lsm_2d, dm, 1, ng_2d);
    cbiom.define(ba_lsm_2d, dm, 1, ng_2d);
    dleaf.define(ba_lsm_2d, dm, 1, ng_2d);
    z0_sfc.define(ba_lsm_2d, dm, 1, ng_2d);
    Khai_L.define(ba_lsm_2d, dm, 1, ng_2d);
    phi_1.define(ba_lsm_2d, dm, 1, ng_2d);
    phi_2.define(ba_lsm_2d, dm, 1, ng_2d);
    IR_emis_vege.define(ba_lsm_2d, dm, 1, ng_2d);
    IR_emis_soil.define(ba_lsm_2d, dm, 1, ng_2d);
    IR_emis_grnd.define(ba_lsm_2d, dm, 1, ng_2d);
    ztop.define(ba_lsm_2d, dm, 1, ng_2d);
    disp_hgt.define(ba_lsm_2d, dm, 1, ng_2d);
    Rgl.define(ba_lsm_2d, dm, 1, ng_2d);
    Rc_min.define(ba_lsm_2d, dm, 1, ng_2d);
    hs_rc.define(ba_lsm_2d, dm, 1, ng_2d);
    rootL.define(ba_lsm_2d, dm, 1, ng_2d);
    root_a.define(ba_lsm_2d, dm, 1, ng_2d);
    root_b.define(ba_lsm_2d, dm, 1, ng_2d);
    precip_extinc.define(ba_lsm_2d, dm, 1, ng_2d);
    mw_mx.define(ba_lsm_2d, dm, 1, ng_2d);
    mws_mx.define(ba_lsm_2d, dm, 1, ng_2d);
    BAI.define(ba_lsm_2d, dm, 1, ng_2d);
    IMPERV.define(ba_lsm_2d, dm, 1, ng_2d);

    cp_vege.setVal(0.0);
    cbiom.setVal(0.02);
    dleaf.setVal(0.04);
    z0_sfc.setVal(0.0);
    Khai_L.setVal(0.0);
    phi_1.setVal(0.0);
    phi_2.setVal(0.0);
    ztop.setVal(0.0);
    disp_hgt.setVal(0.0);
    Rgl.setVal(0.0);
    Rc_min.setVal(0.0);
    hs_rc.setVal(0.0);
    rootL.setVal(0.0);
    root_a.setVal(0.0);
    root_b.setVal(0.0);
    precip_extinc.setVal(0.0);
    mw_mx.setVal(0.0);
    mws_mx.setVal(0.0);
    BAI.setVal(0.0);
    IMPERV.setVal(0.0);

    mw_inc.define(ba_lsm_2d, dm, 1, ng_2d);

    evapo_dry.define(ba_lsm_2d, dm, 1, ng_2d);

    shf_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    shf_soil.define(ba_lsm_2d, dm, 1, ng_2d);
    shf_air.define(ba_lsm_2d, dm, 1, ng_2d);
    lhf_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    lhf_soil.define(ba_lsm_2d, dm, 1, ng_2d);
    lhf_air.define(ba_lsm_2d, dm, 1, ng_2d);
    evp_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    evp_soil.define(ba_lsm_2d, dm, 1, ng_2d);
    evp_air.define(ba_lsm_2d, dm, 1, ng_2d);

    albedovis_v.define(ba_lsm_2d, dm, 1, ng_2d);
    albedonir_v.define(ba_lsm_2d, dm, 1, ng_2d);
    albedovis_s.define(ba_lsm_2d, dm, 1, ng_2d);
    albedonir_s.define(ba_lsm_2d, dm, 1, ng_2d);
    albedovis_v.setVal(0.0);
    albedovis_s.setVal(0.0);
    albedonir_v.setVal(0.0);
    albedonir_s.setVal(0.0);
    IR_emis_vege.setVal(0.0);
    IR_emis_soil.setVal(0.98);
    IR_emis_grnd.setVal(0.0);
    vege_YES.setVal(0.0);

    r_a.define(ba_lsm_2d, dm, 1, ng_2d);
    r_b.define(ba_lsm_2d, dm, 1, ng_2d);
    r_c.define(ba_lsm_2d, dm, 1, ng_2d);
    r_d.define(ba_lsm_2d, dm, 1, ng_2d);
    r_soil.define(ba_lsm_2d, dm, 1, ng_2d);

    net_rad.define(ba_lsm_2d, dm, SLM_NetRad::NumVars, ng_2d);
    wet_canop.define(ba_lsm_2d, dm, 1, ng_2d);

    slm_diag.define(ba_lsm_2d, dm, SLM_Diag::NumVars, ng_2d);
    slm_diag.setVal(0.0);

    zrefxy.define(ba_lsm_2d, dm, 1, ng_2d);
    zrefxy.setVal(0.0);

    r_soil.setVal(0.0);
    lhf_air.setVal(0.0);
    lhf_canop.setVal(0.0);
    shf_air.setVal(0.0);
    shf_canop.setVal(0.0);
    evp_air.setVal(0.0);
    evp_canop.setVal(0.0);
    evp_soil.setVal(0.0);

    t_canop.setVal(0.0);
    t_cas.setVal(0.0);
    t_skin.setVal(0.0);
    t_ground_skin.setVal(0.0);
    q_cas.setVal(0.0);
    
    t_sfc.setVal(0.0);
    q_sfc.setVal(0.0);
    q_gr.setVal(0.0);
    sdew.setVal(0.0);
    wet_canop.setVal(0.0);

    mw.setVal(0.0);
    mws.setVal(0.0);
    mw_inc.setVal(0.0);
    evapo_dry.setVal(0.0);

    r_a.setVal(0.0);
    r_b.setVal(0.0);
    r_c.setVal(0.0);
    r_d.setVal(0.0);
    r_soil.setVal(0.0);
    lsm_fab_vars[LsmVar_SLM::ustar]->setVal(0.1);
    lsm_fab_vars[LsmVar_SLM::tstar]->setVal(0.0);
    lsm_fab_vars[LsmVar_SLM::qstar]->setVal(0.0);

    // Initialize SLM from inputs if specified
    init_from_file();

    landtype.setVal(landtype0);
    nroot.setVal(m_nz_lsm);
    LAI.setVal(LAI0);
    SAI.setVal(0.0);
    sstxy.setVal(0.0);
    mws_mx.setVal(mws_mx0);

    net_rad.setVal(0.0);

    // Initialize NOAHMP radiation state variables
    albold_noahmp.setVal(0.65);  // initial snow albedo
    tauss_noahmp.setVal(0.0);    // initial snow age

    // Set a default vegetation fraction of 1.0.
    // Cells without vegetation (baresoil, urban) are set to 0.0 later on
    // If using WRFInput, this fraction will be overwritten with those values
    lsm_fab_vars[LsmVar_SLM::veg_frac]->setVal(1.0);

    if (!use_wrfinput) {
        // Initialize here if not using wrfinput, otherwise it is done at first time step
        slm_init();
    }

    //Following Noah-MP, zref is modified to become ztop (canopy topheight) + dz0(center height of the atmosphere's lowest grid)
    //First, read in SLM_use_inputs, if true, zref is set from the inputs parameter slm.zref
    // if false, SLM is coupled to ERF atmosphere, therefore, zref is computed  
    pp.query("SLM_use_inputs",set_from_file);

    if (!set_from_file && !use_wrfinput) {
        Real zlo      = m_geom.ProbLo(2);
        Real dz       = m_geom.CellSize(2);
        for ( MFIter mfi(cons_in,TileNoZ()); mfi.isValid(); ++mfi) {
            const Box& xybx      = mfi.growntilebox(0);
            const Array4<const Real>& z_nd_arr = (use_terrain) ? z_phys_nd->const_array(mfi) : Array4<Real>{};
            auto ztop_arr = ztop.array(mfi);
            auto zrefxy_arr = zrefxy.array(mfi);

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                amrex::Real znd = (z_nd_arr) ? Compute_Zrel_AtCellCenter(i, j, 0, z_nd_arr) : zlo + 0.5*dz;
                zrefxy_arr(i,j,0) = ztop_arr(i,j,0) + znd;
            }); 
        }       
    } else {
        zrefxy.setVal(zref);  // set zrefxy with the value read from inputs file
    }

    //slm_to_rad_vars = {Lsm_Data_Ptr(LsmVar_SLM::tsurf), &albedovis_s, &albedovis_v, &albedonir_s, &albedonir_v, &IR_emis_vege, &net_rad};
    slm_to_rad_vars = {Lsm_Data_Ptr(LsmVar_SLM::tsurf),
                       Lsm_Data_Ptr(Lsm_DataIndex("emis_sfc")),
                       Lsm_Data_Ptr(Lsm_DataIndex("alb_vis_sfc")),
                       Lsm_Data_Ptr(Lsm_DataIndex("alb_nir_sfc")),
                       Lsm_Data_Ptr(Lsm_DataIndex("alb_vis_sfc_diff")),
                       Lsm_Data_Ptr(Lsm_DataIndex("alb_nir_sfc_diff"))};
}
/**
 * Initialize SLM from input data - used for testing only
 */
void SLM::init_from_file()
{
    ParmParse pp("slm");
    pp.query("SLM_use_inputs", set_from_file);

    pp.query("landtype0", landtype0);
    pp.query("LAI0", LAI0);

    auto const get_layer_prop = [&pp](std::string name, const int nz, amrex::Vector<amrex::Real> &prop)
    {
        int nval = pp.countval(name.c_str());

        if (nval == 1)
        {
            // read single val and fill prop array over all soil layers
            amrex::Real tmp;
            pp.query(name.c_str(), tmp);
            prop.resize(nz);
            std::fill(prop.begin(), prop.end(), tmp);
        } else {
            pp.queryarr(name.c_str(), prop);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(prop.size() == nz, " Expected " + name + " to have " + std::to_string(nz) + " values!");
        }
    };

    pp.query("soiltnudging", dosoiltnudging);
    pp.query("soilwnudging", dosoilwnudging);
    pp.query("tausoil", tausoil);

    if (!use_wrfinput) {
        get_layer_prop("clay0", m_nz_lsm, clay0);
        get_layer_prop("sand0", m_nz_lsm, sand0);
        get_layer_prop("sw0", m_nz_lsm, sw0);
        get_layer_prop("st0", m_nz_lsm, st0);
    }
    if (dosoiltnudging || dosoilwnudging) {
        get_layer_prop("relax_hgt", m_nz_lsm, relax_hgt);
    } else {
        relax_hgt.resize(m_nz_lsm);
        std::fill(relax_hgt.begin(), relax_hgt.end(), 0.0);
    }

    if (!use_wrfinput) {
        for (int i = 0; i < m_nz_lsm; i++)
        {
            amrex::Print() << " " << i << ": soilt = " << st0[i] << " soilw = " << sw0[i] << " clay = " << clay0[i] << " sand = " << sand0[i] << " relax = " << relax_hgt[i] << std::endl;
        }
    }

    pp.query("tabs_s", tabs_s);
    pp.query("t00", t00);

    pp.query("z0_soil", z0_soil);
    pp.query("mws_mx0", mws_mx0);
    pp.query("Rc_max", Rc_max);
    pp.query("T_opt", T_opt);
    pp.query("zref", zref);

    // Read NoahmpTable.TBL
    pp.query("use_parameter_file", use_param_file);
    pp.query("interpolate_lai", interpolate_lai);
    pp.query("parameter_file", parameter_file);
    pp.query("veg_dataset", veg_dataset);
    pp.query("soil_dataset", soil_dataset);
    pp.query("use_param_tbl", use_wrf_lai);

    std::string radiation_scheme_name = "NoahMP";
    pp.query("radiation_scheme", radiation_scheme_name);
    const std::string radiation_scheme_lower = amrex::toLower(radiation_scheme_name);
    if (radiation_scheme_lower == "slm") {
        radiation_scheme = RadiationScheme::SLM;
        radiation_scheme_name = "SLM";
    } else if (radiation_scheme_lower == "noahmp") {
        radiation_scheme = RadiationScheme::NoahMP;
        radiation_scheme_name = "NoahMP";
        if (!use_param_file) {
            amrex::Abort("slm.radiation_scheme=noahmp requires "
                         "slm.use_parameter_file=true");
        }
    } else {
        amrex::Abort("Invalid slm.radiation_scheme='" + radiation_scheme_lower +
                     "'. Expected 'SLM' or 'NoahMP'");
    }
    amrex::Print() << " SLM radiation scheme: " << radiation_scheme_name << std::endl;

    if (use_param_file) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(use_wrf_lai != use_param_file, "Cannot use parameter file and parameter table, must choose one method");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(parameter_file != "", "Using parameter file, but file path is empty!");
    }

    // Read NoahMP table if provided
    if (use_param_file && parameter_file != "") {
        full_params = ReadParameterFile(parameter_file, param_veg_categories, param_soil_categories);

        std::string soil_param_key, veg_param_key;
        if (soil_dataset == "stas_ruc") {
            soil_param_key = "noahmp_soil_stas_ruc_parameters";
        } else if (soil_dataset == "stas") {
            soil_param_key = "noahmp_soil_stas_parameters";
        } else {
            amrex::Abort("Unrecognized soil dataset type, expected 'stas_ruc' or 'stas'");
        }

        if (veg_dataset == "usgs") {
            veg_param_key = "noahmp_usgs_parameters";
        } else if (veg_dataset == "modis") {
            veg_param_key = "noahmp_modis_parameters";
        } else {
            amrex::Abort("Unrecognized vegetation dataset type, expected 'usgs' or 'modis'");
        }

        // Copy soil and vegetation datasets to GPU
        auto &h_soil_params = full_params.at(soil_param_key);
        auto &h_veg_params = full_params.at(veg_param_key);

        for (int i = 0; i < h_soil_params.size(); i++) {
            const int varsize = h_soil_params[i].second.size();
            if (varsize > 1) {
                //amrex::Print() << " -- copying soil param '" << h_soil_params[i].first << "' to GPU " << std::endl;
                amrex::Gpu::PinnedVector<amrex::Real> *d_var = new amrex::Gpu::PinnedVector<amrex::Real>(varsize);
                d_soil_params.insert(std::make_pair(h_soil_params[i].first, d_var));
                Gpu::copyAsync(Gpu::hostToDevice, h_soil_params[i].second.data(), h_soil_params[i].second.data()+varsize, d_var->data());
            }
        }

        for (int i = 0; i < h_veg_params.size(); i++) {
            const int varsize = h_veg_params[i].second.size();
            if (varsize > 1) {
                //amrex::Print() << " -- copying veg param '" << h_veg_params[i].first << "' to GPU " << std::endl;
                amrex::Gpu::PinnedVector<amrex::Real> *d_var = new amrex::Gpu::PinnedVector<amrex::Real>(varsize);
                d_veg_params.insert(std::make_pair(h_veg_params[i].first, d_var));
                Gpu::copyAsync(Gpu::hostToDevice, h_veg_params[i].second.data(), h_veg_params[i].second.data()+varsize, d_var->data());
            }
        }

        Gpu::streamSynchronize();

        // Read radiation parameters from NoahmpTable.TBL
        // Read noahmp_rad_parameters block
        if (full_params.count("noahmp_rad_parameters") > 0) {
            auto &h_rad_params = full_params.at("noahmp_rad_parameters");
            for (int i = 0; i < h_rad_params.size(); i++) {
                std::string pname = h_rad_params[i].first;
                const int varsize = h_rad_params[i].second.size();

                if (pname == "albsat_vis") {
                    albsat_vis.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_rad_params[i].second.data(),
                                   h_rad_params[i].second.data()+varsize, albsat_vis.data());
                } else if (pname == "albsat_nir") {
                    albsat_nir.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_rad_params[i].second.data(),
                                   h_rad_params[i].second.data()+varsize, albsat_nir.data());
                } else if (pname == "albdry_vis") {
                    albdry_vis.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_rad_params[i].second.data(),
                                   h_rad_params[i].second.data()+varsize, albdry_vis.data());
                } else if (pname == "albdry_nir") {
                    albdry_nir.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_rad_params[i].second.data(),
                                   h_rad_params[i].second.data()+varsize, albdry_nir.data());
                } else if (pname == "alblak") {
                    alblak_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_rad_params[i].second.data(),
                                   h_rad_params[i].second.data()+varsize, alblak_rad.data());
                } else if (pname == "omegas") {
                    omegas_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_rad_params[i].second.data(),
                                   h_rad_params[i].second.data()+varsize, omegas_rad.data());
                } else if (pname == "betads" && varsize == 1) {
                    betads_rad = h_rad_params[i].second[0];
                } else if (pname == "betais" && varsize == 1) {
                    betais_rad = h_rad_params[i].second[0];
                }
            }
        }

        // Read noahmp_global_parameters block for snow albedo parameters
        if (full_params.count("noahmp_global_parameters") > 0) {
            auto &h_global_params = full_params.at("noahmp_global_parameters");
            for (int i = 0; i < h_global_params.size(); i++) {
                std::string pname = h_global_params[i].first;
                const int varsize = h_global_params[i].second.size();

                if (varsize == 1) {
                    if (pname == "tau0") tau0_rad = h_global_params[i].second[0];
                    else if (pname == "grain_growth") grain_growth_rad = h_global_params[i].second[0];
                    else if (pname == "extra_growth") extra_growth_rad = h_global_params[i].second[0];
                    else if (pname == "dirt_soot") dirt_soot_rad = h_global_params[i].second[0];
                    else if (pname == "bats_cosz") bats_cosz_rad = h_global_params[i].second[0];
                    else if (pname == "bats_vis_new") bats_vis_new_rad = h_global_params[i].second[0];
                    else if (pname == "bats_nir_new") bats_nir_new_rad = h_global_params[i].second[0];
                    else if (pname == "bats_vis_age") bats_vis_age_rad = h_global_params[i].second[0];
                    else if (pname == "bats_nir_age") bats_nir_age_rad = h_global_params[i].second[0];
                    else if (pname == "bats_vis_dir") bats_vis_dir_rad = h_global_params[i].second[0];
                    else if (pname == "bats_nir_dir") bats_nir_dir_rad = h_global_params[i].second[0];
                    else if (pname == "swemx") swemx_rad = h_global_params[i].second[0];
                    else if (pname == "snow_emis") snow_emis_rad = h_global_params[i].second[0];
                    else if (pname == "rsurf_exp") rsurf_exp = h_global_params[i].second[0];
                } else if (varsize == 2) {
                    if (pname == "eg") {
                        eg_soil_rad = h_global_params[i].second[0];  // IST=1 (soil)
                        eg_lake_rad = h_global_params[i].second[1];  // IST=2 (lake)
                    }
                }
            }
        }

        // Read vegetation parameters for radiation from modis/usgs parameters
        if (full_params.count(veg_param_key) > 0) {
            auto &h_veg_rad_params = full_params.at(veg_param_key);
            for (int i = 0; i < h_veg_rad_params.size(); i++) {
                std::string pname = h_veg_rad_params[i].first;
                const int varsize = h_veg_rad_params[i].second.size();

                if (pname == "rhol_vis") {
                    rhol_vis_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, rhol_vis_rad.data());
                } else if (pname == "rhol_nir") {
                    rhol_nir_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, rhol_nir_rad.data());
                } else if (pname == "rhos_vis") {
                    rhos_vis_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, rhos_vis_rad.data());
                } else if (pname == "rhos_nir") {
                    rhos_nir_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, rhos_nir_rad.data());
                } else if (pname == "taul_vis") {
                    taul_vis_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, taul_vis_rad.data());
                } else if (pname == "taul_nir") {
                    taul_nir_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, taul_nir_rad.data());
                } else if (pname == "taus_vis") {
                    taus_vis_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, taus_vis_rad.data());
                } else if (pname == "taus_nir") {
                    taus_nir_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, taus_nir_rad.data());
                // NOTE: Skip "xl" - using existing Khai_L from init_landtype()
                } else if (pname == "rc") {
                    rc_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, rc_rad.data());
                // NOTE: Skip "hvt" - using existing ztop from init_landtype()
                } else if (pname == "hvb") {
                    hvb_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, hvb_rad.data());
                } else if (pname == "den") {
                    den_rad.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, den_rad.data());
                } else if (pname == "cwpvt") {
                    cwpvt.resize(varsize);
                    Gpu::copyAsync(Gpu::hostToDevice, h_veg_rad_params[i].second.data(),
                                   h_veg_rad_params[i].second.data()+varsize, cwpvt.data());
                }
            }
        }

        Gpu::streamSynchronize();

        amrex::Print() << " SLM: NOAHMP radiation parameters loaded from " << parameter_file << std::endl;
    }

    if (interpolate_lai) {
        pp.gettable("lai", lai_table);
        pp.gettable("sai", sai_table);

        // validate LAI and SAI tables
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lai_table.size() == 12, "Invalid LAI table size, expected values for all 12 months");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sai_table.size() == 12, "Invalid SAI table size, expected values for all 12 months");
        num_landtypes = lai_table[0].size(); // use first entry as size, check all others against
        for (int i = 0; i < lai_table.size(); i++) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lai_table[i].size() == num_landtypes, "Invalid LAI table - inconsistent number of landtype entries between months");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sai_table[i].size() == num_landtypes, "Invalid SAI table - inconsistent number of landtype entries between months");
        }

        d_lai_curr.resize(num_landtypes);
        d_lai_next.resize(num_landtypes);
        d_sai_curr.resize(num_landtypes);
        d_sai_next.resize(num_landtypes);
    }

    if (use_wrf_lai) {
        pp.gettable("vegparam", param_table);
        amrex::Print() << " param table = " << std::endl;

        int nparam = param_table[0].size();
        for (int t = 0; t < param_table.size(); t++) {
            amrex::Print() << "     LANDTYPE " << t << ": ";
            for (int i = 0; i < param_table[t].size(); i++) {
                amrex::Print() << param_table[t][i] << " ";
            }
            amrex::Print() << std::endl;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(param_table[t].size() == nparam, "Invalid param table, inconsistent number of parameters for landtype");
        }

        // Copy parameter table to GPU
        d_param_table = TableData<Real, 2>({0, 0}, {static_cast<int>(param_table.size()), nparam});
        amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
        TableData<Real, 2> h_param_table({0, 0}, {static_cast<int>(param_table.size()), nparam}, amrex::The_Pinned_Arena());
        auto const &h_tab = h_param_table.table();
#else
        auto const &h_tab = d_param_table.table();
#endif
        for (int t = 0; t < param_table.size(); t++) {
            for (int i = 0; i < param_table[t].size(); i++) {
                h_tab(t, i) = param_table[t][i];
            }
        }
#ifdef AMREX_USE_GPU
        d_param_table.copy(h_param_table);
        Gpu::streamSynchronize();
#endif
    }

    pp.query("rad_input_file", rad_input_file);
    if (rad_input_file != "") {
#ifndef ERF_USE_NETCDF
        amrex::Abort("ERF needs to be compiled with NetCDF to use the SLM rad_input_file option!");
#endif
#ifdef ERF_USE_NETCDF
        ncutils::NCFile rad_forcing_ncf = ncutils::NCFile::open(rad_input_file, NC_NOWRITE);

        ncutils::NCVar rad_time = rad_forcing_ncf.var("time");
        num_rad_times = rad_time.shape()[0];
        rad_times.resize(num_rad_times);
        rad_time.get(rad_times.dataPtr(), {0}, {static_cast<unsigned long>(num_rad_times)});

        amrex::Real start_rad_time = rad_times[0];
        for (int i = 0; i < num_rad_times; i++)
        {
            rad_times[i] = (rad_times[i] - start_rad_time) * 86400.0; // shift relative to first time
        }

        ncutils::NCDim rad_x = rad_forcing_ncf.dim("x");
        ncutils::NCDim rad_y = rad_forcing_ncf.dim("y");

        amrex::Print() << " Read radiation forcing file '" << rad_input_file << "'" << std::endl;
        amrex::Print() << "   rad forcing file has " << num_rad_times << " time values, x = " << rad_x.len() << " y = " << rad_y.len() << std::endl;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_lsm_geom.Domain().length(0) == rad_x.len() && m_lsm_geom.Domain().length(1) == rad_y.len(), "Radiation input file must have same X and Y dimensions!");

        BoxList bl_rad = ba_lsm_2d.boxList();
        for (auto& b : bl_rad) {
            b.setSmall(2, 0);
            b.setBig(2, num_rad_times - 1);
        }
        BoxArray ba_rad = BoxArray(std::move(bl_rad));
        IntVect ng_2d(0, 0, 0);

        // rad_input has 6 components: SWVIS, SWNIR, SWVISD, SWNIRD, COSZRS, LWDS
        rad_input_data.define(ba_rad, landmask.distributionMap, 6, ng_2d);
        rad_input_data.setVal(0.0);

        amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
            Arena_Used = amrex::The_Pinned_Arena();
#endif

        for(int comp = 0; comp < rad_input_data.nComp(); ++comp) {

            ncutils::NCVar rad_var = rad_forcing_ncf.var(rad_names[comp]);
            amrex::Print() << " Reading radiation input var " << rad_var.name() << std::endl;

            for ( MFIter mfi(rad_input_data); mfi.isValid(); ++mfi) {
                const auto& box = mfi.fabbox();

                FArrayBox tmp(box, 1, Arena_Used);
                auto *dataPtr = tmp.dataPtr();
                AMREX_ALWAYS_ASSERT(dataPtr != nullptr);

                std::vector<size_t> starts = {0,
                                              static_cast<unsigned long>(box.loVect()[1]),
                                              static_cast<unsigned long>(box.loVect()[0])};
                std::vector<size_t> counts = {static_cast<unsigned long>(num_rad_times),
                                              static_cast<unsigned long>(box.length()[1]),
                                              static_cast<unsigned long>(box.length()[0])};

                amrex::Print() << "  -- reading {" << counts[0] << "," << counts[1] << "," << counts[2] << "} at index {" << starts[0] << "," << starts[1] << "," << starts[2] << "}" << std::endl;
                rad_var.get(dataPtr, starts, counts);

                auto mf_arr = rad_input_data.array(mfi, comp);
                auto tmp_arr = tmp.array();
                ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    mf_arr(i, j, k) = tmp_arr(i, j, k);
                });
            }
        }

        /*
         for ( MFIter mfi(rad_input_data); mfi.isValid(); ++mfi)
         {
            const auto& box = mfi.fabbox();
            auto rad_arr = rad_input_data.const_array(mfi);
            ParallelFor(box, rad_names.size(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                amrex::Print() << "   i = " << i << " j = " << j << " k = " << k << " n = " << n << " RADVAR = " << rad_arr(i, j, k, n) << std::endl;
            });
         }
         */
#endif
    }

    if (set_from_file)
    {
        pp.query("SLM_num_ref_inputs", num_ref_inputs);
        pp.query("SLM_ref_sounding_file", ref_sounding_file);
        pp.query("SLM_ref_flux_file", ref_flux_file);
        pp.query("SLM_ref_sst_file", ref_sst_file);

        pp.query("start_time", start_time);
        pp.query("time_unit", time_unit);

        // Reads the SLM input flux file assuming the following fields:
        //  t[day], swdn[W/m2], lwdn[W/m2], swup[W/m2], lwup[W/m2]
        fluxes = SLM::read_cols(ref_flux_file, 1);

        // Reads the SLM input SST file assuming the following fields:
        // t[day], sst[K], precip[mm/s]
        sst = SLM::read_cols(ref_sst_file, 1);

        // Reads the SLM input sounding file assuming the following fields:
        //  t[day], pres0[mb], tabs[K], q[g/kg], uvel[m/s], vvel[m/s]
        sounding = SLM::read_cols(ref_sounding_file, 1);

        // setup initial time value and input file indices
        if (start_time != -1.0) {
            // if the user provided a start time, find which row index it occurs at in the file
            auto it = std::find(sst[0].begin(), sst[0].end(), start_time);
            if (it != sst[0].end())
            {
                time_index = it - sst[0].begin();
                time = sst[0][time_index];

            } else {
                amrex::Error("SLM: invalid start_time value: could not find time value in input file!");
            }
        } else {
            // if the user did not provide a starting time, then default to the first time in the file
            time = sst[0][0];
            time_index = 0;
        }

        start_time_index = time_index;

        lsm_fab_vars[LsmVar_SLM::precipref]->setVal(sst[2][time_index]);

        lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->setVal(fluxes[1][time_index]);
        lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->setVal(0.0);
        lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->setVal(0.0);
        lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->setVal(0.0);

        lsm_fab_vars[LsmVar_SLM::lwref]->setVal(fluxes[2][time_index]);
        lsm_fab_vars[LsmVar_SLM::coszrsxy]->setVal(1.0);

        const amrex::Real pres = sounding[1][time_index] * 100.0;
        const amrex::Real qv = sounding[3][time_index] / 1000.0;
        const amrex::Real theta = getThgivenTandP(sounding[2][time_index], pres, R_d / Cp_d);
        lsm_fab_vars[LsmVar_SLM::tref]->setVal(sounding[2][time_index]);
        lsm_fab_vars[LsmVar_SLM::qref]->setVal(qv);
        lsm_fab_vars[LsmVar_SLM::pref]->setVal(pres / 100.0);
        //lsm_fab_vars[LsmVar_SLM::dref]->setVal(getRhogivenThetaPress(theta, pres, R_d / Cp_d, qv));
        lsm_fab_vars[LsmVar_SLM::dref]->setVal(pres / (rair * sounding[2][time_index]));

        lsm_fab_vars[LsmVar_SLM::uref]->setVal(sounding[4][time_index]);
        lsm_fab_vars[LsmVar_SLM::vref]->setVal(sounding[5][time_index]);

        sstxy.setVal(sst[1][time_index]);

    }
    shf_soil.setVal(0.0);
    lhf_soil.setVal(0.0);

    r_a.setVal(1.0e9);
    r_b.setVal(1.0e4);
    r_c.setVal(1.0e4);
    r_d.setVal(1.0e4);

    mw.setVal(0.0);
    mws.setVal(0.0);

    mw_inc.setVal(0.0);
    evapo_dry.setVal(0.0);
}

void SLM::time_interp_from_ref()
{
    // performs time interpolation of SLM inputs from reference data instead of
    // using ERF inputs
    if (set_from_file)
    {
        // Increment time used when reading input from testing files.
        //  Note: file times are in days, so we convert our dt to days
        time += (m_dt * (time_unit / 86400.0));

        amrex::Real t0 = sst[0][time_index];
        amrex::Real t1 = sst[0][time_index+1];
        while (time >= t1)
        {
            int prev_index = time_index;
            // shift time index to next window
            time_index = std::min(time_index + 1, int(sst[0].size() - start_time_index - 2));
            t0 = sst[0][time_index];
            t1 = sst[0][time_index+1];
            if (prev_index == time_index) {
                break;
            }
        }

        // interpolate values from input files based on the current time
        amrex::Real precip_interp = linear_interp(t0, t1, time, sst[2][time_index], sst[2][time_index + 1]);
        amrex::Real sw_interp = linear_interp(t0, t1, time, fluxes[1][time_index], fluxes[1][time_index + 1]);
        amrex::Real lw_interp = linear_interp(t0, t1, time, fluxes[2][time_index], fluxes[2][time_index + 1]);
        amrex::Real sst_interp = linear_interp(t0, t1, time, sst[1][time_index], sst[1][time_index + 1]);

        amrex::Real p_interp = linear_interp(t0, t1, time, sounding[1][time_index], sounding[1][time_index + 1]);
        amrex::Real t_interp = linear_interp(t0, t1, time, sounding[2][time_index], sounding[2][time_index + 1]);
        amrex::Real q_interp = linear_interp(t0, t1, time, sounding[3][time_index], sounding[3][time_index + 1]);
        amrex::Real u_interp = linear_interp(t0, t1, time, sounding[4][time_index], sounding[4][time_index + 1]);
        amrex::Real v_interp = linear_interp(t0, t1, time, sounding[5][time_index], sounding[5][time_index + 1]);
        lsm_fab_vars[LsmVar_SLM::precipref]->setVal(precip_interp);

        // total SW split into:
        //   Diffuse = ~30% SW
        //   Direct  = ~70% SW
        // Visible and NIR = 50/50%
        lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->setVal(0.5*(sw_interp*0.7));
        lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->setVal(0.5*(sw_interp*0.7));
        lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->setVal(0.5*(sw_interp*0.3));
        lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->setVal(0.5*(sw_interp*0.3));

        lsm_fab_vars[LsmVar_SLM::lwref]->setVal(lw_interp);
        lsm_fab_vars[LsmVar_SLM::coszrsxy]->setVal(1.0);

        const amrex::Real pres = p_interp * 100.0;
        const amrex::Real qv = q_interp / 1000.0;
        const amrex::Real theta = getThgivenTandP(t_interp, pres, R_d / Cp_d);
        lsm_fab_vars[LsmVar_SLM::tref]->setVal(t_interp);
        lsm_fab_vars[LsmVar_SLM::qref]->setVal(qv);
        lsm_fab_vars[LsmVar_SLM::pref]->setVal(pres / 100.0);
        //lsm_fab_vars[LsmVar_SLM::dref]->setVal(getRhogivenThetaPress(theta, pres, R_d / Cp_d, qv));

        lsm_fab_vars[LsmVar_SLM::dref]->setVal(pres / (rair * t_interp));
        lsm_fab_vars[LsmVar_SLM::uref]->setVal(u_interp);
        lsm_fab_vars[LsmVar_SLM::vref]->setVal(v_interp);

        sstxy.setVal(sst_interp);
    }
}

/**
 * Initialize SLM
 */
void SLM::slm_init()
{
    // validate the landtype flag
    const int landtype_min = landtype.min(0);
    const int landtype_max = landtype.max(0);
    if (landtype_min < 0 || landtype_max > 16)
    {
        amrex::Abort("SLM: landtype values are outside of valid 0-16 range! min = " +
                     std::to_string(landtype_min) + ", max = " + std::to_string(landtype_max));
    }

    // validate LAI - check that they are within [0,10]
    const double lai_min = LAI.min(0);
    const double lai_max = LAI.max(0);
    if (lai_min < 0.0 || lai_max > 10.0)
    {
      amrex::Abort("SLM: LAI values are outside of valid [0-10] range! min = " +
                   std::to_string(lai_min) + ", max = " + std::to_string(lai_max));
    }

    // Sets model properties based on the landtype
    init_landtype();

    // TODO: VEGYES and set minimum LAI
    vege_YES.setVal(1.0);

    for ( MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();

        auto vege_YES_arr = vege_YES.array(mfi);
        auto vegetype_arr = vegetype.const_array(mfi);
        auto nroot_arr = nroot.array(mfi);
        auto landtype_arr = landtype.const_array(mfi);

        auto IR_emis_vege_arr = IR_emis_vege.array(mfi);
        auto IR_emis_soil_arr = IR_emis_soil.array(mfi);
        auto IR_emis_grnd_arr = IR_emis_grnd.array(mfi);
        auto phi_1_arr = phi_1.array(mfi);
        auto phi_2_arr = phi_2.array(mfi);
        auto precip_extinc_arr = precip_extinc.array(mfi);
        auto mw_mx_arr = mw_mx.array(mfi);
        auto LAI_arr = LAI.array(mfi);
        auto BAI_arr = BAI.const_array(mfi);
        auto IMPERV_arr = IMPERV.const_array(mfi);
        auto ztop_arr = ztop.const_array(mfi);
        auto Khai_L_arr = Khai_L.array(mfi);
        auto landmask_arr = landmask.const_array(mfi);

        auto veg_frac_arr = lsm_fab_vars[LsmVar_SLM::veg_frac]->array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) == 1) {
                if (vegetype_arr(i, j, 0) == 0) {
                    vege_YES_arr(i, j, 0) = 0.0;
                    veg_frac_arr(i, j, 0) = 0.0;
                    nroot_arr(i, j, 0) = 0;
                } else {
                    // set minimum LAI for vegetated land
                    LAI_arr(i, j, 0) = std::max(LAI_arr(i, j, 0), 0.001);
                }

                phi_1_arr(i, j, 0) = 0.5 - 0.633 * Khai_L_arr(i, j, 0) - 0.33 * (std::pow(Khai_L_arr(i, j, 0), 2));
                phi_2_arr(i, j, 0) = 0.877 * (1.0 - 2.0 * phi_1_arr(i, j, 0));
                IR_emis_vege_arr(i, j, 0) = 0.97 * (1.0 - std::exp(-1.0 * (phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0)) * LAI_arr(i, j, 0)));
                precip_extinc_arr(i, j, 0) = phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0);
                //mw_mx_arr(i, j, 0) = 0.1 * LAI_arr(i, j, 0);
                //add the basal area to mw_mx. circumference of trunks per unit area = sqrt(4*PI*basal area), which then is converted from sq.feet/sq acre to m2/m2 by diving by 43560
                mw_mx_arr(i, j, 0) = 0.1 * LAI_arr(i, j, 0) + ztop_arr(i, j, 0) * std::pow(4.0 * PI * BAI_arr(i, j, 0) / 43560., 0.5);

                if(landtype_arr(i, j, 0) == 13)
                {
                    IR_emis_grnd_arr(i, j, 0) = IR_emis_urban * IMPERV_arr(i, j, 0) + IR_emis_soil_arr(i, j, 0) * (1.-IMPERV_arr(i, j, 0)); // urban IR emissivity
                }
                else
                {
                    IR_emis_grnd_arr(i, j, 0) = IR_emis_soil_arr(i, j, 0);
                }
            } else {
                vege_YES_arr(i, j, 0) = 0.0;
                veg_frac_arr(i, j, 0) = 0.0;
            }
        });
    }

    // Initialize soil parameters
    init_soil_tw();

    init_soil_vars();

    // Calculate fraction of root in each soil layer
    vege_root_init();

    // Update any parameters using the values from a parameter file
    init_from_params();
}

/**
 * Helper function to set properties according to the IGBP class
 */
void SLM::init_landtype()
{
    const Real d_z0_soil = z0_soil;

    for ( MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();

        auto landtype_arr = landtype.const_array(mfi);
        auto landmask_arr = landmask.array(mfi);
        auto LAI_arr = LAI.const_array(mfi);

        auto albedovis_v_arr = albedovis_v.array(mfi);
        auto albedonir_v_arr = albedonir_v.array(mfi);
        auto albedovis_s_arr = albedovis_s.array(mfi);
        auto albedonir_s_arr = albedonir_s.array(mfi);
        auto ztop_arr = ztop.array(mfi);
        auto disp_hgt_arr = disp_hgt.array(mfi);
        auto z0_sfc_arr = z0_sfc.array(mfi);
        auto Khai_L_arr = Khai_L.array(mfi);
        auto rootL_arr = rootL.array(mfi);
        auto root_a_arr = root_a.array(mfi);
        auto root_b_arr = root_b.array(mfi);
        auto Rc_min_arr = Rc_min.array(mfi);
        auto Rgl_arr = Rgl.array(mfi);
        auto hs_rc_arr = hs_rc.array(mfi);
        auto BAI_arr = BAI.array(mfi);
        auto IMPERV_arr = IMPERV.array(mfi);
        auto vegetype_arr = vegetype.array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            switch (landtype_arr(i, j, 0))
            {
                case 0: // water
                    // TODO: should this landmask be ERF's lmask_lev variable?
                    landmask_arr(i, j, 0) = 0;
                    vegetype_arr(i, j, 0) = 0;
                    // set albedo for water so that it is populated for RRTMGP
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.06;
                    albedovis_s_arr(i, j, 0) = 0.06;
                    albedonir_s_arr(i, j, 0) = 0.06;
                    break;
                case 1: // evergreen needleleaf forest
                    albedovis_v_arr(i, j, 0) = 0.094;
                    albedonir_v_arr(i, j, 0) = 0.161;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 1.09;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.8;
                    root_a_arr(i, j, 0) = 6.706;
                    root_b_arr(i, j, 0) = 2.175;
                    Rc_min_arr(i, j, 0) = 250.;
                    Rgl_arr(i, j, 0) = 120.;
                    //Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 0.03;
                    BAI_arr(i, j, 0) = 200.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 2: // evergreen broadleaf forest
                    albedovis_v_arr(i, j, 0) = 0.086;
                    albedonir_v_arr(i, j, 0) = 0.146;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 1.00;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 3.0;
                    root_a_arr(i, j, 0) = 7.344;
                    root_b_arr(i, j, 0) = 1.303;
                    Rc_min_arr(i, j, 0) = 250.;
                    Rgl_arr(i, j, 0) = 120.;
                    hs_rc_arr(i, j, 0) = 0.03;
                    BAI_arr(i, j, 0) = 200.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 3: // deciduous needleaf forest
                    albedovis_v_arr(i, j, 0) = 0.102;
                    albedonir_v_arr(i, j, 0) = 0.198;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 1.10;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 2.0;
                    root_a_arr(i, j, 0) = 7.066;
                    root_b_arr(i, j, 0) = 1.953;
                    Rc_min_arr(i, j, 0) = 250.;
                    Rgl_arr(i, j, 0) = 120.;
                    hs_rc_arr(i, j, 0) = 0.03;
                    BAI_arr(i, j, 0) = 200.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 4: // deciduous broadleaf forest
                    albedovis_v_arr(i, j, 0) = 0.056;
                    albedonir_v_arr(i, j, 0) = 0.151;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.8;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 2.0;
                    root_a_arr(i, j, 0) = 5.990;
                    root_b_arr(i, j, 0) = 1.955;
                    Rc_min_arr(i, j, 0) = 250.;
                    Rgl_arr(i, j, 0) = 120.;
                    hs_rc_arr(i, j, 0) = 0.03;
                    BAI_arr(i, j, 0) = 200.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 5: // mixed forest
                    albedovis_v_arr(i, j, 0) = 0.093;
                    albedonir_v_arr(i, j, 0) = 0.179;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.8;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 2.4;
                    root_a_arr(i, j, 0) = 4.453;
                    root_b_arr(i, j, 0) = 1.631;
                    Rc_min_arr(i, j, 0) = 250.;
                    Rgl_arr(i, j, 0) = 120.;
                    hs_rc_arr(i, j, 0) = 0.03;
                    BAI_arr(i, j, 0) = 200.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 6: // closed shrublands
                    albedovis_v_arr(i, j, 0) = 0.069;
                    albedonir_v_arr(i, j, 0) = 0.121;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 1.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.1;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 2.5;
                    root_a_arr(i, j, 0) = 6.326;
                    root_b_arr(i, j, 0) = 1.567;
                    Rc_min_arr(i, j, 0) = 220.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 60.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 7: // open shrublands
                    albedovis_v_arr(i, j, 0) = 0.079;
                    albedonir_v_arr(i, j, 0) = 0.226;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 1.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.1;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 3.1;
                    root_a_arr(i, j, 0) = 7.718;
                    root_b_arr(i, j, 0) = 1.262;
                    Rc_min_arr(i, j, 0) = 220.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 60.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 8: // woody savannas
                    albedovis_v_arr(i, j, 0) = 0.071;
                    albedonir_v_arr(i, j, 0) = 0.085;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 5.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.3;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.7;
                    root_a_arr(i, j, 0) = 7.604;
                    root_b_arr(i, j, 0) = 2.300;
                    Rc_min_arr(i, j, 0) = 180.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.02;
                    BAI_arr(i, j, 0) = 100.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 9: // savannas
                    albedovis_v_arr(i, j, 0) = 0.061;
                    albedonir_v_arr(i, j, 0) = 0.227;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 5.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.3;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 2.4;
                    root_a_arr(i, j, 0) = 8.235;
                    root_b_arr(i, j, 0) = 1.627;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 100.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 10: // grasslands
                    albedovis_v_arr(i, j, 0) = 0.09;
                    albedonir_v_arr(i, j, 0) = 0.269;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.04;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 10.74;
                    root_b_arr(i, j, 0) = 2.608;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 20.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 11: // permanent wetlands
                    albedovis_v_arr(i, j, 0) = 0.081;
                    albedonir_v_arr(i, j, 0) = 0.18;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.2;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.0;
                    Rgl_arr(i, j, 0) = 100.0;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 20.0;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) =1;
                    break;
                case 12: // croplands
                    albedovis_v_arr(i, j, 0) = 0.084;
                    albedonir_v_arr(i, j, 0) = 0.193;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.03;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 60.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 13: // urban
                    albedovis_v_arr(i, j, 0) = 0.102;
                    albedonir_v_arr(i, j, 0) = 0.164;
                    albedovis_s_arr(i, j, 0) = 0.15;
                    albedonir_s_arr(i, j, 0) = 0.25;
                    ztop_arr(i, j, 0) = 10.;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 1.0;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 0.0;
                    Rgl_arr(i, j, 0) = 0.0;
                    hs_rc_arr(i, j, 0) = 0.0;
                    BAI_arr(i, j, 0) = 0.0;
                    IMPERV_arr(i, j, 0) = 0.75;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                case 14: // croplands/natural mozaics
                    albedovis_v_arr(i, j, 0) = 0.071;
                    albedonir_v_arr(i, j, 0) = 0.169;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.04;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 0.01;
                    BAI_arr(i, j, 0) = 20.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 15: // snow/ice
                    albedovis_v_arr(i, j, 0) = 0.0; // actually depends on moisture content
                    albedonir_v_arr(i, j, 0) = 0.0;
                    albedovis_s_arr(i, j, 0) = 0.91;
                    albedonir_s_arr(i, j, 0) = 0.65;
                    ztop_arr(i, j, 0) = 0.0;
                    disp_hgt_arr(i, j, 0) = 0.0;
                    z0_sfc_arr(i, j, 0) = 0.001;
                    Khai_L_arr(i, j, 0) = 0.0;
                    rootL_arr(i, j, 0) = 0.0;
                    root_a_arr(i, j, 0) = 0.0;
                    root_b_arr(i, j, 0) = 0.0;
                    Rc_min_arr(i, j, 0) = 0.0;
                    Rgl_arr(i, j, 0) = 0.0;
                    hs_rc_arr(i, j, 0) = 0.0;
                    BAI_arr(i, j, 0) = 0.0;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                case 16: // baresoil
                    albedovis_v_arr(i, j, 0) = 0.0; // actually depends on moisture content
                    albedonir_v_arr(i, j, 0) = 0.0;
                    albedovis_s_arr(i, j, 0) = 0.208;
                    albedonir_s_arr(i, j, 0) = 0.344;
                    ztop_arr(i, j, 0) = 0.;
                    disp_hgt_arr(i, j, 0) = 0.0;
                    z0_sfc_arr(i, j, 0) = d_z0_soil;
                    Khai_L_arr(i, j, 0) = 0.;
                    rootL_arr(i, j, 0) = 0.;
                    root_a_arr(i, j, 0) = 0.;
                    root_b_arr(i, j, 0) = 0.;
                    Rc_min_arr(i, j, 0) = 0.;
                    Rgl_arr(i, j, 0) = 0.;
                    hs_rc_arr(i, j, 0) = 0.;
                    BAI_arr(i, j, 0) = 0.;
                    IMPERV_arr(i, j, 0) = 0.;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                default:
// TODO: check if AMReX::Abort can be called inside ParFor?
#ifndef AMREX_USE_GPU
                  amrex::Abort(
                    "landtype invalid for i = " + std::to_string(i) +
                    " j = " + std::to_string(j));
#endif
                  break;
            }
#ifndef AMREX_USE_GPU
            if (vegetype_arr(i, j, 0) == 1 && LAI_arr(i, j, 0) == 0.0)
            {
                amrex::Abort("LAI is not set for vegetated land point at i = " + std::to_string(i) + " j = " + std::to_string(j));
            }
#endif
        });
  }
}

/**
 * Reads soil input file, initializes soil wetness and temperature
 */
void SLM::init_soil_tw()
{
    // TODO - read from input file
    const Real d_tabs_s = tabs_s;
    const Real d_t00 = t00;

    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    const bool wrfinput = use_wrfinput;

    amrex::Gpu::DeviceVector<Real> d_dz_lsm_vec(m_dz_lsm.size());
    amrex::Gpu::DeviceVector<Real> d_clay_vec(m_dz_lsm.size());
    amrex::Gpu::DeviceVector<Real> d_sand_vec(m_dz_lsm.size());
    amrex::Gpu::DeviceVector<Real> d_st0_vec(m_dz_lsm.size());
    amrex::Gpu::DeviceVector<Real> d_sw0_vec(m_dz_lsm.size());
    amrex::Gpu::DeviceVector<Real> d_relax_vec(m_dz_lsm.size());

    amrex::Gpu::copy(Gpu::hostToDevice, m_dz_lsm.begin(), m_dz_lsm.end(), d_dz_lsm_vec.begin());
    amrex::Gpu::copy(Gpu::hostToDevice, clay0.begin(), clay0.end(), d_clay_vec.begin());
    amrex::Gpu::copy(Gpu::hostToDevice, sand0.begin(), sand0.end(), d_sand_vec.begin());
    amrex::Gpu::copy(Gpu::hostToDevice, st0.begin(), st0.end(), d_st0_vec.begin());
    amrex::Gpu::copy(Gpu::hostToDevice, sw0.begin(), sw0.end(), d_sw0_vec.begin());
    amrex::Gpu::copy(Gpu::hostToDevice, relax_hgt.begin(), relax_hgt.end(), d_relax_vec.begin());
    Real *d_dz_lsm = d_dz_lsm_vec.data();
    Real *d_clay0 = d_clay_vec.data();
    Real *d_sand0 = d_sand_vec.data();
    Real *d_st0 = d_st0_vec.data();
    Real *d_sw0 = d_sw0_vec.data();
    Real *d_relax = d_relax_vec.data();

    auto tsurf = lsm_fab_vars[LsmVar_SLM::tsurf];
    for ( MFIter mfi(*tsurf, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto landmask_arr = landmask.const_array(mfi);

        auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->array(mfi);
        auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->array(mfi);

        auto tsurf_arr   = lsm_fab_vars[LsmVar_SLM::tsurf]->array(mfi);
        auto lai_arr     = lsm_fab_vars[LsmVar_SLM::lai]->array(mfi);
        auto vegtype_arr = lsm_fab_vars[LsmVar_SLM::vegtype]->array(mfi);

        auto sand_arr = lsm_fab_vars[LsmVar_SLM::sand]->array(mfi);
        auto clay_arr = lsm_fab_vars[LsmVar_SLM::clay]->array(mfi);

        auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->array(mfi);
        auto soil_relax_hgt_arr = lsm_fab_vars[LsmVar_SLM::soil_relax_hgt]->array(mfi);

        auto sstxy_arr = sstxy.array(mfi);

        auto soiltype_arr = lsm_fab_vars[LsmVar_SLM::soiltype]->array(mfi);
        auto landtype_arr = landtype.const_array(mfi);
        auto LAI_local_arr = LAI.const_array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (landmask_arr(i, j, 0) == 1)
            {
                if (!wrfinput) {
                    s_depth_arr(i, j, k) = d_dz_lsm[(k*-1)+d_khi_lsm];
                    clay_arr(i, j, k) = d_clay0[(k*-1)+d_khi_lsm];
                    sand_arr(i, j, k) = d_sand0[(k*-1)+d_khi_lsm];
                    soilt_arr(i, j, k) = d_st0[(k*-1)+d_khi_lsm];
                    soilw_arr(i, j, k) = d_sw0[(k*-1)+d_khi_lsm];
                } else {
                    // initialize sand and clay % based on soil type from WRFInput
                    switch (static_cast<int>(soiltype_arr(i, j, d_khi_lsm)))
                    {
                        case 1: // sand
                            sand_arr(i,j,k) = 0.92;
                            clay_arr(i,j,k) = 0.03;
                            break;
                        case 2: // loamy sand
                            sand_arr(i,j,k) = 0.82;
                            clay_arr(i,j,k) = 0.06;
                            break;
                        case 3: // sandy loam
                            sand_arr(i,j,k) = 0.65;
                            clay_arr(i,j,k) = 0.10;
                            break;
                        case 4: // silt loam
                            sand_arr(i,j,k) = 0.20;
                            clay_arr(i,j,k) = 0.15;
                            break;
                        case 5: // silt
                            sand_arr(i,j,k) = 0.08;
                            clay_arr(i,j,k) = 0.12;
                            break;
                        case 6: // loam
                            sand_arr(i,j,k) = 0.40;
                            clay_arr(i,j,k) = 0.20;
                            break;
                        case 7: // sandy clay loam
                            sand_arr(i,j,k) = 0.60;
                            clay_arr(i,j,k) = 0.30;
                            break;
                        case 8: // clay loam
                            sand_arr(i,j,k) = 0.32;
                            clay_arr(i,j,k) = 0.34;
                            break;
                        case 9: // silty clay loam
                            sand_arr(i,j,k) = 0.20;
                            clay_arr(i,j,k) = 0.40;
                            break;
                        case 10: // sandy clay
                            sand_arr(i,j,k) = 0.52;
                            clay_arr(i,j,k) = 0.42;
                            break;
                        case 11: // silty clay
                            sand_arr(i,j,k) = 0.06;
                            clay_arr(i,j,k) = 0.47;
                            break;
                        case 12: // clay
                            sand_arr(i,j,k) = 0.20;
                            clay_arr(i,j,k) = 0.60;
                            break;
                        case 13: // organic material
                            //sand_arr(i,j,k) = 0.0;
                            //clay_arr(i,j,k) = 0.0;
                            sand_arr(i,j,k) = 0.001;
                            clay_arr(i,j,k) = 0.001;
                            break;
                        case 14: // water
                            sand_arr(i,j,k) = 0.001;
                            clay_arr(i,j,k) = 0.001;
                            break;
                        case 15: // bedrock
                            sand_arr(i,j,k) = 0.001;
                            clay_arr(i,j,k) = 0.001;
                            break;
                        case 16: // other (urban/builtup, etc)
                            sand_arr(i,j,k) = 0.001;
                            clay_arr(i,j,k) = 0.001;
                            break;
                        default:
                            sand_arr(i,j,k) = -1.0;
                            clay_arr(i,j,k) = -1.0;
                            break;
                    }

                    // convert to percentage
                    sand_arr(i,j,k) *= 100.0;
                    clay_arr(i,j,k) *= 100.0;
                }

                soil_relax_hgt_arr(i, j, k) = d_relax[(k*-1)+d_khi_lsm];

                // keep the 3D wrfinput copies consistent with the local SLM values if not using wrfinput
                if (!wrfinput && k == d_khi_lsm) {
                    tsurf_arr(i, j, d_khi_lsm) = d_st0[0];
                    tsurf_arr(i, j, 0) = d_st0[0];
                    lai_arr(i, j, d_khi_lsm) = LAI_local_arr(i, j, 0);
                    lai_arr(i, j, 0) = LAI_local_arr(i, j, 0);
                    vegtype_arr(i, j, d_khi_lsm) = landtype_arr(i, j, 0);
                    vegtype_arr(i, j, 0) = landtype_arr(i, j, 0);
                }

                // TODO: nrestart conditional here
                // TODO: sstxy - should be ERF surface temp array?
                sstxy_arr(i, j, 0) = soilt_arr(i, j, d_khi_lsm) - d_t00;
            } else {
                sstxy_arr(i, j, 0) = d_tabs_s - d_t00;
            }
        });
    }

    // Calculate node_z (depth of the center of each soil layer)
    init_layer_depths();
}

/**
 * Calculates center and face depths for each soil layer
 */
void SLM::init_layer_depths()
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    const bool wrfinput = use_wrfinput;

    for ( MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();

        auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->const_array(mfi);

        auto node_z_arr = lsm_fab_vars[LsmVar_SLM::node_z]->array(mfi);
        auto interface_z_arr = lsm_fab_vars[LsmVar_SLM::interface_z]->array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // k = 0 = reference level
            // k = -1 = top of soil layer = khi_lsm
            // k = -m_nz_lsm = klo_lsm

            if (!wrfinput) {
                node_z_arr(i, j, d_khi_lsm) = 0.5 * s_depth_arr(i, j, d_khi_lsm);
                for (int k = d_khi_lsm - 1; k >= d_klo_lsm; k--)
                {
                    node_z_arr(i, j, k) = 0.5 * s_depth_arr(i, j, k);

                    for (int kk = d_khi_lsm; kk >= k + 1; kk--)
                    {
                        node_z_arr(i, j, k) = node_z_arr(i, j, k) + s_depth_arr(i, j, kk);
                    }
                }
            }

            for (int k = d_khi_lsm; k >= d_klo_lsm; k--)
            {
                interface_z_arr(i, j, k) = node_z_arr(i, j, k) + 0.5*s_depth_arr(i, j, k);
            }
        });
    }
}

/**
 * Initializes soil properties based on clay and sand percentages
 */
void SLM::init_soil_vars()
{
    const bool d_param_updated = params_updated;

    auto tsurf = lsm_fab_vars[LsmVar_SLM::tsurf];
    for ( MFIter mfi(*tsurf, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto landmask_arr = landmask.const_array(mfi);

        auto sand_arr = lsm_fab_vars[LsmVar_SLM::sand]->const_array(mfi);
        auto clay_arr = lsm_fab_vars[LsmVar_SLM::clay]->const_array(mfi);

        auto w_s_FC_arr = lsm_fab_vars[LsmVar_SLM::w_s_FC]->array(mfi);
        auto w_s_WP_arr = lsm_fab_vars[LsmVar_SLM::w_s_WP]->array(mfi);
        auto sst_capa_arr = lsm_fab_vars[LsmVar_SLM::sst_capa]->array(mfi);
        auto sst_cond_arr = lsm_fab_vars[LsmVar_SLM::sst_cond]->array(mfi);
        auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->array(mfi);
        auto theta_FC_arr = lsm_fab_vars[LsmVar_SLM::theta_FC]->array(mfi);
        auto theta_WP_arr = lsm_fab_vars[LsmVar_SLM::theta_WP]->array(mfi);
        auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->array(mfi);
        auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->array(mfi);
        auto ks_arr = lsm_fab_vars[LsmVar_SLM::ks]->array(mfi);

        ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (landmask_arr(i, j, 0) == 1)
            {
                if (!d_param_updated) {
                    // compute the initial values
                    const amrex::Real sand_per = sand_arr(i, j, k); // sand percentage

                    // soil solids thermal conductivity (Johansen 1975)
                    //  quartz (= SAND content) thermal conductivity = 7.7 W/mK
                    const amrex::Real mineral_tcond = sand_per > 20.0 ? 2.0 : 3.0; // thermal conductivity of other minerals [W/mK]
                    sst_cond_arr(i, j, k) = std::pow(7.7, (sand_per * 0.01)) * std::pow(mineral_tcond, (1.0 - (sand_per * 0.01)));

                    // Calculated following Cosby et al. 1984 ( hydraulic properties)
                    // hydraulic conductivity at satuation , mm/s
                    ks_arr(i, j, k) = std::pow(10.0, (0.0153*sand_per) - 0.884) * (25.4 / 3600.0); // [mm/s] from [inch/hr]

                    // constant B
                    Bconst_arr(i, j, k) = 0.159 * clay_arr(i, j, k) + 2.91;

                    // porosity (or saturation volumetric water content)
                    poro_soil_arr(i, j, k) = -0.00126 * sand_per + 0.489; // volume/volume

                    // moisture potential at saturation, [mm]
                    m_pot_sat_arr(i, j, k) = std::min(-150.0, -10.0*(std::pow(10.0, 1.88 - 0.0131*sand_per))); // [mm] from [cm]

                    // soil heat capacity, [J/m^3/K]
                    //  Following de Vries(1963) using SAND=34% CLAY=63%
                    sst_capa_arr(i, j, k) = 1.0e6 * (2.128*sand_per + 2.385*clay_arr(i, j, k)) / (sand_per + clay_arr(i, j, k));

                    // volumetric moisture content at field capacity
                    // field capacity is assumed to be the occasion when hydraulic conductivity is 0.1mm/d
                    theta_FC_arr(i, j, k) = poro_soil_arr(i, j, k) * std::pow((0.1 / 86400.0 / ks_arr(i, j, k)), 1.0 / (2.0 * Bconst_arr(i, j, k) + 3.0));

                    // volumetric moisture content at wilting point
                    theta_WP_arr(i, j, k) = poro_soil_arr(i, j, k) * std::pow((-150000.0 / m_pot_sat_arr(i, j, k)), (-1.0 / Bconst_arr(i, j, k)));
                }

                // soil wetness at field capacity
                w_s_FC_arr(i, j, k) = theta_FC_arr(i, j, k) / poro_soil_arr(i, j, k);

                // soil wetness at wilting point
                w_s_WP_arr(i, j, k) = theta_WP_arr(i, j, k) / poro_soil_arr(i, j, k);
            }
        });
    }
}


/**
 * Assigns root fraction in each soil layer
 * Fraction of total root in each soil layer is determined based on the soil
 * depth and vegetation root parameters
 */
void SLM::vege_root_init()
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    for ( MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box = mfi.tilebox();

        auto landmask_arr = landmask.const_array(mfi);

        auto interface_z_arr = lsm_fab_vars[LsmVar_SLM::interface_z]->const_array(mfi);
        auto rootF_arr = lsm_fab_vars[LsmVar_SLM::rootF]->array(mfi);
        auto rootL_arr = rootL.const_array(mfi);

        auto root_a_arr = root_a.const_array(mfi);
        auto root_b_arr = root_b.const_array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) == 1)
            {
                int nrootind = d_khi_lsm;
                for (int k = d_klo_lsm; k <= d_khi_lsm - 1; k++)
                {
                    if (interface_z_arr(i, j, k) >= rootL_arr(i, j, 0) && interface_z_arr(i, j, k + 1) <= rootL_arr(i, j, 0))
                    {
                        nrootind = k;
                    }
                }

                rootF_arr(i, j, nrootind) = 1.0 - 0.5 * (std::exp(-1.0 * root_a_arr(i, j, 0) * rootL_arr(i, j, 0)) +
                                                         std::exp(-1.0 * root_b_arr(i, j, 0) * rootL_arr(i, j, 0)));

                const amrex::Real tot_root_density = rootF_arr(i, j, nrootind);
                for (int k = d_khi_lsm; k >= nrootind + 1; k--)
                {
                    rootF_arr(i, j, k) = 1.0 - 0.5 * (std::exp(-1.0 * root_a_arr(i, j, 0) * interface_z_arr(i, j, k)) +
                                                      std::exp(-1.0 * root_b_arr(i, j, 0) * interface_z_arr(i, j, k)));
                }

                if (nrootind < d_khi_lsm)
                {
                    for (int k = nrootind; k <= d_khi_lsm - 1; k++)
                    {
                        rootF_arr(i, j, k) = rootF_arr(i, j, k) - rootF_arr(i, j, k + 1);
                    }

                    // to ensure total root density equals 1
                    amrex::Real new_root_density = 0.0;
                    for (int k = d_khi_lsm; k >= nrootind; k--)
                    {
                        rootF_arr(i, j, k) = rootF_arr(i, j, k) / tot_root_density;
                        new_root_density += rootF_arr(i, j, k);
                    }

                    AMREX_ASSERT_WITH_MESSAGE(std::abs(new_root_density - 1.0) < 1.0e-12, "total root density should equal 1.0");
                }
            }
        });
    }
}

void SLM::init_slm_vars()
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    // Initializes SLM quantities from reference values tr, ts, qr
    for ( MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        auto box = mfi.tilebox();

        auto landmask_arr = landmask.const_array(mfi);

        auto t_cas_arr = t_cas.array(mfi);
        auto q_cas_arr = q_cas.array(mfi);

        auto t_canop_arr = t_canop.array(mfi);
        auto t_ground_skin_arr = t_ground_skin.array(mfi);

        auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
        auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
        auto tref_arr  = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
        auto tsurf_arr  = lsm_fab_vars[LsmVar_SLM::tsurf]->const_array(mfi);
        auto qref_arr  = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
        auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);

        //auto landtype_arr = landtype.const_array(mfi);
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) == 1)
            {
                amrex::Real q_gr;

                // soil surface specific humidity
                if (soilt_arr(i, j, d_khi_lsm) > tfriz)
                {
                    erf_qsatw(tsurf_arr(i, j, 0), pref_arr(i, j, 0), q_gr);
                }
                else
                {
                    erf_qsati(tsurf_arr(i, j, 0), pref_arr(i, j, 0), q_gr);
                }
                q_gr *= soilw_arr(i, j, d_khi_lsm);

                // canopy temperature: initialize as (ref level temperature + soil surf temperature)/2
                //t_canop_arr(i, j, 0) = 0.5*(tref_arr(i, j, 0) + tsurf_arr(i, j, 0));
                //t_cas_arr(i, j, 0) = 0.5*(tref_arr(i, j, 0) + tsurf_arr(i, j, 0));

                t_canop_arr(i, j, 0) = tref_arr(i, j, 0);
                t_cas_arr(i, j, 0) = tref_arr(i, j, 0);
                t_ground_skin_arr(i, j, 0) = soilt_arr(i, j, d_khi_lsm);

                // specific humidity in canopy air space
                q_cas_arr(i, j, 0) = 0.5*(qref_arr(i, j, 0) + q_gr);

                //amrex::Print() << " i = " << i << " j = " << j << " k = 0 : t_canop = " << t_canop_arr(i, j, 0) << " t_cas = " << t_cas_arr(i, j, 0) << " tsurf = " << tsurf_arr(i, j, 0) << " qcas = " << q_cas_arr(i, j, 0) << " qref = " << qref_arr(i, j, 0) << std::endl;

            }
        });
    }

    // initialize nudging profiles for soil based on the initial soilt and soilw
    MultiFab::Copy(*lsm_fab_vars[LsmVar_SLM::soilt_obs], *lsm_fab_vars[LsmVar_SLM::soilt], 0, 0, 1, 0);
    MultiFab::Copy(*lsm_fab_vars[LsmVar_SLM::soilw_obs], *lsm_fab_vars[LsmVar_SLM::soilw], 0, 0, 1, 0);
}

/**
 * Loads and parses the given .TBL file (such as NoahMPTable.TBL)
 *
 * Returns an unorderd map containing the name of each parameter block, where
 * each block contains a list of variables, with each variable containing one or more values.
 *
 * Returns the vegetation and soil category names loaded from file as a vector of names.
 */
SLMParameterTable SLM::ReadParameterFile(const std::string &filename, amrex::Vector<std::string> &veg_categories, amrex::Vector<std::string> &soil_categories)
{
    amrex::Print() << " SLM: Reading parameter file '" << filename << "'..." << std::endl;

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(filename, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line;

    std::unordered_map<std::string,
                       amrex::Vector<std::pair<std::string, amrex::Vector<amrex::Real>>>> table;

    const auto trim = [](std::string value) {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) { return std::string{}; }
        const auto last = value.find_last_not_of(" \t\r\n");
        return value.substr(first, last - first + 1);
    };

    std::string block_name;
    while (std::getline(is, line)) {
        const auto comment_pos = line.find('!');
        if (comment_pos != std::string::npos) { line.erase(comment_pos); }
        line = trim(line);
        if (line.empty()) { continue; }

        if (line[0] == '&') {
            block_name = amrex::toLower(trim(line.substr(1)));
            table.try_emplace(block_name);
            continue;
        }
        if (line[0] == '/') {
            block_name.clear();
            continue;
        }
        if (block_name.empty()) { continue; }

        const auto equal_pos = line.find('=');
        if (equal_pos == std::string::npos) { continue; }

        const std::string varname = amrex::toLower(trim(line.substr(0, equal_pos)));
        std::string values = trim(line.substr(equal_pos + 1));
        if (varname == "veg_dataset_description") {
            veg_categories.push_back(values);
            continue;
        }
        if (varname == "sltype") {
            soil_categories.push_back(values);
            continue;
        }

        std::replace(values.begin(), values.end(), ',', ' ');
        std::replace(values.begin(), values.end(), 'D', 'E');
        std::replace(values.begin(), values.end(), 'd', 'e');
        std::istringstream value_stream(values);
        std::string token;
        amrex::Vector<amrex::Real> var_values;
        while (value_stream >> token) {
            try {
                std::size_t parsed = 0;
                const amrex::Real value = std::stod(token, &parsed);
                if (parsed != token.size()) {
                    amrex::Abort("SLM: invalid value '" + token + "' for parameter '" +
                                 varname + "' in block '" + block_name + "'");
                }
                var_values.push_back(value);
            } catch (const std::exception&) {
                amrex::Abort("SLM: invalid value '" + token + "' for parameter '" +
                             varname + "' in block '" + block_name + "'");
            }
        }
        if (!var_values.empty()) {
            table[block_name].push_back(std::make_pair(varname, var_values));
        }
    }


    // Print out all the variables:
    amrex::Print() << "---------------------------------------------------" << std::endl;
    amrex::Print() << " SLM Parameter File '" << filename << "':" << std::endl;
    amrex::Print() << "---------------------------------------------------" << std::endl;
    for (auto &block : table) {
        amrex::Print() << "  '" << block.first << "':" << std::endl;
        for (auto &var : block.second) {
            amrex::Print() << "      - '" << var.first << "': [";
            for (int i = 0; i < var.second.size(); i++) {
                amrex::Print() << " " << var.second[i];
                if (i + 1 < var.second.size()) {
                    amrex::Print() << ",";
                }
            }
            amrex::Print() << " ]" << std::endl;
        }
    }
    amrex::Print() << "---------------------------------------------------" << std::endl;

    return table;
}

/**
 * Overwrites any default SLM variables from values set in the parameter file
 */
void SLM::init_from_params()
{
    // do nothing if not using the parameter file
    if (!use_param_file) {
        return;
    }

    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    // Get pointers to GPU soil parameter values
    const amrex::Real *d_param_poro = d_soil_params["maxsmc"]->data();
    const amrex::Real *d_param_theta_FC = d_soil_params["refsmc"]->data();
	const amrex::Real *d_param_theta_WP = d_soil_params["wltsmc"]->data();
	const amrex::Real *d_param_m_pot_sat = d_soil_params["satpsi"]->data();
    const amrex::Real *d_param_Bconst = d_soil_params["bb"]->data();
    const amrex::Real *d_param_ks = d_soil_params["satdk"]->data();

    // Get pointers to GPU vegetation parameter values from NoahMP table
    const amrex::Real *d_param_rs = nullptr;
    const amrex::Real *d_param_rgl = nullptr;
    const amrex::Real *d_param_xl = nullptr;
    const amrex::Real *d_param_hs = nullptr;
    const amrex::Real *d_param_nroot = nullptr;
    const amrex::Real *d_param_hvt = nullptr;
    const amrex::Real *d_param_z0mvt = nullptr;
    const amrex::Real *d_param_cbiom = nullptr;
    const amrex::Real *d_param_dleaf = nullptr;

    // Check if vegetation parameters are available
    bool has_rs = false, has_rgl = false, has_xl = false, has_hs = false;
    bool has_nroot = false, has_hvt = false, has_z0mvt = false;
    bool has_cbiom = false, has_dleaf = false;
    if (d_veg_params.find("rs") != d_veg_params.end()) {
        d_param_rs = d_veg_params.at("rs")->data();
        has_rs = true;
        amrex::Print() << " SLM: Using RS (minimum stomatal resistance) from NoahMP parameter file" << std::endl;
    }
    if (d_veg_params.find("rgl") != d_veg_params.end()) {
        d_param_rgl = d_veg_params.at("rgl")->data();
        has_rgl = true;
        amrex::Print() << " SLM: Using RGL (radiation stress parameter) from NoahMP parameter file" << std::endl;
    }
    if (d_veg_params.find("xl") != d_veg_params.end()) {
        d_param_xl = d_veg_params.at("xl")->data();
        has_xl = true;
        amrex::Print() << " SLM: Using XL (leaf/stem orientation index) from NoahMP parameter file" << std::endl;
    }
    if (d_veg_params.find("hs") != d_veg_params.end()) {
        d_param_hs = d_veg_params.at("hs")->data();
        has_hs = true;
        amrex::Print() << " SLM: Using HS (VPD sensitivity parameter) from NoahMP parameter file" << std::endl;
        amrex::Print() << " SLM: Note - VPD formula changed to match NoahMP: 1/(1+HS*VPD)" << std::endl;
    }
    if (d_veg_params.find("nroot") != d_veg_params.end()) {
        const auto& nroot_params = *d_veg_params.at("nroot");
        for (const amrex::Real value : nroot_params) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                value >= 0.0 && value <= m_nz_lsm && value == std::floor(value),
                "SLM: NoahMP NROOT must be an integer in [0, slm.nsoil]");
        }
        d_param_nroot = nroot_params.data();
        has_nroot = true;
        amrex::Print() << " SLM: Using NROOT from NoahMP parameter file for the BTR=1 soil-water stress factor" << std::endl;
    }
    if (d_veg_params.find("hvt") != d_veg_params.end()) {
        d_param_hvt = d_veg_params.at("hvt")->data();
        has_hvt = true;
        amrex::Print() << " SLM: Using HVT (canopy top height) from NoahMP parameter file" << std::endl;
    }
    if (d_veg_params.find("z0mvt") != d_veg_params.end()) {
        d_param_z0mvt = d_veg_params.at("z0mvt")->data();
        has_z0mvt = true;
        amrex::Print() << " SLM: Using Z0MVT (vegetation momentum roughness) from NoahMP parameter file" << std::endl;
    }
    if (d_veg_params.find("cbiom") != d_veg_params.end()) {
        d_param_cbiom = d_veg_params.at("cbiom")->data();
        has_cbiom = true;
        amrex::Print() << " SLM: Using CBIOM (canopy biomass heat capacity parameter) from NoahMP parameter file" << std::endl;
    }
    if (d_veg_params.find("dleaf") != d_veg_params.end()) {
        d_param_dleaf = d_veg_params.at("dleaf")->data();
        has_dleaf = true;
        amrex::Print() << " SLM: Using DLEAF (characteristic leaf dimension) from NoahMP parameter file" << std::endl;
    }

    for ( amrex::MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        amrex::Box bx2d = mfi.tilebox();

        auto landmask_arr = landmask.const_array(mfi);
        auto landtype_arr = landtype.const_array(mfi);
        auto vegetype_arr = vegetype.const_array(mfi);
        auto soiltype_arr = lsm_fab_vars[LsmVar_SLM::soiltype]->const_array(mfi);

        auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->array(mfi);
        auto theta_FC_arr = lsm_fab_vars[LsmVar_SLM::theta_FC]->array(mfi);
        auto theta_WP_arr = lsm_fab_vars[LsmVar_SLM::theta_WP]->array(mfi);
        auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->array(mfi);
        auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->array(mfi);
        auto ks_arr = lsm_fab_vars[LsmVar_SLM::ks]->array(mfi);

        // Get vegetation parameter arrays
        auto Rc_min_arr = Rc_min.array(mfi);
        auto Rgl_arr = Rgl.array(mfi);
        auto Khai_L_arr = Khai_L.array(mfi);
        auto hs_rc_arr = hs_rc.array(mfi);
        auto nroot_arr = nroot.array(mfi);
        auto ztop_arr = ztop.array(mfi);
        auto disp_hgt_arr = disp_hgt.array(mfi);
        auto z0_sfc_arr = z0_sfc.array(mfi);
        auto cbiom_arr = cbiom.array(mfi);
        auto dleaf_arr = dleaf.array(mfi);

        amrex::ParallelFor(bx2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
            if (landmask_arr(i, j, 0) == 1) {

                const int ltype = landtype_arr(i,j,0) - 1; // shift by one to match table index (i.e, types 1-20 -> index 0-19)
                const int stype = soiltype_arr(i,j,d_khi_lsm) - 1;

                // Set soil parameters from parameter file
                for (int k = d_khi_lsm; k >= d_klo_lsm; k--) {
                    poro_soil_arr(i, j, k) = d_param_poro[stype];
                    theta_FC_arr(i, j, k) = d_param_theta_FC[stype];
                    theta_WP_arr(i, j, k) = d_param_theta_WP[stype];
                    // NoahMP stores SATPSI as a positive magnitude in meters;
                    // SLM uses a negative matric potential in millimeters.
                    m_pot_sat_arr(i, j, k) = -1000.0 * d_param_m_pot_sat[stype];
                    Bconst_arr(i, j, k) = d_param_Bconst[stype];
                    // NoahMP SATDK is m/s; SLM soil hydraulics use mm/s.
                    ks_arr(i, j, k) = 1000.0 * d_param_ks[stype];
                }

                // Set vegetation parameters from NoahMP table if available
                // Only update for valid vegetation types (ltype >= 0)
                if (ltype >= 0 && vegetype_arr(i,j,0) == 1) {
                    if (has_rs && d_param_rs != nullptr) {
                        // RS = minimum stomatal resistance [s/m]
                        Rc_min_arr(i, j, 0) = d_param_rs[ltype];
                    }
                    if (has_rgl && d_param_rgl != nullptr) {
                        // RGL = radiation stress parameter [W/m2]
                        Rgl_arr(i, j, 0) = d_param_rgl[ltype];
                    }
                    if (has_xl && d_param_xl != nullptr) {
                        // XL = leaf/stem orientation index (dimensionless)
                        // Directly maps to Khai_L in SLM
                        Khai_L_arr(i, j, 0) = d_param_xl[ltype];
                    }
                    if (has_hs && d_param_hs != nullptr) {
                        // HS = VPD sensitivity parameter (dimensionless)
                        // NOTE: This only works with the hyperbolic VPD formula: 1/(1+HS*VPD)
                        // If using exponential formula exp(-hs*VPD), HS values are NOT compatible
                        hs_rc_arr(i, j, 0) = d_param_hs[ltype];
                    }
                    if (has_nroot && d_param_nroot != nullptr) {
                        nroot_arr(i, j, 0) = static_cast<int>(d_param_nroot[ltype]);
                    }
                    if (has_hvt && d_param_hvt != nullptr) {
                        ztop_arr(i, j, 0) = d_param_hvt[ltype];
                        disp_hgt_arr(i, j, 0) = 0.65 * ztop_arr(i, j, 0);
                    }
                    if (has_z0mvt && d_param_z0mvt != nullptr) {
                        z0_sfc_arr(i, j, 0) = d_param_z0mvt[ltype];
                    }
                    if (has_cbiom && d_param_cbiom != nullptr) {
                        cbiom_arr(i, j, 0) = d_param_cbiom[ltype];
                    }
                    if (has_dleaf && d_param_dleaf != nullptr) {
                        dleaf_arr(i, j, 0) = d_param_dleaf[ltype];
                    }
                }
            }
        });
    }

    // set a flag to indicate remaining soil variables should be recomputed
    params_updated = true;

    // recompute soil variables using updated parameters
    init_soil_vars();
}

/**
 * Updates the LAI + SAI based on the current simulation time and monthly values
 * from the LAI and SAI tables
 */
void SLM::UpdateLAI(const amrex::MFIter &mfi)
{
    if (interpolate_lai) {
        Box box = mfi.tilebox();
        box.makeSlab(2, 0);

        auto landmask_arr = landmask.const_array(mfi);
        auto landtype_arr = landtype.const_array(mfi);
        auto vegetype_arr = vegetype.const_array(mfi);
        auto LAI_arr = LAI.array(mfi);
        auto SAI_arr = SAI.array(mfi);

        // Update the day and month
        time_t timestamp = time_t(time + start_time);
        struct tm *timeinfo = gmtime(&timestamp);

        m_orbital_year = timeinfo->tm_year + 1900;
        m_orbital_mon  = timeinfo->tm_mon  + 1;
        m_orbital_day  = timeinfo->tm_mday;
        m_orbital_sec  = timeinfo->tm_hour*3600 + timeinfo->tm_min*60 + timeinfo->tm_sec;
        
        static constexpr double dpy[] = {0.0  ,  31.0,  59.0,  90.0, 120.0, 151.0,
                                        181.0, 212.0, 243.0, 273.0, 304.0, 334.0};
        bool leap = (m_orbital_year % 4 == 0 && (!(m_orbital_year % 100 == 0) || (m_orbital_year % 400 == 0))) ? true : false;
        m_calday = dpy[m_orbital_mon-1] + (m_orbital_day-1.0) + m_orbital_sec/86400.0;
        // add extra day if leap year
        if (leap) { m_calday += 1.0; }

        int curr_mon = m_orbital_mon - 1;
        int next_mon = curr_mon + 1;
        Real t0 = dpy[curr_mon];
        Real t1;
        if (next_mon >= 11) {
            next_mon = 0;
            t1 = 365.0;
        } else {
            t1 = dpy[next_mon];
        }

        const Real d_calday = m_calday;

        // Setup device interpolation points
        Real *d_lai_curr_ptr = d_lai_curr.data();
        Real *d_lai_next_ptr = d_lai_next.data();
        Real *d_sai_curr_ptr = d_sai_curr.data();
        Real *d_sai_next_ptr = d_sai_next.data();
        Gpu::copyAsync(Gpu::hostToDevice, lai_table[curr_mon].data(), lai_table[curr_mon].data()+num_landtypes, d_lai_curr_ptr);
        Gpu::copyAsync(Gpu::hostToDevice, lai_table[next_mon].data(), lai_table[next_mon].data()+num_landtypes, d_lai_next_ptr);
        Gpu::copyAsync(Gpu::hostToDevice, sai_table[curr_mon].data(), sai_table[curr_mon].data()+num_landtypes, d_sai_curr_ptr);
        Gpu::copyAsync(Gpu::hostToDevice, sai_table[next_mon].data(), sai_table[next_mon].data()+num_landtypes, d_sai_next_ptr);
        Gpu::streamSynchronize();

        // Update LAI = LAI + SAI
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) == 1) {
                //if (vegetype_arr(i, j, 0) == 0) { // todo: check, urban cells are not "vegetated", so LAI is incorrect there?
                if (landtype_arr(i,j,0) == 16) {
                    // baresoil point, so no LAI or SAI
                    LAI_arr(i,j,0) = 0.0;
                    SAI_arr(i,j,0) = 0.0;
                //} else if (vegetype_arr(i, j, 0) == 1) {
                } else {
                    // vegetation
                    const int ltype = landtype_arr(i,j,0) - 1;

                    const Real lai_x = d_lai_curr_ptr[ltype];
                    const Real lai_y = d_lai_next_ptr[ltype];

                    const Real sai_x = d_sai_curr_ptr[ltype];
                    const Real sai_y = d_sai_next_ptr[ltype];

                    LAI_arr(i,j,0) = linear_interp(t0, t1, d_calday, lai_x, lai_y);
                    SAI_arr(i,j,0) = linear_interp(t0, t1, d_calday, sai_x, sai_y);

                    LAI_arr(i,j,0) += SAI_arr(i,j,0);
                }
            }
        });

        UpdateLAIParameters(mfi);
    } else if (use_wrf_lai) {

        Box box = mfi.tilebox();
        box.makeSlab(2, 0);

        auto landmask_arr = landmask.const_array(mfi);
        auto landtype_arr = landtype.const_array(mfi);
        auto vegetype_arr = vegetype.const_array(mfi);

        auto veg_frac_arr = lsm_fab_vars[LsmVar_SLM::veg_frac]->array(mfi);
        auto veg_frac_min_arr = lsm_fab_vars[LsmVar_SLM::veg_frac_min]->const_array(mfi);
        auto veg_frac_max_arr = lsm_fab_vars[LsmVar_SLM::veg_frac_max]->const_array(mfi);

        auto IR_emis_veg = IR_emis_vege.array(mfi);

        auto LAI_arr = LAI.array(mfi);
        auto SAI_arr = SAI.array(mfi);

        auto const& d_params = d_param_table.const_table();

        // Update LAI = LAI + SAI
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) == 1) {

                const int ltype = landtype_arr(i,j,0) - 1; // shift by one to match table index (i.e, types 1-20 -> 0-19)

                AMREX_ASSERT(d_params(ltype, 0) == ltype + 1); // debug to make sure landtypes match

                if (ltype + 1 == 13) { // urban cells
                    veg_frac_arr(i,j,0) = d_params(ltype,1);
                }

                //    LAI min = col 8, lai max = col 9
                //  emiss min = col 10,    max = col 11
                // albedo min = col 12,    max = col 13

                if (veg_frac_arr(i,j,0) >= veg_frac_max_arr(i,j,0)) {
                    //emis_sfc_arr(i,j,0) = emissmax;
                    LAI_arr(i,j,0) = d_params(ltype, 9);
                } else if (veg_frac_arr(i,j,0) <= veg_frac_min_arr(i,j,0)) {
                    //emis_sfc_arr(i,j,0) = emissmin;
                    LAI_arr(i,j,0) = d_params(ltype, 8);
                } else {
                    if (veg_frac_max_arr(i,j,0) > veg_frac_min_arr(i,j,0)) {

                        Real interp_frac = (veg_frac_arr(i,j,0) - veg_frac_min_arr(i,j,0)) / (veg_frac_max_arr(i,j,0) - veg_frac_min_arr(i,j,0));
                        interp_frac = std::min(std::max(interp_frac, 0.0), 1.0); // bound between 0.0 and 1.0

                        // Scale emissivitiy and LAI between min/max by interp_frac
                        //emis_sfc_arr(i,j,0) = ( (1.0 - interp_frac) * d_params(ltype, 10)) + interp_frac * d_params(ltype, 11);

                        LAI_arr(i,j,0) = ( (1.0 - interp_frac) * d_params(ltype, 8)) + interp_frac * d_params(ltype, 9);
                    } else {
                        // emis_sfc_arr(i,j,0) = 0.5 * d_params(ltype,10) + 0.5 * d_params(ltype, 11);
                        LAI_arr(i,j,0) = 0.5 * d_params(ltype, 8) + 0.5 * d_params(ltype, 9);
                    }
                }

            }
        });

        UpdateLAIParameters(mfi);
    }
}

/**
 * Updates radiation parameters related to LAI
 */
void SLM::UpdateLAIParameters(const amrex::MFIter &mfi)
{
    Box box = mfi.tilebox();
    box.makeSlab(2, 0);

    auto IR_emis_vege_arr = IR_emis_vege.array(mfi);
    auto phi_1_arr = phi_1.array(mfi);
    auto phi_2_arr = phi_2.array(mfi);
    auto precip_extinc_arr = precip_extinc.array(mfi);
    auto mw_mx_arr = mw_mx.array(mfi);
    auto LAI_arr = LAI.array(mfi);
    auto BAI_arr = BAI.array(mfi);
    auto ztop_arr = ztop.array(mfi); 
    auto Khai_L_arr = Khai_L.array(mfi);
    auto landmask_arr = landmask.const_array(mfi);

    ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) == 1) {
            // set minimum LAI for vegetated land
            LAI_arr(i, j, 0) = std::max(LAI_arr(i, j, 0), 0.001);

            phi_1_arr(i, j, 0) = 0.5 - 0.633 * Khai_L_arr(i, j, 0) - 0.33 * (std::pow(Khai_L_arr(i, j, 0), 2));
            phi_2_arr(i, j, 0) = 0.877 * (1.0 - 2.0 * phi_1_arr(i, j, 0));
            IR_emis_vege_arr(i, j, 0) = 0.97 * (1.0 - std::exp(-1.0 * (phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0)) * LAI_arr(i, j, 0)));
            precip_extinc_arr(i, j, 0) = phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0);
            //mw_mx_arr(i, j, 0) = 0.1 * LAI_arr(i, j, 0);
            mw_mx_arr(i, j, 0) = 0.1 * LAI_arr(i, j, 0) + ztop_arr(i, j, 0) * std::pow(4.0 * PI * BAI_arr(i, j, 0) / 43560., 0.5);
        }
    });
}


/* Advance the solution with a simple explicit update (should use tridiagonal solve) */
void
SLM::AdvanceSLM ()
{
    // Expose for GPU copy
    Real dt = m_dt;
    Real dzInv = m_lsm_geom.InvCellSize(2);
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    net_rad.setVal(0.0, SLM_NetRad::net_sw1, 1, 0);
    net_rad.setVal(0.0, SLM_NetRad::net_sw2, 1, 0);
    net_rad.setVal(0.0, SLM_NetRad::net_rad1, 1, 0);
    net_rad.setVal(0.0, SLM_NetRad::net_rad2, 1, 0);

    // Soil temperature and moisture nudging
    for ( MFIter mfi(*lsm_fab_vars[LsmVar_SLM::tsurf], TileNoZ()); mfi.isValid(); ++mfi) {
        soil_nudging(mfi);
    }

    for ( MFIter mfi(landtype, TileNoZ()); mfi.isValid(); ++mfi) {
        auto box = mfi.tilebox();

        auto landmask_arr = landmask.const_array(mfi);

        auto LAI_arr = LAI.const_array(mfi);
        auto precip_extinc_arr = precip_extinc.const_array(mfi);

        auto mw_arr = mw.array(mfi);
        auto mws_arr = mws.const_array(mfi);
        auto mw_mx_arr = mw_mx.const_array(mfi);

        auto mw_inc_arr = mw_inc.array(mfi);

        auto t_cas_arr = t_cas.array(mfi);
        auto q_cas_arr = q_cas.array(mfi);
        
        auto t_sfc_arr = t_sfc.array(mfi);
        auto q_sfc_arr = q_sfc.array(mfi);
        
        auto q_gr_arr = q_gr.array(mfi);

        auto t_canop_arr = t_canop.array(mfi);
        auto t_skin_arr = t_skin.const_array(mfi);
        auto t_ground_skin_arr = t_ground_skin.const_array(mfi);

        auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
        auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);

        auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
        auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);

        auto ustar_arr = lsm_fab_vars[LsmVar_SLM::ustar]->const_array(mfi);
        auto tstar_arr = lsm_fab_vars[LsmVar_SLM::tstar]->array(mfi);
        auto qstar_arr = lsm_fab_vars[LsmVar_SLM::qstar]->array(mfi);

        auto BAI_arr = BAI.const_array(mfi);
        auto ztop_arr = ztop.const_array(mfi);
        auto SAI_arr = SAI.const_array(mfi);
        auto cbiom_arr = cbiom.const_array(mfi);

        auto shf_canop_arr = shf_canop.array(mfi);
        auto shf_soil_arr = shf_soil.array(mfi);
        auto shf_air_arr = shf_air.array(mfi);

        auto lhf_canop_arr = lhf_canop.array(mfi);
        auto lhf_soil_arr = lhf_soil.array(mfi);
        auto lhf_air_arr = lhf_air.array(mfi);
        auto evp_canop_arr = evp_canop.array(mfi);
        auto evp_air_arr = evp_air.array(mfi);

        auto vegetype_arr = vegetype.const_array(mfi);
        auto vege_YES_arr = vege_YES.const_array(mfi);

        auto cp_vege_arr = cp_vege.array(mfi);
        
        auto sdew_arr = sdew.array(mfi);

        auto tref_arr  = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
        auto ur_arr  = lsm_fab_vars[LsmVar_SLM::uref]->const_array(mfi);
        auto vr_arr  = lsm_fab_vars[LsmVar_SLM::vref]->const_array(mfi);
        auto dref_arr  = lsm_fab_vars[LsmVar_SLM::dref]->const_array(mfi);
        auto qref_arr  = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
        auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);
        auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

        auto tsurf_arr = lsm_fab_vars[LsmVar_SLM::tsurf]->array(mfi);
        auto tveg_arr = lsm_fab_vars[LsmVar_SLM::tv]->array(mfi);
        auto mveg_arr = lsm_fab_vars[LsmVar_SLM::mv]->array(mfi);

        auto r_a_arr = r_a.const_array(mfi);
        auto r_b_arr = r_b.const_array(mfi);
        auto r_c_arr = r_c.const_array(mfi);
        auto r_d_arr = r_d.const_array(mfi);
        auto r_soil_arr = r_soil.const_array(mfi);

        auto flbu_arr  = lsm_fab_vars[LsmVar_SLM::flbu]->array(mfi);
        auto flbv_arr  = lsm_fab_vars[LsmVar_SLM::flbv]->array(mfi);
        auto flbq_arr  = lsm_fab_vars[LsmVar_SLM::flbq]->array(mfi);
        auto flbt_arr  = lsm_fab_vars[LsmVar_SLM::flbt]->array(mfi);
        auto prsfc_arr  = lsm_fab_vars[LsmVar_SLM::prsfc]->array(mfi);

        auto net_rad_arr = net_rad.const_array(mfi);
        auto wet_canop_arr = wet_canop.const_array(mfi);

        // TODO: Copies for MOST
        auto fluxq_arr = lsm_fab_flux[LsmFlux_SLM::q_flux]->array(mfi);
        auto fluxt_arr = lsm_fab_flux[LsmFlux_SLM::t_flux]->array(mfi);
        auto tau13_arr = lsm_fab_flux[LsmFlux_SLM::tau13]->array(mfi);
        auto tau23_arr = lsm_fab_flux[LsmFlux_SLM::tau23]->array(mfi);
        auto olen_arr  = lsm_fab_flux[LsmFlux_SLM::olen]->array(mfi);
    
        auto slm_diag_arr = slm_diag.array(mfi);

        // Update LAI and SAI based on current month
        UpdateLAI(mfi);

        if (radiation_scheme == RadiationScheme::SLM) {
            // Calculate net radiation absorbed by canopy and soil surface
            // Old SLM radiation
            radiative_fluxes(mfi);
        } else {
            // NOAHMP radiation calculation
            radiation_noahmp(mfi);
        }

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) != 1) {
                return;
            }

            amrex::Real precip = 0.0;
            amrex::Real drain = 0.0;
            amrex::Real mws_inc = 0.0;
            amrex::Real cp_vege_tot, fh;
            // precipitation interception rate at canoppy
            // For baresoil, precip = 0, as LAI = 0
            precip = precip_array(i, j, 0)*(1.0 - std::exp(-1.0 * precip_extinc_arr(i, j, 0)*LAI_arr(i, j, 0)));

            if (mw_arr(i, j, 0) < mw_mx_arr(i, j, 0))
            {
                drain = 0.0;
            }
            else if(mw_arr(i, j, 0) > mw_mx_arr(i, j, 0))
            {
                // when water holding storage exceeds its maximum, no precipitation is intercepted
                drain = precip;

                // excess water storage gets drained from mw [kg/m^2]
                drain += std::max(mw_arr(i, j, 0) - mw_mx_arr(i, j, 0), 0.0) / dt; // mm/s
            }

            prsfc_arr(i, j, 0) = precip_array(i, j, 0) - precip + drain;

            // Update output variables
            slm_diag_arr(i, j, 0, SLM_Diag::precip) = precip;
            
            // Original SLM formula:
            // cp_vege_arr(i, j, 0) =
            //     (LAI_arr(i, j, 0) * leaf_thickness * 0.001
            //      + ztop_arr(i, j, 0) * BAI_arr(i, j, 0) / 43560.) * 900. * 2800.;

            // NoahMP dry canopy heat capacity: CBIOM * effective VAI * volumetric water heat capacity.
            const amrex::Real vai = std::min(6.0, LAI_arr(i, j, 0) + SAI_arr(i, j, 0));
            cp_vege_arr(i, j, 0) = cbiom_arr(i, j, 0) * vai * cp_water * vege_YES_arr(i, j, 0);
            cp_vege_tot = cp_vege_arr(i, j, 0) + mw_arr(i, j, 0) * 1.e-3 * cp_water;

            // Add to vegetiation moisture increment from SLM::vapor_fluxes()
            mw_inc_arr(i, j, 0) += dt * (precip - drain);

            // Update vegetation moisture storage
            mw_arr(i, j, 0) += mw_inc_arr(i, j, 0);
            
            // from gSAM-SLM
            // intercepted precip cools the canopy:
            // temperature of intercepted rain is the same as reference level
            // note it cools even when water storage is full as old water on leaves is replaced by new rain water
            if (vegetype_arr(i, j, 0) == 1)
            {
                t_canop_arr(i, j, 0) = (cp_vege_tot * t_canop_arr(i, j, 0) + tref_arr(i, j, 0) * precip * dt * 1.e-3 * cp_water) / (cp_vege_tot + precip * dt * 1.e-3 * cp_water);
            }

            // Note:
            // Assign appropriate "surface level" values for each land type, for the calculation  
            //of surface turbulent fluxes
            //For baresoil,   surface level = soil surface
            //For vegetation, surface level = canopy level
            if (vegetype_arr(i, j, 0) == 1)
            {
                t_sfc_arr(i, j, 0) = t_cas_arr(i, j, 0);
                q_sfc_arr(i, j, 0) = q_cas_arr(i, j, 0);
            }
            else
            {
                t_sfc_arr(i, j, 0) = soilt_arr(i, j, d_khi_lsm);
                // Specific humidity at top soil
                if (soilt_arr(i, j, d_khi_lsm) > tfriz)
                {
                    erf_qsatw(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr_arr(i, j, 0));
                    sdew_arr(i, j, 0) = 1.;
                    if (mws_arr(i, j, 0) == 0.0)
                    {
                        fh = fh_calc(soilt_arr(i, j, d_khi_lsm), m_pot_sat_arr(i, j, d_khi_lsm), soilw_arr(i, j, d_khi_lsm), Bconst_arr(i, j, d_khi_lsm));
                        if ( fh > 0.99)
                        { 
                            sdew_arr(i, j, 0) = 1.;
                        }
                        else
                        { 
                            sdew_arr(i, j, 0) = 0.;
                        }   
                        q_gr_arr(i, j, 0) *= fh;
                    } 
                }
                else
                {
                    erf_qsati(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr_arr(i, j, 0));
                    sdew_arr(i, j, 0) = 1.;
                }
                q_sfc_arr(i, j, 0) = q_gr_arr(i, j, 0);
            }                   
        });

        // Calculate turbulent transfer coeff between reference level and surface
        transfer_coeff(mfi);

        // Calculate aerodynamic resistances + stomatal resistance
        resistances(mfi);

        //Subcycle in time to avoid swings of canopy temperature
        // because of small heat capacity of vegetation in curtain places and seasons
        // Of course, the best would be to use implicit scheme, but for now,
        // subcycling seems like a good fix.
        //from gSAM-SLM 
        fluxes_canopy(mfi); // compute both sHF and LHF of canopy

        solve_ground_skin_temperature(mfi);
        
        // Calculate soil moisture increment
        soil_water(mfi);

        // Advance the finite-volume soil using only the final skin-to-soil flux.
        soil_temperature(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (landmask_arr(i, j, 0) != 1) {
                return;
            }

            amrex::Real cp_vege_tot;
            amrex::Real q_gr;

            // SAM rhow[nz] = air density at vertical velocity levels, kg/m^3
            amrex::Real rhow = dref_arr(i, j, 0); // TODO: double check this
            amrex::Real cond_heat, cond_href, cond_hcnp, cond_hundercnp;
            amrex::Real cond_vapor, cond_vref, cond_vcnp, cond_vundercnp;
            amrex::Real cvh = 2 * LAI_arr(i, j, 0) / r_b_arr(i, j, 0);
            if (vegetype_arr(i, j, 0) == 1) {
                // Compute diagnostic variables at canopy air space level
                //   Calculate heat conductances - non-zero only for canopy land type
                cond_heat = 1.0 / r_a_arr(i, j, 0) + cvh + 1.0 / r_d_arr(i, j, 0);
                cond_href = 1.0 / r_a_arr(i, j, 0) / cond_heat;
                cond_hcnp = cvh / cond_heat;
                cond_hundercnp = 1.0 / r_d_arr(i, j, 0) / cond_heat;

                //  Calculate vapor conductances
                cond_vapor = 1.0 / r_a_arr(i, j, 0) + wet_canop_arr(i, j, 0) * LAI_arr(i, j, 0) / r_b_arr(i, j, 0) + (1.0 - wet_canop_arr(i, j, 0)) * LAI_arr(i, j, 0) /(r_b_arr(i, j, 0) + r_c_arr(i, j, 0)) + 1.0 / (r_d_arr(i, j, 0) + r_soil_arr(i, j, 0) + r_litter);
                cond_vref = 1.0 / r_a_arr(i, j, 0) / cond_vapor; 
                cond_vcnp = (wet_canop_arr(i, j, 0) * LAI_arr(i, j, 0) / r_b_arr(i, j, 0) / cond_vapor + (1.0 - wet_canop_arr(i, j, 0)) * LAI_arr(i, j, 0) / (r_b_arr(i, j, 0) + r_c_arr(i, j, 0)) / cond_vapor);
                cond_vundercnp = 1.0 / (r_d_arr(i, j, 0) + r_soil_arr(i, j, 0) + r_litter) / cond_vapor;
            }
            else 
            {
                cond_heat = 1./r_a_arr(i, j, 0);
                cond_href = 1.;
                cond_hcnp = 0.;
                cond_hundercnp = 0.;

                cond_vapor = 1./(r_a_arr(i, j, 0) + r_soil_arr(i, j, 0));
                cond_vref = 1.;
                cond_vcnp = 0.;
                cond_vundercnp = 0.;
            }   

            t_cas_arr(i, j, 0) = tref_arr(i, j, 0)*cond_href + t_canop_arr(i, j, 0)*cond_hcnp + t_ground_skin_arr(i, j, 0)*cond_hundercnp;
            
            amrex::Real qsat_canop;
            if (t_canop_arr(i, j, 0) >= tfriz)
            {
                erf_qsatw(t_canop_arr(i, j, 0), pref_arr(i, j, 0), qsat_canop);
            }
            else
            {
                erf_qsati(t_canop_arr(i, j, 0), pref_arr(i, j, 0), qsat_canop);
            }

            q_gr = q_gr_arr(i, j, 0);
            q_cas_arr(i, j, 0) = qref_arr(i, j, 0) * cond_vref + qsat_canop*cond_vcnp + q_gr*cond_vundercnp;

            // Output variables
            tveg_arr(i, j, 0) = t_canop_arr(i, j, 0);
            tveg_arr(i, j, d_khi_lsm) = t_canop_arr(i, j, 0);
            mveg_arr(i, j, 0) = mw_arr(i, j, 0);
            mveg_arr(i, j, d_khi_lsm) = mw_arr(i, j, 0);
            tsurf_arr(i, j, d_khi_lsm) = t_skin_arr(i, j, 0); // TODO: ts in SLM is input and output - check how this should be coupled back to ERF
            tsurf_arr(i, j, 0) = tsurf_arr(i, j, d_khi_lsm); // make sure this is set at k=0 for radiation coupling
            flbq_arr(i, j, d_khi_lsm) = evp_air_arr(i, j, 0) / rhow; // kg/kg m/s
            flbt_arr(i, j, d_khi_lsm) = shf_air_arr(i, j, 0) / (Cp_d*rhow); // Km/s

            qstar_arr(i, j, 0) = -1.0 * flbq_arr(i, j, d_khi_lsm) / ustar_arr(i, j, 0);
            qstar_arr(i, j, d_khi_lsm) = qstar_arr(i, j, 0);

            // TODO: fix - copies SLM boundary into flux array for MOST
            fluxq_arr(i,j,0) = rhow * flbq_arr(i, j, d_khi_lsm);
            fluxt_arr(i,j,0) = rhow * flbt_arr(i, j, d_khi_lsm);
            tau13_arr(i,j,0) = rhow * flbu_arr(i, j, d_khi_lsm);
            tau23_arr(i,j,0) = rhow * flbv_arr(i, j, d_khi_lsm);

            amrex::Real tvm = getThgivenRandT(rhow, tref_arr(i,j,0), R_d / Cp_d, qref_arr(i,j,0)) * (1.0 + 0.61 * qref_arr(i,j,0));
            olen_arr(i,j,0) = -ustar_arr(i,j,0) * ustar_arr(i,j,0) * ustar_arr(i,j,0) * tvm / (KAPPA * CONST_GRAV * flbt_arr(i,j,d_khi_lsm));
            
        });
    }

    lsm_fab_vars[LsmVar_SLM::tsurf]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::ustar]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::tstar]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::qstar]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::flbu]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::flbv]->FillBoundary(m_geom.periodicity());

    // TODO: fix;
    lsm_fab_flux[LsmFlux_SLM::q_flux]->FillBoundary(m_geom.periodicity());
    lsm_fab_flux[LsmFlux_SLM::t_flux]->FillBoundary(m_geom.periodicity());
    lsm_fab_flux[LsmFlux_SLM::tau13]->FillBoundary(m_geom.periodicity());
    lsm_fab_flux[LsmFlux_SLM::tau23]->FillBoundary(m_geom.periodicity());
    lsm_fab_flux[LsmFlux_SLM::olen]->FillBoundary(m_geom.periodicity());
}

void SLM::radiative_fluxes(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);

    auto swdsvisxyref_arr  = lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->const_array(mfi);
    auto swdsnirxyref_arr  = lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->const_array(mfi);
    auto swdsvisdxyref_arr = lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->const_array(mfi);
    auto swdsnirdxyref_arr = lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->const_array(mfi);
    auto lwref_arr         = lsm_fab_vars[LsmVar_SLM::lwref]->const_array(mfi);
    auto coszrsxy_arr      = lsm_fab_vars[LsmVar_SLM::coszrsxy]->const_array(mfi);

    auto phi_1_arr = phi_1.const_array(mfi);
    auto phi_2_arr = phi_2.const_array(mfi);
    auto LAI_arr   = LAI.const_array(mfi);

    auto albedovis_v_arr = albedovis_v.const_array(mfi);
    auto albedonir_v_arr = albedonir_v.const_array(mfi);
    auto albedovis_s_arr = albedovis_s.const_array(mfi);
    auto albedonir_s_arr = albedonir_s.const_array(mfi);

    auto IR_emis_vege_arr = IR_emis_vege.const_array(mfi);
    auto IR_emis_soil_arr = IR_emis_soil.const_array(mfi);
    auto IR_emis_grnd_arr = IR_emis_grnd.const_array(mfi);
    auto t_canop_arr      = t_canop.const_array(mfi);
    auto soilw_arr        = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto soilt_arr        = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto veg_frac_arr     = lsm_fab_vars[LsmVar_SLM::veg_frac]->const_array(mfi);

    // Combined (veg+soil) surface emissivity and albedos for coupling to radiation model
    auto emis_sfc_arr     = lsm_fab_vars[LsmVar_SLM::emis_sfc]->array(mfi);
    auto alb_nir_sfc_arr  = lsm_fab_vars[LsmVar_SLM::alb_nir_sfc]->array(mfi);
    auto alb_vis_sfc_arr  = lsm_fab_vars[LsmVar_SLM::alb_vis_sfc]->array(mfi);
    auto alb_nir_sfc_diff_arr  = lsm_fab_vars[LsmVar_SLM::alb_nir_sfc_diff]->array(mfi);
    auto alb_vis_sfc_diff_arr  = lsm_fab_vars[LsmVar_SLM::alb_vis_sfc_diff]->array(mfi);

    auto t_skin_arr  = t_skin.array(mfi);
    auto net_rad_arr = net_rad.array(mfi);

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        amrex::Real fdn1; // downwelling flux on top
        amrex::Real fdn2; // transmitted flux (downwelling flux below)
        amrex::Real fup1; // upwelling flux from top
        amrex::Real fup2; // upwelling flux from below

        amrex::Real ka; // optical depth
        amrex::Real explai, explai0, wetfactor;

        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        // ===================================================
        // Compute shortwave radiation transfer between land surface and reference level
        // ===================================================
        if (coszrsxy_arr(i, j, 0) > 0.0)
        {
            // Optical depth of the direct beam per unit leaf area
            ka = phi_1_arr(i, j, 0) / std::max(0.01, coszrsxy_arr(i, j, 0)) + phi_2_arr(i, j, 0);
            explai = exp(-ka*LAI_arr(i, j, 0)); // for direct radiation

            // Optical depth of the diffuse beam per unit leaf area
            ka = phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0);
            explai0 = exp(-ka*LAI_arr(i, j, 0)); // for diffuse radiation

            // net_rad(1) = net absorbed shortwave radiation by canopy
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) += swdsvisxyref_arr(i, j, 0)*(1.0 - albedovis_v_arr(i, j, 0)*(1.0 - explai)-explai);
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) += swdsvisdxyref_arr(i, j, 0)*(1.0 - albedovis_v_arr(i, j, 0)*(1.0 - explai0)-explai);
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) += swdsnirxyref_arr(i, j, 0)*(1.0 - albedonir_v_arr(i, j, 0)*(1.0 - explai)-explai);
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) += swdsnirdxyref_arr(i, j, 0)*(1.0 - albedonir_v_arr(i, j, 0)*(1.0 - explai0)-explai);

            net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1) = swdsvisxyref_arr(i, j, 0) + swdsvisdxyref_arr(i, j, 0) + swdsnirxyref_arr(i, j, 0) + swdsnirdxyref_arr(i, j, 0);
            net_rad_arr(i, j, 0, SLM_NetRad::net_swup1) = net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) - net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1) * (1. - explai);

            // net_rad(2) = net absorbed shortwave radiation by soil
            wetfactor = 1.0 - 0.5*soilw_arr(i, j, d_khi_lsm); // soil wetness factor: assume that wet soil is twice as dark
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsvisxyref_arr(i, j, 0)*(1.0 - albedovis_s_arr(i, j, 0)*wetfactor)*explai;
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsvisdxyref_arr(i, j, 0)*(1.0 - albedovis_s_arr(i, j, 0)*wetfactor)*explai0;
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsnirxyref_arr(i, j, 0)*(1.0 - albedonir_s_arr(i, j, 0)*wetfactor)*explai;
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsnirdxyref_arr(i, j, 0)*(1.0 - albedonir_s_arr(i, j, 0)*wetfactor)*explai0;

            net_rad_arr(i, j, 0, SLM_NetRad::net_swdn2) = net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1)*explai;
            net_rad_arr(i, j, 0, SLM_NetRad::net_swup2) = net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) - net_rad_arr(i, j, 0, SLM_NetRad::net_swdn2);

            // Albedo is computed as alb = alb_v*(1-exp(-kLAI))+alb_s*exp(-kLAI)
            Real alb_nir_veg_dir = albedonir_v_arr(i, j, 0)*(1.0 - explai);
            Real alb_nir_veg_dif = albedonir_v_arr(i, j, 0)*(1.0 - explai0);
            Real alb_vis_veg_dir = albedovis_v_arr(i, j, 0)*(1.0 - explai);
            Real alb_vis_veg_dif = albedovis_v_arr(i, j, 0)*(1.0 - explai0);

            Real alb_nir_soil_dir = albedonir_s_arr(i, j, 0)*wetfactor*explai;
            Real alb_nir_soil_dif = albedonir_s_arr(i, j, 0)*wetfactor*explai0;
            Real alb_vis_soil_dir = albedovis_s_arr(i, j, 0)*wetfactor*explai;
            Real alb_vis_soil_dif = albedovis_s_arr(i, j, 0)*wetfactor*explai0;

            alb_nir_sfc_arr(i, j, 0) = alb_nir_veg_dir + alb_nir_soil_dir;
            alb_vis_sfc_arr(i, j, 0) = alb_vis_veg_dir + alb_vis_soil_dir;
            alb_nir_sfc_diff_arr(i, j, 0) = alb_nir_veg_dif + alb_nir_soil_dif;
            alb_vis_sfc_diff_arr(i, j, 0) = alb_vis_veg_dif + alb_vis_soil_dif;

            alb_nir_sfc_arr(i, j, d_khi_lsm) = alb_nir_sfc_arr(i, j, 0);
            alb_vis_sfc_arr(i, j, d_khi_lsm) = alb_vis_sfc_arr(i, j, 0);
            alb_nir_sfc_diff_arr(i, j, d_khi_lsm) = alb_nir_sfc_diff_arr(i, j, 0);
            alb_vis_sfc_diff_arr(i, j, d_khi_lsm) = alb_vis_sfc_diff_arr(i, j, 0);


            // Store net absorbed SW
            net_rad_arr(i, j, 0, SLM_NetRad::net_sw1) = net_rad_arr(i, j, 0, SLM_NetRad::net_rad1);
            net_rad_arr(i, j, 0, SLM_NetRad::net_sw2) = net_rad_arr(i, j, 0, SLM_NetRad::net_rad2);
        }

        // ===================================================
        // Longwave radiation
        // ===================================================

        // ===================================================
        // tir: Emitted thermal infrared radiation
        // ===================================================
        //  Note: for no vegetation: IR_trans becomes zero => tir(1) automatically becomes zero
        net_rad_arr(i, j, 0, SLM_NetRad::tir1) = IR_emis_vege_arr(i, j, 0)*sigma*(std::pow(t_canop_arr(i, j, 0), 4));
        net_rad_arr(i, j, 0, SLM_NetRad::tir2) = IR_emis_grnd_arr(i, j, 0)*sigma*(std::pow(soilt_arr(i, j, d_khi_lsm), 4));

        // ===================================================
        // downwelling LW on canopy top: input
        // ===================================================
        fdn1 = lwref_arr(i, j, 0);
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwdn1) = fdn1;

        // ===================================================
        // downwelling LW below canopy layer
        // ===================================================
        //  Note:
        //    below canopy layer: incoming LW (fdn1) that is not absorbed by canopy
        //                        + emitted thermal IR by canopy (tir(1)) toward soil surface
        //    with no canopy: fdn2 is computed to be fdn1
        //    (1 - IR_emis) = area of canopy gap (skyview factor)
        // ===================================================
        fdn2 = (1.0 - IR_emis_vege_arr(i, j, 0))*fdn1 + net_rad_arr(i, j, 0, SLM_NetRad::tir1);
        net_rad_arr(i, j, 0, SLM_NetRad::net_lw1) = fdn1 - fdn2;

        // ===================================================
        //  Note: At this stage,
        //   net_rad(1) = net absorbed SW by canopy
        //           + downwelling LW on canopy top
        //           - transmitted LW through canopy layer
        //           - emitted TIR toward soil surface
        // ===================================================
        net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) += fdn1 - fdn2;

        // ===================================================
        // downwelling LW for soil surface
        // ===================================================
        fdn1 = fdn2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwdn2) = fdn1;

        // no fluxes below topsoil
        fdn2 = 0.0;
        fup2 = 0.0;

        // ===================================================
        // Note:
        //  Emitted LW from topsoil = emitted tir from topsoil + portion of incoming LW that is reflected back toward canopy
        //  IR_emis_soil is set to 1.0
        // ===================================================
        fup1 = (net_rad_arr(i, j, 0, SLM_NetRad::tir2) + (1.0 - IR_emis_grnd_arr(i, j, 0))*fdn1);
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup2) = fup1;

        // ===================================================
        // Note: At this stage,
        //  net_rad(2) = net absorbed SW by soil surface
        //             + net absorbed LW by soil surface
        // ===================================================
        net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += fdn1 - fdn2 - fup1 + fup2;

        // net absorbed LW by soil surface
        net_rad_arr(i, j, 0, SLM_NetRad::net_lw2) = fdn1 - fdn2 - fup1 + fup2;

        // ===================================================
        // Incoming LW from below canopy
        //  Note: Incoming LW from below canopy = upwelling flux at topsoil
        // ===================================================
        fup2 = fup1;

        // ===================================================
        // Upwelling LW from canopy top
        //  Note: fup1 = portion of fup2 that is not absorbed + tir emitted from canopy
        //        for no canopy: fup1 = fup2
        // ===================================================
        fup1 = (1.0 - IR_emis_vege_arr(i,j,0))*fup2 + net_rad_arr(i, j, 0, SLM_NetRad::tir1);

        emis_sfc_arr(i, j, 0) = IR_emis_vege_arr(i,j,0) * veg_frac_arr(i, j, d_khi_lsm) + IR_emis_soil_arr(i,j,0) * (1.0 - veg_frac_arr(i, j, d_khi_lsm));

        emis_sfc_arr(i, j, d_khi_lsm) = emis_sfc_arr(i,j,0);

        // total upward LW from surface (for canopy cover- from canopy top, for no canopy - from soil surface)
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup1) = fup1;
        t_skin_arr(i, j, 0) = std::pow(fup1 / sigma, 0.25);

        net_rad_arr(i, j, 0, SLM_NetRad::net_lw1) += fup2 - fup1;


        // ===================================================
        //  Note: At this stage,
        //   net_rad(1) = net absorbed SW by canopy
        //           + downwelling LW on canopy top
        //           - transmitted LW through canopy layer (downward direction)
        //           - emitted TIR from canopy toward soil surface
        //           - emitted TIR from canopy toward atmosphere
        //           + upwelling LW from topsoil
        //           - transmitted LW through canopy layer (upward direction)
        // ===================================================
        net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) += fup2 - fup1;
    });
}

void SLM::transfer_coeff(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const Real dt = m_dt;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);

    auto t_cas_arr = t_cas.const_array(mfi);
    auto q_cas_arr = q_cas.const_array(mfi);
    
    auto t_sfc_arr = t_sfc.const_array(mfi);
    auto q_sfc_arr = q_sfc.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto vegetype_arr = vegetype.const_array(mfi);

    auto mws_arr = mws.const_array(mfi);

    auto disp_hgt_arr = disp_hgt.const_array(mfi);
    auto z0_sfc_arr = z0_sfc.const_array(mfi);

    auto qr_arr = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
    auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
    auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);

    auto ustar_arr = lsm_fab_vars[LsmVar_SLM::ustar]->array(mfi);
    auto tstar_arr = lsm_fab_vars[LsmVar_SLM::tstar]->array(mfi);

    auto tref_arr  = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
    auto ur_arr  = lsm_fab_vars[LsmVar_SLM::uref]->const_array(mfi);
    auto vr_arr  = lsm_fab_vars[LsmVar_SLM::vref]->const_array(mfi);
    auto dref_arr  = lsm_fab_vars[LsmVar_SLM::dref]->const_array(mfi);
    auto qref_arr  = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
    auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);
    auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

    auto r_a_arr = r_a.array(mfi);
    
    auto zref_arr = zrefxy.array(mfi);

    auto flbu_arr  = lsm_fab_vars[LsmVar_SLM::flbu]->array(mfi);
    auto flbv_arr  = lsm_fab_vars[LsmVar_SLM::flbv]->array(mfi);

    constexpr amrex::Real xsim = -1.574;
    constexpr amrex::Real xsih = -0.465;
    const amrex::Real xm = sqrt(sqrt(1.0 - 16.0*xsim));
    const amrex::Real xh = sqrt(sqrt(1.0 - 16.0*xsih));

    constexpr amrex::Real errormax = 0.01;
    constexpr amrex::Real kk = 0.4;
    constexpr int nitermax = 10;

    auto xx = [] AMREX_GPU_DEVICE (const amrex::Real &yy) -> amrex::Real {
        return sqrt(sqrt((1.0 - 16.0 * yy)));
    };
    // unstable: -1.574<xsi<0
    auto psim1 = [] AMREX_GPU_DEVICE (const amrex::Real &x, const amrex::Real &x0) -> amrex::Real {
        return 2.0 * log((1.0 + x) / (1.0 + x0)) + log((1.0 + x*x) / (1.0 + x0*x0)) - 2.0*(atan(x)-atan(x0));
    };
    auto psih1 = [] AMREX_GPU_DEVICE (const amrex::Real &x, const amrex::Real &x0) -> amrex::Real {
        return 2.0 * log((1.0 + x*x) / (1.0 + x0*x0));
    };
    // very unstable: xsi < -1.574
    auto const psim2 = [=] AMREX_GPU_DEVICE (const amrex::Real &xsi, const amrex::Real &xsim0, const amrex::Real &xm0) -> amrex::Real {
        return log(xsim / xsim0) - psim1(xm, xm0) + 1.14*(std::pow(-xsi, 0.3333) - std::pow(-xsim, 0.3333));
    };
    auto const psih2 = [=] AMREX_GPU_DEVICE (const amrex::Real &xsi, const amrex::Real &xsih0, const amrex::Real &xh0) -> amrex::Real {
        return log(xsih / xsih0) - psih1(xh, xh0) + 0.8*(std::pow(-xsi, 0.3333) - std::pow(-xsih, 0.3333));
    };
    // stable: 0 < xsi < 1
    auto psim3 = [] AMREX_GPU_DEVICE (const amrex::Real &xsi, const amrex::Real &xsim0) -> amrex::Real {
        return -5.0 * (xsi - xsim0);
    };
    auto psih3 = [] AMREX_GPU_DEVICE (const amrex::Real &xsi, const amrex::Real &xsih0) -> amrex::Real {
        return -5.0 * (xsi - xsih0);
    };
    // very stable: 0 < xsi < 1
    auto psim4 = [] AMREX_GPU_DEVICE (const amrex::Real &xsi, const amrex::Real &xsim0) -> amrex::Real {
        return log(std::pow(xsi, 5) / xsim0) + 5.0 * (1.0 - xsim0) + xsi - 1.0;
    };
    auto psih4 = [] AMREX_GPU_DEVICE (const amrex::Real &xsi, const amrex::Real &xsih0) -> amrex::Real {
        return log(std::pow(xsi, 5) / xsih0) + 5.0 * (1.0 - xsih0) + xsi - 1.0;
    };


    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        // Inputs:
        // ts = t_sfc
        // th = tref
        // qh = qr
        // qs = q_sfc
        // h = zref = height of ref level
        // z0 = z0_sfc = surface roughness length
        // disp = disp_hgt

        amrex::Real tsp = t_sfc_arr(i, j, 0) * std::pow(1000.0/pref_arr(i, j, 0), rair / cp);
        amrex::Real thp = tref_arr(i, j, 0) * std::pow(1000.0 / pref_arr(i, j, 0), rair / cp);

        amrex::Real vel;
        // Add additional velocity depending on the stratification
        if ((thp - tsp) >= 0.0)
        {
            vel = std::max(0.5, sqrt(std::pow(ur_arr(i, j, 0), 2) + std::pow(vr_arr(i, j, 0), 2)));
        }
        else
        {
            vel = std::max(0.5, sqrt(std::pow(ur_arr(i, j, 0), 2) + std::pow(vr_arr(i, j, 0), 2)));
        }
        const Real d_zref = zref_arr(i,j,0);

        amrex::Real r = 9.81 / tsp * (thp * (1.0 + epsv * qr_arr(i, j, 0)) - tsp * (1.0 + epsv * q_sfc_arr(i, j, 0))) * (d_zref - disp_hgt_arr(i, j, 0)) / (vel*vel);
        r = std::max(-10.0, std::min(r, 0.5)); // cap r for stability of iterations
        // initial guess
        amrex::Real xsi, fm, fh, xsi1;
        amrex::Real xsim0, xsih0;
        amrex::Real mom_trans_coef, heat_trans_coef;

        //adjust z0 for high wind - gSAM-SLM
        amrex::Real z0_adjusted = z0_sfc_arr(i, j, 0) * std::pow(1. + vel/10., -0.6);
        amrex::Real zt0 = 0.135 * z0_adjusted; // roughness length for scalars: ln(z0/zt0)=2 (from Garatte BL textbook, p.93) 
        
        amrex::Real z0h = z0_adjusted / (d_zref - disp_hgt_arr(i, j, 0));
        // make sure (h - disp) is not negative, otherwise z0dym,z0dyh become nan
        AMREX_ALWAYS_ASSERT(d_zref - disp_hgt_arr(i, j, 0) > 0.0);
        amrex::Real zTh = zt0 / (d_zref - disp_hgt_arr(i, j, 0));
        
        amrex::Real zodym = log(1.0 / z0h);
        amrex::Real zodyh = log(1.0 / zTh);

        // first guess for xsi 
        if (r > 0.0)
        {
            xsi = r * zodym / (1.0 - 5.0 * std::min(0.19,r));
        }
        else
        {
            xsi = r*zodym;
        }

        int niter = 0;
        amrex::Real error = 1000.0;

        while (error > errormax && niter < nitermax)
        {
            xsi1 = xsi;
            niter++;
            xsim0 = z0h * xsi;
            xsih0 = zTh * xsi;

            if (xsi < -0.01)
            {
                if (xsi >= xsim)
                {
                    fm = zodym - psim1(xx(xsi), xx(xsim0));
                }
                else
                {
                    fm = psim2(xsi, xsim0, xx(xsim0));
                }

                if (xsi >= xsih)
                {
                    fh = zodyh - psih1(xx(xsi), xx(xsih0));
                }
                else
                {
                    fh = psih2(xsi, xsih0, xx(xsih0));
                }
            }
            else if (xsi > 0.01)
            {
                if (xsi <= 1.0)
                {
                    fm = zodym - psim3(xsi, xsim0);
                    fh = zodyh - psih3(xsi, xsih0);
                }
                else
                {
                    fm = psim4(xsi, xsim0);
                    fh = psih4(xsi, xsih0);
                }
            }
            else
            {
                fm = zodym;
                fh = zodyh;
            }

            xsi = r * fm * fm / fh;
            error = std::abs(xsi - xsi1);
        }

        // limit fm and fh to avoid too large fluxes especially over large surface roughness.
        // Basically, make the maximum slowdown of the velocity not bigger than 50% in one timestep
        amrex::Real fm0 = sqrt((kk*kk)*vel*dt / 0.5 / d_zref);
        fm = std::max(fm0, fm);
        fh = std::max(fh, fh/fm*fm0);

        // drag coefficient C_D = k**2/fm**2
        // heat transfer coefficient C_H = k**2/fm/fh
        mom_trans_coef = (kk*kk) / (fm*fm);
        heat_trans_coef = (kk*kk) / fm / fh;
        ustar_arr(i, j, 0) = sqrt(mom_trans_coef) * vel;

        // set ustar > 0.2 to avoid too calm conditions at night for the turbulent transfer
        //ustar_arr(i, j, 0) = std::max(0.2, ustar_arr(i, j, 0));
        ustar_arr(i, j, 0) = sqrt(ustar_arr(i, j, 0) * ustar_arr(i, j, 0)+0.05*0.05); // following gSAM-SLM 
        tstar_arr(i, j, 0) = -kk * (thp - tsp) / fh;
        
        // recompute Ch and Cd for onsistency
        mom_trans_coef = std::pow((ustar_arr(i, j, 0)/vel),2); 
        heat_trans_coef = (kk/fh) * (ustar_arr(i, j, 0)/vel);

        // aerodynamic resistance between surface and reference level
        r_a_arr(i, j, 0) = fh / kk / ustar_arr(i, j, 0);


        amrex::Real vel_m = vel;
        // amrex::Real RiB = r; // TODO: not used?
        amrex::Real taux_sfc = -1.0 * mom_trans_coef * vel_m * ur_arr(i, j, 0) * (100.0 * pref_arr(i, j, 0) / 287.0 / tref_arr(i, j, 0));
        amrex::Real tauy_sfc = -1.0 * mom_trans_coef * vel_m * vr_arr(i, j, 0) * (100.0 * pref_arr(i, j, 0) / 287.0 / tref_arr(i, j, 0));

        // Output variables
        flbu_arr(i, j, d_khi_lsm) = taux_sfc;
        flbv_arr(i, j, d_khi_lsm) = tauy_sfc;
        ustar_arr(i, j, d_khi_lsm) = ustar_arr(i, j, 0);
        tstar_arr(i, j, d_khi_lsm) = tstar_arr(i, j, 0);
    });
}

void SLM::resistances(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const int d_nz_lsm = m_nz_lsm;
    const Real d_z0_soil = z0_soil;
    const Real d_T_opt = T_opt;
    const Real d_Rc_max = Rc_max;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);
    auto landtype_arr = landtype.const_array(mfi);
    auto vegetype_arr = vegetype.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->const_array(mfi);
    auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);

    auto LAI_arr = LAI.const_array(mfi);
    auto SAI_arr = SAI.const_array(mfi);
    auto t_cas_arr = t_cas.const_array(mfi);
    auto q_cas_arr = q_cas.const_array(mfi);
    auto ustar_arr = lsm_fab_vars[LsmVar_SLM::ustar]->const_array(mfi);
    auto tstar_arr = lsm_fab_vars[LsmVar_SLM::tstar]->const_array(mfi);

    auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->const_array(mfi);
    auto theta_FC_arr = lsm_fab_vars[LsmVar_SLM::theta_FC]->const_array(mfi);
    auto theta_WP_arr = lsm_fab_vars[LsmVar_SLM::theta_WP]->const_array(mfi);
    auto soil_transp_frac_arr = lsm_fab_vars[LsmVar_SLM::soil_transp_frac]->array(mfi);

    auto ztop_arr = ztop.const_array(mfi);
    auto disp_hgt_arr = disp_hgt.const_array(mfi);
    auto z0_sfc_arr = z0_sfc.const_array(mfi);
    auto Rgl_arr = Rgl.const_array(mfi);
    auto Rc_min_arr = Rc_min.const_array(mfi);
    auto hs_rc_arr = hs_rc.const_array(mfi);
    auto nroot_arr = nroot.const_array(mfi);
    auto dleaf_arr = dleaf.const_array(mfi);
    auto zrefxy_arr = zrefxy.const_array(mfi);
    auto net_rad_arr = net_rad.const_array(mfi);

    auto ur_arr = lsm_fab_vars[LsmVar_SLM::uref]->const_array(mfi);
    auto vr_arr = lsm_fab_vars[LsmVar_SLM::vref]->const_array(mfi);

    auto r_a_arr = r_a.const_array(mfi);
    auto r_b_arr = r_b.array(mfi);
    auto r_c_arr = r_c.array(mfi);
    auto r_d_arr = r_d.array(mfi);

    auto phi_1_arr = phi_1.array(mfi);
    auto phi_2_arr = phi_2.array(mfi);

    // Get CWPVT parameter (canopy wind extinction factor) if available
    const amrex::Real* d_cwpvt = (cwpvt.size() > 0) ? cwpvt.data() : nullptr;
    const amrex::Real d_cwpvt_default = 1.0;  // default value if not from table

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        for (int k = 0; k < d_nz_lsm; k++) {
            soil_transp_frac_arr(i, j, d_khi_lsm - k) = 0.0;
        }

        if (vegetype_arr(i, j, 0) == 0)
        {
            // for baresoil, r_d = r_a
            // under canopy resistance
            r_d_arr(i, j, 0) = r_a_arr(i, j, 0);
        }
        else
        {
            // for vegetated land surfaces
            amrex::Real Cs_dense;
            amrex::Real Cs_bare;
            amrex::Real Cs;
            amrex::Real rc_fac_rad, rc_fac_vpd, rc_fac_t, rc_fac_sw, d_root, tmp_radf;
            amrex::Real k_beer, f_shade, lai_sun, lai_shade, sw_sun, sw_shade, tmp_sun, tmp_shade, rc_fac_sun, rc_fac_shade;

            // ===================================================
            // Aerodynamic resistance for heat and vapor transfer under canopy space : r_d
            // ===================================================
            // NoahMP method: ResistanceLeafToGroundMod.F90
            // Reference: Niu et al. (2011), He et al. (2023)

            const amrex::Real MPE = 1.0e-6;
            const amrex::Real CONST_VON_KARMAN = 0.4;
            amrex::Real temp_diff = t_cas_arr(i, j, 0) - soilt_arr(i, j, d_khi_lsm);

            // Get CWPVT (canopy wind extinction parameter) for this vegetation type
            int veg_idx = landtype_arr(i, j, 0) - 1;  // Convert the vegetation class to a 0-based table index
            amrex::Real CanopyWindExtFac = (d_cwpvt != nullptr && veg_idx >= 0) ? d_cwpvt[veg_idx] : d_cwpvt_default;

            // Stability correction for undercanopy (MoStabCorrShUndCan)
            amrex::Real MoStabCorrShUndCan = 1.0;  // Initialize to neutral
            amrex::Real HeatSenGrdTmp = 0.0;  // Will be computed in flux calculation
            // For now, use simple stability correction based on temperature difference
            if (temp_diff < 0.0) {
                // Unstable: (1 - 15*z/L)^(-0.25)
                amrex::Real zeta = std::min(0.0, -0.1);  // Assume moderately unstable
                MoStabCorrShUndCan = std::pow(1.0 - 15.0 * zeta, -0.25);
            } else {
                // Stable: 1 + 4.7*z/L
                amrex::Real zeta = std::min(1.0, 0.1);  // Assume moderately stable
                MoStabCorrShUndCan = 1.0 + 4.7 * zeta;
            }

            // Wind extinction coefficient
            amrex::Real VegAreaIndEff = std::min(6.0, LAI_arr(i, j, 0) + SAI_arr(i, j, 0));
            amrex::Real CanopyHeight = ztop_arr(i, j, 0);
            amrex::Real WindExtCoeffCanopy = std::sqrt(CanopyWindExtFac * VegAreaIndEff * CanopyHeight * MoStabCorrShUndCan);

            // Roughness lengths
            amrex::Real RoughLenShVegGrd = d_z0_soil;  // Ground roughness for heat under canopy
            amrex::Real RoughLenShCanopy = z0_sfc_arr(i, j, 0);  // Canopy roughness for heat
            amrex::Real ZeroPlaneDispSfc = disp_hgt_arr(i, j, 0);  // Zero plane displacement

            // Exponential wind profile terms
            amrex::Real TMP1 = std::exp(-WindExtCoeffCanopy * RoughLenShVegGrd / CanopyHeight);
            amrex::Real TMP2 = std::exp(-WindExtCoeffCanopy * (RoughLenShCanopy + ZeroPlaneDispSfc) / CanopyHeight);
            amrex::Real TMPRAH2 = CanopyHeight * std::exp(WindExtCoeffCanopy) / WindExtCoeffCanopy * (TMP1 - TMP2);

            // Turbulent transfer coefficient KH
            amrex::Real KH = std::max(CONST_VON_KARMAN * ustar_arr(i, j, 0) * (CanopyHeight - ZeroPlaneDispSfc), MPE);

            // Undercanopy aerodynamic resistance
            r_d_arr(i, j, 0) = TMPRAH2 / KH;

            // Original SLM method (commented out for comparison)
            /*
            // temp_diff > 0 : stable undercanopy
            // temp_diff < 0 : unstable undercanopy
            amrex::Real temp_diff = t_cas_arr(i, j, 0) - soilt_arr(i, j, d_khi_lsm);

            // turbulent transfer coefficient under dense canopy
            if (temp_diff < 0.0)
            {
                Cs_dense = 0.004;
            }
            else
            {
                // rd_correc_fac : undercanopy stability parameter, in effect only for stable undercanopy
                amrex::Real rd_correc_fac = CONST_GRAV * ztop_arr(i, j, 0) * std::max(0.0, temp_diff) / soilt_arr(i, j, d_khi_lsm) / (std::pow(ustar_arr(i, j, 0), 2));
                // stable undercanopy - Cs_dense becomes smaller than 0.004
                Cs_dense = 0.004 / (1.0 + 0.5 * std::min(10.0, rd_correc_fac));
            }

            // turbulent transfer coefficient over the exposed topsoil
            // typical value of Cs_bare ~0.2
            Cs_bare = 0.4 / 0.13 * std::pow((d_z0_soil * ustar_arr(i, j, 0) / (1.5e-5)), -0.45);

            // turbulence transfer coefficient undercanopy
            // LAI-weighed sum of Cs_bare and Cs_dense
            Cs = Cs_bare * std::exp(-1.0 * LAI_arr(i, j, 0)) + Cs_dense * (1.0 - std::exp(-1.0 * LAI_arr(i, j, 0)));

            // Undercanopy aerodynamic resistance depends on the weighed sum of the dense canopy covered soil
            // and baresoil turbulent transfer coefficient and friction velocity
            //   Reference: [Oleson et al., 2004] [Zeng et al., 2005]
            r_d_arr(i, j, 0) = std::min(400., 1.0 / ustar_arr(i, j, 0) / Cs); // prevent r_d from getting too large under stable condition
            */

            // ===================================================
            // Leaf boundary layer resistance : r_b
            // ===================================================
            // Original SLM formula:
            // r_b_arr(i, j, 0) = 0.5 * r_a_arr(i, j, 0);

            // NoahMP leaf boundary resistance from ResistanceLeafToGroundMod.F90.
            const amrex::Real wind_ref = std::max(
                std::sqrt(ur_arr(i, j, 0) * ur_arr(i, j, 0)
                          + vr_arr(i, j, 0) * vr_arr(i, j, 0)), MPE);
            const amrex::Real wind_canopy_top = std::max(
                wind_ref
                    * std::log((CanopyHeight - ZeroPlaneDispSfc + RoughLenShCanopy)
                               / RoughLenShCanopy)
                    / std::log(zrefxy_arr(i, j, 0) / RoughLenShCanopy),
                MPE);
            const amrex::Real rb_factor = WindExtCoeffCanopy * 50.0
                / std::max(1.0 - std::exp(-WindExtCoeffCanopy / 2.0), MPE);
            r_b_arr(i, j, 0) = std::min(50.0, std::max(5.0,
                rb_factor * std::sqrt(dleaf_arr(i, j, 0) / wind_canopy_top)));

            // Original SLM Method 1 (commented out for reference)
            // turbulent transfer coefficient between canopy surface and canopy air : Cv = 0.01m/s^-0.5
            // characteristic dimension of the leaves in the direction of wind flux : d_leaf = 0.04m
            // r_b_arr(i, j, 0) = 1.0 / 0.01 * std::pow( ustar_arr(i, j, 0) / 0.04, -0.5) / std::max(0.1, LAI_arr(i, j, 0));
            // Above equation seems to overestimate LHF 

            // ===================================================
            // Stomatal resistance : r_c
            // ===================================================
            // radiation factor (Noah/NoahMP Jarvis bulk formulation)
            // Reference: Noah LSM CANRES subroutine, Chen et al. (1996, 2001)
            amrex::Real RadFac = 0.55 * net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1) * 2.0 / Rgl_arr(i, j, 0) / LAI_arr(i, j, 0);
            rc_fac_rad = (Rc_min_arr(i, j, 0) / d_Rc_max + RadFac) / (1.0 + RadFac);
            rc_fac_rad = std::max(rc_fac_rad, 0.0001);

            /* Original SLM sunlit/shaded approach (commented out to match Noah bulk approach)
            // TODO: check if this is the correct downwelling SW to use
            //amrex::Real tmp_radf = 0.55 * net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1) * 2.0 / Rgl_arr(i, j, 0) / LAI_arr(i, j, 0);
            //rc_fac_rad = (Rc_min_arr(i, j, 0) / d_Rc_max + tmp_radf) / (1.0 + tmp_radf);
            // modify above to account for the shaded and sunlit part of leaves
            k_beer = phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0); // extinction coef
            f_shade = 0.1 ; // fraction or radiation reaching shaded leaves (empirical)

            // partition LAI into sunlit and shaded components
            lai_sun = std::max(1.e-6, (1. - std::exp(-1.0 * k_beer * LAI_arr(i, j, 0))) / k_beer);
            lai_shade = std::max(0., LAI_arr(i, j, 0) - lai_sun);

            // Radiation reaching to sunlit and shaded leaves
            sw_sun = net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1);
            sw_shade = f_shade * sw_sun;

            // compute radiation factors (NoahMP Jarvis form - no LAI division)
            tmp_sun = 0.55 * sw_sun * 2. / Rgl_arr(i, j, 0);
            tmp_shade = 0.55 * sw_shade * 2. / Rgl_arr(i, j, 0);
            rc_fac_sun   = (Rc_min_arr(i,j,0) / d_Rc_max + tmp_sun) / (1.0 + tmp_sun);
            rc_fac_shade = (Rc_min_arr(i,j,0) / d_Rc_max + tmp_shade) / (1.0 + tmp_shade);

            // combine weighted by LAI
            rc_fac_rad = (lai_sun * rc_fac_sun + lai_shade * rc_fac_shade) / LAI_arr(i, j, 0);
            *///


            // vapor pressure deficit factor (Noah/NoahMP Jarvis formulation)
            // Noah uses VPD in mixing ratio units (kg/kg), not vapor pressure (hPa)
            // Reference: Noah LSM CANRES subroutine, Chen et al. (1996, 2001)
            amrex::Real qsatw;
            erf_qsatw(t_cas_arr(i, j, 0), pref_arr(i, j, 0), qsatw);
            amrex::Real VPD_mixratio = qsatw - q_cas_arr(i, j, 0); // VPD in mixing ratio [kg/kg]
            rc_fac_vpd = 1.0 / (1.0 + hs_rc_arr(i, j, 0) * VPD_mixratio);
            rc_fac_vpd = std::max(rc_fac_vpd, 0.01);

            /* Original SLM VPD calculation using vapor pressure in hPa (commented out)
            // Above is modified following changes in gSAM-SLM
            amrex::Real e_cas, es_cas;
            e_cas = q_cas_arr(i, j, 0) * pref_arr(i, j ,0)/(0.622+0.388*q_cas_arr(i, j, 0)); //vapor pressure in hPa (mb)
            es_cas = erf_esatw(t_cas_arr(i, j, 0));
            // NoahMP Jarvis form: hyperbolic response to VPD
            rc_fac_vpd = 1.0 / (1.0 + hs_rc_arr(i, j, 0) * (es_cas - e_cas));
            *///

            // temperature factor
            rc_fac_t = std::max(0., 1.0 - 0.0016 * std::pow(d_T_opt - t_cas_arr(i, j, 0), 2));

            // NoahMP BTR_OPTION=1 soil-moisture stress factor. All rooted
            // layers remain in the depth denominator, including dry layers.
            const int nroot_cell = nroot_arr(i, j, 0);
            d_root = 0.0;
            for (int k = 0; k < nroot_cell; k++) {
                const int lsm_k = d_khi_lsm - k;
                d_root += s_depth_arr(i, j, lsm_k);
            }

            rc_fac_sw = 0.0;
            for (int k = 0; k < nroot_cell; k++) {
                const int lsm_k = d_khi_lsm - k;
                const amrex::Real theta_liq = soilw_arr(i, j, lsm_k) * poro_soil_arr(i, j, lsm_k);
                const amrex::Real moisture_range = theta_FC_arr(i, j, lsm_k) - theta_WP_arr(i, j, lsm_k);
                const amrex::Real soil_wet_fac = std::min(1.0, std::max(0.0,
                    (theta_liq - theta_WP_arr(i, j, lsm_k)) / moisture_range));
                soil_transp_frac_arr(i, j, lsm_k) = std::max(1.0e-6,
                    s_depth_arr(i, j, lsm_k) / std::max(1.0e-6, d_root) * soil_wet_fac);
                rc_fac_sw += soil_transp_frac_arr(i, j, lsm_k);
            }

            rc_fac_sw = std::max(1.0e-6, rc_fac_sw);
            for (int k = 0; k < nroot_cell; k++) {
                const int lsm_k = d_khi_lsm - k;
                soil_transp_frac_arr(i, j, lsm_k) /= rc_fac_sw;
            }

            tmp_radf = std::max(1.0e-6, rc_fac_rad*rc_fac_vpd*rc_fac_t*rc_fac_sw);
            // Leaf-level stomatal resistance. LAI is applied once when this
            // resistance is converted to canopy transpiration conductance.
            r_c_arr(i, j, 0) = std::min(d_Rc_max, Rc_min_arr(i, j, 0) / tmp_radf);

            // ===================================================
            // r_litter : litter resistance - not used in this version
            // ===================================================
            //r_litter_arr(i, j, 0) = 0.0;
            //   Ref. [Sakaguchi and Zeng 2009]
            //   Set litter LAI as 0.5
            // r_litter_arr(i, j, 0) = 1.0 / 0.004 / ustar_arr(i, j, 0) * (1.0 - std::exp(-0.5));
        }
    });
}

void SLM::fluxes_canopy(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const int d_nz_lsm = m_nz_lsm;
    const Real dt = m_dt;
    const int niter = std::max(1, static_cast<int>(std::round(dt)));
    const Real dt_iter = dt/niter;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);

    auto LAI_arr = LAI.const_array(mfi);
    
    auto q_cas_arr = q_cas.const_array(mfi);
    auto t_sfc_arr = t_sfc.const_array(mfi);
    auto q_sfc_arr = q_sfc.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);

    auto lhf_canop_arr = lhf_canop.array(mfi);
    auto evp_canop_arr = evp_canop.array(mfi);
    auto lhf_soil_arr = lhf_soil.array(mfi);
    auto lhf_air_arr = lhf_air.array(mfi);
    
    auto shf_canop_arr = shf_canop.array(mfi);

    auto vegetype_arr = vegetype.const_array(mfi);
    auto vege_YES_arr = vege_YES.const_array(mfi);

    auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

    auto mw_arr = mw.array(mfi);
    auto mw_mx_arr = mw_mx.const_array(mfi);
    auto mw_inc_arr = mw_inc.array(mfi);
    auto mws_arr = mws.const_array(mfi);
    auto cp_vege_arr = cp_vege.const_array(mfi);

    auto dref_arr = lsm_fab_vars[LsmVar_SLM::dref]->const_array(mfi);
    auto qr_arr = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
    auto tr_arr = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
    auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
    auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);
    auto t_canop_arr = t_canop.array(mfi);

    auto r_a_arr = r_a.const_array(mfi);
    auto r_b_arr = r_b.const_array(mfi);
    auto r_c_arr = r_c.const_array(mfi);
    auto r_d_arr = r_d.const_array(mfi);
    auto r_soil_arr = r_soil.array(mfi);

    auto wet_canop_arr = wet_canop.array(mfi);
    auto evapo_dry_arr = evapo_dry.array(mfi);

    auto slm_diag_arr = slm_diag.array(mfi);
    
    auto net_rad_arr = net_rad.const_array(mfi);
        
    auto prsfc_arr  = lsm_fab_vars[LsmVar_SLM::prsfc]->array(mfi);

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }
        amrex::Real evapo_wet;

        if (vegetype_arr(i, j, 0) == 1)
        {
            amrex::Real shf0 = 0.;
            amrex::Real lhf0 = 0.;
            amrex::Real evp0 = 0.;
            amrex::Real evapo_wet0 = 0.;
            amrex::Real drain0 = 0.;
            // SAM rhow[nz] = air density at vertical velocity levels, kg/m^3
            const amrex::Real rhow = dref_arr(i, j, 0); // TODO: double check this
            amrex::Real qsat_canop;
            amrex::Real cp_vege_tot, t_canop_inc;

            for (int iter = 0; iter < niter; iter++)
            {   
                shf_canop_arr(i, j, 0) = (t_canop_arr(i, j, 0) - t_sfc_arr(i, j, 0)) * rhow * cp * 2.0 * LAI_arr(i, j, 0) / r_b_arr(i, j, 0);
                shf0 += shf_canop_arr(i, j, 0);
                    
                // Evaporation from canopy
                // only treat the case where t_canop_arr > tfriz; 
                AMREX_ALWAYS_ASSERT(t_canop_arr(i, j, 0) > tfriz);
        
                erf_qsatw(t_canop_arr(i, j, 0), pref_arr(i, j, 0), qsat_canop);
        
                // direct evaporation from the water held on canopy
                // evaporation/dew only possible if canopy temperature is above freezing
                evapo_wet = std::min(mw_arr(i, j, 0)/dt_iter, ((qsat_canop - q_sfc_arr(i, j, 0)) * rhow * LAI_arr(i, j, 0) / (r_b_arr(i, j, 0))*vege_YES_arr(i, j, 0)));
        
                // increment/decrement of the water amount held on leaves following the direct evaporation/dew formation
                mw_inc_arr(i, j, 0) = -dt_iter*evapo_wet; // evapo_wet [kg/m2s=mm/s]
                mw_arr(i, j, 0) += mw_inc_arr(i, j, 0);
                wet_canop_arr(i, j, 0) = std::min(1.0, mw_arr(i, j, 0)/mw_mx_arr(i, j, 0));
        
                // Transpiration - only ocurs when qsat_canop > qsfc
                evapo_dry_arr(i, j, 0) = std::max(0.,(qsat_canop - q_sfc_arr(i, j, 0))*rhow*(1.0 - wet_canop_arr(i, j, 0))*LAI_arr(i, j, 0)/(r_b_arr(i, j, 0) + r_c_arr(i, j, 0))*vege_YES_arr(i, j, 0));

                // Convert evaporation (kg/m2/s) to latent heat flux (W/m2)
                lhf_canop_arr(i, j, 0) = lcond*(evapo_wet+evapo_dry_arr(i, j, 0));
                evp_canop_arr(i, j, 0) = evapo_wet + evapo_dry_arr(i, j, 0); 
                lhf0 += lhf_canop_arr(i, j, 0);
                evp0 += evp_canop_arr(i, j, 0);
                evapo_wet0 += evapo_wet;

                // Update vegetation moisture storage
                if (mw_arr(i, j, 0) > mw_mx_arr(i, j, 0)) 
                {
                    drain0 += (mw_arr(i, j, 0) - mw_mx_arr(i, j, 0))/dt_iter; // dripping excess of dew
                    mw_arr(i, j, 0) = mw_mx_arr(i, j, 0); 
                }

                // Update vegetation temeprature
                cp_vege_tot = cp_vege_arr(i, j, 0) + mw_arr(i, j, 0) * 1.e-3 * cp_water;
                t_canop_inc = dt_iter / std::max(1.0e-3, cp_vege_tot)*(net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) - shf_canop_arr(i, j, 0) - lhf_canop_arr(i, j, 0)) * vege_YES_arr(i, j, 0);
                t_canop_arr(i, j, 0) = std::min(343.0, t_canop_arr(i, j, 0) + t_canop_inc);  // 343 K = t_canop_max
            }
            shf_canop_arr(i, j, 0) = shf0 / static_cast<amrex::Real>(niter);
            lhf_canop_arr(i, j, 0) = lhf0 / static_cast<amrex::Real>(niter);
            evp_canop_arr(i, j, 0) = evp0 / static_cast<amrex::Real>(niter);
            prsfc_arr(i, j, 0) += (drain0 / static_cast<amrex::Real>(niter));
            slm_diag_arr(i, j, 0, SLM_Diag::precip_sfc) = prsfc_arr(i, j, 0);
            slm_diag_arr(i, j, 0, SLM_Diag::evapo_wet) = evapo_wet0 / static_cast<amrex::Real>(niter);
            slm_diag_arr(i, j, 0, SLM_Diag::drain) += (drain0/static_cast<amrex::Real>(niter)); 
        }
        else
        {
            shf_canop_arr(i, j, 0) = 0.;
            lhf_canop_arr(i, j, 0) = 0.;
            evp_canop_arr(i, j, 0) = 0.;
            wet_canop_arr(i, j, 0) = 0.;
            mw_arr(i, j, 0) = 0.;
            evapo_wet = 0.;
            evapo_dry_arr(i, j, 0) = 0.;
            t_canop_arr(i, j, 0) = tr_arr(i, j, 0);
        }   
    });
}

void SLM::solve_ground_skin_temperature(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);
    auto landtype_arr = landtype.const_array(mfi);
    
    auto IMPERV_arr = IMPERV.const_array(mfi);

    auto q_cas_arr = q_cas.const_array(mfi);
    auto t_cas_arr = t_cas.const_array(mfi);
    auto t_ground_skin_arr = t_ground_skin.array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);
    auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->const_array(mfi);
    auto sst_cond_arr = lsm_fab_vars[LsmVar_SLM::sst_cond]->const_array(mfi);
    auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->const_array(mfi);

    auto lhf_canop_arr = lhf_canop.array(mfi);
    auto lhf_soil_arr = lhf_soil.array(mfi);
    auto lhf_air_arr = lhf_air.array(mfi);
    
    auto evp_canop_arr = evp_canop.array(mfi);
    auto evp_soil_arr = evp_soil.array(mfi);
    auto evp_air_arr = evp_air.array(mfi);

    auto shf_soil_arr = shf_soil.array(mfi);
    auto shf_canop_arr = shf_canop.array(mfi);
    auto shf_air_arr = shf_air.array(mfi);
    
    auto t_sfc_arr = t_sfc.array(mfi);
    auto q_sfc_arr = q_sfc.array(mfi);
    
    auto q_gr_arr = q_gr.array(mfi);
    auto sdew_arr = sdew.array(mfi);

    auto vegetype_arr = vegetype.const_array(mfi);
    auto vege_YES_arr = vege_YES.const_array(mfi);

    auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

    auto mw_arr = mw.const_array(mfi);
    auto mw_mx_arr = mw_mx.const_array(mfi);
    auto mw_inc_arr = mw_inc.array(mfi);
    auto mws_arr = mws.const_array(mfi);

    auto dref_arr = lsm_fab_vars[LsmVar_SLM::dref]->const_array(mfi);
    auto qr_arr = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
    auto tr_arr = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
    auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
    auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);
    auto w_s_WP_arr = lsm_fab_vars[LsmVar_SLM::w_s_WP]->const_array(mfi);
    auto t_canop_arr = t_canop.const_array(mfi);

    auto r_a_arr = r_a.const_array(mfi);
    auto r_b_arr = r_b.const_array(mfi);
    auto r_c_arr = r_c.const_array(mfi);
    auto r_d_arr = r_d.const_array(mfi);
    auto r_soil_arr = r_soil.array(mfi);

    auto wet_canop_arr = wet_canop.array(mfi);
    auto evapo_dry_arr = evapo_dry.array(mfi);

    auto slm_diag_arr = slm_diag.array(mfi);
    auto net_rad_arr = net_rad.array(mfi);
    auto t_skin_arr = t_skin.array(mfi);
    auto lwref_arr = lsm_fab_vars[LsmVar_SLM::lwref]->const_array(mfi);
    auto LAI_arr = LAI.const_array(mfi);
    auto SAI_arr = SAI.const_array(mfi);

    const amrex::Real d_eg_soil = eg_soil_rad;
    constexpr int max_iterations = 5;
    constexpr amrex::Real derivative_step = 0.1;
    constexpr amrex::Real correction_limit = 5.0;
    constexpr amrex::Real convergence_tolerance = 0.01;
    const amrex::Real d_rsurf_exp = rsurf_exp;


    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        // SAM rhow[nz] = air density at vertical velocity levels, kg/m^3
        const amrex::Real rhow = dref_arr(i, j, 0); // TODO: double check this
        const bool vegetated = vegetype_arr(i, j, 0) == 1;
        const amrex::Real q_air = vegetated ? q_cas_arr(i, j, 0) : qr_arr(i, j, 0);

        /* Legacy SLM soil diffusion and resistance formulation:
        amrex::Real q_ground_guess;
        if (t_ground_skin_arr(i, j, 0) >= tfriz) {
            erf_qsatw(t_ground_skin_arr(i, j, 0), pref_arr(i, j, 0), q_ground_guess);
            if (mws_arr(i, j, 0) == 0.0) {
                q_ground_guess *= fh_calc(t_ground_skin_arr(i, j, 0),
                                          m_pot_sat_arr(i, j, d_khi_lsm),
                                          soilw_arr(i, j, d_khi_lsm),
                                          Bconst_arr(i, j, d_khi_lsm));
            }
        } else {
            erf_qsati(t_ground_skin_arr(i, j, 0), pref_arr(i, j, 0), q_ground_guess);
        }
        amrex::Real soil_diff;
        if (soilw_arr(i, j, d_khi_lsm) >= w_s_FC_arr(i,j,d_khi_lsm) ||
            q_air > q_ground_guess) {
            soil_diff = 1.0;
        } else {
            soil_diff = 0.25 * std::pow(1.0 - std::cos(M_PI *
                std::max(0.01, soilw_arr(i, j, d_khi_lsm)) /
                w_s_FC_arr(i, j, d_khi_lsm)), 2);
        }
        amrex::Real totalR_soil;
        if (!vegetated) {
            r_soil_arr(i, j, 0) = std::min(10000.0, std::max(100.0,
                r_a_arr(i, j, 0) * (1.0 / soil_diff - 1.0)));
            totalR_soil = r_soil_arr(i, j, 0) + r_a_arr(i, j, 0);
        } else {
            r_soil_arr(i, j, 0) = std::min(10000.0, std::max(50.0,
                r_d_arr(i, j, 0) * (1.0 / soil_diff - 1.0)));
            totalR_soil = r_soil_arr(i, j, 0) + r_d_arr(i, j, 0) + r_litter;
        }
        */

        // NoahMP option 1 (Sakaguchi and Zeng, 2009) ground resistance.
        const amrex::Real soil_wetness = std::min(1.0,
            std::max(0.0, soilw_arr(i, j, d_khi_lsm)));
        const amrex::Real dry_soil_thickness = s_depth_arr(i, j, d_khi_lsm)
            * (std::exp(std::pow(1.0 - soil_wetness, d_rsurf_exp)) - 1.0)
            / (std::exp(1.0) - 1.0);
        const amrex::Real soil_b = Bconst_arr(i, j, d_khi_lsm);
        const amrex::Real vapor_diffusivity = soil_b > 0.0
            ? 2.2e-5 * poro_soil_arr(i, j, d_khi_lsm)
              * poro_soil_arr(i, j, d_khi_lsm)
              * std::pow(1.0 - w_s_WP_arr(i, j, d_khi_lsm),
                         2.0 + 3.0 / soil_b)
            : 0.0;
        amrex::Real soil_resistance = dry_soil_thickness
            / std::max(vapor_diffusivity, 1.0e-12);
        const amrex::Real theta_liq = soilw_arr(i, j, d_khi_lsm)
            * poro_soil_arr(i, j, d_khi_lsm);
        if (theta_liq < 0.01 || landtype_arr(i, j, 0) == 13 || soil_b <= 0.0) {
            soil_resistance = 1.0e6;
        }
        r_soil_arr(i, j, 0) = soil_resistance;

        const amrex::Real aerodynamic_resistance = vegetated
            ? r_d_arr(i, j, 0) : r_a_arr(i, j, 0);
        const amrex::Real totalR_soil = aerodynamic_resistance + soil_resistance
            + (vegetated ? r_litter : 0.0);

        amrex::Real k_dry, k_sat, conductivity;
        if (landtype_arr(i, j, 0) == 15) {
            conductivity = 1.6;
        } else {
            const amrex::Real dry_density = (1.0 - poro_soil_arr(i, j, d_khi_lsm)) * 2700.0;
            k_dry = (0.135 * dry_density + 64.7) / (2700.0 - 0.947 * dry_density);
            k_sat = std::pow(sst_cond_arr(i, j, d_khi_lsm), 1.0 - poro_soil_arr(i, j, d_khi_lsm));
            k_sat *= std::pow(soilt_arr(i, j, d_khi_lsm) > tfriz ? 0.57 : 1.60,
                              poro_soil_arr(i, j, d_khi_lsm));
            const amrex::Real ke = std::log10(std::max(0.1, soilw_arr(i, j, d_khi_lsm))) + 1.0;
            conductivity = ke * (k_sat - k_dry) + k_dry;
        }
        const amrex::Real conduction_coeff = 2.0 * conductivity / s_depth_arr(i, j, d_khi_lsm);
        const amrex::Real emv = 1.0 - std::exp(-(LAI_arr(i, j, 0) + SAI_arr(i, j, 0)));
        const amrex::Real emg = d_eg_soil;
        const amrex::Real net_sw = net_rad_arr(i, j, 0, SLM_NetRad::net_sw2);
        const amrex::Real lwdn = lwref_arr(i, j, 0);
        const amrex::Real tv = t_canop_arr(i, j, 0);
        const amrex::Real tsoil = soilt_arr(i, j, d_khi_lsm);
        const amrex::Real potential_factor = std::pow(1000.0 / pref_arr(i, j, 0), rair / cp);

        auto residual = [=] AMREX_GPU_DEVICE (amrex::Real tg) noexcept {
            amrex::Real q_ground;
            amrex::Real dew_factor = 1.0;
            if (tg >= tfriz) {
                erf_qsatw(tg, pref_arr(i, j, 0), q_ground);
                if (mws_arr(i, j, 0) == 0.0) {
                    const amrex::Real humidity_factor = fh_calc(tg,
                        m_pot_sat_arr(i, j, d_khi_lsm), soilw_arr(i, j, d_khi_lsm),
                        Bconst_arr(i, j, d_khi_lsm));
                    dew_factor = humidity_factor > 0.99 ? 1.0 : 0.0;
                    q_ground *= humidity_factor;
                }
            } else {
                erf_qsati(tg, pref_arr(i, j, 0), q_ground);
            }
            const amrex::Real sensible = vegetated
                ? (tg - t_cas_arr(i, j, 0)) * rhow * cp / r_d_arr(i, j, 0)
                : (tg - tr_arr(i, j, 0)) * potential_factor * rhow * cp / r_a_arr(i, j, 0);
            amrex::Real evaporation = (q_ground - q_air) * rhow / totalR_soil
                                    * (1.0 - IMPERV_arr(i, j, 0));
            if (evaporation < 0.0) evaporation *= dew_factor;
            const amrex::Real irg = emg * sigma * std::pow(tg, 4)
                                  - emg * (1.0 - emv) * lwdn
                                  - emg * emv * sigma * std::pow(tv, 4);
            return net_sw - irg - sensible - lcond * evaporation
                   - conduction_coeff * (tg - tsoil);
        };

        amrex::Real tg = t_ground_skin_arr(i, j, 0);
        for (int iter = 0; iter < max_iterations; ++iter) {
            const amrex::Real f = residual(tg);
            const amrex::Real fprime = (residual(tg + derivative_step)
                                      - residual(tg - derivative_step)) / (2.0 * derivative_step);
            AMREX_ALWAYS_ASSERT(std::isfinite(f) && std::isfinite(fprime));
            if (std::abs(fprime) < 1.0e-12) break;
            const amrex::Real delta = std::max(-correction_limit,
                std::min(correction_limit, -f / fprime));
            tg += delta;
            if (std::abs(delta) < convergence_tolerance) {
                break;
            }
        }

        amrex::Real q_ground;
        amrex::Real dew_factor = 1.0;
        if (tg >= tfriz) {
            erf_qsatw(tg, pref_arr(i, j, 0), q_ground);
            if (mws_arr(i, j, 0) == 0.0) {
                const amrex::Real humidity_factor = fh_calc(tg,
                    m_pot_sat_arr(i, j, d_khi_lsm), soilw_arr(i, j, d_khi_lsm),
                    Bconst_arr(i, j, d_khi_lsm));
                dew_factor = humidity_factor > 0.99 ? 1.0 : 0.0;
                q_ground *= humidity_factor;
            }
        } else {
            erf_qsati(tg, pref_arr(i, j, 0), q_ground);
        }
        const amrex::Real sensible = vegetated
            ? (tg - t_cas_arr(i, j, 0)) * rhow * cp / r_d_arr(i, j, 0)
            : (tg - tr_arr(i, j, 0)) * potential_factor * rhow * cp / r_a_arr(i, j, 0);
        amrex::Real evaporation = (q_ground - q_air) * rhow / totalR_soil
                                * (1.0 - IMPERV_arr(i, j, 0));
        if (evaporation < 0.0) evaporation *= dew_factor;
        const amrex::Real irg = emg * sigma * std::pow(tg, 4)
                              - emg * (1.0 - emv) * lwdn
                              - emg * emv * sigma * std::pow(tv, 4);
        const amrex::Real ground_conduction = conduction_coeff * (tg - tsoil);

        t_ground_skin_arr(i, j, 0) = tg;
        q_gr_arr(i, j, 0) = q_ground;
        sdew_arr(i, j, 0) = dew_factor;
        shf_soil_arr(i, j, 0) = sensible;
        evp_soil_arr(i, j, 0) = evaporation;
        lhf_soil_arr(i, j, 0) = lcond * evaporation;
        t_sfc_arr(i, j, 0) = vegetated ? t_cas_arr(i, j, 0) : tg;
        q_sfc_arr(i, j, 0) = vegetated ? q_cas_arr(i, j, 0) : q_ground;
        shf_air_arr(i, j, 0) = vegetated ? shf_canop_arr(i, j, 0) + sensible : sensible;
        evp_air_arr(i, j, 0) = vegetated ? evp_canop_arr(i, j, 0) + evaporation : evaporation;
        lhf_air_arr(i, j, 0) = vegetated ? lhf_canop_arr(i, j, 0) + lcond * evaporation
                                          : lcond * evaporation;

        const amrex::Real tir1 = emv * sigma * std::pow(tv, 4);
        const amrex::Real tir2 = emg * sigma * std::pow(tg, 4);
        const amrex::Real lwdn2 = (1.0 - emv) * lwdn + tir1;
        const amrex::Real lwup2 = tir2 + (1.0 - emg) * lwdn2;
        const amrex::Real irc = -emv * (1.0 + (1.0 - emv) * (1.0 - emg)) * lwdn
                              - emv * emg * sigma * std::pow(tg, 4)
                              + (2.0 - emv * (1.0 - emg)) * emv * sigma * std::pow(tv, 4);
        const amrex::Real lwup1 = lwdn + irc + irg;
        net_rad_arr(i, j, 0, SLM_NetRad::tir2) = tir2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup2) = lwup2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwdn2) = lwdn2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lw2) = -irg;
        net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) = net_sw - irg;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup1) = lwup1;
        const amrex::Real emiss_sfc = emv + emg * (1.0 - emv)
                                    + emv * (1.0 - emv) * (1.0 - emg);
        t_skin_arr(i, j, 0) = std::pow(std::max(1.0e-6,
            (lwup1 - (1.0 - emiss_sfc) * lwdn) / (emiss_sfc * sigma)), 0.25);

        slm_diag_arr(i, j, 0, SLM_Diag::grflux) = -ground_conduction;
    });
}

void SLM::soil_water(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const int d_nz_lsm = m_nz_lsm;
    const Real dt = m_dt;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);
    auto landtype_arr = landtype.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->array(mfi);

    auto sst_capa_arr = lsm_fab_vars[LsmVar_SLM::sst_capa]->const_array(mfi);
    auto sst_cond_arr = lsm_fab_vars[LsmVar_SLM::sst_cond]->const_array(mfi);
    auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->const_array(mfi);
    auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->const_array(mfi);

    auto shf_soil_arr = shf_soil.const_array(mfi);
    auto lhf_soil_arr = lhf_soil.const_array(mfi);
    auto net_rad_arr = net_rad.const_array(mfi);

    auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);
    auto LAI_arr = LAI.const_array(mfi);
    auto IMPERV_arr = IMPERV.const_array(mfi);
    auto precip_extinc_arr = precip_extinc.const_array(mfi);

    auto mw_arr = mw.array(mfi);
    auto mw_mx_arr = mw_mx.const_array(mfi);
    auto mw_inc_arr = mw_inc.array(mfi);
    auto mws_arr = mws.array(mfi);
    auto mws_mx_arr = mws_mx.const_array(mfi);
    
    auto evp_soil_arr = evp_soil.const_array(mfi);
    auto evapo_dry_arr = evapo_dry.const_array(mfi);
    auto prsfc_arr  = lsm_fab_vars[LsmVar_SLM::prsfc]->array(mfi);
    auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
    auto ks_arr = lsm_fab_vars[LsmVar_SLM::ks]->const_array(mfi);
    auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);
    auto w_s_WP_arr = lsm_fab_vars[LsmVar_SLM::w_s_WP]->const_array(mfi);
    auto w_s_FC_arr = lsm_fab_vars[LsmVar_SLM::w_s_FC]->const_array(mfi);
    auto soil_transp_frac_arr = lsm_fab_vars[LsmVar_SLM::soil_transp_frac]->const_array(mfi);

    auto slm_diag_arr = slm_diag.array(mfi);

    auto dsw_vars = soilw_vars.array(mfi);

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        amrex::Real aa, bb, cc, dd;
        amrex::Real precip_in = 0.0;
        amrex::Real precip_sfc = prsfc_arr(i, j, 0);

        if (landtype_arr(i, j, 0) == 15) {
            // ice
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;
                soilw_arr(i, j, lsm_k) = 0.0;
            }
        } else {
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;
                dsw_vars(i, j, lsm_k, SLM_DSW::sdepth_mm) = s_depth_arr(i, j, lsm_k) * 1.0e3;
            }
            // Evaporate standing water
            mws_arr(i, j, 0) = std::max(0., mws_arr(i, j, 0) - evp_soil_arr(i, j, 0) * dt); 
            
            // Calculate precipitation infiltration rate into the first soil layer
            bool any_less_than_one = false;
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;
                if(soilw_arr(i, j, lsm_k) < 1.0) {
                    any_less_than_one = true;
                    break;
                }   
            }
            if (any_less_than_one && soilt_arr(i, j, d_khi_lsm) >= tfriz)
            {
                precip_in = (1.-IMPERV_arr(i, j, 0))*std::min(precip_sfc + mws_arr(i, j, 0)/dt, ks_arr(i, j, d_khi_lsm));
            } else {
                precip_in = 0.0;
            }
            // calculate diffusion coefficient and velocity for soil moisture transfer
            // at each interfacial layer(= between adjacent soil layers)
            for (int k = 0; k < d_nz_lsm - 1; k++) {
                const int lsm_k = d_khi_lsm - k;

                if (soilt_arr(i, j, lsm_k) >= tfriz && soilt_arr(i, j, lsm_k-1) >= tfriz)
                {
                    const amrex::Real sd_lo = dsw_vars(i, j, lsm_k - 1, SLM_DSW::sdepth_mm);
                    const amrex::Real sd    = dsw_vars(i, j, lsm_k    , SLM_DSW::sdepth_mm);
                    const amrex::Real dz = sd + sd_lo;

                    dsw_vars(i, j, lsm_k, SLM_DSW::sh_eff_cond) =
                      (sd * std::pow(
                              soilw_arr(i, j, lsm_k),
                              Bconst_arr(i, j, lsm_k) + 2.0) +
                       sd_lo * std::pow(
                                 soilw_arr(i, j, lsm_k - 1),
                                 Bconst_arr(i, j, lsm_k) + 2.0)) /
                      dz * ks_arr(i, j, lsm_k) * Bconst_arr(i, j, lsm_k) *
                      std::abs(m_pot_sat_arr(i, j, lsm_k)) /
                      poro_soil_arr(i, j, lsm_k);

                    dsw_vars(i, j, lsm_k, SLM_DSW::sh_eff_vel) =
                      (sd * std::pow(
                              soilw_arr(i, j, lsm_k),
                              2.0 * Bconst_arr(i, j, lsm_k) + 2.0) +
                       sd_lo * std::pow(
                                 soilw_arr(i, j, lsm_k - 1),
                                 2.0 * Bconst_arr(i, j, lsm_k) + 2.0)) /
                      dz * ks_arr(i, j, lsm_k) / poro_soil_arr(i, j, lsm_k);
                } else {
                    // no water movement between two layers, one of which is frozen
                    dsw_vars(i, j, lsm_k, SLM_DSW::sh_eff_cond) = 0.0;
                    dsw_vars(i, j, lsm_k, SLM_DSW::sh_eff_vel) = 0.0;
                }
            }

            // make saturated soil for wetlands
            if (landtype_arr(i, j, 0) == 11)
            {
                for (int k = 0; k < d_nz_lsm; k++) {
                    const int lsm_k = d_khi_lsm - k;
                    soilw_arr(i, j, lsm_k) = w_s_FC_arr(i, j, lsm_k);
                }
            }

            // from FDE by the implicit method, Thomas algorithm is applied
            //  aa: terms related with soilw(k-1)
            //  bb: terms related with soilw(k)
            //  cc: terms related with soilw(k+1)
            //  dd: current soil wetness - sink + source
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;
                if (soilw_arr(i, j, lsm_k) < w_s_WP_arr(i, j, lsm_k))
                {
                    dsw_vars(i, j, lsm_k, SLM_DSW::sw_wgt) = 0.0;
                } else {
                    dsw_vars(i, j, lsm_k, SLM_DSW::sw_wgt) = 1.0;
                }
            }
            // when there is a standing water (puddle) on top of soil, no soil top layer evaporation - gSAM + SLM 
            amrex::Real evap = 0.0;
            if (mws_arr(i, j, 0) > 0.) 
            {
                evap = 0.0;
            }
            else
            {   
                evap = 1.0;
            }   

            aa = 0.0;
            cc = -2.0 * dsw_vars(i, j, d_khi_lsm, SLM_DSW::sh_eff_cond)*dt/(dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm)*(dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm) + dsw_vars(i, j, d_khi_lsm - 1, SLM_DSW::sdepth_mm)));
            bb = 1.0 - cc + dsw_vars(i, j, d_khi_lsm, SLM_DSW::sh_eff_vel)*dt/dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm);
            dd = soilw_arr(i, j, d_khi_lsm) - (std::max(0.,evp_soil_arr(i, j, 0))*evap + soil_transp_frac_arr(i, j, d_khi_lsm)*evapo_dry_arr(i, j, 0) - precip_in)*dt / poro_soil_arr(i, j, d_khi_lsm) / dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm);

            dsw_vars(i, j, d_khi_lsm, SLM_DSW::alpha) = cc / bb;
            dsw_vars(i, j, d_khi_lsm, SLM_DSW::beta) = dd / bb;

            // For 2-(nsoil-1) layer:
            for (int k = 1; k < d_nz_lsm - 1; k++) {
                // NOTE: SLM k-1 indices here are switched from original SAM
                // version because our soil level indices are from -1 (top) to
                // -nsoil (bottom)
                // Index mapping:
                //  SAM SLM  -> ERF SLM
                //   1       ->  d_khi_lsm
                //   k - 1   ->  k + 1
                //   k + 1   ->  k - 1
                //   nsoil   ->  d_klo_lsm
                const int lsm_k = d_khi_lsm - k; // lsm_k = -2 to -nsoil+1

                const amrex::Real sd_lo = dsw_vars(i, j, lsm_k - 1, SLM_DSW::sdepth_mm);
                const amrex::Real sd    = dsw_vars(i, j, lsm_k    , SLM_DSW::sdepth_mm);
                const amrex::Real sd_up = dsw_vars(i, j, lsm_k + 1, SLM_DSW::sdepth_mm);

                const amrex::Real sh_eff_cond    = dsw_vars(i, j, lsm_k    , SLM_DSW::sh_eff_cond);
                const amrex::Real sh_eff_cond_up = dsw_vars(i, j, lsm_k + 1, SLM_DSW::sh_eff_cond);

                const amrex::Real sh_eff_vel    = dsw_vars(i, j, lsm_k    , SLM_DSW::sh_eff_vel);
                const amrex::Real sh_eff_vel_up = dsw_vars(i, j, lsm_k + 1, SLM_DSW::sh_eff_vel);

                aa = -1.0 * dt / sd * (sh_eff_vel_up + sh_eff_cond_up * 2.0 / (sd_up + sd));
                cc = -1.0 * dt / sd * 2.0 * sh_eff_cond / (sd + sd_lo);
                bb = 1.0 - cc + dt / sd * sh_eff_vel + 2.0 * dt / sd * sh_eff_cond_up / (sd_up + sd);
                dd = soilw_arr(i, j, lsm_k) - soil_transp_frac_arr(i, j, lsm_k) * evapo_dry_arr(i, j, 0) * dt / poro_soil_arr(i, j, lsm_k) / sd; // current time step

                dsw_vars(i, j, lsm_k, SLM_DSW::alpha) = cc / (bb - aa * dsw_vars(i, j, lsm_k + 1, SLM_DSW::alpha));
                dsw_vars(i, j, lsm_k, SLM_DSW::beta) = (dd - aa * dsw_vars(i, j, lsm_k + 1, SLM_DSW::beta)) / (bb - aa * dsw_vars(i, j, lsm_k + 1, SLM_DSW::alpha));
            }

            // For bottom layer:
            aa = -1.0 * dt / dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm) *
                 (dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::sh_eff_vel) +
                  dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::sh_eff_cond) * 2.0 /
                    (dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::sdepth_mm) +
                     dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm)));
            cc = 0.0;
            bb =
              1.0 + (dt / dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm) * 2.0 *
                     dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::sh_eff_cond) /
                     (dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::sdepth_mm) +
                      dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm)));

            // drainage when it exceeds 1.0 mm/s
            amrex::Real drainage_flux = std::max(soilw_arr(i, j, d_klo_lsm) - 1.0, 0.0)*poro_soil_arr(i, j, d_klo_lsm)*dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm)/dt;
            dd = soilw_arr(i, j, d_klo_lsm) - (soil_transp_frac_arr(i, j, d_klo_lsm)*evapo_dry_arr(i, j, 0) + drainage_flux) * dt / poro_soil_arr(i, j, d_klo_lsm) / dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm);

            slm_diag_arr(i, j, 0, SLM_Diag::drain_flux) = drainage_flux;

            dsw_vars(i, j, d_klo_lsm, SLM_DSW::alpha) = 0.0;
            dsw_vars(i, j, d_klo_lsm, SLM_DSW::beta) = (dd - aa * dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::beta)) / (bb - aa * dsw_vars(i, j, d_klo_lsm + 1, SLM_DSW::alpha));

            // (n+1) time step soil wetness:
            soilw_arr(i, j, d_klo_lsm) = dsw_vars(i, j, d_klo_lsm, SLM_DSW::beta);
            for (int k = d_nz_lsm - 2; k >= 0; k--) {
                const int lsm_k = d_khi_lsm - k;
                soilw_arr(i, j, lsm_k) = std::max(0.0, dsw_vars(i, j, lsm_k, SLM_DSW::beta) - dsw_vars(i, j, lsm_k, SLM_DSW::alpha) * soilw_arr(i, j, lsm_k - 1));
            }

            bool fix = false;
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;
                if (soilw_arr(i, j, lsm_k) > 1.0)
                {
                    fix = true;
                    break;
                }
            }

            amrex::Real excess_water = 0.0;
            if (fix) {
                // fix the levels where wetness exceeds 1 preserving total water:
                // compute the excess
                for (int kk = 0; kk < d_nz_lsm; kk++) {
                    const int lsm_kk = d_khi_lsm - kk;
                    if (soilw_arr(i, j, lsm_kk) > 1.0)
                    {
                        excess_water += (soilw_arr(i, j, lsm_kk) - 1.0)*poro_soil_arr(i, j, lsm_kk)*dsw_vars(i, j, lsm_kk, SLM_DSW::sdepth_mm);
                        soilw_arr(i, j, lsm_kk) = 1.0;
                    }
                }

                // distribute the excess among layer into deepest layers first:
                for (int kk = d_nz_lsm - 1; kk >= 0; kk--) {
                    const int lsm_kk = d_khi_lsm - kk;
                    if (soilw_arr(i, j, lsm_kk) < 1.0)
                    {
                        cc = std::min(excess_water, (1.0 - soilw_arr(i, j, lsm_kk))*poro_soil_arr(i, j, lsm_kk)*dsw_vars(i, j, lsm_kk, SLM_DSW::sdepth_mm));
                        soilw_arr(i, j, lsm_kk) += cc / (poro_soil_arr(i, j, lsm_kk)*dsw_vars(i, j, lsm_kk, SLM_DSW::sdepth_mm));
                        excess_water -= cc;
                        if (excess_water <= 0.0)
                        {
                            break;
                        }
                    }
                }
            }
            // fixing the issue with rain infiltration even if soil is completely saturated
            // move all the access to mws and modify the precip_in accordinally
            if (excess_water > 0.0) {
                // still some water left after saturating all the soil layers
                // modify the precip_in so that the access of water is moved
                // to the surface water
                precip_in -= excess_water/dt;
            }
            slm_diag_arr(i, j, 0, SLM_Diag::precip_in) = precip_in;
            mws_arr(i, j, 0) = std::max(0., mws_arr(i, j, 0) + (precip_sfc - precip_in)*dt);

            amrex::Real drain = 0.0;
            if(mws_arr(i, j, 0) > mws_mx_arr(i, j, 0))
            {
                drain = (mws_arr(i, j, 0) - mws_mx_arr(i, j, 0))/dt;
                mws_arr(i, j, 0) = mws_mx_arr(i, j, 0);
            }
            else
            { 
                drain = 0.;
            }
            slm_diag_arr(i, j, 0, SLM_Diag::drain) = drain;

        }
    });
}

void SLM::soil_temperature(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const int d_nz_lsm = m_nz_lsm;
    const Real dt = m_dt;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);
    auto landtype_arr = landtype.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);

    auto sst_capa_arr = lsm_fab_vars[LsmVar_SLM::sst_capa]->const_array(mfi);
    auto sst_cond_arr = lsm_fab_vars[LsmVar_SLM::sst_cond]->const_array(mfi);
    auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->const_array(mfi);

    auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->const_array(mfi);

    auto slm_diag_arr = slm_diag.array(mfi);

    auto dst_vars = soilt_vars.array(mfi);
        
    auto mws_arr = mws.array(mfi);

    // TODO: Refactor this whole loop for GPU
    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        amrex::Real temp, k_dry, k_sat, Ke;

        amrex::Real grflux0 = slm_diag_arr(i, j, 0, SLM_Diag::grflux);
        slm_diag_arr(i, j, 0, SLM_Diag::grflux) = grflux0;
        if (landtype_arr(i, j, 0) == 15) {
            // ice
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;
                dst_vars(i, j, lsm_k, SLM_DST::st_cond) = 1.6; // thermal conductivity of ice
                dst_vars(i, j, lsm_k, SLM_DST::st_capa) = 917.0 * 2030.0; // ice heat capacity
            }
        } else {
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;

                temp = (1.0 - poro_soil_arr(i, j, lsm_k)) * 2700.0;

                // dry thermal conductivity [W/mK]
                k_dry = (0.135 * temp + 64.7) / (2700.0 - 0.947*temp);

                // saturated soil thermal conductivity
                k_sat = std::pow(sst_cond_arr(i, j, lsm_k), 1.0 - poro_soil_arr(i, j, lsm_k));
                if (soilt_arr(i, j, lsm_k) > tfriz) {
                    k_sat *= std::pow(0.57, poro_soil_arr(i, j, lsm_k)); // 0.57 = cond_water
                }
                else {
                    k_sat *= std::pow(1.60, poro_soil_arr(i, j, lsm_k)); // 1.6 = cond_ice
                }

                // Weighing factor between dry and saturated soil thermal conductivity
                Ke = log10(std::max(0.1, soilw_arr(i, j, lsm_k))) + 1.0;

                // Total soil thermal conductivity at each node_z
                dst_vars(i, j, lsm_k, SLM_DST::st_cond) = Ke * (k_sat - k_dry) + k_dry;

                // Soil volumetric heaet capacity at each node_z depth
                if (soilt_arr(i, j, lsm_k) > tfriz) {
                    dst_vars(i, j, lsm_k, SLM_DST::st_capa) =
                        (1.0 - poro_soil_arr(i, j, lsm_k)) *
                        sst_capa_arr(i, j, lsm_k) +
                        (998.0 * 4182.0) * soilw_arr(i, j, lsm_k) *
                        poro_soil_arr(i, j, lsm_k); // water heat capacity
                } else {
                    dst_vars(i, j, lsm_k, SLM_DST::st_capa) =
                        (1.0 - poro_soil_arr(i, j, lsm_k)) *
                        sst_capa_arr(i, j, lsm_k) +
                        (917.0 * 2030.0) * soilw_arr(i, j, lsm_k) *
                        poro_soil_arr(i, j, lsm_k); // ice heat capacity
                }
            }
        }

        for (int k = 0; k < d_nz_lsm - 1; k++) {
            const int lsm_k = d_khi_lsm - k;

            // calculate effective conductivity at the adjacent soil layer interface
            dst_vars(i, j, lsm_k, SLM_DST::st_eff_cond) =
                (dst_vars(i, j, lsm_k - 1, SLM_DST::st_cond) * s_depth_arr(i, j, lsm_k - 1) +
                 dst_vars(i, j, lsm_k, SLM_DST::st_cond) * s_depth_arr(i, j, lsm_k)) /
                 (s_depth_arr(i, j, lsm_k - 1) + s_depth_arr(i, j, lsm_k));
        }

        // from FDE by the implicit, Thomas algorithm is applied
        //   aa: terms related with T(k-1)
        //   bb: terms related with T(k)
        //   cc: terms related with T(k+1)
        //   dd: current soil temperature + (additional source/sink on soil top)
        // For first layer:
        amrex::Real aa = 0.0;
        amrex::Real cc =
          -2.0 * dst_vars(i, j, d_khi_lsm, SLM_DST::st_eff_cond) * dt /
          dst_vars(i, j, d_khi_lsm, SLM_DST::st_capa) /
          (s_depth_arr(i, j, d_khi_lsm) *
           (s_depth_arr(i, j, d_khi_lsm) + s_depth_arr(i, j, d_khi_lsm - 1)));
        amrex::Real bb = 1.0 - cc;
        amrex::Real dd = soilt_arr(i, j, d_khi_lsm) - grflux0 * dt / s_depth_arr(i, j, d_khi_lsm) / dst_vars(i, j, d_khi_lsm, SLM_DST::st_capa);

        dst_vars(i, j, d_khi_lsm, SLM_DST::alpha) = cc / bb;
        dst_vars(i, j, d_khi_lsm, SLM_DST::beta) = dd / bb;

        // For 2-(nsoil-1) layer:
        for (int k = 1; k < d_nz_lsm - 1; k++) {
            // NOTE: SLM k-1 indices here are switched from original SAM
            // version because our soil level indices are from -1 (top) to
            // -nsoil (bottom)
            // Index mapping:
            //  SAM SLM  -> ERF SLM
            //   1       ->  d_khi_lsm
            //   k - 1   ->  k + 1
            //   k + 1   ->  k - 1
            //   nsoil   ->  d_klo_lsm
            const int lsm_k = d_khi_lsm - k;

            const amrex::Real st_eff_cond_up = dst_vars(i, j, lsm_k + 1, SLM_DST::st_eff_cond);
            const amrex::Real st_eff_cond    = dst_vars(i, j, lsm_k    , SLM_DST::st_eff_cond);
            const amrex::Real st_capa        = dst_vars(i, j, lsm_k,     SLM_DST::st_capa);

            aa = -2.0 * st_eff_cond_up * dt / s_depth_arr(i, j, lsm_k) / st_capa / (s_depth_arr(i, j, lsm_k+1) + s_depth_arr(i, j, lsm_k));
            cc = -2.0 * st_eff_cond * dt / s_depth_arr(i, j, lsm_k) / st_capa / (s_depth_arr(i, j, lsm_k) + s_depth_arr(i, j, lsm_k-1));
            bb = 1.0 - cc - aa;
            dd = soilt_arr(i, j, lsm_k); // current time step

            dst_vars(i, j, lsm_k, SLM_DST::alpha) = cc / (bb - aa * dst_vars(i, j, lsm_k + 1, SLM_DST::alpha));
            dst_vars(i, j, lsm_k, SLM_DST::beta) = (dd - aa * dst_vars(i, j, lsm_k + 1, SLM_DST::beta)) / (bb - aa * dst_vars(i, j, lsm_k + 1, SLM_DST::alpha));
        }

        // For bottom layer:
        aa = -2.0 * dst_vars(i, j, d_klo_lsm + 1, SLM_DST::st_eff_cond) * dt / s_depth_arr(i, j, d_klo_lsm) / dst_vars(i, j, d_klo_lsm, SLM_DST::st_capa) / (s_depth_arr(i, j, d_klo_lsm + 1) + s_depth_arr(i, j, d_klo_lsm));
        cc = 0.0;
        bb = 1.0 - aa;
        dd = soilt_arr(i, j, d_klo_lsm);

        dst_vars(i, j, d_klo_lsm, SLM_DST::alpha)  = 0.0;
        dst_vars(i, j, d_klo_lsm, SLM_DST::beta) = (dd - aa * dst_vars(i, j, d_klo_lsm + 1, SLM_DST::beta)) / (bb - aa * dst_vars(i, j, d_klo_lsm + 1, SLM_DST::alpha));

        // (n+1) time step soil temperature:
        soilt_arr(i, j, d_klo_lsm) = dst_vars(i, j, d_klo_lsm, SLM_DST::beta);
        for (int k = d_nz_lsm - 2; k >= 0; k--) {
            const int lsm_k = d_khi_lsm - k;
            soilt_arr(i, j, lsm_k) = dst_vars(i, j, lsm_k, SLM_DST::beta) - dst_vars(i, j, lsm_k, SLM_DST::alpha) * soilt_arr(i, j, lsm_k - 1);
        }

        // for wetland, surface water is already implied: gSAM+SLM
        if (landtype_arr(i, j, 0) == 11) {
            mws_arr(i, j, 0) = 0.;
        }   
    });
}

void SLM::soil_nudging(const amrex::MFIter &mfi3d)
{
    if (!dosoiltnudging && !dosoilwnudging) return;

    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const Real dt = m_dt;

    auto box = mfi3d.tilebox();

    auto landmask_arr = landmask.const_array(mfi3d);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->array(mfi3d);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->array(mfi3d);

    auto soilt_nudge_arr = lsm_fab_vars[LsmVar_SLM::soilt_nudge]->array(mfi3d);
    auto soilw_nudge_arr = lsm_fab_vars[LsmVar_SLM::soilw_nudge]->array(mfi3d);
    auto soil_relax_hgt_arr = lsm_fab_vars[LsmVar_SLM::soil_relax_hgt]->const_array(mfi3d);

    auto soilt_obs_arr = lsm_fab_vars[LsmVar_SLM::soilt_obs]->const_array(mfi3d);
    auto soilw_obs_arr = lsm_fab_vars[LsmVar_SLM::soilw_obs]->const_array(mfi3d);

    const Real d_tau = tausoil;
    if (dosoiltnudging)
    {
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (landmask_arr(i, j, 0) != 1) {
                return;
            }

            soilt_arr(i, j, k) -= (soilt_arr(i, j, k) - soilt_obs_arr(i, j, k))*dt/d_tau*soil_relax_hgt_arr(i, j, k);
            soilt_nudge_arr(i, j, k) = (soilt_arr(i, j, k) - soilt_obs_arr(i, j, k))*soil_relax_hgt_arr(i, j, k) / d_tau;
        });
    }

    if (dosoilwnudging)
    {
        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (landmask_arr(i, j, 0) != 1) {
                return;
            }

            soilw_arr(i, j, k) -= (soilw_arr(i, j, k) - soilw_obs_arr(i, j, k))*dt/d_tau*soil_relax_hgt_arr(i, j, k);
            soilw_nudge_arr(i, j, k) = (soilw_arr(i, j, k) - soilw_obs_arr(i, j, k))*soil_relax_hgt_arr(i, j, k) / d_tau;
        });
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real SLM::fh_calc(const amrex::Real &t, const amrex::Real &mps, const amrex::Real &sw, const amrex::Real &B)
{
    amrex::Real moist_pot1 = std::max(-100000.0, mps / (std::pow(std::max(0.0001, sw), B)) / 1000.0);
    return std::min(1.0, std::exp(moist_pot1*CONST_GRAV/461.0/t));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
SLM::linear_interp(const amrex::Real t0, const amrex::Real t1, const amrex::Real t,
                   const amrex::Real x, const amrex::Real y)
{
  // returns a value that is linearly interpolated between x and y at time t. x
  // is at t=t0, y is at t=t1.
  if (t0 == t1 || t > t1) {
    return y;
  }
  const amrex::Real dt = (t - t0) / (t1 - t0);
  return x + (y - x) * dt;
}

void SLM::Copy_State_to_Lsm(const MultiFab& cons_in, const MultiFab& u_in, const MultiFab& v_in)
{
    int khi = khi_lsm;

    auto tsurf = lsm_fab_vars[LsmVar_SLM::tsurf];

    // Get the temperature, density, pressure at the reference level
    for ( MFIter mfi(*tsurf, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        // Create a box with the same i,j bounds, but only at z = 0
        amrex::Box b2d = box3d;
        b2d.setRange(2, 0);

        auto states_array = cons_in.array(mfi);
        auto u_array = u_in.array(mfi);
        auto v_array = v_in.array(mfi);

        auto tref_array  = lsm_fab_vars[LsmVar_SLM::tref]->array(mfi);
        auto rho_array   = lsm_fab_vars[LsmVar_SLM::dref]->array(mfi);
        auto pres_array  = lsm_fab_vars[LsmVar_SLM::pref]->array(mfi);
        auto qref_array  = lsm_fab_vars[LsmVar_SLM::qref]->array(mfi);
        auto slm_u       = lsm_fab_vars[LsmVar_SLM::uref]->array(mfi);
        auto slm_v       = lsm_fab_vars[LsmVar_SLM::vref]->array(mfi);

        auto slm_dir_sw_vis  = lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->array(mfi);
        auto slm_dir_sw_nir  = lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->array(mfi);
        auto slm_diff_sw_vis = lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->array(mfi);
        auto slm_diff_sw_nir = lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->array(mfi);
        auto slm_lw          = lsm_fab_vars[LsmVar_SLM::lwref]->array(mfi);
        auto slm_zenith      = lsm_fab_vars[LsmVar_SLM::coszrsxy]->array(mfi);

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real qv = states_array(i,j,k,RhoQ1_comp)/states_array(i,j,k,Rho_comp);
            rho_array(i,j,k)   = states_array(i,j,k,Rho_comp);
            pres_array(i,j,k)  = getPgivenRTh(states_array(i,j,k,RhoTheta_comp), qv)/100.;
            tref_array(i,j,k)  = getTgivenRandRTh(states_array(i,j,k,Rho_comp),
                                                  states_array(i,j,k,RhoTheta_comp),
                                                  qv);

            qref_array(i,j,k) = qv;

            amrex::Real u_cc = 0.5 * (u_array(i, j, k) + u_array(i + 1, j, k));
            amrex::Real v_cc = 0.5 * (v_array(i, j, k) + v_array(i, j + 1, k));
            slm_u(i, j, k) = u_cc;
            slm_v(i, j, k) = v_cc;

            // TODO: this is for plotting purposes.. state arrays are at k=0 which is ghost cell for SLM values
            //  SLM AMREX plotfile does not write ghost cells, but NetCDF does - fix?
            rho_array(i, j, khi) = rho_array(i, j, 0);
            pres_array(i, j, khi) = pres_array(i, j, 0);
            tref_array(i, j, khi) = tref_array(i, j, 0);
            qref_array(i, j, khi) = qref_array(i, j, 0);
            slm_u(i, j, khi) = slm_u(i, j, 0);
            slm_v(i, j, khi) = slm_v(i, j, 0);

            slm_dir_sw_vis(i, j, khi) = slm_dir_sw_vis(i, j, 0);
            slm_dir_sw_nir(i, j, khi) = slm_dir_sw_nir(i, j, 0);
            slm_diff_sw_vis(i, j, khi) = slm_diff_sw_vis(i, j, 0);
            slm_diff_sw_nir(i, j, khi) = slm_diff_sw_nir(i, j, 0);
            slm_lw(i, j, khi) = slm_lw(i, j, 0);
            slm_zenith(i, j, khi) = slm_zenith(i, j, 0);
        });
    }

#ifdef ERF_USE_NETCDF
    if (rad_input_file != "")
    {
        int tindex = 0;
        int tindex_next = 0;
        for (int k = 1; k < num_rad_times; k++)
        {
            if (time > rad_times[k])
            {
                tindex = k;
            }
        }

        if (tindex == num_rad_times - 1)
        {
            tindex_next = tindex;
        } else {
            tindex_next = tindex + 1;
        }

        const amrex::Real t0 = rad_times[tindex];
        const amrex::Real t1 = rad_times[tindex_next];

        amrex::Print() << " SLM: time = " << time << ", interpolating fluxes at index " << tindex << " (t0 = " << t0 << " t1 = " << t1 << ")" << std::endl;

        const amrex::Real d_time = time;
        for ( MFIter mfi(rad_input_data); mfi.isValid(); ++mfi) {
            const auto& box3d = mfi.tilebox();

            // Create a box with the same i,j bounds, but only at z = 0
            amrex::Box b2d = box3d;
            b2d.setRange(2, 0);

            auto rad_arr = rad_input_data.const_array(mfi);

            auto slm_dir_sw_vis  = lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->array(mfi);
            auto slm_dir_sw_nir  = lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->array(mfi);
            auto slm_diff_sw_vis = lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->array(mfi);
            auto slm_diff_sw_nir = lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->array(mfi);

            auto slm_lw          = lsm_fab_vars[LsmVar_SLM::lwref]->array(mfi);
            auto slm_zenith      = lsm_fab_vars[LsmVar_SLM::coszrsxy]->array(mfi);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                slm_dir_sw_vis(i, j, k) = linear_interp(t0, t1, d_time, rad_arr(i, j, tindex, 0), rad_arr(i, j, tindex_next, 0));
                slm_dir_sw_nir(i, j, k) = linear_interp(t0, t1, d_time, rad_arr(i, j, tindex, 1), rad_arr(i, j, tindex_next, 1));

                slm_diff_sw_vis(i, j, k) = linear_interp(t0, t1, d_time, rad_arr(i, j, tindex, 2), rad_arr(i, j, tindex_next, 2));
                slm_diff_sw_nir(i, j, k) = linear_interp(t0, t1, d_time, rad_arr(i, j, tindex, 3), rad_arr(i, j, tindex_next, 3));

                slm_zenith(i, j, k) = linear_interp(t0, t1, d_time, rad_arr(i, j, tindex, 4), rad_arr(i, j, tindex_next, 4));
                slm_lw(i, j, k) = linear_interp(t0, t1, d_time, rad_arr(i, j, tindex, 5), rad_arr(i, j, tindex_next, 5));


                // TODO: this is for plotting purposes.. state arrays are at k=0 which is ghost cell for SLM values
                //  SLM AMREX plotfile does not write ghost cells, but NetCDF does - fix?
                slm_dir_sw_vis(i, j, khi) = slm_dir_sw_vis(i, j, 0);
                slm_dir_sw_nir(i, j, khi) = slm_dir_sw_nir(i, j, 0);
                slm_diff_sw_vis(i, j, khi) = slm_diff_sw_vis(i, j, 0);
                slm_diff_sw_nir(i, j, khi) = slm_diff_sw_nir(i, j, 0);
                slm_lw(i, j, khi) = slm_lw(i, j, 0);
                slm_zenith(i, j, khi) = slm_zenith(i, j, 0);
            });
        }
        return;
    }
#endif
}


void
SLM::set_precip_input(const amrex::MultiFab* precip_in)
{
    int khi = khi_lsm;

    auto tsurf = lsm_fab_vars[LsmVar_SLM::tsurf];

    for ( MFIter mfi(*tsurf, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        // Create a box with the same i,j bounds, but only at z = 0
        amrex::Box b2d = box3d;
        b2d.setRange(2, 0);

        auto precip_array = precip_in->const_array(mfi);

        auto slm_precip   = lsm_fab_vars[LsmVar_SLM::precipref]->array(mfi);

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            slm_precip(i, j, k) = precip_array(i, j, k, 0);

            // TODO: this is for plotting purposes.. state arrays are at k=0 which is ghost cell for SLM values
            //  SLM AMREX plotfile does not write ghost cells, but NetCDF does - fix?
            slm_precip(i, j, khi) = slm_precip(i, j, 0);
        });
    }
}

void
SLM::set_terrain_inputs(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& sst_in,
                        const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& lmask_in)
{
    if (!first_step || use_wrfinput)
    {
        return;
    }
    auto tsurf = lsm_fab_vars[LsmVar_SLM::tsurf];
    const int d_khi_lsm = khi_lsm;

    if (sst_in[0] && lmask_in[0]) {
        // Set SLM SST and land mask input from ERF
        for ( MFIter mfi(*tsurf, TileNoZ()); mfi.isValid(); ++mfi) {
            const auto& box3d = mfi.tilebox();

            // Create a box with the same i,j bounds, but only at z = 0
            amrex::Box b2d = box3d;
            b2d.setRange(2, 0);

            auto sst_array = sst_in[0]->array(mfi);
            auto lmask_array = lmask_in[0]->array(mfi);

            auto slm_sst   = sstxy.array(mfi);
            auto slm_lmask   = landmask.array(mfi);

            auto slm_tsk   = lsm_fab_vars[LsmVar_SLM::tsurf]->array(mfi);
            auto slm_tskin = t_skin.array(mfi);
            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                slm_sst(i, j, k) = sst_array(i, j, k, 0);
                slm_tsk(i, j, k) = sst_array(i, j, k, 0);
                slm_tskin(i, j, k) = sst_array(i, j, k, 0);
                slm_lmask(i, j, k) = lmask_array(i, j, k, 0);
            });
        }
    } else {
        // If no terrain is setup, initialize SST to reference temperature and set default land mask
        landmask.setVal(1);
        lsm_fab_vars[LsmVar_SLM::tsurf]->setVal(st0[0]);
    }
}


void SLM::Copy_Lsm_to_State(MultiFab& cons_in)
{
    for ( amrex::MFIter mfi(cons_in,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        auto states_arr = cons_in.array(mfi);

        auto tref_array  = lsm_fab_vars[LsmVar_SLM::tref]->array(mfi);
        auto rho_array   = lsm_fab_vars[LsmVar_SLM::dref]->array(mfi);
        auto pres_array  = lsm_fab_vars[LsmVar_SLM::pref]->array(mfi);

        // get potential total density, temperature, qt, qp
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // TODO
            /*
            states_arr(i,j,k,RhoTheta_comp) = rho_arr(i,j,k)*theta_arr(i,j,k);
            states_arr(i,j,k,RhoQ1_comp)    = rho_arr(i,j,k)*qv_arr(i,j,k);
            states_arr(i,j,k,RhoQ2_comp)    = rho_arr(i,j,k)*qc_arr(i,j,k);
            states_arr(i,j,k,RhoQ3_comp)    = rho_arr(i,j,k)*qp_arr(i,j,k);
            */
        });
    }

    // Fill interior ghost cells and periodic boundaries
    cons_in.FillBoundary(m_geom.periodicity());
}

/**
 * Read columns of data from a file, returning each column in a vector.
 */
std::vector<std::vector<amrex::Real>> SLM::read_cols(const std::string &fname, const int skip_nlines)
{
    std::ifstream ifs(fname);
    if (!ifs.is_open())
    {
        amrex::Error("Error opening input file " + fname);
    }

    std::vector<std::vector<amrex::Real>> datasets;
    std::string line;
    int i = 0;
    int ncols = -1;

    while (std::getline(ifs, line))
    {
        i++;
        if (i <= skip_nlines) continue;

        std::istringstream iss(line);

        amrex::Real tmp;
        // Get the number of columns in the file
        if (ncols == -1) {
            int j = 0;
            while (iss >> tmp) {
                datasets.push_back(std::vector<amrex::Real>());
                j+= 1;
            }

            ncols = j;
            iss = std::istringstream(line);

            amrex::Print() << "-> got " << std::to_string(j) << " columns\n";
        }

        int j = 0;
        while (iss >> tmp) {
            // verify each line has the same number of columns
            if (j > ncols) {
              amrex::Error(
                "Error reading file '" + fname + "': expected line " +
                std::to_string(i) + " to have " + std::to_string(ncols) +
                " columns, but got " + std::to_string(j));
            }
            datasets[j].push_back(tmp);
            j+= 1;
        }
    }

    ifs.close();

    return datasets;
}

void SLM::writeSLM_Data(const PlotFileType plotfile_type, const amrex::Real time, const std::string plot_prefix, const int level_step, const int lev, const int finest_lev, amrex::MultiFab &fab, amrex::Geometry &geom, amrex::Vector<std::string> &varnames)
{
    geom.define(amrex::makeSlab(m_lsm_geom.Domain(), 2, 0), m_lsm_geom.ProbDomain(), m_lsm_geom.Coord(), m_lsm_geom.isPeriodic());

    amrex::Vector<amrex::MultiFab*> mf_data;

    mf_data.push_back(&net_rad);

    mf_data.push_back(&mw);
    mf_data.push_back(&mws);
    mf_data.push_back(&t_canop);
    mf_data.push_back(&t_skin);
    mf_data.push_back(&t_ground_skin);
    mf_data.push_back(&t_cas);
    mf_data.push_back(&q_cas);
    mf_data.push_back(&mw_inc);

    mf_data.push_back(&evapo_dry);
    mf_data.push_back(&shf_air);
    mf_data.push_back(&shf_canop);
    mf_data.push_back(&shf_soil);
    mf_data.push_back(&lhf_air);
    mf_data.push_back(&lhf_canop);
    mf_data.push_back(&lhf_soil);

    //mf_data.push_back(&ustar);
    //mf_data.push_back(&tstar);

    mf_data.push_back(&r_a);
    mf_data.push_back(&r_b);
    mf_data.push_back(&r_c);
    mf_data.push_back(&r_d);
    mf_data.push_back(&r_soil);

    mf_data.push_back(&wet_canop);

    mf_data.push_back(&albedovis_v);
    mf_data.push_back(&albedovis_s);
    mf_data.push_back(&albedonir_v);
    mf_data.push_back(&albedonir_s);
    mf_data.push_back(&IR_emis_vege);
    mf_data.push_back(&IR_emis_soil);
    mf_data.push_back(&zrefxy);
    mf_data.push_back(&vege_YES);

    mf_data.push_back(&LAI);
    mf_data.push_back(&SAI);

    MultiFab tmp_landtype = amrex::ToMultiFab(landtype);
    MultiFab tmp_landmask = amrex::ToMultiFab(landmask);

    mf_data.push_back(&tmp_landtype);
    mf_data.push_back(&tmp_landmask);

    IntVect ng(0, 0, 0);

    // Total number of output MFs: net_rad components + mf_data size - 1 + diag vars + olen
    const int output_size = SLM_NetRad::NumVars + mf_data.size() - 1 + SLM_Diag::NumVars + 1;
    fab.define(ba_lsm_2d, net_rad.DistributionMap(), output_size, ng);
    MultiFab::Copy(fab, *(mf_data[0]), 0, 0, SLM_NetRad::NumVars, 0);
    for (int i = 1; i < mf_data.size(); i++)
    {
        MultiFab::Copy(fab, *(mf_data[i]), 0, i + SLM_NetRad::NumVars - 1, 1, 0);
    }
    MultiFab::Copy(fab, slm_diag, 0, output_size - SLM_Diag::NumVars - 1, SLM_Diag::NumVars, 0);
    MultiFab::Copy(fab, *(lsm_fab_flux[LsmFlux_SLM::olen]), 0, output_size - 1, 1, 0);


    varnames = amrex::Vector<std::string>();
    // net_rad component names:
    varnames.push_back("net_swup1");
    varnames.push_back("net_swup2");
    varnames.push_back("net_swdn1");
    varnames.push_back("net_swdn2");
    varnames.push_back("net_sw1");
    varnames.push_back("net_sw2");
    varnames.push_back("tir1");
    varnames.push_back("tir2");
    varnames.push_back("net_lwup1");
    varnames.push_back("net_lwup2");
    varnames.push_back("net_lwdn1");
    varnames.push_back("net_lwdn2");
    varnames.push_back("net_lw1");
    varnames.push_back("net_lw2");
    varnames.push_back("net_rad1");
    varnames.push_back("net_rad2");
    // ----------------------------

    varnames.push_back("mw");
    varnames.push_back("mws");
    varnames.push_back("t_canop");
    varnames.push_back("t_skin");
    varnames.push_back("t_ground_skin");
    varnames.push_back("t_cas");
    varnames.push_back("q_cas");
    varnames.push_back("mw_inc");

    varnames.push_back("evapo_dry");
    varnames.push_back("shf_air");
    varnames.push_back("shf_canop");
    varnames.push_back("shf_soil");
    varnames.push_back("lhf_air");
    varnames.push_back("lhf_canop");
    varnames.push_back("lhf_soil");

    //varnames.push_back("ustar");
    //varnames.push_back("tstar");

    varnames.push_back("r_a");
    varnames.push_back("r_b");
    varnames.push_back("r_c");
    varnames.push_back("r_d");
    varnames.push_back("r_soil");

    varnames.push_back("wet_canop");

    varnames.push_back("albedovis_veg");
    varnames.push_back("albedovis_soil");
    varnames.push_back("albedonir_veg");
    varnames.push_back("albedonir_soil");
    varnames.push_back("IR_emis_veg");
    varnames.push_back("IR_emis_soil");
    varnames.push_back("zrefxy");
    varnames.push_back("veg_flag");

    varnames.push_back("LAI");
    varnames.push_back("SAI");

    varnames.push_back("vegtype");
    varnames.push_back("landmask");

    for (int i = 0; i < diag_names.size(); i++)
    {
        varnames.push_back(diag_names[i]);
    }

    varnames.push_back("olen");

    AMREX_ALWAYS_ASSERT(varnames.size() == output_size);

    if (plotfile_type == PlotFileType::Amrex) {
        //amrex::WriteSingleLevelPlotfile(plotfilename, fab, varnames, lsm_2d_geom, time, level_step);
#ifdef ERF_USE_NETCDF
        // Temporarily write NetCDF always
        //writeSLM_NetCDF(fab, varnames, time, plot_prefix, level_step);
#endif
#ifdef ERF_USE_NETCDF
    } else if (plotfile_type == PlotFileType::Netcdf) {
        writeSLM_NetCDF(fab, varnames, time, plot_prefix, level_step);
#endif
    } else {
        Abort("Dont know this plot_filetype");
    }
}

//==============================================================================
// NOAHMP RADIATION SUBROUTINES
//==============================================================================

// ----------------------------------------------------------------------
// SUBROUTINE SNOW_AGE
// ----------------------------------------------------------------------
// from BATS
// ----------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SLM::snow_age_noahmp(amrex::Real dt, amrex::Real tg, amrex::Real sneqvo, amrex::Real sneqv,
                          amrex::Real tau0, amrex::Real grain_growth, amrex::Real extra_growth,
                          amrex::Real dirt_soot, amrex::Real swemx,
                          amrex::Real& tauss, amrex::Real& fage)
{
    //input
    //  DT        !main time step (s)
    //  TG        !ground temperature (k)
    //  SNEQVO    !snow mass at last time step(mm)
    //  SNEQV     !snow water per unit ground area (mm)

    //output
    //  FAGE     !snow age

    //input/output
    //  TAUSS      !non-dimensional snow age

    //local
    amrex::Real tage;       //total aging effects
    amrex::Real age1;       //effects of grain growth due to vapor diffusion
    amrex::Real age2;       //effects of grain growth at freezing of melt water
    amrex::Real age3;       //effects of soot
    amrex::Real dela;       //temporary variable
    amrex::Real sge;        //temporary variable
    amrex::Real dels;       //temporary variable
    amrex::Real dela0;      //temporary variable
    amrex::Real arg;        //temporary variable
    // See Yang et al. (1997) J.of Climate for detail.

    constexpr amrex::Real TFRZ = 273.16; // freezing/melting point (k)

    if(sneqv <= 0.0) {
        tauss = 0.0;
    } else {
        dela0 = dt/tau0;
        arg   = grain_growth*(1.0/TFRZ-1.0/tg);
        age1  = std::exp(arg);
        age2  = std::exp(std::min(0.0, extra_growth*arg));
        age3  = dirt_soot;
        tage  = age1+age2+age3;
        dela  = dela0*tage;
        dels  = std::max(0.0, sneqv-sneqvo) / swemx;
        sge   = (tauss+dela)*(1.0-dels);
        tauss = std::max(0.0,sge);
    }

    fage= tauss/(tauss+1.0);
}

// --------------------------------------------------------------------------------------------------
// SUBROUTINE SNOWALB_BATS
// --------------------------------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SLM::snowalb_bats_noahmp(int nband, amrex::Real fsno, amrex::Real cosz, amrex::Real fage,
                              amrex::Real bats_cosz, amrex::Real bats_vis_new, amrex::Real bats_nir_new,
                              amrex::Real bats_vis_age, amrex::Real bats_nir_age,
                              amrex::Real bats_vis_dir, amrex::Real bats_nir_dir,
                              amrex::Real* albsnd, amrex::Real* albsni)
{
    // --------------------------------------------------------------------------------------------------
    // input
    //  NBAND  !number of waveband classes
    //  COSZ    !cosine solar zenith angle
    //  FSNO    !snow cover fraction (-)
    //  FAGE    !snow age correction

    // output
    //  ALBSND !snow albedo for direct(1=vis, 2=nir)
    //  ALBSNI !snow albedo for diffuse
    // ---------------------------------------------------------------------------------------------

    // ------------------------ local variables ----------------------------------------------------
    amrex::Real fzen;                 //zenith angle correction
    amrex::Real cf1;                  //temperary variable
    amrex::Real sl2;                  //2.*SL
    amrex::Real sl1;                  //1/SL
    amrex::Real sl;                   //adjustable parameter
    //  REAL, PARAMETER :: C1 = 0.2  !default in BATS
    //  REAL, PARAMETER :: C2 = 0.5  !default in BATS
    //  REAL, PARAMETER :: C1 = 0.2 * 2. ! double the default to match Sleepers River's
    //  REAL, PARAMETER :: C2 = 0.5 * 2. ! snow surface albedo (double aging effects)
    // ---------------------------------------------------------------------------------------------
    // zero albedos for all points

    albsnd[0] = 0.0;
    albsnd[1] = 0.0;
    albsni[0] = 0.0;
    albsni[1] = 0.0;

    // when cosz > 0

    sl=bats_cosz;
    sl1=1.0/sl;
    sl2=2.0*sl;
    cf1=((1.0+sl1)/(1.0+sl2*cosz)-sl1);
    fzen=std::max(cf1,0.0);

    albsni[0]=bats_vis_new*(1.0-bats_vis_age*fage);
    albsni[1]=bats_nir_new*(1.0-bats_nir_age*fage);

    albsnd[0]=albsni[0]+bats_vis_dir*fzen*(1.0-albsni[0]);    //  vis direct
    albsnd[1]=albsni[1]+bats_nir_dir*fzen*(1.0-albsni[1]);    //  nir direct
}

// --------------------------------------------------------------------------------------------------
// SUBROUTINE GROUNDALB
// --------------------------------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SLM::groundalb_noahmp(int nsoil, int nband, int ice, int ist, amrex::Real fsno,
                           const amrex::Real* smc, const amrex::Real* albsnd, const amrex::Real* albsni,
                           amrex::Real cosz, amrex::Real tg,
                           const amrex::Real* albsat, const amrex::Real* albdry, const amrex::Real* alblak,
                           amrex::Real* albgrd, amrex::Real* albgri)
{
    // --------------------------------------------------------------------------------------------------
    //input
    //  ILOC   !grid index
    //  JLOC   !grid index
    //  NSOIL  !number of soil layers
    //  NBAND  !number of solar radiation waveband classes
    //  ICE    !value of ist for land ice
    //  IST    !surface type
    //  FSNO   !fraction of surface covered with snow (-)
    //  TG     !ground temperature (k)
    //  COSZ   !cosine solar zenith angle (0-1)
    //  SMC    !volumetric soil water content (m3/m3)
    //  ALBSND !direct beam snow albedo (vis, nir)
    //  ALBSNI !diffuse snow albedo (vis, nir)

    //output
    //  ALBGRD !ground albedo (direct beam: vis, nir)
    //  ALBGRI !ground albedo (diffuse: vis, nir)

    //local
    amrex::Real inc;    //soil water correction factor for soil albedo
    amrex::Real albsod; //soil albedo (direct)
    amrex::Real albsoi; //soil albedo (diffuse)
    // --------------------------------------------------------------------------------------------------
    constexpr amrex::Real TFRZ = 273.16; // freezing/melting point (k)

    for (int ib = 0; ib < nband; ib++) {
        inc = std::max(0.11-0.40*smc[0], 0.0);
        if (ist == 1) {                     //soil
            albsod = std::min(albsat[ib]+inc, albdry[ib]);
            albsoi = albsod;
        } else if (tg > TFRZ) {               //unfrozen lake, wetland
            albsod = 0.06/(std::max(0.01,cosz)*std::max(0.01,cosz)*std::max(0.01,cosz)*
                           std::max(0.01,cosz)*std::max(0.01,cosz)*std::max(0.01,cosz)*
                           std::max(0.01,cosz) + 0.15);
            albsoi = 0.06;
        } else {                                      //frozen lake, wetland
            albsod = alblak[ib];
            albsoi = albsod;
        }

        albgrd[ib] = albsod*(1.0-fsno) + albsnd[ib]*fsno;
        albgri[ib] = albsoi*(1.0-fsno) + albsni[ib]*fsno;
    }
}

// --------------------------------------------------------------------------------------------------
// SUBROUTINE TWOSTREAM
// --------------------------------------------------------------------------------------------------
// use two-stream approximation of Dickinson (1983) Adv Geophysics
// 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
// to calculate fluxes absorbed by vegetation, reflected by vegetation,
// and transmitted through vegetation for unit incoming direct or diffuse
// flux given an underlying surface with known albedo.
// --------------------------------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SLM::twostream_noahmp(int ib, int ic, int vegtyp, amrex::Real cosz, amrex::Real vai,
                           amrex::Real fwet, amrex::Real t, const amrex::Real* albgrd, const amrex::Real* albgri,
                           const amrex::Real* rho, const amrex::Real* tau, amrex::Real fveg, int ist,
                           amrex::Real xl, amrex::Real omegas_param, amrex::Real betads, amrex::Real betais,
                           int opt_rad, amrex::Real rc, amrex::Real hvt, amrex::Real hvb, amrex::Real den,
                           amrex::Real* fab, amrex::Real* fre, amrex::Real* ftd, amrex::Real* fti,
                           amrex::Real& gdir, amrex::Real* frev, amrex::Real* freg,
                           amrex::Real& bgap, amrex::Real& wgap,
                           amrex::Real& xl_out, amrex::Real& chil_out, amrex::Real& phi1_out, amrex::Real& phi2_out)
{
    // --------------------------------------------------------------------------------------------------
    // input
    //   IST     !surface type
    //   IB      !waveband number
    //   IC      !0=unit incoming direct; 1=unit incoming diffuse
    //   VEGTYP  !vegetation type
    //   COSZ    !cosine of direct zenith angle (0-1)
    //   VAI     !one-sided leaf+stem area index (m2/m2)
    //   FWET    !fraction of lai, sai that is wetted (-)
    //   T       !surface temperature (k)
    //   ALBGRD  !direct  albedo of underlying surface (-)
    //   ALBGRI  !diffuse albedo of underlying surface (-)
    //   RHO     !leaf+stem reflectance
    //   TAU     !leaf+stem transmittance
    //   FVEG    !green vegetation fraction [0.0-1.0]

    // output
    //   FAB     !flux abs by veg layer (per unit incoming flux)
    //   FRE     !flux refl above veg layer (per unit incoming flux)
    //   FTD     !down dir flux below veg layer (per unit in flux)
    //   FTI     !down dif flux below veg layer (per unit in flux)
    //   GDIR    !projected leaf+stem area in solar direction
    //   FREV    !flux reflected by veg layer   (per unit incoming flux)
    //   FREG    !flux reflected by ground (per unit incoming flux)
    //   BGAP    !between canopy gap fraction for beam (-)
    //   WGAP    !within canopy gap fraction for beam (-)

    // local
    amrex::Real omega;   //fraction of intercepted radiation that is scattered
    amrex::Real omegal;  //omega for leaves
    amrex::Real betai;   //upscatter parameter for diffuse radiation
    amrex::Real betail;  //betai for leaves
    amrex::Real betad;   //upscatter parameter for direct beam radiation
    amrex::Real betadl;  //betad for leaves
    amrex::Real ext;     //optical depth of direct beam per unit leaf area
    amrex::Real avmu;    //average diffuse optical depth
    amrex::Real coszi;   //0.001 <= cosz <= 1.000
    amrex::Real asu;     //single scattering albedo
    amrex::Real chil;    // -0.4 <= xl <= 0.6

    amrex::Real tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
    amrex::Real p1,p2,p3,p4,s1,s2,u1,u2,u3;
    amrex::Real b,c,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10;
    amrex::Real phi1,phi2,sigma;
    amrex::Real ftds,ftis,fres;
    amrex::Real denfveg;
    amrex::Real vai_spread;
    //jref:start
    amrex::Real freveg,frebar,ftdveg,ftiveg,ftdbar,ftibar;
    amrex::Real thetaz;
    //jref:end

    //  variables for the modified two-stream scheme
    //  Niu and Yang (2004), JGR
    constexpr amrex::Real PAI = 3.14159265;
    constexpr amrex::Real TFRZ = 273.16; // freezing/melting point (k)
    amrex::Real hd;       //crown depth (m)
    amrex::Real bb;       //vertical crown radius (m)
    amrex::Real thetap;   //angle conversion from SZA
    amrex::Real fa;       //foliage volume density (m-1)
    amrex::Real newvai;   //effective LSAI (-)
    amrex::Real kopen;    //gap fraction for diffue light (-)
    amrex::Real gap;      //total gap fraction for beam ( <=1-shafac )

    // -----------------------------------------------------------------
    // compute within and between gaps
    vai_spread = vai;
    if(vai == 0.0) {
        gap     = 1.0;
        kopen   = 1.0;
    } else {
        if(opt_rad == 1) {
            denfveg = -std::log(std::max(1.0-fveg,0.01))/(PAI*rc*rc);
            hd      = hvt - hvb;
            bb      = 0.5 * hd;
            thetap  = std::atan(bb/rc * std::tan(std::acos(std::max(0.01,cosz))) );
            // BGAP    = EXP(-parameters%DEN * PAI * parameters%RC**2/COS(THETAP) )
            bgap    = std::exp(-denfveg * PAI * rc*rc/std::cos(thetap) );
            fa      = vai/(1.33 * PAI * rc*rc*rc *(bb/rc)*denfveg);
            newvai  = hd*fa;
            wgap    = (1.0-bgap) * std::exp(-0.5*newvai/cosz);
            gap     = std::min(1.0-fveg, bgap+wgap);

            kopen   = 0.05;
        }

        if(opt_rad == 2) {
            gap     = 0.0;
            kopen   = 0.0;
        }

        if(opt_rad == 3) {
            gap     = 1.0-fveg;
            kopen   = 1.0-fveg;
        }
    }

    // calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
    // OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
    // and OMEGA*BETAI are calculated and then divided by the new OMEGA
    // because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
    // also, the transmittances and reflectances (TAU, RHO) are linear
    // weights of leaf and stem values.

    coszi  = std::max(0.001, cosz);
    chil   = std::min( std::max(xl, -0.4), 0.6);
    if (std::abs(chil) <= 0.01) chil = 0.01;
    phi1   = 0.5 - 0.633*chil - 0.330*chil*chil;
    phi2   = 0.877 * (1.0-2.0*phi1);
    gdir   = phi1 + phi2*coszi;

    // Output diagnostics
    xl_out = xl;
    chil_out = chil;
    phi1_out = phi1;
    phi2_out = phi2;
    ext    = gdir/coszi;
    avmu   = ( 1.0 - phi1/phi2 * std::log((phi1+phi2)/phi1) ) / phi2;
    omegal = rho[ib] + tau[ib];
    tmp0   = gdir + phi2*coszi;
    tmp1   = phi1*coszi;
    asu    = 0.5*omegal*gdir/tmp0 * ( 1.0-tmp1/tmp0*std::log((tmp1+tmp0)/tmp1) );
    betadl = (1.0+avmu*ext)/(omegal*avmu*ext)*asu;
    betail = 0.5 * ( rho[ib]+tau[ib] + (rho[ib]-tau[ib])
                    * ((1.0+chil)/2.0)*((1.0+chil)/2.0) ) / omegal;

    // adjust omega, betad, and betai for intercepted snow

    if (t > TFRZ) {                                //no snow
        tmp0 = omegal;
        tmp1 = betadl;
        tmp2 = betail;
    } else {
        tmp0 =   (1.0-fwet)*omegal        + fwet*omegas_param;
        tmp1 = ( (1.0-fwet)*omegal*betadl + fwet*omegas_param*betads ) / tmp0;
        tmp2 = ( (1.0-fwet)*omegal*betail + fwet*omegas_param*betais ) / tmp0;
    }

    omega = tmp0;
    betad = tmp1;
    betai = tmp2;

    // absorbed, reflected, transmitted fluxes per unit incoming radiation

    b = 1.0 - omega + omega*betai;
    c = omega*betai;
    tmp0 = avmu*ext;
    d = tmp0 * omega*betad;
    f = tmp0 * omega*(1.0-betad);
    tmp1 = b*b - c*c;
    h = std::sqrt(tmp1) / avmu;
    sigma = tmp0*tmp0 - tmp1;
    if ( std::abs(sigma) < 1.0e-6 ) sigma = (sigma >= 0) ? 1.0e-6 : -1.0e-6;
    p1 = b + avmu*h;
    p2 = b - avmu*h;
    p3 = b + tmp0;
    p4 = b - tmp0;
    s1 = std::exp(-h*vai);
    s2 = std::exp(-ext*vai);
    if (ic == 0) {
        u1 = b - c/albgrd[ib];
        u2 = b - c*albgrd[ib];
        u3 = f + c*albgrd[ib];
    } else {
        u1 = b - c/albgri[ib];
        u2 = b - c*albgri[ib];
        u3 = f + c*albgri[ib];
    }
    tmp2 = u1 - avmu*h;
    tmp3 = u1 + avmu*h;
    d1 = p1*tmp2/s1 - p2*tmp3*s1;
    tmp4 = u2 + avmu*h;
    tmp5 = u2 - avmu*h;
    d2 = tmp4/s1 - tmp5*s1;
    h1 = -d*p4 - c*f;
    tmp6 = d - h1*p3/sigma;
    tmp7 = ( d - c - h1/sigma*(u1+tmp0) ) * s2;
    h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1;
    h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1;
    h4 = -f*p3 - c*d;
    tmp8 = h4/sigma;
    tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2;
    h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2;
    h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2;
    h7 = (c*tmp2) / (d1*s1);
    h8 = (-c*tmp3*s1) / d1;
    h9 = tmp4 / (d2*s1);
    h10 = (-tmp5*s1) / d2;

    // downward direct and diffuse fluxes below vegetation
    // Niu and Yang (2004), JGR.

    if (ic == 0) {
        ftds = s2                           *(1.0-gap) + gap;
        ftis = (h4*s2/sigma + h5*s1 + h6/s1)*(1.0-gap);
    } else {
        ftds = 0.0;
        ftis = (h9*s1 + h10/s1)*(1.0-kopen) + kopen;
    }
    ftd[ib] = ftds;
    fti[ib] = ftis;

    // flux reflected by the surface (veg. and ground)

    if (ic == 0) {
        fres   = (h1/sigma + h2 + h3)*(1.0-gap  ) + albgrd[ib]*gap;
        freveg = (h1/sigma + h2 + h3)*(1.0-gap  );
        frebar = albgrd[ib]*gap;                   //jref - separate veg. and ground reflection
    } else {
        fres   = (h7 + h8) *(1.0-kopen) + albgri[ib]*kopen;
        freveg = (h7 + h8) *(1.0-kopen) + albgri[ib]*kopen;
        frebar = 0.0;                                //jref - separate veg. and ground reflection
    }
    fre[ib] = fres;

    frev[ib] = freveg;
    freg[ib] = frebar;

    // flux absorbed by vegetation

    fab[ib] = 1.0 - fre[ib] - (1.0-albgrd[ib])*ftd[ib]
                            - (1.0-albgri[ib])*fti[ib];
}

// --------------------------------------------------------------------------------------------------
// SUBROUTINE SURRAD
// --------------------------------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SLM::surrad_noahmp(amrex::Real mpe, amrex::Real fsun, amrex::Real fsha, amrex::Real elai, amrex::Real vai,
                        amrex::Real laisun, amrex::Real laisha, const amrex::Real* solad, const amrex::Real* solai,
                        const amrex::Real* fabd, const amrex::Real* fabi, const amrex::Real* ftdd,
                        const amrex::Real* ftid, const amrex::Real* ftii, const amrex::Real* albgrd,
                        const amrex::Real* albgri, const amrex::Real* albd, const amrex::Real* albi,
                        const amrex::Real* frevd, const amrex::Real* frevi, const amrex::Real* fregd, const amrex::Real* fregi,
                        amrex::Real& parsun, amrex::Real& parsha, amrex::Real& sav, amrex::Real& sag,
                        amrex::Real& fsa, amrex::Real& fsr, amrex::Real& fsrv, amrex::Real& fsrg)
{
    // --------------------------------------------------------------------------------------------------
    // input
    //  MPE     !prevents underflow errors if division by zero
    //  FSUN    !sunlit fraction of canopy
    //  FSHA    !shaded fraction of canopy
    //  ELAI    !leaf area, one-sided
    //  VAI     !leaf + stem area, one-sided
    //  LAISUN  !sunlit leaf area index, one-sided
    //  LAISHA  !shaded leaf area index, one-sided
    //  SOLAD   !incoming direct solar radiation (w/m2)
    //  SOLAI   !incoming diffuse solar radiation (w/m2)
    //  FABD    !flux abs by veg (per unit incoming direct flux)
    //  FABI    !flux abs by veg (per unit incoming diffuse flux)
    //  FTDD    !down dir flux below veg (per incoming dir flux)
    //  FTID    !down dif flux below veg (per incoming dir flux)
    //  FTII    !down dif flux below veg (per incoming dif flux)
    //  ALBGRD  !ground albedo (direct)
    //  ALBGRI  !ground albedo (diffuse)
    //  ALBD    !overall surface albedo (direct)
    //  ALBI    !overall surface albedo (diffuse)
    //  FREVD    !overall surface albedo veg (direct)
    //  FREVI    !overall surface albedo veg (diffuse)
    //  FREGD    !overall surface albedo grd (direct)
    //  FREGI    !overall surface albedo grd (diffuse)

    // output
    //  PARSUN  !average absorbed par for sunlit leaves (w/m2)
    //  PARSHA  !average absorbed par for shaded leaves (w/m2)
    //  SAV     !solar radiation absorbed by vegetation (w/m2)
    //  SAG     !solar radiation absorbed by ground (w/m2)
    //  FSA     !total absorbed solar radiation (w/m2)
    //  FSR     !total reflected solar radiation (w/m2)
    //  FSRV    !reflected solar radiation by vegetation
    //  FSRG    !reflected solar radiation by ground

    // ------------------------ local variables ----------------------------------------------------
    constexpr int NBAND = 2;   //number of solar radiation waveband classes

    amrex::Real abs;     //absorbed solar radiation (w/m2)
    amrex::Real rnir;    //reflected solar radiation [nir] (w/m2)
    amrex::Real rvis;    //reflected solar radiation [vis] (w/m2)
    amrex::Real laifra;  //leaf area fraction of canopy
    amrex::Real trd;     //transmitted solar radiation: direct (w/m2)
    amrex::Real tri;     //transmitted solar radiation: diffuse (w/m2)
    amrex::Real cad[2];     //direct beam absorbed by canopy (w/m2)
    amrex::Real cai[2];     //diffuse radiation absorbed by canopy (w/m2)
    // ---------------------------------------------------------------------------------------------

    // zero summed solar fluxes

    sag = 0.0;
    sav = 0.0;
    fsa = 0.0;

    // loop over nband wavebands

    for (int ib = 0; ib < NBAND; ib++) {

        // absorbed by canopy

        cad[ib] = solad[ib]*fabd[ib];
        cai[ib] = solai[ib]*fabi[ib];
        sav     = sav + cad[ib] + cai[ib];
        fsa     = fsa + cad[ib] + cai[ib];

        // transmitted solar fluxes incident on ground

        trd = solad[ib]*ftdd[ib];
        tri = solad[ib]*ftid[ib] + solai[ib]*ftii[ib];

        // solar radiation absorbed by ground surface

        abs = trd*(1.0-albgrd[ib]) + tri*(1.0-albgri[ib]);
        sag = sag + abs;
        fsa = fsa + abs;
    }

    // partition visible canopy absorption to sunlit and shaded fractions
    // to get average absorbed par for sunlit and shaded leaves

    laifra = elai / std::max(vai,mpe);
    if (fsun > 0.0) {
        parsun = (cad[0]+fsun*cai[0]) * laifra / std::max(laisun,mpe);
        parsha = (fsha*cai[0])*laifra / std::max(laisha,mpe);
    } else {
        parsun = 0.0;
        parsha = (cad[0]+cai[0])*laifra / std::max(laisha,mpe);
    }

    // reflected solar radiation

    rvis = albd[0]*solad[0] + albi[0]*solai[0];
    rnir = albd[1]*solad[1] + albi[1]*solai[1];
    fsr  = rvis + rnir;

    // reflected solar radiation of veg. and ground (combined ground)
    fsrv = frevd[0]*solad[0]+frevi[0]*solai[0]+frevd[1]*solad[1]+frevi[1]*solai[1];
    fsrg = fregd[0]*solad[0]+fregi[0]*solai[0]+fregd[1]*solad[1]+fregi[1]*solai[1];
}

// --------------------------------------------------------------------------------------------------
// SUBROUTINE ALBEDO
// --------------------------------------------------------------------------------------------------
// surface albedos. also fluxes (per unit incoming direct and diffuse
// radiation) reflected, transmitted, and absorbed by vegetation.
// also sunlit fraction of the canopy.
// --------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------
// SUBROUTINE RADIATION_NOAHMP - Main NOAHMP radiation routine
// --------------------------------------------------------------------------------------------------
void SLM::radiation_noahmp(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const int d_nz_lsm = m_nz_lsm;
    const int d_opt_rad = opt_rad;
    const int d_opt_alb = opt_alb;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);
    auto vegtype_arr = lsm_fab_vars[LsmVar_SLM::vegtype]->const_array(mfi);

    // Input arrays
    auto swdsvisxyref_arr  = lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->const_array(mfi);
    auto swdsnirxyref_arr  = lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->const_array(mfi);
    auto swdsvisdxyref_arr = lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->const_array(mfi);
    auto swdsnirdxyref_arr = lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->const_array(mfi);
    auto coszrsxy_arr      = lsm_fab_vars[LsmVar_SLM::coszrsxy]->const_array(mfi);
    auto precipref_arr     = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

    auto LAI_arr   = LAI.const_array(mfi);
    auto SAI_arr   = SAI.const_array(mfi);
    auto t_canop_arr = t_canop.const_array(mfi);
    auto t_ground_skin_arr = t_ground_skin.const_array(mfi);
    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->const_array(mfi);
    auto veg_frac_arr = lsm_fab_vars[LsmVar_SLM::veg_frac]->const_array(mfi);

    // Existing SLM vegetation structure parameters (used instead of reading from NoahmpTable.TBL)
    // NOTE: Khai_L is equivalent to NOAHMP's XL (leaf/stem orientation index)
    // NOTE: ztop is equivalent to NOAHMP's HVT (top of canopy height in meters)
    auto Khai_L_arr = Khai_L.const_array(mfi);  // leaf/stem orientation (NOAHMP name: XL)
    auto ztop_arr = ztop.const_array(mfi);      // canopy height (NOAHMP name: HVT)

    // State variables for NOAHMP radiation
    auto albold_arr = albold_noahmp.array(mfi);
    auto tauss_arr  = tauss_noahmp.array(mfi);

    // Net radiation arrays for output
    auto net_rad_arr = net_rad.array(mfi);
    auto t_skin_arr  = t_skin.array(mfi);

    // Incoming longwave radiation
    auto lwref_arr = lsm_fab_vars[LsmVar_SLM::lwref]->const_array(mfi);

    // Outputs
    auto emis_sfc_arr     = lsm_fab_vars[LsmVar_SLM::emis_sfc]->array(mfi);
    auto alb_nir_sfc_arr  = lsm_fab_vars[LsmVar_SLM::alb_nir_sfc]->array(mfi);
    auto alb_vis_sfc_arr  = lsm_fab_vars[LsmVar_SLM::alb_vis_sfc]->array(mfi);
    auto alb_nir_sfc_diff_arr  = lsm_fab_vars[LsmVar_SLM::alb_nir_sfc_diff]->array(mfi);
    auto alb_vis_sfc_diff_arr  = lsm_fab_vars[LsmVar_SLM::alb_vis_sfc_diff]->array(mfi);

    amrex::Real dt_loc = m_dt;

    // Get radiation parameters on device
    const amrex::Real* d_rhol_vis = rhol_vis_rad.data();
    const amrex::Real* d_rhol_nir = rhol_nir_rad.data();
    const amrex::Real* d_rhos_vis = rhos_vis_rad.data();
    const amrex::Real* d_rhos_nir = rhos_nir_rad.data();
    const amrex::Real* d_taul_vis = taul_vis_rad.data();
    const amrex::Real* d_taul_nir = taul_nir_rad.data();
    const amrex::Real* d_taus_vis = taus_vis_rad.data();
    const amrex::Real* d_taus_nir = taus_nir_rad.data();
    // NOTE: d_xl uses existing Khai_L_arr, d_hvt uses existing ztop_arr
    const amrex::Real* d_rc = rc_rad.data();
    const amrex::Real* d_hvb = hvb_rad.data();
    const amrex::Real* d_den = den_rad.data();

    const amrex::Real* d_albsat_vis = albsat_vis.data();
    const amrex::Real* d_albsat_nir = albsat_nir.data();
    const amrex::Real* d_albdry_vis = albdry_vis.data();
    const amrex::Real* d_albdry_nir = albdry_nir.data();
    const amrex::Real* d_alblak = alblak_rad.data();
    const amrex::Real* d_omegas = omegas_rad.data();

    const amrex::Real d_betads = betads_rad;
    const amrex::Real d_betais = betais_rad;
    const amrex::Real d_tau0 = tau0_rad;
    const amrex::Real d_grain_growth = grain_growth_rad;
    const amrex::Real d_extra_growth = extra_growth_rad;
    const amrex::Real d_dirt_soot = dirt_soot_rad;
    const amrex::Real d_swemx = swemx_rad;
    const amrex::Real d_bats_cosz = bats_cosz_rad;
    const amrex::Real d_bats_vis_new = bats_vis_new_rad;
    const amrex::Real d_bats_nir_new = bats_nir_new_rad;
    const amrex::Real d_bats_vis_age = bats_vis_age_rad;
    const amrex::Real d_bats_nir_age = bats_nir_age_rad;
    const amrex::Real d_bats_vis_dir = bats_vis_dir_rad;
    const amrex::Real d_bats_nir_dir = bats_nir_dir_rad;

    // Emissivity parameters for longwave radiation
    const amrex::Real d_snow_emis = snow_emis_rad;
    const amrex::Real d_eg_soil = eg_soil_rad;
    const amrex::Real d_eg_lake = eg_lake_rad;

    // Stefan-Boltzmann constant
    constexpr amrex::Real SB = 5.67e-08;  // W/m²/K⁴

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        // --------------------------------------------------------------------------------------------------
        // Local variables for radiation calculation (all arrays from NOAHMP RADIATION subroutine)
        // --------------------------------------------------------------------------------------------------
        constexpr amrex::Real MPE = 1.0E-6;
        constexpr int NBAND = 2;
        constexpr amrex::Real TFRZ = 273.16;

        amrex::Real fage;   //snow age function (0 - new snow)
        amrex::Real albgrd[2]; //ground albedo (direct)
        amrex::Real albgri[2]; //ground albedo (diffuse)
        amrex::Real albd[2];   //surface albedo (direct)
        amrex::Real albi[2];   //surface albedo (diffuse)
        amrex::Real fabd[2];   //flux abs by veg (per unit direct flux)
        amrex::Real fabi[2];   //flux abs by veg (per unit diffuse flux)
        amrex::Real ftdd[2];   //down direct flux below veg (per unit dir flux)
        amrex::Real ftid[2];   //down diffuse flux below veg (per unit dir flux)
        amrex::Real ftii[2];   //down diffuse flux below veg (per unit dif flux)
        amrex::Real frevd[2], frevi[2], fregd[2], fregi[2];
        amrex::Real albsnd[2];   //snow albedo (direct)
        amrex::Real albsni[2];   //snow albedo (diffuse)
        amrex::Real ftdi[2];     //down direct flux below veg per unit dif flux = 0

        amrex::Real fsha;   //shaded fraction of canopy
        amrex::Real vai;    //total LAI + stem area index, one sided
        amrex::Real fsun;   //sunlit fraction of canopy (-)
        amrex::Real laisun; //sunlit leaf area (-)
        amrex::Real laisha; //shaded leaf area (-)
        amrex::Real parsun; //average absorbed par for sunlit leaves (w/m2)
        amrex::Real parsha; //average absorbed par for shaded leaves (w/m2)
        amrex::Real sav;    //solar radiation absorbed by vegetation (w/m2)
        amrex::Real sag;    //solar radiation absorbed by ground (w/m2)
        amrex::Real fsa;    //total absorbed solar radiation (w/m2)
        amrex::Real fsr;    //total reflected solar radiation (w/m2)
        amrex::Real fsrv;   //veg. reflected solar radiation (w/m2)
        amrex::Real fsrg;   //ground reflected solar radiation (w/m2)
        amrex::Real bgap, wgap;
        amrex::Real gdir;   //average projected leaf/stem area in solar direction
        amrex::Real ext;    //optical depth direct beam per unit leaf + stem area
        amrex::Real xl_diag_val, chil_diag_val, phi1_diag_val, phi2_diag_val; //diagnostic outputs

        amrex::Real rho[2];      //leaf/stem reflectance weighted by fraction LAI and SAI
        amrex::Real tau[2];      //leaf/stem transmittance weighted by fraction LAI and SAI
        amrex::Real wl, ws;      //fraction of LAI+SAI that is LAI/SAI

        // Get input values
        amrex::Real cosz = coszrsxy_arr(i, j, 0);
        amrex::Real elai = LAI_arr(i, j, 0);
        amrex::Real esai = SAI_arr(i, j, 0);
        amrex::Real tv = t_canop_arr(i, j, 0);
        amrex::Real tg = t_ground_skin_arr(i, j, 0);
        amrex::Real fveg = veg_frac_arr(i, j, d_khi_lsm);
        int vegtyp = static_cast<int>(vegtype_arr(i, j, d_khi_lsm));
        int veg_idx = vegtyp - 1;  // convert to 0-based index
        if (veg_idx < 0) veg_idx = 0;

        // Incoming solar radiation
        amrex::Real solad[2], solai[2];
        solad[0] = swdsvisxyref_arr(i, j, 0);   // direct visible
        solad[1] = swdsnirxyref_arr(i, j, 0);   // direct NIR
        solai[0] = swdsvisdxyref_arr(i, j, 0);  // diffuse visible
        solai[1] = swdsnirdxyref_arr(i, j, 0);  // diffuse NIR

        // Snow variables (set to zero for now - can be added later)
        amrex::Real fsno = 0.0;     // snow cover fraction
        amrex::Real snowh = 0.0;    // snow height (mm)
        amrex::Real sneqvo = 0.0;   // snow mass at last time step (mm)
        amrex::Real sneqv = 0.0;    // snow mass (mm)
        amrex::Real qsnow = 0.0;    // snowfall (mm/s)

        // Wetness fraction (simplified - using precipitation as indicator)
        amrex::Real fwet = (precipref_arr(i, j, 0) > 0.0) ? 0.1 : 0.0;

        // Soil moisture for top layer
        amrex::Real smc[1];
        smc[0] = soilw_arr(i, j, d_khi_lsm) * poro_soil_arr(i, j, d_khi_lsm);

        int nsoil = 1;  // using top layer only for albedo calc
        int ice = 0;    // not ice
        int ist = 1;    // soil surface type

        // Initialize outputs
        bgap = 0.0;
        wgap = 0.0;
        for (int ib = 0; ib < NBAND; ib++) {
            albgrd[ib] = 0.0;
            albgri[ib] = 0.0;
            albd[ib] = 0.0;
            albi[ib] = 0.0;
            fabd[ib] = 0.0;
            fabi[ib] = 0.0;
            ftdd[ib] = 0.0;
            ftid[ib] = 0.0;
            ftii[ib] = 0.0;
            albsnd[ib] = 0.0;
            albsni[ib] = 0.0;
            frevd[ib] = 0.0;
            frevi[ib] = 0.0;
            fregd[ib] = 0.0;
            fregi[ib] = 0.0;
            ftdi[ib] = 0.0;
        }
        fsun = 0.0;

        // --------------------------------------------------------------------------------------------------
        // ALBEDO CALCULATION (inline to make it device-callable)
        // --------------------------------------------------------------------------------------------------

        // snow age (allow nighttime aging)
        snow_age_noahmp(dt_loc, tg, sneqvo, sneqv, d_tau0, d_grain_growth, d_extra_growth,
                        d_dirt_soot, d_swemx, tauss_arr(i,j,0), fage);

        if (cosz > 0.0) {

            // weight reflectance/transmittance by LAI and SAI
            vai = elai + esai;
            wl  = elai / std::max(vai,MPE);
            ws  = esai / std::max(vai,MPE);

            // Get parameters for this vegetation type
            for (int ib = 0; ib < NBAND; ib++) {
                amrex::Real rhol_val = (ib == 0) ? d_rhol_vis[veg_idx] : d_rhol_nir[veg_idx];
                amrex::Real rhos_val = (ib == 0) ? d_rhos_vis[veg_idx] : d_rhos_nir[veg_idx];
                amrex::Real taul_val = (ib == 0) ? d_taul_vis[veg_idx] : d_taul_nir[veg_idx];
                amrex::Real taus_val = (ib == 0) ? d_taus_vis[veg_idx] : d_taus_nir[veg_idx];
                rho[ib] = std::max(rhol_val*wl + rhos_val*ws, MPE);
                tau[ib] = std::max(taul_val*wl + taus_val*ws, MPE);
            }

            // snow albedos: only if COSZ > 0 and FSNO > 0
            if(d_opt_alb == 1) {
                snowalb_bats_noahmp(NBAND, fsno, cosz, fage,
                                   d_bats_cosz, d_bats_vis_new, d_bats_nir_new,
                                   d_bats_vis_age, d_bats_nir_age,
                                   d_bats_vis_dir, d_bats_nir_dir,
                                   albsnd, albsni);
            }

            // ground surface albedo
            // Get soil albedo parameters (using soil color index = 4 as default for now)
            int soil_color_idx = 3;  // index into albedo arrays (0-based for index 4 soil color)
            amrex::Real albsat[2], albdry[2], alblak[2];
            albsat[0] = d_albsat_vis[soil_color_idx];
            albsat[1] = d_albsat_nir[soil_color_idx];
            albdry[0] = d_albdry_vis[soil_color_idx];
            albdry[1] = d_albdry_nir[soil_color_idx];
            alblak[0] = d_alblak[0];
            alblak[1] = d_alblak[1];

            groundalb_noahmp(nsoil, NBAND, ice, ist, fsno, smc, albsnd, albsni, cosz, tg,
                            albsat, albdry, alblak, albgrd, albgri);

            // loop over NBAND wavebands to calculate surface albedos and solar
            // fluxes for unit incoming direct (IC=0) and diffuse flux (IC=1)

            for (int ib = 0; ib < NBAND; ib++) {
                // direct (IC=0)
                twostream_noahmp(ib, 0, vegtyp, cosz, vai, fwet, tv, albgrd, albgri,
                                rho, tau, fveg, ist,
                                Khai_L_arr(i,j,0), d_omegas[ib], d_betads, d_betais,
                                d_opt_rad, d_rc[veg_idx], ztop_arr(i,j,0), d_hvb[veg_idx], d_den[veg_idx],
                                fabd, albd, ftdd, ftid, gdir, frevd, fregd, bgap, wgap,
                                xl_diag_val, chil_diag_val, phi1_diag_val, phi2_diag_val);

                // diffuse (IC=1)
                twostream_noahmp(ib, 1, vegtyp, cosz, vai, fwet, tv, albgrd, albgri,
                                rho, tau, fveg, ist,
                                Khai_L_arr(i,j,0), d_omegas[ib], d_betads, d_betais,
                                d_opt_rad, d_rc[veg_idx], ztop_arr(i,j,0), d_hvb[veg_idx], d_den[veg_idx],
                                fabi, albi, ftdi, ftii, gdir, frevi, fregi, bgap, wgap,
                                xl_diag_val, chil_diag_val, phi1_diag_val, phi2_diag_val);
            }
            // sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01.
            ext = gdir/cosz * std::sqrt(1.0-rho[0]-tau[0]);
            fsun = (1.0-std::exp(-ext*vai)) / std::max(ext*vai,MPE);
            ext = fsun;

            if (ext < 0.01) {
                wl = 0.0;
            } else {
                wl = ext;
            }
            fsun = wl;

        } // cosz>0

        // --------------------------------------------------------------------------------------------------
        // SURRAD CALCULATION
        // --------------------------------------------------------------------------------------------------

        fsha = 1.0-fsun;
        laisun = elai*fsun;
        laisha = elai*fsha;
        vai = elai+ esai;

        surrad_noahmp(MPE, fsun, fsha, elai, vai, laisun, laisha, solad, solai,
                     fabd, fabi, ftdd, ftid, ftii, albgrd, albgri, albd, albi,
                     frevd, frevi, fregd, fregi,
                     parsun, parsha, sav, sag, fsa, fsr, fsrv, fsrg);

        // --------------------------------------------------------------------------------------------------
        // LONGWAVE RADIATION CALCULATION (NOAHMP formulation)
        // --------------------------------------------------------------------------------------------------
        // NOAHMP longwave scheme includes multiple reflections between canopy and ground
        // References: NOAHMP module_sf_noahmplsm.F lines 2137-2144 (emissivity)
        //                                               lines 3874-3875, 4041 (canopy LW)
        //                                               lines 4086-4087, 4104 (ground LW)
        // --------------------------------------------------------------------------------------------------

        // Get incoming longwave radiation from atmosphere (W/m²)
        // lwdn is positive downward (standard atmospheric convention)
        amrex::Real lwdn = lwref_arr(i, j, 0);

        // Calculate emissivities
        // Vegetation emissivity (NOAHMP line 2139)
        // EMV depends on total vegetation area index (LAI + SAI)
        // EMV = 1 for dense canopy, EMV → 0 for sparse/no vegetation
        amrex::Real emv = 1.0 - std::exp(-(elai + esai));

        // Ground emissivity (NOAHMP lines 2140-2144)
        // EMG is weighted by snow cover fraction (fsno)
        // Note: ice = 1 if IST = 2 (lake/ice), for now assume ice = 0 (soil)
        // Note: fsno and ist are already declared earlier (lines 5159, 5174)
        amrex::Real emg;
        if (ist == 1) {
            emg = d_eg_soil * (1.0 - fsno) + d_snow_emis * fsno;
        } else {
            emg = d_eg_lake * (1.0 - fsno) + d_snow_emis * fsno;
        }

        // --------------------------------------------------------------------------------------------------
        // Canopy net longwave radiation (IRC) - NOAHMP lines 3874-3875, 4041
        // --------------------------------------------------------------------------------------------------
        // NOAHMP SIGN CONVENTION: IRC is NET UPWARD longwave flux (positive = energy loss)
        //   IRC = emitted_upward - absorbed_from_atmosphere - absorbed_from_ground
        // NOAHMP energy balance: SAV - IRC - SHC - EVC = 0
        //
        // Physical interpretation of AIR_C term:
        //   -EMV*LWDN: canopy absorbs atmospheric LW (gain, hence negative)
        //   -EMV*(1-EMV)*(1-EMG)*LWDN: multiple reflection of atmospheric LW
        //                               (transmitted through canopy, reflected by ground, absorbed by canopy)
        //   -EMV*EMG*SB*TG^4: canopy absorbs ground emission (gain, hence negative)
        //
        // Physical interpretation of CIR_C*TV^4 term:
        //   2*EMV*SB*TV^4: canopy emits in both directions (up and down)
        //   -EMV^2*(1-EMG)*SB*TV^4: multiple reflection of canopy downward emission
        //                            (reflected by ground and re-absorbed by canopy)
        // --------------------------------------------------------------------------------------------------
        amrex::Real air_c = -emv * (1.0 + (1.0 - emv) * (1.0 - emg)) * lwdn
                          - emv * emg * SB * std::pow(tg, 4);
        amrex::Real cir_c = (2.0 - emv * (1.0 - emg)) * emv * SB;
        amrex::Real irc = air_c + cir_c * std::pow(tv, 4);

        // --------------------------------------------------------------------------------------------------
        // Ground net longwave radiation (IRG) - NOAHMP lines 4086-4087, 4104
        // --------------------------------------------------------------------------------------------------
        // NOAHMP SIGN CONVENTION: IRG is NET UPWARD longwave flux (positive = energy loss)
        //   IRG = emitted_upward - absorbed_from_atmosphere - absorbed_from_canopy
        // NOAHMP energy balance: SAG - IRG - SHG - EVG = 0
        //
        // Physical interpretation of AIR_G term:
        //   -EMG*(1-EMV)*LWDN: ground absorbs transmitted atmospheric LW (gain, hence negative)
        //   -EMG*EMV*SB*TV^4: ground absorbs canopy downward emission (gain, hence negative)
        //
        // Physical interpretation of CIR_G*TG^4 term:
        //   EMG*SB*TG^4: ground emits upward
        // --------------------------------------------------------------------------------------------------
        amrex::Real air_g = -emg * (1.0 - emv) * lwdn - emg * emv * SB * std::pow(tv, 4);
        amrex::Real cir_g = emg * SB;
        amrex::Real irg = cir_g * std::pow(tg, 4) + air_g;

        amrex::Real tir1 = emv * SB * std::pow(tv, 4);
        amrex::Real tir2 = emg * SB * std::pow(tg, 4);
        amrex::Real lwdn1 = lwdn;
        amrex::Real lwdn2 = (1.0 - emv) * lwdn + tir1;
        amrex::Real lwup2 = tir2 + (1.0 - emg) * lwdn2;
        amrex::Real lwup1 = lwdn + irc + irg;

        // --------------------------------------------------------------------------------------------------
        // Store shortwave and longwave radiation to net_rad arrays
        // --------------------------------------------------------------------------------------------------
        // SLM SIGN CONVENTION: net_rad is NET ABSORBED radiation (positive = energy gain)
        //   net_rad = absorbed_radiation (all sources) - emitted_radiation
        // SLM energy balance: net_rad - SH - LH = heat storage change
        //
        // SIGN CONVERSION: IRC/IRG (net upward) → net_lw (net absorbed)
        //   net_absorbed = -net_upward
        //   When IRC > 0 (net upward flux, energy loss), net_lw1 < 0 (energy loss)
        //   When IRC < 0 (net downward flux, energy gain), net_lw1 > 0 (energy gain)
        // --------------------------------------------------------------------------------------------------

        // Shortwave absorbed (W/m²)
        net_rad_arr(i, j, 0, SLM_NetRad::net_sw1) = sav;  // canopy absorbed solar
        net_rad_arr(i, j, 0, SLM_NetRad::net_sw2) = sag;  // ground absorbed solar

        // Longwave absorbed (W/m²)
        // Sign flip converts NOAHMP convention (net upward) to SLM convention (net absorbed)
        net_rad_arr(i, j, 0, SLM_NetRad::tir1) = tir1;
        net_rad_arr(i, j, 0, SLM_NetRad::tir2) = tir2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup1) = lwup1;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup2) = lwup2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwdn1) = lwdn1;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwdn2) = lwdn2;
        net_rad_arr(i, j, 0, SLM_NetRad::net_lw1) = -irc;  // canopy absorbed longwave
        net_rad_arr(i, j, 0, SLM_NetRad::net_lw2) = -irg;  // ground absorbed longwave

        // Total net radiation (W/m²)
        // net_rad = shortwave_absorbed + longwave_absorbed
        //         = sav + (-irc) = sav - irc
        net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) = sav - irc;  // canopy: SW + LW
        net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) = sag - irg;  // ground: SW + LW
        
        // Downwelling shortwave for stomatal resistance calculation
        // Total incoming SW = direct + diffuse for both visible and NIR
        net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1) = solad[0] + solad[1] + solai[0] + solai[1];
        net_rad_arr(i, j, 0, SLM_NetRad::net_swup1) = (solad[0] + solad[1] + solai[0] + solai[1]) - sav;  // incoming - absorbed = reflected

        // Downwelling SW reaching ground (transmitted through canopy)
        // This could be ftdd*solad + ftid*solad + ftii*solai, but simpler approximation:
        net_rad_arr(i, j, 0, SLM_NetRad::net_swdn2) = (solad[0] + solai[0]) * (ftdd[0] + ftid[0] + ftii[0]) +
                                                      (solad[1] + solai[1]) * (ftdd[1] + ftid[1] + ftii[1]);
        net_rad_arr(i, j, 0, SLM_NetRad::net_swup2) = net_rad_arr(i, j, 0, SLM_NetRad::net_swdn2) - sag;
        // --------------------------------------------------------------------------------------------------
        // Store albedo outputs
        // --------------------------------------------------------------------------------------------------

        alb_vis_sfc_arr(i, j, 0) = albd[0];
        alb_nir_sfc_arr(i, j, 0) = albd[1];
        alb_vis_sfc_diff_arr(i, j, 0) = albi[0];
        alb_nir_sfc_diff_arr(i, j, 0) = albi[1];

        alb_vis_sfc_arr(i, j, d_khi_lsm) = albd[0];
        alb_nir_sfc_arr(i, j, d_khi_lsm) = albd[1];
        alb_vis_sfc_diff_arr(i, j, d_khi_lsm) = albi[0];
        alb_nir_sfc_diff_arr(i, j, d_khi_lsm) = albi[1];

        // Surface emissivity for the vegetated canopy-ground system, including
        // the Noah-MP multiple-reflection term.
        amrex::Real emiss_sfc = emv + emg * (1.0 - emv)
                              + emv * (1.0 - emv) * (1.0 - emg);
        emis_sfc_arr(i, j, 0) = emiss_sfc;
        emis_sfc_arr(i, j, d_khi_lsm) = emis_sfc_arr(i, j, 0);

        t_skin_arr(i, j, 0) = std::pow(
            std::max(MPE, (lwup1 - (1.0 - emiss_sfc) * lwdn) / (emiss_sfc * SB)), 0.25);

    }); // End ParallelFor
}
