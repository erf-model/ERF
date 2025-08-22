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
           const Real& dt,
           std::unique_ptr<amrex::MultiFab>& z_phys_cc_in)
{
    m_dt = dt;
    m_geom = geom;
    z_phys_cc = z_phys_cc_in.get();

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

    LsmVarMap.resize(m_lsm_size);
    LsmVarMap = {
      LsmVar_SLM::tsurf,         LsmVar_SLM::ustar,        LsmVar_SLM::tstar,       LsmVar_SLM::qstar,
      LsmVar_SLM::tv,            LsmVar_SLM::mv,           LsmVar_SLM::soilt,       LsmVar_SLM::soilw,
      LsmVar_SLM::sand,          LsmVar_SLM::clay,         LsmVar_SLM::s_depth,
      LsmVar_SLM::flbu,          LsmVar_SLM::flbv,         LsmVar_SLM::flbq,
      LsmVar_SLM::flbt,          LsmVar_SLM::prsfc,        LsmVar_SLM::precipref,
      LsmVar_SLM::swdsvisxyref,  LsmVar_SLM::swdsnirxyref, LsmVar_SLM::swdsvisdxyref,
      LsmVar_SLM::swdsnirdxyref, LsmVar_SLM::lwref,        LsmVar_SLM::tref,
      LsmVar_SLM::uref,          LsmVar_SLM::vref,         LsmVar_SLM::dref,
      LsmVar_SLM::qref,          LsmVar_SLM::pref,         LsmVar_SLM::node_z,
      LsmVar_SLM::soilt_nudge,   LsmVar_SLM::soilw_nudge,  LsmVar_SLM::lai,
      LsmVar_SLM::vegtype,       LsmVar_SLM::soiltype};

    LsmVarName.resize(m_lsm_size);
    LsmVarName = {"tsurf",         "ustar",          "tstar",
                  "qstar",         "tveg",           "mv",
                  "tsoil",         "wsoil",          "sand",
                  "clay",          "soil_thickness", "surface_u",
                  "surface_v",     "surface_vapor",  "surface_heat",
                  "precip_soil",   "ref_precip",     "SW_dw_dir_vis",
                  "SW_dw_dir_nir", "SW_dw_dif_vis",  "SW_dw_dif_nir",
                  "LW_dw",         "ref_t",          "ref_u",
                  "ref_v",         "ref_d",          "ref_q",
                  "ref_p",         "node_z",         "soilt_nudge",
                  "soilw_nudge",   "lai",            "vegtype",
                  "soiltype"};

    AMREX_ALWAYS_ASSERT(LsmVarMap.size() == LsmVarName.size());
    AMREX_ALWAYS_ASSERT(LsmVarMap.size() == m_lsm_size);

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
    m_lsm_geom.define( ba_lsm.minimalBox(), lsm_rb, m_geom.Coord(), m_geom.isPeriodic() );

    // Create the data and fluxes
    for (auto ivar = 0; ivar < LsmVar_SLM::NumVars; ++ivar) {
        // State vars are CC
        Real theta_0 = m_theta_dir;
        lsm_fab_vars[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_vars[ivar]->setVal(0.0);

        // Fluxes are nodal in z
        //lsm_fab_flux[ivar] = std::make_shared<MultiFab>(convert(ba_lsm, IntVect(0,0,1)), dm, 1, IntVect(0,0,0));
        //lsm_fab_flux[ivar]->setVal(0.);
    }

    // packed temporary 1D arrays in soil water and soil temperature
    soilt_vars.define(ba_lsm, dm, SLM_DST::NumVars, 0);
    soilw_vars.define(ba_lsm, dm, SLM_DSW::NumVars, 0);
    soilt_vars.setVal(0.0);
    soilw_vars.setVal(0.0);

    // Create local 2D data
    BoxList bl_lsm_2d = ba_lsm.boxList();
    for (auto& b : bl_lsm_2d) {
        b.setRange(2, 0, 1);
    }
    ba_lsm_2d = BoxArray(std::move(bl_lsm_2d));
    IntVect ng_2d(0, 0, 0);

    // TODO: Placeholder landmask array - fix!
    landmask.define(ba_lsm_2d, dm, 1, ng_2d);
    landmask.setVal(1);

    landtype.define(ba_lsm_2d, dm, 1, ng_2d);
    LAI.define(ba_lsm_2d, dm, 1, ng_2d);
    sstxy.define(ba_lsm_2d, dm, 1, ng_2d);

    t_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    mw.define(ba_lsm_2d, dm, 1, ng_2d);
    mws.define(ba_lsm_2d, dm, 1, ng_2d);
    t_skin.define(ba_lsm_2d, dm, 1, ng_2d);
    t_cas.define(ba_lsm_2d, dm, 1, ng_2d);
    q_cas.define(ba_lsm_2d, dm, 1, ng_2d);

    vegetype.define(ba_lsm_2d, dm, 1, ng_2d);
    vege_YES.define(ba_lsm_2d, dm, 1, ng_2d);
    cp_vege.define(ba_lsm_2d, dm, 1, ng_2d);
    z0_sfc.define(ba_lsm_2d, dm, 1, ng_2d);
    Khai_L.define(ba_lsm_2d, dm, 1, ng_2d);
    phi_1.define(ba_lsm_2d, dm, 1, ng_2d);
    phi_2.define(ba_lsm_2d, dm, 1, ng_2d);
    IR_emis_vege.define(ba_lsm_2d, dm, 1, ng_2d);
    IR_emis_soil.define(ba_lsm_2d, dm, 1, ng_2d);
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

    mw_inc.define(ba_lsm_2d, dm, 1, ng_2d);

    evapo_dry.define(ba_lsm_2d, dm, 1, ng_2d);

    shf_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    shf_soil.define(ba_lsm_2d, dm, 1, ng_2d);
    shf_air.define(ba_lsm_2d, dm, 1, ng_2d);
    lhf_canop.define(ba_lsm_2d, dm, 1, ng_2d);
    lhf_soil.define(ba_lsm_2d, dm, 1, ng_2d);
    lhf_air.define(ba_lsm_2d, dm, 1, ng_2d);

    albedovis_v.define(ba_lsm_2d, dm, 1, ng_2d);
    albedonir_v.define(ba_lsm_2d, dm, 1, ng_2d);
    albedovis_s.define(ba_lsm_2d, dm, 1, ng_2d);
    albedonir_s.define(ba_lsm_2d, dm, 1, ng_2d);

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

    r_soil.setVal(0.0);
    lhf_air.setVal(0.0);
    lhf_canop.setVal(0.0);
    shf_air.setVal(0.0);
    shf_canop.setVal(0.0);
    t_canop.setVal(0.0);
    t_cas.setVal(0.0);
    t_skin.setVal(0.0);
    q_cas.setVal(0.0);
    wet_canop.setVal(0.0);

    // Initialize 1D arrays
    soilw_inc.resize({klo_lsm},  {khi_lsm});
    alpha.resize({klo_lsm},  {khi_lsm});
    beta.resize({klo_lsm},  {khi_lsm});

    r_a.setVal(0.0);
    r_b.setVal(0.0);
    r_c.setVal(0.0);
    r_d.setVal(0.0);
    lsm_fab_vars[LsmVar_SLM::ustar]->setVal(0.1);
    lsm_fab_vars[LsmVar_SLM::tstar]->setVal(0.0);
    lsm_fab_vars[LsmVar_SLM::qstar]->setVal(0.0);

    // Initialize SLM from inputs if specified
    init_from_file();

    landtype.setVal(landtype0);
    LAI.setVal(LAI0);
    mws_mx.setVal(mws_mx0);

    net_rad.setVal(0.0);

    if (!use_wrfinput) {
        // Initialize here if not using wrfinput, otherwise it is done at first time step
        slm_init();
    }

	//Following Noah-MP, zref is modified to become ztop (canopy topheight) + dz0(center height of the atmosphere's lowest grid)
    //First, read in SLM_use_inputs, if true, zref is set from the inputs parameter slm.zref
	// if false, SLM is coupled to ERF atmosphere, therefore, zref is computed  
	pp.query("SLM_use_inputs",set_from_file);

    if (!set_from_file) {
        ParmParse pp_erf("erf");
        pp_erf.query("use_terrain", use_terrain);

        Real zlo      = m_geom.ProbLo(2);
        Real dz       = m_geom.CellSize(2);
        for ( MFIter mfi(cons_in,TileNoZ()); mfi.isValid(); ++mfi) {
            const Box& xybx      = mfi.growntilebox(0);
            const Array4<const Real>& z_cc_arr = (use_terrain) ? z_phys_cc->const_array(mfi) : Array4<Real>{};
            auto ztop_arr = ztop.array(mfi);
			auto zrefxy_arr = zrefxy.array(mfi);

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
                Real zcc = (z_cc_arr) ? z_cc_arr(i,j,0) : zlo + 0.5*dz;
			    zrefxy_arr(i,j,0) = ztop_arr(i,j,0) + zcc;
		    });	
	    }		
    } else {
	    zrefxy.setVal(zref);  // set zrefxy with the value read from inputs file
    }

    //slm_to_rad_vars = {Lsm_Data_Ptr(LsmVar_SLM::tsurf), &albedovis_s, &albedovis_v, &albedonir_s, &albedonir_v, &IR_emis_vege, &net_rad};
    slm_to_rad_vars = {Lsm_Data_Ptr(LsmVar_SLM::tsurf), &IR_emis_vege, &albedovis_v, &albedonir_v, &albedovis_v, &albedonir_v};
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

    get_layer_prop("clay0", m_nz_lsm, clay0);
    get_layer_prop("sand0", m_nz_lsm, sand0);
    if (!use_wrfinput) {
        get_layer_prop("sw0", m_nz_lsm, sw0);
        get_layer_prop("st0", m_nz_lsm, st0);
    }
    if (dosoiltnudging || dosoilwnudging) {
        get_layer_prop("relax_hgt", m_nz_lsm, relax_hgt);
    } else {
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

        auto IR_emis_vege_arr = IR_emis_vege.array(mfi);
        auto phi_1_arr = phi_1.array(mfi);
        auto phi_2_arr = phi_2.array(mfi);
        auto precip_extinc_arr = precip_extinc.array(mfi);
        auto mw_mx_arr = mw_mx.array(mfi);
        auto LAI_arr = LAI.array(mfi);
        auto Khai_L_arr = Khai_L.array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (vegetype_arr(i, j, 0) == 0) {
                vege_YES_arr(i, j, 0) = 0.0;
            } else {
                // set minimum LAI for vegetated land
                LAI_arr(i, j, 0) = std::max(LAI_arr(i, j, 0), 0.001);
            }

            IR_emis_vege_arr(i, j, 0) = 0.97 * (1.0 - std::exp(-1.0 * LAI_arr(i, j, 0)));
            phi_1_arr(i, j, 0) = 0.5 - 0.633 * Khai_L_arr(i, j, 0) - 0.33 * (std::pow(Khai_L_arr(i, j, 0), 2));
            phi_2_arr(i, j, 0) = 0.877 * (1.0 - 2.0 * phi_1_arr(i, j, 0));
            precip_extinc_arr(i, j, 0) = phi_1_arr(i, j, 0) + phi_2_arr(i, j, 0);
            mw_mx_arr(i, j, 0) = 0.1 * LAI_arr(i, j, 0);
        });
    }

    // Initialize soil parameters
    init_soil_tw();

    IR_emis_soil.setVal(0.98);

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

                // soil wetness at field capacity
                w_s_FC_arr(i, j, k) = theta_FC_arr(i, j, k) / poro_soil_arr(i, j, k);

                // soil wetness at wilting point
                w_s_WP_arr(i, j, k) = theta_WP_arr(i, j, k) / poro_soil_arr(i, j, k);
            }
        });
    }

    // Calculate fraction of root in each soil layer
    vege_root_init();
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
        auto vegetype_arr = vegetype.array(mfi);

        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            switch (landtype_arr(i, j, 0))
            {
                case 0: // water
                    // TODO: should this landmask be ERF's lmask_lev variable?
                    landmask_arr(i, j, 0) = 0;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                case 1: // evergreen needleaf forest
                    albedovis_v_arr(i, j, 0) = 0.04;
                    albedonir_v_arr(i, j, 0) = 0.20;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 1.09;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 100.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 2: // evergreen broadleaf forest
                    albedovis_v_arr(i, j, 0) = 0.04;
                    albedonir_v_arr(i, j, 0) = 0.20;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 2.65;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 150.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 100.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 3: // deciduous needleaf forest
                    albedovis_v_arr(i, j, 0) = 0.05;
                    albedonir_v_arr(i, j, 0) = 0.23;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.85;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 100.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 4: // deciduous broadleaf forest
                    albedovis_v_arr(i, j, 0) = 0.07;
                    albedonir_v_arr(i, j, 0) = 0.24;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.8;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 100.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 5: // mixed forest
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.24;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 20.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.8;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 100.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 6: // closed shrublands
                    albedovis_v_arr(i, j, 0) = 0.07;
                    albedonir_v_arr(i, j, 0) = 0.26;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 1.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.03;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 30.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 7: // open shrublands
                    albedovis_v_arr(i, j, 0) = 0.14;
                    albedonir_v_arr(i, j, 0) = 0.32;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 1.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.03;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 50.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 30.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 8: // woody savannas
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.21;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 5.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.86;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 50.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 9: // savannas
                    albedovis_v_arr(i, j, 0) = 0.07;
                    albedonir_v_arr(i, j, 0) = 0.26;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 5.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.86;
                    Khai_L_arr(i, j, 0) = 0.25;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 30.;
                    hs_rc_arr(i, j, 0) = 54.53;
                    BAI_arr(i, j, 0) = 50.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 10: // grasslands
                    albedovis_v_arr(i, j, 0) = 0.07;
                    albedonir_v_arr(i, j, 0) = 0.25;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.08;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 170.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 39.18;
                    BAI_arr(i, j, 0) = 30.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 11: // permanent wetlands
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.18;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.04;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 0.0;
                    Rgl_arr(i, j, 0) = 0.0;
                    hs_rc_arr(i, j, 0) = 0.0;
                    BAI_arr(i, j, 0) = 0.0;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                case 12: // croplands
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.24;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.075;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 39.18;
                    BAI_arr(i, j, 0) = 50.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 13: // urban
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.22;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 10.;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 1.0;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 0.0;
                    Rgl_arr(i, j, 0) = 0.0;
                    hs_rc_arr(i, j, 0) = 0.0;
                    BAI_arr(i, j, 0) = 0.0;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                case 14: // croplands/natural mozaics
                    albedovis_v_arr(i, j, 0) = 0.06;
                    albedonir_v_arr(i, j, 0) = 0.22;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
                    ztop_arr(i, j, 0) = 0.5;
                    disp_hgt_arr(i, j, 0) = 0.68 * ztop_arr(i, j, 0);
                    z0_sfc_arr(i, j, 0) = 0.07;
                    Khai_L_arr(i, j, 0) = -0.3;
                    rootL_arr(i, j, 0) = 1.5;
                    root_a_arr(i, j, 0) = 5.558;
                    root_b_arr(i, j, 0) = 2.614;
                    Rc_min_arr(i, j, 0) = 100.;
                    Rgl_arr(i, j, 0) = 100.;
                    hs_rc_arr(i, j, 0) = 39.18;
                    BAI_arr(i, j, 0) = 50.;
                    vegetype_arr(i, j, 0) = 1;
                    break;
                case 15: // snow/ice
                    albedovis_v_arr(i, j, 0) = 0.0; // actually depends on moisture content
                    albedonir_v_arr(i, j, 0) = 0.0;
                    albedovis_s_arr(i, j, 0) = 0.95;
                    albedonir_s_arr(i, j, 0) = 0.65;
                    ztop_arr(i, j, 0) = 0.0;
                    disp_hgt_arr(i, j, 0) = 0.0;
                    z0_sfc_arr(i, j, 0) = 0.01;
                    Khai_L_arr(i, j, 0) = 0.0;
                    rootL_arr(i, j, 0) = 0.0;
                    root_a_arr(i, j, 0) = 0.0;
                    root_b_arr(i, j, 0) = 0.0;
                    Rc_min_arr(i, j, 0) = 0.0;
                    Rgl_arr(i, j, 0) = 0.0;
                    hs_rc_arr(i, j, 0) = 0.0;
                    BAI_arr(i, j, 0) = 0.0;
                    vegetype_arr(i, j, 0) = 0;
                    break;
                case 16: // baresoil
                    albedovis_v_arr(i, j, 0) = 0.0; // actually depends on moisture content
                    albedonir_v_arr(i, j, 0) = 0.0;
                    albedovis_s_arr(i, j, 0) = 0.19;
                    albedonir_s_arr(i, j, 0) = 0.38;
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
        auto soilt_obs_arr = lsm_fab_vars[LsmVar_SLM::soilt_obs]->array(mfi);
        auto soilw_obs_arr = lsm_fab_vars[LsmVar_SLM::soilw_obs]->array(mfi);

        auto sand_arr = lsm_fab_vars[LsmVar_SLM::sand]->array(mfi);
        auto clay_arr = lsm_fab_vars[LsmVar_SLM::clay]->array(mfi);

        auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->array(mfi);
        auto soil_relax_hgt_arr = lsm_fab_vars[LsmVar_SLM::soil_relax_hgt]->array(mfi);

        auto sstxy_arr = sstxy.array(mfi);

        auto soiltype_arr = lsm_fab_vars[LsmVar_SLM::soiltype]->array(mfi);

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
                        default:
                            sand_arr(i,j,k) = -1.0;
                            clay_arr(i,j,k) = -1.0;
                            break;
                    }

                }

                // initialize nudging profiles for soil based on the initial soilt and soilw
                soilt_obs_arr(i, j, k) = soilt_arr(i, j, k);
                soilw_obs_arr(i, j, k) = soilw_arr(i, j, k);
                soil_relax_hgt_arr(i, j, k) = d_relax[(k*-1)+d_khi_lsm];

                // TODO: nrestart conditional here
                // TODO: sstxy - should be ERF surface temp array?
                sstxy_arr(i, j, 0) = soilt_arr(i, j, d_khi_lsm) - d_t00;
            } else {
                sstxy_arr(i, j, 0) = d_tabs_s - d_t00;
            }
        });
    }

    // TODO: is it worth storing s_depth_mm, or should we always convert?
    // s_depth_mm[] = s_depth[] * 1000.0

    // Calculate node_z (depth of the center of each soil layer)
    // TODO: this can probably be done natively via AMReX?
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

        auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
        auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
        auto tref_arr  = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
        auto tsurf_arr  = lsm_fab_vars[LsmVar_SLM::tsurf]->const_array(mfi);
        auto qref_arr  = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
        auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);

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
                t_canop_arr(i, j, 0) = 0.5*(tref_arr(i, j, 0) + tsurf_arr(i, j, 0));
                t_cas_arr(i, j, 0) = 0.5*(tref_arr(i, j, 0) + tsurf_arr(i, j, 0));

                // specific humidity in canopy air space
                q_cas_arr(i, j, 0) = 0.5*(qref_arr(i, j, 0) + q_gr);

                //amrex::Print() << " i = " << i << " j = " << j << " k = 0 : t_canop = " << t_canop_arr(i, j, 0) << " t_cas = " << t_cas_arr(i, j, 0) << " tsurf = " << tsurf_arr(i, j, 0) << " qcas = " << q_cas_arr(i, j, 0) << " qref = " << qref_arr(i, j, 0) << std::endl;

            }
        });
    }
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

        auto mw_arr = mw.const_array(mfi);
        auto mws_arr = mws.const_array(mfi);
        auto mw_mx_arr = mw_mx.const_array(mfi);

        auto mw_inc_arr = mw_inc.const_array(mfi);

        auto t_cas_arr = t_cas.array(mfi);
        auto q_cas_arr = q_cas.array(mfi);

        auto t_canop_arr = t_canop.array(mfi);
        auto t_skin_arr = t_skin.const_array(mfi);

        auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
        auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);

        auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
        auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);

        auto ustar_arr = lsm_fab_vars[LsmVar_SLM::ustar]->const_array(mfi);
        auto tstar_arr = lsm_fab_vars[LsmVar_SLM::tstar]->array(mfi);
        auto qstar_arr = lsm_fab_vars[LsmVar_SLM::qstar]->array(mfi);

        auto BAI_arr = BAI.const_array(mfi);
        auto ztop_arr = ztop.const_array(mfi);

        auto shf_canop_arr = shf_canop.array(mfi);
        auto shf_soil_arr = shf_soil.array(mfi);
        auto shf_air_arr = shf_air.array(mfi);

        auto lhf_canop_arr = lhf_canop.array(mfi);
        auto lhf_soil_arr = lhf_soil.array(mfi);
        auto lhf_air_arr = lhf_air.array(mfi);

        auto vegetype_arr = vegetype.const_array(mfi);
        auto vege_YES_arr = vege_YES.const_array(mfi);

        auto cp_vege_arr = cp_vege.array(mfi);

        auto tref_arr  = lsm_fab_vars[LsmVar_SLM::tref]->const_array(mfi);
        auto ur_arr  = lsm_fab_vars[LsmVar_SLM::uref]->const_array(mfi);
        auto vr_arr  = lsm_fab_vars[LsmVar_SLM::vref]->const_array(mfi);
        auto dref_arr  = lsm_fab_vars[LsmVar_SLM::dref]->const_array(mfi);
        auto qref_arr  = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
        auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);
        auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

        auto tsurf_arr = lsm_fab_vars[LsmVar_SLM::tsurf]->array(mfi);

        auto r_a_arr = r_a.const_array(mfi);
        auto r_b_arr = r_b.const_array(mfi);
        auto r_c_arr = r_c.const_array(mfi);
        auto r_d_arr = r_d.const_array(mfi);
        auto r_soil_arr = r_soil.const_array(mfi);

        //auto flbu_arr  = lsm_fab_vars[LsmVar_SLM::flbu]->array(mfi);
        //auto flbv_arr  = lsm_fab_vars[LsmVar_SLM::flbv]->array(mfi);
        auto flbq_arr  = lsm_fab_vars[LsmVar_SLM::flbq]->array(mfi);
        auto flbt_arr  = lsm_fab_vars[LsmVar_SLM::flbt]->array(mfi);
        auto prsfc_arr  = lsm_fab_vars[LsmVar_SLM::prsfc]->const_array(mfi);

        auto net_rad_arr = net_rad.const_array(mfi);
        auto wet_canop_arr = wet_canop.const_array(mfi);

        // Calculate net radiation absorbed by canopy and soil surface
        radiative_fluxes(mfi);

        // Calculate turbulent transfer coeff between reference level and surface
        transfer_coeff(mfi);

        // Calculate aerodynamic resistances + stomatal resistance
        resistances(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            //amrex::Real precip = 0.0; // precipitation interception rate at canopy
            //amrex::Real precip_sfc = 0.0; // precipitation reaching to soil surface
            //amrex::Real drain = 0.0; // drainage rate from canopy
            //amrex::Real cnp_mw_drip = 0.0; // dripping from canopy water storage when the storage exceeds mw_mx - mm/s

            amrex::Real t_sfc;
            // SAM rhow[nz] = air density at vertical velocity levels, kg/m^3
            amrex::Real rhow = dref_arr(i, j, 0); // TODO: double check this

            if (landmask_arr(i, j, 0) == 1)
            {

                // Calculate net radiation absorbed by canopy and soil surface
                //radiative_fluxes(i, j);

                if (vegetype_arr(i, j, 0) != 0)
                {
                    // vegetated surfaces
                    t_sfc = t_cas_arr(i, j, 0);
                    //q_sfc = q_cas_arr(i, j, 0);
                }
                else
                {
                    // baresoil
                    t_sfc = soilt_arr(i, j, d_khi_lsm);
                    //q_sfc = q_gr;
                }

                // Calculate turbulent transfer coeff between reference level and surface
                //transfer_coeff(i, j);

                // Calculate surface momentum fluxes
                //taux_sfc = -1.0 * mom_trans_coef * vel_m * ur_arr(i, j, 0) * (pres0 * 100.0 / 287.0 / tref_arr(i, j, 0));
                //tauy_sfc = -1.0 * mom_trans_coef * vel_m * vr_arr(i, j, 0) * (pres0 * 100.0 / 287.0 / tref_arr(i, j, 0));

                // Calculate aerodynamic resistances + stomatal resistance
                //resistances(i, j);

                // Sensible heat fluxes
                shf_canop_arr(i, j, 0) = (t_canop_arr(i, j, 0) - t_sfc) * rhow * cp / r_b_arr(i, j, 0) * vege_YES_arr(i, j, 0);
                shf_soil_arr(i, j, 0) = (soilt_arr(i, j, d_khi_lsm) - t_sfc) * rhow * cp / r_d_arr(i, j, 0);
                shf_air_arr(i, j, 0) = (t_sfc - tref_arr(i, j, 0)) * rhow * cp / r_a_arr(i, j, 0);

                if (vegetype_arr(i, j, 0) == 0)
                {
                    shf_soil_arr(i, j, 0) = shf_air_arr(i, j, 0);
                }
                else
                {
                    shf_air_arr(i, j, 0) = shf_canop_arr(i, j, 0) + shf_soil_arr(i, j, 0);
                }

                // Calculate temperature scale for z0hsfc
                //tstar_arr(i, j, 0) = -1.0 * shf_air_arr(i, j, 0) / rhow / cp / ustar_arr(i, j, 0);
                tstar_arr(i, j, 0) = shf_air_arr(i, j, 0) / rhow / cp;
                tstar_arr(i, j, d_khi_lsm) = tstar_arr(i, j, 0);

                // Calculate latent heat fluxes
                //vapor_fluxes(i, j); // -- computes r_soil

                // Calculate soil moisture increment
                //soil_water(mfi);

                // Update vegetation moisture storage
                //mw_arr(i, j, 0) += mw_inc_arr(i, j, 0);

                // Output variables
                //flbu_arr(i, j, 0) = taux_sfc;
                //flbv_arr(i, j, 0) = tauy_sfc;
                //prsfc_arr(i, j, 0) = precip_sfc;
            }
        });

        // Calculate latent heat fluxes
        vapor_fluxes(mfi); // -- computes r_soil

        // Calculate soil moisture increment
        soil_water(mfi);

        // Calculate soil temperature increment
        //  amrex::Real grflux0 = -1.0 * (net_rad[1] - shf_soil_arr(i, j, 0) - lhf_soil_arr(i, j, 0));
        soil_temperature(mfi);

        ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            amrex::Real cp_vege_tot;
            amrex::Real q_gr;

            // SAM rhow[nz] = air density at vertical velocity levels, kg/m^3
            amrex::Real rhow = dref_arr(i, j, 0); // TODO: double check this

            if (landmask_arr(i, j, 0) == 1) {
                // Update vegetation temperature
                cp_vege_arr(i, j, 0) = (LAI_arr(i, j, 0)*0.001 + ztop_arr(i, j, 0)*BAI_arr(i, j, 0) / 43560.0) * 900.0 * 2800.0;
                cp_vege_tot = cp_vege_arr(i, j, 0) + mw_arr(i, j, 0) * 1.0e-3 * cp_water;

                amrex::Real t_canop_inc = dt / std::max(1.0e-3, cp_vege_tot)*(net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) - shf_canop_arr(i, j, 0) - lhf_canop_arr(i, j, 0)) * vege_YES_arr(i, j, 0);
                t_canop_arr(i, j, 0) += t_canop_inc;

                // Compute diagnostic variables at canopy air space level
                //   Calculate heat conductances - non-zero only for canopy land type
                amrex::Real cond_heat = 1.0 / r_a_arr(i, j, 0) + 1.0 / r_b_arr(i, j, 0) + 1.0 / r_d_arr(i, j, 0);
                amrex::Real cond_href = 1.0 / r_a_arr(i, j, 0) / cond_heat * vege_YES_arr(i, j, 0);
                amrex::Real cond_hcnp = 1.0 / r_b_arr(i, j, 0) / cond_heat * vege_YES_arr(i, j, 0);
                amrex::Real cond_hundercnp = 1.0 / r_d_arr(i, j, 0) / cond_heat * vege_YES_arr(i, j, 0);

                //  Calculate vapor conductances
                amrex::Real cond_vapor = 1.0 / r_a_arr(i, j, 0) + wet_canop_arr(i, j, 0) / (2.0 * r_b_arr(i, j, 0)) + (1.0 - wet_canop_arr(i, j, 0))/(2.0 * r_b_arr(i, j, 0) + r_c_arr(i, j, 0)) + 1.0 / (r_d_arr(i, j, 0) + r_soil_arr(i, j, 0) + r_litter);
                amrex::Real cond_vref = 1.0 / r_a_arr(i, j, 0) / cond_vapor * vege_YES_arr(i, j, 0);
                amrex::Real cond_vcnp = (wet_canop_arr(i, j, 0) / (2.0 * r_b_arr(i, j, 0)) / cond_vapor + (1.0 - wet_canop_arr(i, j, 0)) / (2.0 * r_b_arr(i, j, 0) + r_c_arr(i, j, 0)) / cond_vapor) * vege_YES_arr(i, j, 0);
                amrex::Real cond_vundercnp = 1.0 / (r_d_arr(i, j, 0) + r_soil_arr(i, j, 0) + r_litter) / cond_vapor * vege_YES_arr(i, j, 0);

                // Canopy air space
                if (vegetype_arr(i, j, 0) == 0)
                {
                    cond_hcnp = 1.0;
                    cond_vundercnp = 1.0;
                }

                t_cas_arr(i, j, 0) = tref_arr(i, j, 0)*cond_href + t_canop_arr(i, j, 0)*cond_hcnp + soilt_arr(i, j, d_khi_lsm)*cond_hundercnp;
                amrex::Real qsat_canop;
                if (t_canop_arr(i, j, 0) >= tfriz)
                {
                    erf_qsatw(t_canop_arr(i, j, 0), pref_arr(i, j, 0), qsat_canop);
                }
                else
                {
                    erf_qsati(t_canop_arr(i, j, 0), pref_arr(i, j, 0), qsat_canop);
                }

                if (soilt_arr(i, j, d_khi_lsm) >= tfriz)
                {
                    erf_qsatw(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
                    q_gr *= fh_calc(soilt_arr(i, j, d_khi_lsm), m_pot_sat_arr(i, j, d_khi_lsm), soilw_arr(i, j, d_khi_lsm), Bconst_arr(i, j, d_khi_lsm));
                }
                else
                {
                    erf_qsati(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
                }

                q_cas_arr(i, j, 0) = qref_arr(i, j, 0) * cond_vref + qsat_canop*cond_vcnp + q_gr*cond_vundercnp;

                // Output variables
                tsurf_arr(i, j, d_khi_lsm) = t_skin_arr(i, j, 0); // TODO: ts in SLM is input and output - check how this should be coupled back to ERF
                flbq_arr(i, j, d_khi_lsm) = lhf_air_arr(i, j, 0) / (lcond*rhow);
                flbt_arr(i, j, d_khi_lsm) = shf_air_arr(i, j, 0) / (Cp_d*rhow);

                //qstar_arr(i, j, 0) = -1.0 * lhf_air_arr(i, j, 0) / (lcond*rhow) / ustar_arr(i, j, 0);
                qstar_arr(i, j, 0) = lhf_air_arr(i, j, 0) / (lcond*rhow);
                qstar_arr(i, j, d_khi_lsm) = qstar_arr(i, j, 0);
                // Collect 2D stat variables
                // collect_2D_stat_vars(i, j);
            }
        });
    }

    lsm_fab_vars[LsmVar_SLM::tsurf]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::ustar]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::tstar]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::qstar]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::flbu]->FillBoundary(m_geom.periodicity());
    lsm_fab_vars[LsmVar_SLM::flbv]->FillBoundary(m_geom.periodicity());
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
    auto t_canop_arr      = t_canop.const_array(mfi);
    auto soilw_arr        = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto soilt_arr        = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);

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
            ka = phi_1_arr(i, j, 0) / coszrsxy_arr(i, j, 0) + phi_2_arr(i, j, 0);
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
            net_rad_arr(i, j, 0, SLM_NetRad::net_swup1) = net_rad_arr(i, j, 0, SLM_NetRad::net_rad1) - net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1);

            // net_rad(2) = net absorbed shortwave radiation by soil
            wetfactor = 1.0 - 0.5*soilw_arr(i, j, d_khi_lsm); // soil wetness factor: assume that wet soil is twice as dark
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsvisxyref_arr(i, j, 0)*(1.0 - albedovis_s_arr(i, j, 0)*wetfactor)*explai;
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsvisdxyref_arr(i, j, 0)*(1.0 - albedovis_s_arr(i, j, 0)*wetfactor)*explai0;
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsnirxyref_arr(i, j, 0)*(1.0 - albedonir_s_arr(i, j, 0)*wetfactor)*explai;
            net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) += swdsnirdxyref_arr(i, j, 0)*(1.0 - albedonir_s_arr(i, j, 0)*wetfactor)*explai0;

            net_rad_arr(i, j, 0, SLM_NetRad::net_swdn2) = net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1)*explai;
            net_rad_arr(i, j, 0, SLM_NetRad::net_swup2) = net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) - net_rad_arr(i, j, 0, SLM_NetRad::net_swdn2);

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
        net_rad_arr(i, j, 0, SLM_NetRad::tir2) = IR_emis_soil_arr(i, j, 0)*sigma*(std::pow(soilt_arr(i, j, d_khi_lsm), 4));

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
        //  Emitted LW from topsoil = emitted tir from topspoil + portion of incoming LW that is reflected back toward canopy
        //  IR_emis_soil is set to 1.0
        // ===================================================
        fup1 = (net_rad_arr(i, j, 0, SLM_NetRad::tir2) + (1.0 - IR_emis_soil_arr(i, j, 0))*fdn1);
        net_rad_arr(i, j, 0, SLM_NetRad::net_lwup2) = fup1;

        // ===================================================
        // Note: At this stage,
        //  net_rad(2) = net absorbed SW by soil surface
        //          + net absorbed LW by soil surface
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
    auto tstar_arr = lsm_fab_vars[LsmVar_SLM::tstar]->const_array(mfi);

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
    constexpr amrex::Real kk = 0.35;
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

        amrex::Real t_sfc, q_sfc, q_gr;
        // TODO: This code is repeated across several functions - store q_sfc,t_sfc for reuse?
        // determine input q_sfc and t_sfc
        if (vegetype_arr(i, j, 0) != 0)
        {
            // vegetated surfaces
            q_sfc = q_cas_arr(i, j, 0);
            t_sfc = t_cas_arr(i, j, 0);
        }
        else
        {
            // baresoil
            t_sfc = soilt_arr(i, j, d_khi_lsm);

            // Specific humidity at top soil
            if (soilt_arr(i, j, d_khi_lsm) > tfriz)
            {
                erf_qsatw(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
                if (mws_arr(i, j, 0) == 0.0)
                {
                    q_gr *= fh_calc(soilt_arr(i, j, d_khi_lsm), m_pot_sat_arr(i, j, d_khi_lsm), soilw_arr(i, j, d_khi_lsm), Bconst_arr(i, j, d_khi_lsm));
                }
            }
            else
            {
                erf_qsati(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
            }

            q_sfc = q_gr;
        }
        // Inputs:
        // ts = t_sfc
        // th = tref
        // qh = qr
        // qs = q_sfc
        // h = zref = height of ref level
        // z0 = z0_sfc = surface roughness length
        // disp = disp_hgt

        amrex::Real tsp = t_sfc * std::pow(1000.0/pref_arr(i, j, 0), rair / cp);
        amrex::Real thp = tref_arr(i, j, 0) * std::pow(1000.0 / pref_arr(i, j, 0), rair / cp);

        amrex::Real vel;
        // Add additional velocity depending on the stratification
        if ((thp - tsp) >= 0.0)
        {
            vel = sqrt(std::pow(ur_arr(i, j, 0), 2) + std::pow(vr_arr(i, j, 0), 2) + 0.1*0.1);
        }
        else
        {
            vel = sqrt(std::pow(ur_arr(i, j, 0), 2) + std::pow(vr_arr(i, j, 0), 2) + 1.0);
        }
        const Real d_zref = zref_arr(i,j,0);

        amrex::Real r = 9.81 / tsp * (thp * (1.0 + epsv * qr_arr(i, j, 0)) - tsp * (1.0 + epsv * q_sfc)) * (d_zref - disp_hgt_arr(i, j, 0)) / (vel*vel);
        r = std::max(-10.0, std::min(r, 0.19));
        // initial guess
        amrex::Real xsi, fm, fh, xsi1;
        amrex::Real xsim0, xsih0;
        amrex::Real mom_trans_coef, heat_trans_coef;
        amrex::Real z0h = z0_sfc_arr(i, j, 0) / (d_zref - disp_hgt_arr(i, j, 0));

        // make sure (h - disp) is not negative, otherwise z0dym,z0dyh become nan
        AMREX_ALWAYS_ASSERT(d_zref - disp_hgt_arr(i, j, 0) > 0.0);

        amrex::Real tstar_in = -1.0 * tstar_arr(i, j, 0) / ustar_arr(i, j, 0);
        amrex::Real zt0 = std::max(0.0001, (70.0*1.5e-5 / ustar_arr(i, j, 0)) * std::exp(-7.2*sqrt(ustar_arr(i, j, 0))*(std::pow(std::abs(tstar_in), 0.25))));

        amrex::Real zTh = zt0 / (d_zref - disp_hgt_arr(i, j, 0));
        amrex::Real zodym = log(1.0 / z0h);
        amrex::Real zodyh = log(1.0 / zTh);

        if (r > 0.0)
        {
            xsi = r * zodym / (1.0 - 5.0 * r);
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
        ustar_arr(i, j, 0) = std::max(0.2, ustar_arr(i, j, 0));

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
    auto vegetype_arr = vegetype.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto s_depth_arr = lsm_fab_vars[LsmVar_SLM::s_depth]->const_array(mfi);
    auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);

    auto LAI_arr = LAI.const_array(mfi);
    auto t_cas_arr = t_cas.const_array(mfi);
    auto q_cas_arr = q_cas.const_array(mfi);
    auto ustar_arr = lsm_fab_vars[LsmVar_SLM::ustar]->const_array(mfi);
    auto tstar_arr = lsm_fab_vars[LsmVar_SLM::tstar]->const_array(mfi);

    auto w_s_WP_arr = lsm_fab_vars[LsmVar_SLM::w_s_WP]->const_array(mfi);
    auto w_s_FC_arr = lsm_fab_vars[LsmVar_SLM::w_s_FC]->const_array(mfi);
    auto rootF_arr = lsm_fab_vars[LsmVar_SLM::rootF]->const_array(mfi);
    auto poro_soil_arr = lsm_fab_vars[LsmVar_SLM::poro_soil]->const_array(mfi);
    auto theta_FC_arr = lsm_fab_vars[LsmVar_SLM::theta_FC]->const_array(mfi);
    auto theta_WP_arr = lsm_fab_vars[LsmVar_SLM::theta_WP]->const_array(mfi);

    auto ztop_arr = ztop.const_array(mfi);
    auto Rgl_arr = Rgl.const_array(mfi);
    auto Rc_min_arr = Rc_min.const_array(mfi);
    auto hs_rc_arr = hs_rc.const_array(mfi);
    auto net_rad_arr = net_rad.const_array(mfi);

    auto r_a_arr = r_a.const_array(mfi);
    auto r_b_arr = r_b.array(mfi);
    auto r_c_arr = r_c.array(mfi);
    auto r_d_arr = r_d.array(mfi);

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
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
            amrex::Real rc_fac_rad, rc_fac_vpd, rc_fac_t, rc_fac_sw, d_root;

            // Aerodynamic resistance for heat and vapor transfer under canopy space : r_d
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
            r_d_arr(i, j, 0) = 1.0 / ustar_arr(i, j, 0) / Cs;

            // ===================================================
            // Leaf boundary layer resistance : r_b
            // ===================================================
            // turbulent transfer coefficient between canopy surface and canopy air : Cv = 0.01m/s^-0.5
            // characteristic dimension of the elaves in the direction of wind flux : d_leaf = 0.04m
            r_b_arr(i, j, 0) = 1.0 / 0.01 * std::pow( ustar_arr(i, j, 0) / 0.04, -0.5) / std::max(0.1, LAI_arr(i, j, 0));

            // ===================================================
            // Stomatal resistance : r_c
            // ===================================================
            // radiation factor
            // TODO: check if this is the correct downwelling SW to use
            amrex::Real tmp_radf = 0.55 * net_rad_arr(i, j, 0, SLM_NetRad::net_swdn1) * 2.0 / Rgl_arr(i, j, 0) / LAI_arr(i, j, 0);
            rc_fac_rad = (Rc_min_arr(i, j, 0) / d_Rc_max + tmp_radf) / (1.0 + tmp_radf);

            // vapor pressure deficit factor
            amrex::Real qsatw;
            erf_qsatw(t_cas_arr(i, j, 0), pref_arr(i, j, 0), qsatw);
            rc_fac_vpd = 1.0 / (1.0 + hs_rc_arr(i, j, 0) * (qsatw - q_cas_arr(i, j, 0)));

            // temperature factor
            rc_fac_t = 1.0 - 0.0016 * std::pow(d_T_opt - t_cas_arr(i, j, 0), 2);

            // rootzone soil moisture factor
            rc_fac_sw = 0.0;
            d_root = 0.0;
            for (int k = 0; k < d_nz_lsm; k++) {
                const int lsm_k = d_khi_lsm - k;

                if (rootF_arr(i, j, lsm_k) > 0.0)
                {
                    amrex::Real tmp;
                    // soil layer with the root
                    if (soilw_arr(i, j, lsm_k) > w_s_FC_arr(i, j, lsm_k))
                    {
                        // no water stree if water level exceeds the value at field capacity
                        tmp = 1.0 * s_depth_arr(i, j, lsm_k);
                        d_root += s_depth_arr(i, j, lsm_k);
                    }
                    else if (soilw_arr(i, j, lsm_k) < w_s_WP_arr(i, j, lsm_k))
                    {
                        // below wilting point
                        tmp = 0.0;
                    }
                    else
                    {
                        // otherwise
                        tmp = s_depth_arr(i, j, lsm_k) * (soilw_arr(i, j, lsm_k) * poro_soil_arr(i, j, lsm_k) - theta_WP_arr(i, j, lsm_k)) / (theta_FC_arr(i, j, lsm_k) - theta_WP_arr(i, j, lsm_k));
                        d_root += s_depth_arr(i, j, lsm_k);
                    }
                    rc_fac_sw += tmp;
                }
            }

            rc_fac_sw = rc_fac_sw / std::max(1.0e-6, d_root);

            tmp_radf = std::max(1.0e-6, rc_fac_rad*rc_fac_vpd*rc_fac_t*rc_fac_sw);
            r_c_arr(i, j, 0) = std::min(d_Rc_max, Rc_min_arr(i, j, 0) / LAI_arr(i, j, 0) / tmp_radf);

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

void SLM::vapor_fluxes(const amrex::MFIter &mfi)
{
    const int d_khi_lsm = khi_lsm;
    const int d_klo_lsm = klo_lsm;
    const int d_nz_lsm = m_nz_lsm;
    const Real dt = m_dt;

    auto box = mfi.tilebox();

    auto landmask_arr = landmask.const_array(mfi);

    auto q_cas_arr = q_cas.const_array(mfi);

    auto soilt_arr = lsm_fab_vars[LsmVar_SLM::soilt]->const_array(mfi);
    auto soilw_arr = lsm_fab_vars[LsmVar_SLM::soilw]->const_array(mfi);
    auto pref_arr  = lsm_fab_vars[LsmVar_SLM::pref]->const_array(mfi);

    auto lhf_canop_arr = lhf_canop.array(mfi);
    auto lhf_soil_arr = lhf_soil.array(mfi);
    auto lhf_air_arr = lhf_air.array(mfi);

    auto vegetype_arr = vegetype.const_array(mfi);
    auto vege_YES_arr = vege_YES.const_array(mfi);

    auto precip_array  = lsm_fab_vars[LsmVar_SLM::precipref]->const_array(mfi);

    auto mw_arr = mw.const_array(mfi);
    auto mw_mx_arr = mw_mx.const_array(mfi);
    auto mw_inc_arr = mw_inc.array(mfi);
    auto mws_arr = mws.const_array(mfi);

    auto dref_arr = lsm_fab_vars[LsmVar_SLM::dref]->const_array(mfi);
    auto qr_arr = lsm_fab_vars[LsmVar_SLM::qref]->const_array(mfi);
    auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
    auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);
    auto rootF_arr = lsm_fab_vars[LsmVar_SLM::rootF]->const_array(mfi);
    auto w_s_FC_arr = lsm_fab_vars[LsmVar_SLM::w_s_FC]->const_array(mfi);
    auto t_canop_arr = t_canop.const_array(mfi);

    auto r_a_arr = r_a.const_array(mfi);
    auto r_b_arr = r_b.const_array(mfi);
    auto r_c_arr = r_c.const_array(mfi);
    auto r_d_arr = r_d.const_array(mfi);
    auto r_soil_arr = r_soil.array(mfi);

    auto wet_canop_arr = wet_canop.array(mfi);
    auto evapo_dry_arr = evapo_dry.array(mfi);

    auto slm_diag_arr = slm_diag.array(mfi);

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        amrex::Real q_sfc;
        amrex::Real q_gr = 0.0;
        amrex::Real qref_tmp;

        // SAM rhow[nz] = air density at vertical velocity levels, kg/m^3
        const amrex::Real rhow = dref_arr(i, j, 0); // TODO: double check this

        // Specific humidity at top soil
        if (soilt_arr(i, j, d_khi_lsm) > tfriz)
        {
            if (mws_arr(i, j, 0) > 0.0)
            {
                erf_qsatw(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
            } else {
                erf_qsatw(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
                q_gr *= fh_calc(soilt_arr(i, j, d_khi_lsm), m_pot_sat_arr(i, j, d_khi_lsm), soilw_arr(i, j, d_khi_lsm), Bconst_arr(i, j, d_khi_lsm));
            }
        }
        else
        {
            erf_qsati(soilt_arr(i, j, d_khi_lsm), pref_arr(i, j, 0), q_gr);
        }

        // determine input q_sfc
        if (vegetype_arr(i, j, 0) != 0)
        {
            // vegetated surfaces
            q_sfc = q_cas_arr(i, j, 0);

            qref_tmp = q_sfc;
        }
        else
        {
            // baresoil
            q_sfc = q_gr;

            qref_tmp = qr_arr(i, j, 0);
        }

        amrex::Real soil_diff, totalR_soil, evapo_s, qsat_canop, evapo_wet;

        // Soil diffusion factor
        if (soilw_arr(i, j, d_khi_lsm) >= w_s_FC_arr(i,j,d_khi_lsm) || (qref_tmp > q_gr))
        {
            soil_diff = 1.0;
        } else {
            soil_diff = 0.25 * (std::pow(1.0 - cos(M_PI*std::max(0.01, soilw_arr(i, j, d_khi_lsm))/w_s_FC_arr(i, j, d_khi_lsm)), 2));
        }

        // total resistance for the topsoil
        if (vegetype_arr(i, j, 0) == 0)
        {
            // for baresoil
            r_soil_arr(i, j, 0) = std::min(10000.0, std::max(100.0, r_a_arr(i, j, 0)*(1.0 / soil_diff - 1.0)));
            totalR_soil = r_soil_arr(i, j, 0) + r_a_arr(i, j, 0) + r_litter;
        } else {
            r_soil_arr(i, j, 0) = std::min(10000.0, std::max(50.0, r_d_arr(i, j, 0)*(1.0 / soil_diff - 1.0)));
            totalR_soil = r_soil_arr(i, j, 0) + r_d_arr(i, j, 0) + r_litter;
        }

        // Evaporation from soil top underneath the canopy
        evapo_s = vege_YES_arr(i, j, 0) * (q_gr - q_sfc) * rhow / totalR_soil;

        // Evaporation from canopy
        wet_canop_arr(i, j, 0) = 1.0;
        erf_qsatw(t_canop_arr(i, j, 0), pref_arr(i, j, 0), qsat_canop);
        if (vegetype_arr(i, j, 0) != 0 && qsat_canop > q_sfc)
        {
            // wet portion of canopy
            wet_canop_arr(i, j, 0) = min(1.0, mw_arr(i, j, 0)/mw_mx_arr(i, j, 0));
        }

        if (vegetype_arr(i, j, 0) == 0)
        {
            wet_canop_arr(i, j, 0) = 0.0;
        }

        // direct evaporation from the water held on canopy
        evapo_wet = (qsat_canop - q_sfc)*rhow*wet_canop_arr(i, j, 0)/(2.0 * r_b_arr(i, j, 0))*vege_YES_arr(i, j, 0);
        slm_diag_arr(i, j, 0, SLM_Diag::evapo_wet) = evapo_wet;

        // increment/decrement of the water amount held on leaves following the direct evaporation/dew formation
        mw_inc_arr(i, j, 0) = -dt*evapo_wet; // evapo_wet [kg/m2s=mm/s]

        // Transpiration
        // For dew formation situation (qsat_canop < q_sfc), wet_canop = 1. automatically makes evapo_dry = 0
        evapo_dry_arr(i, j, 0) = (qsat_canop - q_sfc)*rhow*(1.0 - wet_canop_arr(i, j, 0))/(2.0 * r_b_arr(i, j, 0) + r_c_arr(i, j, 0))*vege_YES_arr(i, j, 0);

        // Check soil moisture availability for transpiration
        for (int k = 0; k < d_nz_lsm; k++) { // TODO: switch to regular lsm k loop
            const int lsm_k = d_khi_lsm - k;
            if (soilw_arr(i, j, lsm_k) < 0.1)
            {
                evapo_dry_arr(i, j, 0) -= evapo_dry_arr(i, j, 0) * rootF_arr(i, j, lsm_k);
            }
        }

        // Convert evaporation (kg/m2/s) to latent heat flux (W/m2)
        lhf_canop_arr(i, j, 0) = lcond*(evapo_wet+evapo_dry_arr(i, j, 0));
        lhf_soil_arr(i, j, 0) = lcond*evapo_s;
        lhf_air_arr(i, j, 0) = lcond*(q_sfc-qr_arr(i, j, 0))*rhow / r_a_arr(i, j, 0);

        if (vegetype_arr(i, j, 0) == 0)
        {
            // For baresoil case, lhf_air is equal to lhf_soil, lhf_air is updated to include soil diffusion factor
            lhf_air_arr(i, j, 0) *= r_a_arr(i, j, 0) / (r_soil_arr(i, j, 0) + r_a_arr(i, j, 0));
            evapo_s = lhf_air_arr(i, j, 0) / lcond;
            lhf_soil_arr(i, j, 0) = lhf_air_arr(i, j, 0);
        } else {
            lhf_air_arr(i, j, 0) = lhf_canop_arr(i, j, 0) + lhf_soil_arr(i, j, 0);
        }
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
    auto precip_extinc_arr = precip_extinc.const_array(mfi);

    auto mw_arr = mw.array(mfi);
    auto mw_mx_arr = mw_mx.const_array(mfi);
    auto mw_inc_arr = mw_inc.array(mfi);
    auto mws_arr = mws.array(mfi);
    auto mws_mx_arr = mws_mx.const_array(mfi);

    auto evapo_dry_arr = evapo_dry.const_array(mfi);
    auto prsfc_arr  = lsm_fab_vars[LsmVar_SLM::prsfc]->array(mfi);
    auto m_pot_sat_arr = lsm_fab_vars[LsmVar_SLM::m_pot_sat]->const_array(mfi);
    auto ks_arr = lsm_fab_vars[LsmVar_SLM::ks]->const_array(mfi);
    auto Bconst_arr = lsm_fab_vars[LsmVar_SLM::Bconst]->const_array(mfi);
    auto w_s_WP_arr = lsm_fab_vars[LsmVar_SLM::w_s_WP]->const_array(mfi);
    auto rootF_arr = lsm_fab_vars[LsmVar_SLM::rootF]->const_array(mfi);

    auto slm_diag_arr = slm_diag.array(mfi);

    auto dsw_vars = soilw_vars.array(mfi);

    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        amrex::Real drain = 0.0, puddle = 0.0;
        amrex::Real aa, bb, cc, dd;
        amrex::Real precip = 0.0;
        amrex::Real precip_in = 0.0;
        amrex::Real mws_inc = 0.0;

        // precipitation interception rate at canoppy
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

        amrex::Real precip_sfc = precip_array(i, j, 0) - precip + drain;

        // Update output variables
        prsfc_arr(i, j, 0) = precip_sfc;
        slm_diag_arr(i, j, 0, SLM_Diag::precip_sfc) = precip_sfc;
        slm_diag_arr(i, j, 0, SLM_Diag::drain) = drain;
        slm_diag_arr(i, j, 0, SLM_Diag::precip) = precip;

        // Add to vegetiation moisture increment from SLM::vapor_fluxes()
        mw_inc_arr(i, j, 0) += dt * (precip - drain);

        // Update vegetation moisture storage
        mw_arr(i, j, 0) += mw_inc_arr(i, j, 0);

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

            // Calculate precipitation infiltration rate into the first soil layer
            puddle = mws_arr(i, j, 0) / dt;
            if (soilw_arr(i, j, d_khi_lsm) < 1.0 && soilt_arr(i, j, d_khi_lsm) >= tfriz)
            {
                precip_in = std::min(precip_sfc + puddle, ks_arr(i, j, d_khi_lsm));
            } else {
                precip_in = 0.0;
            }

            puddle += precip_sfc - precip_in;
            drain = std::max(0.0, puddle - mws_mx_arr(i, j, 0) / dt);
            puddle -= drain;

            mws_inc = puddle*dt - mws_arr(i, j, 0);

            // Update puddle water storage
            mws_arr(i, j, 0) += mws_inc;

            //run_off_sfc = drain;
            slm_diag_arr(i, j, 0, SLM_Diag::run_off_sfc) = drain;
            slm_diag_arr(i, j, 0, SLM_Diag::precip_in) = precip_in;

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

            aa = 0.0;
            cc = -2.0 * dsw_vars(i, j, d_khi_lsm, SLM_DSW::sh_eff_cond)*dt/(dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm)*(dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm) + dsw_vars(i, j, d_khi_lsm - 1, SLM_DSW::sdepth_mm)));
            bb = 1.0 - cc + dsw_vars(i, j, d_khi_lsm, SLM_DSW::sh_eff_vel)*dt/dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm);
            dd = soilw_arr(i, j, d_khi_lsm) - (lhf_soil_arr(i, j, 0)/lcond + rootF_arr(i, j, d_khi_lsm)*evapo_dry_arr(i, j, 0)*dsw_vars(i, j, d_khi_lsm, SLM_DSW::sw_wgt) - precip_in)*dt / poro_soil_arr(i, j, d_khi_lsm) / dsw_vars(i, j, d_khi_lsm, SLM_DSW::sdepth_mm);

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
                dd = soilw_arr(i, j, lsm_k) - rootF_arr(i, j, lsm_k) * evapo_dry_arr(i, j, 0) * dsw_vars(i, j, lsm_k, SLM_DSW::sw_wgt) * dt / poro_soil_arr(i, j, lsm_k) / sd; // current time step

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
            dd = soilw_arr(i, j, d_klo_lsm) - (rootF_arr(i, j, d_klo_lsm)*evapo_dry_arr(i, j, 0)*dsw_vars(i, j, d_klo_lsm, SLM_DSW::sw_wgt) + drainage_flux) * dt / poro_soil_arr(i, j, d_klo_lsm) / dsw_vars(i, j, d_klo_lsm, SLM_DSW::sdepth_mm);

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

            if (fix) {
                // fix the levels where wetness exceeds 1 preserving total water:
                dd = 0.0;

                // compute the excess
                for (int kk = 0; kk < d_nz_lsm; kk++) {
                    const int lsm_kk = d_khi_lsm - kk;
                    if (soilw_arr(i, j, lsm_kk) > 1.0)
                    {
                        dd = dd + (soilw_arr(i, j, lsm_kk) - 1.0)*poro_soil_arr(i, j, lsm_kk)*dsw_vars(i, j, lsm_kk, SLM_DSW::sdepth_mm);
                        soilw_arr(i, j, lsm_kk) = 1.0;
                    }
                }

                // distribute the excess among layer into deepest layers first:
                for (int kk = d_nz_lsm - 1; kk >= 0; kk--) {
                    const int lsm_kk = d_khi_lsm - kk;
                    if (soilw_arr(i, j, lsm_kk) < 1.0)
                    {
                        cc = std::min(dd, (1.0 - soilw_arr(i, j, lsm_kk))*poro_soil_arr(i, j, lsm_kk)*dsw_vars(i, j, lsm_kk, SLM_DSW::sdepth_mm));
                        soilw_arr(i, j, lsm_kk) += cc / (poro_soil_arr(i, j, lsm_kk)*dsw_vars(i, j, lsm_kk, SLM_DSW::sdepth_mm));
                        dd -= cc;
                        if (dd < 0.0)
                        {
                            break;
                        }
                    }
                }
            }
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

    auto shf_soil_arr = shf_soil.const_array(mfi);
    auto lhf_soil_arr = lhf_soil.const_array(mfi);

    auto net_rad_arr = net_rad.const_array(mfi);
    auto slm_diag_arr = slm_diag.array(mfi);

    auto dst_vars = soilt_vars.array(mfi);

    // TODO: Refactor this whole loop for GPU
    ParallelFor( box, [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
        if (landmask_arr(i, j, 0) != 1) {
            return;
        }

        amrex::Real temp, k_dry, k_sat, Ke;

        amrex::Real grflux0 = -1.0 * (net_rad_arr(i, j, 0, SLM_NetRad::net_rad2) - shf_soil_arr(i, j, 0) - lhf_soil_arr(i, j, 0));
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
                    k_sat *= std::pow(0.57, poro_soil_arr(i, j, lsm_k));
                }
                else {
                    k_sat *= std::pow(1.60, poro_soil_arr(i, j, lsm_k));
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
SLM::set_flux_inputs(const amrex::MultiFab* sw_lw_fluxes_in,
                     const amrex::MultiFab* zenith_in)
{
    int khi = khi_lsm;

    auto tsurf = lsm_fab_vars[LsmVar_SLM::tsurf];

    for ( MFIter mfi(*tsurf, TileNoZ()); mfi.isValid(); ++mfi) {
        const auto& box3d = mfi.tilebox();

        // Create a box with the same i,j bounds, but only at z = 0
        amrex::Box b2d = box3d;
        b2d.setRange(2, 0);

        auto sw_lw_fluxes_arr = sw_lw_fluxes_in->const_array(mfi);
        auto zenith_array     = zenith_in->const_array(mfi);

        auto slm_dir_sw_vis  = lsm_fab_vars[LsmVar_SLM::swdsvisxyref]->array(mfi);
        auto slm_dir_sw_nir  = lsm_fab_vars[LsmVar_SLM::swdsnirxyref]->array(mfi);
        auto slm_diff_sw_vis = lsm_fab_vars[LsmVar_SLM::swdsvisdxyref]->array(mfi);
        auto slm_diff_sw_nir = lsm_fab_vars[LsmVar_SLM::swdsnirdxyref]->array(mfi);

        auto slm_lw          = lsm_fab_vars[LsmVar_SLM::lwref]->array(mfi);
        auto slm_zenith      = lsm_fab_vars[LsmVar_SLM::coszrsxy]->array(mfi);

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            slm_dir_sw_vis(i, j, k) = sw_lw_fluxes_arr(i, j, k, 0);
            slm_dir_sw_nir(i, j, k) = sw_lw_fluxes_arr(i, j, k, 1);

            slm_diff_sw_vis(i, j, k) = sw_lw_fluxes_arr(i, j, k, 2);
            slm_diff_sw_nir(i, j, k) = sw_lw_fluxes_arr(i, j, k, 3);

            slm_lw(i, j, k) = sw_lw_fluxes_arr(i, j, k, 5);
            slm_zenith(i, j, k) = zenith_array(i, j, k, 0);
            //slm_zenith(i, j, k) = 1.0;


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

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                slm_sst(i, j, k) = sst_array(i, j, k, 0);
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

void SLM::writeSLM_Data(const PlotFileType plotfile_type, const amrex::Real time, const std::string plot_prefix, const int level_step)
{
    std::string plotfilename = amrex::Concatenate(plot_prefix + "2D_", level_step, 5);

    amrex::Geometry lsm_2d_geom;
    lsm_2d_geom.define( ba_lsm_2d.minimalBox(), m_lsm_geom.ProbDomain(), m_lsm_geom.Coord(), m_lsm_geom.isPeriodic());

    amrex::Vector<amrex::MultiFab*> mf_data;

    mf_data.push_back(&net_rad);

    mf_data.push_back(&mw);
    mf_data.push_back(&mws);
    mf_data.push_back(&t_canop);
    mf_data.push_back(&t_skin);
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

    IntVect ng(0, 0, 0);

    // Total number of output MFs: net_rad components + mf_data size - 1
    const int output_size = SLM_NetRad::NumVars + mf_data.size() - 1;
    MultiFab fab(ba_lsm_2d, net_rad.DistributionMap(), output_size, ng);
    MultiFab::Copy(fab, *(mf_data[0]), 0, 0, SLM_NetRad::NumVars, 0);
    for (int i = 1; i < mf_data.size(); i++)
    {
        MultiFab::Copy(fab, *(mf_data[i]), 0, i + SLM_NetRad::NumVars - 1, 1, 0);
    }


    amrex::Vector<std::string> varnames;
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

    AMREX_ALWAYS_ASSERT(varnames.size() == output_size);

    if (plotfile_type == PlotFileType::Amrex) {
        amrex::WriteSingleLevelPlotfile(plotfilename, fab, varnames, lsm_2d_geom, time, level_step);
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
