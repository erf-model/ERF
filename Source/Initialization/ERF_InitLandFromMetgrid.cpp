/**
 * \file ERF_InitLandFromMetgrid.cpp
 *
 * Fill the SLM land surface model's initial state from a met_em file.
 *
 * SLM's real-data initialization was written against the fields real.exe puts in a
 * wrfinput file (TSLB, SMOIS, DZS, ZS, IVGTYP, ISLTYP, TSK, VEGFRA, SHDMIN, SHDMAX, LAI);
 * see SLM::wrfinput_map and SLM::init_wrfinput_vars.  A met_em file holds the raw
 * ingredients real.exe builds those from, so here we do the same derivation and hand SLM
 * the result in the same form, so that it takes exactly the path it takes for a wrfinput
 * run.
 */
#include "ERF.H"
#include "ERF_Constants.H"
#include "ERF_LandSurface.H"

#include <algorithm>
#include <cmath>
#include <map>
#include <regex>

using namespace amrex;

#ifdef ERF_USE_NETCDF

#include "ERF_NCWpsFile.H"

namespace {

// The Noah soil layers real.exe writes to wrfinput (DZS, ZS) for a 4-layer LSM.
constexpr int nsoil_noah = 4;
constexpr Real dzs_noah[nsoil_noah] = {Real(0.1), Real(0.3), Real(0.6), Real(1.0)};
constexpr Real zs_noah [nsoil_noah] = {Real(0.05), Real(0.25), Real(0.7), Real(1.5)};

// One soil layer of a met_em file, e.g. ST010040 is the layer from 10 to 40 cm.
struct MetgridSoilLayer {
    int top_cm;
    int bot_cm;
};

// Name of the met_em variable holding the given soil layer, e.g. ("ST", {10,40}) -> "ST010040"
std::string
soil_layer_name (const std::string& prefix, const MetgridSoilLayer& layer)
{
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%s%03d%03d", prefix.c_str(), layer.top_cm, layer.bot_cm);
    return std::string(buf);
}

//
// Find the soil layers in a met_em file.  Which layers are present depends on the source
// dataset (GFS gives 0-10, 10-40, 40-100 and 100-200 cm; ERA5 gives 0-7, 7-28, 28-100 and
// 100-289 cm), so they are discovered from the variable names rather than assumed.  The
// temperature and moisture layers must be the same.  Read on the I/O rank and broadcast.
//
Vector<MetgridSoilLayer>
find_metgrid_soil_layers (const std::string& fname)
{
    Vector<int> depths;
    if (ParallelDescriptor::IOProcessor()) {
        std::map<std::pair<int,int>, int> found; // (top,bot) -> bit 1 for ST, bit 2 for SM
        const std::regex layer_re("^(ST|SM)([0-9]{3})([0-9]{3})$");
        auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
        for (const auto& var : ncf.all_vars()) {
            std::smatch m;
            const std::string name = var.name();
            if (std::regex_match(name, m, layer_re)) {
                const auto key = std::make_pair(std::stoi(m[2].str()), std::stoi(m[3].str()));
                found[key] |= (m[1].str() == "ST") ? 1 : 2;
            }
        }
        ncf.close();
        for (const auto& [key, bits] : found) {
            if (bits != 3) {
                Abort("init_lsm_from_metgrid: " + fname + " has soil " +
                      std::string((bits == 1) ? "temperature" : "moisture") + " for the layer " +
                      std::to_string(key.first) + "-" + std::to_string(key.second) +
                      " cm but not soil " + std::string((bits == 1) ? "moisture" : "temperature"));
            }
            depths.push_back(key.first);
            depths.push_back(key.second);
        }
    }

    int n = depths.size();
    ParallelDescriptor::Bcast(&n, 1, ParallelDescriptor::IOProcessorNumber());
    depths.resize(n);
    if (n > 0) {
        ParallelDescriptor::Bcast(depths.data(), n, ParallelDescriptor::IOProcessorNumber());
    }

    // std::map iterates in order of increasing top depth
    Vector<MetgridSoilLayer> layers;
    for (int i = 0; i < n; i += 2) {
        layers.push_back({depths[i], depths[i+1]});
    }
    return layers;
}

//
// Linear interpolation in depth from the source layer midpoints to z, held constant
// beyond the shallowest and deepest midpoints.
//
Real
interp_in_depth (const Vector<Real>& z_src, const Vector<Real>& v_src, Real z)
{
    const int n = z_src.size();
    if (z <= z_src[0])   { return v_src[0]; }
    if (z >= z_src[n-1]) { return v_src[n-1]; }
    int m = 0;
    while (z > z_src[m+1]) { ++m; }
    const Real w = (z - z_src[m]) / (z_src[m+1] - z_src[m]);
    return (one - w) * v_src[m] + w * v_src[m+1];
}

//
// The two months whose mid-month values bracket the given date, and the weight of the
// second, as in real.exe's monthly_interp_to_date: each monthly value is taken to apply on
// the 15th of its month, and values in between are interpolated linearly in time.
//
void
monthly_interp_weights (const std::string& date, int& m0, int& m1, Real& w1)
{
    // date is "YYYY-MM-DD_hh:mm:ss"
    const int year  = std::stoi(date.substr(0,4));
    const int month = std::stoi(date.substr(5,2));
    const int day   = std::stoi(date.substr(8,2));
    const int hour  = std::stoi(date.substr(11,2));

    const bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    int days_in_month[12] = {31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    // Day of the year (starting from 0) of the date and of the middle of each month
    Real mid[12];
    int start = 0;
    for (int m = 0; m < 12; ++m) {
        mid[m] = Real(start + 14);
        start += days_in_month[m];
    }
    const int days_in_year = start;
    int doy = day - 1;
    for (int m = 0; m < month-1; ++m) { doy += days_in_month[m]; }
    const Real t = Real(doy) + Real(hour) / Real(24.0);

    if (t < mid[0] || t >= mid[11]) {
        // Between mid-December and mid-January, wrapping around the end of the year
        m0 = 11;
        m1 = 0;
        const Real t0 = mid[11];
        const Real t1 = mid[0] + Real(days_in_year);
        const Real tt = (t < mid[0]) ? t + Real(days_in_year) : t;
        w1 = (tt - t0) / (t1 - t0);
    } else {
        m0 = 0;
        while (t >= mid[m0+1]) { ++m0; }
        m1 = m0 + 1;
        w1 = (t - mid[m0]) / (mid[m1] - mid[m0]);
    }
}

} // namespace

/**
 * Fill the SLM inputs at this level from a met_em file, building the fields real.exe would
 * write to wrfinput:
 *
 *   TSK            SKINTEMP, replaced by SST over water where SST is given
 *   IVGTYP         LU_INDEX
 *   ISLTYP         the dominant category of SOILCTOP over land (soil water, ISOILWATER, is
 *                  excluded there), ISOILWATER over water
 *   TSLB, SMOIS    the ST* and SM* soil layers interpolated in depth to the Noah layers;
 *                  over water, TSK and 1
 *   DZS, ZS        the thicknesses and depths of the Noah layers: 0.1, 0.3, 0.6 and 1.0 m
 *   VEGFRA         GREENFRAC interpolated to the date, in percent
 *   SHDMIN/SHDMAX  the smallest and largest monthly GREENFRAC, in percent
 *   LAI            LAI12M interpolated to the date
 *
 * SLM then makes its own conversions from these in SLM::init_wrfinput_vars, as it does for a
 * wrfinput run.
 *
 * @param lev   Integer specifying the current level
 * @param fname Path to the met_em file for the initial time
 * @param date  The date of that file, as "YYYY-MM-DD_hh:mm:ss" (needed on the I/O rank only)
 */
void
ERF::init_lsm_from_metgrid (int lev, const std::string& fname, const std::string& date)
{
    AMREX_ALWAYS_ASSERT(solverChoice.lsm_type == LandSurfaceType::SLM);

    Print() << "Initializing SLM from the land fields in " << fname << std::endl;

    // SLM's layers must be the Noah layers we interpolate to
    {
        const int nsoil = lsm.Get_Lsm_Geom(lev).Domain().length(2);
        const auto& soil_dz = lsm.get_model_lev<SLM>(lev)->get_soil_dz();
        bool match = (nsoil == nsoil_noah) && (static_cast<int>(soil_dz.size()) == nsoil_noah);
        for (int k = 0; match && k < nsoil_noah; ++k) {
            match = std::fabs(soil_dz[k] - dzs_noah[k]) < Real(1.e-6);
        }
        if (!match) {
            Abort("init_lsm_from_metgrid: SLM initialized from met_em files uses the 4 Noah "
                  "soil layers that real.exe writes to wrfinput; set slm.nsoil = 4 and "
                  "slm.soil_dz = 0.1 0.3 0.6 1.0");
        }
    }

    const Vector<MetgridSoilLayer> layers = find_metgrid_soil_layers(fname);
    if (layers.empty()) {
        Abort("init_lsm_from_metgrid: " + fname + " has no soil temperature or moisture "
              "layers (variables ST<top><bottom> and SM<top><bottom>, depths in cm)");
    }
    const int nlay = layers.size();
    Vector<Real> z_src(nlay);
    Print() << "  met_em soil layers (cm):";
    for (int m = 0; m < nlay; ++m) {
        z_src[m] = Real(0.005) * Real(layers[m].top_cm + layers[m].bot_cm);
        Print() << " " << layers[m].top_cm << "-" << layers[m].bot_cm;
    }
    Print() << std::endl;

    // Category numbers for water from the file's global attributes, with the defaults for
    // the MODIS land use and STATSGO soil categories that SLM's tables assume
    int isoilwater = 14;
    std::string mminlu;
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
        if (ncf.has_attr("ISOILWATER")) {
            std::vector<int> attr;
            ncf.get_attr("ISOILWATER", attr);
            isoilwater = attr[0];
        }
        if (ncf.has_attr("MMINLU")) { mminlu = ncf.get_attr("MMINLU"); }
        ncf.close();
    }
    ParallelDescriptor::Bcast(&isoilwater, 1, ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor() && !mminlu.empty() &&
        mminlu.find("MODIS") == std::string::npos) {
        Warning("init_lsm_from_metgrid: " + fname + " uses the land use categories " + mminlu +
                "; SLM's vegetation tables assume the MODIS categories");
    }

    //
    // Read the land fields.  Every rank gets the whole domain, kept on the host since the
    // derivation below is done there.  The z extent of the domain box bounds how many
    // levels are read, so make it large enough for the 12 months and 16 soil categories.
    //
    Box domain = boxes_at_level[lev][0];
    domain.setRange(2, 0, 64);

    FArrayBox NC_lmask_fab, NC_tsk_fab, NC_sst_fab, NC_luindex_fab;
    FArrayBox NC_soilctop_fab, NC_greenfrac_fab, NC_lai12m_fab;
    Vector<FArrayBox> NC_st_fab(nlay), NC_sm_fab(nlay);

    Vector<FArrayBox*>  NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;
    auto add = [&] (FArrayBox* fab, const std::string& name, NC_Data_Dims_Type type)
    {
        NC_fabs.push_back(fab); NC_names.push_back(name); NC_dim_types.push_back(type);
    };
    add(&NC_lmask_fab,     "LANDMASK",  NC_Data_Dims_Type::Time_SN_WE);
    add(&NC_tsk_fab,       "SKINTEMP",  NC_Data_Dims_Type::Time_SN_WE);
    add(&NC_luindex_fab,   "LU_INDEX",  NC_Data_Dims_Type::Time_SN_WE);
    add(&NC_soilctop_fab,  "SOILCTOP",  NC_Data_Dims_Type::Time_BT_SN_WE);
    add(&NC_greenfrac_fab, "GREENFRAC", NC_Data_Dims_Type::Time_BT_SN_WE);
    add(&NC_lai12m_fab,    "LAI12M",    NC_Data_Dims_Type::Time_BT_SN_WE);
    for (int m = 0; m < nlay; ++m) {
        add(&NC_st_fab[m], soil_layer_name("ST", layers[m]), NC_Data_Dims_Type::Time_SN_WE);
        add(&NC_sm_fab[m], soil_layer_name("SM", layers[m]), NC_Data_Dims_Type::Time_SN_WE);
    }
    const int n_required = NC_fabs.size();
    add(&NC_sst_fab,       "SST",       NC_Data_Dims_Type::Time_SN_WE);

    Vector<int> success(NC_fabs.size());
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs,
                                            success, /*keep_on_host=*/ true);
    for (int i = 0; i < n_required; ++i) {
        if (success[i] != 1) {
            Abort("init_lsm_from_metgrid: " + NC_names[i] + " was not found in " + fname +
                  "; it is required to initialize SLM from met_em files");
        }
    }
    const bool have_sst = (success[n_required] == 1);

    const int nmonths = NC_greenfrac_fab.box().length(2);
    const int ncat    = NC_soilctop_fab.box().length(2);
    if (nmonths != 12 || NC_lai12m_fab.box().length(2) != 12) {
        Abort("init_lsm_from_metgrid: GREENFRAC and LAI12M in " + fname + " must have 12 months");
    }

    // read_from_metgrid sets the date on the I/O rank only
    int m0 = 0, m1 = 0;
    Real w1 = zero;
    if (ParallelDescriptor::IOProcessor()) {
        monthly_interp_weights(date, m0, m1, w1);
    }
    ParallelDescriptor::Bcast(&m0, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&m1, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&w1, 1, ParallelDescriptor::IOProcessorNumber());
    Print() << "  monthly fields interpolated between months " << m0+1 << " and " << m1+1
            << " with weight " << w1 << " on the later" << std::endl;

    //
    // Build the wrfinput fields on the host.  The soil fields have one component per layer,
    // with k = 0 the top layer as in wrfinput.
    //
    Arena* arena = The_Arena();
#ifdef AMREX_USE_GPU
    // Make sure these live on the CPU, where they are filled, and on the GPU, where they
    // are copied from
    arena = The_Pinned_Arena();
#endif
    const Box bx2d  = makeSlab(NC_tsk_fab.box(), 2, 0);
    Box bxsoil = bx2d; bxsoil.setRange(2, 0, nsoil_noah);

    std::map<std::string, FArrayBox> wrf;
    for (const std::string name : {"TSK", "IVGTYP", "ISLTYP", "VEGFRA", "SHDMIN", "SHDMAX", "LAI"}) {
        wrf[name].resize(bx2d, 1, arena);
    }
    for (const std::string name : {"TSLB", "SMOIS", "DZS", "ZS"}) {
        wrf[name].resize(bxsoil, 1, arena);
    }

    const auto lmask    = NC_lmask_fab.const_array();
    const auto skin     = NC_tsk_fab.const_array();
    const auto luindex  = NC_luindex_fab.const_array();
    const auto soilctop = NC_soilctop_fab.const_array();
    const auto green    = NC_greenfrac_fab.const_array();
    const auto lai12m   = NC_lai12m_fab.const_array();
    const auto sst      = have_sst ? NC_sst_fab.const_array() : Array4<const Real>{};

    auto tsk    = wrf["TSK"].array();
    auto ivgtyp = wrf["IVGTYP"].array();
    auto isltyp = wrf["ISLTYP"].array();
    auto vegfra = wrf["VEGFRA"].array();
    auto shdmin = wrf["SHDMIN"].array();
    auto shdmax = wrf["SHDMAX"].array();
    auto lai    = wrf["LAI"].array();
    auto tslb   = wrf["TSLB"].array();
    auto smois  = wrf["SMOIS"].array();
    auto dzs    = wrf["DZS"].array();
    auto zs     = wrf["ZS"].array();

    Long n_bad_soil = 0;
    Long n_no_soil_cat = 0;
    Vector<Real> st(nlay), sm(nlay);

    const auto lo = lbound(bx2d);
    const auto hi = ubound(bx2d);
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        const bool land = (lmask(i,j,0) > Real(0.5));

        tsk(i,j,0) = skin(i,j,0);
        if (!land && have_sst && sst(i,j,0) > Real(170.)) {
            tsk(i,j,0) = sst(i,j,0);
        }

        ivgtyp(i,j,0) = std::round(luindex(i,j,0));

        if (land) {
            int dominant = -1;
            Real max_frac = zero;
            for (int c = 0; c < ncat; ++c) {
                if (c+1 == isoilwater) { continue; }
                if (soilctop(i,j,c) > max_frac) {
                    max_frac = soilctop(i,j,c);
                    dominant = c+1;
                }
            }
            if (dominant < 0) {
                // No soil data at this land point: use clay loam, as real.exe does
                dominant = 8;
                ++n_no_soil_cat;
            }
            isltyp(i,j,0) = dominant;
        } else {
            isltyp(i,j,0) = isoilwater;
        }

        Real gmin = green(i,j,0), gmax = green(i,j,0);
        for (int m = 1; m < nmonths; ++m) {
            gmin = std::min(gmin, green(i,j,m));
            gmax = std::max(gmax, green(i,j,m));
        }
        vegfra(i,j,0) = Real(100.) * ((one - w1) * green(i,j,m0) + w1 * green(i,j,m1));
        shdmin(i,j,0) = Real(100.) * gmin;
        shdmax(i,j,0) = Real(100.) * gmax;
        lai(i,j,0)    = (one - w1) * lai12m(i,j,m0) + w1 * lai12m(i,j,m1);

        for (int m = 0; m < nlay; ++m) {
            st[m] = NC_st_fab[m].const_array()(i,j,0);
            sm[m] = NC_sm_fab[m].const_array()(i,j,0);
        }
        for (int k = 0; k < nsoil_noah; ++k) {
            dzs(i,j,k) = dzs_noah[k];
            zs (i,j,k) = zs_noah[k];
            if (land) {
                tslb (i,j,k) = interp_in_depth(z_src, st, zs_noah[k]);
                smois(i,j,k) = interp_in_depth(z_src, sm, zs_noah[k]);
                if (!std::isfinite(tslb(i,j,k)) || tslb(i,j,k) < Real(170.) || tslb(i,j,k) > Real(400.) ||
                    !std::isfinite(smois(i,j,k)) || smois(i,j,k) < zero || smois(i,j,k) > one) {
                    ++n_bad_soil;
                }
            } else {
                tslb (i,j,k) = tsk(i,j,0);
                smois(i,j,k) = one;
            }
        }
    }
    }

    // Every rank has worked through the whole domain, so these counts agree on all of them
    if (n_bad_soil > 0) {
        Abort("init_lsm_from_metgrid: " + std::to_string(n_bad_soil) + " land soil values in " +
              fname + " are missing or out of range (temperature outside 170-400 K or moisture "
              "outside 0-1)");
    }
    if (n_no_soil_cat > 0) {
        Print() << "  " << n_no_soil_cat << " land points have no soil category in SOILCTOP; "
                << "using category 8 (clay loam) there" << std::endl;
    }

    //
    // Copy into the LSM's input fields, as init_from_wrfinput does: the soil layers go to
    // LSM levels k = -1 (top) down to -nsoil, and the 2D fields to the top level.
    //
    int i_lo = bx2d.smallEnd(0); int i_hi = bx2d.bigEnd(0);
    int j_lo = bx2d.smallEnd(1); int j_hi = bx2d.bigEnd(1);
    for (const auto& [wrfname, lsmname] : lsm.Get_WRFInputNames()) {
        auto it = wrf.find(wrfname);
        if (it == wrf.end()) {
            Abort("init_lsm_from_metgrid: SLM expects the wrfinput field " + wrfname +
                  ", which is not built from met_em files");
        }
        const FArrayBox& src = it->second;
        const int nk = src.box().length(2);

        int lsm_idx = lsm.Get_DataIdx(lev, lsmname);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lsm_idx != -1, "LSM variable mapping invalid!");
        AMREX_ALWAYS_ASSERT(lsm_data[lev][lsm_idx]);

        for (MFIter mfi(*lsm_data[lev][lsm_idx], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box tbx = mfi.tilebox();
            const int lsm_khi = tbx.bigEnd(2);
            tbx.setRange(2, 0, nk);
            const Array4<      Real>& dst_arr = lsm_data[lev][lsm_idx]->array(mfi);
            const Array4<const Real>& src_arr = src.const_array();
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int li = amrex::min(amrex::max(i, i_lo), i_hi);
                int lj = amrex::min(amrex::max(j, j_lo), j_hi);
                dst_arr(i,j,lsm_khi-k) = src_arr(li,lj,k);
            });
        }
        lsm_data[lev][lsm_idx]->FillBoundary(geom[lev].periodicity());
    }
    Gpu::streamSynchronize();
}

#endif // ERF_USE_NETCDF
