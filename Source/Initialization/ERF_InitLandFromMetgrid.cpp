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
#include <cfenv>
#include <cmath>
#include <map>
#include <regex>

using namespace amrex;

#ifdef ERF_USE_NETCDF

#include "ERF_NCInterface.H"
#include "ERF_MetgridLandUtils.H"

namespace {

// The depths of the Noah soil layers real.exe writes to wrfinput (ZS); their thicknesses
// (DZS) are metgrid_land::dzs_noah
using metgrid_land::nsoil_noah;
using metgrid_land::dzs_noah;
constexpr Real zs_noah [nsoil_noah] = {Real(0.05), Real(0.25), Real(0.7), Real(1.5)};

// A land soil value that is missing or out of range is filled in (see init_lsm_from_metgrid),
// unless there are more of them than this fraction of the land points
constexpr Real max_bad_soil_fraction = Real(0.01);

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

// The shape of a met_em variable: (Time, south_north, west_east) for a 2D field, or
// (Time, nz, south_north, west_east) with nz months or categories
struct MetgridVar {
    int found = 0;
    int nz = 1;
    int ny = 0;
    int nx = 0;
};

// The shapes of the given met_em variables, read on the I/O rank and broadcast
Vector<MetgridVar>
metgrid_var_info (const std::string& fname, const Vector<std::string>& names)
{
    const int n = names.size();
    Vector<int> info(4*n, 0);
    if (ParallelDescriptor::IOProcessor()) {
        auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
        for (int i = 0; i < n; ++i) {
            if (!ncf.has_var(names[i])) { continue; }
            const std::vector<size_t> shape = ncf.var(names[i]).shape();
            if (shape.size() == 3) {
                info[4*i]   = 1;
                info[4*i+1] = 1;
                info[4*i+2] = static_cast<int>(shape[1]);
                info[4*i+3] = static_cast<int>(shape[2]);
            } else if (shape.size() == 4) {
                info[4*i]   = 1;
                info[4*i+1] = static_cast<int>(shape[1]);
                info[4*i+2] = static_cast<int>(shape[2]);
                info[4*i+3] = static_cast<int>(shape[3]);
            } else {
                Abort("init_lsm_from_metgrid: " + names[i] + " in " + fname + " has " +
                      std::to_string(shape.size()) + " dimensions; expected 3 or 4");
            }
        }
        ncf.close();
    }
    ParallelDescriptor::Bcast(info.data(), 4*n, ParallelDescriptor::IOProcessorNumber());

    Vector<MetgridVar> vars(n);
    for (int i = 0; i < n; ++i) {
        vars[i] = {info[4*i], info[4*i+1], info[4*i+2], info[4*i+3]};
    }
    return vars;
}

//
// Read the columns of a met_em variable over region (a box in the level's index space, one
// cell in z) into fab, whose k index is the month or category.  The file's (0,0) is the
// cell file_lo.  At most max_readers ranks read the file at the same time, as in
// read_region_from_wrfinput.  A rank with an empty region reads nothing, but every rank
// must call this.
//
void
read_metgrid_region (const std::string& fname, const std::string& name, const MetgridVar& var,
                     const IntVect& file_lo, const Box& region, FArrayBox& fab, int max_readers)
{
    const int nprocs = ParallelDescriptor::NProcs();
    const int myproc = ParallelDescriptor::MyProc();
    const int nwaves = (max_readers > 0) ? (nprocs + max_readers - 1) / max_readers : 1;

    for (int wave = 0; wave < nwaves; ++wave) {
        if ((myproc % nwaves) == wave && region.ok()) {
            const int ilo = region.smallEnd(0) - file_lo[0];
            const int jlo = region.smallEnd(1) - file_lo[1];
            const int ni  = region.length(0);
            const int nj  = region.length(1);
            const int nk  = var.nz;

            // Only the first time is ever used
            std::vector<size_t> start, count;
            if (var.nz == 1) {
                start = {0, size_t(jlo), size_t(ilo)};
                count = {1, size_t(nj), size_t(ni)};
            } else {
                start = {0, 0, size_t(jlo), size_t(ilo)};
                count = {1, size_t(nk), size_t(nj), size_t(ni)};
            }

            std::vector<float> buf(size_t(ni) * nj * nk);
            {
                auto prev_fpe_except = amrex::disableFPExcept(amrex::FPExcept::all);
                auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
                ncf.var(name).get(buf.data(), start, count);
                ncf.close();
                std::feclearexcept(FE_ALL_EXCEPT);
                amrex::setFPExcept(prev_fpe_except);
            }

            Box bx = region;
            bx.setRange(2, 0, nk);
#ifdef AMREX_USE_GPU
            fab.resize(bx, 1, The_Pinned_Arena());
#else
            fab.resize(bx, 1);
#endif
            const Array4<Real>& arr = fab.array();
            const Dim3 lo = lbound(bx);
            for (int k = 0; k < nk; ++k) {
                for (int j = 0; j < nj; ++j) {
                    for (int i = 0; i < ni; ++i) {
                        arr(lo.x+i,lo.y+j,k) = static_cast<Real>(buf[(size_t(k)*nj + j)*ni + i]);
                    }
                }
            }
        }
        if (nwaves > 1) { ParallelDescriptor::Barrier(); }
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
 * A land soil temperature or moisture that is missing or out of range (metgrid can leave
 * such points where its source data had none, typically along coastlines) is replaced by
 * TSK, or by the mean moisture of that layer over the valid land points; the run stops if
 * more than max_bad_soil_fraction of the land points need it.
 *
 * Each rank reads and builds only the columns of its own boxes.  SLM then makes its own
 * conversions from these fields in SLM::init_wrfinput_vars, as it does for a wrfinput run.
 *
 * @param lev   Integer specifying the current level
 * @param fname Path to this level's met_em file for the initial time
 * @param date  The date of that file, as "YYYY-MM-DD_hh:mm:ss"
 */
void
ERF::init_lsm_from_metgrid (int lev, const std::string& fname, const std::string& date)
{
    AMREX_ALWAYS_ASSERT(solverChoice.lsm_type == LandSurfaceType::SLM);

    Print() << "Initializing SLM at level " << lev << " from the land fields in " << fname << std::endl;

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
    // The met_em fields, with their shapes.  SST is optional; the others are required.
    //
    enum { i_lmask = 0, i_skin, i_luindex, i_soilctop, i_green, i_lai12m, i_first_layer };
    Vector<std::string> names = {"LANDMASK", "SKINTEMP", "LU_INDEX", "SOILCTOP", "GREENFRAC", "LAI12M"};
    for (int m = 0; m < nlay; ++m) {
        names.push_back(soil_layer_name("ST", layers[m]));
        names.push_back(soil_layer_name("SM", layers[m]));
    }
    const int i_sst = names.size();
    names.push_back("SST");

    const Vector<MetgridVar> vars = metgrid_var_info(fname, names);
    for (int i = 0; i < i_sst; ++i) {
        if (!vars[i].found) {
            Abort("init_lsm_from_metgrid: " + names[i] + " was not found in " + fname +
                  "; it is required to initialize SLM from met_em files");
        }
    }
    const bool have_sst = vars[i_sst].found;

    const int nmonths = vars[i_green].nz;
    const int ncat    = vars[i_soilctop].nz;
    if (nmonths != 12 || vars[i_lai12m].nz != 12) {
        Abort("init_lsm_from_metgrid: GREENFRAC and LAI12M in " + fname + " must have 12 months");
    }

    // The file covers this level's domain, its (0,0) at the domain's low corner, as for the
    // atmospheric fields in init_from_metgrid
    const Box& domain = boxes_at_level[lev][0];
    const IntVect file_lo = domain.smallEnd();
    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
        if (vars[i].found && (vars[i].nx < domain.length(0) || vars[i].ny < domain.length(1))) {
            Abort("init_lsm_from_metgrid: " + names[i] + " in " + fname + " is " +
                  std::to_string(vars[i].nx) + " x " + std::to_string(vars[i].ny) +
                  ", smaller than level " + std::to_string(lev) + "'s " +
                  std::to_string(domain.length(0)) + " x " + std::to_string(domain.length(1)));
        }
    }

    //
    // The columns this rank needs: those of its boxes of SLM's input fields (all on the
    // same grids), and the bounding box of them, which is what it reads
    //
    const MultiFab* lsm_mf0 = nullptr;
    for (const auto& [wrfname, lsmname] : lsm.Get_WRFInputNames()) {
        int lsm_idx = lsm.Get_DataIdx(lev, lsmname);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lsm_idx != -1, "LSM variable mapping invalid!");
        AMREX_ALWAYS_ASSERT(lsm_data[lev][lsm_idx]);
        lsm_mf0 = lsm_data[lev][lsm_idx];
        break;
    }
    AMREX_ALWAYS_ASSERT(lsm_mf0 != nullptr);

    Vector<Box> my_columns;
    Box region;
    for (MFIter mfi(*lsm_mf0); mfi.isValid(); ++mfi) {
        Box b = mfi.validbox();
        b.setRange(2, 0, 1);
        b &= makeSlab(domain, 2, 0);
        if (!b.ok()) { continue; }
        my_columns.push_back(b);
        region = region.ok() ? region.minBox(b) : b;
    }

    int max_readers = 256;
    {
        ParmParse pp("erf");
        pp.query("metgrid_max_readers", max_readers);
    }

    Vector<FArrayBox> NC_fab(names.size());
    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
        if (vars[i].found) {
            read_metgrid_region(fname, names[i], vars[i], file_lo, region, NC_fab[i], max_readers);
        }
    }

    // read_from_metgrid gives every rank the date
    int m0 = 0, m1 = 0;
    Real w1 = zero;
    metgrid_land::monthly_interp_weights(date, m0, m1, w1);
    Print() << "  monthly fields interpolated between months " << m0+1 << " and " << m1+1
            << " with weight " << w1 << " on the later" << std::endl;

    //
    // Build the wrfinput fields on the host over this rank's region.  The soil fields have
    // one component per layer, with k = 0 the top layer as in wrfinput.  A land soil value
    // that is missing or out of range is marked with -1 here and filled in below.
    //
    Arena* arena = The_Arena();
#ifdef AMREX_USE_GPU
    // Make sure these live on the CPU, where they are filled, and on the GPU, where they
    // are copied from
    arena = The_Pinned_Arena();
#endif
    const Box bx2d = region;
    Box bxsoil = bx2d; bxsoil.setRange(2, 0, nsoil_noah);

    // A rank with no boxes here has empty fields, which it never reads
    const bool have_region = region.ok();
    std::map<std::string, FArrayBox> wrf;
    for (const char* name : {"TSK", "IVGTYP", "ISLTYP", "VEGFRA", "SHDMIN", "SHDMAX", "LAI"}) {
        if (have_region) { wrf[name].resize(bx2d, 1, arena); } else { wrf[name]; }
    }
    for (const char* name : {"TSLB", "SMOIS", "DZS", "ZS"}) {
        if (have_region) { wrf[name].resize(bxsoil, 1, arena); } else { wrf[name]; }
    }

    const auto lmask    = have_region ? NC_fab[i_lmask   ].const_array() : Array4<const Real>{};
    const auto skin     = have_region ? NC_fab[i_skin    ].const_array() : Array4<const Real>{};
    const auto luindex  = have_region ? NC_fab[i_luindex ].const_array() : Array4<const Real>{};
    const auto soilctop = have_region ? NC_fab[i_soilctop].const_array() : Array4<const Real>{};
    const auto green    = have_region ? NC_fab[i_green   ].const_array() : Array4<const Real>{};
    const auto lai12m   = have_region ? NC_fab[i_lai12m  ].const_array() : Array4<const Real>{};
    const auto sst      = (have_region && have_sst) ? NC_fab[i_sst].const_array() : Array4<const Real>{};

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

    constexpr Real bad = Real(-1.0);
    auto is_land = [&] (int i, int j) { return lmask(i,j,0) > Real(0.5); };

    // The dominant soil category of a land point in SOILCTOP, soil water excluded, or -1 if
    // SOILCTOP gives it none
    auto dominant_soil = [&] (int i, int j) -> int
    {
        int dominant = -1;
        Real max_frac = zero;
        for (int c = 0; c < ncat; ++c) {
            if (c+1 == isoilwater) { continue; }
            if (soilctop(i,j,c) > max_frac) {
                max_frac = soilctop(i,j,c);
                dominant = c+1;
            }
        }
        return dominant;
    };

    Vector<Real> st(nlay), sm(nlay);
    if (have_region) {
        const auto lo = lbound(bx2d);
        const auto hi = ubound(bx2d);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            const bool land = is_land(i,j);

            tsk(i,j,0) = skin(i,j,0);
            if (!land && have_sst && sst(i,j,0) > Real(170.)) {
                tsk(i,j,0) = sst(i,j,0);
            }

            ivgtyp(i,j,0) = std::round(luindex(i,j,0));

            if (land) {
                // No soil data at this land point: use clay loam, as real.exe does
                const int dominant = dominant_soil(i,j);
                isltyp(i,j,0) = (dominant < 0) ? 8 : dominant;
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
                st[m] = NC_fab[i_first_layer + 2*m    ].const_array()(i,j,0);
                sm[m] = NC_fab[i_first_layer + 2*m + 1].const_array()(i,j,0);
            }
            for (int k = 0; k < nsoil_noah; ++k) {
                dzs(i,j,k) = dzs_noah[k];
                zs (i,j,k) = zs_noah[k];
                if (land) {
                    const Real t = metgrid_land::interp_in_depth(z_src, st, zs_noah[k]);
                    const Real w = metgrid_land::interp_in_depth(z_src, sm, zs_noah[k]);
                    const bool t_ok = std::isfinite(t) && t >= Real(170.) && t <= Real(400.);
                    const bool w_ok = std::isfinite(w) && w >= zero && w <= one;
                    tslb (i,j,k) = t_ok ? t : bad;
                    smois(i,j,k) = w_ok ? w : bad;
                } else {
                    tslb (i,j,k) = tsk(i,j,0);
                    smois(i,j,k) = one;
                }
            }
        }
        }
    }

    //
    // Count, over the columns this rank owns, the land points and those with a missing or
    // out-of-range soil value, and sum the valid moisture of each layer
    //
    Long n_land = 0, n_bad = 0, n_no_soil_cat = 0;
    Vector<Real> sm_sum(nsoil_noah, zero);
    Vector<Long> sm_cnt(nsoil_noah, 0);
    for (const Box& b : my_columns) {
        const auto lo = lbound(b);
        const auto hi = ubound(b);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (!is_land(i,j)) { continue; }
            ++n_land;
            if (dominant_soil(i,j) < 0) { ++n_no_soil_cat; }
            bool any_bad = false;
            for (int k = 0; k < nsoil_noah; ++k) {
                if (tslb(i,j,k) == bad || smois(i,j,k) == bad) { any_bad = true; }
                if (smois(i,j,k) != bad) {
                    sm_sum[k] += smois(i,j,k);
                    ++sm_cnt[k];
                }
            }
            if (any_bad) { ++n_bad; }
        }
        }
    }
    {
        Long counts[3] = {n_land, n_bad, n_no_soil_cat};
        ParallelDescriptor::ReduceLongSum(counts, 3);
        n_land = counts[0]; n_bad = counts[1]; n_no_soil_cat = counts[2];
    }
    ParallelDescriptor::ReduceRealSum(sm_sum.data(), nsoil_noah);
    ParallelDescriptor::ReduceLongSum(sm_cnt.data(), nsoil_noah);

    if (n_bad > max_bad_soil_fraction * Real(n_land)) {
        Abort("init_lsm_from_metgrid: " + std::to_string(n_bad) + " of the " +
              std::to_string(n_land) + " land points in " + fname + " have a soil value that "
              "is missing or out of range (temperature outside 170-400 K or moisture outside "
              "0-1); at most " +
              std::to_string(std::lround(max_bad_soil_fraction * Real(100.))) +
              "% of them are filled in");
    }
    if (n_bad > 0) {
        Print() << "  " << n_bad << " of the " << n_land << " land points have a soil value "
                << "that is missing or out of range; using the skin temperature there for "
                << "temperature, and the mean over the other land points for moisture" << std::endl;
    }
    if (n_no_soil_cat > 0) {
        Print() << "  " << n_no_soil_cat << " land points have no soil category in SOILCTOP; "
                << "using category 8 (clay loam) there" << std::endl;
    }

    if (n_bad > 0 && have_region) {
        Vector<Real> sm_mean(nsoil_noah, zero);
        for (int k = 0; k < nsoil_noah; ++k) {
            if (sm_cnt[k] > 0) { sm_mean[k] = sm_sum[k] / Real(sm_cnt[k]); }
        }
        const auto lo = lbound(bx2d);
        const auto hi = ubound(bx2d);
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            for (int k = 0; k < nsoil_noah; ++k) {
                if (tslb (i,j,k) == bad) { tslb (i,j,k) = tsk(i,j,0); }
                if (smois(i,j,k) == bad) { smois(i,j,k) = sm_mean[k]; }
            }
        }
        }
    }

    //
    // Copy into the LSM's input fields, as init_from_wrfinput does: the soil layers go to
    // LSM levels k = -1 (top) down to -nsoil, and the 2D fields to the top level.
    //
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
            AMREX_ALWAYS_ASSERT(src.box().contains(tbx));
            const Array4<      Real>& dst_arr = lsm_data[lev][lsm_idx]->array(mfi);
            const Array4<const Real>& src_arr = src.const_array();
            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dst_arr(i,j,lsm_khi-k) = src_arr(i,j,k);
            });
        }
        lsm_data[lev][lsm_idx]->FillBoundary(geom[lev].periodicity());
    }
    Gpu::streamSynchronize();
}

#endif // ERF_USE_NETCDF
