/**
 * \file ERF_InitGeowind.cpp
 */
#include <ERF.H>

//#include <ERF_Interpolation_1D.H>

using namespace amrex;

/**
 * Read and interpolate a geostrophic wind profile onto ERF vertical levels.
 *
 * @param input_file Path to the geostrophic wind profile
 * @param u_geos Host vector storing interpolated geostrophic x-wind
 * @param u_geos_d Device vector storing interpolated geostrophic x-wind
 * @param v_geos Host vector storing interpolated geostrophic y-wind
 * @param v_geos_d Device vector storing interpolated geostrophic y-wind
 * @param lgeom Geometry defining the vertical domain
 * @param zlev_stag Optional stretched-grid staggered height levels
 */
void ERF::init_geo_wind_profile (const std::string input_file,
                                 Vector<Real>& u_geos,
                                 Gpu::DeviceVector<Real>& u_geos_d,
                                 Vector<Real>& v_geos,
                                 Gpu::DeviceVector<Real>& v_geos_d,
                                 const Geometry& lgeom,
                                 const Vector<Real>& zlev_stag)
{
    const int klo = 0;
    const int khi = lgeom.Domain().bigEnd()[AMREX_SPACEDIM-1];
    const amrex::Real dz = lgeom.CellSize()[AMREX_SPACEDIM-1];

    // zlevels_stag is filled on every mesh, so the ProbLo/ProbHi branch below is reached only
    // by a caller that passes no z levels.  The table's heights are physical heights, as the
    // z levels are: they start at geometry.prob_lo[2] on every mesh.
    const bool grid_stretch = (zlev_stag.size() > 0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!grid_stretch || zlev_stag.size() == khi+2,
        "init_geo_wind_profile: the z levels must be those of the level being filled");
    const Real zbot = (grid_stretch) ? zlev_stag[klo]   : lgeom.ProbLo(AMREX_SPACEDIM-1);
    const Real ztop = (grid_stretch) ? zlev_stag[khi+1] : lgeom.ProbHi(AMREX_SPACEDIM-1);

    amrex::Print() << "Reading geostrophic wind profile from " << input_file << std::endl;
    std::ifstream profile_reader(input_file);
    if(!profile_reader.is_open()) {
        amrex::Error("Error opening the abl_geo_wind_table\n");
    }

    // First, read the input data into temp vectors. Blank lines and lines starting
    // with '#' are skipped; any other line must hold exactly three numbers, z Ug Vg,
    // optionally followed by a '#' comment (the same rule as Exec/Tools/erf_input_checker).
    std::string line;
    Vector<Real> z_inp, Ug_inp, Vg_inp;
    Real z, Ug, Vg;
    int line_no = 0;
    amrex::Print() << "z  Ug  Vg" << std::endl;
    while(std::getline(profile_reader, line)) {
        ++line_no;
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') { continue; }
        std::istringstream iss(line);
        std::string rest;
        if (!(iss >> z >> Ug >> Vg) || ((iss >> rest) && rest[0] != '#')) {
            amrex::Abort("init_geo_wind_profile: line " + std::to_string(line_no) + " of " + input_file
                         + " (erf.abl_geo_wind_table) does not hold exactly three numbers z Ug Vg");
        }
        if (!z_inp.empty() && !(z > z_inp.back())) {
            amrex::Abort("init_geo_wind_profile: the heights in " + input_file
                         + " (erf.abl_geo_wind_table) must increase; see line " + std::to_string(line_no));
        }
        amrex::Print() << z << " " << Ug << " " << Vg << std::endl;
        z_inp.push_back(z);
        Ug_inp.push_back(Ug);
        Vg_inp.push_back(Vg);
        if (z >= ztop) break;
    }

    const int Ninp = z_inp.size();
    if (Ninp == 0) {
        amrex::Abort("init_geo_wind_profile: " + input_file + " (erf.abl_geo_wind_table) has no z Ug Vg rows");
    }
    if (z_inp[0] > zbot || z_inp[Ninp-1] < ztop) {
        amrex::Abort("init_geo_wind_profile: the rows of " + input_file + " (erf.abl_geo_wind_table) span "
                     + std::to_string(z_inp[0]) + " to " + std::to_string(z_inp[Ninp-1])
                     + " m and must cover the domain, " + std::to_string(zbot) + " to " + std::to_string(ztop) + " m");
    }

    // Now, interpolate vectors to the cell centers
    for (int k = 0; k <= khi; k++) {
        z = (grid_stretch) ? myhalf * (zlev_stag[k] + zlev_stag[k+1])
                           : zbot + (k + myhalf) * dz;
        u_geos[k] = interpolate_1d(z_inp.dataPtr(), Ug_inp.dataPtr(), z, Ninp);
        v_geos[k] = interpolate_1d(z_inp.dataPtr(), Vg_inp.dataPtr(), z, Ninp);
    }

    // Copy from host version to device version
    Gpu::copy(Gpu::hostToDevice, u_geos.begin(), u_geos.end(), u_geos_d.begin());
    Gpu::copy(Gpu::hostToDevice, v_geos.begin(), v_geos.end(), v_geos_d.begin());

    profile_reader.close();
}
