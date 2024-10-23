/**
 * \file ERF_init_geowind.cpp
 */
//#include <ERF_EOS.H>
#include <ERF.H>
//#include <ERF_TileNoZ.H>
//#include <ERF_prob_common.H>
//#include <ERF_ParFunctions.H>
//#include <ERF_Utils.H>

//#include <ERF_Interpolation_1D.H>

using namespace amrex;

void ERF::init_geo_wind_profile(const std::string input_file,
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

    const bool grid_stretch = (zlev_stag.size() > 0);
    const Real zbot = (grid_stretch) ? zlev_stag[klo]   : lgeom.ProbLo(AMREX_SPACEDIM-1);
    const Real ztop = (grid_stretch) ? zlev_stag[khi+1] : lgeom.ProbHi(AMREX_SPACEDIM-1);

    amrex::Print() << "Reading geostrophic wind profile from " << input_file << std::endl;
    std::ifstream profile_reader(input_file);
    if(!profile_reader.is_open()) {
        amrex::Error("Error opening the abl_geo_wind_table\n");
    }

    // First, read the input data into temp vectors
    std::string line;
    Vector<Real> z_inp, Ug_inp, Vg_inp;
    Real z, Ug, Vg;
    amrex::Print() << "z  Ug  Vg" << std::endl;
    while(std::getline(profile_reader, line)) {
        std::istringstream iss(line);
        iss >> z >> Ug >> Vg;
        amrex::Print() << z << " " << Ug << " " << Vg << std::endl;
        z_inp.push_back(z);
        Ug_inp.push_back(Ug);
        Vg_inp.push_back(Vg);
        if (z >= ztop) break;
    }

    const int Ninp = z_inp.size();
    AMREX_ALWAYS_ASSERT(z_inp[0] <= zbot);
    AMREX_ALWAYS_ASSERT(z_inp[Ninp-1] >= ztop);

    // Now, interpolate vectors to the cell centers
    for (int k = 0; k <= khi; k++) {
        z = (grid_stretch) ? 0.5 * (zlev_stag[k] + zlev_stag[k+1])
                           : zbot + (k + 0.5) * dz;
        u_geos[k] = interpolate_1d(z_inp.dataPtr(), Ug_inp.dataPtr(), z, Ninp);
        v_geos[k] = interpolate_1d(z_inp.dataPtr(), Vg_inp.dataPtr(), z, Ninp);
    }

    // Copy from host version to device version
    Gpu::copy(Gpu::hostToDevice, u_geos.begin(), u_geos.end(), u_geos_d.begin());
    Gpu::copy(Gpu::hostToDevice, v_geos.begin(), v_geos.end(), v_geos_d.begin());

    profile_reader.close();
}
