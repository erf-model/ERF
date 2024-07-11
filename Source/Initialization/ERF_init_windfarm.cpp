/**
 * \file ERF_init_windfarm.cpp
 */

#include <ERF.H>
#include <WindFarm.H>

using namespace amrex;

/**
 * Read in the turbine locations in latitude-longitude from windturbines.txt
 * and convert it into x and y coordinates in metres
 *
 * @param lev Integer specifying the current level
 */

void
ERF::init_windfarm (int lev)
{
    Vector<Real> xloc, yloc;
    if(solverChoice.windfarm_loc_type == WindFarmLocType::lat_lon) {
        init_windfarm_lat_lon(lev, solverChoice.windfarm_loc_table, xloc, yloc);
    }
    else if(solverChoice.windfarm_loc_type == WindFarmLocType::x_y) {
        init_windfarm_x_y(lev, solverChoice.windfarm_loc_table, xloc, yloc);
    }

    else {
        amrex::Abort("Are you using windfarms? For windfarm simulations, the inputs need to have an"
                     " entry erf.windfarm_loc_table which should not be None. \n");
    }

    // Write out a vtk file for turbine locations
    if (ParallelDescriptor::IOProcessor()){
        FILE* file_turbloc_vtk;
        file_turbloc_vtk = fopen("turbine_locations.vtk","w");
        fprintf(file_turbloc_vtk, "%s\n","# vtk DataFile Version 3.0");
        fprintf(file_turbloc_vtk, "%s\n","Wind turbine locations");
        fprintf(file_turbloc_vtk, "%s\n","ASCII");
        fprintf(file_turbloc_vtk, "%s\n","DATASET POLYDATA");
        fprintf(file_turbloc_vtk, "%s %ld %s\n", "POINTS", xloc.size(), "float");
        for(int it=0; it<xloc.size(); it++){
            fprintf(file_turbloc_vtk, "%0.15g %0.15g %0.15g\n", xloc[it], yloc[it], 1e-12);
        }
        fclose(file_turbloc_vtk);
    }

    amrex::Gpu::DeviceVector<Real> d_xloc(xloc.size());
    amrex::Gpu::DeviceVector<Real> d_yloc(yloc.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xloc.begin(), xloc.end(), d_xloc.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yloc.begin(), yloc.end(), d_yloc.begin());

    Real* d_xloc_ptr     = d_xloc.data();
    Real* d_yloc_ptr     = d_yloc.data();

    Nturb[lev].setVal(0);

    int i_lo = geom[lev].Domain().smallEnd(0); int i_hi = geom[lev].Domain().bigEnd(0);
    int j_lo = geom[lev].Domain().smallEnd(1); int j_hi = geom[lev].Domain().bigEnd(1);
    auto dx = geom[lev].CellSizeArray();
    int num_turb = xloc.size();

     // Initialize wind farm
    for ( MFIter mfi(Nturb[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx     = mfi.tilebox();
        auto  Nturb_array = Nturb[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int li = amrex::min(amrex::max(i, i_lo), i_hi);
            int lj = amrex::min(amrex::max(j, j_lo), j_hi);

            Real x1 = li*dx[0], x2 = (li+1)*dx[0];
            Real y1 = lj*dx[1], y2 = (lj+1)*dx[1];

            for(int it=0; it<num_turb; it++){
                if( d_xloc_ptr[it]+1e-12 > x1 and d_xloc_ptr[it]+1e-12 < x2 and
                    d_yloc_ptr[it]+1e-12 > y1 and d_yloc_ptr[it]+1e-12 < y2){
                       Nturb_array(i,j,k,0) = Nturb_array(i,j,k,0) + 1;
                }
            }
        });
    }

    read_in_table(solverChoice.windfarm_spec_table);
}

void
ERF::init_windfarm_lat_lon (const int lev,
                            const std::string windfarm_loc_table,
                            Vector<Real>& xloc,
                            Vector<Real>& yloc)
{

    if(solverChoice.latitude_lo == -1e10) {
        amrex::Error("Are you using latitude-longitude for wind turbine locations? There needs to be"
                     " entry for erf.latitude_lo in the inputs for"
                     " the lower bottom corner of the domain");
    }

    if(solverChoice.longitude_lo == -1e10) {
        amrex::Error("Are you using latitude-longitude for wind turbine locations? There needs to be"
                     " entry erf.longitude_lo in the inputs for"
                     " the lower bottom corner of the domain");
    }

    if(std::fabs(solverChoice.latitude_lo) > 90.0) {
        amrex::Error("The value of erf.latitude_lo in the input file should be within -90 and 90");
    }

    if(std::fabs(solverChoice.longitude_lo) > 180.0) {
        amrex::Error("The value of erf.longitude_lo in the input file should be within -180 and 180");
    }

    // Read turbine locations from windturbines.txt
    std::ifstream file(windfarm_loc_table);
    if (!file.is_open()) {
        amrex::Error("Wind turbines location table not found. Either the inputs is missing the"
                     " erf.windfarm_loc_table entry or the file specified in the entry " + windfarm_loc_table + " is missing.");
    }
    // Vector of vectors to store the matrix
    Vector<Real> lat, lon;
    Real value1, value2, value3;

    while (file >> value1 >> value2 >> value3) {

        if(std::fabs(value1) > 90.0) {
            amrex::Error("The value of latitude for entry " + std::to_string(lat.size() + 1) +
                         " in " + windfarm_loc_table + " should be within -90 and 90");
        }

        if(std::fabs(value2) > 180.0) {
            amrex::Error("The value of longitude for entry " + std::to_string(lat.size() + 1) +
                         " in " + windfarm_loc_table + " should be within -180 and 180");
        }

        if(std::fabs(solverChoice.longitude_lo) > 180.0) {
            amrex::Error("The values of erf.longitude_lo should be within -180 and 180");
        }
        lat.push_back(value1);
        lon.push_back(value2);
    }
    file.close();

    Real rad_earth = 6371.0e3; // Radius of the earth

    Real lat_lo  = solverChoice.latitude_lo*M_PI/180.0;
    Real lon_lo  = solverChoice.longitude_lo*M_PI/180.0;

    // (lat_lo, lon_lo) is mapped to (0,0)

    for(int it=0;it<lat.size();it++){
        lat[it] = lat[it]*M_PI/180.0;
        lon[it] = lon[it]*M_PI/180.0;
        Real delta_lat = (lat[it] - lat_lo);
        Real delta_lon = (lon[it] - lon_lo);

        Real term1 = std::pow(sin(delta_lat/2.0),2);
        Real term2 = cos(lat[it])*cos(lat_lo)*std::pow(sin(delta_lon/2.0),2);
        Real dist =  2.0*rad_earth*std::asin(std::sqrt(term1 + term2));
        Real dy_turb = (lat[it] - lat_lo) * 111000.0 * 180.0/M_PI ;
        yloc.push_back(dy_turb);
        Real dx_turb = std::sqrt(std::pow(dist,2) - std::pow(dy_turb,2));
        xloc.push_back(dx_turb);
    }
}

void
ERF::init_windfarm_x_y (const int lev,
                            const std::string windfarm_loc_table,
                            Vector<Real>& xloc,
                            Vector<Real>& yloc)
{
    // Read turbine locations from windturbines.txt
    std::ifstream file(windfarm_loc_table);
    if (!file.is_open()) {
        amrex::Error("Wind turbines location table not found. Either the inputs is missing the"
                     " erf.windfarm_loc_table entry or the file specified in the entry " + windfarm_loc_table + " is missing.");
    }
    // Vector of vectors to store the matrix
    Real value1, value2;

    while (file >> value1 >> value2) {
        xloc.push_back(value1);
        yloc.push_back(value2);
    }
    file.close();
}
