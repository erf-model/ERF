/**
 * \file ERF_init_windfarm.cpp
 */

#include <ERF.H>

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
    // Read turbine locations from windturbines.txt
    FILE *file_windturbines;
    file_windturbines = fopen("windturbines.txt","r");
    std::ifstream file("windturbines.txt");
    if (!file.is_open()) {
        amrex::Error("Wind turbines location file windturbines.txt not found");
    }
    // Vector of vectors to store the matrix
    std::vector<Real> lat, lon, xloc, yloc;
    Real value1, value2, value3;
    while (file >> value1 >> value2 >> value3) {
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

    Nturb[lev].setVal(0);

    int i_lo = geom[lev].Domain().smallEnd(0); int i_hi = geom[lev].Domain().bigEnd(0);
    int j_lo = geom[lev].Domain().smallEnd(1); int j_hi = geom[lev].Domain().bigEnd(1);

     // Initialize wind farm
    for ( MFIter mfi(Nturb[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx     = mfi.tilebox();
        auto  Nturb_array = Nturb[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int li = amrex::min(amrex::max(i, i_lo), i_hi);
            int lj = amrex::min(amrex::max(j, j_lo), j_hi);

            auto dx = geom[lev].CellSizeArray();
            Real x1 = li*dx[0], x2 = (li+1)*dx[0];
            Real y1 = lj*dx[1], y2 = (lj+1)*dx[1];

            for(int it=0; it<xloc.size(); it++){
                if( xloc[it]+1e-12 > x1 and xloc[it]+1e-12 < x2 and
                    yloc[it]+1e-12 > y1 and yloc[it]+1e-12 < y2){
                    Nturb_array(i,j,k,0) = Nturb_array(i,j,k,0) + 1;
                }
            }
        });
    }
}

