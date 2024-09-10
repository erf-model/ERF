/**
 * \file InitWindFarm.cpp
 */

#include <WindFarm.H>

using namespace amrex;

/**
 * Read in the turbine locations in latitude-longitude from windturbines.txt
 * and convert it into x and y coordinates in metres
 *
 * @param lev Integer specifying the current level
 */
void
WindFarm::read_tables (std::string windfarm_loc_table,
                         std::string windfarm_spec_table,
                         bool x_y, bool lat_lon,
                         const Real latitude_lo,
                         const Real longitude_lo)
{
    amrex::Print() << "Reading wind turbine locations table" << "\n";
    read_windfarm_locations_table(windfarm_loc_table,
                                   x_y, lat_lon,
                                  latitude_lo,longitude_lo);

    amrex::Print() << "Reading wind turbine specifications table" << "\n";
    read_windfarm_spec_table(windfarm_spec_table);
}

void
WindFarm::read_windfarm_locations_table(const std::string windfarm_loc_table,
                                           bool x_y, bool lat_lon,
                                           const Real latitude_lo,
                                           const Real longitude_lo)
{
    if(x_y) {
        init_windfarm_x_y(windfarm_loc_table);
    }
    else if(lat_lon) {
        init_windfarm_lat_lon(windfarm_loc_table, latitude_lo, longitude_lo);
    }
    else {
        amrex::Abort("Are you using windfarms? For windfarm simulations, the inputs need to have an"
                     " entry erf.windfarm_loc_table which should not be None. \n");
    }

    set_turb_loc(xloc, yloc);
}

void
WindFarm::init_windfarm_lat_lon (const std::string windfarm_loc_table,
                                 const Real latitude_lo, const Real longitude_lo)
{

    if(latitude_lo == -1e10) {
        amrex::Error("Are you using latitude-longitude for wind turbine locations? There needs to be"
                     " entry for erf.latitude_lo in the inputs for"
                     " the lower bottom corner of the domain");
    }

    if(longitude_lo == -1e10) {
        amrex::Error("Are you using latitude-longitude for wind turbine locations? There needs to be"
                     " entry erf.longitude_lo in the inputs for"
                     " the lower bottom corner of the domain");
    }

    if(std::fabs(latitude_lo) > 90.0) {
        amrex::Error("The value of erf.latitude_lo in the input file should be within -90 and 90");
    }

    if(std::fabs(longitude_lo) > 180.0) {
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

        if(std::fabs(longitude_lo) > 180.0) {
            amrex::Error("The values of erf.longitude_lo should be within -180 and 180");
        }
        lat.push_back(value1);
        lon.push_back(value2);
    }
    file.close();

    Real rad_earth = 6371.0e3; // Radius of the earth

    Real lat_lo  = latitude_lo*M_PI/180.0;
    Real lon_lo  = longitude_lo*M_PI/180.0;

    // Find the coordinates of average of min and max of the farm
    // Rotate about that point

    Real lat_min = *std::min_element(lat.begin(), lat.end());
    Real lat_max = *std::max_element(lat.begin(), lat.end());
    Real lon_min = *std::min_element(lon.begin(), lon.end());
    Real lon_max = *std::max_element(lon.begin(), lon.end());

    Real lat_cen = 0.5*(lat_min+lat_max)*M_PI/180.0;
    Real lon_cen = 0.5*(lon_min+lon_max)*M_PI/180.0;

    // (lat_lo, lon_lo) is mapped to (0,0)


    for(int it=0;it<lat.size();it++){
        lat[it] = lat[it]*M_PI/180.0;
        lon[it] = lon[it]*M_PI/180.0;
        Real delta_lat = (lat[it] - lat_cen);
        Real delta_lon = (lon[it] - lon_cen);

        Real term1 = std::pow(sin(delta_lat/2.0),2);
        Real term2 = cos(lat[it])*cos(lat_cen)*std::pow(sin(delta_lon/2.0),2);
        Real dist =  2.0*rad_earth*std::asin(std::sqrt(term1 + term2));
        Real dy_turb = (lat[it] - lat_cen) * 111000.0 * 180.0/M_PI ;
        Real dx_turb = std::sqrt(std::pow(dist,2) - std::pow(dy_turb,2));
        if(delta_lon >= 0.0) {
            xloc.push_back(dx_turb);
        }
        else {
            xloc.push_back(-dx_turb);
        }
        yloc.push_back(dy_turb);
    }

    Real xloc_min = *std::min_element(xloc.begin(),xloc.end());
    Real yloc_min = *std::min_element(yloc.begin(),yloc.end());

    for(int it = 0;it<xloc.size(); it++){
        xloc[it] = xloc[it] - xloc_min + 1000;
        yloc[it] = yloc[it] - yloc_min + 1000;
    }

    Real xmin = *std::min_element(xloc.begin(), xloc.end());
    Real xmax = *std::max_element(xloc.begin(), xloc.end());
    Real ymin = *std::min_element(yloc.begin(), yloc.end());
    Real ymax = *std::max_element(yloc.begin(), yloc.end());

    Real xcen = 0.5*(xmin+xmax);
    Real ycen = 0.5*(ymin+ymax);

    Real theta = 0.25*M_PI;

    for(int it = 0;it<xloc.size(); it++){
        Real xnew = ( (xloc[it]-xcen)*std::cos(theta) + (yloc[it]-ycen)*std::sin(theta));
        Real ynew = (-(xloc[it]-xcen)*std::sin(theta) + (yloc[it]-ycen)*std::cos(theta));

        xloc[it] = xcen + xnew;
        yloc[it] = ycen + ynew;
    }


    Real xloc_min1 = *std::min_element(xloc.begin(),xloc.end());
    Real yloc_min1 = *std::min_element(yloc.begin(),yloc.end());

    for(int it = 0;it<xloc.size(); it++){
        xloc[it] = xloc[it] - xloc_min1 + 1000;
        yloc[it] = yloc[it] - yloc_min1 + 1000;
    }
}

void
WindFarm::init_windfarm_x_y (const std::string windfarm_loc_table)
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


void
WindFarm::read_windfarm_spec_table(const std::string windfarm_spec_table)
{
    //The first line is the number of pairs entries for the power curve and thrust coefficient.
    //The second line gives first the height in meters of the turbine hub, second, the diameter in
    //meters of the rotor, third the standing thrust coefficient, and fourth the nominal power of
    //the turbine in MW.
    //The remaining lines contain the three values of: wind speed, thrust coefficient, and power production in kW.

     // Read turbine data from wind-turbine-1.tbl
    std::ifstream file_turb_table(windfarm_spec_table);
    if (!file_turb_table.is_open()) {
        Error("Wind farm specifications table not found. Either the inputs is missing the "
                      "erf.windfarm_spec_table entry or the file specified in the entry - " + windfarm_spec_table + " is missing.");
    }
    else {
        Print() << "Reading in wind farm specifications table: " << windfarm_spec_table << "\n";
    }

    int nlines;
    file_turb_table >> nlines;
    wind_speed.resize(nlines);
    thrust_coeff.resize(nlines);
    power.resize(nlines);

    Real rotor_dia;
    file_turb_table >> hub_height >> rotor_dia >> thrust_coeff_standing >> nominal_power;
    rotor_rad = rotor_dia*0.5;
    if(rotor_rad > hub_height) {
        Abort("The blade length is more than the hub height. Check the second line in wind-turbine-1.tbl. Aborting.....");
    }
    if(thrust_coeff_standing > 1.0) {
        Abort("The standing thrust coefficient is greater than 1. Check the second line in wind-turbine-1.tbl. Aborting.....");
    }

    for(int iline=0;iline<nlines;iline++){
        file_turb_table >> wind_speed[iline] >> thrust_coeff[iline] >> power[iline];
        if(thrust_coeff[iline] > 1.0) {
            Abort("The thrust coefficient is greater than 1. Check wind-turbine-1.tbl. Aborting.....");
        }
    }
    file_turb_table.close();

    set_turb_spec(rotor_rad, hub_height, thrust_coeff_standing,
                  wind_speed, thrust_coeff, power);

}

void
WindFarm::fill_Nturb_multifab(const Geometry& geom,
                              MultiFab& mf_Nturb)
{

    amrex::Gpu::DeviceVector<Real> d_xloc(xloc.size());
    amrex::Gpu::DeviceVector<Real> d_yloc(yloc.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xloc.begin(), xloc.end(), d_xloc.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yloc.begin(), yloc.end(), d_yloc.begin());

    Real* d_xloc_ptr     = d_xloc.data();
    Real* d_yloc_ptr     = d_yloc.data();

    mf_Nturb.setVal(0);

    int i_lo = geom.Domain().smallEnd(0); int i_hi = geom.Domain().bigEnd(0);
    int j_lo = geom.Domain().smallEnd(1); int j_hi = geom.Domain().bigEnd(1);
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    int num_turb = xloc.size();

     // Initialize wind farm
    for ( MFIter mfi(mf_Nturb,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx     = mfi.tilebox();
        auto  Nturb_array = mf_Nturb.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int li = amrex::min(amrex::max(i, i_lo), i_hi);
            int lj = amrex::min(amrex::max(j, j_lo), j_hi);

            Real x1 = ProbLoArr[0] + li*dx[0];
            Real x2 = ProbLoArr[0] + (li+1)*dx[0];
            Real y1 = ProbLoArr[1] + lj*dx[1];
            Real y2 = ProbLoArr[1] + (lj+1)*dx[1];

            for(int it=0; it<num_turb; it++){
                if( d_xloc_ptr[it]+1e-12 > x1 and d_xloc_ptr[it]+1e-12 < x2 and
                    d_yloc_ptr[it]+1e-12 > y1 and d_yloc_ptr[it]+1e-12 < y2){
                       Nturb_array(i,j,k,0) = Nturb_array(i,j,k,0) + 1;
                }
            }
        });
    }
}


void
WindFarm::fill_SMark_multifab(const Geometry& geom,
                              MultiFab& mf_SMark,
                              const Real& sampling_distance_by_D)
{
    amrex::Gpu::DeviceVector<Real> d_xloc(xloc.size());
    amrex::Gpu::DeviceVector<Real> d_yloc(yloc.size());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xloc.begin(), xloc.end(), d_xloc.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yloc.begin(), yloc.end(), d_yloc.begin());

    Real d_rotor_rad = rotor_rad;
    Real d_hub_height = hub_height;
    Real d_sampling_distance = sampling_distance_by_D*2.0*rotor_rad;

    Real* d_xloc_ptr     = d_xloc.data();
    Real* d_yloc_ptr     = d_yloc.data();

    mf_SMark.setVal(-1.0);

    int i_lo = geom.Domain().smallEnd(0); int i_hi = geom.Domain().bigEnd(0);
    int j_lo = geom.Domain().smallEnd(1); int j_hi = geom.Domain().bigEnd(1);
    int k_lo = geom.Domain().smallEnd(2); int k_hi = geom.Domain().bigEnd(2);
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    int num_turb = xloc.size();

     // Initialize wind farm
    for ( MFIter mfi(mf_SMark,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx     = mfi.tilebox();
        auto  SMark_array = mf_SMark.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = amrex::min(amrex::max(i, i_lo), i_hi);
            int jj = amrex::min(amrex::max(j, j_lo), j_hi);
            int kk = amrex::min(amrex::max(k, k_lo), k_hi);

            Real x1 = ProbLoArr[0] + ii*dx[0];
            Real x2 = ProbLoArr[0] + (ii+1)*dx[0];

            Real y = ProbLoArr[1] + (jj+0.5) * dx[1];
            Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

            for(int it=0; it<num_turb; it++){
                if(d_xloc_ptr[it]-d_sampling_distance+1e-12 > x1 and d_xloc_ptr[it]-d_sampling_distance+1e-12 < x2) {
                   if(std::pow((y-d_yloc_ptr[it])*(y-d_yloc_ptr[it]) + (z-d_hub_height)*(z-d_hub_height),0.5) < d_rotor_rad) {
                       SMark_array(i,j,k,0) = it;
                    }
                }
            }
        });
    }
}

void
WindFarm::write_turbine_locations_vtk()
{
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
}


void
WindFarm::write_actuator_disks_vtk(const Geometry& geom)
{

    if (ParallelDescriptor::IOProcessor()){
        FILE *file_actuator_disks_all, *file_actuator_disks_in_dom;
        file_actuator_disks_all = fopen("actuator_disks_all.vtk","w");
        fprintf(file_actuator_disks_all, "%s\n","# vtk DataFile Version 3.0");
        fprintf(file_actuator_disks_all, "%s\n","Actuator Disks");
        fprintf(file_actuator_disks_all, "%s\n","ASCII");
        fprintf(file_actuator_disks_all, "%s\n","DATASET POLYDATA");

        file_actuator_disks_in_dom = fopen("actuator_disks_in_dom.vtk","w");
        fprintf(file_actuator_disks_in_dom, "%s\n","# vtk DataFile Version 3.0");
        fprintf(file_actuator_disks_in_dom, "%s\n","Actuator Disks");
        fprintf(file_actuator_disks_in_dom, "%s\n","ASCII");
        fprintf(file_actuator_disks_in_dom, "%s\n","DATASET POLYDATA");

        int npts = 100;
        fprintf(file_actuator_disks_all, "%s %ld %s\n", "POINTS", xloc.size()*npts, "float");
        auto ProbLoArr = geom.ProbLoArray();
        auto ProbHiArr = geom.ProbHiArray();
        int num_turb_in_dom = 0;

        // Find the number of turbines inside the specified computational domain

        for(int it=0; it<xloc.size(); it++){
            Real x = xloc[it];
            Real y = yloc[it];
            if(x > ProbLoArr[0] and x < ProbHiArr[0] and y > ProbLoArr[1] and y < ProbHiArr[1]) {
                num_turb_in_dom++;
            }
        }
        fprintf(file_actuator_disks_in_dom, "%s %ld %s\n", "POINTS", num_turb_in_dom*npts, "float");

        for(int it=0; it<xloc.size(); it++){
            Real x;
            x = xloc[it]+1e-12;
            for(int pt=0;pt<100;pt++){
                Real y, z;
                Real theta = 2.0*M_PI/npts*pt;
                y = yloc[it]+rotor_rad*cos(theta);
                z = hub_height+rotor_rad*sin(theta);
                fprintf(file_actuator_disks_all, "%0.15g %0.15g %0.15g\n", x, y, z);
                if(xloc[it] > ProbLoArr[0] and xloc[it] < ProbHiArr[0] and yloc[it] > ProbLoArr[1] and yloc[it] < ProbHiArr[1]) {
                    fprintf(file_actuator_disks_in_dom, "%0.15g %0.15g %0.15g\n", x, y, z);
                }
            }
        }
        fprintf(file_actuator_disks_all, "%s %ld %ld\n", "LINES", xloc.size()*(npts-1), static_cast<long int>(xloc.size()*(npts-1)*3));
        fprintf(file_actuator_disks_in_dom, "%s %ld %ld\n", "LINES", num_turb_in_dom*(npts-1), static_cast<long int>(num_turb_in_dom*(npts-1)*3));
        for(int it=0; it<xloc.size(); it++){
            for(int pt=0;pt<99;pt++){
                fprintf(file_actuator_disks_all, "%ld %ld %ld\n",
                                             static_cast<long int>(2),
                                             static_cast<long int>(it*npts+pt),
                                             static_cast<long int>(it*npts+pt+1));
            }
        }
         for(int it=0; it<num_turb_in_dom; it++){
            for(int pt=0;pt<99;pt++){
                fprintf(file_actuator_disks_in_dom, "%ld %ld %ld\n",
                                             static_cast<long int>(2),
                                             static_cast<long int>(it*npts+pt),
                                             static_cast<long int>(it*npts+pt+1));
            }
        }

        fclose(file_actuator_disks_all);
        fclose(file_actuator_disks_in_dom);
    }
}


