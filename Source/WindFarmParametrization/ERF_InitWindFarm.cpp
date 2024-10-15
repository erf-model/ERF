/**
 * \file ERF_InitWindFarm.cpp
 */

#include <ERF_WindFarm.H>
#include <filesystem>
#include <dirent.h>   // For POSIX directory handling

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
                       const Real windfarm_x_shift,
                       const Real windfarm_y_shift)
{
    amrex::Print() << "Reading wind turbine locations table" << "\n";
    read_windfarm_locations_table(windfarm_loc_table,
                                  x_y, lat_lon,
                                  windfarm_x_shift, windfarm_y_shift);

    amrex::Print() << "Reading wind turbine specifications table" << "\n";
    read_windfarm_spec_table(windfarm_spec_table);
}

void
WindFarm::read_windfarm_locations_table(const std::string windfarm_loc_table,
                                        bool x_y, bool lat_lon,
                                        const Real windfarm_x_shift,
                                        const Real windfarm_y_shift)
{
    if(x_y) {
        init_windfarm_x_y(windfarm_loc_table);
    }
    else if(lat_lon) {
        init_windfarm_lat_lon(windfarm_loc_table, windfarm_x_shift, windfarm_y_shift);
    }
    else {
        amrex::Abort("Are you using windfarms? For windfarm simulations, the inputs need to have an"
                     " entry erf.windfarm_loc_type which should be either lat_lon or x_y. \n");
    }

    set_turb_loc(xloc, yloc);
}

void
WindFarm::init_windfarm_lat_lon (const std::string windfarm_loc_table,
                                 const Real windfarm_x_shift,
                                 const Real windfarm_y_shift)
{

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
        lat.push_back(value1);
        lon.push_back(value2);
    }
    file.close();

    Real rad_earth = 6371.0e3; // Radius of the earth

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
        xloc[it] = xloc[it] - xloc_min + windfarm_x_shift;
        yloc[it] = yloc[it] - yloc_min + windfarm_y_shift;
    }

    FILE* file_xy_loc;
    file_xy_loc = fopen("file_xy_loc_KingPlains.txt","w");

    for(int it = 0;it<xloc.size(); it++){
        fprintf(file_xy_loc,"%0.15g %0.15g %0.15g\n", xloc[it], yloc[it], 89.0);
    }
    fclose(file_xy_loc);
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
        value1 = value1 + 1e-3;
        value2 = value2 + 1e-3;
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
WindFarm::read_windfarm_blade_table(const std::string windfarm_blade_table)
{
    std::ifstream filename(windfarm_blade_table);
    std::string line;
    Real temp, var1, var2, var3;
    if (!filename.is_open()) {
        Error("You are using a generalized actuator disk model based on blade element theory. This needs info of blades."
                      " An entry erf.windfarm_blade_table is needed. Either the entry is missing or the file specified"
                      " in the entry - " + windfarm_blade_table + " is missing.");
    }
    else {
        Print() << "Reading in wind farm blade table: " << windfarm_blade_table << "\n";

        // First 6 lines are comments

        for (int i = 0; i < 6; ++i) {
            if (std::getline(filename, line)) {  // Read one line into the array
            }
        }

        while(filename >> var1 >> temp >> temp >> temp >> var2 >> var3 >> temp) {
            bld_rad_loc.push_back(var1);
            bld_twist.push_back(var2);
            bld_chord.push_back(var3);
            //int idx = bld_rad_loc.size()-1;
            //printf("Values are = %0.15g %0.15g %0.15g\n", bld_rad_loc[idx], bld_twist[idx], bld_chord[idx]);
        }
        set_blade_spec(bld_rad_loc, bld_twist, bld_chord);
        n_bld_sections = bld_rad_loc.size();
    }
}

void
WindFarm::read_windfarm_spec_table_extra(const std::string windfarm_spec_table_extra)
{
    // Open the file
    std::ifstream file(windfarm_spec_table_extra);

    // Check if file opened successfully
    if (!file.is_open()) {
        Abort("Error: You are using generalized wind farms option. This requires an input file erf.windfarm_spec_table_extra."
              " Either this entry is missing in the inputs or the file specified -" + windfarm_spec_table_extra + " does"
              " not exist. Exiting...");
    } else {
        printf("Reading in windfarm_spec_table_extra %s", windfarm_spec_table_extra.c_str());
    }

    // Ignore the first line (header)
    std::string header;
    std::getline(file, header);

    // Variables to hold each row's values
    double V, Cp, Ct, rpm, pitch, temp;

    // Read the file row by row
    while (file >> V) {
        char comma;  // To ignore the commas
        file >> comma >> Cp >> comma >> Ct >> comma >> temp >> comma >> temp >> comma
             >> temp >> comma >> rpm >> comma >> pitch >> comma >> temp;

        velocity.push_back(V);
        C_P.push_back(Cp);
        C_T.push_back(Ct);
        rotor_RPM.push_back(rpm);
        blade_pitch.push_back(pitch);
    }

    set_turb_spec_extra(velocity, C_P, C_T, rotor_RPM, blade_pitch);
}


void
WindFarm::read_windfarm_airfoil_tables(const std::string windfarm_airfoil_tables,
                                       const std::string windfarm_blade_table)
{
    DIR* dir;
    struct dirent* entry;
    std::vector<std::string> files;

    // Check if directory exists
    if ((dir = opendir(windfarm_airfoil_tables.c_str())) == nullptr) {
       Abort("You are using a generalized actuator disk model based on blade element theory. This needs info of airfoil"
             " cross sections over the span of the blade. There needs to be an entry erf.airfoil_tables which is the directory that"
             " contains the angle of attack, Cl, Cd data for each airfoil cross-section. Either the entry is missing or the directory specified"
             " in the entry - " + windfarm_airfoil_tables + " is missing. Exiting...");
    }

    // Loop through directory entries and collect filenames
    while ((entry = readdir(dir)) != nullptr) {
        // Skip special directory entries "." and ".."
        if (std::string(entry->d_name) == "." || std::string(entry->d_name) == "..") {
            continue;
        }
        files.emplace_back(windfarm_airfoil_tables + "/" + entry->d_name);  // Add file path to vector
    }

    // Close the directory
    closedir(dir);

    if (files.empty()) {
        Abort("It seems the directory containing the info of airfoil cross sections of the blades - " + windfarm_airfoil_tables +
              " is empty. Exiting...");
    }

    if(files.size() != static_cast<long double>(n_bld_sections)) {
        printf("There are %d airfoil sections in the last column of %s. But the number"
               " of files in %s is only %ld.\n", n_bld_sections, windfarm_blade_table.c_str(),
                windfarm_airfoil_tables.c_str(), files.size());
        Abort("The number of blade sections from " + windfarm_blade_table + " should match the number of"
              " files in " + windfarm_airfoil_tables + ". Exiting...");
    }

    // Sort filenames in lexicographical (alphabetical) order
    std::sort(files.begin(), files.end());

    // Process each file
    int count = 0;
    bld_airfoil_aoa.resize(n_bld_sections);
    bld_airfoil_Cl.resize(n_bld_sections);
    bld_airfoil_Cd.resize(n_bld_sections);
    for (const auto& filePath : files) {
        std::ifstream filename(filePath.c_str());

        if (!filename.is_open()) {
            std::cerr << "Failed to open file: " << filePath << std::endl;
            continue;  // Move on to the next file
        }

           std::cout << "Reading file: " << filePath << std::endl;

        std::string line;
        for (int i = 0; i < 54; ++i) {
            if (std::getline(filename, line)) {  // Read one line into the array
            }
        }

        Real var1, var2, var3, temp;

        while(filename >> var1 >> var2 >> var3 >> temp) {
            bld_airfoil_aoa[count].push_back(var1);
            bld_airfoil_Cl[count].push_back(var2);
            bld_airfoil_Cd[count].push_back(var3);
            //int idx = bld_airfoil_aoa.size()-1;
            //printf("Values are = %0.15g %0.15g %0.15g\n", bld_airfoil_aoa[idx], bld_airfoil_Cl[idx], bld_airfoil_Cd[idx]);
        }
        count++;
    }

    set_blade_airfoil_spec(bld_airfoil_aoa, bld_airfoil_Cl, bld_airfoil_Cd);
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
    if(dx[0]<= 1e-3 or dx[1]<=1e-3 or dx[2]<= 1e-3) {
        Abort("The value of mesh spacing for wind farm parametrization cannot be less than 1e-3 m. "
              "It should be usually of order 1 m");
    }
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
                if( d_xloc_ptr[it]+1e-3 > x1 and d_xloc_ptr[it]+1e-3 < x2 and
                    d_yloc_ptr[it]+1e-3 > y1 and d_yloc_ptr[it]+1e-3 < y2){
                       Nturb_array(i,j,k,0) = Nturb_array(i,j,k,0) + 1;
                }
            }
        });
    }
}

void
WindFarm::fill_SMark_multifab(const Geometry& geom,
                              MultiFab& mf_SMark,
                              const Real& sampling_distance_by_D,
                              const Real& turb_disk_angle)
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

    Real theta = turb_disk_angle*M_PI/180.0-0.5*M_PI;

    set_turb_disk_angle(theta);
    my_turb_disk_angle = theta;

    Real nx = -std::cos(theta);
    Real ny = -std::sin(theta);

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
            Real y1 = ProbLoArr[1] + jj*dx[1];
            Real y2 = ProbLoArr[1] + (jj+1)*dx[1];

            Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

            for(int it=0; it<num_turb; it++){
                Real x0 = d_xloc_ptr[it] + d_sampling_distance*nx;
                Real y0 = d_yloc_ptr[it] + d_sampling_distance*ny;

                bool is_cell_marked = find_if_marked(x1, x2, y1, y2, x0, y0,
                                                     nx, ny, d_hub_height, d_rotor_rad, z);
                if(is_cell_marked) {
                    SMark_array(i,j,k,0) = it;
                }
                x0 = d_xloc_ptr[it];
                y0 = d_yloc_ptr[it];

                is_cell_marked = find_if_marked(x1, x2, y1, y2, x0, y0,
                                                nx, ny, d_hub_height, d_rotor_rad, z);
                if(is_cell_marked) {
                    SMark_array(i,j,k,1) = it;
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
            fprintf(file_turbloc_vtk, "%0.15g %0.15g %0.15g\n", xloc[it], yloc[it], hub_height);
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
        fprintf(file_actuator_disks_in_dom, "%s %ld %s\n", "POINTS", static_cast<long int>(num_turb_in_dom*npts), "float");

        Real nx = std::cos(my_turb_disk_angle+0.5*M_PI);
        Real ny = std::sin(my_turb_disk_angle+0.5*M_PI);

        for(int it=0; it<xloc.size(); it++){
            for(int pt=0;pt<100;pt++){
                Real x, y, z;
                Real theta = 2.0*M_PI/npts*pt;
                x = xloc[it] + rotor_rad*cos(theta)*nx;
                y = yloc[it] + rotor_rad*cos(theta)*ny;
                z = hub_height + rotor_rad*sin(theta);
                fprintf(file_actuator_disks_all, "%0.15g %0.15g %0.15g\n", x, y, z);
                if(xloc[it] > ProbLoArr[0] and xloc[it] < ProbHiArr[0] and yloc[it] > ProbLoArr[1] and yloc[it] < ProbHiArr[1]) {
                    fprintf(file_actuator_disks_in_dom, "%0.15g %0.15g %0.15g\n", x, y, z);
                }
            }
        }
        fprintf(file_actuator_disks_all, "%s %ld %ld\n", "LINES", xloc.size()*(npts-1), static_cast<long int>(xloc.size()*(npts-1)*3));
        fprintf(file_actuator_disks_in_dom, "%s %ld %ld\n", "LINES", static_cast<long int>(num_turb_in_dom*(npts-1)), static_cast<long int>(num_turb_in_dom*(npts-1)*3));
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


