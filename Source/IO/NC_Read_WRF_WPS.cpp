#include "ERF.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
ERF::read_from_wrfinput(int lev)
{
    amrex::Print() << "Loading initial data from NetCDF file " << std::endl;
    Box input_box;

    NC_xvel_fab[lev].resize(num_boxes_at_level[lev]);
    NC_yvel_fab[lev].resize(num_boxes_at_level[lev]);
    NC_zvel_fab[lev].resize(num_boxes_at_level[lev]);
    NC_rho_fab[lev].resize(num_boxes_at_level[lev]);
    NC_rhotheta_fab[lev].resize(num_boxes_at_level[lev]);

#ifdef ERF_USE_TERRAIN
    NC_PH_fab[lev].resize(num_boxes_at_level[lev]);
    NC_PHB_fab[lev].resize(num_boxes_at_level[lev]);
#endif

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        if (lev == 0) {
            input_box = geom[0].Domain();
        } else {
            input_box = boxes_at_level[lev][idx];
        }

        // We allocate these here so they exist on all ranks
        Box ubx(input_box); ubx.surroundingNodes(0);
        Box vbx(input_box); vbx.surroundingNodes(1);
        Box wbx(input_box); wbx.surroundingNodes(2);

        NC_xvel_fab[lev][idx].resize(ubx,1);
        NC_yvel_fab[lev][idx].resize(vbx,1);
        NC_zvel_fab[lev][idx].resize(wbx,1);
        NC_rho_fab[lev][idx].resize(input_box,1);
        NC_rhotheta_fab[lev][idx].resize(input_box,1);

#ifdef ERF_USE_TERRAIN
        NC_PH_fab[lev][idx].resize(wbx,1);
        NC_PHB_fab[lev][idx].resize(wbx,1);
#endif

#ifdef AMREX_USE_GPU
        FArrayBox host_NC_xvel_fab   (NC_xvel_fab[lev][idx].box(),     NC_xvel_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_yvel_fab   (NC_yvel_fab[lev][idx].box(),     NC_yvel_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_zvel_fab   (NC_zvel_fab[lev][idx].box(),     NC_zvel_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_rho_fab    (NC_rho_fab[lev][idx].box(),      NC_rho_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
       FArrayBox host_NC_rhotheta_fab(NC_rhotheta_fab[lev][idx].box(), NC_rhotheta_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
#ifdef ERF_USE_TERRAIN
        FArrayBox host_NC_PH_fab (NC_PH_fab[lev][idx].box(),  NC_PH_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_PHB_fab(NC_PHB_fab[lev][idx].box(), NC_PHB_fab[lev][idx].nComp(),amrex::The_Pinned_Arena());
#endif
#else
        FArrayBox host_NC_xvel_fab    (NC_xvel_fab[lev][idx]    , amrex::make_alias, 0, NC_xvel_fab[lev][idx].nComp());
        FArrayBox host_NC_yvel_fab    (NC_yvel_fab[lev][idx]    , amrex::make_alias, 0, NC_yvel_fab[lev][idx].nComp());
        FArrayBox host_NC_zvel_fab    (NC_zvel_fab[lev][idx]    , amrex::make_alias, 0, NC_zvel_fab[lev][idx].nComp());
        FArrayBox host_NC_rho_fab     (NC_rho_fab[lev][idx]     , amrex::make_alias, 0, NC_rho_fab[lev][idx].nComp());
        FArrayBox host_NC_rhotheta_fab(NC_rhotheta_fab[lev][idx], amrex::make_alias, 0, NC_rhotheta_fab[lev][idx].nComp());
#ifdef ERF_USE_TERRAIN
        FArrayBox host_NC_PH_fab      (NC_PH_fab[lev][idx]      , amrex::make_alias, 0, NC_PH_fab[lev][idx].nComp());
        FArrayBox host_NC_PHB_fab     (NC_PHB_fab[lev][idx]     , amrex::make_alias, 0, NC_PHB_fab[lev][idx].nComp());
#endif
#endif

        if (ParallelDescriptor::IOProcessor())
        {
            Vector<FArrayBox*> NC_fabs;
            Vector<std::string> NC_names;
            Vector<enum NC_Data_Dims_Type> NC_dim_types;

            NC_fabs.push_back(&host_NC_xvel_fab);     NC_names.push_back("U");
            NC_fabs.push_back(&host_NC_yvel_fab);     NC_names.push_back("V");
            NC_fabs.push_back(&host_NC_zvel_fab);     NC_names.push_back("W");
            NC_fabs.push_back(&host_NC_rho_fab);      NC_names.push_back("ALB");
            NC_fabs.push_back(&host_NC_rhotheta_fab); NC_names.push_back("T_INIT");

            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);

#ifdef ERF_USE_TERRAIN
            NC_fabs.push_back(&host_NC_PH_fab);  NC_names.push_back("PH");
            NC_fabs.push_back(&host_NC_PHB_fab); NC_names.push_back("PHB");

            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
#endif

            // Read the netcdf file and fill these FABs
            // NOTE: right now we are hard-wired to one "domain" per level -- but that can be generalized
            //       once we know how to determine the level for each input file
            BuildFABsFromWRFInputFile(nc_init_file[lev][idx], NC_names, NC_fabs, NC_dim_types);

        } // if ParalleDescriptor::IOProcessor()

        // We put a barrier here so the rest of the processors wait to do anything until they have the data
        amrex::ParallelDescriptor::Barrier();

        // When an FArrayBox is built, space is allocated on every rank.  However, we only
        //    filled the data in these FABs on the IOProcessor.  So here we broadcast
        //    the data to every rank.

        int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
        ParallelDescriptor::Bcast(host_NC_xvel_fab.dataPtr(),NC_xvel_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_yvel_fab.dataPtr(),NC_yvel_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_zvel_fab.dataPtr(),NC_zvel_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_rho_fab.dataPtr(),NC_rho_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_rhotheta_fab.dataPtr(),NC_rhotheta_fab[lev][idx].box().numPts(),ioproc);
#ifdef ERF_USE_TERRAIN
        ParallelDescriptor::Bcast(host_NC_PHB_fab.dataPtr(),NC_PHB_fab[lev][idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_PH_fab.dataPtr() ,NC_PH_fab[lev][idx].box().numPts() ,ioproc);
#endif

#ifdef AMREX_USE_GPU
         Gpu::copy(Gpu::hostToDevice, host_NC_xvel_fab.dataPtr(), host_NC_xvel_fab.dataPtr()+host_NC_xvel_fab.size(),
                                           NC_xvel_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_yvel_fab.dataPtr(), host_NC_yvel_fab.dataPtr()+host_NC_yvel_fab.size(),
                                           NC_yvel_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_zvel_fab.dataPtr(), host_NC_zvel_fab.dataPtr()+host_NC_zvel_fab.size(),
                                           NC_zvel_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rho_fab.dataPtr(), host_NC_rho_fab.dataPtr()+host_NC_rho_fab.size(),
                                           NC_rho_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rhotheta_fab.dataPtr(), host_NC_rhotheta_fab.dataPtr()+host_NC_rhotheta_fab.size(),
                                           NC_rhotheta_fab[lev][idx].dataPtr());
#ifdef ERF_USE_TERRAIN
         Gpu::copy(Gpu::hostToDevice, host_NC_PH_fab.dataPtr(), host_NC_PH_fab.dataPtr()+host_NC_PH_fab.size(),
                                           NC_PH_fab[lev][idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_PHB_fab.dataPtr(), host_NC_PHB_fab.dataPtr()+host_NC_PHB_fab.size(),
                                           NC_PHB_fab[lev][idx].dataPtr());
#endif
#endif

        // Convert to rho by inverting
        NC_rho_fab[lev][idx].template invert<RunOn::Device>(1.0);

        // The ideal.exe NetCDF file has this ref value subtracted from theta or T_INIT. Need to add in ERF.
        const Real theta_ref = 300.0;
        NC_rhotheta_fab[lev][idx].template plus<RunOn::Device>(theta_ref);

        // Now multiply by rho to get (rho theta) instead of theta
        NC_rhotheta_fab[lev][idx].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);

        amrex::Print() <<
          "Successfully loaded data from the wrfinput (output of 'ideal.exe' / 'real.exe') NetCDF file at level " << lev << std::endl;
    } // idx
}

void
ERF::read_from_wrfbdy()
{
    amrex::Print() << "Loading boundary data from NetCDF file " << std::endl;
    // *********************************************************
    // Allocate space for all of the boundary planes we may need
    // Here we make only one enough space for one time -- we will
    //    add space for the later time slices later
    // *********************************************************

    const amrex::Box& domain = geom[0].Domain();

    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

    amrex::IntVect plo(lo);
    amrex::IntVect phi(hi);

    /*
     bdy_data_xlo[time] (etc) contain 4 different variables U_BXS, V_BXS, W_BXS, T_BXS.
     The dimensions of these on xlo and xhi are:
     U_BXS -> dim -> 1     NY       NZ
     V_BXS -> dim -> 1    (NY+1)    NZ
     W_BXS -> dim -> 1     NY      (NZ+1)
     T_BXS -> dim -> 1     NY       NZ

     The dimensions of these on ylo and yhi are:
     U_BYS -> dim -> 1    (NX+1)    NZ
     V_BYS -> dim -> 1     NX       NZ
     W_BYS -> dim -> 1     NX      (NZ+1)
     T_BYS -> dim -> 1     NX       NZ
    */

    // *******************************************************************************
    // First allocate space for just one time
    // *******************************************************************************
    bdy_data_xlo.resize(1);
    bdy_data_xhi.resize(1);
    bdy_data_ylo.resize(1);
    bdy_data_yhi.resize(1);

    // *******************************************************************************
    // xlo bdy
    // *******************************************************************************
    plo[0] = -1; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = -1; phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xlo(plo, phi);

    Box xlo_plane_no_stag(pbx_xlo);
    Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});
    Box xlo_plane_z_stag = convert(pbx_xlo, {0, 0, 1});

    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_no_stag, 1)); // U
    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_y_stag , 1)); // V
    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_z_stag , 1)); // W
    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_no_stag, 1)); // T

    // *******************************************************************************
    // xhi bdy
    // *******************************************************************************
    plo[0] = hi[0] + 1; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0] + 1; phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xhi(plo, phi);

    Box xhi_plane_no_stag(pbx_xhi);
    Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});
    Box xhi_plane_z_stag = convert(pbx_xhi, {0, 0, 1});

    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_no_stag, 1)); // U
    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_y_stag , 1)); // V
    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_z_stag , 1)); // W
    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_no_stag, 1)); // T

    // *******************************************************************************
    // ylo bdy
    // *******************************************************************************
    plo[1] = -1; plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = -1; phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_ylo(plo, phi);

    Box ylo_plane_no_stag(pbx_ylo);
    Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
    Box ylo_plane_z_stag = convert(pbx_ylo, {0, 0, 1});

    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_x_stag , 1)); // U
    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_no_stag, 1)); // V
    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_z_stag , 1)); // W
    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_no_stag, 1)); // T

    // *******************************************************************************
    // yhi bdy
    // *******************************************************************************
    plo[1] = hi[1] + 1; plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = hi[1] + 1; phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_yhi(plo, phi);

    Box yhi_plane_no_stag(pbx_yhi);
    Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
    Box yhi_plane_z_stag = convert(pbx_yhi, {0, 0, 1});

    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_x_stag , 1)); // U
    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_no_stag, 1)); // V
    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_z_stag , 1)); // W
    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_no_stag, 1)); // T

    // *******************************************************************************

    int ntimes;
    if (ParallelDescriptor::IOProcessor())
    {
        // Read the netcdf file and fill these FABs
        ntimes = BuildFABsFromWRFBdyFile(nc_bdy_file, bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi);

    } // if ParalleDescriptor::IOProcessor()

    amrex::Print() << "NTIMES " << ntimes << std::endl;

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // Make sure all processors know how many times are stored
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int nvars = bdy_data_xlo[0].size();
    for (int nt = 0; nt < ntimes; nt++)
    {
        for (int i = 0; i < nvars; i++)
        {
            ParallelDescriptor::Bcast(bdy_data_xlo[nt][i].dataPtr(),bdy_data_xlo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_xhi[nt][i].dataPtr(),bdy_data_xhi[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_ylo[nt][i].dataPtr(),bdy_data_ylo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_yhi[nt][i].dataPtr(),bdy_data_yhi[nt][i].box().numPts(),ioproc);
        }

        // CONVERT (THETA - 300) to (RHO THETA)
        amrex::Real theta_ref = 300.;
        bdy_data_xlo[nt][3].template plus<RunOn::Device>(theta_ref);
        bdy_data_xhi[nt][3].template plus<RunOn::Device>(theta_ref);
        bdy_data_ylo[nt][3].template plus<RunOn::Device>(theta_ref);
        bdy_data_yhi[nt][3].template plus<RunOn::Device>(theta_ref);

        // Now multiply by rho to get (rho theta) instead of theta
        // bdy_data_xlo[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_xhi[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_ylo[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_yhi[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
    }

    amrex::Print() << "Successfully loaded data from the wrfbdy (output of 'real.exe') NetCDF file" << std::endl << std::endl;
}
#endif // ERF_USE_NETCDF
