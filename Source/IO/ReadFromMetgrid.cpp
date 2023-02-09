#include "ERF.H"
#include "AMReX_FArrayBox.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
ERF::read_from_metgrid(int lev, int idx,
                       Vector<FArrayBox>& NC_xvel_fab, Vector<FArrayBox>& NC_yvel_fab,
                       Vector<FArrayBox>& NC_zvel_fab, Vector<FArrayBox>& NC_rho_fab,
                       Vector<FArrayBox>& NC_rhop_fab, Vector<FArrayBox>& NC_rhotheta_fab,
                       Vector<FArrayBox>& NC_MUB_fab ,
                       Vector<FArrayBox>& NC_MSFU_fab, Vector<FArrayBox>& NC_MSFV_fab,
                       Vector<FArrayBox>& NC_MSFM_fab,
                       Vector<FArrayBox>& NC_SST_fab,
                       Vector<FArrayBox>& NC_C1H_fab , Vector<FArrayBox>& NC_C2H_fab,
                       Vector<FArrayBox>& NC_RDNW_fab,
                       Vector<FArrayBox>& NC_PH_fab, Vector<FArrayBox>& NC_PHB_fab,
                       Vector<FArrayBox>& NC_ALB_fab, Vector<FArrayBox>& NC_PB_fab)
{
    amrex::Print() << "Loading initial data from NetCDF file at level " << lev << std::endl;
    Box input_box;

    if (lev == 0) {
        input_box = geom[0].Domain();
    } else {
        input_box = boxes_at_level[lev][idx];
    }

    // We allocate these here so they exist on all ranks
    Box ubx(input_box); ubx.surroundingNodes(0);
    Box vbx(input_box); vbx.surroundingNodes(1);
    Box wbx(input_box); wbx.surroundingNodes(2);

    Box  zbx(input_box);   zbx.setRange(0,0); zbx.setRange(1,0);
    Box mubx(input_box);   mubx.setRange(2,0);
    Box msfubx(ubx);       msfubx.setRange(2,0);
    Box msfvbx(vbx);       msfvbx.setRange(2,0);
    Box msfmbx(input_box); msfmbx.setRange(2,0);

    Box sstbx(input_box); msfmbx.setRange(2,0); //same size as msfmbx

    // These are all 3D arrays
    NC_xvel_fab[idx].resize(ubx,1);
    NC_yvel_fab[idx].resize(vbx,1);
    NC_zvel_fab[idx].resize(wbx,1);
    NC_rho_fab[idx].resize(input_box,1);
    NC_rhop_fab[idx].resize(input_box, 1);
    NC_rhotheta_fab[idx].resize(input_box,1);

    // Height on w-faces
    NC_PH_fab[idx].resize(wbx,1);
    NC_PHB_fab[idx].resize(wbx,1);

    // Pressure on cell centers
    NC_PB_fab[idx].resize(input_box,1);

    // Reference density on cell centers
    NC_ALB_fab[idx].resize(input_box,1);

    // These are 2D (x-y) arrays
    NC_MUB_fab[idx].resize(mubx,1);
    NC_MSFU_fab[idx].resize(msfubx,1);
    NC_MSFV_fab[idx].resize(msfvbx,1);
    NC_MSFM_fab[idx].resize(msfmbx,1);
    NC_SST_fab[idx].resize(msfmbx,1);

    // These are 1D (z) arrays
    NC_C1H_fab[idx].resize(zbx,1);
    NC_C2H_fab[idx].resize(zbx,1);
    NC_RDNW_fab[idx].resize(zbx,1);

#ifdef AMREX_USE_GPU
        FArrayBox host_NC_xvel_fab    (NC_xvel_fab[idx].box(),     NC_xvel_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_yvel_fab    (NC_yvel_fab[idx].box(),     NC_yvel_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_zvel_fab    (NC_zvel_fab[idx].box(),     NC_zvel_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_rho_fab     (NC_rho_fab[idx].box(),      NC_rho_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_rhop_fab(NC_rhop_fab[idx].box(), NC_rhop_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_rhotheta_fab(NC_rhotheta_fab[idx].box(), NC_rhotheta_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_PH_fab (NC_PH_fab[idx].box(),  NC_PH_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_PHB_fab(NC_PHB_fab[idx].box(), NC_PHB_fab[idx].nComp(),amrex::The_Pinned_Arena());

        FArrayBox host_NC_PB_fab (NC_PB_fab[idx].box() , NC_PB_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_ALB_fab(NC_ALB_fab[idx].box(), NC_ALB_fab[idx].nComp(),amrex::The_Pinned_Arena());

        FArrayBox host_NC_MUB_fab (NC_MUB_fab[idx].box(),  NC_MUB_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_MSFU_fab(NC_MSFU_fab[idx].box(), NC_MSFU_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_MSFV_fab(NC_MSFV_fab[idx].box(), NC_MSFV_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_MSFM_fab(NC_MSFM_fab[idx].box(), NC_MSFM_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_SST_fab(NC_SST_fab[idx].box(), NC_SST_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_C1H_fab (NC_C1H_fab[idx].box(),  NC_C1H_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_C2H_fab (NC_C2H_fab[idx].box(),  NC_C2H_fab[idx].nComp(),amrex::The_Pinned_Arena());
        FArrayBox host_NC_RDNW_fab (NC_RDNW_fab[idx].box(),  NC_RDNW_fab[idx].nComp(),amrex::The_Pinned_Arena());
#else
        FArrayBox host_NC_xvel_fab    (NC_xvel_fab[idx]    , amrex::make_alias, 0, NC_xvel_fab[idx].nComp());
        FArrayBox host_NC_yvel_fab    (NC_yvel_fab[idx]    , amrex::make_alias, 0, NC_yvel_fab[idx].nComp());
        FArrayBox host_NC_zvel_fab    (NC_zvel_fab[idx]    , amrex::make_alias, 0, NC_zvel_fab[idx].nComp());
        FArrayBox host_NC_rho_fab     (NC_rho_fab[idx]     , amrex::make_alias, 0, NC_rho_fab[idx].nComp());
        FArrayBox host_NC_rhop_fab(NC_rhop_fab[idx], amrex::make_alias, 0, NC_rhop_fab[idx].nComp());
        FArrayBox host_NC_rhotheta_fab(NC_rhotheta_fab[idx], amrex::make_alias, 0, NC_rhotheta_fab[idx].nComp());
        FArrayBox host_NC_PH_fab      (NC_PH_fab[idx]      , amrex::make_alias, 0, NC_PH_fab[idx].nComp());
        FArrayBox host_NC_PHB_fab     (NC_PHB_fab[idx]     , amrex::make_alias, 0, NC_PHB_fab[idx].nComp());
        FArrayBox host_NC_PB_fab      (NC_PB_fab[idx]      , amrex::make_alias, 0, NC_PB_fab[idx].nComp());
        FArrayBox host_NC_ALB_fab     (NC_ALB_fab[idx]     , amrex::make_alias, 0, NC_ALB_fab[idx].nComp());
        FArrayBox host_NC_MUB_fab     (NC_MUB_fab[idx]     , amrex::make_alias, 0, NC_MUB_fab[idx].nComp());
        FArrayBox host_NC_MSFU_fab    (NC_MSFU_fab[idx]    , amrex::make_alias, 0, NC_MSFU_fab[idx].nComp());
        FArrayBox host_NC_MSFV_fab    (NC_MSFV_fab[idx]    , amrex::make_alias, 0, NC_MSFV_fab[idx].nComp());
        FArrayBox host_NC_MSFM_fab    (NC_MSFM_fab[idx]    , amrex::make_alias, 0, NC_MSFM_fab[idx].nComp());
        FArrayBox host_NC_SST_fab     (NC_SST_fab[idx]     , amrex::make_alias, 0, NC_SST_fab[idx].nComp());
        FArrayBox host_NC_C1H_fab     (NC_C1H_fab[idx]     , amrex::make_alias, 0, NC_C1H_fab[idx].nComp());
        FArrayBox host_NC_C2H_fab     (NC_C2H_fab[idx]     , amrex::make_alias, 0, NC_C2H_fab[idx].nComp());
        FArrayBox host_NC_RDNW_fab    (NC_RDNW_fab[idx]    , amrex::make_alias, 0, NC_RDNW_fab[idx].nComp());
#endif

        if (ParallelDescriptor::IOProcessor())
        {
            Vector<FArrayBox*> NC_fabs;
            Vector<std::string> NC_names;
            Vector<enum NC_Data_Dims_Type> NC_dim_types;

            NC_fabs.push_back(&host_NC_xvel_fab);      NC_names.push_back("U");    NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_yvel_fab);      NC_names.push_back("V");    NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_zvel_fab);      NC_names.push_back("W");    NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_rho_fab);       NC_names.push_back("ALB");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_rhop_fab),      NC_names.push_back("AL");   NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_rhotheta_fab);  NC_names.push_back("T");    NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_PH_fab);        NC_names.push_back("PH");   NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_PHB_fab);       NC_names.push_back("PHB");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_PB_fab);        NC_names.push_back("PB");   NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_ALB_fab);       NC_names.push_back("ALB");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
            NC_fabs.push_back(&host_NC_MUB_fab);       NC_names.push_back("MUB");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
            NC_fabs.push_back(&host_NC_MSFU_fab);      NC_names.push_back("MAPFAC_UY"); NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
            NC_fabs.push_back(&host_NC_MSFV_fab);      NC_names.push_back("MAPFAC_VY"); NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
            NC_fabs.push_back(&host_NC_MSFM_fab);      NC_names.push_back("MAPFAC_MY"); NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE);
            NC_fabs.push_back(&host_NC_SST_fab);       NC_names.push_back("SST");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); //SST NOTE: Data dims type?
            NC_fabs.push_back(&host_NC_C1H_fab);       NC_names.push_back("C1H");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT);
            NC_fabs.push_back(&host_NC_C2H_fab);       NC_names.push_back("C2H");  NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT);
            NC_fabs.push_back(&host_NC_RDNW_fab);      NC_names.push_back("RDNW"); NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT);

            // Read the netcdf file and fill these FABs
            amrex::Print() << "Building initial FABS from file " << nc_init_file[lev][idx] << std::endl;
            BuildFABsFromWRFInputFile(nc_init_file[lev][idx], NC_names, NC_dim_types, NC_fabs);

        } // if ParalleDescriptor::IOProcessor()

        // We put a barrier here so the rest of the processors wait to do anything until they have the data
        amrex::ParallelDescriptor::Barrier();

        // When an FArrayBox is built, space is allocated on every rank.  However, we only
        //    filled the data in these FABs on the IOProcessor.  So here we broadcast
        //    the data to every rank.

        int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
        ParallelDescriptor::Bcast(host_NC_xvel_fab.dataPtr(),NC_xvel_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_yvel_fab.dataPtr(),NC_yvel_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_zvel_fab.dataPtr(),NC_zvel_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_rho_fab.dataPtr(),NC_rho_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_rhop_fab.dataPtr(), NC_rhop_fab[idx].box().numPts(), ioproc);
        ParallelDescriptor::Bcast(host_NC_rhotheta_fab.dataPtr(),NC_rhotheta_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_PHB_fab.dataPtr() ,NC_PHB_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_PH_fab.dataPtr()  ,NC_PH_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_PB_fab.dataPtr()  ,NC_PB_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_ALB_fab.dataPtr() ,NC_ALB_fab[idx].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(host_NC_MUB_fab.dataPtr() ,NC_MUB_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_MSFU_fab.dataPtr(),NC_MSFU_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_MSFV_fab.dataPtr(),NC_MSFV_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_MSFM_fab.dataPtr(),NC_MSFM_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_SST_fab.dataPtr() ,NC_SST_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_C1H_fab.dataPtr() ,NC_C1H_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_C2H_fab.dataPtr() ,NC_C2H_fab[idx].box().numPts() ,ioproc);
        ParallelDescriptor::Bcast(host_NC_RDNW_fab.dataPtr() ,NC_RDNW_fab[idx].box().numPts() ,ioproc);

#ifdef AMREX_USE_GPU
         Gpu::copy(Gpu::hostToDevice, host_NC_xvel_fab.dataPtr(), host_NC_xvel_fab.dataPtr()+host_NC_xvel_fab.size(),
                                           NC_xvel_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_yvel_fab.dataPtr(), host_NC_yvel_fab.dataPtr()+host_NC_yvel_fab.size(),
                                           NC_yvel_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_zvel_fab.dataPtr(), host_NC_zvel_fab.dataPtr()+host_NC_zvel_fab.size(),
                                           NC_zvel_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rho_fab.dataPtr(), host_NC_rho_fab.dataPtr()+host_NC_rho_fab.size(),
                                           NC_rho_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rhop_fab.dataPtr(), host_NC_rhop_fab.dataPtr()+host_NC_rhop_fab.size(),
                                           NC_rhop_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_rhotheta_fab.dataPtr(), host_NC_rhotheta_fab.dataPtr()+host_NC_rhotheta_fab.size(),
                                           NC_rhotheta_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_PH_fab.dataPtr(), host_NC_PH_fab.dataPtr()+host_NC_PH_fab.size(),
                                           NC_PH_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_PHB_fab.dataPtr(), host_NC_PHB_fab.dataPtr()+host_NC_PHB_fab.size(),
                                           NC_PHB_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_PB_fab.dataPtr(), host_NC_PB_fab.dataPtr()+host_NC_PB_fab.size(),
                                           NC_PB_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_ALB_fab.dataPtr(), host_NC_ALB_fab.dataPtr()+host_NC_ALB_fab.size(),
                                           NC_ALB_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_MUB_fab.dataPtr(), host_NC_MUB_fab.dataPtr()+host_NC_MUB_fab.size(),
                                           NC_MUB_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_MSFU_fab.dataPtr(), host_NC_MSFU_fab.dataPtr()+host_NC_MSFU_fab.size(),
                                           NC_MSFU_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_MSFV_fab.dataPtr(), host_NC_MSFV_fab.dataPtr()+host_NC_MSFV_fab.size(),
                                           NC_MSFV_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_MSFM_fab.dataPtr(), host_NC_MSFM_fab.dataPtr()+host_NC_MSFM_fab.size(),
                                           NC_MSFM_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_SST_fab.dataPtr(), host_NC_SST_fab.dataPtr()+host_NC_SST_fab.size(),
                                           NC_SST_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_C1H_fab.dataPtr(), host_NC_C1H_fab.dataPtr()+host_NC_C1H_fab.size(),
                                           NC_C1H_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_C2H_fab.dataPtr(), host_NC_C2H_fab.dataPtr()+host_NC_C2H_fab.size(),
                                           NC_C2H_fab[idx].dataPtr());
         Gpu::copy(Gpu::hostToDevice, host_NC_RDNW_fab.dataPtr(), host_NC_RDNW_fab.dataPtr()+host_NC_RDNW_fab.size(),
                                           NC_RDNW_fab[idx].dataPtr());
#endif

        //
        // Convert the velocities using the map factors
        //
        const Box& uubx = NC_xvel_fab[idx].box();
        const Array4<Real>    u_arr = NC_xvel_fab[idx].array();
        const Array4<Real> msfu_arr = NC_MSFU_fab[idx].array();
        ParallelFor(uubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            u_arr(i,j,k) /= msfu_arr(i,j,0);
        });

        const Box& vvbx = NC_yvel_fab[idx].box();
        const Array4<Real>    v_arr = NC_yvel_fab[idx].array();
        const Array4<Real> msfv_arr = NC_MSFV_fab[idx].array();
        ParallelFor(vvbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            v_arr(i,j,k) /= msfv_arr(i,j,0);
        });

        const Box& wwbx = NC_zvel_fab[idx].box();
        const Array4<Real>    w_arr = NC_zvel_fab[idx].array();
        const Array4<Real> msfw_arr = NC_MSFM_fab[idx].array();
        ParallelFor(wwbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            w_arr(i,j,k) /= msfw_arr(i,j,0);
        });

        //
        // WRF decomposes (1/rho) rather than rho so rho = 1/(ALB + AL)
        //
        NC_rho_fab[idx].template plus<RunOn::Device>(NC_rhop_fab[idx], 0, 0, 1);
        NC_rho_fab[idx].template invert<RunOn::Device>(1.0);

        const Real theta_ref = 300.0;
        NC_rhotheta_fab[idx].template plus<RunOn::Device>(theta_ref);

        // Now multiply by rho to get (rho theta) instead of theta
        NC_rhotheta_fab[idx].template mult<RunOn::Device>(NC_rho_fab[idx],0,0,1);
}
#endif // ERF_USE_NETCDF
