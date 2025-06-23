/**
 * \file ERF_InitFromNCFile.cpp
 */
#include <ERF.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_ProbCommon.H>
#include <ERF_DataStruct.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

#include <ERF_NCWpsFile.H>

/**
 * Initializes ERF data using data supplied by an external NetCDF file.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_ncfile (int lev)
{
    if (nc_init_file.empty()) {
        Error("NetCDF initialization file name must be provided via input");
    }

    if (nc_init_file[lev].empty()) {
        Error("NetCDF initialization file name must be provided via input");
    }

    // ***********************************************************
    // Initialize base state to be non-zero so we don't divide by zero anywhere
    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    r_hse.setVal(1.0); p_hse.setVal(1.0); pi_hse.setVal(1.); th_hse.setVal(1.0); qv_hse.setVal(0.0);

    // ***********************************************************

    const std::string fname = nc_init_file[lev][0];

    FArrayBox NC_rho_fab;
    FArrayBox NC_theta_fab;
    FArrayBox NC_xvel_fab;
    FArrayBox NC_yvel_fab;
    FArrayBox NC_zvel_fab;

    FArrayBox NC_r_hse_fab;
    FArrayBox NC_p_hse_fab;
    FArrayBox NC_th_hse_fab;

    Print() << "Loading data from NetCDF file " << fname << " at level " << lev << std::endl;

    Vector<FArrayBox*>  NC_fabs;
    Vector<std::string> NC_fnames;
    Vector<enum NC_Data_Dims_Type> NC_fdim_types;

    NC_fabs.push_back(&NC_rho_fab)  ; NC_fnames.push_back("RHO"); NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_theta_fab); NC_fnames.push_back("T")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_xvel_fab) ; NC_fnames.push_back("U")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_yvel_fab) ; NC_fnames.push_back("V")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_zvel_fab) ; NC_fnames.push_back("W")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);

    NC_fabs.push_back(&NC_r_hse_fab)  ; NC_fnames.push_back("RHO_HSE"); NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_p_hse_fab)    ; NC_fnames.push_back("P_HSE")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);
    NC_fabs.push_back(&NC_th_hse_fab)    ; NC_fnames.push_back("T_HSE")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);

    Vector<int> success; success.resize(NC_fabs.size());

    Print() << "Building initial FABS from file " << fname << std::endl;

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(geom[lev].Domain(), fname, NC_fnames, NC_fdim_types, NC_fabs, success);

    // ***********************************************************

    auto& lev_new = vars_new[lev];

    // This defines all the z(i,j,k) values given z(i,j,0)
    make_terrain_fitted_coords(lev, geom[lev], *z_phys_nd[lev], zlevels_stag[lev], phys_bc_type);

    // Default all cell-centered variables to 0
    lev_new[Vars::cons].setVal(0.0);

    // Default density to 1
    Real den_ref = 1.0;
    lev_new[Vars::cons].setVal(den_ref,0,1);

    // Default theta to 300; multiply by rho below
    Real theta_ref = 300.0;
    lev_new[Vars::cons].setVal(theta_ref,1,1);

    // Default xvel to 0
    lev_new[Vars::xvel].setVal(0.0,0,0,1);

    // Default yvel to 0
    lev_new[Vars::yvel].setVal(0.0,0,0,1);

    // Default zvel to 0
    lev_new[Vars::zvel].setVal(0.0,0,0,1);

    int src_comp = 0;
    int num_comp = 1;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
        FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
        FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];

        FArrayBox &r_hse_fab  = r_hse[mfi];
        FArrayBox &p_hse_fab  = p_hse[mfi];
        FArrayBox &th_hse_fab = th_hse[mfi];

        // Copy on intersect...
        if (success[0]) {
            int dest_comp = 0;
            cons_fab.template copy<RunOn::Device>(NC_rho_fab  , src_comp, dest_comp, num_comp);
        }
        if (success[1]) {
            int dest_comp = 1;
            cons_fab.template copy<RunOn::Device>(NC_theta_fab, src_comp, dest_comp, num_comp);
            // cons_fab.template plus<RunOn::Device>(theta_ref, 1);
        }

        if (success[2]) {
            int dest_comp = 0;
            xvel_fab.template copy<RunOn::Device>(NC_xvel_fab , src_comp, dest_comp, num_comp);
        }
        if (success[3]) {
            int dest_comp = 0;
            yvel_fab.template copy<RunOn::Device>(NC_yvel_fab , src_comp, dest_comp, num_comp);
        }
        if (success[4]) {
            int dest_comp = 0;
            zvel_fab.template copy<RunOn::Device>(NC_zvel_fab , src_comp, dest_comp, num_comp);
        }

        if (success[5]) {
            int dest_comp = 0;
            r_hse_fab.template copy<RunOn::Device>(NC_r_hse_fab , src_comp, dest_comp, num_comp);
        }
        if (success[6]) {
            int dest_comp = 0;
            p_hse_fab.template copy<RunOn::Device>(NC_p_hse_fab , src_comp, dest_comp, num_comp);
        }
        if (success[7]) {
            int dest_comp = 0;
            th_hse_fab.template copy<RunOn::Device>(NC_th_hse_fab , src_comp, dest_comp, num_comp);
        }

    } // mf

    MultiFab::Multiply(lev_new[Vars::cons], lev_new[Vars::cons], 0, 1, 1, 0);

    // Populate the base state ghost cells and derived data structures if needed
    if (success[6]) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(p_hse, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Box vbx  = mfi.validbox();
            Box gtbx = mfi.growntilebox();

            int ilo = vbx.smallEnd(0); int ihi = vbx.bigEnd(0);
            int jlo = vbx.smallEnd(1); int jhi = vbx.bigEnd(1);
            int klo = vbx.smallEnd(2); int khi = vbx.bigEnd(2);

            Array4<Real>  r_hse_arr = r_hse.array(mfi);
            Array4<Real>  p_hse_arr = p_hse.array(mfi);
            Array4<Real> th_hse_arr = th_hse.array(mfi);
            Array4<Real> pi_hse_arr = pi_hse.array(mfi);

            ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Base state needs ghost cells filled, protect FAB access
                int ii = std::max(i , ilo);
                    ii = std::min(ii, ihi);
                int jj = std::max(j , jlo);
                    jj = std::min(jj, jhi);
                int kk = std::max(k , klo);
                    kk = std::min(kk, khi);

                 r_hse_arr(i,j,k) =  r_hse_arr(ii,jj,kk);
                 p_hse_arr(i,j,k) =  p_hse_arr(ii,jj,kk);
                th_hse_arr(i,j,k) = th_hse_arr(ii,jj,kk);
                pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(ii,jj,kk), R_d/Cp_d);
            });
        }// mfi

        // Fill boundary if periodic
        r_hse.FillBoundary(geom[lev].periodicity());
        p_hse.FillBoundary(geom[lev].periodicity());
        th_hse.FillBoundary(geom[lev].periodicity());
        pi_hse.FillBoundary(geom[lev].periodicity());

    }// success[6]
}
#endif // ERF_USE_NETCDF
