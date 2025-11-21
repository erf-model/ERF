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

    r_hse.setVal(1.0); p_hse.setVal(1.e5); pi_hse.setVal(1.); th_hse.setVal(300.0); qv_hse.setVal(0.0);

    // ***********************************************************

    const std::string fname = nc_init_file[lev][0];

    FArrayBox NC_rho_fab;
    FArrayBox NC_theta_fab;
    FArrayBox NC_scalar_fab;
    FArrayBox NC_xvel_fab;
    FArrayBox NC_yvel_fab;
    FArrayBox NC_zvel_fab;

    FArrayBox NC_r_hse_fab;
    FArrayBox NC_p_hse_fab;
    FArrayBox NC_th_hse_fab;

    FArrayBox NC_qv_fab;

    Print() << "Loading data from NetCDF file " << fname << " at level " << lev << std::endl;

    Vector<FArrayBox*>  NC_fabs;
    Vector<std::string> NC_fnames;
    Vector<enum NC_Data_Dims_Type> NC_fdim_types;

    NC_fabs.push_back(&NC_rho_fab)   ; NC_fnames.push_back("RHO") ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 0
    NC_fabs.push_back(&NC_theta_fab) ; NC_fnames.push_back("T")   ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 1
    NC_fabs.push_back(&NC_scalar_fab); NC_fnames.push_back("SCAL"); NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 2
    NC_fabs.push_back(&NC_xvel_fab)  ; NC_fnames.push_back("U")   ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 3
    NC_fabs.push_back(&NC_yvel_fab)  ; NC_fnames.push_back("V")   ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 4
    NC_fabs.push_back(&NC_zvel_fab)  ; NC_fnames.push_back("W")   ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 5

    NC_fabs.push_back(&NC_r_hse_fab) ; NC_fnames.push_back("RHO_HSE"); NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);   // 6
    NC_fabs.push_back(&NC_p_hse_fab) ; NC_fnames.push_back("P_HSE")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);   // 7
    NC_fabs.push_back(&NC_th_hse_fab); NC_fnames.push_back("T_HSE")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);   // 8

    NC_fabs.push_back(&NC_qv_fab)    ; NC_fnames.push_back("QV")  ; NC_fdim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE);      // 9

    Vector<int> success; success.resize(NC_fabs.size());

    Print() << "Building initial FABS from file " << fname << std::endl;

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(geom[lev].Domain(), fname, NC_fnames, NC_fdim_types, NC_fabs, success);

    // ***********************************************************

    auto& lev_new = vars_new[lev];

    auto have_moisture = (lev_new[Vars::cons].nComp() > RhoQ1_comp);

    // This defines all the z(i,j,k) values given z(i,j,0)
    make_terrain_fitted_coords(lev, geom[lev], *z_phys_nd[lev], zlevels_stag[lev], phys_bc_type);

    // Default all cell-centered variables to 0
    lev_new[Vars::cons].setVal(0.0);

    // Default density to 1
    Real den_ref = 1.0;
    lev_new[Vars::cons].setVal(den_ref,Rho_comp,1);

    // Default theta to 300; multiply by rho below
    Real theta_ref = 300.0;
    lev_new[Vars::cons].setVal(theta_ref,RhoTheta_comp,1);

    // Default scalar to 0; multiply by rho below
    Real scal_ref = 0.0;
    lev_new[Vars::cons].setVal(scal_ref,RhoScalar_comp,1);

    // Default xvel to 0
    lev_new[Vars::xvel].setVal(0.0,0,0,1);

    // Default yvel to 0
    lev_new[Vars::yvel].setVal(0.0,0,0,1);

    // Default zvel to 0
    lev_new[Vars::zvel].setVal(0.0,0,0,1);

    const int src_comp = 0;
    const int num_comp = 1;

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
            int dest_comp = Rho_comp;
            cons_fab.template copy<RunOn::Device>(NC_rho_fab  , src_comp, dest_comp, num_comp);
        }
        if (success[1]) {
            int dest_comp = RhoTheta_comp;
            cons_fab.template copy<RunOn::Device>(NC_theta_fab, src_comp, dest_comp, num_comp);
        }
        if (success[2]) {
            int dest_comp = RhoScalar_comp;
            cons_fab.template copy<RunOn::Device>(NC_scalar_fab, src_comp, dest_comp, num_comp);
        }

        if (success[3]) {
            int dest_comp = 0;
            xvel_fab.template copy<RunOn::Device>(NC_xvel_fab , src_comp, dest_comp, num_comp);
        }
        if (success[4]) {
            int dest_comp = 0;
            yvel_fab.template copy<RunOn::Device>(NC_yvel_fab , src_comp, dest_comp, num_comp);
        }
        if (success[5]) {
            int dest_comp = 0;
            zvel_fab.template copy<RunOn::Device>(NC_zvel_fab , src_comp, dest_comp, num_comp);
        }

        // HSE vars
        if (success[6]) {
            int dest_comp = 0;
            r_hse_fab.template copy<RunOn::Device>(NC_r_hse_fab , src_comp, dest_comp, num_comp);
        }
        if (success[7]) {
            int dest_comp = 0;
            p_hse_fab.template copy<RunOn::Device>(NC_p_hse_fab , src_comp, dest_comp, num_comp);
        }
        if (success[8]) {
            int dest_comp = 0;
            th_hse_fab.template copy<RunOn::Device>(NC_th_hse_fab , src_comp, dest_comp, num_comp);
        }

        // Moisture vars
        if (success[9]) {
            int dest_comp = RhoQ1_comp;
            if (have_moisture) {
                cons_fab.template copy<RunOn::Device>(NC_qv_fab , src_comp, dest_comp, num_comp);
            } else {
                Print() << "QV was read, but no moisture model is active" << std::endl;
            }
        } else {
            Print() << "QV not found, defaulting to dry" << std::endl;
        }
    } // mf

    lev_new[Vars::xvel].setVal(0.0,0,0,1);
    lev_new[Vars::yvel].setVal(0.0,0,0,1);
    lev_new[Vars::zvel].setVal(0.0,0,0,1);

    MultiFab::Multiply(lev_new[Vars::cons], lev_new[Vars::cons], Rho_comp, RhoTheta_comp, 1, 0);
    MultiFab::Multiply(lev_new[Vars::cons], lev_new[Vars::cons], Rho_comp, RhoScalar_comp, 1, 0);

    // Populate the base state ghost cells and derived data structures if needed
    if (success[6] && success[7] && success[8]) { // have RHO_HSE, P_HSE, and T_HSE
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

            const Array4<const Real>& con_arr = lev_new[Vars::cons].const_array(mfi);

            const Array4<Real>&  r_hse_arr = r_hse.array(mfi);
            const Array4<Real>&  p_hse_arr = p_hse.array(mfi);
            const Array4<Real>& th_hse_arr = th_hse.array(mfi);
            const Array4<Real>& pi_hse_arr = pi_hse.array(mfi);
            const Array4<Real>& qv_hse_arr = (have_moisture) ? qv_hse.array(mfi) : Array4<Real>{};

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

                // qv_hse == qv
                if (have_moisture) qv_hse_arr(i,j,k) = con_arr(ii,jj,kk,RhoQ1_comp);
            });
        }// mfi

    } else {
        // Same code as Initialization/ERF_InitFromWRFInput.cpp
        Print() << "Calculating HSE state!\n";

        int k_dom_lo = geom[lev].Domain().smallEnd(2);
        int k_dom_hi = geom[lev].Domain().bigEnd(2);
        Real tol = 1.0e-10;
        Real grav = CONST_GRAV;
        for ( MFIter mfi(lev_new[Vars::cons],TileNoZ()); mfi.isValid(); ++mfi ) {
            Box bx  = mfi.tilebox();
            int klo = bx.smallEnd(2);
            int khi = bx.bigEnd(2);
            AMREX_ALWAYS_ASSERT((klo == k_dom_lo) && (khi == k_dom_hi));
            bx.makeSlab(2,klo);

            const Array4<const Real>& con_arr    = lev_new[Vars::cons].const_array(mfi);
            const Array4<const Real>& z_arr      = z_phys_nd[lev]->const_array(mfi);
            const Array4<      Real>&  r_hse_arr = r_hse.array(mfi);
            const Array4<      Real>&  p_hse_arr = p_hse.array(mfi);
            const Array4<      Real>& th_hse_arr = th_hse.array(mfi);
            const Array4<      Real>& pi_hse_arr = pi_hse.array(mfi);
            const Array4<      Real>& qv_hse_arr = (have_moisture) ? qv_hse.array(mfi) : Array4<Real>{};

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int /*k*/) noexcept
            {
                // integrate from surface to domain top
                Real dz, F, C;
                Real rho_tot_hi, rho_tot_lo;
                Real z_lo, z_hi;
                Real R_lo, R_hi;
                Real qv_lo, qv_hi;
                Real Th_lo, Th_hi;
                Real P_lo, P_hi;

                // First integrate from sea level to the height at klo
                {
                    // Vertical grid spacing
                    z_lo = 0.0; // corresponding to p_0
                    z_hi = 0.125 * (z_arr(i,j,klo  ) + z_arr(i+1,j,klo  ) + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  )
                                   +z_arr(i,j,klo+1) + z_arr(i+1,j,klo+1) + z_arr(i,j+1,klo+1) + z_arr(i+1,j+1,klo+1));
                    dz = z_hi - z_lo;

                    // Establish known constant
                    qv_lo = (have_moisture) ? con_arr(i,j,klo,RhoQ1_comp) / con_arr(i,j,klo,Rho_comp) : 0.0;
                    Th_lo = con_arr(i,j,klo,RhoTheta_comp) / con_arr(i,j,klo,Rho_comp);
                    P_lo  = p_0;
                    R_lo  = getRhogivenThetaPress(Th_lo, P_lo, R_d/Cp_d, qv_lo);
                    rho_tot_lo = R_lo * (1. + qv_lo);
                    C  = -P_lo + 0.5*rho_tot_lo*grav*dz;

                    // Initial guess and residual
                    qv_hi = (have_moisture) ? con_arr(i,j,klo,RhoQ1_comp) / con_arr(i,j,klo,Rho_comp) : 0.0;
                    Th_hi = con_arr(i,j,klo,RhoTheta_comp) / con_arr(i,j,klo,Rho_comp);
                    P_hi  = p_0;
                    R_hi  = getRhogivenThetaPress(Th_hi, P_hi, R_d/Cp_d, qv_hi);
                    rho_tot_hi = R_hi * (1. + qv_hi);
                    F = P_hi + 0.5*rho_tot_hi*grav*dz + C;

                    // Do iterations
                    HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                                 grav, C, Th_hi,
                                                 qv_hi, qv_hi,
                                                 P_hi, R_hi, F);

                    // Assign data
                     r_hse_arr(i,j,klo) = R_hi;
                     p_hse_arr(i,j,klo) = P_hi;
                    th_hse_arr(i,j,klo) = Th_hi;
                    pi_hse_arr(i,j,klo) = getExnergivenP(p_hse_arr(i,j,klo), R_d/Cp_d);
                    if (have_moisture) qv_hse_arr(i,j,klo) = qv_hi;
                    P_lo = P_hi;
                    z_lo = z_hi;
                }

                for (int k(klo+1); k<=khi; ++k) {
                    // Vertical grid spacing
                  z_hi = 0.125 * (z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                 +z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1));
                  dz   = z_hi - z_lo;

                  // Establish known constant
                  qv_lo = (have_moisture) ? con_arr(i,j,k,RhoQ1_comp) / con_arr(i,j,k,Rho_comp) : 0.0;
                  Th_lo = con_arr(i,j,k,RhoTheta_comp) / con_arr(i,j,k,Rho_comp);
                  R_lo  = getRhogivenThetaPress(Th_lo, P_lo, R_d/Cp_d, qv_lo);
                  rho_tot_lo = R_lo * (1. + qv_lo);
                  C  = -P_lo + 0.5*rho_tot_lo*grav*dz;

                  // Initial guess and residual
                  qv_hi = (have_moisture) ? con_arr(i,j,k,RhoQ1_comp) / con_arr(i,j,k,Rho_comp) : 0.0;
                  Th_hi = con_arr(i,j,k,RhoTheta_comp) / con_arr(i,j,k,Rho_comp);
                  R_hi  = getRhogivenThetaPress(Th_hi, P_hi, R_d/Cp_d, qv_hi);
                  rho_tot_hi = R_hi * (1. + qv_hi);
                  F = P_hi + 0.5*rho_tot_hi*grav*dz + C;

                  // Do iterations
                  HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                               grav, C, Th_hi,
                                               qv_hi, qv_hi,
                                               P_hi, R_hi, F);

                  // Assign data
                   r_hse_arr(i,j,k) = R_hi;
                   p_hse_arr(i,j,k) = P_hi;
                  th_hse_arr(i,j,k) = Th_hi;
                  pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), R_d/Cp_d);
                  if (have_moisture) qv_hse_arr(i,j,k) = qv_hi;
                  P_lo = P_hi;
                  z_lo = z_hi;
                }
            });
        } // mfi
    }

    // Fill boundary if periodic
    r_hse.FillBoundary(geom[lev].periodicity());
    p_hse.FillBoundary(geom[lev].periodicity());
    th_hse.FillBoundary(geom[lev].periodicity());
    pi_hse.FillBoundary(geom[lev].periodicity());
    qv_hse.FillBoundary(geom[lev].periodicity());
}
#endif // ERF_USE_NETCDF
