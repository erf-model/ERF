#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <ERF_Constants.H>
#include <SpatialStencils.H>
#include <TimeIntegration.H>
#include <EOS.H>

using namespace amrex;

void erf_fast_rhs (int level,
                   Vector<MultiFab >& S_rhs,
                   const Vector<MultiFab >& S_stage_data,
                   const Vector<MultiFab >& S_data,
                   std::array< MultiFab, AMREX_SPACEDIM>&  advflux,
                   std::array< MultiFab, AMREX_SPACEDIM>& diffflux,
                   const amrex::Geometry geom, const amrex::Real dt,
                         amrex::InterpFaceRegister* ifr,
                   const SolverChoice& solverChoice,
                   const amrex::Real* dptr_dens_hse, const amrex::Real* dptr_pres_hse,
                   const amrex::Real* dptr_rayleigh_tau, const amrex::Real* dptr_rayleigh_ubar,
                   const amrex::Real* dptr_rayleigh_vbar, const amrex::Real* dptr_rayleigh_thetabar)
{
    BL_PROFILE_VAR("erf_rhs()",erf_rhs);

    int klo = geom.Domain().smallEnd()[2];
    int khi = geom.Domain().bigEnd()[2];

    const GpuArray<Real, AMREX_SPACEDIM> dxInv = geom.InvCellSizeArray();
    const auto& ba = S_data[IntVar::cons].boxArray();
    const auto& dm = S_data[IntVar::cons].DistributionMap();

    amrex::MultiFab Delta_rho_u(    convert(ba,IntVect(1,0,0)), dm, 1, 1);
    amrex::MultiFab Delta_rho_v(    convert(ba,IntVect(0,1,0)), dm, 1, 1);
    amrex::MultiFab Delta_rho_w(    convert(ba,IntVect(0,0,1)), dm, 1, 1);
    amrex::MultiFab Delta_rho  (            ba                , dm, 1, 1);
    amrex::MultiFab Delta_rho_theta(        ba                , dm, 1, 1);

    // Create delta_rho_u/v/w/theta  = U'', V'', W'', Theta'' in the docs
    MultiFab::Copy(Delta_rho_u    , S_data[IntVar::xmom], 0, 0, 1, 1);
    MultiFab::Copy(Delta_rho_v    , S_data[IntVar::ymom], 0, 0, 1, 1);
    MultiFab::Copy(Delta_rho_w    , S_data[IntVar::zmom], 0, 0, 1, 1);
    MultiFab::Copy(Delta_rho      , S_data[IntVar::cons], Rho_comp     , 0, 1,1);
    MultiFab::Copy(Delta_rho_theta, S_data[IntVar::cons], RhoTheta_comp, 0, 1,1);

    MultiFab::Subtract(Delta_rho_u    , S_stage_data[IntVar::xmom], 0, 0, 1, 1);
    MultiFab::Subtract(Delta_rho_v    , S_stage_data[IntVar::ymom], 0, 0, 1, 1);
    MultiFab::Subtract(Delta_rho_w    , S_stage_data[IntVar::zmom], 0, 0, 1, 1);
    MultiFab::Subtract(Delta_rho      , S_stage_data[IntVar::cons], Rho_comp     , 0, 1, 1);
    MultiFab::Subtract(Delta_rho_theta, S_stage_data[IntVar::cons], RhoTheta_comp, 0, 1, 1);

    // Not sure if we need these
    Delta_rho_u.FillBoundary(geom.periodicity());
    Delta_rho_v.FillBoundary(geom.periodicity());
    Delta_rho_w.FillBoundary(geom.periodicity());
    Delta_rho_theta.FillBoundary(geom.periodicity());

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    const iMultiFab *mlo_mf_x, *mhi_mf_x;
    const iMultiFab *mlo_mf_y, *mhi_mf_y;
    const iMultiFab *mlo_mf_z, *mhi_mf_z;

    if (level > 0)
    {
        mlo_mf_x = &(ifr->mask(Orientation(0,Orientation::low)));
        mhi_mf_x = &(ifr->mask(Orientation(0,Orientation::high)));
        mlo_mf_y = &(ifr->mask(Orientation(1,Orientation::low)));
        mhi_mf_y = &(ifr->mask(Orientation(1,Orientation::high)));
        mlo_mf_z = &(ifr->mask(Orientation(2,Orientation::low)));
        mhi_mf_z = &(ifr->mask(Orientation(2,Orientation::high)));
    }

    // *************************************************************************
    // Define updates in the current RK stage, fluxes are computed here itself
    //TODO: Benchmarking of performance. We are computing the fluxes on the fly.
    //If the performance slows, consider saving all the fluxes apriori and accessing them here.

    // *************************************************************************
    for ( MFIter mfi(S_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        Box const& valid_bx = mfi.validbox();
        int vlo_x = valid_bx.smallEnd(0);
        int vhi_x = valid_bx.bigEnd(0);
        int vlo_y = valid_bx.smallEnd(1);
        int vhi_y = valid_bx.bigEnd(1);
        int vlo_z = valid_bx.smallEnd(2);
        int vhi_z = valid_bx.bigEnd(2);

        auto mlo_x = (level > 0) ? mlo_mf_x->const_array(mfi) : Array4<const int>{};
        auto mhi_x = (level > 0) ? mhi_mf_x->const_array(mfi) : Array4<const int>{};
        auto mlo_y = (level > 0) ? mlo_mf_y->const_array(mfi) : Array4<const int>{};
        auto mhi_y = (level > 0) ? mhi_mf_y->const_array(mfi) : Array4<const int>{};
        auto mlo_z = (level > 0) ? mlo_mf_z->const_array(mfi) : Array4<const int>{};
        auto mhi_z = (level > 0) ? mhi_mf_z->const_array(mfi) : Array4<const int>{};

        const Array4<const Real> & cell_data        = S_data[IntVar::cons].array(mfi);
        const Array4<Real> & cell_rhs         = S_rhs[IntVar::cons].array(mfi);
        const Array4<const Real> & cell_stage_data  = S_stage_data[IntVar::cons].array(mfi);

        const Array4<Real>& delta_rho_u     = Delta_rho_u.array(mfi);
        const Array4<Real>& delta_rho_v     = Delta_rho_v.array(mfi);
        const Array4<Real>& delta_rho_w     = Delta_rho_w.array(mfi);
        const Array4<Real>& delta_rho       = Delta_rho.array(mfi);
        const Array4<Real>& delta_rho_theta = Delta_rho_theta.array(mfi);

        const Array4<Real>& rho_u_rhs = S_rhs[IntVar::xmom].array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVar::ymom].array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVar::zmom].array(mfi);

        const Array4<Real>& xflux_rhs = S_rhs[IntVar::xflux].array(mfi);
        const Array4<Real>& yflux_rhs = S_rhs[IntVar::yflux].array(mfi);
        const Array4<Real>& zflux_rhs = S_rhs[IntVar::zflux].array(mfi);

        // These are temporaries we use to add to the S_rhs for the fluxes
        const Array4<Real>& advflux_x = advflux[0].array(mfi);
        const Array4<Real>& advflux_y = advflux[1].array(mfi);
        const Array4<Real>& advflux_z = advflux[2].array(mfi);

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        amrex::ParallelFor(bx, S_data[IntVar::cons].nComp(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {

            if (n == Rho_comp)
            {
                // We use the most current momenta here to update rho
                cell_rhs(i, j, k, n) += -AdvectionContributionForState(i, j, k, delta_rho_u, delta_rho_v, delta_rho_w,
                                         cell_data, Rho_comp,
                                         advflux_x, advflux_y, advflux_z, dxInv, solverChoice.spatial_order);
            } else if (n == RhoTheta_comp) {
                // We use the most current momenta but the "lagged" theta here to update (rho theta)
                cell_rhs(i, j, k, n) += -AdvectionContributionForState(i, j, k, delta_rho_u, delta_rho_v, delta_rho_w,
                                         cell_stage_data, RhoTheta_comp,
                                         advflux_x, advflux_y, advflux_z, dxInv, solverChoice.spatial_order);
            } else {
                cell_rhs(i, j, k, n) = 0.0;
            }

        }
        );

        // Compute the RHS for the flux terms from this stage --
        //     we do it this way so we don't double count
        amrex::ParallelFor(tbx, S_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == Rho_comp && n == RhoTheta_comp)
                xflux_rhs(i,j,k,n) = advflux_x(i,j,k,n);
            else
                xflux_rhs(i,j,k,n) = 0.0;
        });
        amrex::ParallelFor(tby, S_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == Rho_comp && n == RhoTheta_comp)
                yflux_rhs(i,j,k,n) = advflux_y(i,j,k,n);
            else
                yflux_rhs(i,j,k,n) = 0.0;
        });
        amrex::ParallelFor(tbz, S_data[IntVar::cons].nComp(),
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (n == Rho_comp && n == RhoTheta_comp)
                zflux_rhs(i,j,k,n) = advflux_z(i,j,k,n);
            else
                zflux_rhs(i,j,k,n) = 0.0;
        });

        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // x-momentum equation

            rho_u_rhs(i, j, k) = 0.0; // Initialize the updated x-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real exner_pi_l = std::pow(getPgivenRTh(cell_stage_data(i-1,j,k,RhoTheta_comp))/p_0,R_d/c_p);
                Real exner_pi_r = std::pow(getPgivenRTh(cell_stage_data(i  ,j,k,RhoTheta_comp))/p_0,R_d/c_p);
                Real exner_pi_c =  0.5 * (exner_pi_l + exner_pi_r);
                rho_u_rhs(i, j, k) += (-dxInv[0]) * Gamma * R_d * exner_pi_c *
                  (delta_rho_theta(i  ,j,k,0) - delta_rho_theta(i-1,j,k,0));

            } // not on coarse-fine boundary
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // y-momentum equation

            rho_v_rhs(i, j, k) = 0.0; // Initialize the updated y-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (j == vlo_y && mlo_y(i,j,k)) || (j == vhi_y+1 && mhi_y(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                Real exner_pi_l = std::pow(getPgivenRTh(cell_stage_data(i,j-1,k,RhoTheta_comp))/p_0,R_d/c_p);
                Real exner_pi_r = std::pow(getPgivenRTh(cell_stage_data(i,j  ,k,RhoTheta_comp))/p_0,R_d/c_p);
                Real exner_pi_c =  0.5 * (exner_pi_l + exner_pi_r);
                rho_v_rhs(i, j, k) += (-dxInv[1]) * Gamma * R_d * exner_pi_c *
                  (delta_rho_theta(i,j  ,k,0) - delta_rho_theta(i,j-1,k,0));

            } // not on coarse-fine boundary
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation

            rho_w_rhs(i, j, k) = 0.0; // Initialize the updated z-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (k == vlo_z && mlo_z(i,j,k)) || (k == vhi_z+1 && mhi_z(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {
                Real exner_pi_l = std::pow(getPgivenRTh(cell_stage_data(i,j,k-1,RhoTheta_comp))/p_0,R_d/c_p);
                Real exner_pi_r = std::pow(getPgivenRTh(cell_stage_data(i,j,k  ,RhoTheta_comp))/p_0,R_d/c_p);
                Real exner_pi_c =  0.5 * (exner_pi_l + exner_pi_r);

                // Add (negative) gradient of (rho theta) multiplied by lagged "pi"
                rho_w_rhs(i, j, k) += -1.0 * Gamma * R_d * exner_pi_c *
                  (delta_rho_theta(i,j,k,0) - delta_rho_theta(i,j,k-1,0))*dxInv[2];

                // Add gravity term
                Real rhobar = 0.5 * (dptr_dens_hse[k] + dptr_dens_hse[k-1]);
                Real  pibar = std::pow(0.5 * (dptr_pres_hse[k] + dptr_pres_hse[k-1]),R_d/c_p);
                Real wadv = 0.0; // We will not use this since we are using 2nd order interpolation
                rho_w_rhs(i, j, k) += grav_gpu[2] * (
                    InterpolateFromCellOrFace(i, j, k, delta_rho, 0, wadv, Coord::z, 2)
                    - (R_d / (c_p - R_d)) * exner_pi_c * rhobar / pibar *
                    InterpolateFromCellOrFace(i, j, k, delta_rho_theta,             0, wadv, Coord::z, 2) /
                    InterpolateFromCellOrFace(i, j, k, cell_stage_data, RhoTheta_comp, wadv, Coord::z, 2) );

            } // not on coarse-fine boundary
        });
    }
}
