#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
//#include <AMReX_VisMF.H>

#include <Constants.H>

#include <ERF.H>
#include <RK3.H>
#include <EOS.H>

using namespace amrex;

void RK3_stage  (MultiFab& cons_old,  MultiFab& cons_upd,
                 MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old, 
                 MultiFab& xmom_upd, MultiFab& ymom_upd, MultiFab& zmom_upd, 
                 MultiFab& xvel    , MultiFab& yvel    , MultiFab& zvel    ,  
                 MultiFab& prim    , MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt,
                 const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("RK3_stage()",RK3_stage);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray(); // TODO: Remove this once all the flux computations are modularized

    // ************************************************************************************
    // Fill the ghost cells/faces of the MultiFabs we will need
    // Apply BC on state data at cells
    // ************************************************************************************ 
    cons_old.FillBoundary(geom.periodicity());

    xmom_old.FillBoundary(geom.periodicity());
    ymom_old.FillBoundary(geom.periodicity());
    zmom_old.FillBoundary(geom.periodicity());

    // Apply BC on velocity data on faces
    // Note that in RK3_advance, the BC was applied on momentum
    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    amrex::Vector<MultiFab*> vars{&cons_old, &xvel, &yvel, &zvel};

    ERF::applyBCs(geom, vars);

    // **************************************************************************************
    // Deal with gravity
    Real gravity = solverChoice.use_gravity? CONST_GRAV: 0.0;
    // CONST_GRAV is a positive constant, but application of grav_gpu to the vertical momentum
    // tendency assumes this quantity is negative.
    //    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // **************************************************************************************
    // Calculate face-based fluxes to update cell-centered quantities, and 
    //           edge-based and cell-based fluxes to update face-centered quantities
    // **************************************************************************************
    // TODO: No need for 'CalcAdvFlux'. Remove/clean this
//    CalcAdvFlux(cons_old, xmom_old, ymom_old, zmom_old,
//                xvel    , yvel    , zvel    ,
//                faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux,
//                geom, dxp, dt,
//                solverChoice);
//#if 0
//    CalcDiffFlux(cons_old, xmom_old, ymom_old, zmom_old,
//                 xvel    , yvel    , zvel    ,
//                 eta, zeta, kappa,
//                 faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux,
//                 geom, dxp, dt,
//                 solverChoice);
//#endif

    // **************************************************************************************
    // Define updates in the current RK stage, fluxes are computed here itself
    //TODO: Benchmarking of performance. We are computing the fluxes on the fly.
    //If the performance slows, consider saving all the fluxes apriori and accessing them here.

    // ************************************************************************************** 
    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cell_data_old     = cons_old.array(mfi);
        const Array4<Real> & cell_data_upd     = cons_upd.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real> & u = xvel.array(mfi);
        const Array4<Real> & v = yvel.array(mfi);
        const Array4<Real> & w = zvel.array(mfi);

        const Array4<Real>& rho_u = xmom_old.array(mfi);
        const Array4<Real>& rho_v = ymom_old.array(mfi);
        const Array4<Real>& rho_w = zmom_old.array(mfi);

        const Array4<Real>& rho_u_upd = xmom_upd.array(mfi);
        const Array4<Real>& rho_v_upd = ymom_upd.array(mfi);
        const Array4<Real>& rho_w_upd = zmom_upd.array(mfi);

        // TODO: We won't need faceflux, edgeflux, and centflux when using the new code architecture.Remove them. Fluxes are computed here itself.
//        Array4<Real const> const& xflux = faceflux[0].array(mfi);
//        Array4<Real const> const& yflux = faceflux[1].array(mfi);
//        Array4<Real const> const& zflux = faceflux[2].array(mfi);

//        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
//        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);
//
//        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
//        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);
//
//        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
//        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);
//
//        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
//        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
//        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

//        amrex::ParallelFor(bx, cons_old.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
//        {
//            cell_data_upd(i,j,k,n) = - dt *
//                ( AMREX_D_TERM(  (xflux(i+1,j,k,n) - xflux(i,j,k,n)) / dx[0],
//                               + (yflux(i,j+1,k,n) - yflux(i,j,k,n)) / dx[1],
//                               + (zflux(i,j,k+1,n) - zflux(i,j,k,n)) / dx[2])
//                                                                                       )
//                + dt*source_fab(i,j,k,n);
//        });

        // **************************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************************
        amrex::ParallelFor(bx, cons_old.nComp(),
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            cell_data_upd(i, j, k, n) = 0.0; // Initialize the updated state eqn term to zero.

            // Add advection terms.
            if (solverChoice.use_state_advection)
                cell_data_upd(i, j, k, n) += (-dt) * AdvectionContributionForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, n, geom, solverChoice.spatial_order);

            // Add diffusive terms.
            if (solverChoice.use_thermal_diffusion && n == RhoTheta_comp)
                cell_data_upd(i, j, k, n) += dt * DiffusionContributionForState(i, j, k,cell_data_old, RhoTheta_comp, geom, solverChoice);
            if (solverChoice.use_scalar_diffusion && n == RhoScalar_comp)
                cell_data_upd(i, j, k, n) += dt * DiffusionContributionForState(i, j, k,cell_data_old, RhoScalar_comp, geom, solverChoice);

            // Add source terms. TODO: Put this under a if condition when we implement source term
            cell_data_upd(i, j, k, n) += dt * source_fab(i, j, k, n);
        }
        );

        // TODO: Fine-tune dealing with kinematic and eddy viscosity
        Array4<Real> nut;
//        if (solverChoice.use_smagorinsky) // Compute nut, otherwise remains whatever it's initialized to
//            ComputeTurbulentViscosity(u, v, w, geom,nut);

        // **************************************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // **************************************************************************************
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // x-momentum equation
//            rho_u_upd(i,j,k) =
//                    -dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
//                    -dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
//                    -dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
//                    +0.5*dt*grav_gpu[0]*(cell_data_old(i-1,j,k,0)+cell_data_old(i,j,k,0));

            rho_u_upd(i, j, k) = 0.0; // Initialize the updated x-mom eqn term to zero

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_u_upd(i, j, k) += (-dt) * AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::x, geom, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_u_upd(i, j, k) += dt * DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::x, geom, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_u_upd(i, j, k) += (-dt / dx[0]) * (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i - 1, j, k, RhoTheta_comp)));

            // Add gravity term
            rho_u_upd(i, j, k) += dt * grav_gpu[0] * InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::x, solverChoice.spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // y-momentum equation
//            rho_v_upd(i,j,k) =
//                    -dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
//                    -dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
//                    -dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
//                    +0.5*dt*grav_gpu[1]*(cell_data_old(i,j-1,k,0)+cell_data_old(i,j,k,0));

            rho_v_upd(i, j, k) = 0.0; // Initialize the updated y-mom eqn term to zero

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_v_upd(i, j, k) += (-dt) * AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::y, geom, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_v_upd(i, j, k) += dt * DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::y, geom, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_v_upd(i, j, k) += (-dt / dx[1]) * (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i, j - 1, k, RhoTheta_comp)));

            // Add gravity term
            rho_v_upd(i, j, k) += dt * grav_gpu[1] * InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::y, solverChoice.spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation
//            rho_w_upd(i,j,k) =
//                    -dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
//                    -dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
//                    -dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
//                    +0.5*dt*grav_gpu[2]*(cell_data_old(i,j,k-1,0)+cell_data_old(i,j,k,0));

            rho_w_upd(i, j, k) = 0.0; // Initialize the updated z-mom eqn term to zero

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_w_upd(i, j, k) += (-dt) * AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::z, geom, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_w_upd(i, j, k) += dt * DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::z, geom, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_w_upd(i, j, k) += (-dt / dx[2]) * (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i, j, k - 1, RhoTheta_comp)));

            // Add gravity term
            rho_w_upd(i, j, k) += dt * grav_gpu[2] * InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::z, solverChoice.spatial_order);
        });
    }
}
