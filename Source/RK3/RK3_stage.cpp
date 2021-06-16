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

    int nvars = cons_old.nComp(); 

    // ************************************************************************************
    // Fill the ghost cells/faces of the MultiFabs we will need
    // Apply BC on state data at cells
    // ************************************************************************************ 
    cons_old.FillBoundary(geom.periodicity());

    // Apply BC on velocity data on faces
    // Note that in RK3_advance, the BC was applied on momentum
    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // **************************************************************************************
    // Deal with gravity
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    Real gravity = solverChoice.use_gravity? CONST_GRAV: 0.0;
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // ************************************************************************************** 
    // 
    // Calculate face-based fluxes to update cell-centered quantities, and 
    //           edge-based and cell-based fluxes to update face-centered quantities
    // 
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
    //TODO: Benchmarking ofperformance. We are computing the fluxes on the fly. 
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

//        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
        amrex::ParallelFor(bx, nvars,
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            cell_data_upd(i, j, k, n) = 0.0; // Initialize the updated state eqn term to zero.

            // Compute advection terms. TODO: Put these under if condition. No need to compute these for pure diffusion problems
            Real xFaceFluxNext, xFaceFluxPrev, yFaceFluxNext, yFaceFluxPrev, zFaceFluxNext, zFaceFluxPrev;
            switch(n) {
            case Rho_comp: // Continuity
                xFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::unity, AdvectingQuantity::rho_u, solverChoice.spatial_order);
                xFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::unity, AdvectingQuantity::rho_u, solverChoice.spatial_order);
                yFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::unity, AdvectingQuantity::rho_v, solverChoice.spatial_order);
                yFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::unity, AdvectingQuantity::rho_v, solverChoice.spatial_order);
                zFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::unity, AdvectingQuantity::rho_w, solverChoice.spatial_order);
                zFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::unity, AdvectingQuantity::rho_w, solverChoice.spatial_order);
                break;
            case RhoTheta_comp: // Temperature
                xFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::theta, AdvectingQuantity::rho_u, solverChoice.spatial_order);
                xFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::theta, AdvectingQuantity::rho_u, solverChoice.spatial_order);
                yFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::theta, AdvectingQuantity::rho_v, solverChoice.spatial_order);
                yFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::theta, AdvectingQuantity::rho_v, solverChoice.spatial_order);
                zFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::theta, AdvectingQuantity::rho_w, solverChoice.spatial_order);
                zFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::theta, AdvectingQuantity::rho_w, solverChoice.spatial_order);
                break;
            case RhoScalar_comp: // Scalar
                xFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::scalar, AdvectingQuantity::rho_u, solverChoice.spatial_order);
                xFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::scalar, AdvectingQuantity::rho_u, solverChoice.spatial_order);
                yFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::scalar, AdvectingQuantity::rho_v, solverChoice.spatial_order);
                yFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::scalar, AdvectingQuantity::rho_v, solverChoice.spatial_order);
                zFaceFluxNext = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::next, AdvectedQuantity::scalar, AdvectingQuantity::rho_w, solverChoice.spatial_order);
                zFaceFluxPrev = ComputeAdvectedQuantityForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, NextOrPrev::prev, AdvectedQuantity::scalar, AdvectingQuantity::rho_w, solverChoice.spatial_order);
                break;
            default:
                amrex::Abort("Error: Conserved quantity index is unrecognized");
            }
            // Add advective terms. TODO: Put this under a if condition
            cell_data_upd(i, j, k, n) += (-dt) *
                    (  (xFaceFluxNext - xFaceFluxPrev) / dx[0]
                     + (yFaceFluxNext - yFaceFluxPrev) / dx[1]
                     + (zFaceFluxNext - zFaceFluxPrev) / dx[2]
                    );

            // Add diffusive terms. TODO: Put this under a if condition

            // Add source terms. TODO: Put this under a if condition
            cell_data_upd(i, j, k, n) += dt * source_fab(i, j, k, n);
        }
        );

        // TODO: Fine-tune dealing with kinematic and eddy viscosity
        Real nu = 0.0; // Obtained from solver choice
        Array4<Real> nut;
        if (solverChoice.use_smagorinsky) // Compute nut, otherwise remains whatever it's initialized to
            ComputeTurbulentViscosity(u, v, w, geom,nut);

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
            Real centFluxXXNext, centFluxXXPrev, edgeFluxXYNext, edgeFluxXYPrev, edgeFluxXZNext, edgeFluxXZPrev;
            centFluxXXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            centFluxXXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxXYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxXYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxXZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::u, AdvectingQuantity::rho_w, solverChoice.spatial_order);
            edgeFluxXZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::u, AdvectingQuantity::rho_w, solverChoice.spatial_order);

            rho_u_upd(i, j, k) += (-dt) * (
                         (centFluxXXNext - centFluxXXPrev) / dx[0] // Contribution to x-mom eqn from advective flux in x-dir
                        +(edgeFluxXYNext - edgeFluxXYPrev) / dx[1] // Contribution to x-mom eqn from advective flux in y-dir
                        +(edgeFluxXZNext - edgeFluxXZPrev) / dx[2] // Contribution to x-mom eqn from advective flux in z-dir
                    );

            // Add diffusive terms
            Real tau11Next, tau11Prev, tau12Next, tau12Prev, tau13Next, tau13Prev;
            tau11Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::x,DiffusionDir::x, geom, nut, solverChoice, nu);
            tau11Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::x,DiffusionDir::x, geom, nut, solverChoice, nu);
            tau12Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::x,DiffusionDir::y, geom, nut, solverChoice, nu);
            tau12Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::x,DiffusionDir::y, geom, nut, solverChoice, nu);
            tau13Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::x,DiffusionDir::z, geom, nut, solverChoice, nu);
            tau13Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::x,DiffusionDir::z, geom, nut, solverChoice, nu);

            rho_u_upd(i, j, k) += dt * (
                    (tau11Next - tau11Prev) / dx[0] // Contribution to x-mom eqn from diffusive flux in x-dir
                   +(tau12Next - tau12Prev) / dx[1] // Contribution to x-mom eqn from diffusive flux in y-dir
                   +(tau13Next - tau13Prev) / dx[2] // Contribution to x-mom eqn from diffusive flux in z-dir
            );

            // Add pressure gradient
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
            Real centFluxYYNext, centFluxYYPrev, edgeFluxYXNext, edgeFluxYXPrev, edgeFluxYZNext, edgeFluxYZPrev;
            edgeFluxYXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxYXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            centFluxYYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            centFluxYYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxYZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::v, AdvectingQuantity::rho_w, solverChoice.spatial_order);
            edgeFluxYZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::v, AdvectingQuantity::rho_w, solverChoice.spatial_order);

            rho_v_upd(i, j, k) += (-dt) * (
                         (edgeFluxYXNext - edgeFluxYXPrev) / dx[0] // Contribution to y-mom eqn from advective flux in x-dir
                        +(centFluxYYNext - centFluxYYPrev) / dx[1] // Contribution to y-mom eqn from advective flux in y-dir
                        +(edgeFluxYZNext - edgeFluxYZPrev) / dx[2] // Contribution to y-mom eqn from advective flux in z-dir
            );

            // Add diffusive terms
            Real tau21Next, tau21Prev, tau22Next, tau22Prev, tau23Next, tau23Prev;
            tau21Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::y, DiffusionDir::x, geom, nut, solverChoice, nu);
            tau21Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::y, DiffusionDir::x, geom, nut, solverChoice, nu);
            tau22Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::y, DiffusionDir::y, geom, nut, solverChoice, nu);
            tau22Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::y, DiffusionDir::y, geom, nut, solverChoice, nu);
            tau23Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::y, DiffusionDir::z, geom, nut, solverChoice, nu);
            tau23Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::y, DiffusionDir::z, geom, nut, solverChoice, nu);

            rho_v_upd(i, j, k) += dt * (
                      (tau21Next - tau21Prev) / dx[0] // Contribution to y-mom eqn from diffusive flux in x-dir
                    + (tau22Next - tau22Prev) / dx[1] // Contribution to y-mom eqn from diffusive flux in y-dir
                    + (tau23Next - tau23Prev) / dx[2] // Contribution to y-mom eqn from diffusive flux in z-dir
            );

            // Add pressure gradient
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
            Real centFluxZZNext, centFluxZZPrev, edgeFluxZXNext, edgeFluxZXPrev, edgeFluxZYNext, edgeFluxZYPrev;
            edgeFluxZXNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxZXPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_u, solverChoice.spatial_order);
            edgeFluxZYNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            edgeFluxZYPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_v, solverChoice.spatial_order);
            centFluxZZNext = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::next, AdvectedQuantity::w, AdvectingQuantity::rho_w, solverChoice.spatial_order);
            centFluxZZPrev = ComputeAdvectedQuantityForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, NextOrPrev::prev, AdvectedQuantity::w, AdvectingQuantity::rho_w, solverChoice.spatial_order);

            rho_w_upd(i, j, k) += (-dt) * (
                         (edgeFluxZXNext - edgeFluxZXPrev) / dx[0] // Contribution to z-mom eqn from advective flux in x-dir
                        +(edgeFluxZYNext - edgeFluxZYPrev) / dx[1] // Contribution to z-mom eqn from advective flux in y-dir
                        +(centFluxZZNext - centFluxZZPrev) / dx[2] // Contribution to z-mom eqn from advective flux in z-dir
            );

            // Add diffusive terms
            Real tau31Next, tau31Prev, tau32Next, tau32Prev, tau33Next, tau33Prev;
            tau31Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::z, DiffusionDir::x, geom, nut, solverChoice, nu);
            tau31Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::z, DiffusionDir::x, geom, nut, solverChoice, nu);
            tau32Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::z, DiffusionDir::y, geom, nut, solverChoice, nu);
            tau32Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::z, DiffusionDir::y, geom, nut, solverChoice, nu);
            tau33Next = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::next, MomentumEqn::z, DiffusionDir::z, geom, nut, solverChoice, nu);
            tau33Prev = ComputeStressTerm(i, j, k, u, v, w, NextOrPrev::prev, MomentumEqn::z, DiffusionDir::z, geom, nut, solverChoice, nu);

            rho_w_upd(i, j, k) += dt * (
                      (tau31Next - tau31Prev) / dx[0] // Contribution to y-mom eqn from diffusive flux in x-dir
                    + (tau32Next - tau32Prev) / dx[1] // Contribution to y-mom eqn from diffusive flux in y-dir
                    + (tau33Next - tau33Prev) / dx[2] // Contribution to y-mom eqn from diffusive flux in z-dir
            );

            // Add pressure gradient
            rho_w_upd(i, j, k) += (-dt / dx[2]) * (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i, j, k - 1, RhoTheta_comp)));

            // Add gravity term
            rho_w_upd(i, j, k) += dt * grav_gpu[2] * InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::z, solverChoice.spatial_order);
        });
    }
}
