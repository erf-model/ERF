#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
//#include <AMReX_InterpFaceRegister.H>

#include <Constants.H>

#include <RK3.H>
#include <EOS.H>

using namespace amrex;

void RK3_stage  (int level,
                 MultiFab& cons_old,  MultiFab& cons_upd,
                 MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old,
                 MultiFab& xmom_upd, MultiFab& ymom_upd, MultiFab& zmom_upd,
                 MultiFab& xvel    , MultiFab& yvel    , MultiFab& zvel    ,
                 MultiFab& source,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt,
                       amrex::InterpFaceRegister* ifr,
                 const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("RK3_stage()",RK3_stage);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // *************************************************************************
    // Fill the ghost cells/faces of the MultiFabs we will need
    // *************************************************************************
    cons_old.FillBoundary(geom.periodicity());

    xmom_old.FillBoundary(geom.periodicity());
    ymom_old.FillBoundary(geom.periodicity());
    zmom_old.FillBoundary(geom.periodicity());

    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // Apply BC on state data at cells
    // Apply BC on velocity data on faces
    // Note that in RK3_advance, the BC was applied on momentum
    amrex::Vector<MultiFab*> vars{&cons_old, &xmom_old, &ymom_old, &zmom_old, &xvel, &yvel, &zvel};
    ERF::applyBCs(geom, vars);

    // *************************************************************************
    // Deal with gravity
    Real gravity = solverChoice.use_gravity? CONST_GRAV: 0.0;
    // CONST_GRAV is a positive constant, but application of grav_gpu to the vertical momentum
    // tendency assumes this quantity is negative.
    //    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, gravity};
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -gravity};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // *************************************************************************
    // Calculate cell-centered eddy viscosity
    //
    // Notes:
    // 1. Need to apply BC on velocity field so that ComputeStrainRate works
    //    properly
    // 2. Need to call FillBoundary and applyBCs to set ghost values on all
    //    boundaries so that InterpolateTurbulentViscosity works properly
    // *************************************************************************
    MultiFab eddyViscosity(cons_old.boxArray(),cons_old.DistributionMap(),1,1);
    if (solverChoice.use_smagorinsky)
    {
        ComputeTurbulentViscosity(xvel, yvel, zvel, cons_old, eddyViscosity, dx, solverChoice);
        eddyViscosity.FillBoundary(geom.periodicity());
        amrex::Vector<MultiFab*> eddyvisc_update{&eddyViscosity};
        ERF::applyBCs(geom, eddyvisc_update);
    }

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
    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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

        const Array4<Real>& nut = eddyViscosity.array(mfi);

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        amrex::ParallelFor(bx, cons_old.nComp(),
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            cell_data_upd(i, j, k, n) = 0.0; // Initialize the updated state eqn term to zero.

            // Add advection terms.
            if (solverChoice.use_state_advection)
                cell_data_upd(i, j, k, n) += (-dt) * AdvectionContributionForState(i, j, k, rho_u, rho_v, rho_w, cell_data_old, n, dx, solverChoice.spatial_order);

            // Add diffusive terms.
            if (solverChoice.use_thermal_diffusion && n == RhoTheta_comp)
                cell_data_upd(i, j, k, n) += dt * DiffusionContributionForState(i, j, k,cell_data_old, RhoTheta_comp, dx, solverChoice);
            if (solverChoice.use_scalar_diffusion && n == RhoScalar_comp)
                cell_data_upd(i, j, k, n) += dt * DiffusionContributionForState(i, j, k,cell_data_old, RhoScalar_comp, dx, solverChoice);

            // Add source terms. TODO: Put this under a if condition when we implement source term
            cell_data_upd(i, j, k, n) += dt * source_fab(i, j, k, n);
        }
        );
        // *********************************************************************
        // Define updates in the RHS of {x, y, z}-momentum equations
        // *********************************************************************
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // x-momentum equation

            rho_u_upd(i, j, k) = 0.0; // Initialize the updated x-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_u_upd(i, j, k) += (-dt) *
                  AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::x, dx, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_u_upd(i, j, k) += dt *
                  DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::x, dx, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_u_upd(i, j, k) += (-dt / dx[0]) *
                  (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i - 1, j, k, RhoTheta_comp)));

            // Add gravity term
            if (solverChoice.use_gravity)
                rho_u_upd(i, j, k) += dt * grav_gpu[0] *
                  InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::x, solverChoice.spatial_order);

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_u_upd(i, j, k) += (-dt) * solverChoice.abl_pressure_grad[0];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_v_loc = 0.25 * (rho_u(i,j+1,k) + rho_u(i,j,k) + rho_u(i-1,j+1,k) + rho_u(i-1,j,k));
                Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                rho_u_upd(i, j, k) += dt * solverChoice.coriolis_factor *
                        (rho_v_loc * solverChoice.sinphi - rho_w_loc * solverChoice.cosphi);
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_u_upd(i, j, k) += dt * solverChoice.abl_geo_forcing[0];

            } // not on coarse-fine boundary
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // y-momentum equation

            rho_v_upd(i, j, k) = 0.0; // Initialize the updated y-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (j == vlo_y && mlo_y(i,j,k)) || (j == vhi_y+1 && mhi_y(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_v_upd(i, j, k) += (-dt) * AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::y, dx, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_v_upd(i, j, k) += dt *
                  DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::y, dx, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_v_upd(i, j, k) += (-dt / dx[1]) *
                  (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i, j - 1, k, RhoTheta_comp)));

            // Add gravity term
            if (solverChoice.use_gravity)
               rho_v_upd(i, j, k) += dt * grav_gpu[1] *
                  InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::y, solverChoice.spatial_order);

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_v_upd(i, j, k) += (-dt) * solverChoice.abl_pressure_grad[1];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                rho_v_upd(i, j, k) += (-dt) * solverChoice.coriolis_factor * rho_u_loc * solverChoice.sinphi;
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_v_upd(i, j, k) += dt * solverChoice.abl_geo_forcing[1];

            } // not on coarse-fine boundary
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) { // z-momentum equation

            rho_w_upd(i, j, k) = 0.0; // Initialize the updated z-mom eqn term to zero

            bool on_coarse_fine_boundary = false;
            if (level > 0)
            {
               on_coarse_fine_boundary =
                 ( (k == vlo_z && mlo_z(i,j,k)) || (k == vhi_z+1 && mhi_z(i,j,k)) );
            }

            if (!on_coarse_fine_boundary)
            {

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_w_upd(i, j, k) += (-dt) *
                    AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::z, dx, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_w_upd(i, j, k) += dt *
                    DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::z, dx, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_w_upd(i, j, k) += (-dt / dx[2]) *
                    (getPgivenRTh(cell_data_old(i, j, k, RhoTheta_comp)) - getPgivenRTh(cell_data_old(i, j, k - 1, RhoTheta_comp)));

            // Add gravity term
            if (solverChoice.use_gravity)
               rho_w_upd(i, j, k) += dt * grav_gpu[2] *
                   InterpolateDensityFromCellToFace(i, j, k, cell_data_old, NextOrPrev::prev, Coord::z, solverChoice.spatial_order);

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_w_upd(i, j, k) += (-dt) * solverChoice.abl_pressure_grad[2];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                rho_w_upd(i, j, k) += dt * solverChoice.coriolis_factor * rho_u_loc * solverChoice.cosphi;
            }
            } // not on coarse-fine boundary
        });
    }
}
