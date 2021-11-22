#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <ERF_Constants.H>
#include <TimeIntegration.H>
#include <EOS.H>

using namespace amrex;

void erf_rhs (int level,
              Vector<std::unique_ptr<MultiFab> >& S_rhs,
              const Vector<std::unique_ptr<MultiFab> >& S_data,
              MultiFab& source,
              std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
              const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt,
                    amrex::InterpFaceRegister* ifr,
              const SolverChoice& solverChoice,
              const amrex::Real* dptr_dens_hse, const amrex::Real* dptr_pres_hse)
{
    BL_PROFILE_VAR("erf_rhs()",erf_rhs);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    const auto& ba = S_data[IntVar::cons]->boxArray();
    const auto& dm = S_data[IntVar::cons]->DistributionMap();

    amrex::MultiFab xvel(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    amrex::MultiFab yvel(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    amrex::MultiFab zvel(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    MomentumToVelocity(xvel, yvel, zvel,
                       *S_data[IntVar::cons],
                       *S_data[IntVar::xmom],
                       *S_data[IntVar::ymom],
                       *S_data[IntVar::zmom],
                       1, solverChoice);

    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // Apply BC on velocity data on faces
    // Note that the BC was already applied on momentum
    amrex::Vector<MultiFab*> vel_vars{&xvel, &yvel, &zvel};
    ERF::applyBCs(geom, vel_vars);

    // *************************************************************************
    // Set gravity as a vector
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -solverChoice.gravity};
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
    MultiFab eddyViscosity(S_data[IntVar::cons]->boxArray(),S_data[IntVar::cons]->DistributionMap(),1,1);
    if (solverChoice.les_type == LESType::Smagorinsky)
    {
        ComputeTurbulentViscosity(xvel, yvel, zvel, *S_data[IntVar::cons], eddyViscosity, dx, solverChoice);
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
    for ( MFIter mfi(*S_data[IntVar::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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

        const Array4<Real> & cell_data  = S_data[IntVar::cons]->array(mfi);
        const Array4<Real> & cell_rhs   = S_rhs[IntVar::cons]->array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real> & u = xvel.array(mfi);
        const Array4<Real> & v = yvel.array(mfi);
        const Array4<Real> & w = zvel.array(mfi);

        const Array4<Real>& rho_u = S_data[IntVar::xmom]->array(mfi);
        const Array4<Real>& rho_v = S_data[IntVar::ymom]->array(mfi);
        const Array4<Real>& rho_w = S_data[IntVar::zmom]->array(mfi);

        const Array4<Real>& rho_u_rhs = S_rhs[IntVar::xmom]->array(mfi);
        const Array4<Real>& rho_v_rhs = S_rhs[IntVar::ymom]->array(mfi);
        const Array4<Real>& rho_w_rhs = S_rhs[IntVar::zmom]->array(mfi);

        const Array4<Real>& nut = eddyViscosity.array(mfi);

        // **************************************************************************
        // Define updates in the RHS of continuity, temperature, and scalar equations
        // **************************************************************************
        amrex::ParallelFor(bx, S_data[IntVar::cons]->nComp(),
       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
            cell_rhs(i, j, k, n) = 0.0; // Initialize the updated state eqn term to zero.

            // Add advection terms.
            if (solverChoice.use_state_advection)
                cell_rhs(i, j, k, n) += -AdvectionContributionForState(i, j, k, rho_u, rho_v, rho_w, cell_data, n, dx, solverChoice.spatial_order);

            // Add diffusive terms.
            if (solverChoice.use_thermal_diffusion && n == RhoTheta_comp)
                cell_rhs(i, j, k, n) += DiffusionContributionForState(i, j, k,cell_data, RhoTheta_comp, dx, solverChoice);
            if (solverChoice.use_scalar_diffusion && n == RhoScalar_comp)
                cell_rhs(i, j, k, n) += DiffusionContributionForState(i, j, k,cell_data, RhoScalar_comp, dx, solverChoice);

            // Add source terms. TODO: Put this under a if condition when we implement source term
            cell_rhs(i, j, k, n) += source_fab(i, j, k, n);
        }
        );
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

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_u_rhs(i, j, k) += -AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::x, dx, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_u_rhs(i, j, k) += DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::x, dx, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
            {
                rho_u_rhs(i, j, k) += (-1.0_rt / dx[0]) *
                  (getPprimegivenRTh(cell_data(i    , j, k, RhoTheta_comp),dptr_pres_hse[k]) -
                   getPprimegivenRTh(cell_data(i - 1, j, k, RhoTheta_comp),dptr_pres_hse[k]));
            }

            // Add gravity term
            if (solverChoice.use_gravity)
                rho_u_rhs(i, j, k) += grav_gpu[0] *
                  InterpolateDensityPertFromCellToFace(i, j, k, cell_data, NextOrPrev::prev,
                                                       Coord::x, solverChoice.spatial_order, dptr_dens_hse);

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_u_rhs(i, j, k) += -solverChoice.abl_pressure_grad[0];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_v_loc = 0.25 * (rho_u(i,j+1,k) + rho_u(i,j,k) + rho_u(i-1,j+1,k) + rho_u(i-1,j,k));
                Real rho_w_loc = 0.25 * (rho_w(i,j,k+1) + rho_w(i,j,k) + rho_w(i,j-1,k+1) + rho_w(i,j-1,k));
                rho_u_rhs(i, j, k) += solverChoice.coriolis_factor *
                        (rho_v_loc * solverChoice.sinphi - rho_w_loc * solverChoice.cosphi);
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_u_rhs(i, j, k) += solverChoice.abl_geo_forcing[0];

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

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_v_rhs(i, j, k) += -AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::y, dx, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_v_rhs(i, j, k) += DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::y, dx, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_v_rhs(i, j, k) += (-1.0_rt / dx[1]) *
                  (getPprimegivenRTh(cell_data(i, j    , k, RhoTheta_comp),dptr_pres_hse[k]) -
                   getPprimegivenRTh(cell_data(i, j - 1, k, RhoTheta_comp),dptr_pres_hse[k]));

            // Add gravity term
            if (solverChoice.use_gravity)
               rho_v_rhs(i, j, k) += grav_gpu[1] *
                  InterpolateDensityPertFromCellToFace(i, j, k, cell_data, NextOrPrev::prev,
                                                       Coord::y, solverChoice.spatial_order, dptr_dens_hse);

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_v_rhs(i, j, k) += -solverChoice.abl_pressure_grad[1];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j-1,k) + rho_u(i,j-1,k));
                rho_v_rhs(i, j, k) += -solverChoice.coriolis_factor * rho_u_loc * solverChoice.sinphi;
            }

            // Add geostrophic forcing
            if (solverChoice.abl_driver_type == ABLDriverType::GeostrophicWind)
                rho_v_rhs(i, j, k) += solverChoice.abl_geo_forcing[1];

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

            // Add advective terms
            if (solverChoice.use_momentum_advection)
                rho_w_rhs(i, j, k) += -AdvectionContributionForMom(i, j, k, rho_u, rho_v, rho_w, u, v, w, MomentumEqn::z, dx, solverChoice);

            // Add diffusive terms
            if (solverChoice.use_momentum_diffusion)
                rho_w_rhs(i, j, k) += DiffusionContributionForMom(i, j, k, u, v, w, MomentumEqn::z, dx, nut, solverChoice);

            // Add pressure gradient
            if (solverChoice.use_pressure)
                rho_w_rhs(i, j, k) += (-1.0_rt / dx[2]) *
                    (getPprimegivenRTh(cell_data(i, j, k    , RhoTheta_comp),dptr_pres_hse[k  ]) -
                     getPprimegivenRTh(cell_data(i, j, k - 1, RhoTheta_comp),dptr_pres_hse[k-1]));

            // Add gravity term
            if (solverChoice.use_gravity)
               rho_w_rhs(i, j, k) += grav_gpu[2] *
                   InterpolateDensityPertFromCellToFace(i, j, k, cell_data, NextOrPrev::prev,
                                                       Coord::z, solverChoice.spatial_order, dptr_dens_hse);

            // Add driving pressure gradient
            if (solverChoice.abl_driver_type == ABLDriverType::PressureGradient)
                rho_w_rhs(i, j, k) += -solverChoice.abl_pressure_grad[2];

            // Add Coriolis forcing (that assumes east is +x, north is +y)
            if (solverChoice.use_coriolis)
            {
                Real rho_u_loc = 0.25 * (rho_u(i+1,j,k) + rho_u(i,j,k) + rho_u(i+1,j,k-1) + rho_u(i,j,k-1));
                rho_w_rhs(i, j, k) += solverChoice.coriolis_factor * rho_u_loc * solverChoice.cosphi;
            }
            } // not on coarse-fine boundary
        });
    }
}
