#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>

#include <ERF.H>
#include <RK3.H>

using namespace amrex;

void RK3_advance(MultiFab& cons_old,  MultiFab& cons_new,
                 MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old, 
                 MultiFab& xmom_new, MultiFab& ymom_new, MultiFab& zmom_new, 
                 MultiFab& xvel    , MultiFab& yvel    , MultiFab& zvel    , 
                 MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt)
{
    BL_PROFILE_VAR("RK3stepStag()",RK3stepStag);

    int nvars = cons_old.nComp(); 
    int ngc = 1;

    // ************************************************************************************** 
    // 
    // Allocate temporary MultiFab to hold the primitive variables
    // 
    // ************************************************************************************** 
    MultiFab prim(cons_old.boxArray(),cons_old.DistributionMap(),nvars,2); 

    // ************************************************************************************** 
    // 
    // Allocate temporary MultiFabs to hold the intermediate updates
    // 
    // ************************************************************************************** 

    const BoxArray& ba            = cons_old.boxArray();
    const DistributionMapping& dm = cons_old.DistributionMap();

    // At cell centers
    MultiFab cons_update_1(ba,dm,nvars,ngc);
    MultiFab cons_update_2(ba,dm,nvars,ngc);
    cons_update_1.setVal(0.0,0,nvars,ngc);
    cons_update_2.setVal(0.0,0,nvars,ngc);

    // On faces
    MultiFab xmom_update_1(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab ymom_update_1(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab zmom_update_1(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    MultiFab xmom_update_2(convert(ba,IntVect(1,0,0)), dm, 1, 1);
    MultiFab ymom_update_2(convert(ba,IntVect(0,1,0)), dm, 1, 1);
    MultiFab zmom_update_2(convert(ba,IntVect(0,0,1)), dm, 1, 1);

    // ************************************************************************************** 
    // 
    // Convert old conservative variables cons_old --> prim
    //     and old face-based momentum --> velocity
    // 
    // ************************************************************************************** 
    ConservedToPrimitive(prim, xvel, yvel, zvel, cons_old, xmom_old, ymom_old, zmom_old);

    // ************************************************************************************ 
    // 
    // Fill the ghost cells/faces of the MultiFabs we have just filled
    // 
    // ************************************************************************************ 
    prim.FillBoundary(geom.periodicity());

    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // ************************************************************************************** 

    Real rho0 = 1.0;
    cons_update_1.setVal(rho0,0,1,ngc);
    cons_update_2.setVal(rho0,0,1,ngc);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, 0.0};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // ************************************************************************************** 
    // 
    // Calculate face-based fluxes to update cell-centered quantities, and 
    //           edge-based and cell-based fluxes to update face-centered quantities
    // 
    // ************************************************************************************** 
    calculateFluxStag(cons_old, xmom_old, ymom_old, zmom_old, 
                      prim    , xvel    , yvel    , zvel    , 
                      eta, zeta, kappa,
                      faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                      geom, dxp, dt);

    // ************************************************************************************** 
    // 
    // Define updates in the first RK stage
    // 
    // ************************************************************************************** 
    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab     = cons_old.array(mfi);
        const Array4<Real> & cup_fab    = cons_update_1.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real>& momx = xmom_old.array(mfi);
        const Array4<Real>& momy = ymom_old.array(mfi);
        const Array4<Real>& momz = zmom_old.array(mfi);

        const Array4<Real>& mompx = xmom_update_1.array(mfi);
        const Array4<Real>& mompy = ymom_update_1.array(mfi);
        const Array4<Real>& mompz = zmom_update_1.array(mfi);

        Array4<Real const> const& xflux_fab = faceflux[0].array(mfi);
        Array4<Real const> const& yflux_fab = faceflux[1].array(mfi);
        Array4<Real const> const& zflux_fab = faceflux[2].array(mfi);

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);

        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);

        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup_fab(i,j,k,n) = cu_fab(i,j,k,n) - dt *
                ( AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                               + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                               + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                       )
                + dt*source_fab(i,j,k,n);
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompx(i,j,k) = momx(i,j,k)
                    -dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                    -dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                    -dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                    +0.5*dt*grav_gpu[0]*(cu_fab(i-1,j,k,0)+cu_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompy(i,j,k) = momy(i,j,k)
                    -dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                    -dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*dt*grav_gpu[1]*(cu_fab(i,j-1,k,0)+cu_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompz(i,j,k) = momz(i,j,k)
                    -dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2] 
                    +0.5*dt*grav_gpu[2]*(cu_fab(i,j,k-1,0)+cu_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup_fab(i,j,k,4) += 0.5 * dt * (  grav_gpu[0]*(momx(i+1,j,k)+momx(i,j,k))
                                            + grav_gpu[1]*(momy(i,j+1,k)+momy(i,j,k))
                                            + grav_gpu[2]*(momz(i,j,k+1)+momz(i,j,k)) );
        });
    }

    xmom_update_1.FillBoundary(geom.periodicity());
    ymom_update_1.FillBoundary(geom.periodicity());
    zmom_update_1.FillBoundary(geom.periodicity());

    cons_update_1.FillBoundary(geom.periodicity());

    // ************************************************************************************** 
    // 
    // Convert intermediate conservative variables cons_update_1 --> prim
    //     and intermediate face-based momentum                  --> velocity
    // 
    // ************************************************************************************** 
    ConservedToPrimitive(prim, xvel, yvel, zvel, cons_update_1, xmom_update_1, ymom_update_1, zmom_update_1);

    // ************************************************************************************ 
    // 
    // Fill the ghost cells/faces of the MultiFabs we have just filled
    // 
    // ************************************************************************************ 
    prim.FillBoundary(geom.periodicity());

    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // ************************************************************************************** 
    // 
    // Calculate face-based fluxes to update cell-centered quantities, and 
    //           edge-based and cell-based fluxes to update face-centered quantities
    // 
    // ************************************************************************************** 
    calculateFluxStag(cons_update_1, xmom_update_1, ymom_update_1, zmom_update_1, 
                      prim         , xvel         , yvel         , zvel         , 
                      eta, zeta, kappa,
                      faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                      geom, dxp, dt);

    // ************************************************************************************** 
    // 
    // Define updates in the second RK stage
    // 
    // ************************************************************************************** 

    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab     = cons_old.array(mfi);
        const Array4<Real> & cup_fab    = cons_update_1.array(mfi);
        const Array4<Real> & cup2_fab   = cons_update_2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real>& momx = xmom_old.array(mfi);
        const Array4<Real>& momy = ymom_old.array(mfi);
        const Array4<Real>& momz = zmom_old.array(mfi);

        const Array4<Real>& mompx = xmom_update_1.array(mfi);
        const Array4<Real>& mompy = ymom_update_1.array(mfi);
        const Array4<Real>& mompz = zmom_update_1.array(mfi);

        const Array4<Real>& momp2x = xmom_update_2.array(mfi);
        const Array4<Real>& momp2y = ymom_update_2.array(mfi);
        const Array4<Real>& momp2z = zmom_update_2.array(mfi);

        Array4<Real const> const& xflux_fab = faceflux[0].array(mfi);
        Array4<Real const> const& yflux_fab = faceflux[1].array(mfi);
        Array4<Real const> const& zflux_fab = faceflux[2].array(mfi);

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);
        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);
        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cup2_fab(i,j,k,n) = 0.25*( 3.0* cu_fab(i,j,k,n) + cup_fab(i,j,k,n) - dt *
                                       ( AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                                                      + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                                                      + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2])
                                                                                                                )
                                       +dt*source_fab(i,j,k,n)  );
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2x(i,j,k) = 0.25*3.0*momx(i,j,k) + 0.25*mompx(i,j,k)
                    -0.25*dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                    -0.25*dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                    -0.25*dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                    +0.5*0.25*dt*grav_gpu[0]*(cup_fab(i-1,j,k,0)+cup_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2y(i,j,k) = 0.25*3.0*momy(i,j,k) + 0.25*mompy(i,j,k)
                    -0.25*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -0.25*dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                    -0.25*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*0.25*dt*grav_gpu[1]*(cup_fab(i,j-1,k,0)+cup_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momp2z(i,j,k) = 0.25*3.0*momz(i,j,k) + 0.25*mompz(i,j,k)
                    -0.25*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -0.25*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -0.25*dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
                    +0.5*0.25*dt*grav_gpu[2]*(cup_fab(i,j,k-1,0)+cup_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cup2_fab(i,j,k,4) += 0.5 * 0.25 * dt * (  grav_gpu[0]*(mompx(i+1,j,k)+mompx(i,j,k))
                                                    + grav_gpu[1]*(mompy(i,j+1,k)+mompy(i,j,k))
                                                    + grav_gpu[2]*(mompz(i,j,k+1)+mompz(i,j,k)) );
        });
    }
        
    xmom_update_2.FillBoundary(geom.periodicity());
    ymom_update_2.FillBoundary(geom.periodicity());
    zmom_update_2.FillBoundary(geom.periodicity());

    cons_update_2.FillBoundary(geom.periodicity());

    // ************************************************************************************** 
    // 
    // Convert intermediate conservative variables cons_update_2 --> prim
    //     and intermediate face-based momentum                  --> velocity
    // 
    // ************************************************************************************** 
    ConservedToPrimitive(prim, xvel, yvel, zvel, cons_update_2, xmom_update_2, ymom_update_2, zmom_update_2);

    // ************************************************************************************ 
    // 
    // Fill the ghost cells/faces of the MultiFabs we have just filled
    // 
    // ************************************************************************************ 
    prim.FillBoundary(geom.periodicity());

    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // ************************************************************************************** 
    // 
    // Calculate face-based fluxes to update cell-centered quantities, and 
    //           edge-based and cell-based fluxes to update face-centered quantities
    // 
    // ************************************************************************************** 
    calculateFluxStag(cons_update_2, xmom_update_2, ymom_update_2, zmom_update_2, 
                      prim         , xvel         , yvel         , zvel         , 
                      eta, zeta, kappa,
                      faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                      geom, dxp, dt);

    // ************************************************************************************** 
    // 
    // Define the final update in the third RK stage
    // 
    // ************************************************************************************** 

    for ( MFIter mfi(cons_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab     = cons_old.array(mfi);
        const Array4<Real> & cup2_fab   = cons_update_2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real> & cu_new_fab = cons_new.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = xmom_old.array(mfi);,
                     const Array4<Real>& momy = ymom_old.array(mfi);,
                     const Array4<Real>& momz = zmom_old.array(mfi););

        AMREX_D_TERM(const Array4<Real>& newmomx = xmom_new.array(mfi);,
                     const Array4<Real>& newmomy = ymom_new.array(mfi);,
                     const Array4<Real>& newmomz = zmom_new.array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = xmom_update_2.array(mfi);,
                     const Array4<Real>& momp2y = ymom_update_2.array(mfi);,
                     const Array4<Real>& momp2z = zmom_update_2.array(mfi););

        AMREX_D_TERM(Array4<Real const> const& xflux_fab = faceflux[0].array(mfi);,
                     Array4<Real const> const& yflux_fab = faceflux[1].array(mfi);,
                     Array4<Real const> const& zflux_fab = faceflux[2].array(mfi););

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);
        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);
        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cu_new_fab(i,j,k,n) = (2./3.) * ( 0.5* cu_fab(i,j,k,n) + cup2_fab(i,j,k,n) - dt * (
                                    (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0] 
                                   +(yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1] 
                                   +(zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2] )
                                    + dt*source_fab(i,j,k,n) );
            
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newmomx(i,j,k) = (2./3.)*(0.5*momx(i,j,k) + momp2x(i,j,k))
                  -(2./3.)*dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                  -(2./3.)*dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                  +0.5*(2./3.)*dt*grav_gpu[0]*(cup2_fab(i-1,j,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newmomy(i,j,k) = (2./3.)*(0.5*momy(i,j,k) + momp2y(i,j,k))
                  -(2./3.)*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                  -(2./3.)*dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                  -(2./3.)*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                  +0.5*(2/3.)*dt*grav_gpu[1]*(cup2_fab(i,j-1,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            newmomz(i,j,k) = (2./3.)*(0.5*momz(i,j,k) + momp2z(i,j,k))
                  -(2./3.)*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                  -(2./3.)*dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
                  +0.5*(2./3.)*dt*grav_gpu[2]*(cup2_fab(i,j,k-1,0)+cup2_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cu_new_fab(i,j,k,4) += 0.5 * (2./3.) * dt * ( grav_gpu[0]*(momp2x(i+1,j,k)+momp2x(i,j,k))
                                                        + grav_gpu[1]*(momp2y(i,j+1,k)+momp2y(i,j,k))
                                                        + grav_gpu[2]*(momp2z(i,j,k+1)+momp2z(i,j,k)) );
        });
    }

    xmom_new.FillBoundary(geom.periodicity());
    ymom_new.FillBoundary(geom.periodicity());
    zmom_new.FillBoundary(geom.periodicity());

    cons_new.FillBoundary(geom.periodicity());

    // ************************************************************************************** 
    // 
    // Convert new momentum to new velocity on faces
    // 
    // ************************************************************************************** 
    MomentumToVelocity(xvel, yvel, zvel, cons_new, xmom_new, ymom_new, zmom_new);

    // ************************************************************************************ 
    // 
    // Fill the ghost cells/faces of the MultiFabs we have just filled
    // 
    // ************************************************************************************ 
    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());
}
