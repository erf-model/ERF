#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>

#include <Constants.H>

#include <ERF.H>
#include <RK3.H>

using namespace amrex;

void RK3_stage  (MultiFab& cons_old,  MultiFab& cons_upd,
                 MultiFab& xmom_old, MultiFab& ymom_old, MultiFab& zmom_old, 
                 MultiFab& xmom_upd, MultiFab& ymom_upd, MultiFab& zmom_upd, 
                 MultiFab& xvel    , MultiFab& yvel    , MultiFab& zvel    ,  
                 MultiFab& prim    , MultiFab& source,
                 MultiFab& eta, MultiFab& zeta, MultiFab& kappa, MultiFab& alpha,
                 std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, AMREX_SPACEDIM>& facegrad_scalar,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom, const amrex::Real* dxp, const amrex::Real dt)
{
    BL_PROFILE_VAR("RK3_stage()",RK3_stage);

    int nvars = cons_old.nComp(); 

    // ************************************************************************************ 
    // 
    // Fill the ghost cells/faces of the MultiFabs we will need
    // 
    // ************************************************************************************ 
    cons_old.FillBoundary(geom.periodicity());

    xvel.FillBoundary(geom.periodicity());
    yvel.FillBoundary(geom.periodicity());
    zvel.FillBoundary(geom.periodicity());

    // ************************************************************************************** 

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    // const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, CONST_GRAV};
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, 0.0};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    // ************************************************************************************** 
    // 
    // Calculate face-based fluxes to update cell-centered quantities, and 
    //           edge-based and cell-based fluxes to update face-centered quantities
    // 
    // ************************************************************************************** 
    CalcAdvFlux(cons_old, xmom_old, ymom_old, zmom_old, 
                xvel    , yvel    , zvel    , 
                faceflux, facegrad_scalar, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                geom, dxp, dt);
#if 0
    CalcDiffFlux(cons_old, xmom_old, ymom_old, zmom_old, 
                 xvel    , yvel    , zvel    , 
                 eta, zeta, kappa,
                 faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                 geom, dxp, dt);
#endif

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

        const Array4<Real> & cu_old     = cons_old.array(mfi);
        const Array4<Real> & cu_upd     = cons_upd.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real>& momx = xmom_old.array(mfi);
        const Array4<Real>& momy = ymom_old.array(mfi);
        const Array4<Real>& momz = zmom_old.array(mfi);

        const Array4<Real>& mompx = xmom_upd.array(mfi);
        const Array4<Real>& mompy = ymom_upd.array(mfi);
        const Array4<Real>& mompz = zmom_upd.array(mfi);

        Array4<Real const> const& xflux = faceflux[0].array(mfi);
        Array4<Real const> const& yflux = faceflux[1].array(mfi);
        Array4<Real const> const& zflux = faceflux[2].array(mfi);

        Array4<Real const> const& xgrad_scalar = facegrad_scalar[0].array(mfi);
        Array4<Real const> const& ygrad_scalar = facegrad_scalar[1].array(mfi);
        Array4<Real const> const& zgrad_scalar = facegrad_scalar[2].array(mfi);

        Array4<Real const> const& edgex_v = edgeflux_x[0].array(mfi);
        Array4<Real const> const& edgex_w = edgeflux_x[1].array(mfi);

        Array4<Real const> const& edgey_u = edgeflux_y[0].array(mfi);
        Array4<Real const> const& edgey_w = edgeflux_y[1].array(mfi);

        Array4<Real const> const& edgez_u = edgeflux_z[0].array(mfi);
        Array4<Real const> const& edgez_v = edgeflux_z[1].array(mfi);

        Array4<Real const> const& cenx_u = cenflux[0].array(mfi);
        Array4<Real const> const& ceny_v = cenflux[1].array(mfi);
        Array4<Real const> const& cenz_w = cenflux[2].array(mfi);

        const Array4<const Real> alpha_s = alpha.array(mfi);
        const Array4<const Real> kappa_t = kappa.array(mfi);
        amrex::ParallelFor(bx, nvars, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            //could unwarp the n loop to accelarate
            Real coef = (n==0) ? 0.0 : (n==1) ? kappa_t(i,j,k)/(cu_old(i,j,k,0)*c_p) : alpha_s(i,j,k);

            cu_upd(i,j,k,n) = 
                - dt * 
                    ( (xflux(i+1,j,k,n) - xflux(i,j,k,n)) / dx[0]
                    + (yflux(i,j+1,k,n) - yflux(i,j,k,n)) / dx[1]
                    + (zflux(i,j,k+1,n) - zflux(i,j,k,n)) / dx[2] )
                + dt * coef *
                    ( (xgrad_scalar(i+1,j,k,n) - xgrad_scalar(i,j,k,n)) / dx[0]
                    + (ygrad_scalar(i,j+1,k,n) - ygrad_scalar(i,j,k,n)) / dx[1]
                    + (zgrad_scalar(i,j,k+1,n) - zgrad_scalar(i,j,k,n)) / dx[2] )
                + dt * source_fab(i,j,k,n);
        });

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompx(i,j,k) = 
                    -dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                    -dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                    -dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                    +0.5*dt*grav_gpu[0]*(cu_old(i-1,j,k,0)+cu_old(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompy(i,j,k) = 
                    -dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                    -dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                    -dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                    +0.5*dt*grav_gpu[1]*(cu_old(i,j-1,k,0)+cu_old(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            mompz(i,j,k) = 
                    -dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                    -dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                    -dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2] 
                    +0.5*dt*grav_gpu[2]*(cu_old(i,j,k-1,0)+cu_old(i,j,k,0));
        });
    }
}
