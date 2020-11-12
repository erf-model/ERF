#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>

#include <ERF.H>
#include <RK3.H>

using namespace amrex;

void RK3_advance(MultiFab& cu,  MultiFab& cu_new,
                 MultiFab& cu_x, MultiFab& cu_y, MultiFab& cu_z, 
                 MultiFab& cunew_x, MultiFab& cunew_y, MultiFab& cunew_z, 
                 MultiFab& prim, 
                 MultiFab&  u_x, MultiFab&  v_y, MultiFab& w_z, 
                 MultiFab&  unew_x, MultiFab&  vnew_y, MultiFab& wnew_z, 
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

    int nvars = cu.nComp(); 
    int ngc = 1;

    prim.setVal(0.);
    conservedToPrimitiveStag(prim, u_x, v_y, w_z, cu, cu_x, cu_y, cu_z);

    u_x.FillBoundary(geom.periodicity());
    v_y.FillBoundary(geom.periodicity());
    w_z.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    MultiFab cup (cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    MultiFab cup2(cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    cup.setVal(0.0,0,nvars,ngc);
    cup2.setVal(0.0,0,nvars,ngc);

    Real rho0 = 1.0;
    cup.setVal(rho0,0,1,ngc);
    cup2.setVal(rho0,0,1,ngc);

    amrex::Print() << "cu.boxArray() "  << cu.boxArray()  << std::endl;

    MultiFab cupmom_x(convert(cu.boxArray(),IntVect(1,0,0)), cu.DistributionMap(), 1, 1);
    MultiFab cupmom_y(convert(cu.boxArray(),IntVect(0,1,0)), cu.DistributionMap(), 1, 1);
    MultiFab cupmom_z(convert(cu.boxArray(),IntVect(0,0,1)), cu.DistributionMap(), 1, 1);

    MultiFab cup2mom_x(convert(cu.boxArray(),IntVect(1,0,0)), cu.DistributionMap(), 1, 1);
    MultiFab cup2mom_y(convert(cu.boxArray(),IntVect(0,1,0)), cu.DistributionMap(), 1, 1);
    MultiFab cup2mom_z(convert(cu.boxArray(),IntVect(0,0,1)), cu.DistributionMap(), 1, 1);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, 0.0};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    calculateFluxStag(cu, cu_x, cu_y, cu_z, 
                      prim, u_x, v_y, w_z, 
                      eta, zeta, kappa,
                      faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                      geom, dxp, dt);

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cenflux[d].FillBoundary(geom.periodicity());
        faceflux[d].FillBoundary(geom.periodicity());
    }

    for (int d=0; d<2; d++) {
        edgeflux_x[d].FillBoundary(geom.periodicity());
        edgeflux_y[d].FillBoundary(geom.periodicity());
        edgeflux_z[d].FillBoundary(geom.periodicity());
    }

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup_fab = cup.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real>& momx = cu_x.array(mfi);
        const Array4<Real>& momy = cu_y.array(mfi);
        const Array4<Real>& momz = cu_z.array(mfi);

        const Array4<Real>& mompx = cupmom_x.array(mfi);
        const Array4<Real>& mompy = cupmom_y.array(mfi);
        const Array4<Real>& mompz = cupmom_z.array(mfi);

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

    cupmom_x.FillBoundary(geom.periodicity());
    cupmom_y.FillBoundary(geom.periodicity());
    cupmom_z.FillBoundary(geom.periodicity());

    cup.FillBoundary(geom.periodicity());

    conservedToPrimitiveStag(prim, u_x, v_y, w_z, cup, cupmom_x, cupmom_y, cupmom_z);

    u_x.FillBoundary(geom.periodicity());
    v_y.FillBoundary(geom.periodicity());
    w_z.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    // setBC(prim, cup);

    ///////////////////////////////////////////////////////////

    calculateFluxStag(cup, cupmom_x, cupmom_y, cupmom_z, 
                      prim, u_x, v_y, w_z, 
                      eta, zeta, kappa,
                      faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                      geom, dxp, dt);

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cenflux[d].FillBoundary(geom.periodicity());
        faceflux[d].FillBoundary(geom.periodicity());
    }

    for (int d=0; d<2; d++) {
        edgeflux_x[d].FillBoundary(geom.periodicity());
        edgeflux_y[d].FillBoundary(geom.periodicity());
        edgeflux_z[d].FillBoundary(geom.periodicity());
    }

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) 
    {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup_fab = cup.array(mfi);
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = cu_x.array(mfi);,
                     const Array4<Real>& momy = cu_y.array(mfi);,
                     const Array4<Real>& momz = cu_z.array(mfi););

        AMREX_D_TERM(const Array4<Real>& mompx = cupmom_x.array(mfi);,
                     const Array4<Real>& mompy = cupmom_y.array(mfi);,
                     const Array4<Real>& mompz = cupmom_z.array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = cup2mom_x.array(mfi);,
                     const Array4<Real>& momp2y = cup2mom_y.array(mfi);,
                     const Array4<Real>& momp2z = cup2mom_z.array(mfi););

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
        
    cup2mom_x.FillBoundary(geom.periodicity());
    cup2mom_y.FillBoundary(geom.periodicity());
    cup2mom_z.FillBoundary(geom.periodicity());

    cup2.FillBoundary(geom.periodicity());

    conservedToPrimitiveStag(prim, u_x, v_y, w_z, cup2, cu_x, cu_y, cu_z);

    u_x.FillBoundary(geom.periodicity());
    v_y.FillBoundary(geom.periodicity());
    w_z.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

    // setBC(prim, cup2);

    ///////////////////////////////////////////////////////////

    calculateFluxStag(cup2, cup2mom_x, cup2mom_y, cup2mom_z, 
                      prim, u_x, v_y, w_z, 
                      eta, zeta, kappa,
                      faceflux, edgeflux_x, edgeflux_y, edgeflux_z, cenflux, 
                      geom, dxp, dt);


    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cenflux[d].FillBoundary(geom.periodicity());
        faceflux[d].FillBoundary(geom.periodicity());
    }

    for (int d=0; d<2; d++) {
        edgeflux_x[d].FillBoundary(geom.periodicity());
        edgeflux_y[d].FillBoundary(geom.periodicity());
        edgeflux_z[d].FillBoundary(geom.periodicity());
    }

    for ( MFIter mfi(cu,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<Real> & cu_fab = cu.array(mfi);
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        const Array4<Real> & cu_new_fab = cu_new.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = cu_x.array(mfi);,
                     const Array4<Real>& momy = cu_y.array(mfi);,
                     const Array4<Real>& momz = cu_z.array(mfi););

        AMREX_D_TERM(const Array4<Real>& newmomx = cunew_x.array(mfi);,
                     const Array4<Real>& newmomy = cunew_y.array(mfi);,
                     const Array4<Real>& newmomz = cunew_z.array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = cup2mom_x.array(mfi);,
                     const Array4<Real>& momp2y = cup2mom_y.array(mfi);,
                     const Array4<Real>& momp2z = cup2mom_z.array(mfi););

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

    cunew_x.FillBoundary(geom.periodicity());
    cunew_y.FillBoundary(geom.periodicity());
    cunew_z.FillBoundary(geom.periodicity());

    cu_new.FillBoundary(geom.periodicity());

    conservedToPrimitiveStag(prim, unew_x, vnew_y, wnew_z, cu_new, cunew_x, cunew_y, cunew_z); 

    unew_x.FillBoundary(geom.periodicity());
    vnew_y.FillBoundary(geom.periodicity());
    wnew_z.FillBoundary(geom.periodicity()); 

    prim.FillBoundary(geom.periodicity());

    // setBC(prim, cu);
}
