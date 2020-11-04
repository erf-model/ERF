#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_VisMF.H>

#include <RK3.H>

using namespace amrex;

void RK3stepStag(MultiFab& cu, 
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                 MultiFab& prim, std::array< MultiFab, AMREX_SPACEDIM >& vel,
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

    MultiFab cup (cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    MultiFab cup2(cu.boxArray(),cu.DistributionMap(),nvars,ngc);
    cup.setVal(0.0,0,nvars,ngc);
    cup2.setVal(0.0,0,nvars,ngc);

    Real rho0 = 1.0;
    cup.setVal(rho0,0,1,ngc);
    cup2.setVal(rho0,0,1,ngc);

    Vector<amrex::IntVect> nodal_flag_dir;
    nodal_flag_dir.resize(AMREX_SPACEDIM);

    IntVect nodal_flag_x;
    IntVect nodal_flag_y;
    IntVect nodal_flag_z;

    std::array< MultiFab, AMREX_SPACEDIM > cupmom;
    std::array< MultiFab, AMREX_SPACEDIM > cup2mom;
    for (int d=0; d<AMREX_SPACEDIM; d++) 
    {
        // Designates data on faces
        nodal_flag_x[d] = int(d==0);
        nodal_flag_y[d] = int(d==1);
        nodal_flag_z[d] = int(d==2);

        // Enable indexing flags above in loops
        nodal_flag_dir[0][d] = nodal_flag_x[d];
        nodal_flag_dir[1][d] = nodal_flag_y[d];
        nodal_flag_dir[2][d] = nodal_flag_z[d];
    }

    for (int d=0; d<AMREX_SPACEDIM; d++) 
    {
        cupmom[d].define(convert(cu.boxArray(),nodal_flag_dir[d]), cu.DistributionMap(), 1, 1);
        cupmom[d].setVal(0.);
        cup2mom[d].define(convert(cu.boxArray(),nodal_flag_dir[d]), cu.DistributionMap(), 1, 1);
        cup2mom[d].setVal(0.);
    }

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    
    const    Array<Real,AMREX_SPACEDIM> grav{0.0, 0.0, -9.8};
    const GpuArray<Real,AMREX_SPACEDIM> grav_gpu{grav[0], grav[1], grav[2]};

    /////////////////////////////////////////////////////

    // Set transport coefs 
    eta.setVal(0.);
    zeta.setVal(0.);
    kappa.setVal(0.);

    calculateFluxStag(cu, cumom, prim, vel, eta, zeta, kappa,
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

        AMREX_D_TERM(const Array4<Real>& momx = cumom[0].array(mfi);,
                     const Array4<Real>& momy = cumom[1].array(mfi);,
                     const Array4<Real>& momz = cumom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& mompx = cupmom[0].array(mfi);,
                     const Array4<Real>& mompy = cupmom[1].array(mfi);,
                     const Array4<Real>& mompz = cupmom[2].array(mfi););

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

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cupmom[d].FillBoundary(geom.periodicity());
    }
    cup.FillBoundary(geom.periodicity());
    conservedToPrimitiveStag(prim, vel, cup, cupmom);
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
    // setBC(prim, cup);

    // Set transport coeffs
    eta.setVal(0.);
    zeta.setVal(0.);
    kappa.setVal(0.);

    ///////////////////////////////////////////////////////////

    calculateFluxStag(cup, cupmom, prim, vel, eta, zeta, kappa,
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
        const Array4<Real> & cup2_fab = cup2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        AMREX_D_TERM(const Array4<Real>& momx = cumom[0].array(mfi);,
                     const Array4<Real>& momy = cumom[1].array(mfi);,
                     const Array4<Real>& momz = cumom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& mompx = cupmom[0].array(mfi);,
                     const Array4<Real>& mompy = cupmom[1].array(mfi);,
                     const Array4<Real>& mompz = cupmom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = cup2mom[0].array(mfi);,
                     const Array4<Real>& momp2y = cup2mom[1].array(mfi);,
                     const Array4<Real>& momp2z = cup2mom[2].array(mfi););

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
        
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cup2mom[d].FillBoundary(geom.periodicity());
    }
    cup2.FillBoundary(geom.periodicity());
    conservedToPrimitiveStag(prim, vel, cup2, cup2mom);
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
    // setBC(prim, cup2);

    // Set transport coeffs
    eta.setVal(0.);
    zeta.setVal(0.);
    kappa.setVal(0.);

    ///////////////////////////////////////////////////////////

    calculateFluxStag(cup2, cup2mom, prim, vel, eta, zeta, kappa,
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

        AMREX_D_TERM(const Array4<Real>& momx = cumom[0].array(mfi);,
                     const Array4<Real>& momy = cumom[1].array(mfi);,
                     const Array4<Real>& momz = cumom[2].array(mfi););

        AMREX_D_TERM(const Array4<Real>& momp2x = cup2mom[0].array(mfi);,
                     const Array4<Real>& momp2y = cup2mom[1].array(mfi);,
                     const Array4<Real>& momp2z = cup2mom[2].array(mfi););

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
            cu_fab(i,j,k,n) = (2./3.) *( 0.5* cu_fab(i,j,k,n) + cup2_fab(i,j,k,n) - dt *
                                    (   AMREX_D_TERM(  (xflux_fab(i+1,j,k,n) - xflux_fab(i,j,k,n)) / dx[0],
                                                     + (yflux_fab(i,j+1,k,n) - yflux_fab(i,j,k,n)) / dx[1],
                                                     + (zflux_fab(i,j,k+1,n) - zflux_fab(i,j,k,n)) / dx[2]) 
                                                                                                            )
                                    + dt*source_fab(i,j,k,n) );
            
        }); // [1:3 indices are not valuable -- momentum flux]

        // momentum flux
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momx(i,j,k) = (2./3.)*(0.5*momx(i,j,k) + momp2x(i,j,k))
                  -(2./3.)*dt*(cenx_u(i,j,k) - cenx_u(i-1,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_u(i,j+1,k) - edgey_u(i,j,k))/dx[1]
                  -(2./3.)*dt*(edgez_u(i,j,k+1) - edgez_u(i,j,k))/dx[2]
                  +0.5*(2./3.)*dt*grav_gpu[0]*(cup2_fab(i-1,j,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momy(i,j,k) = (2./3.)*(0.5*momy(i,j,k) + momp2y(i,j,k))
                  -(2./3.)*dt*(edgex_v(i+1,j,k) - edgex_v(i,j,k))/dx[0]
                  -(2./3.)*dt*(ceny_v(i,j,k) - ceny_v(i,j-1,k))/dx[1]
                  -(2./3.)*dt*(edgez_v(i,j,k+1) - edgez_v(i,j,k))/dx[2]
                  +0.5*(2/3.)*dt*grav_gpu[1]*(cup2_fab(i,j-1,k,0)+cup2_fab(i,j,k,0));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            momz(i,j,k) = (2./3.)*(0.5*momz(i,j,k) + momp2z(i,j,k))
                  -(2./3.)*dt*(edgex_w(i+1,j,k) - edgex_w(i,j,k))/dx[0]
                  -(2./3.)*dt*(edgey_w(i,j+1,k) - edgey_w(i,j,k))/dx[1]
                  -(2./3.)*dt*(cenz_w(i,j,k) - cenz_w(i,j,k-1))/dx[2]
                  +0.5*(2./3.)*dt*grav_gpu[2]*(cup2_fab(i,j,k-1,0)+cup2_fab(i,j,k,0));
        });
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cu_fab(i,j,k,4) += 0.5 * (2./3.) * dt * (  grav_gpu[0]*(momp2x(i+1,j,k)+momp2x(i,j,k))
                                                    + grav_gpu[1]*(momp2y(i,j+1,k)+momp2y(i,j,k))
                                                    + grav_gpu[2]*(momp2z(i,j,k+1)+momp2z(i,j,k)) );
        });
    }

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cumom[d].FillBoundary(geom.periodicity());
    }
    cu.FillBoundary(geom.periodicity());
    conservedToPrimitiveStag(prim, vel, cu, cumom);
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
    // setBC(prim, cu);
}
