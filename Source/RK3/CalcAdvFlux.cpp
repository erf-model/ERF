#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include <EOS.H>
#include <RK3.H>

using namespace amrex;

void CalcAdvFlux(const MultiFab& cons_in, 
                 const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in, 
                 const MultiFab& xvel_in, const MultiFab& yvel_in, const MultiFab& zvel_in, 
                 std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                 std::array< MultiFab, 2 >& edgeflux_x,
                 std::array< MultiFab, 2 >& edgeflux_y,
                 std::array< MultiFab, 2 >& edgeflux_z,
                 std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                 const amrex::Geometry geom,
                 const amrex::Real* dx, const amrex::Real dt,
                 const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("CalcAdvFlux()",CalcAdvFlux);

    GpuArray<Real,AMREX_SPACEDIM> dx_gpu;
    for (int n=0; n<AMREX_SPACEDIM; ++n) {
        dx_gpu[n] = dx[n];
    }
    
    faceflux[0].setVal(0.0);
    faceflux[1].setVal(0.0);
    faceflux[2].setVal(0.0);

    edgeflux_x[0].setVal(0.0);
    edgeflux_x[1].setVal(0.0);

    edgeflux_y[0].setVal(0.0);
    edgeflux_y[1].setVal(0.0);

    edgeflux_z[0].setVal(0.0);
    edgeflux_z[1].setVal(0.0);

    cenflux[0].setVal(0.0);
    cenflux[1].setVal(0.0);
    cenflux[2].setVal(0.0);

    // Loop over boxes
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) 
    {
        const Box & bx_xy = mfi.tilebox(IntVect(1,1,0));
        const Box & bx_xz = mfi.tilebox(IntVect(1,0,1));
        const Box & bx_yz = mfi.tilebox(IntVect(0,1,1));

        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& yflux = faceflux[1].array(mfi);
        const Array4<Real>& zflux = faceflux[2].array(mfi);

        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);
        const Array4<Real>& edgey_u = edgeflux_y[0].array(mfi);
        const Array4<Real>& edgey_w = edgeflux_y[1].array(mfi);
        const Array4<Real>& edgez_u = edgeflux_z[0].array(mfi);
        const Array4<Real>& edgez_v = edgeflux_z[1].array(mfi);

        const Array4<Real>& cenx_u = cenflux[0].array(mfi);
        const Array4<Real>& ceny_v = cenflux[1].array(mfi);
        const Array4<Real>& cenz_w = cenflux[2].array(mfi);

        Array4<Real const> const& momx = xmom_in.array(mfi);
        Array4<Real const> const& momy = ymom_in.array(mfi);
        Array4<Real const> const& momz = zmom_in.array(mfi);

        Array4<Real const> const& velx = xvel_in.array(mfi);
        Array4<Real const> const& vely = yvel_in.array(mfi);
        Array4<Real const> const& velz = zvel_in.array(mfi);

        const Array4<const Real> cons = cons_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Box& bx = mfi.tilebox();

        // Loop over the cells and compute fluxes
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // Define average values on faces
            Real rho      = 0.5 * (cons(i,j,k, Density_comp) + cons(i-1,j,k, Density_comp));
            Real rhotheta = 0.5 * (cons(i,j,k,RhoTheta_comp) + cons(i-1,j,k,RhoTheta_comp));
            Real rhoscalar   = 0.5 * (cons(i, j, k, RhoScalar_comp) + cons(i - 1, j, k, RhoScalar_comp));

            // Density flux
            xflux(i,j,k,Density_comp)  = rho * velx(i,j,k);

            // Theta: conservative flux is (rho u theta)
            xflux(i,j,k,RhoTheta_comp) = rhotheta * velx(i,j,k);

            // Scalar: conservative flux is (rho u s)
            xflux(i, j, k, RhoScalar_comp)   = rhoscalar * velx(i, j, k);
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // Define average values on faces
            Real rho      = 0.5 * (cons(i,j,k, Density_comp) + cons(i,j-1,k, Density_comp));
            Real rhotheta = 0.5 * (cons(i,j,k,RhoTheta_comp) + cons(i,j-1,k,RhoTheta_comp));
            Real rhoscalar   = 0.5 * (cons(i, j, k, RhoScalar_comp) + cons(i, j - 1, k, RhoScalar_comp));

            // Density
            yflux(i,j,k,Density_comp)  = rho * vely(i,j,k);

            // Theta: conservative flux is (rho u theta)
            yflux(i,j,k,RhoTheta_comp) = rhotheta * vely(i,j,k);

            // Scalar: conservative flux is (rho u s)
            yflux(i, j, k, RhoScalar_comp)   = rhoscalar * vely(i, j, k);
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // Define average values on faces
            Real rho      = 0.5 * (cons(i,j,k, Density_comp) + cons(i,j,k-1, Density_comp));
            Real rhotheta = 0.5 * (cons(i,j,k,RhoTheta_comp) + cons(i,j,k-1,RhoTheta_comp));
            Real rhoscalar   = 0.5 * (cons(i, j, k, RhoScalar_comp) + cons(i, j, k - 1, RhoScalar_comp));

            // Density
            zflux(i,j,k,Density_comp)  = rho * velz(i,j,k);

            // Theta: conservative flux is (rho u theta)
            zflux(i,j,k,RhoTheta_comp) = rhotheta * velz(i,j,k);

            // Scalar: conservative flux is (rho u s)
            zflux(i, j, k, RhoScalar_comp)   = rhoscalar * velz(i, j, k);
        }
        ); // 'ParallelFor' ends here

            amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgey_u(i,j,k) = 0.25*(momx(i,j-1,k)+momx(i,j,k))*(vely(i-1,j,k)+vely(i,j,k));
                edgex_v(i,j,k) = 0.25*(momy(i-1,j,k)+momy(i,j,k))*(velx(i,j-1,k)+velx(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgez_u(i,j,k) = 0.25*(momx(i,j,k-1)+momx(i,j,k))*(velz(i-1,j,k)+velz(i,j,k));
                edgex_w(i,j,k) = 0.25*(momz(i-1,j,k)+momz(i,j,k))*(velx(i,j,k-1)+velx(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgez_v(i,j,k) = 0.25*(momy(i,j,k-1)+momy(i,j,k))*(velz(i,j-1,k)+velz(i,j,k));
                edgey_w(i,j,k) = 0.25*(momz(i,j-1,k)+momz(i,j,k))*(vely(i,j,k-1)+vely(i,j,k));
            });

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real rhotheta = cons(i,j,k, RhoTheta_comp);
                Real pressure = getPgivenRTh(rhotheta);
                cenx_u(i,j,k) = 0.25*(momx(i,j,k)+momx(i+1,j,k))*(velx(i,j,k)+velx(i+1,j,k)) + pressure;
                ceny_v(i,j,k) = 0.25*(momy(i,j,k)+momy(i,j+1,k))*(vely(i,j,k)+vely(i,j+1,k)) + pressure;
                cenz_w(i,j,k) = 0.25*(momz(i,j,k)+momz(i,j,k+1))*(velz(i,j,k)+velz(i,j,k+1)) + pressure;
            });

    } // end mfi

    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cenflux[d].FillBoundary(geom.periodicity());
        faceflux[d].FillBoundary(geom.periodicity());
    }

    for (int d=0; d<2; d++) {
        edgeflux_x[d].FillBoundary(geom.periodicity());
        edgeflux_y[d].FillBoundary(geom.periodicity());
        edgeflux_z[d].FillBoundary(geom.periodicity());
    }
}
