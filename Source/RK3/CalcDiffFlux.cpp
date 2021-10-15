#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include <RK3.H>

using namespace amrex;

void CalcDiffFlux(const MultiFab& cons_in,
                       const MultiFab& xmom_in, const MultiFab& ymom_in, const MultiFab& zmom_in,
                       const MultiFab& prim_in,
                       const MultiFab& xvel_in, const MultiFab& yvel_in, const MultiFab& zvel_in,
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                       std::array< MultiFab, 2 >& edgeflux_x,
                       std::array< MultiFab, 2 >& edgeflux_y,
                       std::array< MultiFab, 2 >& edgeflux_z,
                       std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                       const amrex::Geometry geom,
                       const amrex::Real* dx, const amrex::Real dt,
                       const SolverChoice& solverChoice)
{
    BL_PROFILE_VAR("CalcDiffFlux()",CalcDiffFlux);

//    int visc_type_gpu = 3;  // Include bulk viscosity

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

    int ngc = 1;

    std::array< MultiFab, AMREX_SPACEDIM > tau_diag; // diagonal stress (defined at cell centers)
    tau_diag[0].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);
    tau_diag[1].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);
    tau_diag[2].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);

    std::array< MultiFab, AMREX_SPACEDIM > tau_diagoff; // off diagonal stress at edges
    tau_diagoff[0].define(convert(cons_in.boxArray(),IntVect(1,1,0)),cons_in.DistributionMap(),1,ngc);
    tau_diagoff[1].define(convert(cons_in.boxArray(),IntVect(0,1,1)),cons_in.DistributionMap(),1,ngc);
    tau_diagoff[2].define(convert(cons_in.boxArray(),IntVect(1,0,1)),cons_in.DistributionMap(),1,ngc);

    tau_diag[0].setVal(0.0);
    tau_diag[1].setVal(0.0);
    tau_diag[2].setVal(0.0);

    tau_diagoff[0].setVal(0.0);
    tau_diagoff[1].setVal(0.0);
    tau_diagoff[2].setVal(0.0);

    ////////////////////
    // diffusive fluxes
    ////////////////////

    /*
    // Loop over boxes
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        const Box & bx_xy = mfi.tilebox(IntVect(1,1,0));
        const Box & bx_xz = mfi.tilebox(IntVect(1,0,1));
        const Box & bx_yz = mfi.tilebox(IntVect(0,1,1));

        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& yflux = faceflux[1].array(mfi);
        const Array4<Real>& zflux = faceflux[2].array(mfi);

        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);
        const Array4<Real>& edgey_u = edgeflux_y[0].array(mfi);
        const Array4<Real>& edgey_w = edgeflux_y_in[1].array(mfi);
        const Array4<Real>& edgez_u = edgeflux_z_in[0].array(mfi);
        const Array4<Real>& edgez_v = edgeflux_z_in[1].array(mfi);

        const Array4<Real>& cenx_u = cenflux[0].array(mfi);
        const Array4<Real>& ceny_v = cenflux[1].array(mfi);
        const Array4<Real>& cenz_w = cenflux[2].array(mfi);

        const Array4<Real> tauxx = tau_diag[0].array(mfi);
        const Array4<Real> tauyy = tau_diag[1].array(mfi);
        const Array4<Real> tauzz = tau_diag[2].array(mfi);

        const Array4<Real> tauxy = tau_diagoff[0].array(mfi);
        const Array4<Real> tauyz = tau_diagoff[1].array(mfi);
        const Array4<Real> tauxz = tau_diagoff[2].array(mfi);

        Array4<Real const> const& momx = cu_x.array(mfi);
        Array4<Real const> const& momy = cu_y.array(mfi);
        Array4<Real const> const& momz = cu_z.array(mfi);

        Array4<Real const> const& velx = u_x.array(mfi);
        Array4<Real const> const& vely = v_y.array(mfi);
        Array4<Real const> const& velz = w_z.array(mfi);

        const Array4<const Real> prim = prim_in.array(mfi);
        const Array4<const Real> cons = cons_in.array(mfi);

        const Array4<const Real> eta   = eta_in.array(mfi);
        const Array4<const Real> zeta  = zeta_in.array(mfi);
        const Array4<const Real> kappa = kappa_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        IntVect nd(AMREX_D_DECL(1,1,1));
        const Box& tbn = mfi.tilebox(nd);

        const Box& bx = mfi.tilebox();

        Real half = 0.5;

        // Populate diagonal stress
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_x, v_y, w_z; // velocity gradients
            u_x = (velx(i+1,j,k) - velx(i,j,k))/dx_gpu[0];
            v_y = (vely(i,j+1,k) - vely(i,j,k))/dx_gpu[1];
            w_z = (velz(i,j,k+1) - velz(i,j,k))/dx_gpu[2];

            Real div = u_x + v_y + w_z;
            if (amrex::Math::abs(visc_type_gpu) == 3) {
              tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
              tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
              tauzz(i,j,k) = 2*eta(i,j,k)*w_z + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
            }
            else {
              tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (0.0 - 2*eta(i,j,k)/3.)*div;
              tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (0.0 - 2*eta(i,j,k)/3.)*div;
              tauzz(i,j,k) = 2*eta(i,j,k)*w_z + (0.0 - 2*eta(i,j,k)/3.)*div;
            }

        });

        // Populate off-diagonal stress
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_y, v_x; // velocity gradients
            u_y = (velx(i,j,k) - velx(i,j-1,k))/dx_gpu[1];
            v_x = (vely(i,j,k) - vely(i-1,j,k))/dx_gpu[0];
            tauxy(i,j,k) = 0.25*(eta(i-1,j-1,k)+eta(i-1,j,k)+eta(i,j-1,k)+eta(i,j,k))*(u_y+v_x);

        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_z, w_x; // velocity gradients
            u_z = (velx(i,j,k) - velx(i,j,k-1))/dx_gpu[2];
            w_x = (velz(i,j,k) - velz(i-1,j,k))/dx_gpu[0];
            tauxz(i,j,k) = 0.25*(eta(i-1,j,k-1)+eta(i-1,j,k)+eta(i,j,k-1)+eta(i,j,k))*(u_z+w_x);
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real v_z, w_y; // velocity gradients
            v_z = (vely(i,j,k) - vely(i,j,k-1))/dx_gpu[2];
            w_y = (velz(i,j,k) - velz(i,j-1,k))/dx_gpu[1];
            tauyz(i,j,k) = 0.25*(eta(i,j-1,k-1)+eta(i,j-1,k)+eta(i,j,k-1)+eta(i,j,k))*(v_z+w_y);

        });

        // Loop over faces for flux calculations (4:5+ns)
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            xflux(i,j,k,4) -= 0.5*velx(i,j,k)*(tauxx(i-1,j,k)+tauxx(i,j,k));
            xflux(i,j,k,4) -= 0.25*((vely(i,j+1,k)+vely(i-1,j+1,k))*tauxy(i,j+1,k) + (vely(i,j,k)+vely(i-1,j,k))*tauxy(i,j,k));
            xflux(i,j,k,4) -= 0.25*((velz(i,j,k+1)+velz(i-1,j,k+1))*tauxz(i,j,k+1) + (velz(i,j,k)+velz(i-1,j,k))*tauxz(i,j,k));

            Real kxp = 0.5*(kappa(i,j,k) + kappa(i-1,j,k));
            xflux(i,j,k,4) -= kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/dx_gpu[0];
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            yflux(i,j,k,4) -= 0.25*((velx(i+1,j,k)+velx(i+1,j-1,k))*tauxy(i+1,j,k) + (velx(i,j,k)+velx(i,j-1,k))*tauxy(i,j,k));
            yflux(i,j,k,4) -= 0.5*vely(i,j,k)*(tauyy(i,j-1,k)+tauyy(i,j,k));
            yflux(i,j,k,4) -= 0.25*((velz(i,j,k+1)+velz(i,j-1,k+1))*tauyz(i,j,k+1) + (velz(i,j,k)+velz(i,j-1,k))*tauyz(i,j,k));

            Real kyp = 0.5*(kappa(i,j,k) + kappa(i,j-1,k));
            yflux(i,j,k,4) -= kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/dx_gpu[1];
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            zflux(i,j,k,4) -= 0.25*((velx(i+1,j,k-1)+velx(i+1,j,k))*tauxz(i+1,j,k) + (velx(i,j,k)+velx(i,j,k-1))*tauxz(i,j,k));
            zflux(i,j,k,4) -= 0.25*((vely(i,j+1,k-1)+vely(i,j+1,k))*tauyz(i,j+1,k) + (vely(i,j,k)+vely(i,j,k-1))*tauyz(i,j,k));
            zflux(i,j,k,4) -= 0.5*velz(i,j,k)*(tauzz(i,j,k-1)+tauzz(i,j,k));

            Real kzp = 0.5*(kappa(i,j,k) + kappa(i,j,k-1));
            zflux(i,j,k,4) -= kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/dx_gpu[2];
        });

        // Loop over edges for momemntum flux calculations [1:3]
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgey_u(i,j,k) -= tauxy(i,j,k);
            edgex_v(i,j,k) -= tauxy(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgez_u(i,j,k) -= tauxz(i,j,k);
            edgex_w(i,j,k) -= tauxz(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgez_v(i,j,k) -= tauyz(i,j,k);
            edgey_w(i,j,k) -= tauyz(i,j,k);
        });

        // Loop over the center cells and compute fluxes (diagonal momentum terms)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            cenx_u(i,j,k) -= tauxx(i,j,k);
            ceny_v(i,j,k) -= tauyy(i,j,k);
            cenz_w(i,j,k) -= tauzz(i,j,k);
        });
    }
    */

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
