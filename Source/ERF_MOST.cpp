#include <ERF_PhysBCFunct.H>

using namespace amrex;

void
ERFPhysBCFunct::impose_most_bcs(const Box& bx,
                                const Array4<Real>& dest_arr,
                                const Array4<Real>& cons_arr,
                                const Array4<Real>& velx_arr,
                                const Array4<Real>& vely_arr,
                                const Array4<Real>&  eta_arr,
                                const int idx, const int icomp, const int zlo)
{
    // check idx to distinguish between Vars::Dvel and Vars::Dmom
    // (in Legacy this was controlled by the is_derived flag)
    bool var_is_derived = false;
    if (idx == Vars::xvel || idx == Vars::yvel) {
        var_is_derived = true;
    }

    /**
    * NOTE: the number of ghost zone for state variables are different from face centered
    *       variables in the new version.
    */

    if (idx == Vars::cons) {
        amrex::Box b2d = bx; // Copy constructor
        b2d.setBig(2,zlo-1);
        int n = Cons::RhoTheta;
        ParallelFor(b2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                Real velx, vely, rho, theta, eta;
                int ix, jx, iy, jy, ie, je;

                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                jx = j < lbound(velx_arr).y ? lbound(velx_arr).y : j;
                ix = ix > ubound(velx_arr).x-1 ? ubound(velx_arr).x-1 : ix;
                jx = jx > ubound(velx_arr).y ? ubound(velx_arr).y : jx;

                iy = i < lbound(vely_arr).x ? lbound(vely_arr).x : i;
                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                iy = iy > ubound(vely_arr).x ? ubound(vely_arr).x : iy;
                jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
                je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;

                velx  = 0.5*(velx_arr(ix,jx,zlo)+velx_arr(ix+1,jx,zlo));
                vely  = 0.5*(vely_arr(iy,jy,zlo)+vely_arr(iy,jy+1,zlo));
                rho   = cons_arr(ie,je,zlo,Rho_comp);
                theta = cons_arr(ie,je,zlo,RhoTheta_comp)/rho;

                // TODO: Verify this is the correct Diff component
                eta   = eta_arr(ie,je,zlo,EddyDiff::Mom_h);

                Real vmag    = sqrt(velx*velx+vely*vely);
                Real num1    = (theta-m_most.theta_mean)*m_most.vmag_mean;
                Real num2    = (m_most.theta_mean-m_most.surf_temp)*vmag;
                Real motheta = (num1+num2)*m_most.utau*m_most.kappa/m_most.phi_h();

                if (!var_is_derived) {
                    dest_arr(i,j,k,icomp+n) = rho*(m_most.surf_temp + motheta*rho/eta);
                } else {
                    dest_arr(i,j,k,icomp+n) = m_most.surf_temp + motheta/eta;
                }
            });

    } else if (idx == Vars::xvel || idx == Vars::xmom) { //for velx

        amrex::Box xb2d = surroundingNodes(bx,0); // Copy constructor
        xb2d.setBig(2,zlo-1);

        ParallelFor(xb2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                Real velx, vely, rho, eta;
                int jy, ie, je;

                int iylo = i <= lbound(vely_arr).x ? lbound(vely_arr).x : i-1;
                int iyhi = i >  ubound(vely_arr).x ? ubound(vely_arr).x : i;

                jy = j < lbound(vely_arr).y ? lbound(vely_arr).y : j;
                jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                ie = ie > ubound(eta_arr).x-1 ? ubound(eta_arr).x-1 : ie;
                je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;

                velx  = velx_arr(i,j,zlo);
                vely  = 0.25*( vely_arr(iyhi,jy,zlo)+vely_arr(iyhi,jy+1,zlo)
                              +vely_arr(iylo,jy,zlo)+vely_arr(iylo,jy+1,zlo));
                rho   = 0.5*(cons_arr(ie-1,je,zlo,Rho_comp)+
                             cons_arr(ie  ,je,zlo,Rho_comp));
                eta   = 0.5*( eta_arr(ie-1,je,zlo,EddyDiff::Mom_h)+
                              eta_arr(ie  ,je,zlo,EddyDiff::Mom_h));

                Real vmag  = sqrt(velx*velx+vely*vely);
                Real vgx   = ((velx-m_most.vel_mean[0])*m_most.vmag_mean + vmag*m_most.vel_mean[0])/
                              (m_most.vmag_mean*m_most.vmag_mean) * m_most.utau*m_most.utau;

                if (!var_is_derived) {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - vgx*rho/eta;
                } else {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - vgx/eta;
                }
        });

    } else if (idx == Vars::yvel || idx == Vars::ymom) { //for vely

        amrex::Box yb2d = surroundingNodes(bx,1); // Copy constructor
        yb2d.setBig(2,zlo-1);

        ParallelFor(yb2d, [=,m_most=m_most] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                Real velx, vely, rho, eta;
                int ix, ie, je;

                ix = i < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                ix = ix > ubound(velx_arr).x ? ubound(velx_arr).x : ix;

                int jxlo = j <= lbound(velx_arr).y ? lbound(velx_arr).y : j-1;
                int jxhi = j >  ubound(velx_arr).y ? ubound(velx_arr).y : j;

                ie = i < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                je = j < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
                je = je > ubound(eta_arr).y-1 ? ubound(eta_arr).y-1 : je;

                velx  = 0.25*( velx_arr(ix,jxhi,zlo)+velx_arr(ix+1,jxhi,zlo)
                              +velx_arr(ix,jxlo,zlo)+velx_arr(ix+1,jxlo,zlo));
                vely  = vely_arr(i,j,zlo);
                rho   = 0.5*(cons_arr(ie,je-1,zlo,Rho_comp)+
                             cons_arr(ie,je  ,zlo,Rho_comp));
                eta   = 0.5*(eta_arr(ie,je-1,zlo,EddyDiff::Mom_h)+
                             eta_arr(ie,je  ,zlo,EddyDiff::Mom_h));
                Real vmag  = sqrt(velx*velx+vely*vely);
                Real vgy   = ((vely-m_most.vel_mean[1])*m_most.vmag_mean + vmag*m_most.vel_mean[1]) /
                             (m_most.vmag_mean*m_most.vmag_mean)*m_most.utau*m_most.utau;

                if (!var_is_derived) {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - vgy*rho/eta;
                } else {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - vgy/eta;
                }
        });
    }
}
