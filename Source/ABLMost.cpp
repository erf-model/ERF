#include <ABLMost.H>
#include <PlaneAverage.H>
#include <VelPlaneAverage.H>

using namespace amrex;

void ABLMost::update_fluxes(int lev,
                            amrex::MultiFab& Theta_new, amrex::MultiFab& U_new,
                            amrex::MultiFab& V_new, amrex::MultiFab& W_new,
                            int max_iters)
{
    PlaneAverage Thave(&Theta_new, m_geom[lev], 2, true);
    PlaneAverage vxave(&U_new, m_geom[lev], 2, true);
    PlaneAverage vyave(&V_new, m_geom[lev], 2, true);
    PlaneAverage vzave(&W_new, m_geom[lev], 2, true);
    VelPlaneAverage vmagave({&U_new,&V_new,&W_new}, m_geom[lev], 2, true);

    Thave.compute_averages(ZDir(), Thave.field());
    vxave.compute_averages(ZDir(), vxave.field());
    vyave.compute_averages(ZDir(), vyave.field());
    vzave.compute_averages(ZDir(), vzave.field());
    vmagave.compute_hvelmag_averages(ZDir(), 0, 1, vmagave.field());

    vel_mean[0] = vxave.line_average_interpolated(zref, 0);
    vel_mean[1] = vyave.line_average_interpolated(zref, 0);
    vel_mean[2] = vzave.line_average_interpolated(zref, 0);
    vmag_mean   = vmagave.line_hvelmag_average_interpolated(zref);
    theta_mean  = Thave.line_average_interpolated(zref, 0);

    constexpr amrex::Real eps = 1.0e-16;
    amrex::Real zeta = 0.0;
    amrex::Real utau_iter = 0.0;

    // Initialize variables
    amrex::Real psi_m = 0.0;
    amrex::Real psi_h = 0.0;
    utau = kappa * vmag_mean / (std::log(zref / z0_const));

    int iter = 0;
    do {
            utau_iter = utau;
            switch (alg_type) {
            case HEAT_FLUX:
                surf_temp = surf_temp_flux * (std::log(zref / z0_const) - psi_h) /
                                (utau * kappa) + theta_mean;
                break;

            case SURFACE_TEMPERATURE:
                surf_temp_flux = -(theta_mean - surf_temp) * utau * kappa /
                                (std::log(zref / z0_const) - psi_h);
                break;
            }

            if (std::abs(surf_temp_flux) > eps) {
                // Stable and unstable ABL conditions
                obukhov_len = -utau * utau * utau * theta_mean /
                            (kappa * gravity * surf_temp_flux);
                zeta = zref / obukhov_len;
            } else {
                // Neutral conditions
                obukhov_len = std::numeric_limits<amrex::Real>::max();
                zeta = 0.0;
            }
            psi_m = calc_psi_m(zeta);
            psi_h = calc_psi_h(zeta);
            utau = kappa * vmag_mean / (std::log(zref / z0_const) - psi_m);
            ++iter;
    } while ((std::abs(utau_iter - utau) > 1e-5) && iter <= max_iters);

    if (iter >= max_iters) {
        amrex::Print()
            << "MOData::update_fluxes: Convergence criteria not met after "
            << max_iters << " iterations"
            << "\nObuhov length = " << obukhov_len << " zeta = " << zeta
            << "\npsi_m = " << psi_m << " psi_h = " << psi_h
            << "\nutau = " << utau << " Tsurf = " << surf_temp
            << " q = " << surf_temp_flux << std::endl;
    }
}

void
ABLMost::impose_most_bcs(const int lev, const Box& bx,
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

    // Define temporaries so we can access these on GPU
    Real d_kappa = kappa;
    Real d_utau  = utau;
    Real d_sfcT  = surf_temp;
    Real d_thM   = theta_mean;
    Real d_vmM   = vmag_mean;
    Real d_vxM   = vel_mean[0];
    Real d_vyM   = vel_mean[1];
    Real d_dz    = m_geom[lev].CellSize(2);
    ABLMostData d_most = get_most_data();

    const Array4<Real> z0_arr = z_0[lev].array();

    if (idx == Vars::cons) {
        amrex::Box b2d = bx; // Copy constructor
        b2d.setBig(2,zlo-1);
        int n = Cons::RhoTheta;
        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                int k0 = 0;
                Real d_phi_h = d_most.phi_h(z0_arr(i,j,k0));
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
                eta   = eta_arr(ie,je,zlo,EddyDiff::Theta_v);

                Real vmag    = sqrt(velx*velx+vely*vely);
                Real num1    = (theta-d_thM)*d_vmM;
                Real num2    = (d_thM-d_sfcT)*vmag;
                Real moflux  = (num1+num2)*d_utau*d_kappa/(d_phi_h*d_vmM);
                Real deltaz  = d_dz * (zlo - k);

                if (!var_is_derived) {
                    dest_arr(i,j,k,icomp+n) = rho*(theta - moflux*rho/eta*deltaz);
                } else {
                    dest_arr(i,j,k,icomp+n) = theta - moflux/eta*deltaz;
                }
            });

    } else if (idx == Vars::xvel || idx == Vars::xmom) { //for velx

        amrex::Box xb2d = surroundingNodes(bx,0); // Copy constructor
        xb2d.setBig(2,zlo-1);

        ParallelFor(xb2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
                eta   = 0.5*( eta_arr(ie-1,je,zlo,EddyDiff::Mom_v)+
                              eta_arr(ie  ,je,zlo,EddyDiff::Mom_v));

                Real vmag  = sqrt(velx*velx+vely*vely);
                Real stressx = ((velx-d_vxM)*d_vmM + vmag*d_vxM)/
                               (d_vmM*d_vmM) * d_utau*d_utau;
                Real deltaz  = d_dz * (zlo - k);

                if (!var_is_derived) {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressx*rho*rho/eta*deltaz;
                } else {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressx*rho/eta*deltaz;
                }
        });

    } else if (idx == Vars::yvel || idx == Vars::ymom) { //for vely

        amrex::Box yb2d = surroundingNodes(bx,1); // Copy constructor
        yb2d.setBig(2,zlo-1);

        ParallelFor(yb2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
                eta   = 0.5*(eta_arr(ie,je-1,zlo,EddyDiff::Mom_v)+
                             eta_arr(ie,je  ,zlo,EddyDiff::Mom_v));
                Real vmag  = sqrt(velx*velx+vely*vely);
                Real stressy = ((vely-d_vyM)*d_vmM + vmag*d_vyM) /
                               (d_vmM*d_vmM)*d_utau*d_utau;
                Real deltaz  = d_dz * (zlo - k);

                if (!var_is_derived) {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressy*rho*rho/eta*deltaz;
                } else {
                    dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressy*rho/eta*deltaz;
                }
        });
    }
}
