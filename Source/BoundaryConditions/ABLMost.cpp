#include <ABLMost.H>
#include <MOSTAverage.H>

using namespace amrex;

/**
 * Function to update the fluxs (u^star and t^star) for Monin Obukhov similarity theory.
 *
 * @param[in] lev Current level
 * @param[in] max_iters maximum iterations to use
 */
void ABLMost::update_fluxes(int lev, int max_iters)
{
    // Compute plane averages for all vars
    m_ma.compute_averages(lev);

    // Pointers to the computed averages
    const auto *const tm_ptr  = m_ma.get_average(lev,2);
    const auto *const umm_ptr = m_ma.get_average(lev,3);

    // GPU device captures
    amrex::Real d_kappa = kappa;
    amrex::Real d_zref  = m_ma.get_zref();
    amrex::Real d_gravity = gravity;
    amrex::Real d_surf_temp_flux = surf_temp_flux;
    ABLMostData d_most = get_most_data();

    // Tolerances
    constexpr amrex::Real eps = std::numeric_limits<Real>::epsilon();
    constexpr amrex::Real tol = 1.0e-5;

    // Ghost cells for CC var
    amrex::IntVect ng = u_star[lev]->nGrowVect(); ng[2]=0;

    // Initialize to the adiabatic q=0 case
    if (mom_type == MomCalcType::CHARNOCK) {
        for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(ng);

            auto t_surf_arr = t_surf[lev]->array(mfi);
            auto t_star_arr = t_star[lev]->array(mfi);
            auto u_star_arr = u_star[lev]->array(mfi);

            const auto tm_arr  = tm_ptr->array(mfi);
            const auto umm_arr = umm_ptr->array(mfi);
            const auto z0_arr  = z_0[lev].array();

            auto olen_arr = olen[lev]->array(mfi);

            amrex::Real d_a = d_most.Cnk_a;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int iter = 0;
                u_star_arr(i,j,k) = umm_arr(i,j,k);
                amrex::Real ustar = 0.0;
                amrex::Real z0    = 0.0;
                do {
                    ustar = u_star_arr(i,j,k);
                    z0    = (d_a / d_gravity) * ustar * ustar;
                    u_star_arr(i,j,k) = d_kappa * umm_arr(i,j,k) / std::log(d_zref / z0);
                    ++iter;
                } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
                olen_arr(i,j,k) = std::numeric_limits<Real>::max();
            });
        }
    } else if (mom_type == MomCalcType::MODIFIED_CHARNOCK) {
        for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(ng);

            auto t_surf_arr = t_surf[lev]->array(mfi);
            auto t_star_arr = t_star[lev]->array(mfi);
            auto u_star_arr = u_star[lev]->array(mfi);

            const auto tm_arr  = tm_ptr->array(mfi);
            const auto umm_arr = umm_ptr->array(mfi);
            const auto z0_arr  = z_0[lev].array();

            auto olen_arr = olen[lev]->array(mfi);

            amrex::Real d_b = d_most.Cnk_b1 * std::log(d_most.Cnk_b2 / d_most.Cnk_d);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int iter = 0;
                u_star_arr(i,j,k) = umm_arr(i,j,k);
                amrex::Real ustar = 0.0;
                amrex::Real z0    = 0.0;
                do {
                    ustar = u_star_arr(i,j,k);
                    z0    = std::exp( (2.7*ustar - 1.8/d_b) / (ustar + 0.17/d_b) );
                    u_star_arr(i,j,k) = d_kappa * umm_arr(i,j,k) / std::log(d_zref / z0);
                    ++iter;
                } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
                olen_arr(i,j,k) = std::numeric_limits<Real>::max();
            });
        }
    } else {
        for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(ng);

            auto t_star_arr = t_star[lev]->array(mfi);
            auto u_star_arr = u_star[lev]->array(mfi);
            auto olen_arr   = olen[lev]->array(mfi);

            const auto umm_arr = umm_ptr->array(mfi);
            const auto z0_arr  = z_0[lev].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                t_star_arr(i,j,k) = 0.0;
                olen_arr(i,j,k)   = std::numeric_limits<Real>::max();
                u_star_arr(i,j,k) = d_kappa * umm_arr(i,j,k) / std::log(d_zref / z0_arr(i,j,k));
            });
        }
    }

    // Specified finite heat flux
    if ( (theta_type == ThetaCalcType::HEAT_FLUX) && (std::abs(surf_temp_flux) > eps) ) {
        for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(ng);

            auto t_surf_arr = t_surf[lev]->array(mfi);
            auto t_star_arr = t_star[lev]->array(mfi);
            auto u_star_arr = u_star[lev]->array(mfi);

            const auto tm_arr  = tm_ptr->array(mfi);
            const auto umm_arr = umm_ptr->array(mfi);
            const auto z0_arr  = z_0[lev].array();

            auto olen_arr = olen[lev]->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                int iter = 0;
                amrex::Real ustar = 0.0;
                amrex::Real zeta  = 0.0;
                amrex::Real psi_m = 0.0;
                amrex::Real psi_h = 0.0;
                amrex::Real Olen  = 0.0;
                do {
                    ustar = u_star_arr(i,j,k);
                    Olen = -ustar * ustar * ustar * tm_arr(i,j,k) /
                           (d_kappa * d_gravity * d_surf_temp_flux);
                    zeta  = d_zref / Olen;
                    psi_m = d_most.calc_psi_m(zeta);
                    psi_h = d_most.calc_psi_h(zeta);
                    u_star_arr(i,j,k) = d_kappa * umm_arr(i,j,k) / (std::log(d_zref / z0_arr(i,j,k)) - psi_m);
                    ++iter;
                } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

                t_surf_arr(i,j,k) = d_surf_temp_flux * (std::log(d_zref / z0_arr(i,j,k)) - psi_h) /
                                    (u_star_arr(i,j,k) * d_kappa) + tm_arr(i,j,k);
                t_star_arr(i,j,k) = -d_surf_temp_flux / u_star_arr(i,j,k);
                olen_arr(i,j,k)   = Olen;
            });
        }
    // Specified surface temperature
    } else if ( theta_type == ThetaCalcType::SURFACE_TEMPERATURE ) {
        for (MFIter mfi(*u_star[lev]); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(ng);

            auto t_surf_arr = t_surf[lev]->array(mfi);
            auto t_star_arr = t_star[lev]->array(mfi);
            auto u_star_arr = u_star[lev]->array(mfi);

            const auto tm_arr  = tm_ptr->array(mfi);
            const auto umm_arr = umm_ptr->array(mfi);
            const auto z0_arr  = z_0[lev].array();

            auto olen_arr = olen[lev]->array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                // Nothing to do unless the flux != 0
                if (std::abs(t_surf_arr(i,j,k)-tm_arr(i,j,k)) > eps) {
                    int iter = 0;
                    amrex::Real ustar = 0.0;
                    amrex::Real tflux = 0.0;
                    amrex::Real zeta  = 0.0;
                    amrex::Real psi_m = 0.0;
                    amrex::Real psi_h = 0.0;
                    amrex::Real Olen  = 0.0;
                    do {
                        ustar = u_star_arr(i,j,k);
                        tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * d_kappa /
                                 (std::log(d_zref / z0_arr(i,j,k)) - psi_h);
                        Olen = -ustar * ustar * ustar * tm_arr(i,j,k) /
                               (d_kappa * d_gravity * tflux);
                        zeta  = d_zref / Olen;
                        psi_m = d_most.calc_psi_m(zeta);
                        psi_h = d_most.calc_psi_h(zeta);
                        u_star_arr(i,j,k) = d_kappa * umm_arr(i,j,k) / (std::log(d_zref / z0_arr(i,j,k)) - psi_m);
                        ++iter;
                    } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

                    t_star_arr(i,j,k) = d_kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                                        (std::log(d_zref / z0_arr(i,j,k)) - psi_h);
                    olen_arr(i,j,k)   = Olen;
                }
            });
        }
    }
}

/**
 * Function to impose Monin Obukhov similarity theory fluxes by populating ghost cells.
 *
 * @param[in] lev Current level
 * @param[in,out] mfs Multifabs to populate
 * @param[in] eddyDiffs Diffusion coefficients from turbulence model
 */
void
ABLMost::impose_most_bcs(const int lev,
                         const Vector<MultiFab*>& mfs,
                         MultiFab* eddyDiffs)
{

    const int icomp = 0;
    for (MFIter mfi(*mfs[0]); mfi.isValid(); ++mfi)
    {
        // Get field arrays
        const auto cons_arr = mfs[Vars::cons]->array(mfi);
        const auto velx_arr = mfs[Vars::xvel]->array(mfi);
        const auto vely_arr = mfs[Vars::yvel]->array(mfi);
        const auto  eta_arr = eddyDiffs->array(mfi);

        // Get average arrays
        const auto *const u_mean     = m_ma.get_average(lev,0);
        const auto *const v_mean     = m_ma.get_average(lev,1);
        const auto *const t_mean     = m_ma.get_average(lev,2);
        const auto *const u_mag_mean = m_ma.get_average(lev,3);

        const auto um_arr  = u_mean->array(mfi);
        const auto vm_arr  = v_mean->array(mfi);
        const auto tm_arr  = t_mean->array(mfi);
        const auto umm_arr = u_mag_mean->array(mfi);

        // Get derived arrays
        const auto u_star_arr = u_star[lev]->array(mfi);
        const auto t_star_arr = t_star[lev]->array(mfi);
        const auto t_surf_arr = t_surf[lev]->array(mfi);

        // Define temporaries so we can access these on GPU
        Real d_dz = m_geom[lev].CellSize(2);

        for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
        {
            const Box& bx = (*mfs[var_idx])[mfi].box();
            auto dest_arr = (*mfs[var_idx])[mfi].array();
            int zlo = 0;

            if (var_idx == Vars::cons) {
                amrex::Box b2d = bx; // Copy constructor
                b2d.setBig(2,zlo-1);
                int n = Cons::RhoTheta;

                ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real velx, vely, rho, theta, eta;
                    int ix, jx, iy, jy, ie, je, ic, jc;

                    ix = i < lbound(velx_arr).x    ? lbound(velx_arr).x   : i;
                    jx = j < lbound(velx_arr).y    ? lbound(velx_arr).y   : j;
                    ix = ix > ubound(velx_arr).x-1 ? ubound(velx_arr).x-1 : ix;
                    jx = jx > ubound(velx_arr).y   ? ubound(velx_arr).y   : jx;

                    iy = i  < lbound(vely_arr).x   ? lbound(vely_arr).x   : i;
                    jy = j  < lbound(vely_arr).y   ? lbound(vely_arr).y   : j;
                    iy = iy > ubound(vely_arr).x   ? ubound(vely_arr).x   : iy;
                    jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

                    ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
                    je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
                    ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
                    je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;

                    ic = i  < lbound(cons_arr).x ? lbound(cons_arr).x : i;
                    jc = j  < lbound(cons_arr).y ? lbound(cons_arr).y : j;
                    ic = ic > ubound(cons_arr).x ? ubound(cons_arr).x : ic;
                    jc = jc > ubound(cons_arr).y ? ubound(cons_arr).y : jc;

                    velx  = 0.5*(velx_arr(ix,jx,zlo)+velx_arr(ix+1,jx  ,zlo));
                    vely  = 0.5*(vely_arr(iy,jy,zlo)+vely_arr(iy  ,jy+1,zlo));
                    rho   = cons_arr(ic,jc,zlo,Rho_comp);
                    theta = cons_arr(ic,jc,zlo,RhoTheta_comp) / rho;
                    eta   = eta_arr(ie,je,zlo,EddyDiff::Theta_v); // == rho * alpha [kg/m^3 * m^2/s]

                    Real d_thM  =  tm_arr(ic,jc,zlo);
                    Real d_vmM  = umm_arr(ic,jc,zlo);
                    Real d_utau = u_star_arr(ic,jc,zlo);
                    Real d_ttau = t_star_arr(ic,jc,zlo);
                    Real d_sfcT = t_surf_arr(ic,jc,zlo);

                    Real vmag    = sqrt(velx*velx+vely*vely);
                    Real num1    = (theta-d_thM)*d_vmM;
                    Real num2    = (d_thM-d_sfcT)*vmag;
                    Real moflux  = d_ttau*d_utau*(num1+num2)/((d_thM-d_sfcT)*d_vmM);
                    Real deltaz  = d_dz * (zlo - k);

                    dest_arr(i,j,k,icomp+n) = rho*(theta - moflux*rho/eta*deltaz);
                });

            } else if (var_idx == Vars::xvel || var_idx == Vars::xmom) { //for velx

                amrex::Box xb2d = surroundingNodes(bx,0); // Copy constructor
                xb2d.setBig(2,zlo-1);

                ParallelFor(xb2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real velx, vely, rho, eta;
                    int jy, ie, je, ic, jc;

                    int iylo = i <= lbound(vely_arr).x ? lbound(vely_arr).x : i-1;
                    int iyhi = i >  ubound(vely_arr).x ? ubound(vely_arr).x : i;

                    jy = j  < lbound(vely_arr).y   ? lbound(vely_arr).y   : j;
                    jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

                    ie = i  < lbound(eta_arr).x+1 ? lbound(eta_arr).x+1 : i;
                    je = j  < lbound(eta_arr).y   ? lbound(eta_arr).y   : j;
                    ie = ie > ubound(eta_arr).x   ? ubound(eta_arr).x   : ie;
                    je = je > ubound(eta_arr).y   ? ubound(eta_arr).y   : je;

                    ic = i  < lbound(cons_arr).x+1 ? lbound(cons_arr).x+1 : i;
                    jc = j  < lbound(cons_arr).y   ? lbound(cons_arr).y   : j;
                    ic = ic > ubound(cons_arr).x   ? ubound(cons_arr).x   : ic;
                    jc = jc > ubound(cons_arr).y   ? ubound(cons_arr).y   : jc;

                    velx  = velx_arr(i,j,zlo);
                    vely  = 0.25*( vely_arr(iyhi,jy,zlo)+vely_arr(iyhi,jy+1,zlo)
                                  +vely_arr(iylo,jy,zlo)+vely_arr(iylo,jy+1,zlo));
                    rho   = 0.5 *( cons_arr(ic-1,jc,zlo,Rho_comp)+
                                   cons_arr(ic  ,jc,zlo,Rho_comp));
                    eta   = 0.5 *( eta_arr(ie-1,je,zlo,EddyDiff::Mom_v)+
                                   eta_arr(ie  ,je,zlo,EddyDiff::Mom_v));

                    Real d_vxM  = um_arr(i,j,zlo);
                    Real d_vmM  = 0.5 * ( umm_arr(ic-1,jc,zlo) + umm_arr(ic,jc,zlo) );
                    Real d_utau = 0.5 * ( u_star_arr(ic-1,jc,zlo) + u_star_arr(ic,jc,zlo) );

                    Real vmag    = sqrt(velx*velx+vely*vely);
                    Real stressx = ( (velx-d_vxM)*d_vmM + vmag*d_vxM )/
                                   (d_vmM*d_vmM) * d_utau*d_utau;
                    Real deltaz  = d_dz * (zlo - k);

                    if (var_idx == Vars::xmom) {
                        dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressx*rho*rho/eta*deltaz;
                    } else {
                        AMREX_ALWAYS_ASSERT(var_idx == Vars::xvel);
                        dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressx*rho/eta*deltaz;
                    }
                });

            } else if (var_idx == Vars::yvel || var_idx == Vars::ymom) { //for vely

                amrex::Box yb2d = surroundingNodes(bx,1); // Copy constructor
                yb2d.setBig(2,zlo-1);

                ParallelFor(yb2d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real velx, vely, rho, eta;
                    int ix, ie, je, ic, jc;

                    ix = i  < lbound(velx_arr).x ? lbound(velx_arr).x : i;
                    ix = ix > ubound(velx_arr).x ? ubound(velx_arr).x : ix;

                    int jxlo = j <= lbound(velx_arr).y ? lbound(velx_arr).y : j-1;
                    int jxhi = j >  ubound(velx_arr).y ? ubound(velx_arr).y : j;

                    ie = i  < lbound(eta_arr).x   ? lbound(eta_arr).x   : i;
                    je = j  < lbound(eta_arr).y+1 ? lbound(eta_arr).y+1 : j;
                    ie = ie > ubound(eta_arr).x   ? ubound(eta_arr).x   : ie;
                    je = je > ubound(eta_arr).y   ? ubound(eta_arr).y   : je;

                    ic = i  < lbound(cons_arr).x   ? lbound(cons_arr).x   : i;
                    jc = j  < lbound(cons_arr).y+1 ? lbound(cons_arr).y+1 : j;
                    ic = ic > ubound(cons_arr).x   ? ubound(cons_arr).x   : ic;
                    jc = jc > ubound(cons_arr).y   ? ubound(cons_arr).y   : jc;

                    velx  = 0.25*( velx_arr(ix,jxhi,zlo)+velx_arr(ix+1,jxhi,zlo)
                                  +velx_arr(ix,jxlo,zlo)+velx_arr(ix+1,jxlo,zlo));
                    vely  = vely_arr(i,j,zlo);
                    rho   = 0.5*(cons_arr(ic,jc-1,zlo,Rho_comp)+
                                 cons_arr(ic,jc  ,zlo,Rho_comp));
                    eta   = 0.5*(eta_arr(ie,je-1,zlo,EddyDiff::Mom_v)+
                                 eta_arr(ie,je  ,zlo,EddyDiff::Mom_v));

                    Real d_vyM  = vm_arr(i,j,zlo);
                    Real d_vmM  = 0.5 * ( umm_arr(ic,jc-1,zlo) + umm_arr(ic,jc,zlo) );
                    Real d_utau = 0.5 * ( u_star_arr(ic,jc-1,zlo) + u_star_arr(ic,jc,zlo) );

                    Real vmag    = sqrt(velx*velx+vely*vely);
                    Real stressy = ( (vely-d_vyM)*d_vmM + vmag*d_vyM ) /
                                    (d_vmM*d_vmM)*d_utau*d_utau;
                    Real deltaz  = d_dz * (zlo - k);

                    if (var_idx == Vars::ymom) {
                        dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressy*rho*rho/eta*deltaz;
                    } else {
                        AMREX_ALWAYS_ASSERT(var_idx == Vars::yvel);
                        dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressy*rho/eta*deltaz;
                    }
                });
            }
        } // var_idx
    } // mf
}
