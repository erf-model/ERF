/*!
 * positively definite monotonic advection with non-oscillatory option
 * and gravitational sedimentation
 * NOTES: this code is modified from SAMXX (C++ version of SAM code)
 */

#include "Microphysics.H"
#include "ERF_Constants.H"

using namespace amrex;

void Microphysics::PrecipFall(int hydro_type) {

  amrex::Real constexpr eps = 1.e-10;
  bool constexpr nonos = true;

  const auto& box3d = m_geom.Domain();
  const auto& lo = amrex::lbound(box3d);
  const auto& hi = amrex::ubound(box3d);

  const auto nx = hi.x - lo.x + 1;
  const auto ny = hi.y - lo.y + 1;
  const auto nz = hi.z - lo.z + 1;

  auto qt    = mic_fab_vars[MicVar::qt];
  auto omega = mic_fab_vars[MicVar::omega];
  auto tabs  = mic_fab_vars[MicVar::tabs];
  auto theta = mic_fab_vars[MicVar::theta];

  MultiFab mx;
  mx.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  mx.setVal(0.);

  MultiFab mn;
  mn.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  mn.setVal(0.);

  MultiFab lfac;
  lfac.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  lfac.setVal(0.);

  MultiFab www;
  www.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  www.setVal(0.);

  MultiFab fz;
  fz.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  fz.setVal(0.);

  MultiFab wp;
  wp.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  wp.setVal(0.);

  MultiFab tmp_qp;
  tmp_qp.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());
  tmp_qp.setVal(0.);

  TableData<Real, 1> irhoadz;
  TableData<Real, 1> iwmax;
  TableData<Real, 1> rhofac;
  TableData<Real, 1> adz;

  irhoadz.resize({0},{nz});
  iwmax.resize({0},{nz});
  rhofac.resize({0},{nz});
  adz.resize({0},{nz});

  auto irhoadz_t = irhoadz.table();
  auto iwmax_t   = iwmax.table();
  auto rhofac_t  = rhofac.table();
  auto adz_t     = adz.table();
  auto rho1d_t   = rho1d.table();

  auto dz = m_geom.CellSize(2);

  // for (int k=0; k<nzm; k++) {
  //  for (int icrm=0; icrm<ncrms; icrm++) {
  //parallel_for( SimpleBounds<2>(nzm,ncrms) , YAKL_LAMBDA (int k, int icrm) {
  ParallelFor(nz, [=] AMREX_GPU_DEVICE (int k) noexcept {
    rhofac_t(k)  = std::sqrt(1.29/rho1d_t(k));
    adz_t(k)     = 1.0;
    irhoadz_t(k) = 1.0/(rho1d_t(k)*adz_t(k));
    int kb       = std::max(0,k-1);
    Real wmax    = dz*adz_t(kb)/dt;   // Velocity equivalent to a cfl of 1.0.
    iwmax_t(k)   = 1.0/wmax;
  });

  //  Add sedimentation of precipitation field to the vert. vel.
  Real prec_cfl = 0.0;
  TableData<Real, 3> prec_cfl_arr;
  prec_cfl_arr.resize({0,0,0},{nx,ny,nz});
  auto prec_cfl_arr_t = prec_cfl_arr.table();

  // for (int k=0; k<nzm; k++) {
  //   for (int j=0; j<ny; j++) {
  //     for (int i=0; i<nx; i++) {
  //       for (int icrm=0; icrm<ncrms; icrm++) {
  for ( MFIter mfi(lfac, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto lfac_array  = lfac.array(mfi);
     auto omega_array = omega->array(mfi);
     auto qt_array    = qt->array(mfi);
     auto tabs_array  = tabs->array(mfi);
     auto wp_array    = wp.array(mfi);
     auto www_array   = www.array(mfi);
     auto fz_array    = fz.array(mfi);

     ParallelFor(box3d, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
       if (hydro_type == 0) {
         lfac_array(i,j,k) = fac_cond;
       }
       else if (hydro_type == 1) {
         lfac_array(i,j,k) = fac_sub;
       }
       else if (hydro_type == 2) {
         lfac_array(i,j,k) = fac_cond + (1.0-omega_array(i,j,k))*fac_fus;
       }
       else if (hydro_type == 3) {
         lfac_array(i,j,k) = 0.0;
       }
       Real tmp = term_vel_qp(i,j,k,qt_array(i,j,k),
                  vrain, vsnow, vgrau, crain, csnow, cgrau, rho1d_t(k),
                  tabs_array(i,j,k), a_pr, a_gr);
       wp_array(i,j,k)=rhofac_t(k)*tmp;
       tmp = wp_array(i,j,k)*iwmax_t(k);
       prec_cfl_arr_t(i,j,k) = tmp;
       wp_array(i,j,k) = -wp_array(i,j,k)*rho1d_t(k)*dt/dz;
       if (k == 0) {
         fz_array(i,j,nz-1)   = 0.0;
         www_array(i,j,nz-1)  = 0.0;
         lfac_array(i,j,nz-1) = 0.0;
       }
    });
  }

  Real prec_cfl_loc  = 0.0;
  ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
     amrex::Gpu::Atomic::Min(&prec_cfl_arr_t(k,j,i), prec_cfl_loc);
  });
  prec_cfl = std::max(prec_cfl, prec_cfl_loc);

  // If maximum CFL due to precipitation velocity is greater than 0.9,
  // take more than one advection step to maintain stability.
  int nprec;
  if (prec_cfl > 0.9) {
    nprec = std::ceil(prec_cfl/0.9);
    // for (int k=0; k<nzm; k++) {
    //   for (int j=0; j<ny; j++) {
    //     for (int i=0; i<nx; i++) {
    //       for (int icrm=0; icrm<ncrms; icrm++) {
    for (MFIter mfi(wp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      auto wp_array = wp.array(mfi);

      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int k, int j, int i) {
        // wp already includes factor of dt, so reduce it by a
        // factor equal to the number of precipitation steps.
        wp_array(i,j,k) = wp_array(i,j,k)/nprec;
      });
    }
  } else {
    nprec = 1;
  }

  //if(nprec > 1){std::cout << nprec << std::endl;}

//#ifdef MMF_FIXED_SUBCYCLE
    nprec = 4;
//#endif

  for(int iprec = 1; iprec<=nprec; iprec++) {
    //parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
    for ( MFIter mfi(tmp_qp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       auto qt_array     = qt->array(mfi);
       auto tabs_array   = tabs->array(mfi);
       auto theta_array  = theta->array(mfi);
       auto tmp_qp_array = tmp_qp.array(mfi);
       auto mx_array     = mx.array(mfi);
       auto mn_array     = mn.array(mfi);
       auto fz_array     = fz.array(mfi);
       auto wp_array     = wp.array(mfi);
       auto lfac_array   = lfac.array(mfi);
       auto www_array    = www.array(mfi);

       ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
         tmp_qp_array(i,j,k) = qt_array(i,j,k); // Temporary array for qp in this column
       });

      //parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        if (nonos) {
          int kc=min(nz-1,k+1);
          int kb=max(0,k-1);
          mx_array(i,j,k) = max(tmp_qp_array(i,j,kb), max(tmp_qp_array(i,j,kc), tmp_qp_array(i,j,k)));
          mn_array(i,j,k) = min(tmp_qp_array(i,j,kb), min(tmp_qp_array(i,j,kc), tmp_qp_array(i,j,k)));
        }
        // Define upwind precipitation flux
        fz_array(i,j,k) = tmp_qp_array(i,j,k)*wp_array(i,j,k);
      });

      //parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        int kc=k+1;
        if (kc < nz)
          tmp_qp_array(i,j,k) = tmp_qp_array(i,j,k)-(fz_array(i,j,kc)-fz_array(i,j,k))*irhoadz_t(k); //Update temporary qp
      });

      //parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        // Also, compute anti-diffusive correction to previous
        // (upwind) approximation to the flux
        int kb=max(0,k-1);
        // The precipitation velocity is a cell-centered quantity,
        // since it is computed from the cell-centered
        // precipitation mass fraction.  Therefore, a reformulated
        // anti-diffusive flux is used here which accounts for
        // this and results in reduced numerical diffusion.
        www_array(i,j,k) = 0.5*(1.0+wp_array(i,j,k)*irhoadz_t(k))*(tmp_qp_array(i,j,kb)*wp_array(i,j,kb) -
                           tmp_qp_array(i,j,k)*wp_array(i,j,k)); // works for wp(k)<0
      });

      if (nonos) {
        //parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
          int kc=min(nz-1,k+1);
          int kb=max(0,k-1);
          mx_array(i,j,k) = max(tmp_qp_array(i,j,kb),max(tmp_qp_array(i,j,kc), max(tmp_qp_array(i,j,k), mx_array(i,j,k))));
          mn_array(i,j,k) = min(tmp_qp_array(i,j,kb),min(tmp_qp_array(i,j,kc), min(tmp_qp_array(i,j,k), mn_array(i,j,k))));
          kc = min(nz-1,k+1);
          mx_array(i,j,k) = rho1d_t(k)*adz_t(k)*(mx_array(i,j,k)-tmp_qp_array(i,j,k))/(pn(www_array(i,j,kc)) +
                                                                                       pp(www_array(i,j,k))+eps);
          mn_array(i,j,k) = rho1d_t(k)*adz_t(k)*(tmp_qp_array(i,j,k)-mn_array(i,j,k))/(pp(www_array(i,j,kc)) +
                                                                                       pn(www_array(i,j,k))+eps);
        });

        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
          int kb=max(0,k-1);
          // Add limited flux correction to fz(k).
          fz_array(i,j,k) = fz_array(i,j,k) + pp(www_array(i,j,k))*std::min(1.0,std::min(mx_array(i,j,k), mn_array(i,j,kb))) -
                                              pn(www_array(i,j,k))*std::min(1.0,std::min(mx_array(i,j,kb),mn_array(i,j,k))); // Anti-diffusive flux
        });
      }

      // Update precipitation mass fraction and liquid-ice static
      // energy using precipitation fluxes computed in this column.

      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        int kc=min(k+1, nz-1);
        // Update precipitation mass fraction.
        // Note that fz is the total flux, including both the
        // upwind flux and the anti-diffusive correction.
        Real flagstat = 1.0;
        qt_array(i,j,k) = qt_array(i,j,k)-(fz_array(i,j,kc) - fz_array(i,j,k))*irhoadz_t(k);
        Real tmp  = -(fz_array(i,j,kc)-fz_array(i,j,k))*irhoadz_t(k)*flagstat;  // For qp budget
//
// NOTE qpfall, tlat, and precflux,...are output diagnostic variables, not sure whether we need to calculate these variables here?
// Please correct me!!! by xyuan@anl.gov
//
//      amrex::Gpu::Atomic::Add(&qpfall_t(k), tmp);
        Real lat_heat = -(lfac_array(i,j,kc)*fz_array(i,j,kc)-lfac_array(i,j,k)*fz_array(i,j,k))*irhoadz_t(k);
        theta_array(i,j,k) = theta_array(i,j,k)-lat_heat*inv_cp;
//      amrex::Gpu::Atomic::Add(&tlat_t(k), -lat_heat);
//      tmp = fz_t(k,j,i)*flagstat;
//      amrex::Gpu::Atomic::Add(&precflux_t(k), -tmp);

//if(qt_t(k,j,i) < 0.08) {
//std::cout << "precipfall: " << i << "; " << j << "; " << k << "; lat " << lat_heat << "; theta " << t_t(k,j,i) << "; qt " << qt_t(k,j,i)
//          << "; fz " << fz_t(k,j,i) << "; fz-1 " << fz_t(kc,j,i) << "; irhoadz " << irhoadz_t(k) << std::endl;
//}

//      yakl::atomicAdd(precflux(k,icrm),-tmp);
        if (k == 0) {
//        precsfc_t(j,i)  = precsfc_t(j,i)  - fz_t(0,j,i)*flagstat;
//        precssfc_t(j,i) = precssfc_t(j,i) - fz_t(0,j,i)*(1.0-omega_t(0,j,i))*flagstat;
//        prec_xy_t(j,i)  = prec_xy_t(j,i)  - fz_t(0,j,i)*flagstat;
        }
      });

      if (iprec < nprec) {
        // Re-compute precipitation velocity using new value of qp.
        //parallel_for( SimpleBounds<4>(nzm,ny,nx,ncrms) , YAKL_LAMBDA (int k, int j, int i, int icrm) {
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
          Real tmp = term_vel_qp(i,j,k,qt_array(i,j,k),
                     vrain, vsnow, vgrau, crain, csnow, cgrau, rho1d_t(k),
                     tabs_array(i,j,k), a_pr, a_gr);
          wp_array(i,j,k) = rhofac_t(k)*tmp;
          // Decrease precipitation velocity by factor of nprec
          wp_array(i,j,k) = -wp_array(i,j,k)*rho1d_t(k)*dt/dz/nprec;
          // Note: Don't bother checking CFL condition at each
          // substep since it's unlikely that the CFL will
          // increase very much between substeps when using
          // monotonic advection schemes.
          if (k == 0) {
            fz_array(i,j,nz-1)   = 0.0;
            www_array(i,j,nz-1)  = 0.0;
            lfac_array(i,j,nz-1) = 0.0;
          }
        });
      }
    }
  } // iprec loop
}


void Microphysics::MicroPrecipFall() {

  crain = b_rain / 4.0;
  csnow = b_snow / 4.0;
  cgrau = b_grau / 4.0;
  vrain = a_rain * gamr3 / 6.0 / pow((PI * rhor * nzeror), crain);
  vsnow = a_snow * gams3 / 6.0 / pow((PI * rhos * nzeros), csnow);
  vgrau = a_grau * gamg3 / 6.0 / pow((PI * rhog * nzerog), cgrau);

  const auto& box3d = m_geom.Domain();

  for ( MFIter mfi(*(mic_fab_vars[MicVar::omega]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto omega_array = mic_fab_vars[MicVar::omega]->array(mfi);
     auto tabs_array  = mic_fab_vars[MicVar::tabs]->array(mfi);

     ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
       omega_array(i,j,k) = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
     });
  }
  PrecipFall(2);
}
