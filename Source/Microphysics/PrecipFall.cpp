/*!
 * positively definite monotonic advection with non-oscillatory option
 * and gravitational sedimentation
 * NOTES: this code is modified from SAMXX (C++ version of SAM code)
 */

#include "Microphysics.H"

using namespace amrex;

void Microphysics::PrecipFall(int hydro_type) {

  Real constexpr eps = 1.e-10;
  bool constexpr nonos = true;

  Real gam3  = erf_gammafff(3.0             );
  Real gamr1 = erf_gammafff(3.0+b_rain      );
  Real gamr2 = erf_gammafff((5.0+b_rain)/2.0);
  Real gamr3 = erf_gammafff(4.0+b_rain      );
  Real gams1 = erf_gammafff(3.0+b_snow      );
  Real gams2 = erf_gammafff((5.0+b_snow)/2.0);
  Real gams3 = erf_gammafff(4.0+b_snow      );
  Real gamg1 = erf_gammafff(3.0+b_grau      );
  Real gamg2 = erf_gammafff((5.0+b_grau)/2.0);
  Real gamg3 = erf_gammafff(4.0+b_grau      );

  Real vrain = a_rain*gamr3/6.0/pow((PI*rhor*nzeror),crain);
  Real vsnow = a_snow*gams3/6.0/pow((PI*rhos*nzeros),csnow);
  Real vgrau = a_grau*gamg3/6.0/pow((PI*rhog*nzerog),cgrau);

  Real dt_advance = dt;
  int nz = nlev;

  auto qp    = mic_fab_vars[MicVar::qp];
  auto omega = mic_fab_vars[MicVar::omega];
  auto tabs  = mic_fab_vars[MicVar::tabs];
  auto theta = mic_fab_vars[MicVar::theta];

  auto ba    = tabs->boxArray();
  auto dm    = tabs->DistributionMap();
  auto ngrow = tabs->nGrowVect();

  MultiFab mx;
  MultiFab mn;
  MultiFab lfac;
  MultiFab www;
  MultiFab fz;
  MultiFab wp;
  MultiFab tmp_qp;

  mx.define(ba,dm, 1, ngrow);
  mn.define(ba,dm, 1, ngrow);
  lfac.define(ba, dm, 1, ngrow);
  www.define(ba, dm, 1, ngrow);
  fz.define(ba, dm, 1, ngrow);
  wp.define(ba, dm, 1, ngrow);
  tmp_qp.define(ba, dm, 1, ngrow);

  TableData<Real, 1> irho;
  TableData<Real, 1> iwmax;
  TableData<Real, 1> rhofac;

  irho.resize({zlo},{zhi});
  iwmax.resize({zlo},{zhi});
  rhofac.resize({zlo},{zhi});

  auto irho_t    = irho.table();
  auto iwmax_t   = iwmax.table();
  auto rhofac_t  = rhofac.table();
  auto rho1d_t   = rho1d.table();

  auto dz = m_geom.CellSize(2);

  ParallelFor(nz, [=] AMREX_GPU_DEVICE (int k) noexcept {
    rhofac_t(k)  = std::sqrt(1.29/rho1d_t(k));
    irho_t(k)    = 1.0/rho1d_t(k);
    Real wmax    = dz/dt_advance;   // Velocity equivalent to a cfl of 1.0.
    iwmax_t(k)   = 1.0/wmax;
  });

  //  Add sedimentation of precipitation field to the vert. vel.
  MultiFab prec_cfl_fab;
  prec_cfl_fab.define(tabs->boxArray(),tabs->DistributionMap(), 1, tabs->nGrowVect());

  for ( MFIter mfi(lfac, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto lfac_array     = lfac.array(mfi);
     auto omega_array    = omega->array(mfi);
     auto qp_array       = qp->array(mfi);
     auto tabs_array     = tabs->array(mfi);
     auto wp_array       = wp.array(mfi);
     auto www_array      = www.array(mfi);
     auto fz_array       = fz.array(mfi);
     auto prec_cfl_array = prec_cfl_fab.array(mfi);

     const auto& box3d = mfi.tilebox();

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
       Real tmp = term_vel_qp(i,j,k,qp_array(i,j,k),
                  vrain, vsnow, vgrau, crain, csnow, cgrau, rho1d_t(k),
                  tabs_array(i,j,k), a_pr, a_gr);
       wp_array(i,j,k)=rhofac_t(k)*tmp;
       tmp = wp_array(i,j,k)*iwmax_t(k);
       prec_cfl_array(i,j,k) = tmp;
       wp_array(i,j,k) = -wp_array(i,j,k)*rho1d_t(k)*dt_advance/dz;
       if (k == 0) {
         fz_array(i,j,nz-1)   = 0.0;
         www_array(i,j,nz-1)  = 0.0;
         lfac_array(i,j,nz-1) = 0.0;
       }
    });
  }

  auto const& cfl_arrays = prec_cfl_fab.const_arrays();
  Real prec_cfl = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
          prec_cfl_fab, IntVect(0),
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        ->GpuTuple<Real>
        {
          return { cfl_arrays[box_no](i,j,k) };
        });

  // If maximum CFL due to precipitation velocity is greater than 0.9,
  // take more than one advection step to maintain stability.
  int nprec;
  if (prec_cfl > 0.9) {
    nprec = std::ceil(prec_cfl/0.9);
    for (MFIter mfi(wp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      auto wp_array = wp.array(mfi);
      const auto& box3d = mfi.tilebox();

      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int k, int j, int i) {
        // wp already includes factor of dt, so reduce it by a
        // factor equal to the number of precipitation steps.
        wp_array(i,j,k) = wp_array(i,j,k)/Real(nprec);
      });
    }
  } else {
    nprec = 1;
  }

std::cout << "precipfall: nprec= " << nprec << std::endl;

#ifdef ERF_FIXED_SUBCYCLE
    nprec = 4;
#endif

  for(int iprec = 1; iprec<=nprec; iprec++) {
    for ( MFIter mfi(tmp_qp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
       auto qp_array     = qp->array(mfi);
       auto tabs_array   = tabs->array(mfi);
       auto theta_array  = theta->array(mfi);
       auto tmp_qp_array = tmp_qp.array(mfi);
       auto mx_array     = mx.array(mfi);
       auto mn_array     = mn.array(mfi);
       auto fz_array     = fz.array(mfi);
       auto wp_array     = wp.array(mfi);
       auto lfac_array   = lfac.array(mfi);
       auto www_array    = www.array(mfi);

       const auto& box3d = mfi.tilebox();

       ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
         tmp_qp_array(i,j,k) = qp_array(i,j,k); // Temporary array for qp in this column
       });

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

      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        int kc = min(k+1, nz-1);
        tmp_qp_array(i,j,k) = tmp_qp_array(i,j,k)-(fz_array(i,j,kc)-fz_array(i,j,k))*irho_t(k); //Update temporary qp
      });

      ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        // Also, compute anti-diffusive correction to previous
        // (upwind) approximation to the flux
        int kb=max(0,k-1);
        // The precipitation velocity is a cell-centered quantity,
        // since it is computed from the cell-centered
        // precipitation mass fraction.  Therefore, a reformulated
        // anti-diffusive flux is used here which accounts for
        // this and results in reduced numerical diffusion.
        www_array(i,j,k) = 0.5*(1.0+wp_array(i,j,k)*irho_t(k))*(tmp_qp_array(i,j,kb)*wp_array(i,j,kb) -
                           tmp_qp_array(i,j,k)*wp_array(i,j,k)); // works for wp(k)<0
      });

      if (nonos) {
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
          int kc=min(nz-1,k+1);
          int kb=max(0,k-1);
          mx_array(i,j,k) = max(tmp_qp_array(i,j,kb),max(tmp_qp_array(i,j,kc), max(tmp_qp_array(i,j,k), mx_array(i,j,k))));
          mn_array(i,j,k) = min(tmp_qp_array(i,j,kb),min(tmp_qp_array(i,j,kc), min(tmp_qp_array(i,j,k), mn_array(i,j,k))));
          kc = min(nz-1,k+1);
          mx_array(i,j,k) = rho1d_t(k)*(mx_array(i,j,k)-tmp_qp_array(i,j,k))/(pn(www_array(i,j,kc)) +
                                                                              pp(www_array(i,j,k))+eps);
          mn_array(i,j,k) = rho1d_t(k)*(tmp_qp_array(i,j,k)-mn_array(i,j,k))/(pp(www_array(i,j,kc)) +
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
      ParallelFor(box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
        int kc=min(k+1, nz-1);
        // Update precipitation mass fraction.
        // Note that fz is the total flux, including both the
        // upwind flux and the anti-diffusive correction.
        Real flagstat = 1.0;
//        qp_array(i,j,k) = qp_array(i,j,k) - (fz_array(i,j,kc) - fz_array(i,j,k))*irho_t(k);
        Real tmp  = -(fz_array(i,j,kc)-fz_array(i,j,k))*irho_t(k)*flagstat;  // For qp budget
//
// NOTE qpfall, tlat, and precflux,...are output diagnostic variables, not sure whether we need to calculate these variables here?
// Please correct me!!! by xyuan@anl.gov
//
//      amrex::Gpu::Atomic::Add(&qpfall_t(k), tmp);
        Real lat_heat = -(lfac_array(i,j,kc)*fz_array(i,j,kc)-lfac_array(i,j,k)*fz_array(i,j,k))*irho_t(k);
        amrex::Gpu::Atomic::Add(&theta_array(i,j,k), -lat_heat);
//        amrex::Gpu::Atomic::Add(&tlat_t(k), -lat_heat);
        tmp = fz_array(i,j,k)*flagstat;
//      amrex::Gpu::Atomic::Add(&precflux_t(k), -tmp);
//      yakl::atomicAdd(precflux(k,icrm),-tmp);
        if (k == 0) {
//        precsfc_t(i,i)  = precsfc_t(i,j)  - fz_array(i,j,0)*flagstat;
//        precssfc_t(i,j) = precssfc_t(i,j) - fz_array(i,j,0)*(1.0-omega_array(i,j,0))*flagstat;
//        prec_xy_t(i,j)  = prec_xy_t(i,j)  - fz_array(i,j,0)*flagstat;
        }
      });

      if (iprec < nprec) {
        // Re-compute precipitation velocity using new value of qp.
        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
          Real tmp = term_vel_qp(i,j,k,qp_array(i,j,k),
                     vrain, vsnow, vgrau, crain, csnow, cgrau, rho1d_t(k),
                     tabs_array(i,j,k), a_pr, a_gr);
          wp_array(i,j,k) = rhofac_t(k)*tmp;
          // Decrease precipitation velocity by factor of nprec
          wp_array(i,j,k) = -wp_array(i,j,k)*rho1d_t(k)*dt_advance/dz/nprec;
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

  for ( MFIter mfi(*(mic_fab_vars[MicVar::omega]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {
     auto omega_array = mic_fab_vars[MicVar::omega]->array(mfi);
     auto tabs_array  = mic_fab_vars[MicVar::tabs]->array(mfi);

     const auto& box3d = mfi.tilebox();

     ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
       omega_array(i,j,k) = std::max(0.0,std::min(1.0,(tabs_array(i,j,k)-tprmin)*a_pr));
     });
  }
  PrecipFall(2);
}
