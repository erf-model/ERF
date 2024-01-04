#include "prob.H"
#include "EOS.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);

  init_base_parms(parms.rho_0, parms.T_0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_temperature(amrex::Real z, amrex::Real z_0, amrex::Real z_trop, amrex::Real z_top,
                           amrex::Real T_0, amrex::Real T_trop, amrex::Real T_top)
{
    if (z <= z_trop) {
      amrex::Real lapse = - (T_trop - T_0) / (z_trop - z_0);
      return T_0 - lapse * (z - z_0);
    } else {
      amrex::Real lapse = - (T_top - T_trop) / (z_top - z_trop);
      return T_trop - lapse * (z - z_trop);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_pressure(amrex::Real z, amrex::Real z_0,
                        amrex::Real z_trop, amrex::Real z_top,
                        amrex::Real T_0, amrex::Real T_trop, amrex::Real T_top)
{
    if (z <= z_trop) {
      amrex::Real lapse = - (T_trop - T_0) / (z_trop - z_0);
      amrex::Real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
      return p_0 * std::pow( T / T_0 , CONST_GRAV/(R_d*lapse) );
    } else {
      // Get pressure at the tropopause
      amrex::Real lapse = - (T_trop - T_0) / (z_trop - z_0);
      amrex::Real p_trop = p_0 * std::pow( T_trop / T_0 , CONST_GRAV/(R_d*lapse) );
      // Get pressure at requested height
      lapse = - (T_top - T_trop) / (z_top - z_trop);
      if (lapse != 0) {
        amrex::Real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
        return p_trop * std::pow( T / T_trop , CONST_GRAV/(R_d*lapse) );
      } else {
        return p_trop * std::exp(-CONST_GRAV*(z-z_trop)/(R_d*T_trop));
      }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_sat_mix(amrex::Real press, amrex::Real T ) {
    return 380./(press) * std::exp(17.27*(T-273.)/(T-36.));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
init_supercell_relhum(amrex::Real z, amrex::Real z_trop)
{
  if (z <= z_trop) {
      return 1. - 0.75 * pow(z / z_trop , 1.25);
  } else {
      return 0.25;
  }
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& state,
    Array4<Real      > const& x_vel,
    Array4<Real      > const& y_vel,
    Array4<Real      > const& z_vel,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    const int khi = geomdata.Domain().bigEnd()[2];

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real& dz        = geomdata.CellSize()[2];
  const amrex::Real& prob_lo_z = geomdata.ProbLo()[2];
  const amrex::Real& prob_hi_z = geomdata.ProbHi()[2];

  const amrex::Real rdOcp   = sc.rdOcp;

  // const amrex::Real thetatr = 343.0;
  // const amrex::Real theta0  = 300.0;
  const amrex::Real ztr     = 12000.0;
  const amrex::Real Ttr     = 213.0;
  const amrex::Real Ttop    = 213.0;
  const amrex::Real deltaz  = 1000.*0.0;
  const amrex::Real zs      = 5000.;
  const amrex::Real us      = 30.;
  const amrex::Real uc      = 15.;

  amrex::ParallelForRNG(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
  {
    // Geometry (note we must include these here to get the data on device)
    const auto prob_lo         = geomdata.ProbLo();
    const auto dx              = geomdata.CellSize();
    const amrex::Real z        = prob_lo[2] + (k + 0.5) * dx[2];

    amrex::Real relhum = init_supercell_relhum(z, ztr);
    amrex::Real temp   = init_supercell_temperature(z, prob_lo_z, ztr, prob_hi_z, parms.T_0, Ttr, Ttop);
    amrex::Real press  = init_supercell_pressure(z, prob_lo_z, ztr, prob_hi_z, parms.T_0, Ttr, Ttop);
    amrex::Real qvs    = init_supercell_sat_mix(press, temp);

#if 0
    amrex::Real thetaeq;
    amrex::Real faceq;
    if (z <= ztr) {
       thetaeq = theta0 + (thetatr - theta0)*std::pow(z/ztr, 5./4.);
    }
    else {
       thetaeq = thetatr*std::exp(CONST_GRAV*(z-ztr)/c_p*Ttr);
    }

    faceq = 1.;
    for (int km = 0; km < k; ++km) {
       amrex::Real zloc = prob_lo[2]+(km+0.5)*dx[2];
       amrex::Real theq;
      if (zloc <= ztr) {
       theq = theta0 + (thetatr - theta0)*std::pow(zloc/ztr, 5./4.);
      } else {
       theq = thetatr*std::exp(CONST_GRAV*(zloc-ztr)/c_p*Ttr);
      }
       faceq -= CONST_GRAV/(c_p*theq)*dz;
    }

    temp = faceq*thetaeq;
#endif

    if (relhum*qvs > 0.014) relhum = 0.014/qvs;
    amrex::Real qvapor  = std::min(0.014, qvs*relhum);
    amrex::Real rho     = press/(R_d+qvapor*R_v)/temp;

    // perturb theta
    amrex::Real rand_double = amrex::Random(engine) - 1.0;        // Random number in [-1,1]
    amrex::Real scaling = (khi-static_cast<amrex::Real>(k))/khi;  // Less effect at higher levels
    amrex::Real deltaT = parms.T_pert*scaling*rand_double;

    amrex::Real theta = getThgivenRandT(rho, temp+deltaT, rdOcp);

    // This version perturbs rho but not p
    state(i, j, k, RhoTheta_comp) = rho*theta - getRhoThetagivenP(p_hse(i,j,k));
    state(i, j, k, Rho_comp)      = rho - r_hse(i,j,k);

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;

    // mean states
    if (use_moisture) {
        state(i, j, k, RhoQ1_comp) = rho*qvapor;
        state(i, j, k, RhoQ2_comp) = 0.0;
      }
  });

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real z = prob_lo_z + (k+0.5) * dz;
    if (z < zs-deltaz) {
      x_vel(i, j, k) = us*(z/zs) - uc;
    } else if (std::abs(z-zs) < deltaz) {
      x_vel(i, j, k) = (-0.8+3.*(z/zs)-1.25*(z/zs)*(z/zs))*us-uc;
    } else {
      x_vel(i, j, k) = us-uc;
    }
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_terrain(
    const Geometry& /*geom*/,
    MultiFab& z_phys_nd,
    const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        amrex::Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}

void
Problem::erf_init_rayleigh(
    amrex::Vector<amrex::Real>& tau,
    amrex::Vector<amrex::Real>& ubar,
    amrex::Vector<amrex::Real>& vbar,
    amrex::Vector<amrex::Real>& wbar,
    amrex::Vector<amrex::Real>& thetabar,
    amrex::Geometry      const& geom)
{
  const int khi = geom.Domain().bigEnd()[2];

  // We just use these values to test the Rayleigh damping
  for (int k = 0; k <= khi; k++)
  {
      tau[k]  = 1.0;
      ubar[k] = 2.0;
      vbar[k] = 1.0;
      wbar[k] = 0.0;
      thetabar[k] = parms.T_0;
  }
}
