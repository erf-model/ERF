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
Real
init_supercell_temperature(Real z, Real z_0, Real z_trop, Real z_top,
                           Real T_0, Real T_trop, Real T_top)
{
    if (z <= z_trop) {
        Real lapse = - (T_trop - T_0) / (z_trop - z_0);
        return T_0 - lapse * (z - z_0);
    } else {
        Real lapse = - (T_top - T_trop) / (z_top - z_trop);
        return T_trop - lapse * (z - z_trop);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
init_supercell_pressure(Real z, Real z_0,
                        Real z_trop, Real z_top,
                        Real T_0, Real T_trop, Real T_top)
{
    if (z <= z_trop) {
        Real lapse = - (T_trop - T_0) / (z_trop - z_0);
        Real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
        return p_0 * std::pow( T / T_0 , CONST_GRAV/(R_d*lapse) );
    } else {
        // Get pressure at the tropopause
        Real lapse = - (T_trop - T_0) / (z_trop - z_0);
        Real p_trop = p_0 * std::pow( T_trop / T_0 , CONST_GRAV/(R_d*lapse) );
        // Get pressure at requested height
        lapse = - (T_top - T_trop) / (z_top - z_trop);
        if (lapse != 0) {
            Real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
            return p_trop * std::pow( T / T_trop , CONST_GRAV/(R_d*lapse) );
        } else {
            return p_trop * std::exp(-CONST_GRAV*(z-z_trop)/(R_d*T_trop));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
init_supercell_sat_mix(Real press, Real T ) {
    return 380./(press) * std::exp(17.27*(T-273.)/(T-36.));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
init_supercell_relhum(Real z, Real z_trop)
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
    amrex::Array4<amrex::Real const> const& /*state*/,
    amrex::Array4<amrex::Real      > const& state_pert,
    amrex::Array4<amrex::Real      > const& x_vel_pert,
    amrex::Array4<amrex::Real      > const& y_vel_pert,
    amrex::Array4<amrex::Real      > const& z_vel_pert,
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
    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
    const Real& dz        = geomdata.CellSize()[2];
    const Real& prob_lo_z = geomdata.ProbLo()[2];
    const Real& prob_hi_z = geomdata.ProbHi()[2];

    const Real rdOcp   = sc.rdOcp;

    // const Real thetatr = 343.0;
    // const Real theta0  = 300.0;
    const Real ztr     = 12000.0;
    const Real Ttr     = 213.0;
    const Real Ttop    = 213.0;
    const Real deltaz  = 1000.*0.0;
    const Real zs      = 5000.;
    const Real us      = 30.;
    const Real uc      = 15.;

    ParallelForRNG(bx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
    {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo         = geomdata.ProbLo();
        const auto dx              = geomdata.CellSize();
        const Real z        = prob_lo[2] + (k + 0.5) * dx[2];

        Real relhum = init_supercell_relhum(z, ztr);
        Real temp   = init_supercell_temperature(z, prob_lo_z, ztr, prob_hi_z, parms_d.T_0, Ttr, Ttop);
        Real press  = init_supercell_pressure(z, prob_lo_z, ztr, prob_hi_z, parms_d.T_0, Ttr, Ttop);
        Real qvs    = init_supercell_sat_mix(press, temp);

#if 0
        Real thetaeq;
        Real faceq;
        if (z <= ztr) {
           thetaeq = theta0 + (thetatr - theta0)*std::pow(z/ztr, 5./4.);
        } else {
           thetaeq = thetatr*std::exp(CONST_GRAV*(z-ztr)/c_p*Ttr);
        }

        faceq = 1.;
        for (int km = 0; km < k; ++km) {
           Real zloc = prob_lo[2]+(km+0.5)*dx[2];
           Real theq;
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
        Real qvapor  = std::min(0.014, qvs*relhum);
        Real rho     = press/(R_d+qvapor*R_v)/temp;

        // perturb theta
        Real rand_double = amrex::Random(engine) - 1.0;        // Random number in [-1,1]
        Real scaling = (khi-static_cast<Real>(k))/khi;  // Less effect at higher levels
        Real deltaT = parms_d.T_pert*scaling*rand_double;

        Real theta = getThgivenRandT(rho, temp+deltaT, rdOcp);

        // This version perturbs rho but not p
        state_pert(i, j, k, RhoTheta_comp) = rho*theta - getRhoThetagivenP(p_hse(i,j,k));
        state_pert(i, j, k, Rho_comp)      = rho - r_hse(i,j,k);

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        // mean states
        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = rho*qvapor;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    // Set the x-velocity
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real z = prob_lo_z + (k+0.5) * dz;
        if (z < zs-deltaz) {
            x_vel_pert(i, j, k) = us*(z/zs) - uc;
        } else if (std::abs(z-zs) < deltaz) {
            x_vel_pert(i, j, k) = (-0.8+3.*(z/zs)-1.25*(z/zs)*(z/zs))*us-uc;
        } else {
            x_vel_pert(i, j, k) = us-uc;
        }
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        y_vel_pert(i, j, k) = 0.0;
    });

    // Set the z-velocity
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        z_vel_pert(i, j, k) = 0.0;
    });

    amrex::Gpu::streamSynchronize();
}

void
Problem::erf_init_rayleigh(
    amrex::Vector<amrex::Vector<Real> >& rayleigh_ptrs,
    amrex::Geometry      const& geom,
    std::unique_ptr<MultiFab>& /*z_phys_cc*/)
{
    const int khi = geom.Domain().bigEnd()[2];

    // We just use these values to test the Rayleigh damping
    for (int k = 0; k <= khi; k++)
    {
        rayleigh_ptrs[Rayleigh::tau][k]      = 1.0;
        rayleigh_ptrs[Rayleigh::ubar][k]     = 2.0;
        rayleigh_ptrs[Rayleigh::vbar][k]     = 1.0;
        rayleigh_ptrs[Rayleigh::wbar][k]     = 0.0;
        rayleigh_ptrs[Rayleigh::thetabar][k] = parms.T_0;
    }
}
