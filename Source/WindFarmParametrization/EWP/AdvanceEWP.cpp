#include <EWP.H>
#include <IndexDefines.H>
#include <ERF_Constants.H>
#include <Interpolation_1D.H>

using namespace amrex;

void
EWP::advance (const Geometry& geom,
              const Real& dt_advance,
              MultiFab& cons_in,
              MultiFab& mf_vars_ewp,
              MultiFab& U_old,
              MultiFab& V_old,
              MultiFab& W_old,
              const MultiFab& mf_Nturb)
 {
    source_terms_cellcentered(geom, cons_in, mf_vars_ewp, U_old, V_old, W_old, mf_Nturb);
    update(dt_advance, cons_in, U_old, V_old, mf_vars_ewp);
}


void
EWP::update (const Real& dt_advance,
             MultiFab& cons_in,
             MultiFab& U_old, MultiFab& V_old,
             const MultiFab& mf_vars_ewp)
{

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx  = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        auto cons_array  = cons_in.array(mfi);
        auto ewp_array = mf_vars_ewp.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);

        ParallelFor(tbx, tby, bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_vel(i,j,k) = u_vel(i,j,k) + (ewp_array(i-1,j,k,0) + ewp_array(i,j,k,0))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            v_vel(i,j,k) = v_vel(i,j,k) + (ewp_array(i,j-1,k,1) + ewp_array(i,j,k,1))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            cons_array(i,j,k,RhoQKE_comp) = cons_array(i,j,k,RhoQKE_comp) + ewp_array(i,j,k,2)*dt_advance;
        });
    }
}

void
EWP::source_terms_cellcentered (const Geometry& geom,
                                const MultiFab& cons_in,
                                MultiFab& mf_vars_ewp,
                                const MultiFab& U_old,
                                const MultiFab& V_old,
                                const MultiFab& W_old,
                                const MultiFab& mf_Nturb)
{

  get_turb_spec(rotor_rad, hub_height, thrust_coeff_standing,
                  wind_speed, thrust_coeff, power);

  auto dx = geom.CellSizeArray();
  auto ProbLoArr = geom.ProbLoArray();
  Real sigma_0 = 1.7*rotor_rad;

  Real d_rotor_rad = rotor_rad;
  Real d_hub_height = hub_height;

  Gpu::DeviceVector<Real> d_wind_speed(wind_speed.size());
  Gpu::DeviceVector<Real> d_thrust_coeff(thrust_coeff.size());

  // Copy data from host vectors to device vectors
  Gpu::copy(Gpu::hostToDevice, wind_speed.begin(), wind_speed.end(), d_wind_speed.begin());
  Gpu::copy(Gpu::hostToDevice, thrust_coeff.begin(), thrust_coeff.end(), d_thrust_coeff.begin());


  // Domain valid box
  const amrex::Box& domain = geom.Domain();
  int domlo_z = domain.smallEnd(2);
  int domhi_z = domain.bigEnd(2) + 1;

  // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
  mf_vars_ewp.setVal(0.0);

  for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& gbx = mfi.growntilebox(1);
        auto ewp_array = mf_vars_ewp.array(mfi);
        auto Nturb_array = mf_Nturb.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);
        auto w_vel       = W_old.array(mfi);

        const Real* wind_speed_d     = d_wind_speed.dataPtr();
        const Real* thrust_coeff_d   = d_thrust_coeff.dataPtr();
        const int n_spec_table = d_wind_speed.size();

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);
            Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

            // Compute Fitch source terms

            Real Vabs = std::pow(u_vel(i,j,k)*u_vel(i,j,k) +
                                 v_vel(i,j,k)*v_vel(i,j,k) +
                                 w_vel(i,j,kk)*w_vel(i,j,kk), 0.5);

            Real C_T = interpolate_1d(wind_speed_d, thrust_coeff_d, Vabs, n_spec_table);

            Real C_TKE = 0.0;
            Real K_turb = 1.0;

            Real L_wake = std::pow(dx[0]*dx[1],0.5)/2.0;
            Real sigma_e = Vabs/(3.0*K_turb*L_wake)*
                           (std::pow(2.0*K_turb*L_wake/Vabs + std::pow(sigma_0,2),3.0/2.0) - std::pow(sigma_0,3));

            Real phi     = std::atan2(v_vel(i,j,k),u_vel(i,j,k)); // Wind direction w.r.t the x-dreiction
            Real fac = -std::pow(PI/8.0,0.5)*C_T*std::pow(d_rotor_rad,2)*
                        std::pow(Vabs,2)/(dx[0]*dx[1]*sigma_e)*
                        std::exp(-0.5*std::pow((z - d_hub_height)/sigma_e,2));
            ewp_array(i,j,k,0) = fac*std::cos(phi)*Nturb_array(i,j,k);
            ewp_array(i,j,k,1) = fac*std::sin(phi)*Nturb_array(i,j,k);
            ewp_array(i,j,k,2) = C_TKE*0.0;
         });
    }
}
