#include <ERF_Fitch.H>
#include <ERF_IndexDefines.H>
#include <ERF_Constants.H>
#include <ERF_Interpolation_1D.H>

using namespace amrex;

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
Real compute_A(const Real z,
               const Real hub_height,
               const Real rotor_rad)
{

    Real d  = std::min(std::fabs(z - hub_height), rotor_rad);
    Real theta = std::acos(d/rotor_rad);
    Real A_s = rotor_rad*rotor_rad*theta - d*std::pow(std::max(rotor_rad*rotor_rad - d*d,0.0), 0.5);
    Real A = PI*rotor_rad*rotor_rad/2.0 - A_s;

    return A;
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
Real compute_Aijk(const Real z_k,
                  const Real z_kp1,
                  const Real hub_height,
                  const Real rotor_rad)
{

    Real A_k   = compute_A(z_k, hub_height, rotor_rad);
    Real A_kp1 = compute_A(z_kp1, hub_height, rotor_rad);

    Real check = (z_k - hub_height)*(z_kp1 - hub_height);
    Real A_ijk;
    if(check > 0){
        A_ijk = std::fabs(A_k -A_kp1);
    }
    else{
        A_ijk = A_k + A_kp1;
    }

    return A_ijk;
}


void
Fitch::advance (const Geometry& geom,
                const Real& dt_advance,
                MultiFab& cons_in,
                MultiFab& mf_vars_fitch,
                MultiFab& U_old,
                MultiFab& V_old,
                MultiFab& W_old,
                const MultiFab& mf_Nturb,
                const MultiFab& mf_SMark)
{
    AMREX_ALWAYS_ASSERT(W_old.nComp() > 0);
    AMREX_ALWAYS_ASSERT(mf_SMark.nComp() > 0);
    source_terms_cellcentered(geom, cons_in, mf_vars_fitch, U_old, V_old, W_old, mf_Nturb);
    update(dt_advance, cons_in, U_old, V_old, mf_vars_fitch);
}


void
Fitch::update (const Real& dt_advance,
               MultiFab& cons_in,
               MultiFab& U_old, MultiFab& V_old,
               const MultiFab& mf_vars_fitch)
{

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx  = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        auto cons_array  = cons_in.array(mfi);
        auto fitch_array = mf_vars_fitch.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);

        ParallelFor(tbx, tby, bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_vel(i,j,k) = u_vel(i,j,k) + (fitch_array(i-1,j,k,2) + fitch_array(i,j,k,2))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            v_vel(i,j,k) = v_vel(i,j,k) + (fitch_array(i,j-1,k,3) + fitch_array(i,j,k,3))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            cons_array(i,j,k,RhoQKE_comp) = cons_array(i,j,k,RhoQKE_comp) + fitch_array(i,j,k,4)*dt_advance;
        });
    }
}

void
Fitch::source_terms_cellcentered (const Geometry& geom,
                                  const MultiFab& cons_in,
                                  MultiFab& mf_vars_fitch,
                                  const MultiFab& U_old,
                                  const MultiFab& V_old,
                                  const MultiFab& W_old,
                                  const MultiFab& mf_Nturb)
{

  get_turb_spec(rotor_rad, hub_height, thrust_coeff_standing,
                  wind_speed, thrust_coeff, power);

  auto dx = geom.CellSizeArray();

  // Domain valid box
  const amrex::Box& domain = geom.Domain();
  int domlo_z = domain.smallEnd(2);
  int domhi_z = domain.bigEnd(2) + 1;


  //Real sum = 0.0;
  //Real *sum_area = &sum;

  // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
  mf_vars_fitch.setVal(0.0);
  Real d_hub_height = hub_height;
  Real d_rotor_rad = rotor_rad;
  Gpu::DeviceVector<Real> d_wind_speed(wind_speed.size());
  Gpu::DeviceVector<Real> d_thrust_coeff(thrust_coeff.size());

  // Copy data from host vectors to device vectors
  Gpu::copy(Gpu::hostToDevice, wind_speed.begin(), wind_speed.end(), d_wind_speed.begin());
  Gpu::copy(Gpu::hostToDevice, thrust_coeff.begin(), thrust_coeff.end(), d_thrust_coeff.begin());

  for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& gbx = mfi.growntilebox(1);
        auto fitch_array = mf_vars_fitch.array(mfi);
        auto Nturb_array = mf_Nturb.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);
        auto w_vel       = W_old.array(mfi);

        const Real* wind_speed_d     = d_wind_speed.dataPtr();
        const Real* thrust_coeff_d   = d_thrust_coeff.dataPtr();
        const int n_spec_table = d_wind_speed.size();

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

            Real z_k   = kk*dx[2];
            Real z_kp1 = (kk+1)*dx[2];

            Real A_ijk = compute_Aijk(z_k, z_kp1, d_hub_height, d_rotor_rad);

            // Compute Fitch source terms

            Real Vabs = std::pow(u_vel(i,j,k)*u_vel(i,j,k) +
                                 v_vel(i,j,k)*v_vel(i,j,k) +
                                 w_vel(i,j,kk)*w_vel(i,j,kk), 0.5);

            Real C_T = interpolate_1d(wind_speed_d, thrust_coeff_d, Vabs, n_spec_table);
            Real C_TKE = 0.0;

            fitch_array(i,j,k,0) = Vabs;
            fitch_array(i,j,k,1) =  -0.5*Nturb_array(i,j,k)/(dx[0]*dx[1])*C_T*Vabs*Vabs*A_ijk/(z_kp1 - z_k);
            fitch_array(i,j,k,2) = u_vel(i,j,k)/Vabs*fitch_array(i,j,k,1);
            fitch_array(i,j,k,3) = v_vel(i,j,k)/Vabs*fitch_array(i,j,k,1);
            fitch_array(i,j,k,4) = 0.5*Nturb_array(i,j,k)/(dx[0]*dx[1])*C_TKE*std::pow(Vabs,3)*A_ijk/(z_kp1 - z_k);

                 //amrex::Gpu::Atomic::Add(sum_area, A_ijk);
        });
    }
        //std::cout << "Checking sum here...." <<"\n";
        //printf("%0.15g, %0.15g\n", *sum_area , PI*R*R);
        //exit(0);
}
