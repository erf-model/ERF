#include <ERF.H>
#include <WindFarm.H>
#include <Fitch.H>
#include <IndexDefines.H>

using namespace amrex;

Real C_TKE = 0.0;

Real compute_A(Real z)
{

    Real d  = std::min(std::fabs(z - hub_height), rotor_rad);
    Real theta = std::acos(d/rotor_rad);
    Real A_s = rotor_rad*rotor_rad*theta - d*std::pow(rotor_rad*rotor_rad - d*d, 0.5);
    Real A = M_PI*rotor_rad*rotor_rad/2.0 - A_s;

    return A;
}

Real compute_Aijk(Real z_k, Real z_kp1)
{

    Real A_k   = compute_A(z_k);
    Real A_kp1 = compute_A(z_kp1);

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


void fitch_advance (int lev,
                    const Geometry& geom,
                    const Real& dt_advance,
                    MultiFab& cons_in,
                    MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                    MultiFab& mf_vars_fitch, const amrex::MultiFab& mf_Nturb)
{
    fitch_source_terms_cellcentered(geom, cons_in, U_old, V_old, W_old, mf_vars_fitch, mf_Nturb);
    fitch_update(dt_advance, cons_in, U_old, V_old, mf_vars_fitch);
}


void fitch_update (const Real& dt_advance,
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

void fitch_source_terms_cellcentered (const Geometry& geom,
                                      const MultiFab& cons_in,
                                      const MultiFab& U_old, const MultiFab& V_old, const MultiFab& W_old,
                                      MultiFab& mf_vars_fitch, const amrex::MultiFab& mf_Nturb)
{

  auto dx = geom.CellSizeArray();
  auto ProbHiArr = geom.ProbHiArray();
  auto ProbLoArr = geom.ProbLoArray();

  // Domain valid box
  const amrex::Box& domain = geom.Domain();
  int domlo_x = domain.smallEnd(0);
  int domhi_x = domain.bigEnd(0) + 1;
  int domlo_y = domain.smallEnd(1);
  int domhi_y = domain.bigEnd(1) + 1;
  int domlo_z = domain.smallEnd(2);
  int domhi_z = domain.bigEnd(2) + 1;


  Real sum = 0.0;
  Real *sum_area = &sum;

  // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
  mf_vars_fitch.setVal(0.0);

  for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox(1);
        auto cons_array  = cons_in.array(mfi);
        auto fitch_array = mf_vars_fitch.array(mfi);
        auto Nturb_array = mf_Nturb.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);
        auto w_vel       = W_old.array(mfi);

        amrex::IntVect lo = bx.smallEnd();

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
            int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);


            Real x = (ii+0.5) * dx[0];
            Real y = (jj+0.5) * dx[1];
            Real z = (kk+0.5) * dx[2];

            Real z_k   = kk*dx[2];
            Real z_kp1 = (kk+1)*dx[2];

            Real A_ijk = compute_Aijk(z_k, z_kp1);

            // Compute Fitch source terms

            Real Vabs = std::pow(u_vel(i,j,k)*u_vel(i,j,k) +
                                 v_vel(i,j,k)*v_vel(i,j,k) +
                                 w_vel(i,j,kk)*w_vel(i,j,kk), 0.5);

            Real C_T = interpolate_1d(wind_speed.dataPtr(), thrust_coeff.dataPtr(), z, wind_speed.size());

            fitch_array(i,j,k,0) = Vabs;
            fitch_array(i,j,k,1) =  -0.5*Nturb_array(i,j,k)/(dx[0]*dx[1])*C_T*Vabs*Vabs*A_ijk/(z_kp1 - z_k);
            fitch_array(i,j,k,2) = u_vel(i,j,k)/Vabs*fitch_array(i,j,k,1);
            fitch_array(i,j,k,3) = v_vel(i,j,k)/Vabs*fitch_array(i,j,k,1);
            fitch_array(i,j,k,4) = 0.5*Nturb_array(i,j,k)/(dx[0]*dx[1])*C_TKE*std::pow(Vabs,3)*A_ijk/(z_kp1 - z_k);

                 //amrex::Gpu::Atomic::Add(sum_area, A_ijk);
        });
    }
        //std::cout << "Checking sum here...." <<"\n";
        //printf("%0.15g, %0.15g\n", *sum_area , M_PI*R*R);
        //exit(0);
}
