#include <EWP.H>
#include <IndexDefines.H>
using namespace amrex;

Real R_ewp = 30.0;
Real z_c_ewp = 100.0;
Real C_T_ewp = 0.8, C_TKE_ewp = 0.0;

void ewp_advance (int lev,
                  const Geometry& geom,
                  const Real& dt_advance,
                  MultiFab& cons_in,
                  MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                  MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb)
 {
    ewp_source_terms_cellcentered(geom, cons_in, U_old, V_old, W_old, mf_vars_ewp, mf_Nturb);
    ewp_update(dt_advance, cons_in, U_old, V_old, mf_vars_ewp);
}


void ewp_update (const Real& dt_advance,
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

void ewp_source_terms_cellcentered (const Geometry& geom,
                                      const MultiFab& cons_in,
                                      const MultiFab& U_old, const MultiFab& V_old, const MultiFab& W_old,
                                      MultiFab& mf_vars_ewp, const amrex::MultiFab& mf_Nturb)
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

  // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
  mf_vars_ewp.setVal(0.0);

  for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox(1);
        auto cons_array  = cons_in.array(mfi);
        auto ewp_array = mf_vars_ewp.array(mfi);
        auto Nturb_array = mf_Nturb.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);
        auto w_vel       = W_old.array(mfi);

        amrex::IntVect lo = bx.smallEnd();

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
            int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);


            Real x = ProbLoArr[0] + (ii+0.5) * dx[0];
            Real y = ProbLoArr[1] + (jj+0.5) * dx[1];
            Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

            // Compute Fitch source terms

            Real Vabs = std::pow(u_vel(i,j,k)*u_vel(i,j,k) +
                                 v_vel(i,j,k)*v_vel(i,j,k) +
                                 w_vel(i,j,kk)*w_vel(i,j,kk), 0.5);

            Real L_wake = std::pow(dx[0]*dx[1],0.5);
            Real K_turb = 100.0;
            Real sigma_0 = 60.0;
            Real sigma_e = Vabs/(3.0*K_turb*L_wake)*(std::pow(2.0*K_turb*L_wake/Vabs + std::pow(sigma_0,2),3.0/2.0) - std::pow(sigma_0,3));

            Real phi     = std::atan2(v_vel(i,j,k),u_vel(i,j,k)); // Wind direction w.r.t the x-dreiction
            Real fac = -std::pow(M_PI/8.0,0.5)*C_T_ewp*std::pow(R_ewp,2)*std::pow(Vabs,2)/(dx[0]*dx[1]*sigma_e)*std::exp(-0.5*std::pow((z - z_c_ewp)/sigma_e,2));
            ewp_array(i,j,k,0) = fac*std::cos(phi)*Nturb_array(i,j,k);
            ewp_array(i,j,k,1) = fac*std::sin(phi)*Nturb_array(i,j,k);
            ewp_array(i,j,k,2) = 0.0;
         });
    }
}
