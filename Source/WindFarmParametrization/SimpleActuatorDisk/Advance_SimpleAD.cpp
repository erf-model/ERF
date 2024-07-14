#include <SimpleAD.H>
#include <IndexDefines.H>

using namespace amrex;

void SimpleAD::advance (int lev,
                  const Geometry& geom,
                  const Real& dt_advance,
                  MultiFab& cons_in,
                  MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                  MultiFab& mf_vars_simpleAD, const amrex::MultiFab& mf_Nturb)
{

	std::cout << "Values of " << this->rotor_rad << " " << this->hub_height << "\n";

	exit(0);

    source_terms_cellcentered(geom, cons_in, U_old, V_old, W_old, mf_vars_simpleAD, mf_Nturb);
    update(dt_advance, cons_in, U_old, V_old, mf_vars_simpleAD);
}

void SimpleAD::update (const Real& dt_advance,
                  MultiFab& cons_in,
                  MultiFab& U_old, MultiFab& V_old,
                  const MultiFab& mf_vars_simpleAD)
{

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx  = mfi.tilebox();
        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        auto cons_array  = cons_in.array(mfi);
        auto simpleAD_array = mf_vars_simpleAD.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_vel(i,j,k) = u_vel(i,j,k) + (simpleAD_array(i-1,j,k,0) + simpleAD_array(i,j,k,0))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            v_vel(i,j,k) = v_vel(i,j,k) + (simpleAD_array(i,j-1,k,1) + simpleAD_array(i,j,k,1))/2.0*dt_advance;
        });
    }
}

void SimpleAD::source_terms_cellcentered (const Geometry& geom,
                                          const MultiFab& cons_in,
                                          const MultiFab& U_old, const MultiFab& V_old, const MultiFab& W_old,
                                          MultiFab& mf_vars_simpleAD, const amrex::MultiFab& mf_Nturb)
{

  auto dx = geom.CellSizeArray();
  auto ProbHiArr = geom.ProbHiArray();
  auto ProbLoArr = geom.ProbLoArray();
  Real d_rotor_rad = rotor_rad;
  Real d_hub_height = hub_height;

  // Domain valid box
  const amrex::Box& domain = geom.Domain();
  int domlo_x = domain.smallEnd(0);
  int domhi_x = domain.bigEnd(0) + 1;
  int domlo_y = domain.smallEnd(1);
  int domhi_y = domain.bigEnd(1) + 1;
  int domlo_z = domain.smallEnd(2);
  int domhi_z = domain.bigEnd(2) + 1;

  // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
  mf_vars_simpleAD.setVal(0.0);

  for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx       = mfi.tilebox();
        const Box& gbx      = mfi.growntilebox(1);
        auto cons_array     = cons_in.array(mfi);
        auto simpleAD_array = mf_vars_simpleAD.array(mfi);
        auto u_vel          = U_old.array(mfi);
        auto v_vel          = V_old.array(mfi);
        auto w_vel          = W_old.array(mfi);

        amrex::IntVect lo = bx.smallEnd();

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
            int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

            Real x1 = ProbLoArr[0] + ii     * dx[0];
            Real x2 = ProbLoArr[0] + (ii+1) * dx[0];

            Real x = ProbLoArr[0] + (ii+0.5) * dx[0];
            Real y = ProbLoArr[1] + (jj+0.5) * dx[1];
            Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

            // Compute Simple AD source terms

            Real Vabs = std::pow(u_vel(i,j,k)*u_vel(i,j,k) +
                                 v_vel(i,j,k)*v_vel(i,j,k) +
                                 w_vel(i,j,kk)*w_vel(i,j,kk), 0.5);

            Real phi = std::atan2(v_vel(i,j,k),u_vel(i,j,k)); // Wind direction w.r.t the x-dreiction


            Real fac = 0.0;
            int check_int = 0;
            for(int it=0;it<xloc.size();it++){
                if(xloc[it]+1e-12 > x1 and xloc[it]+1e-12 < x2) {
                   if(std::pow((y-yloc[it])*(y-yloc[it]) + (z-d_hub_height)*(z-d_hub_height),0.5) < d_rotor_rad) {
                        check_int++;
                        fac = -2.0*std::pow(u_vel(i,j,k)*std::cos(phi) + v_vel(i,j,k)*std::sin(phi), 2.0)*0.5*(1.0-0.5);
                    }
                }
            }
            if(check_int > 1){
                amrex::Error("Actuator disks are overlapping. Visualize actuator_disks.vtk "
                             "and check the windturbine locations input file. Exiting..");
            }

            simpleAD_array(i,j,k,0) = fac*std::cos(phi);
            simpleAD_array(i,j,k,1) = fac*std::sin(phi);
         });
    }
}
