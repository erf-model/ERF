#include <ERF_GeneralAD.H>
#include <ERF_IndexDefines.H>
#include <ERF_Interpolation_1D.H>

using namespace amrex;

void
GeneralAD::advance (const Geometry& geom,
                  const Real& dt_advance,
                  MultiFab& cons_in,
                  MultiFab& mf_vars_generalAD,
                  MultiFab& U_old,
                  MultiFab& V_old,
                  MultiFab& W_old,
                  const MultiFab& mf_Nturb,
                  const MultiFab& mf_SMark)
{
    AMREX_ALWAYS_ASSERT(W_old.nComp() > 0);
    AMREX_ALWAYS_ASSERT(mf_Nturb.nComp() > 0);
    AMREX_ALWAYS_ASSERT(mf_vars_generalAD.nComp() > 0);
    source_terms_cellcentered(geom, cons_in, mf_SMark, mf_vars_generalAD);
    update(dt_advance, cons_in, U_old, V_old, mf_vars_generalAD);
}

void
GeneralAD::update (const Real& dt_advance,
                  MultiFab& cons_in,
                  MultiFab& U_old, MultiFab& V_old,
                  const MultiFab& mf_vars_generalAD)
{

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);

        auto generalAD_array = mf_vars_generalAD.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);

        ParallelFor(tbx, tby,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_vel(i,j,k) = u_vel(i,j,k) + (generalAD_array(i-1,j,k,0) + generalAD_array(i,j,k,0))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            v_vel(i,j,k) = v_vel(i,j,k) + (generalAD_array(i,j-1,k,1) + generalAD_array(i,j,k,1))/2.0*dt_advance;
        });
    }
}

void
GeneralAD::source_terms_cellcentered (const Geometry& geom,
                                     const MultiFab& cons_in,
                                     const MultiFab& mf_SMark,
                                     MultiFab& mf_vars_generalAD)
{

    get_turb_loc(xloc, yloc);
    get_turb_spec(rotor_rad, hub_height, thrust_coeff_standing,
                  wind_speed, thrust_coeff, power);

    Gpu::DeviceVector<Real> d_xloc(xloc.size());
    Gpu::DeviceVector<Real> d_yloc(yloc.size());
    Gpu::copy(Gpu::hostToDevice, xloc.begin(), xloc.end(), d_xloc.begin());
    Gpu::copy(Gpu::hostToDevice, yloc.begin(), yloc.end(), d_yloc.begin());

      auto dx = geom.CellSizeArray();

  // Domain valid box
      const amrex::Box& domain = geom.Domain();
      int domlo_x = domain.smallEnd(0);
      int domhi_x = domain.bigEnd(0) + 1;
      int domlo_y = domain.smallEnd(1);
      int domhi_y = domain.bigEnd(1) + 1;
      int domlo_z = domain.smallEnd(2);
      int domhi_z = domain.bigEnd(2) + 1;

      // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
      mf_vars_generalAD.setVal(0.0);

      long unsigned int nturbs = xloc.size();

    get_turb_disk_angle(turb_disk_angle);
    Real nx = -std::cos(turb_disk_angle);
    Real ny = -std::sin(turb_disk_angle);
    Real d_turb_disk_angle = turb_disk_angle;

    Gpu::DeviceVector<Real> d_wind_speed(wind_speed.size());
    Gpu::DeviceVector<Real> d_thrust_coeff(thrust_coeff.size());

    // Copy data from host vectors to device vectors
    Gpu::copy(Gpu::hostToDevice, wind_speed.begin(), wind_speed.end(), d_wind_speed.begin());
    Gpu::copy(Gpu::hostToDevice, thrust_coeff.begin(), thrust_coeff.end(), d_thrust_coeff.begin());

    const Real* wind_speed_d     = d_wind_speed.dataPtr();
    const Real* thrust_coeff_d   = d_thrust_coeff.dataPtr();
    const int n_spec_table = d_wind_speed.size();

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& gbx      = mfi.growntilebox(1);
        auto SMark_array    = mf_SMark.array(mfi);
        auto generalAD_array = mf_vars_generalAD.array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
            int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

            int check_int = 0;

            Real source_x = 0.0;
            Real source_y = 0.0;

            for(long unsigned int it=0;it<nturbs;it++) {
                Real avg_vel  = 10.0;
                Real phi      = 0.0;

                Real C_T = interpolate_1d(wind_speed_d, thrust_coeff_d, avg_vel, n_spec_table);
                Real a;
                if(C_T <= 1) {
                    a = 0.5 - 0.5*std::pow(1.0-C_T,0.5);
                }
                Real Uinfty_dot_nhat = avg_vel*(std::cos(phi)*nx + std::sin(phi)*ny);
                // Add the source terms in the location of the actuator disk
                if(SMark_array(ii,jj,kk,1) == static_cast<double>(it)) {
                    check_int++;
                    if(C_T <= 1) {
                        source_x = -2.0*std::pow(Uinfty_dot_nhat, 2.0)*a*(1.0-a)*dx[1]*dx[2]*std::cos(d_turb_disk_angle)/(dx[0]*dx[1]*dx[2])*std::cos(phi);
                        source_y = -2.0*std::pow(Uinfty_dot_nhat, 2.0)*a*(1.0-a)*dx[1]*dx[2]*std::cos(d_turb_disk_angle)/(dx[0]*dx[1]*dx[2])*std::sin(phi);
                    }
                    else {
                        source_x = -0.5*C_T*std::pow(Uinfty_dot_nhat, 2.0)*dx[1]*dx[2]*std::cos(d_turb_disk_angle)/(dx[0]*dx[1]*dx[2])*std::cos(phi);
                        source_y = -0.5*C_T*std::pow(Uinfty_dot_nhat, 2.0)*dx[1]*dx[2]*std::cos(d_turb_disk_angle)/(dx[0]*dx[1]*dx[2])*std::sin(phi);
                    }
                }
            }
            if(check_int > 1){
                amrex::Error("Actuator disks are overlapping. Visualize actuator_disks.vtk "
                             "and check the windturbine locations input file. Exiting..");
            }

            generalAD_array(i,j,k,0) = source_x;
            generalAD_array(i,j,k,1) = source_y;
         });
    }
}
