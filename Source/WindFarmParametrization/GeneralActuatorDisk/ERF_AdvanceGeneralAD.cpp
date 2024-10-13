#include <ERF_GeneralAD.H>
#include <ERF_IndexDefines.H>
#include <ERF_Interpolation_1D.H>
#include <ERF_Constants.H>

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
    compute_freestream_velocity(cons_in, U_old, V_old, mf_SMark);
    source_terms_cellcentered(geom, cons_in, mf_SMark, mf_vars_generalAD);
    update(dt_advance, cons_in, U_old, V_old, W_old, mf_vars_generalAD);
}

void
GeneralAD::update (const Real& dt_advance,
                  MultiFab& cons_in,
                  MultiFab& U_old,
                  MultiFab& V_old,
                  MultiFab& W_old,
                  const MultiFab& mf_vars_generalAD)
{

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box tbx = mfi.nodaltilebox(0);
        Box tby = mfi.nodaltilebox(1);
        Box tbz = mfi.nodaltilebox(2);

        auto generalAD_array = mf_vars_generalAD.array(mfi);
        auto u_vel       = U_old.array(mfi);
        auto v_vel       = V_old.array(mfi);
        auto w_vel       = W_old.array(mfi);

        ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            u_vel(i,j,k) = u_vel(i,j,k) + (generalAD_array(i-1,j,k,0) + generalAD_array(i,j,k,0))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            v_vel(i,j,k) = v_vel(i,j,k) + (generalAD_array(i,j-1,k,1) + generalAD_array(i,j,k,1))/2.0*dt_advance;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            w_vel(i,j,k) = w_vel(i,j,k) + (generalAD_array(i,j-1,k,2) + generalAD_array(i,j,k,2))/2.0*dt_advance;
        });

    }
}

void GeneralAD::compute_freestream_velocity(const MultiFab& cons_in,
                                           const MultiFab& U_old,
                                           const MultiFab& V_old,
                                           const MultiFab& mf_SMark)
{
     get_turb_loc(xloc, yloc);
     freestream_velocity.clear();
     freestream_phi.clear();
     disk_cell_count.clear();
     freestream_velocity.resize(xloc.size(),0.0);
     freestream_phi.resize(xloc.size(),0.0);
     disk_cell_count.resize(xloc.size(),0.0);

     Gpu::DeviceVector<Real> d_freestream_velocity(xloc.size());
     Gpu::DeviceVector<Real> d_freestream_phi(yloc.size());
     Gpu::DeviceVector<Real> d_disk_cell_count(yloc.size());
     Gpu::copy(Gpu::hostToDevice, freestream_velocity.begin(), freestream_velocity.end(), d_freestream_velocity.begin());
     Gpu::copy(Gpu::hostToDevice, freestream_phi.begin(), freestream_phi.end(), d_freestream_phi.begin());
     Gpu::copy(Gpu::hostToDevice, disk_cell_count.begin(), disk_cell_count.end(), d_disk_cell_count.begin());

     Real* d_freestream_velocity_ptr = d_freestream_velocity.data();
     Real* d_freestream_phi_ptr = d_freestream_phi.data();
     Real* d_disk_cell_count_ptr     = d_disk_cell_count.data();


     for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        auto SMark_array    = mf_SMark.array(mfi);
        auto u_vel          = U_old.array(mfi);
        auto v_vel          = V_old.array(mfi);
        Box tbx = mfi.nodaltilebox(0);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            if(SMark_array(i,j,k,0) != -1.0) {
                int turb_index = static_cast<int>(SMark_array(i,j,k,0));
                Real phi = std::atan2(v_vel(i,j,k),u_vel(i,j,k)); // Wind direction w.r.t the x-dreiction
                Gpu::Atomic::Add(&d_freestream_velocity_ptr[turb_index],std::pow(u_vel(i,j,k)*u_vel(i,j,k) + v_vel(i,j,k)*v_vel(i,j,k),0.5));
                Gpu::Atomic::Add(&d_disk_cell_count_ptr[turb_index],1.0);
                Gpu::Atomic::Add(&d_freestream_phi_ptr[turb_index],phi);
            }
        });
    }

    // Copy back to host
    Gpu::copy(Gpu::deviceToHost, d_freestream_velocity.begin(), d_freestream_velocity.end(), freestream_velocity.begin());
    Gpu::copy(Gpu::deviceToHost, d_freestream_phi.begin(), d_freestream_phi.end(), freestream_phi.begin());
    Gpu::copy(Gpu::deviceToHost, d_disk_cell_count.begin(), d_disk_cell_count.end(), disk_cell_count.begin());

    // Reduce the data on every processor
    amrex::ParallelAllReduce::Sum(freestream_velocity.data(),
                                  freestream_velocity.size(),
                                  amrex::ParallelContext::CommunicatorAll());

    amrex::ParallelAllReduce::Sum(freestream_phi.data(),
                                  freestream_phi.size(),
                                  amrex::ParallelContext::CommunicatorAll());


   amrex::ParallelAllReduce::Sum(disk_cell_count.data(),
                                 disk_cell_count.size(),
                                 amrex::ParallelContext::CommunicatorAll());

    get_turb_loc(xloc, yloc);


    if (ParallelDescriptor::IOProcessor()){
        for(int it=0; it<xloc.size(); it++){
            //std::cout << "turbine index, freestream velocity is " << it << " " << freestream_velocity[it] << " " <<
            //                                                 disk_cell_count[it]  <<  " " <<
            //                                                freestream_velocity[it]/(disk_cell_count[it] + 1e-10) << " " <<
            //                                              freestream_phi[it]/(disk_cell_count[it] + 1e-10) << "\n";
        }
    }
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
int find_rad_loc_index(const Real rad,
                       const Real* bld_rad_loc,
                       const int n_bld_sections)
{
      // Find the index of the radial location
    int index=-1;
    Real rhub = 2.0;
    Real rad_from_hub = rad - rhub;
    if(rad_from_hub < 0.0) {
        index = 0;
    }
    else {
        for(int i=0;i<n_bld_sections;i++){
            if(bld_rad_loc[i] > rad) {
                index = i;
                break;
            }
        }
    }
    if(index == -1 and rad > bld_rad_loc[n_bld_sections-1]) {
        index = n_bld_sections-1;
    }
    if(index == -1) {
        //printf("The radial section is at %0.15g m\n",rad);
        Abort("Could not find index of the radial section.");
    }

    return index;
}

AMREX_FORCE_INLINE
AMREX_GPU_DEVICE
std::array<Real,2>
compute_source_terms_Fn_Ft (const Real rad,
                            const Real avg_vel,
                            const Real* bld_rad_loc,
                            const Real* bld_twist,
                            const Real* bld_chord,
                            int n_bld_sections,
                            const Real* bld_airfoil_aoa,
                            const Real* bld_airfoil_Cl,
                            const Real* bld_airfoil_Cd,
                            const int n_pts_airfoil)
{
    Real Omega = 9.0/60.0*2.0*PI;
    Real rho = 1.226;

    Real B = 3.0;
    Real rhub = 2.0;
    Real rtip = 63.5;

    Real twist = interpolate_1d(bld_rad_loc, bld_twist, rad, n_bld_sections);
    Real c = interpolate_1d(bld_rad_loc, bld_chord, rad, n_bld_sections);

    // Iteration procedure

    Real s = 0.5*c*B/(PI*rad);

    Real at, an, V1, Vt, Vr, psi, L, D, Cn, Ct;
    Real ftip, fhub, F, Cl, Cd, at_new, an_new;

    at = 0.1;
    an = 0.1;

    bool is_converged = false;

    for(int i=0;i<100;i++) {
        V1 = avg_vel*(1-an);
        Vt = Omega*(1.0+at)*rad;
        Vr = std::pow(V1*V1+Vt*Vt,0.5);

        psi = std::atan2(V1,Vt);

        Real aoa = psi*180.0/PI - twist;

        Cl = interpolate_1d(bld_airfoil_aoa, bld_airfoil_Cl, aoa, n_pts_airfoil);
        Cd = interpolate_1d(bld_airfoil_aoa, bld_airfoil_Cd, aoa, n_pts_airfoil);

        //Cl = 1.37;
        //Cd = 0.014;

        //printf("rad, aoa, Cl, Cd    = %0.15g %0.15g %0.15g %0.15g\n", rad, aoa, Cl, Cd);

        Cn = Cl*std::cos(psi) + Cd*std::sin(psi);
        Ct = Cl*std::sin(psi) - Cd*std::cos(psi);

        ftip = B*(rtip-rad)/(2.0*rad*std::sin(psi)+1e-10);
        fhub = B*(rad-rhub)/(2.0*rad*std::sin(psi)+1e-10);

        AMREX_ALWAYS_ASSERT(std::fabs(std::exp(-fhub))<=1.0);
        AMREX_ALWAYS_ASSERT(std::fabs(std::exp(-ftip))<=1.0);

        F = 1.0;//2.0/PI*(std::acos(std::exp(-ftip)) + std::acos(std::exp(-fhub)) );

        at_new = 1.0/ ( 4.0*F*std::sin(psi)*std::cos(psi)/(s*Ct+1e-10) - 1.0 );
        an_new = 1.0/ ( 1.0 + 4.0*F*std::pow(std::sin(psi),2)/(s*Cn + 1e-10) );
        at_new = std::max(0.0, at_new);

        if(std::fabs(at_new-at) < 1e-5 and std::fabs(an_new-an) < 1e-5) {
            //printf("Converged at, an = %d %0.15g %0.15g %0.15g\n",i, at, an, psi);
            at = at_new;
            an = an_new;
            is_converged = true;
            break;
        }
        at = at_new;
        an = an_new;
        //printf("Iteration, at, an = %0.15g %0.15g %0.15g\n",at, an, psi);
    }

    if(!is_converged) {
        Abort("The iteration procedure for the generalized actuator disk did not converge. Exiting...");
    }

    // Iterations converged. Now compute Ft, Fn

    L = 0.5*rho*Vr*Vr*c*Cl;
    D = 0.5*rho*Vr*Vr*c*Cd;

    Real Fn = L*std::cos(psi) + D*std::sin(psi);
    Real Ft = L*std::sin(psi) - D*std::cos(psi);

     //printf("Fn and Ft %0.15g %0.15g %0.15g %0.15g\n", L, D, std::cos(psi), std::sin(psi));

    std::array<Real, 2> Fn_and_Ft;
    Fn_and_Ft[0] = Fn;
    Fn_and_Ft[1] = Ft;

    return Fn_and_Ft;

    //exit(0);
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

    get_blade_spec(bld_rad_loc,bld_twist,bld_chord);

    get_blade_airfoil_spec(bld_airfoil_aoa, bld_airfoil_Cl, bld_airfoil_Cd);

    for(int i=0;i<bld_rad_loc.size();i++) {
        //std::cout << "I am here ...." << bld_rad_loc[i] << " " << bld_twist[i] << " " << bld_chord[i] << "\n";
    }


    Real d_hub_height = hub_height;
    Real d_rotor_rad = rotor_rad;

    Gpu::DeviceVector<Real> d_xloc(xloc.size());
    Gpu::DeviceVector<Real> d_yloc(yloc.size());
    Gpu::copy(Gpu::hostToDevice, xloc.begin(), xloc.end(), d_xloc.begin());
    Gpu::copy(Gpu::hostToDevice, yloc.begin(), yloc.end(), d_yloc.begin());

      auto dx = geom.CellSizeArray();

  // Domain valid box
      const amrex::Box& domain = geom.Domain();
      auto ProbLoArr = geom.ProbLoArray();
      int domlo_x = domain.smallEnd(0);
      int domhi_x = domain.bigEnd(0) + 1;
      int domlo_y = domain.smallEnd(1);
      int domhi_y = domain.bigEnd(1) + 1;
      int domlo_z = domain.smallEnd(2);
      int domhi_z = domain.bigEnd(2) + 1;

      // The order of variables are - Vabs dVabsdt, dudt, dvdt, dTKEdt
      mf_vars_generalAD.setVal(0.0);

     long unsigned int nturbs = xloc.size();

    // This is the angle phi in Fig. 10 in Mirocha et. al. 2014
    // set_turb_disk angle in ERF_InitWindFarm.cpp sets this phi as
    // the turb_disk_angle
    get_turb_disk_angle(turb_disk_angle);
    Real d_turb_disk_angle = turb_disk_angle;

    Gpu::DeviceVector<Real> d_freestream_velocity(nturbs);
    Gpu::DeviceVector<Real> d_disk_cell_count(nturbs);
    Gpu::copy(Gpu::hostToDevice, freestream_velocity.begin(), freestream_velocity.end(), d_freestream_velocity.begin());
    Gpu::copy(Gpu::hostToDevice, disk_cell_count.begin(), disk_cell_count.end(), d_disk_cell_count.begin());

    Real* d_xloc_ptr = d_xloc.data();
    Real* d_yloc_ptr = d_yloc.data();
    Real* d_freestream_velocity_ptr = d_freestream_velocity.data();
    Real* d_disk_cell_count_ptr     = d_disk_cell_count.data();

    int n_bld_sections = bld_rad_loc.size();

    Gpu::DeviceVector<Real>    d_bld_rad_loc(n_bld_sections);
    Gpu::DeviceVector<Real>    d_bld_twist(n_bld_sections);
    Gpu::DeviceVector<Real>    d_bld_chord(n_bld_sections);

    Gpu::copy(Gpu::hostToDevice, bld_rad_loc.begin(), bld_rad_loc.end(), d_bld_rad_loc.begin());
    Gpu::copy(Gpu::hostToDevice, bld_twist.begin(), bld_twist.end(), d_bld_twist.begin());
    Gpu::copy(Gpu::hostToDevice, bld_chord.begin(), bld_chord.end(), d_bld_chord.begin());

    Real* bld_rad_loc_ptr = d_bld_rad_loc.data();
    Real* bld_twist_ptr = d_bld_twist.data();
    Real* bld_chord_ptr = d_bld_chord.data();

    Vector<Gpu::DeviceVector<Real>> d_bld_airfoil_aoa(n_bld_sections);
    Vector<Gpu::DeviceVector<Real>> d_bld_airfoil_Cl(n_bld_sections);
    Vector<Gpu::DeviceVector<Real>> d_bld_airfoil_Cd(n_bld_sections);

    int n_pts_airfoil = bld_airfoil_aoa[0].size();

    for(int i=0;i<n_bld_sections;i++){
        d_bld_airfoil_aoa[i].resize(n_pts_airfoil);
        d_bld_airfoil_Cl[i].resize(n_pts_airfoil);
        d_bld_airfoil_Cd[i].resize(n_pts_airfoil);
        Gpu::copy(Gpu::hostToDevice, bld_airfoil_aoa[i].begin(), bld_airfoil_aoa[i].end(), d_bld_airfoil_aoa[i].begin());
        Gpu::copy(Gpu::hostToDevice, bld_airfoil_Cl[i].begin(), bld_airfoil_Cl[i].end(), d_bld_airfoil_Cl[i].begin());
        Gpu::copy(Gpu::hostToDevice, bld_airfoil_Cd[i].begin(), bld_airfoil_Cd[i].end(), d_bld_airfoil_Cd[i].begin());
    }

    Vector<Real*> hp_bld_airfoil_aoa, hp_bld_airfoil_Cl, hp_bld_airfoil_Cd;
    for (auto & v :d_bld_airfoil_aoa) {
        hp_bld_airfoil_aoa.push_back(v.data());
    }
    for (auto & v :d_bld_airfoil_Cl) {
        hp_bld_airfoil_Cl.push_back(v.data());
    }
    for (auto & v :d_bld_airfoil_Cd) {
        hp_bld_airfoil_Cd.push_back(v.data());
    }

    Gpu::AsyncArray<Real*> aoa(hp_bld_airfoil_aoa.data(), n_bld_sections);
    Gpu::AsyncArray<Real*> Cl(hp_bld_airfoil_Cl.data(), n_bld_sections);
    Gpu::AsyncArray<Real*> Cd(hp_bld_airfoil_Cd.data(), n_bld_sections);

    auto d_bld_airfoil_aoa_ptr = aoa.data();
    auto d_bld_airfoil_Cl_ptr  = Cl.data();
    auto d_bld_airfoil_Cd_ptr  = Cd.data();

    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& gbx      = mfi.growntilebox(1);
        auto SMark_array    = mf_SMark.array(mfi);
        auto generalAD_array = mf_vars_generalAD.array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
            int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
            int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

            Real x   = ProbLoArr[0] + (ii+0.5)*dx[0];
            Real y   = ProbLoArr[1] + (jj+0.5)*dx[1];
            Real z   = ProbLoArr[2] + (kk+0.5)*dx[2];
            // ?? Density needed here
            Real inv_dens_vol = 1.0/(1.0*dx[0]*dx[1]*dx[2]);

            int check_int = 0;

            Real source_x = 0.0, source_y = 0.0, source_z = 0.0;
            std::array<Real,2> Fn_and_Ft;

            for(long unsigned int it=0;it<nturbs;it++) {
                 Real avg_vel  = d_freestream_velocity_ptr[it]/(d_disk_cell_count_ptr[it] + 1e-10);
                 Real phi = d_turb_disk_angle;

                // This if check makes sure it is a point on the actuator disk
                if(SMark_array(ii,jj,kk,1) == static_cast<double>(it)) {
                    check_int++;

                    // Find radial distance of the point and the zeta angle
                    Real rad = std::pow( (x-d_xloc_ptr[it])*(x-d_xloc_ptr[it]) +
                                         (y-d_yloc_ptr[it])*(y-d_yloc_ptr[it]) +
                                         (z-d_hub_height)*(z-d_hub_height), 0.5 );

                    int index = find_rad_loc_index(rad, bld_rad_loc_ptr, n_bld_sections);

                    // This if check makes sure it is a point with radial distance
                    // between the hub radius and the rotor radius.
                    // ?? hub radius needed here
                    if(rad >= 2.0 and rad <= d_rotor_rad) {
                        //AMREX_ASSERT( (z-d_hub_height) <= rad );
                        // Conside the vector that joines the point and the turbine center.
                        // Dot it on to the vector that joins the turbine center and along
                        // the plane of the disk. See fig. 10 in Mirocha et. al. 2014.

                        Real vec_proj = (x-d_xloc_ptr[it])*(std::sin(phi)) +
                                        (y-d_yloc_ptr[it])*(-std::cos(phi));


                        Real zeta = std::atan2(z-d_hub_height, vec_proj);
                        //printf("zeta val is %0.15g\n", zeta*180.0/PI);
                        Fn_and_Ft = compute_source_terms_Fn_Ft(rad, avg_vel,
                                                               bld_rad_loc_ptr,
                                                               bld_twist_ptr,
                                                               bld_chord_ptr,
                                                               n_bld_sections,
                                                               d_bld_airfoil_aoa_ptr[index],
                                                               d_bld_airfoil_Cl_ptr[index],
                                                               d_bld_airfoil_Cd_ptr[index],
                                                               n_pts_airfoil);

                        Real Fn = Fn_and_Ft[0];
                        Real Ft = Fn_and_Ft[1];
                        // Compute the source terms - pass in radial distance, free stream velocity

                        Real Fx = Fn*std::cos(phi) + Ft*std::sin(zeta)*std::sin(phi);
                        Real Fy = Fn*std::sin(phi) - Ft*std::sin(zeta)*std::cos(phi);
                        Real Fz = -Ft*std::cos(zeta);

                        source_x = -Fx*inv_dens_vol;
                        source_y = -Fy*inv_dens_vol;
                        source_z = -Fz*inv_dens_vol;


                        //printf("Val source_x, is %0.15g, %0.15g, %0.15g %0.15g %0.15g %0.15g\n", rad, Fn, Ft, source_x, source_y, source_z);
                    }
                }
            }

            if(check_int > 1){
                amrex::Error("Actuator disks are overlapping. Visualize actuator_disks.vtk "
                             "and check the windturbine locations input file. Exiting..");
            }

            generalAD_array(i,j,k,0) = source_x;
            generalAD_array(i,j,k,1) = source_y;
            generalAD_array(i,j,k,2) = source_z;
         });
    }
}
