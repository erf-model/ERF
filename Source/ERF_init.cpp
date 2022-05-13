/**
 * \file ERF_init.cpp
 */

#include <ERF.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF
void
ERF::init_from_wrfinput(const amrex::Box& bx, FArrayBox& state_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab, FArrayBox& z_vel_fab
#ifdef ERF_USE_TERRAIN
                       ,FArrayBox& z_phys
#endif // ERF_USE_TERRAIN
                                       )
{
    //
    // FArrayBox to FArrayBox copy does "copy on intersection"
    // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
    //
    // This copies x-vel
    x_vel_fab.template copy<RunOn::Device>(NC_xvel_fab);

    // This copies y-vel
    y_vel_fab.template copy<RunOn::Device>(NC_yvel_fab);

    // This copies z-vel
    z_vel_fab.template copy<RunOn::Device>(NC_xvel_fab);

    // We first initialize all state_fab variables to zero
    state_fab.template setVal<RunOn::Device>(0.);

    // This copies the density
    state_fab.template copy<RunOn::Device>(NC_rho_fab, 0, Rho_comp, 1);

    // This copies (rho*theta)
    state_fab.template copy<RunOn::Device>(NC_rhotheta_fab, 0, RhoTheta_comp, 1);

#ifdef ERF_USE_TERRAIN
    // This copies from NC_zphys on z-faces to z_phys_nd on nodes
    Array4<Real>    z_arr   = z_phys.array();
    Array4<Real> nc_phb_arr = NC_PHB_fab.array();
    Array4<Real> nc_ph_arr  = NC_PH_fab.array();

    const Box& z_phys_box(z_phys.box());

    Box nodal_dom = amrex::surroundingNodes(geom[0].Domain());

    int ilo = nodal_dom.smallEnd()[0];
    int ihi = nodal_dom.bigEnd()[0];
    int jlo = nodal_dom.smallEnd()[1];
    int jhi = nodal_dom.bigEnd()[1];
    int klo = nodal_dom.smallEnd()[2];
    int khi = nodal_dom.bigEnd()[2];

    //
    // We must be careful not to read out of bounds of the WPS data
    //
    amrex::ParallelFor(z_phys_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int ii = std::max(std::min(i,ihi-1),ilo+1);
        int jj = std::max(std::min(j,jhi-1),jlo+1);
        if (!z_phys_box.contains(ii,jj,k) || !z_phys_box.contains(ii-1,jj-1,k)) amrex::Print() << "BAD " << IntVect(ii,jj,k) << std::endl;
        if (k < 0) {
            Real z_klo   =  0.25 * ( nc_ph_arr (ii,jj,klo  ) +  nc_ph_arr(ii-1,jj,klo  ) + nc_ph_arr (ii,jj-1,klo  ) + nc_ph_arr (ii-1,jj-1,klo) +
                                     nc_phb_arr(ii,jj,klo  ) + nc_phb_arr(ii-1,jj,klo  ) + nc_phb_arr(ii,jj-1,klo  ) + nc_phb_arr(ii-1,jj-1,klo) ) / CONST_GRAV;
            Real z_klop1 =  0.25 * ( nc_ph_arr (ii,jj,klo+1) +  nc_ph_arr(ii-1,jj,klo+1) + nc_ph_arr (ii,jj-1,klo+1) + nc_ph_arr (ii-1,jj-1,klo+1) +
                                     nc_phb_arr(ii,jj,klo+1) + nc_phb_arr(ii-1,jj,klo+1) + nc_phb_arr(ii,jj-1,klo+1) + nc_phb_arr(ii-1,jj-1,klo+1) ) / CONST_GRAV;
            z_arr(i, j, k) = 2.0 * z_klo - z_klop1;
        } else if (k > khi) {
            Real z_khi   =  0.25 * ( nc_ph_arr (ii,jj,khi  ) + nc_ph_arr (ii-1,jj,khi  ) + nc_ph_arr (ii,jj-1,khi  ) + nc_ph_arr (ii-1,jj-1,khi) +
                                     nc_phb_arr(ii,jj,khi  ) + nc_phb_arr(ii-1,jj,khi  ) + nc_phb_arr(ii,jj-1,khi  ) + nc_phb_arr(ii-1,jj-1,khi) ) / CONST_GRAV;
            Real z_khim1 =  0.25 * ( nc_ph_arr (ii,jj,khi-1) + nc_ph_arr (ii-1,jj,khi-1) + nc_ph_arr (ii,jj-1,khi-1) + nc_ph_arr (ii-1,jj-1,khi-1) +
                                     nc_phb_arr(ii,jj,khi-1) + nc_phb_arr(ii-1,jj,khi-1) + nc_phb_arr(ii,jj-1,khi-1) + nc_phb_arr(ii-1,jj-1,khi-1) ) / CONST_GRAV;
            z_arr(i, j, k) = 2.0 * z_khi - z_khim1;
          } else {
            z_arr(i, j, k) = 0.25 * ( nc_ph_arr (ii,jj,k) +  nc_ph_arr(ii-1,jj,k) +  nc_ph_arr(ii,jj-1,k) +  nc_ph_arr(ii-1,jj-1,k) +
                                      nc_phb_arr(ii,jj,k) + nc_phb_arr(ii-1,jj,k) + nc_phb_arr(ii,jj-1,k) + nc_phb_arr(ii-1,jj-1,k) ) / CONST_GRAV;
        } // k
    });
#endif
}
#endif // ERF_USE_NETCDF

void ERF::init_from_input_sounding(
        const amrex::Box &bx,
        amrex::Array4<amrex::Real> const &state,
        amrex::Array4<amrex::Real> const &x_vel,
        amrex::Array4<amrex::Real> const &y_vel,
        amrex::Array4<amrex::Real> const &z_vel,
        amrex::GeometryData const &geomdata,
        InputSoundingData const &inputSoundingData) {

    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d.dataPtr();
    const Real* U_inp_sound     = inputSoundingData.U_inp_sound_d.dataPtr();
    const Real* V_inp_sound     = inputSoundingData.V_inp_sound_d.dataPtr();
    const int   inp_sound_size  = inputSoundingData.size();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Geometry
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // TODO: Read this from file, the way we do for custom problems
        // Or provide rho = rho (z) as applicable or computer rho = rho(z) as WRF does
        Real rho_0 = 1.0;

        // Set the density
        state(i, j, k, Rho_comp) = rho_0;

        // Initial Rho0*Theta0
        state(i, j, k, RhoTheta_comp) = rho_0 * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state(i, j, k, RhoScalar_comp) = 0;
    });

    // Construct a box that is on x-faces
    const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
    // Construct a box that is on y-faces
    const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
    // Construct a box that is on z-faces
    const amrex::Box& zbx = amrex::surroundingNodes(bx,2);

    // Set the x,y,z-velocities
    amrex::ParallelFor(xbx, ybx, zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of x-faces
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel(i, j, k) = interpolate_1d(z_inp_sound, U_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of y-faces
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel(i, j, k) = interpolate_1d(z_inp_sound, V_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of z-faces
        // Set the z-velocity
        z_vel(i, j, k) = 0.0;
    });
}
