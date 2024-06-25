#include "prob.H"
#include "EOS.H"
#include "TerrainMetrics.H"

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
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("U_0", parms.U_0);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& z_nd,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    // Geometry (note we must include these here to get the data on device)
    ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const auto prob_lo  = geomdata.ProbLo();
        const auto dx       = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = z_cc(i,j,k);

        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    // Set the x-velocity
    ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = 0.0;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = 0.0;
    });

    // Set the z-velocity from impenetrable condition
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = 0.0;
    });

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_terrain(
    const Geometry& geom,
    MultiFab& z_phys_nd,
    const Real& /*time*/)
{
    // Check if a valid csv file exists for the terrain
    std::string fname;
    amrex::ParmParse pp("erf");
    auto valid_fname = pp.query("terrain_file_name",fname);
    if (valid_fname) {
        this->read_custom_terrain(fname,geom,z_phys_nd,time);
    } else {

        // Domain cell size and real bounds
        auto dx = geom.CellSizeArray();
        auto ProbLoArr = geom.ProbLoArray();
        auto ProbHiArr = geom.ProbHiArray();

        // Domain valid box (z_nd is nodal)
        const amrex::Box& domain = geom.Domain();
        int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
        // int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
        int domlo_z = domain.smallEnd(2);

        // User function parameters
        Real a    = 0.5;
        Real num  = 8 * a * a * a;
        Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]);
        Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

        // Number of ghost cells
        int ngrow = z_phys_nd.nGrow();

        // Populate bottom plane
        int k0 = domlo_z;

        for ( amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // Grown box with no z range
            amrex::Box xybx = mfi.growntilebox(ngrow);
            xybx.setRange(2,0);

            amrex::Array4<Real> const& z_arr = z_phys_nd.array(mfi);

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {

                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                // Location of nodes
                Real x = (ii  * dx[0] - xcen);

                if(fabs(x)<a){
                    z_arr(i,j,k0) = pow(a*a - x*x, 0.5);
                }
                else{
                    z_arr(i,j,k0) = 0.0;
                }

            });
        }
    }
}
