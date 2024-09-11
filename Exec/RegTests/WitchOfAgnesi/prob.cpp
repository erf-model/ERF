#include "prob.H"
#include "TerrainMetrics.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit (
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem ()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("dampcoef", parms.dampcoef);
  pp.query("zdamp", parms.zdamp);

  pp.query("hmax", parms.hmax);
  pp.query("L", parms.L);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert (
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
    Array4<Real const> const& /*z_cc*/,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Set scalar = 0 everywhere
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    // Set the x-velocity
    ParallelFor(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        x_vel_pert(i, j, k) = parms_d.U_0;
    });

    // Set the y-velocity
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = 0.0;
    });

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0];
    dxInv[1] = 1. / dx[1];
    dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0, x_vel_pert, y_vel_pert, z_nd, dxInv);
    });

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_terrain (
    const Geometry& geom,
    MultiFab& z_phys_nd,
    const Real& time)
{

    // Check if a valid text file exists for the terrain
    std::string fname;
    ParmParse pp("erf");
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
        // Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]);

        // if hm is nonzero, then use alternate hill definition
        Real hm = parms.hmax;
        Real L = parms.L;

        // Number of ghost cells
        int ngrow = z_phys_nd.nGrow();

        // Populate bottom plane
        int k0 = domlo_z;

        for ( amrex::MFIter mfi(z_phys_nd,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            amrex::Box zbx = mfi.nodaltilebox(2);
            if (zbx.smallEnd(2) > k0) continue;

            // Grown box with no z range
            amrex::Box xybx = mfi.growntilebox(ngrow);
            xybx.setRange(2,0);

            amrex::Array4<Real> const& z_arr = z_phys_nd.array(mfi);

            ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                // int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen);
                // Real y = (jj  * dx[1] - ycen);

                // WoA Hill in x-direction
                if (hm==0) {
                    z_arr(i,j,k0) = num / (x*x + 4 * a * a);
                } else {
                    Real x_L = x / L;
                    z_arr(i,j,k0) = hm / (1 + x_L*x_L);
                }
            });
        }
    }
}
