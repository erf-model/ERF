#include "ERF_Prob.H"
#include "ERF_TerrainMetrics.H"

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
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);

  pp.query("hmax", parms.hmax);
  pp.query("L", parms.L);
  pp.query("z_offset", parms.z_offset);

  pp.query("dir", parms.dir);

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
    Array4<Real const> const& mf_u,
    Array4<Real const> const& mf_v,
    const SolverChoice& sc,
    const int /*lev*/)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

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
    ParallelFor(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        y_vel_pert(i, j, k) = parms_d.V_0;
    });

    const auto dx = geomdata.CellSize();
    amrex::GpuArray<Real, AMREX_SPACEDIM> dxInv;
    dxInv[0] = 1. / dx[0];
    dxInv[1] = 1. / dx[1];
    dxInv[2] = 1. / dx[2];

    // Set the z-velocity from impenetrable condition
    if (sc.terrain_type == TerrainType::StaticFittedMesh) {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = WFromOmega(i, j, k, 0.0,
                                             x_vel_pert, y_vel_pert,
                                             mf_u, mf_v, z_nd, dxInv);
        });
    } else {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = 0.0;
        });
    }

    amrex::Gpu::streamSynchronize();
}

void
Problem::init_custom_terrain (
    const Geometry& geom,
    FArrayBox& terrain_fab,
    const Real& /*time*/)
{
    //
    // We put this here as a convenience for testing the map factor implementation
    // Note that these factors must match those in Source/ERF_MakeNewArrays.cpp
    //
    ParmParse pp("erf");
    bool test_mapfactor = false;
    pp.query("test_mapfactor",test_mapfactor);

    Real mf_m;
    if (test_mapfactor) {
        mf_m = 0.5;
    } else {
        mf_m = 1.;
    }

    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbLoArr = geom.ProbLoArray();
    auto ProbHiArr = geom.ProbHiArray();

    const amrex::Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2);

    // User function parameters
    Real a    = 0.5;
    Real num  = 8. * a * a * a;
    Real xcen = 0.5 * (ProbLoArr[0] + ProbHiArr[0]) / mf_m;
    Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]) / mf_m;

    // if hm is nonzero, then use alternate hill definition
    Real hm = parms.hmax;
    Real L = parms.L;
    Real z_offset = parms.z_offset;

    // Populate bottom plane
    int k0 = domlo_z;

    amrex::Box zbx = terrain_fab.box();
    if (zbx.smallEnd(2) <= k0)
    {
        amrex::Array4<Real> const& z_arr = terrain_fab.array();

#if 0
        Real mf_y;
        if (test_mapfactor) {
            mf_y = 0.25;
        } else {
            mf_y = 1.;
        }

        // This is a 3D hill with variation in both the x- and y-directions
        int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1) + 1;
        Real ycen = 0.5 * (ProbLoArr[1] + ProbHiArr[1]) / mf_y;
        num  = 8000. * a * a * a;

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Clip indices for ghost-cells
            int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
            int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

            // Location of nodes
            Real x = (ProbLoArr[0] + ii * dx[0]) / mf_m - xcen;
            Real y = (ProbLoArr[1] + jj * dx[1]) / mf_m - ycen;

            // WoA Hill in x-direction
            z_arr(i,j,k0) = num / (x*x + y*y + 4 * a * a);
        });
#else
        // This is a 2D hill with variation in only the x-direction
        if (parms.dir == 0) {
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen) * mf_m;

                // WoA Hill in x-direction
                if (hm==0) {
                    z_arr(i,j,k0) = num / (x*x + 4 * a * a);
                } else {
                    Real x_L = x / L;
                    z_arr(i,j,k0) = hm / (1 + x_L*x_L) + z_offset;
                }
            });
        } else if (parms.dir == 1) {
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                // Location of nodes
                Real y = (ProbLoArr[1] + jj * dx[1] - ycen) * mf_m;

                // WoA Hill in x-direction
                if (hm==0) {
                    z_arr(i,j,k0) = num / (y*y + 4.0 * a * a);
                } else {
                    Real y_L = y / L;
                    z_arr(i,j,k0) = hm / (1.0 + y_L*y_L) + z_offset;
                }
            });
        } else if (parms.dir == 2) {
            ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                // Clip indices for ghost-cells
                int ii = amrex::min(amrex::max(i,domlo_x),domhi_x);
                int jj = amrex::min(amrex::max(j,domlo_y),domhi_y);

                // Location of nodes
                Real x = (ProbLoArr[0] + ii * dx[0] - xcen) * mf_m;
                Real y = (ProbLoArr[1] + jj * dx[1] - ycen) * mf_m;
                Real r = std::sqrt(x*x + y*y);

                // WoA Hill in x-direction
                if (hm==0) {
                    z_arr(i,j,k0) = num / (r*r + 4.0 * a * a);
                } else {
                    Real r_L = r / L;
                    z_arr(i,j,k0) = hm / (1.0 + r_L*r_L) + z_offset;
                }
            });
        } else {
            amrex::Abort("Unknown dir in ERF_Prob.cpp");
        }
#endif
    }
}
