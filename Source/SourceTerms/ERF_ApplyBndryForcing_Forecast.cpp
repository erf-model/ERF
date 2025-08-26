#include <AMReX_MultiFab.H>
#include <ERF_SrcHeaders.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
ApplyBndryForcing_Forecast (
  const SolverChoice& solverChoice,
  const Geometry geom,
  const Box& tbx,
  const Box& tby,
  const Box& tbz,
  const Array4<const Real>& z_phys_nd,
  const Array4<Real>& rho_u_rhs,
  const Array4<Real>& rho_v_rhs,
  const Array4<Real>& rho_w_rhs,
  const Array4<const Real>& rho_u,
  const Array4<const Real>& rho_v,
  const Array4<const Real>& rho_w,
  const Array4<const Real>& rho_u_initial_state,
  const Array4<const Real>& rho_v_initial_state,
  const Array4<const Real>& rho_w_initial_state,
  const Array4<const Real>& cons_initial_state)

{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();

    amrex::ignore_unused(rho_w_initial_state);

    // Domain valid box
    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1) + 1;

    Real hindcast_lateral_sponge_length   = solverChoice.hindcast_lateral_sponge_length;
    Real hindcast_lateral_sponge_strength = solverChoice.hindcast_lateral_sponge_strength;

    bool hindcast_zhi_sponge_damping      = solverChoice.hindcast_zhi_sponge_damping;
    Real hindcast_zhi_sponge_length       = solverChoice.hindcast_zhi_sponge_length;
    Real hindcast_zhi_sponge_strength     = solverChoice.hindcast_zhi_sponge_strength;

    Real xlo_sponge_end   = ProbLoArr[0] + hindcast_lateral_sponge_length;
    Real xhi_sponge_start = ProbHiArr[0] - hindcast_lateral_sponge_length;
    Real ylo_sponge_end   = ProbLoArr[1] + hindcast_lateral_sponge_length;
    Real yhi_sponge_start = ProbHiArr[1] - hindcast_lateral_sponge_length;
    Real zhi_sponge_start = ProbHiArr[2] - hindcast_zhi_sponge_length;

    AMREX_ALWAYS_ASSERT(xlo_sponge_end   > ProbLoArr[0]);
    AMREX_ALWAYS_ASSERT(xhi_sponge_start < ProbHiArr[0]);
    AMREX_ALWAYS_ASSERT(ylo_sponge_end   > ProbLoArr[1]);
    AMREX_ALWAYS_ASSERT(yhi_sponge_start < ProbHiArr[1]);

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);

        Real x = ProbLoArr[0] + ii * dx[0];
        Real y = ProbLoArr[1] + (jj+0.5) * dx[1];

        Real rho_u_sponge = rho_u_initial_state(i,j,k)*cons_initial_state(i,j,k,0);
        // x lo sponge
            if (x < xlo_sponge_end) {
                Real xi = (xlo_sponge_end - x) / hindcast_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // x hi sponge
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / hindcast_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // y lo sponge
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / hindcast_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // x right sponge
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / hindcast_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
    });


    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);

        Real x = ProbLoArr[0] + (ii+0.5) * dx[0];
        Real y = ProbLoArr[1] + jj * dx[1];

        Real rho_v_sponge    = rho_v_initial_state(i,j,k)*cons_initial_state(i,j,k,0);

        // x lo sponge
            if (x < xlo_sponge_end) {
                Real xi = (xlo_sponge_end - x) / hindcast_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
        // x hi sponge
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / hindcast_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }

        // y lo sponge
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / hindcast_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
        // x right sponge
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / hindcast_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
    });

    ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        Real z = z_phys_nd(i,j,k);

        if(hindcast_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / hindcast_zhi_sponge_length;
                rho_w_rhs(i, j, k) -= hindcast_zhi_sponge_strength * xi * xi * (rho_w(i, j, k) - 0.0);
            }
        }
    });
}
