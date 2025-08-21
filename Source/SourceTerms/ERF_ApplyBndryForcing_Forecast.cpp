#include <AMReX_MultiFab.H>
#include <ERF_SrcHeaders.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
ApplyBndryForcing_Forecast (
  const Geometry geom,
  const Box& tbx,
  const Box& tby,
  const Box& tbz,
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

    Real sponge_strength = -1.0, sponge_length = -1.0;

    amrex::ParmParse pp("erf");

    pp.query("hurricane_sponge_strength", sponge_strength);
    pp.query("hurricane_sponge_length", sponge_length);

    if (sponge_strength == -1.0 || sponge_length == -1.0) {
        amrex::Abort("ERROR: Missing required parameters 'erf.hurricane_sponge_strength' or 'erf.hurricane_sponge_length'");
    }

    // Domain valid box
    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1) + 1;

    Real xlo_sponge_end   = ProbLoArr[0] + sponge_length;
    Real xhi_sponge_start = ProbHiArr[0] - sponge_length;
    Real ylo_sponge_end   = ProbLoArr[1] + sponge_length;
    Real yhi_sponge_start = ProbHiArr[1] - sponge_length;

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
                /*printf("Inside this routine rhou............%0.15g, %0.15g, %0.15g, %0.15g\n", rho_u(i, j, k),
                        rho_u_initial_state(i,j,k), cons_initial_state(i,j,k,0), rho_u_sponge);
                printf("Inside this routine rhov............%0.15g, %0.15g\n", rho_v(i, j, k), rho_v_initial_state(i,j,k));
                printf("Inside this routine rhow............%0.15g, %0.15g\n", rho_w(i, j, k), rho_w_initial_state(i,j,k));*/
                Real xi = (xlo_sponge_end - x) / sponge_length;
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // x hi sponge
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / sponge_length;
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // y lo sponge
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / sponge_length;
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // x right sponge
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / sponge_length;
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
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
                Real xi = (xlo_sponge_end - x) / sponge_length;
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
        // x hi sponge
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / sponge_length;
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }

        // y lo sponge
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / sponge_length;
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
        // x right sponge
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / sponge_length;
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
    });
}
