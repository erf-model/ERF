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

    Real met_lateral_sponge_strength = -1.0, met_lateral_sponge_length = -1.0;
    Real met_zhi_sponge_strength = -1.0, met_zhi_sponge_length = -1.0;
    bool use_met_zhi_sponge_damping = false;

    amrex::ParmParse pp("erf");

    pp.query("met_lateral_sponge_strength", met_lateral_sponge_strength);
    pp.query("met_lateral_sponge_length", met_lateral_sponge_length);

    pp.query("met_zhi_sponge_length", met_zhi_sponge_length);
    pp.query("met_zhi_sponge_strength", met_zhi_sponge_strength);

    pp.query("use_met_zhi_sponge_damping", use_met_zhi_sponge_damping);

    if (met_lateral_sponge_strength < 0.0) {
        amrex::Abort("ERROR: Missing input parameter 'erf.met_lateral_sponge_strength' or it is specified to be less than zero");
    }

    if (met_lateral_sponge_length < 0.0) {
        amrex::Abort("ERROR: Missing input parameter 'erf.met_lateral_sponge_length' or it is specified to be less than zero");
    }

    if (met_zhi_sponge_strength < 0.0) {
        amrex::Abort("ERROR: Missing input parameter 'erf.met_zhi_sponge_strength' or it is specified to be less than zero");
    }

    if (met_zhi_sponge_strength < 0.0) {
        amrex::Abort("ERROR: Missing input parameter 'erf.met_zhi_sponge_strength' or it is specified to be less than zero");
    }

    // Domain valid box
    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1) + 1;

    Real xlo_sponge_end   = ProbLoArr[0] + met_lateral_sponge_length;
    Real xhi_sponge_start = ProbHiArr[0] - met_lateral_sponge_length;
    Real ylo_sponge_end   = ProbLoArr[1] + met_lateral_sponge_length;
    Real yhi_sponge_start = ProbHiArr[1] - met_lateral_sponge_length;
    Real zhi_sponge_start = ProbHiArr[2] - met_zhi_sponge_length;

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
                Real xi = (xlo_sponge_end - x) / met_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // x hi sponge
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / met_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // y lo sponge
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / met_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
            }
        // x right sponge
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / met_lateral_sponge_length;
                rho_u_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
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
                Real xi = (xlo_sponge_end - x) / met_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
        // x hi sponge
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / met_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }

        // y lo sponge
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / met_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
        // x right sponge
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / met_lateral_sponge_length;
                rho_v_rhs(i, j, k) -= met_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
            }
    });

    ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);

        Real z = z_phys_nd(i,j,k);

        if(use_met_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / met_zhi_sponge_length;
                rho_w_rhs(i, j, k) -= met_zhi_sponge_strength * xi * xi * (rho_w(i, j, k) - 0.0);
            }
        }
    });
}
