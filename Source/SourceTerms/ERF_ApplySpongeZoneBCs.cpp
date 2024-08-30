#include <AMReX_MultiFab.H>
#include <Src_headers.H>

//#include <TerrainMetrics.H>
//#include <IndexDefines.H>
//#include <PlaneAverage.H>

using namespace amrex;

void
ApplySpongeZoneBCsForCC (
  const SpongeChoice& spongeChoice,
  const Geometry geom,
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data)
{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();

    const Real sponge_strength = spongeChoice.sponge_strength;
    const int use_xlo_sponge_damping = spongeChoice.use_xlo_sponge_damping;
    const int use_xhi_sponge_damping = spongeChoice.use_xhi_sponge_damping;
    const int use_ylo_sponge_damping = spongeChoice.use_ylo_sponge_damping;
    const int use_yhi_sponge_damping = spongeChoice.use_yhi_sponge_damping;
    const int use_zlo_sponge_damping = spongeChoice.use_zlo_sponge_damping;
    const int use_zhi_sponge_damping = spongeChoice.use_zhi_sponge_damping;

    const Real xlo_sponge_end   = spongeChoice.xlo_sponge_end;
    const Real xhi_sponge_start = spongeChoice.xhi_sponge_start;
    const Real ylo_sponge_end   = spongeChoice.ylo_sponge_end;
    const Real yhi_sponge_start = spongeChoice.yhi_sponge_start;
    const Real zlo_sponge_end   = spongeChoice.zlo_sponge_end;
    const Real zhi_sponge_start = spongeChoice.zhi_sponge_start;

    const Real sponge_density = spongeChoice.sponge_density;

    // Domain valid box
    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2);
    int domhi_z = domain.bigEnd(2) + 1;

    if(use_xlo_sponge_damping)AMREX_ALWAYS_ASSERT(xlo_sponge_end   > ProbLoArr[0]);
    if(use_xhi_sponge_damping)AMREX_ALWAYS_ASSERT(xhi_sponge_start < ProbHiArr[0]);
    if(use_ylo_sponge_damping)AMREX_ALWAYS_ASSERT(ylo_sponge_end   > ProbLoArr[1]);
    if(use_yhi_sponge_damping)AMREX_ALWAYS_ASSERT(yhi_sponge_start < ProbHiArr[1]);
    if(use_zlo_sponge_damping)AMREX_ALWAYS_ASSERT(zlo_sponge_end   > ProbLoArr[2]);
    if(use_zhi_sponge_damping)AMREX_ALWAYS_ASSERT(zhi_sponge_start < ProbHiArr[2]);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
        int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

        Real x = ProbLoArr[0] + (ii+0.5) * dx[0];
        Real y = ProbLoArr[1] + (jj+0.5) * dx[1];
        Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

        // x left sponge
        if(use_xlo_sponge_damping){
            if (x < xlo_sponge_end) {
                Real xi = (xlo_sponge_end - x) / (xlo_sponge_end - ProbLoArr[0]);
                cell_rhs(i, j, k, 0) -= sponge_strength * xi * xi * (cell_data(i, j, k, 0) - sponge_density);
            }
        }
        // x right sponge
        if(use_xhi_sponge_damping){
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / (ProbHiArr[0] - xhi_sponge_start);
                cell_rhs(i, j, k, 0) -= sponge_strength * xi * xi * (cell_data(i, j, k, 0) - sponge_density);
            }
        }

        // y left sponge
        if(use_ylo_sponge_damping){
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / (ylo_sponge_end - ProbLoArr[1]);
                cell_rhs(i, j, k, 0) -= sponge_strength * xi * xi * (cell_data(i, j, k, 0) - sponge_density);
            }
        }
        // x right sponge
        if(use_yhi_sponge_damping){
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / (ProbHiArr[1] - yhi_sponge_start);
                cell_rhs(i, j, k, 0) -= sponge_strength * xi * xi * (cell_data(i, j, k, 0) - sponge_density);
            }
        }

        // x left sponge
        if(use_zlo_sponge_damping){
            if (z < zlo_sponge_end) {
                Real xi = (zlo_sponge_end - z) / (zlo_sponge_end - ProbLoArr[2]);
                cell_rhs(i, j, k, 0) -= sponge_strength * xi * xi * (cell_data(i, j, k, 0) - sponge_density);
            }
        }
        // x right sponge
        if(use_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / (ProbHiArr[2] - zhi_sponge_start);
                cell_rhs(i, j, k, 0) -= sponge_strength * xi * xi * (cell_data(i, j, k, 0) - sponge_density);
            }
        }
    });
}

void
ApplySpongeZoneBCsForMom (
  const SpongeChoice& spongeChoice,
  const Geometry geom,
  const Box& tbx,
  const Box& tby,
  const Box& tbz,
  const Array4<Real>& rho_u_rhs,
  const Array4<Real>& rho_v_rhs,
  const Array4<Real>& rho_w_rhs,
  const Array4<const Real>& rho_u,
  const Array4<const Real>& rho_v,
  const Array4<const Real>& rho_w)
{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();

    const Real sponge_strength = spongeChoice.sponge_strength;
    const int use_xlo_sponge_damping = spongeChoice.use_xlo_sponge_damping;
    const int use_xhi_sponge_damping = spongeChoice.use_xhi_sponge_damping;
    const int use_ylo_sponge_damping = spongeChoice.use_ylo_sponge_damping;
    const int use_yhi_sponge_damping = spongeChoice.use_yhi_sponge_damping;
    const int use_zlo_sponge_damping = spongeChoice.use_zlo_sponge_damping;
    const int use_zhi_sponge_damping = spongeChoice.use_zhi_sponge_damping;

    const Real xlo_sponge_end   = spongeChoice.xlo_sponge_end;
    const Real xhi_sponge_start = spongeChoice.xhi_sponge_start;
    const Real ylo_sponge_end   = spongeChoice.ylo_sponge_end;
    const Real yhi_sponge_start = spongeChoice.yhi_sponge_start;
    const Real zlo_sponge_end   = spongeChoice.zlo_sponge_end;
    const Real zhi_sponge_start = spongeChoice.zhi_sponge_start;

    const Real sponge_density = spongeChoice.sponge_density;
    const Real sponge_x_velocity = spongeChoice.sponge_x_velocity;
    const Real sponge_y_velocity = spongeChoice.sponge_y_velocity;
    const Real sponge_z_velocity = spongeChoice.sponge_z_velocity;

    // Domain valid box
    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1) + 1;
    int domlo_z = domain.smallEnd(2);
    int domhi_z = domain.bigEnd(2) + 1;

    if(use_xlo_sponge_damping)AMREX_ALWAYS_ASSERT(xlo_sponge_end   > ProbLoArr[0]);
    if(use_xhi_sponge_damping)AMREX_ALWAYS_ASSERT(xhi_sponge_start < ProbHiArr[0]);
    if(use_ylo_sponge_damping)AMREX_ALWAYS_ASSERT(ylo_sponge_end   > ProbLoArr[1]);
    if(use_yhi_sponge_damping)AMREX_ALWAYS_ASSERT(yhi_sponge_start < ProbHiArr[1]);
    if(use_zlo_sponge_damping)AMREX_ALWAYS_ASSERT(zlo_sponge_end   > ProbLoArr[2]);
    if(use_zhi_sponge_damping)AMREX_ALWAYS_ASSERT(zhi_sponge_start < ProbHiArr[2]);

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
        int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

        Real x = ProbLoArr[0] + ii * dx[0];
        Real y = ProbLoArr[1] + (jj+0.5) * dx[1];
        Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

        // x lo sponge
        if(use_xlo_sponge_damping){
            if (x < xlo_sponge_end) {
                Real xi = (xlo_sponge_end - x) / (xlo_sponge_end - ProbLoArr[0]);
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
            }
        }
        // x hi sponge
        if(use_xhi_sponge_damping){
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / (ProbHiArr[0] - xhi_sponge_start);
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
            }
        }

        // y lo sponge
        if(use_ylo_sponge_damping){
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / (ylo_sponge_end - ProbLoArr[1]);
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
            }
        }
        // x right sponge
        if(use_yhi_sponge_damping){
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / (ProbHiArr[1] - yhi_sponge_start);
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
            }
        }

        // z lo sponge
        if(use_zlo_sponge_damping){
            if (z < zlo_sponge_end) {
                Real xi = (zlo_sponge_end - z) / (zlo_sponge_end - ProbLoArr[2]);
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
            }
        }


        // z hi sponge
        if(use_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / (ProbHiArr[2] - zhi_sponge_start);
                rho_u_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_u(i, j, k) - sponge_density*sponge_x_velocity);
            }
        }
    });


    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
        int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

        Real x = ProbLoArr[0] + (ii+0.5) * dx[0];
        Real y = ProbLoArr[1] + jj * dx[1];
        Real z = ProbLoArr[2] + (kk+0.5) * dx[2];

        // x lo sponge
        if(use_xlo_sponge_damping){
            if (x < xlo_sponge_end) {
                Real xi = (xlo_sponge_end - x) / (xlo_sponge_end - ProbLoArr[0]);
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - sponge_density*sponge_y_velocity);
            }
        }
        // x hi sponge
        if(use_xhi_sponge_damping){
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / (ProbHiArr[0] - xhi_sponge_start);
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - sponge_density*sponge_y_velocity);
            }
        }

        // y lo sponge
        if(use_ylo_sponge_damping){
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / (ylo_sponge_end - ProbLoArr[1]);
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - sponge_density*sponge_y_velocity);
            }
        }
        // x right sponge
        if(use_yhi_sponge_damping){
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / (ProbHiArr[1] - yhi_sponge_start);
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - sponge_density*sponge_y_velocity);
            }
        }

        // z lo sponge
        if(use_zlo_sponge_damping){
            if (z < zlo_sponge_end) {
                Real xi = (zlo_sponge_end - z) / (zlo_sponge_end - ProbLoArr[2]);
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - sponge_density*sponge_y_velocity);
            }
        }


        // z hi sponge
        if(use_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / (ProbHiArr[2] - zhi_sponge_start);
                rho_v_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_v(i, j, k) - sponge_density*sponge_y_velocity);
            }
        }
    });


    ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
        int kk = amrex::min(amrex::max(k, domlo_z), domhi_z);

        Real x = ProbLoArr[0] + (ii+0.5) * dx[0];
        Real y = ProbLoArr[1] + (jj+0.5) * dx[1];
        Real z = ProbLoArr[2] + kk * dx[2];

        // x left sponge
        if(use_xlo_sponge_damping){
            if (x < xlo_sponge_end) {
                Real xi = (xlo_sponge_end - x) / (xlo_sponge_end - ProbLoArr[0]);
                rho_w_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_w(i, j, k) - sponge_density*sponge_z_velocity);
            }
        }
        // x right sponge
        if(use_xhi_sponge_damping){
            if (x > xhi_sponge_start) {
                Real xi = (x - xhi_sponge_start) / (ProbHiArr[0] - xhi_sponge_start);
                rho_w_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_w(i, j, k) - sponge_density*sponge_z_velocity);
            }
        }

        // y lo sponge
        if(use_ylo_sponge_damping){
            if (y < ylo_sponge_end) {
                Real xi = (ylo_sponge_end - y) / (ylo_sponge_end - ProbLoArr[1]);
                rho_w_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_w(i, j, k) - sponge_density*sponge_z_velocity);
            }
        }
        // x right sponge
        if(use_yhi_sponge_damping){
            if (y > yhi_sponge_start) {
                Real xi = (y - yhi_sponge_start) / (ProbHiArr[1] - yhi_sponge_start);
                rho_w_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_w(i, j, k) - sponge_density*sponge_z_velocity);
            }
        }

        // z lo sponge
        if(use_zlo_sponge_damping){
            if (z < zlo_sponge_end) {
                Real xi = (zlo_sponge_end - z) / (zlo_sponge_end - ProbLoArr[2]);
                rho_w_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_w(i, j, k) - sponge_density*sponge_z_velocity);
            }
        }


        // z top sponge
        if(use_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / (ProbHiArr[2] - zhi_sponge_start);
                rho_w_rhs(i, j, k) -= sponge_strength * xi * xi * (rho_w(i, j, k) - sponge_density*sponge_z_velocity);
            }
        }
    });
}






