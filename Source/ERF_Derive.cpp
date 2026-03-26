#include "ERF_Derive.H"
#include "ERF_EOS.H"
#include "ERF_StormDiagnostics.H"
#include "ERF_IndexDefines.H"

using namespace amrex;

namespace derived {

/**
 * Function to define a derived quantity by dividing by density
 * (analogous to cons_to_prim)
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity
 * @params[in] datfab array of data used to construct derived quantity
 * @params[in] scalar_index index of quantity to be divided by density
*/
void erf_derrhodivide (const Box& bx,
                       FArrayBox& derfab,
                       const FArrayBox& datfab,
                       const int scalar_index)
{
    // This routine divides any cell-centered conserved quantity by density
    auto const dat = datfab.array();
    auto primitive  = derfab.array();

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real rho       = dat(i, j, k, Rho_comp);
        const Real conserved = dat(i, j, k, scalar_index);
        primitive(i,j,k) = conserved / rho;
    });
}

/**
 * Placeholder function that does nothing
*/
void
erf_dernull (const Box& /*bx*/,
             FArrayBox& /*derfab*/,
             int /*dcomp*/,
             int /*ncomp*/,
             const FArrayBox& /*datfab*/,
             const FArrayBox& /*zcc_fab*/,
             const Geometry& /*geomdata*/,
             Real /*time*/,
             const int* /*bcrec*/,
             const int /*level*/)
{ }

/**
 * Function to define the sound speed by calling an EOS routine
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds pressure
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_dersoundspeed (const Box& bx,
                   FArrayBox& derfab,
                   int /*dcomp*/,
                   int /*ncomp*/,
                   const FArrayBox& datfab,
                   const FArrayBox& /*zcc_fab*/,
                   const Geometry& /*geomdata*/,
                   Real /*time*/,
                   const int* /*bcrec*/,
                   const int /*level*/)
{
    auto const dat = datfab.array();
    auto cfab      = derfab.array();

    // NOTE: we compute the soundspeed of dry air -- we do not account for any moisture effects here
    Real qv = Real(0.0);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real rhotheta = dat(i, j, k, RhoTheta_comp);
        const Real rho      = dat(i, j, k, Rho_comp);
        AMREX_ALWAYS_ASSERT(rhotheta > 0);
        cfab(i,j,k) = std::sqrt(Gamma * getPgivenRTh(rhotheta,qv) / rho);
    });
}

/**
 * Function to define the temperature by calling an EOS routine
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds pressure
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_dertemp (const Box& bx,
             FArrayBox& derfab,
             int /*dcomp*/,
             int /*ncomp*/,
             const FArrayBox& datfab,
             const FArrayBox& /*zcc_fab*/,
             const Geometry& /*geomdata*/,
             Real /*time*/,
             const int* /*bcrec*/,
             const int /*level*/)
{
    auto const dat = datfab.array();
    auto tfab      = derfab.array();

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real rho = dat(i, j, k, Rho_comp);
        const Real rhotheta = dat(i, j, k, RhoTheta_comp);
        AMREX_ALWAYS_ASSERT(rhotheta > Real(0.0));
        tfab(i,j,k) = getTgivenRandRTh(rho,rhotheta);
    });
}
void
erf_dermoisttemp (const Box& bx,
             FArrayBox& derfab,
             int /*dcomp*/,
             int /*ncomp*/,
             const FArrayBox& datfab,
             const FArrayBox& /*zcc_fab*/,
             const Geometry& /*geomdata*/,
             Real /*time*/,
             const int* /*bcrec*/,
             const int /*level*/)
{
    auto const dat = datfab.array();
    auto tfab      = derfab.array();

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real rho = dat(i, j, k, Rho_comp);
        const Real rhotheta = dat(i, j, k, RhoTheta_comp);
        AMREX_ALWAYS_ASSERT(rhotheta > Real(0.0));
        const Real qv = dat(i, j, k, RhoQ1_comp) / rho;
        tfab(i,j,k) = getTgivenRandRTh(rho,rhotheta,qv);
    });
}

/**
 * Function to define the potential temperature by calling an EOS routine
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds pressure
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_dertheta (const Box& bx,
              FArrayBox& derfab,
              int /*dcomp*/,
              int /*ncomp*/,
              const FArrayBox& datfab,
              const FArrayBox& /*zcc_fab*/,
              const Geometry& /*geomdata*/,
              Real /*time*/,
              const int* /*bcrec*/,
              const int /*level*/)
{
    erf_derrhodivide(bx, derfab, datfab, RhoTheta_comp);
}

/**
 * Function to define a scalar s by dividing (rho s) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds scalar s
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derscalar (const Box& bx,
               FArrayBox& derfab,
               int /*dcomp*/,
               int /*ncomp*/,
               const FArrayBox& datfab,
               const FArrayBox& /*zcc_fab*/,
               const Geometry& /*geomdata*/,
               Real /*time*/,
               const int* /*bcrec*/,
               const int /*level*/)
{
    erf_derrhodivide(bx, derfab, datfab, RhoScalar_comp);
}

/**
 * Function to define the kinetic energy KE by dividing (rho KE) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds KE
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derKE (const Box& bx,
           FArrayBox& derfab,
           int /*dcomp*/,
           int /*ncomp*/,
           const FArrayBox& datfab,
           const FArrayBox& /*zcc_fab*/,
           const Geometry& /*geomdata*/,
           Real /*time*/,
           const int* /*bcrec*/,
           const int /*level*/)
{
    erf_derrhodivide(bx, derfab, datfab, RhoKE_comp);
}

void
erf_dervortx ( const Box& bx,
               FArrayBox& derfab,
               int dcomp,
               int ncomp,
               const FArrayBox& datfab,
               const FArrayBox& zcc_fab,
               const Geometry& geomdata,
               Real /*time*/,
               const int* /*bcrec*/,
               const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity x-component
    auto z_arr     = zcc_fab.array(); // cell-centered height z

    const Real two_dy = two * geomdata.CellSize(1);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real two_dz = z_arr(i,j,k+1) - z_arr(i,j,k-1);
        tfab(i,j,k,dcomp) = (dat(i,j+1,k,2) - dat(i,j-1,k,2)) / two_dy  // dw/dy
                          - (dat(i,j,k+1,1) - dat(i,j,k-1,1)) / two_dz; // dv/dz
    });
}

void
erf_dervorty ( const Box& bx,
               FArrayBox& derfab,
               int dcomp,
               int ncomp,
               const FArrayBox& datfab,
               const FArrayBox& zcc_fab,
               const Geometry& geomdata,
               Real /*time*/,
               const int* /*bcrec*/,
               const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity y-component
    auto z_arr     = zcc_fab.array(); // cell-centered height z

    const Real two_dx = two * geomdata.CellSize(0);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real two_dz = z_arr(i,j,k+1) - z_arr(i,j,k-1);
        tfab(i,j,k,dcomp) = (dat(i,j,k+1,0) - dat(i,j,k-1,0)) / two_dz  // du/dz
                          - (dat(i+1,j,k,2) - dat(i-1,j,k,2)) / two_dx; // dw/dx
    });
}

void
erf_dervortz ( const Box& bx,
               FArrayBox& derfab,
               int dcomp,
               int ncomp,
               const FArrayBox& datfab,
               const FArrayBox& /*zcc_fab*/,
               const Geometry& geomdata,
               Real /*time*/,
               const int* /*bcrec*/,
               const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity z-component

    const Real dx = geomdata.CellSize(0);
    const Real dy = geomdata.CellSize(2);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        tfab(i,j,k,dcomp) = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / (two*dx)  // dv/dx
                          - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / (two*dy); // du/dy
    });
}

void
erf_derenstrophysq ( const Box& bx,
                     FArrayBox& derfab,
                     int dcomp,
                     int ncomp,
                     const FArrayBox& datfab,
                     const FArrayBox& zcc_fab,
                     const Geometry& geomdata,
                     Real /*time*/,
                     const int* /*bcrec*/,
                     const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity x-component
    auto z_arr     = zcc_fab.array(); // cell-centered height z

    const Real two_dx = two * geomdata.CellSize(0);
    const Real two_dy = two * geomdata.CellSize(1);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real two_dz = z_arr(i,j,k+1) - z_arr(i,j,k-1);

        Real vortx = (dat(i,j+1,k,2) - dat(i,j-1,k,2)) / two_dy  // dw/dy
                    -(dat(i,j,k+1,1) - dat(i,j,k-1,1)) / two_dz; // dv/dz
        Real vorty = (dat(i,j,k+1,0) - dat(i,j,k-1,0)) / two_dz  // du/dz
                    -(dat(i+1,j,k,2) - dat(i-1,j,k,2)) / two_dx; // dw/dx
        Real vortz = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / two_dx  // dv/dx
                    -(dat(i,j+1,k,0) - dat(i,j-1,k,0)) / two_dy; // du/dy

        tfab(i,j,k,dcomp) = vortx*vortx + vorty*vorty + vortz*vortz;
    });
}

void
erf_dermagvel ( const Box& bx,
                FArrayBox& derfab,
                int dcomp,
                int ncomp,
                const FArrayBox& datfab,
                const FArrayBox& /*zcc_fab*/,
                const Geometry& /*geomdata*/,
                Real /*time*/,
                const int* /*bcrec*/,
                const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered magvel

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real u = dat(i,j,k,0);
        Real v = dat(i,j,k,1);
        Real w = dat(i,j,k,2);
        tfab(i,j,k,dcomp) = std::sqrt(u*u + v*v + w*w);
    });
}

void
erf_dermagvelsq ( const Box& bx,
                  FArrayBox& derfab,
                  int dcomp,
                  int ncomp,
                  const FArrayBox& datfab,
                  const FArrayBox& /*zcc_fab*/,
                  const Geometry& /*geomdata*/,
                  Real /*time*/,
                  const int* /*bcrec*/,
                  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered magvel

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real u = dat(i,j,k,0);
        Real v = dat(i,j,k,1);
        Real w = dat(i,j,k,2);
        tfab(i,j,k,dcomp) = u*u + v*v + w*w;
    });
}

void
erf_derreflectivity ( const Box& bx,
                      FArrayBox& derfab,
                      int dcomp,
                      int /*ncomp*/,
                      const FArrayBox& datfab,
                      const FArrayBox& /*zcc_fab*/,
                      const Geometry& /*geomdata*/,
                      Real /*time*/,
                      const int* /*bcrec*/,
                      const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);

    auto const dat = datfab.array(); // cell-centered state vector
    auto rfab      = derfab.array(); // cell-centered reflectivity

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real rho = dat(i,j,k,Rho_comp);
        Real qv  = std::max(Real(0.0),dat(i,j,k,RhoQ1_comp)/rho);
        Real qpr = std::max(Real(0.0),dat(i,j,k,RhoQ4_comp)/rho);
        Real qps = std::max(Real(0.0),dat(i,j,k,RhoQ5_comp)/rho);
        Real qpg = std::max(Real(0.0),dat(i,j,k,RhoQ6_comp)/rho);

        Real temp  = getTgivenRandRTh(rho, dat(i,j,k,RhoTheta_comp), qv);

        rfab(i, j, k, dcomp) = compute_max_reflectivity_dbz(rho, temp, qpr, qps, qpg,
                                                            1, 1, 1, 1);
    });
}

void
erf_dermaxreflectivity (
  const Box& bx,
                  FArrayBox& derfab,
                  int dcomp,
                  int /*ncomp*/,
                  const FArrayBox& datfab,
                  const FArrayBox& /*zcc_fab*/,
                  const Geometry& /*geomdata*/,
                  Real /*time*/,
                  const int* /*bcrec*/,
                  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);

    auto const dat = datfab.array(); // cell-centered state vector
    auto rfab      = derfab.array(); // cell-centered max reflectivity

    // Collapse to i,j box (ignore vertical for now)
    Box b2d = bx;
    b2d.setSmall(2,0);
    b2d.setBig(2,0);

    ParallelFor(b2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {

        Real max_dbz = -1.0e30;

        // find max reflectivity over k
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {

            Real rho = dat(i,j,k,Rho_comp);
            Real qv  = std::max(Real(0.0),dat(i,j,k,RhoQ1_comp)/rho);
            Real qpr = std::max(Real(0.0),dat(i,j,k,RhoQ4_comp)/rho);
            Real qps = std::max(Real(0.0),dat(i,j,k,RhoQ5_comp)/rho);
            Real qpg = std::max(Real(0.0),dat(i,j,k,RhoQ6_comp)/rho);

            Real temp = getTgivenRandRTh(rho, dat(i,j,k,RhoTheta_comp), qv);

            Real dbz = compute_max_reflectivity_dbz(rho, temp, qpr, qps, qpg,
                                                    1, 1, 1, 1);
            max_dbz = amrex::max(max_dbz, dbz);
        }

        // Store max_dbz into *all* levels for this (i,j)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            rfab(i, j, k, dcomp) = max_dbz;
        }
    });
}

void
erf_derlocalhelicity (
  const Box& bx,
                  FArrayBox& derfab,
                  int dcomp,
                  int /*ncomp*/,
                  const FArrayBox& datfab,
                  const FArrayBox& /*zcc_fab*/,
                  const Geometry& geomdata,
                  Real /*time*/,
                  const int* /*bcrec*/,
                  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);

    auto const dat = datfab.array(); // cell-centered velocity
    auto dfab      = derfab.array(); // cell-centered local helicity

    const Real two_dx = Real(2.0)*geomdata.CellSize(0);
    const Real two_dy = Real(2.0)*geomdata.CellSize(1);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        Real vortz = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / two_dx  // dv/dx
                   - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / two_dy; // du/dy
        Real w     = dat(i,j,k,2);

        // Helicity
        dfab(i,j,k,dcomp) = vortz * w;
    });
}

void
erf_derhelicity ( const Box& bx,
                  FArrayBox& derfab,
                  int dcomp,
                  int /*ncomp*/,
                  const FArrayBox& datfab,
                  const FArrayBox& zcc_fab,
                  const Geometry& geomdata,
                  Real /*time*/,
                  const int* /*bcrec*/,
                  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);

    auto const dat = datfab.array(); // cell-centered velocity
    auto dfab      = derfab.array(); // cell-centered local helicity
    auto z_arr     = zcc_fab.array(); // cell-centered height z

    const Real dx = geomdata.CellSize(0);
    const Real dy = geomdata.CellSize(1);

    // Collapse to i,j box (ignore vertical for now)
    Box b2d = bx;
    b2d.setSmall(2,0);
    b2d.setBig(2,0);

    ParallelFor(b2d, [=] AMREX_GPU_DEVICE(int i, int j, int ) noexcept
    {
        Real int_hel = Real(0.0);
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
        {
            Real z = z_arr(i,j,k);

            // Helicity is defined as integral from 2km to 5km in vertical
            if (z > Real(2000.0) && z < Real(5000.0)) {

                Real z_hi = Real(0.5) * (z_arr(i,j,k) + z_arr(i,j,k+1));
                Real z_lo = Real(0.5) * (z_arr(i,j,k) + z_arr(i,j,k-1));
                Real dz = z_hi - z_lo;

                Real vortz = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / (2.0*dx)  // dv/dx
                           - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / (2.0*dy); // du/dy
                Real w     = dat(i,j,k,2); // vertical velocity

                int_hel += vortz * w * dz;
            }
        }

        // Store vertical integral into *all* levels for this (i,j)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            dfab(i, j, k, dcomp) = int_hel;
        }
    });
}

} // namespace
