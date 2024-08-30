#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

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
                   const Geometry& /*geomdata*/,
                   Real /*time*/,
                   const int* /*bcrec*/,
                   const int /*level*/)
{
    auto const dat = datfab.array();
    auto cfab      = derfab.array();

    // NOTE: we compute the soundspeed of dry air -- we do not account for any moisture effects here
    Real qv = 0.;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real rhotheta = dat(i, j, k, RhoTheta_comp);
        const Real rho      = dat(i, j, k, Rho_comp);
        AMREX_ALWAYS_ASSERT(rhotheta > 0.);
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
        AMREX_ALWAYS_ASSERT(rhotheta > 0.);
        tfab(i,j,k) = getTgivenRandRTh(rho,rhotheta);
    });
}
void
erf_dermoisttemp (const Box& bx,
             FArrayBox& derfab,
             int /*dcomp*/,
             int /*ncomp*/,
             const FArrayBox& datfab,
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
        AMREX_ALWAYS_ASSERT(rhotheta > 0.);
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
           const Geometry& /*geomdata*/,
           Real /*time*/,
           const int* /*bcrec*/,
           const int /*level*/)
{
    erf_derrhodivide(bx, derfab, datfab, RhoKE_comp);
}

/**
 * Function to define QKE by dividing (rho QKE) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds QKE
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derQKE (const Box& bx,
            FArrayBox& derfab,
            int /*dcomp*/,
            int /*ncomp*/,
            const FArrayBox& datfab,
            const Geometry& /*geomdata*/,
            Real /*time*/,
            const int* /*bcrec*/,
            const int /*level*/)
{
    erf_derrhodivide(bx, derfab, datfab, RhoQKE_comp);
}

void
erf_dervortx (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity x-component

    const Real dy = geomdata.CellSize(1);
    const Real dz = geomdata.CellSize(2);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        tfab(i,j,k,dcomp) = (dat(i,j+1,k,2) - dat(i,j-1,k,2)) / (2.0*dy)  // dw/dy
                          - (dat(i,j,k+1,1) - dat(i,j,k-1,1)) / (2.0*dz); // dv/dz
    });
}

void
erf_dervorty (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array(); // cell-centered velocity
    auto tfab      = derfab.array(); // cell-centered vorticity y-component

    const Real dx = geomdata.CellSize(0);
    const Real dz = geomdata.CellSize(2);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        tfab(i,j,k,dcomp) = (dat(i,j,k+1,0) - dat(i,j,k-1,0)) / (2.0*dz)  // du/dz
                          - (dat(i+1,j,k,2) - dat(i-1,j,k,2)) / (2.0*dx); // dw/dx
    });
}

void
erf_dervortz (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& geomdata,
  amrex::Real /*time*/,
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
        tfab(i,j,k,dcomp) = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / (2.0*dx)  // dv/dx
                          - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / (2.0*dy); // du/dy
    });
}

void
erf_dermagvel (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
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

} // namespace
