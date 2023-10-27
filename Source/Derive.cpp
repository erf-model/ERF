#include "Derive.H"
#include "EOS.H"
#include "IndexDefines.H"

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
void erf_derrhodivide (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  const amrex::FArrayBox& datfab,
  const int scalar_index)
{
  // This routine divides any cell-centered conserved quantity by density
  auto const dat = datfab.array();
  auto primitive  = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho       = dat(i, j, k, Rho_comp);
    const amrex::Real conserved = dat(i, j, k, scalar_index);
    primitive(i,j,k) = conserved / rho;
  });
}

/**
 * Placeholder function that does nothing
*/
void
erf_dernull (
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*datfab*/,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{ }

/**
 * Function to define pressure by calling an EOS routine
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds pressure
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derpres (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto pfab      = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    AMREX_ALWAYS_ASSERT(rhotheta > 0.);
#if defined(ERF_USE_MOISTURE)
    // HACK HACK HACK -- FOR NOW THIS IS ONLY THE PARTIAL PRESSURE OF THE DRY AIR
    amrex::Real qv = 0.;
    pfab(i,j,k) = getPgivenRTh(rhotheta,qv);
#elif defined(ERF_USE_WARM_NO_PRECIP)
    amrex::Real qv = dat(i,j,k,RhoQv_comp) / dat(i,j,k,Rho_comp);
    pfab(i,j,k) = getPgivenRTh(rhotheta,qv);
#else
    pfab(i,j,k) = getPgivenRTh(rhotheta);
#endif
  });
}

/**
 * Function to define the sound speed by calling an EOS routine
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds pressure
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_dersoundspeed (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto cfab      = derfab.array();

  // NOTE: we compute the soundspeed of dry air -- we do not account for any moisture effects here
  amrex::Real qv = 0.;

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    const amrex::Real rho      = dat(i, j, k, Rho_comp);
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
erf_dertemp (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  auto const dat = datfab.array();
  auto tfab      = derfab.array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    const amrex::Real rho = dat(i, j, k, Rho_comp);
    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
    AMREX_ALWAYS_ASSERT(rhotheta > 0.);
    tfab(i,j,k) = getTgivenRandRTh(rho,rhotheta);
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
erf_dertheta (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
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
erf_derscalar (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
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
erf_derKE (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
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
erf_derQKE (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  erf_derrhodivide(bx, derfab, datfab, RhoQKE_comp);
}

#if defined(ERF_USE_MOISTURE)
/**
 * Function to define total water Qt by dividing (rho Qt) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds Qt
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derQt (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  erf_derrhodivide(bx, derfab, datfab, RhoQt_comp);
}

/**
 * Function to define precipitating water Qp by dividing (rho Qp) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds Qp
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derQp (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  erf_derrhodivide(bx, derfab, datfab, RhoQp_comp);
}
#elif defined(ERF_USE_WARM_NO_PRECIP)
/**
 * Function to define water vapor Qv by dividing (rho Qv) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds Qv
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derQv (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  erf_derrhodivide(bx, derfab, datfab, RhoQv_comp);
}
/**
 * Function to define cloud water Qc by dividing (rho Qc) by rho
 *
 * @params[in] bx box on which to divide by density
 * @params[out] derfab array of derived quantity -- here it holds Qc
 * @params[in] datfab array of data used to construct derived quantity
*/
void
erf_derQc (
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& datfab,
  const amrex::Geometry& /*geomdata*/,
  amrex::Real /*time*/,
  const int* /*bcrec*/,
  const int /*level*/)
{
  erf_derrhodivide(bx, derfab, datfab, RhoQc_comp);
}
#endif

}
