#include "ERF_Derive.H"
#include "ERF_EOS.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_StormDiagnostics.H"
#include "ERF_IndexDefines.H"

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr Real
mucape_search_depth_pa ()
{
    return Real(30000.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr Real
mucape_min_pressure_pa ()
{
    return Real(100.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr Real
mucape_min_temperature ()
{
    return Real(180.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr Real
mucape_min_qv ()
{
    return Real(1.0e-12);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr int
mucape_max_substeps ()
{
    return 64;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_virtual_temperature (Real T, Real qv)
{
    return T * (Real(1) + Real(0.61) * amrex::max(qv, Real(0)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_vapor_pressure_pa (Real p, Real qv)
{
    Real qv_clamped = amrex::max(qv, Real(0));
    return p * qv_clamped / amrex::max(Rd_on_Rv + qv_clamped, mucape_min_qv());
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_dewpoint_temperature (Real p, Real qv, Real T)
{
    Real e_mb = amrex::max(mucape_vapor_pressure_pa(p, qv) * Real(0.01), Real(1.0e-6));
    Real log_ratio = std::log(e_mb / Real(6.112));
    Real Td = Real(243.5) * log_ratio / (Real(17.67) - log_ratio) + Real(273.15);
    return amrex::max(mucape_min_temperature(), amrex::min(T, Td));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_lcl_temperature (Real T, Real Td)
{
    Real Td_safe = amrex::max(mucape_min_temperature(), amrex::min(T, Td));
    Real inv_tlcl = Real(1) / (Td_safe - Real(56.0)) + std::log(T / Td_safe) / Real(800.0);
    Real Tlcl = Real(1) / inv_tlcl + Real(56.0);
    return amrex::max(mucape_min_temperature(), amrex::min(T, Tlcl));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_qsat (Real T, Real p_pa)
{
    Real qsat = Real(0);
    erf_qsatw(T, amrex::max(p_pa, mucape_min_pressure_pa()) * Real(0.01), qsat);
    return amrex::max(qsat, Real(0));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_saturated_dTdp (Real T, Real p_pa)
{
    Real p_safe = amrex::max(p_pa, mucape_min_pressure_pa());
    Real T_safe = amrex::max(T, mucape_min_temperature());
    Real qsat = mucape_qsat(T_safe, p_safe);
    Real num = R_d * T_safe * (Real(1) + Real(0.61) * qsat) *
               (Real(1) + L_v * qsat / (R_d * T_safe));
    Real denom = p_safe * (Cp_d + (L_v * L_v * qsat * Rd_on_Rv) / (R_d * T_safe * T_safe));
    return num / amrex::max(denom, Real(1.0e-12));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_integrate_saturated_temperature (Real T_start, Real p_start, Real p_target)
{
    if (p_target >= p_start) {
        return T_start;
    }

    Real T = T_start;
    Real p = p_start;
    int nsub = amrex::min(mucape_max_substeps(),
                          amrex::max(1, static_cast<int>(std::abs(p_start - p_target) / Real(500.0)) + 1));
    Real dp = (p_target - p_start) / nsub;

    for (int n = 0; n < nsub; ++n) {
        Real dTdp_1 = mucape_saturated_dTdp(T, p);
        Real p_next = p + dp;
        Real T_pred = amrex::max(mucape_min_temperature(), T + dTdp_1 * dp);
        Real dTdp_2 = mucape_saturated_dTdp(T_pred, p_next);
        T = amrex::max(mucape_min_temperature(), T + myhalf * (dTdp_1 + dTdp_2) * dp);
        p = p_next;
    }

    return T;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
mucape_positive_area (Real b_lo, Real b_hi, Real dz)
{
    if (dz <= Real(0)) {
        return Real(0);
    }

    if (b_lo > Real(0) && b_hi > Real(0)) {
        return myhalf * (b_lo + b_hi) * dz;
    }

    if (b_lo <= Real(0) && b_hi <= Real(0)) {
        return Real(0);
    }

    Real frac = b_lo / (b_lo - b_hi);
    frac = amrex::max(Real(0), amrex::min(Real(1), frac));

    if (b_lo > Real(0)) {
        return myhalf * b_lo * frac * dz;
    } else {
        return myhalf * b_hi * (Real(1) - frac) * dz;
    }
}

} // namespace

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
erf_dermaxreflectivity ( const Box& bx,
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

        Real max_dbz = Real(-1.0e30);

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
erf_derlocalhelicity ( const Box& bx,
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

    auto const dat = datfab.array();  // cell-centered velocity
    auto dfab      = derfab.array();  // integral of local helicity
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

                Real z_hi = myhalf * (z_arr(i,j,k) + z_arr(i,j,k+1));
                Real z_lo = myhalf * (z_arr(i,j,k) + z_arr(i,j,k-1));
                Real dz = z_hi - z_lo;

                Real vortz = (dat(i+1,j,k,1) - dat(i-1,j,k,1)) / (two*dx)  // dv/dx
                           - (dat(i,j+1,k,0) - dat(i,j-1,k,0)) / (two*dy); // du/dy
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

void
erf_derprecipitable ( const Box& bx,
                      FArrayBox& derfab,
                      int dcomp,
                      int /*ncomp*/,
                      const FArrayBox& datfab,
                      const FArrayBox& zcc_fab,
                      const Geometry& /*geomdata*/,
                      Real /*time*/,
                      const int* /*bcrec*/,
                      const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);

    auto const dat = datfab.array(); // cell-centered state vector
    auto dfab      = derfab.array(); // integral of qv to define precipitable water

    // Collapse to i,j box (ignore vertical for now)
    Box b2d = bx;
    b2d.setSmall(2,0);
    b2d.setBig(2,0);

    auto z_arr     = zcc_fab.array(); // cell-centered height z

    ParallelFor(b2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
    {

        Real int_qv = Real(0.0);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k)
        {
            Real z_hi = myhalf * (z_arr(i,j,k) + z_arr(i,j,k+1));
            Real z_lo = myhalf * (z_arr(i,j,k) + z_arr(i,j,k-1));
            Real dz = z_hi - z_lo;

            Real rhoQ1 = dat(i, j, k, RhoQ1_comp);

            int_qv += rhoQ1 * dz;
        }

        // Store vertical integral into *all* levels for this (i,j)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            dfab(i, j, k, dcomp) = int_qv;
        }
    });
}

void
erf_dermucape ( const Box& bx,
                FArrayBox& derfab,
                int dcomp,
                int ncomp,
                const FArrayBox& datfab,
                const FArrayBox& zcc_fab,
                const Geometry& /*geomdata*/,
                Real /*time*/,
                const int* /*bcrec*/,
                const int /*level*/)
{
    AMREX_ALWAYS_ASSERT(dcomp == 0);
    AMREX_ALWAYS_ASSERT(ncomp == 1);

    auto const dat = datfab.array();
    auto dfab      = derfab.array();
    auto z_arr     = zcc_fab.array();
    const int ncons = datfab.nComp();

    Box b2d = bx;
    b2d.setSmall(2,0);
    b2d.setBig(2,0);

    ParallelFor(b2d, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
    {
        Real mucape = Real(0);
        int klo = bx.smallEnd(2);
        int khi = bx.bigEnd(2);

        if (ncons > RhoQ1_comp) {
            Real p_sfc = Real(0);
            for (int k = klo; k <= khi; ++k) {
                Real rho = dat(i,j,k,Rho_comp);
                if (rho <= Real(0)) { continue; }
                Real qv = amrex::max(Real(0), dat(i,j,k,RhoQ1_comp) / rho);
                p_sfc = amrex::max(p_sfc, getPgivenRTh(dat(i,j,k,RhoTheta_comp), qv));
            }

            if (p_sfc > mucape_search_depth_pa()) {
                Real p_search_min = p_sfc - mucape_search_depth_pa();

                for (int ks = klo; ks < khi; ++ks) {
                    Real rho_src = dat(i,j,ks,Rho_comp);
                    if (rho_src <= Real(0)) { continue; }

                    Real qv_src = amrex::max(Real(0), dat(i,j,ks,RhoQ1_comp) / rho_src);
                    Real rt_src = dat(i,j,ks,RhoTheta_comp);
                    Real p_src  = getPgivenRTh(rt_src, qv_src);

                    if (p_src < p_search_min) { continue; }

                    Real T_src = getTgivenRandRTh(rho_src, rt_src, qv_src);
                    if (T_src <= Real(0)) { continue; }

                    Real Td_src = mucape_dewpoint_temperature(p_src, qv_src, T_src);
                    Real Tlcl   = mucape_lcl_temperature(T_src, Td_src);
                    Real plcl   = p_src * std::pow(Tlcl / T_src, Cp_d / R_d);
                    plcl = amrex::min(p_src, amrex::max(plcl, mucape_min_pressure_pa()));

                    Real theta_src = getThgivenTandP(T_src, p_src, R_d / Cp_d);

                    Real candidate_cape = Real(0);
                    Real z_prev = z_arr(i,j,ks);
                    Real b_prev = Real(0);

                    bool saturated = false;
                    Real T_sat_prev = Tlcl;
                    Real p_sat_prev = plcl;

                    for (int k = ks + 1; k <= khi; ++k) {
                        Real rho_env = dat(i,j,k,Rho_comp);
                        if (rho_env <= Real(0)) { continue; }

                        Real qv_env = amrex::max(Real(0), dat(i,j,k,RhoQ1_comp) / rho_env);
                        Real rt_env = dat(i,j,k,RhoTheta_comp);
                        Real p_env  = getPgivenRTh(rt_env, qv_env);
                        Real T_env  = getTgivenRandRTh(rho_env, rt_env, qv_env);
                        Real Tv_env = mucape_virtual_temperature(T_env, qv_env);
                        Real z_env  = z_arr(i,j,k);

                        Real T_parcel;
                        Real qv_parcel;

                        if (p_env >= plcl) {
                            T_parcel  = getTgivenPandTh(p_env, theta_src, R_d / Cp_d);
                            qv_parcel = qv_src;
                        } else {
                            if (!saturated) {
                                saturated = true;
                                T_sat_prev = Tlcl;
                                p_sat_prev = plcl;
                            }

                            if (p_env < p_sat_prev) {
                                T_sat_prev = mucape_integrate_saturated_temperature(T_sat_prev, p_sat_prev, p_env);
                                p_sat_prev = p_env;
                            }

                            T_parcel  = T_sat_prev;
                            qv_parcel = mucape_qsat(T_parcel, p_env);
                        }

                        Real Tv_parcel = mucape_virtual_temperature(T_parcel, qv_parcel);
                        Real buoyancy = CONST_GRAV * (Tv_parcel - Tv_env) /
                                        amrex::max(Tv_env, mucape_min_temperature());

                        candidate_cape += mucape_positive_area(b_prev, buoyancy, z_env - z_prev);
                        b_prev = buoyancy;
                        z_prev = z_env;
                    }

                    mucape = amrex::max(mucape, candidate_cape);
                }
            }
        }

        for (int k = klo; k <= khi; ++k) {
            dfab(i,j,k,dcomp) = mucape;
        }
    });
}

} // namespace
