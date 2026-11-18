#ifndef ERF_GTEST_SURFACE_LAYER_COMMON_H_
#define ERF_GTEST_SURFACE_LAYER_COMMON_H_

#include <AMReX_Geometry.H>

#include <ERF_Constants.H>
#include <ERF_SurfaceLayer.H>

#include <cmath>
#include <algorithm>
#include <limits>

namespace erf_surface_layer_test {

using namespace amrex;

constexpr Real test_rho = Real(1.0);
constexpr Real test_rho_theta = Real(1.0);
constexpr Real test_u = Real(3.0);
constexpr Real test_v = Real(4.0);
constexpr Real test_w = Real(5.0);
constexpr Real test_qv = Real(0.01);
constexpr Real test_surface_temperature = Real(300.0);
constexpr Real tau_sentinel = Real(-987654.0);

// Match ERF's default non-EB allocation: ComputeGhostCells() is two,
// conserved state gets one additional ghost cell, and primitive fields have
// one vertical ghost cell.
constexpr int test_state_ng = 3;
constexpr int test_velocity_ng = 2;
constexpr int test_primitive_z_ng = 1;

// Independent qsat oracle for the flat fixture. Horizontal faces use the
// adjacent cell-centre pressure; the vertical faces use the signed distance
// from the cell centre to their boundary face.
inline Real
expected_surface_pressure (const Geometry& geom, const Orientation face,
                           const Real delta_z_override = std::numeric_limits<Real>::quiet_NaN())
{
    Real pressure = getPgivenRTh(test_rho_theta, test_qv);
    Real delta_z = Real(0.0);
    if (face.coordDir() == 2) {
        delta_z = face.isLow() ? myhalf * geom.CellSize(2)
                               : -myhalf * geom.CellSize(2);
        if (std::isfinite(delta_z_override)) { delta_z = delta_z_override; }
    }
    pressure += test_rho * CONST_GRAV * delta_z;
    return pressure;
}

inline Real
expected_qsat (const Geometry& geom, const Orientation face,
               const Real delta_z_override = std::numeric_limits<Real>::quiet_NaN())
{
    Real pressure = expected_surface_pressure(geom, face, delta_z_override);
    pressure *= Real(0.01); // Pa -> mbar, as in fill_qsurf_with_qsat.
    Real qsat = 0;
    erf_qsatw(test_surface_temperature, pressure, qsat);
    return qsat;
}

inline Real
qsat_tolerance (const Real expected)
{
    return Real(128.0) * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.0), std::abs(expected));
}

inline bool
stress_has_expected_sign (const Real value, const bool low_face)
{
    return low_face ? value < Real(0) : value > Real(0);
}

inline bool
stress_is_antisymmetric (const Real low_value, const Real high_value,
                         const Real tolerance_value)
{
    return std::abs(high_value + low_value) <= tolerance_value;
}

} // namespace erf_surface_layer_test

#endif
