// Member drag 0.5 rho Cd w L |U_n| U_n and the ASCE 7 lattice force coefficient.

#include "ERF_MemberDrag.H"

#include <algorithm>
#include <cmath>
#include <string>

#include <AMReX.H>
#include <AMReX_BLassert.H>

using amrex::Real;

namespace erf_towers {

std::array<Real,3> member_drag (const MemberNode& n, const std::array<Real,3>& wind, const std::array<Real,3>& velocity, Real rho)
{
    std::array<Real,3> u{{wind[0] - velocity[0], wind[1] - velocity[1], wind[2] - velocity[2]}};
    const Real along = u[0] * n.axis[0] + u[1] * n.axis[1] + u[2] * n.axis[2];
    for (int d = 0; d < 3; ++d) { u[d] -= along * n.axis[d]; }
    const Real un = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    Real k = Real(0.5) * rho * n.drag_coefficient * n.drag_width * n.length * un;
    if (n.diagonal_gain > Real(0.0)) {
        // the direction factor from the angle of the normal flow's horizontal part to the face normal
        const Real uh = std::sqrt(u[0] * u[0] + u[1] * u[1]);
        if (uh > Real(0.0)) {
            const Real c = (u[0] * n.face[0] + u[1] * n.face[1]) / uh;
            const Real s = (n.face[0] * u[1] - n.face[1] * u[0]) / uh;
            k *= Real(1.0) + n.diagonal_gain * (Real(2.0) * s * c) * (Real(2.0) * s * c);
        }
    }
    return {{k * u[0], k * u[1], k * u[2]}};
}

Real lattice_force_coefficient (Real phi)
{
    return Real(4.0) * phi * phi - Real(5.9) * phi + Real(4.0);
}

Real lattice_diagonal_gain (Real phi)
{
    return std::min(Real(0.75) * phi, Real(0.2));
}

MemberDrag::MemberDrag (Real rho)
    : m_rho(rho)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(rho) && rho > 0.0, "MemberDrag: the air density must be finite and positive (kg/m^3)");
}

void MemberDrag::loads (const std::vector<MemberNode>& nodes, const std::vector<Real>& wind,
                        const std::vector<Real>& velocity, std::vector<Real>& force) const
{
    if (wind.size() != 3 * nodes.size() || velocity.size() != wind.size()) {
        amrex::Abort("MemberDrag::loads: 3 wind and 3 velocity components per drag node are needed (" + std::to_string(nodes.size()) +
                     " nodes, " + std::to_string(wind.size()) + " wind and " + std::to_string(velocity.size()) + " velocity values given)");
    }
    force.assign(3 * nodes.size(), 0.0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const std::array<Real,3> u{{wind[3*i], wind[3*i+1], wind[3*i+2]}};
        const std::array<Real,3> v{{velocity[3*i], velocity[3*i+1], velocity[3*i+2]}};
        const auto f = member_drag(nodes[i], u, v, m_rho);
        for (int d = 0; d < 3; ++d) { force[3*i+d] = f[d]; }
    }
}

} // namespace erf_towers
