/**
 * Terrain-following inflow profiles (xlo.inflow_profile = log_law or file).
 *
 * These run after the lateral boundary fills. On an Inflow or InflowOutflow
 * face with a profile they overwrite the ghost cells, and the boundary face
 * itself for the normal velocity, with the profile evaluated at each cell's own
 * height above the ground beneath it. On an InflowOutflow face only the cells
 * where the flow enters are overwritten.
 */
#include "ERF_PhysBCFunct.H"

using namespace amrex;

namespace {

// The BC type of face `ori` for BC component `bccomp` when that face has an
// active profile and prescribes the component (ext_dir or ext_dir_upwind);
// -1 otherwise.
int
profile_bc_type (const Array<InflowProfile,2*AMREX_SPACEDIM>* profiles,
                 const Vector<BCRec>& bcs, const int bccomp, const Orientation ori)
{
    if (profiles == nullptr || !(*profiles)[ori].active) { return -1; }
    const int dir = ori.coordDir();
    const int t = ori.isLow() ? bcs[bccomp].lo(dir) : bcs[bccomp].hi(dir);
    return (t == ERFBCType::ext_dir || t == ERFBCType::ext_dir_upwind) ? t : -1;
}

} // namespace

void
ERFPhysBCFunct_u::impose_inflow_profile_xvel (const Array4<Real>& dest_arr,
                                             const Array4<Real const>& xvel_arr,
                                             const Array4<Real const>& yvel_arr,
                                             const Array4<Real const>& z_nd,
                                             const Box& bx, const Box& domain,
                                             const int bccomp)
{
    if (m_inflow_profiles == nullptr) { return; }
    BL_PROFILE("impose_inflow_profile_xvel()");

    const auto dom_lo  = lbound(domain);
    const auto dom_hi  = ubound(domain);
    const bool terrain = static_cast<bool>(z_nd);
    const Real dz      = m_geom.CellSize(2);

    for (int dir = 0; dir < 2; ++dir) {
        if (m_geom.isPeriodic(dir)) { continue; }
        for (const auto side : {Orientation::low, Orientation::high}) {
            const Orientation ori(dir, side);
            const int bct = profile_bc_type(m_inflow_profiles, m_domain_bcs_type, bccomp, ori);
            if (bct < 0) { continue; }

            const auto& prof  = (*m_inflow_profiles)[ori];
            const Real* zt    = prof.z_d.data();
            const Real* ft    = prof.u_d.data();
            const int   nt    = static_cast<int>(prof.z_d.size());
            const bool upwind = (bct == ERFBCType::ext_dir_upwind);
            const bool low    = ori.isLow();

            // x-velocity is normal to an x-face (the face and the ghost faces
            // beyond it) and tangential to a y-face (the ghost rows)
            Box b(bx);
            if (dir == 0) {
                if (low) { b.setBig(0, dom_lo.x); } else { b.setSmall(0, dom_hi.x+1); }
            } else {
                if (low) { b.setBig(1, dom_lo.y-1); } else { b.setSmall(1, dom_hi.y+1); }
            }
            if (!b.ok()) { continue; }

            const int iface = low ? dom_lo.x : dom_hi.x+1;
            const int jcol  = low ? dom_lo.y : dom_hi.y;
            const int jface = low ? dom_lo.y : dom_hi.y+1;

            ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int kk = amrex::max(dom_lo.z, amrex::min(k, dom_hi.z));
                if (upwind) {
                    const Real un = (dir == 0) ? xvel_arr(iface,j,kk) : yvel_arr(i,jface,kk);
                    if (low ? (un < Real(0.0)) : (un > Real(0.0))) { return; }
                }
                Real zq;
                if (dir == 0) {
                    const int jj = amrex::max(dom_lo.y, amrex::min(j, dom_hi.y));
                    zq = (terrain) ? height_above_ground(z_nd, iface, iface, jj, jj+1, kk, dom_lo.z)
                                   : (static_cast<Real>(kk - dom_lo.z) + Real(0.5)) * dz;
                } else {
                    const int ii = amrex::max(dom_lo.x, amrex::min(i, dom_hi.x+1));
                    zq = (terrain) ? height_above_ground(z_nd, ii, ii, jcol, jcol+1, kk, dom_lo.z)
                                   : (static_cast<Real>(kk - dom_lo.z) + Real(0.5)) * dz;
                }
                dest_arr(i,j,k) = inflow_profile_value(zt, ft, nt, zq);
            });
        }
    }
}

void
ERFPhysBCFunct_v::impose_inflow_profile_yvel (const Array4<Real>& dest_arr,
                                             const Array4<Real const>& xvel_arr,
                                             const Array4<Real const>& yvel_arr,
                                             const Array4<Real const>& z_nd,
                                             const Box& bx, const Box& domain,
                                             const int bccomp)
{
    if (m_inflow_profiles == nullptr) { return; }
    BL_PROFILE("impose_inflow_profile_yvel()");

    const auto dom_lo  = lbound(domain);
    const auto dom_hi  = ubound(domain);
    const bool terrain = static_cast<bool>(z_nd);
    const Real dz      = m_geom.CellSize(2);

    for (int dir = 0; dir < 2; ++dir) {
        if (m_geom.isPeriodic(dir)) { continue; }
        for (const auto side : {Orientation::low, Orientation::high}) {
            const Orientation ori(dir, side);
            const int bct = profile_bc_type(m_inflow_profiles, m_domain_bcs_type, bccomp, ori);
            if (bct < 0) { continue; }

            const auto& prof  = (*m_inflow_profiles)[ori];
            const Real* zt    = prof.z_d.data();
            const Real* ft    = prof.v_d.data();
            const int   nt    = static_cast<int>(prof.z_d.size());
            const bool upwind = (bct == ERFBCType::ext_dir_upwind);
            const bool low    = ori.isLow();

            // y-velocity is tangential to an x-face (the ghost columns) and
            // normal to a y-face (the face and the ghost faces beyond it)
            Box b(bx);
            if (dir == 0) {
                if (low) { b.setBig(0, dom_lo.x-1); } else { b.setSmall(0, dom_hi.x+1); }
            } else {
                if (low) { b.setBig(1, dom_lo.y); } else { b.setSmall(1, dom_hi.y+1); }
            }
            if (!b.ok()) { continue; }

            const int icol  = low ? dom_lo.x : dom_hi.x;
            const int iface = low ? dom_lo.x : dom_hi.x+1;
            const int jface = low ? dom_lo.y : dom_hi.y+1;

            ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int kk = amrex::max(dom_lo.z, amrex::min(k, dom_hi.z));
                if (upwind) {
                    const Real un = (dir == 0) ? xvel_arr(iface,j,kk) : yvel_arr(i,jface,kk);
                    if (low ? (un < Real(0.0)) : (un > Real(0.0))) { return; }
                }
                Real zq;
                if (dir == 0) {
                    const int jj = amrex::max(dom_lo.y, amrex::min(j, dom_hi.y+1));
                    zq = (terrain) ? height_above_ground(z_nd, icol, icol+1, jj, jj, kk, dom_lo.z)
                                   : (static_cast<Real>(kk - dom_lo.z) + Real(0.5)) * dz;
                } else {
                    const int ii = amrex::max(dom_lo.x, amrex::min(i, dom_hi.x));
                    zq = (terrain) ? height_above_ground(z_nd, ii, ii+1, jface, jface, kk, dom_lo.z)
                                   : (static_cast<Real>(kk - dom_lo.z) + Real(0.5)) * dz;
                }
                dest_arr(i,j,k) = inflow_profile_value(zt, ft, nt, zq);
            });
        }
    }
}

void
ERFPhysBCFunct_cons::impose_inflow_profile_cons (const Array4<Real>& dest_arr,
                                               const Array4<Real const>& xvel_arr,
                                               const Array4<Real const>& yvel_arr,
                                               const Array4<Real const>& z_nd,
                                               const Box& bx, const Box& domain,
                                               const int icomp, const int ncomp)
{
    if (m_inflow_profiles == nullptr) { return; }
    BL_PROFILE("impose_inflow_profile_cons()");

    const auto dom_lo  = lbound(domain);
    const auto dom_hi  = ubound(domain);
    const bool terrain = static_cast<bool>(z_nd);
    const Real dz      = m_geom.CellSize(2);

    for (int dir = 0; dir < 2; ++dir) {
        if (m_geom.isPeriodic(dir)) { continue; }
        for (const auto side : {Orientation::low, Orientation::high}) {
            const Orientation ori(dir, side);
            const auto& prof = (*m_inflow_profiles)[ori];
            if (!prof.active) { continue; }

            // The profile prescribes theta and tke where it has them, as
            // rho times the tabulated value with rho from the adjacent cell;
            // density itself stays extrapolated from the interior
            for (const int comp : {RhoTheta_comp, RhoKE_comp}) {
                if (comp < icomp || comp >= icomp + ncomp) { continue; }
                const bool is_ke = (comp == RhoKE_comp);
                if (is_ke ? !prof.has_tke : !prof.has_theta) { continue; }
                const int bcvar = is_ke ? BCVars::RhoKE_bc_comp : BCVars::RhoTheta_bc_comp;
                const int bct = profile_bc_type(m_inflow_profiles, m_domain_bcs_type, bcvar, ori);
                if (bct < 0) { continue; }

                const Real* zt    = prof.z_d.data();
                const Real* ft    = is_ke ? prof.tke_d.data() : prof.theta_d.data();
                const int   nt    = static_cast<int>(prof.z_d.size());
                const bool upwind = (bct == ERFBCType::ext_dir_upwind);
                const bool low    = ori.isLow();

                Box b(bx);
                if (dir == 0) {
                    if (low) { b.setBig(0, dom_lo.x-1); } else { b.setSmall(0, dom_hi.x+1); }
                } else {
                    if (low) { b.setBig(1, dom_lo.y-1); } else { b.setSmall(1, dom_hi.y+1); }
                }
                if (!b.ok()) { continue; }

                const int icol  = low ? dom_lo.x : dom_hi.x;
                const int iface = low ? dom_lo.x : dom_hi.x+1;
                const int jcol  = low ? dom_lo.y : dom_hi.y;
                const int jface = low ? dom_lo.y : dom_hi.y+1;

                ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int kk = amrex::max(dom_lo.z, amrex::min(k, dom_hi.z));
                    const int ii = (dir == 0) ? icol : amrex::max(dom_lo.x, amrex::min(i, dom_hi.x));
                    const int jj = (dir == 1) ? jcol : amrex::max(dom_lo.y, amrex::min(j, dom_hi.y));
                    if (upwind) {
                        const Real un = (dir == 0) ? xvel_arr(iface,jj,kk) : yvel_arr(ii,jface,kk);
                        if (low ? (un < Real(0.0)) : (un > Real(0.0))) { return; }
                    }
                    const Real zq = (terrain) ? height_above_ground(z_nd, ii, ii+1, jj, jj+1, kk, dom_lo.z)
                                              : (static_cast<Real>(kk - dom_lo.z) + Real(0.5)) * dz;
                    dest_arr(i,j,k,comp) = dest_arr(ii,jj,kk,Rho_comp) * inflow_profile_value(zt, ft, nt, zq);
                });
            }
        }
    }
}
