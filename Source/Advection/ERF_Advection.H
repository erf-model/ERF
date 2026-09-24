#ifndef ERF_ADVECTION_H_
#define ERF_ADVECTION_H_

#include "AMReX.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_BCRec.H"
#include "AMReX_YAFluxRegister.H"

#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"
#include "ERF_EB.H"

/** Compute advection tendency for density and averaged momentum fluxes */
void AdvectionSrcForRho (const amrex::Box& bx,
                         const amrex::Array4<amrex::Real>& src,
                         const amrex::Array4<const amrex::Real>& rho_u,    // These are being used
                         const amrex::Array4<const amrex::Real>& rho_v,    //  to define the fluxes
                         const amrex::Array4<const amrex::Real>& omega,
                         const amrex::Array4<      amrex::Real>& avg_xmom, // These are being defined
                         const amrex::Array4<      amrex::Real>& avg_ymom, //  from the rho fluxes
                         const amrex::Array4<      amrex::Real>& avg_zmom,
                         const amrex::Array4<const amrex::Real>& ax_arr,
                         const amrex::Array4<const amrex::Real>& ay_arr,
                         const amrex::Array4<const amrex::Real>& az_arr,
                         const amrex::Array4<const amrex::Real>& detJ,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const amrex::Array4<const amrex::Real>& mf_mx,
                         const amrex::Array4<const amrex::Real>& mf_my,
                         const amrex::Array4<const amrex::Real>& mf_uy,
                         const amrex::Array4<const amrex::Real>& mf_vx,
                         const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                         const bool fixed_rho);

/** Compute advection tendencies for cell-centered scalar components */
void AdvectionSrcForScalars (const amrex::Box& bx,
                             const int icomp, const int ncomp,
                             const amrex::Array4<const amrex::Real>& avg_xmom,
                             const amrex::Array4<const amrex::Real>& avg_ymom,
                             const amrex::Array4<const amrex::Real>& avg_zmom,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<amrex::Real>& src,
                             const amrex::Array4<const amrex::Real>& vf_arr,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                             const amrex::Array4<const amrex::Real>& mf_mx,
                             const amrex::Array4<const amrex::Real>& mf_my,
                             const AdvType horiz_adv_type, const AdvType vert_adv_type,
                             const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                             const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                             const amrex::Box& domain,
                             const amrex::BCRec* bc_ptr_h);

/** Compute advection tendencies for all components of momentum */
void AdvectionSrcForMom (const amrex::MFIter& mfi,
                         const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                         const amrex::Vector<amrex::Box>& bxx_grown,
                         const amrex::Vector<amrex::Box>& bxy_grown,
                         const amrex::Vector<amrex::Box>& bxz_grown,
                         const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                         const amrex::Array4<      amrex::Real>& rho_w_rhs,
                         const amrex::Array4<const amrex::Real>& rho,
                         const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                         const amrex::Array4<const amrex::Real>& w        ,
                         const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                         const amrex::Array4<const amrex::Real>& Omega    ,
                         const amrex::Array4<const amrex::Real>& z_nd,
                         const amrex::Array4<const amrex::Real>& ax,
                         const amrex::Array4<const amrex::Real>& ay,
                         const amrex::Array4<const amrex::Real>& az,
                         const amrex::Array4<const amrex::Real>& detJ,
                               amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                         const amrex::Array4<const amrex::Real>& mf_mx,
                         const amrex::Array4<const amrex::Real>& mf_ux,
                         const amrex::Array4<const amrex::Real>& mf_vx,
                         const amrex::Array4<const amrex::Real>& mf_my,
                         const amrex::Array4<const amrex::Real>& mf_uy,
                         const amrex::Array4<const amrex::Real>& mf_vy,
                         const AdvType horiz_adv_type, const AdvType vert_adv_type,
                         const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                         MeshType& mesh_type, TerrainType& terrain_type,
                         const eb_& ebfact,
                               amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_u_arr,
                               amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_v_arr,
                               amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_w_arr,
                         const amrex::Vector<amrex::iMultiFab>& physbnd_mask,
                         const bool already_on_centroids,
                         const int lo_z_face, const int hi_z_face,
                         const amrex::Box& domain,
                         const amrex::BCRec* bc_ptr_h);

/** Compute advection tendencies for all components of momentum when dz is constant */
void AdvectionSrcForMom_ConstantDz (
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            TerrainType& terrain_type,
                            const int lo_z_face, const int hi_z_face);

/** Compute advection tendencies for all components of momentum when dz is stretched */
void AdvectionSrcForMom_StretchedDz (
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Gpu::DeviceVector<amrex::Real>& stretched_dz_d,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const int lo_z_face, const int hi_z_face);

/** Compute advection tendencies for all components of momentum for TERRAIN-FITTED COORDS only */
void AdvectionSrcForMom_TF (
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::Array4<const amrex::Real>& z_nd,
                            const amrex::Array4<const amrex::Real>& ax,
                            const amrex::Array4<const amrex::Real>& ay,
                            const amrex::Array4<const amrex::Real>& az,
                            const amrex::Array4<const amrex::Real>& detJ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const int lo_z_face, const int hi_z_face);

/** Compute advection tendencies for all components of momentum when using EB for terrain */
void AdvectionSrcForMom_EB (const amrex::MFIter& mfi,
                            const amrex::Box& bxx, const amrex::Box& bxy, const amrex::Box& bxz,
                            const amrex::Vector<amrex::Box>& bxx_grown,
                            const amrex::Vector<amrex::Box>& bxy_grown,
                            const amrex::Vector<amrex::Box>& bxz_grown,
                            const amrex::Array4<      amrex::Real>& rho_u_rhs, const amrex::Array4<      amrex::Real>& rho_v_rhs,
                            const amrex::Array4<      amrex::Real>& rho_w_rhs,
                            const amrex::Array4<const amrex::Real>& u        , const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w        ,
                            const amrex::Array4<const amrex::Real>& rho_u    , const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& Omega    ,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_mx,
                            const amrex::Array4<const amrex::Real>& mf_ux,
                            const amrex::Array4<const amrex::Real>& mf_vx,
                            const amrex::Array4<const amrex::Real>& mf_my,
                            const amrex::Array4<const amrex::Real>& mf_uy,
                            const amrex::Array4<const amrex::Real>& mf_vy,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const eb_& ebfact,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_u_arr,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_v_arr,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_w_arr,
                            const amrex::Vector<amrex::iMultiFab>& physbnd_mask,
                            const bool already_on_centroids,
                            const int lo_z_face, const int hi_z_face,
                            const amrex::Box& domain);

/**
 * Shrink a z-momentum box away from the bottom and top of the domain, where no
 * z-momentum source term is computed.
 *
 * The shrink must be applied only at an end that actually coincides with the domain,
 * which is why this takes `domain` rather than shrinking unconditionally: on a grid
 * decomposed in z, an unconditional shrink also removes the w-face at every internal
 * split, and any operator run on the result silently skips those faces.
 *
 * `b` is z-nodal; `domain` is cell-centered.
 */
inline amrex::Box
ShrinkZmomBoxAtDomainEnds (amrex::Box b, const amrex::Box& domain)
{
    if (b.smallEnd(2) == domain.smallEnd(2))   { b.growLo(2,-1); }
    if (b.bigEnd(2)   == domain.bigEnd(2) + 1) { b.growHi(2,-1); }
    return b;
}

/**
 * Which side of the domain is open in one coordinate direction, as seen by the
 * cells an open-boundary kernel is about to update.
 */
enum struct OpenSide : int { lo = -1, none = 0, hi = 1 };

/** One box on an open lateral boundary, tagged with the side(s) that are open there. */
struct OpenBCPatch {
    amrex::Box  box;
    OpenSide x_side = OpenSide::none;
    OpenSide y_side = OpenSide::none;
};

/**
 * Remove from an open-boundary box the cells that also sit on a perpendicular
 * open boundary, so that no cell is updated by two different open-boundary kernels.
 *
 * `b` is trimmed in dimension `dim` only at an end that actually coincides with the
 * corresponding end of `domain`. The domain-level open flags by themselves are not a
 * sufficient test: `b` is a per-box (or per-tile) slab whose ends generally lie in
 * the interior, and trimming those would silently drop interior cells from the
 * open-boundary treatment and make the answer depend on the grid decomposition.
 *
 * `domain` is cell-centered; the index of the boundary is taken in the index space
 * of `b`, which may be node-centered in `dim`.
 */
inline amrex::Box
TrimOpenBCCorner (amrex::Box b,
                  const int dim,
                  const amrex::Box& domain,
                  const bool lo_open,
                  const bool hi_open)
{
    if (!b.ok()) { return b; }

    const int dom_lo = domain.smallEnd(dim);
    const int dom_hi = domain.bigEnd(dim) + (b.ixType().nodeCentered(dim) ? 1 : 0);

    if (lo_open && b.smallEnd(dim) == dom_lo) { b.growLo(dim,-1); }
    if (hi_open && b.bigEnd(dim)   == dom_hi) { b.growHi(dim,-1); }

    return b;
}

/**
 * Build the one-cell-wide box where the `x_side` and `y_side` open boundaries meet,
 * or an empty box if `b` does not reach both of them.
 */
inline amrex::Box
OpenBCCornerBox (amrex::Box b,
                 const amrex::Box& domain,
                 const OpenSide x_side,
                 const OpenSide y_side)
{
    AMREX_ALWAYS_ASSERT(x_side != OpenSide::none && y_side != OpenSide::none);

    if (!b.ok()) { return amrex::Box(); }

    const OpenSide side[2] = {x_side, y_side};
    for (int dim(0); dim < 2; ++dim) {
        const int idx = (side[dim] == OpenSide::lo) ?
                         domain.smallEnd(dim) :
                         domain.bigEnd(dim) + (b.ixType().nodeCentered(dim) ? 1 : 0);
        if (idx < b.smallEnd(dim) || idx > b.bigEnd(dim)) { return amrex::Box(); }
        b = amrex::makeSlab(b, dim, idx);
    }

    return b;
}

/**
 * Enumerate the boxes on which the open-boundary tangential operator must be applied
 * for a quantity that has no boundary-normal owner at a lateral corner -- that is,
 * w and the cell-centered state, both of which are tangent to every lateral boundary.
 *
 * The returned boxes are disjoint and together cover every cell of `b` that lies on an
 * open lateral boundary: the four edges, trimmed so they do not reach into a corner,
 * plus the corners themselves, which carry both sides so the kernel differences across
 * neither open boundary there.
 */
inline amrex::Vector<OpenBCPatch>
OpenBCTangentPatches (const amrex::Box& b,
                      const amrex::Box& domain,
                      const bool xlo_open,
                      const bool xhi_open,
                      const bool ylo_open,
                      const bool yhi_open)
{
    amrex::Vector<OpenBCPatch> patches;
    if (!b.ok()) { return patches; }
    patches.reserve(8);

    const int dom_xhi = domain.bigEnd(0) + (b.ixType().nodeCentered(0) ? 1 : 0);
    const int dom_yhi = domain.bigEnd(1) + (b.ixType().nodeCentered(1) ? 1 : 0);

    // Edges: each is trimmed by the open boundaries perpendicular to it, which own
    // the corner cells by way of the corner patches added below.
    if (xlo_open && b.smallEnd(0) == domain.smallEnd(0)) {
        amrex::Box e = TrimOpenBCCorner(amrex::makeSlab(b,0,domain.smallEnd(0)), 1, domain, ylo_open, yhi_open);
        if (e.ok()) { patches.push_back({e, OpenSide::lo, OpenSide::none}); }
    }
    if (xhi_open && b.bigEnd(0) == dom_xhi) {
        amrex::Box e = TrimOpenBCCorner(amrex::makeSlab(b,0,dom_xhi), 1, domain, ylo_open, yhi_open);
        if (e.ok()) { patches.push_back({e, OpenSide::hi, OpenSide::none}); }
    }
    if (ylo_open && b.smallEnd(1) == domain.smallEnd(1)) {
        amrex::Box e = TrimOpenBCCorner(amrex::makeSlab(b,1,domain.smallEnd(1)), 0, domain, xlo_open, xhi_open);
        if (e.ok()) { patches.push_back({e, OpenSide::none, OpenSide::lo}); }
    }
    if (yhi_open && b.bigEnd(1) == dom_yhi) {
        amrex::Box e = TrimOpenBCCorner(amrex::makeSlab(b,1,dom_yhi), 0, domain, xlo_open, xhi_open);
        if (e.ok()) { patches.push_back({e, OpenSide::none, OpenSide::hi}); }
    }

    // Corners
    for (OpenSide x_side : {OpenSide::lo, OpenSide::hi}) {
        if (!((x_side == OpenSide::lo) ? xlo_open : xhi_open)) { continue; }
        for (OpenSide y_side : {OpenSide::lo, OpenSide::hi}) {
            if (!((y_side == OpenSide::lo) ? ylo_open : yhi_open)) { continue; }
            amrex::Box c = OpenBCCornerBox(b, domain, x_side, y_side);
            if (c.ok()) { patches.push_back({c, x_side, y_side}); }
        }
    }

    return patches;
}

/** Compute advection tendencies for all normal components of momentum */
void
AdvectionSrcForOpenBC_Normal (const amrex::Box& bx,
                              const int& dir,
                              const amrex::Array4<      amrex::Real>& rhs_arr,
                              const amrex::Array4<const amrex::Real>& vel_norm_arr,
                              const amrex::Array4<const amrex::Real>& cell_data_arr,
                              const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
                              const bool do_lo=false);

/** Compute advection tendencies for u momentum tangent to open BC */
void
AdvectionSrcForOpenBC_Tangent_Xmom (const amrex::Box& bxx,
                                    const int& dir,
                                    const amrex::Array4<      amrex::Real>& rho_u_rhs,
                                    const amrex::Array4<const amrex::Real>& u,
                                    const amrex::Array4<const amrex::Real>& rho_u,
                                    const amrex::Array4<const amrex::Real>& rho_v,
                                    const amrex::Array4<const amrex::Real>& Omega,
                                    const amrex::Array4<const amrex::Real>& ax,
                                    const amrex::Array4<const amrex::Real>& az,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo=false);

/** Compute advection tendencies for v momentum tangent to open BC */
void
AdvectionSrcForOpenBC_Tangent_Ymom (const amrex::Box& bxy,
                                    const int& dir,
                                    const amrex::Array4<      amrex::Real>& rho_v_rhs,
                                    const amrex::Array4<const amrex::Real>& v,
                                    const amrex::Array4<const amrex::Real>& rho_u,
                                    const amrex::Array4<const amrex::Real>& rho_v,
                                    const amrex::Array4<const amrex::Real>& Omega,
                                    const amrex::Array4<const amrex::Real>& ay,
                                    const amrex::Array4<const amrex::Real>& az,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo=false);

/** Compute advection tendencies for w momentum tangent to open BC
 *  `x_side` and `y_side` say which lateral boundaries the cells being updated lie on;
 *  both are set at a corner, where neither direction may be differenced across. */
void
AdvectionSrcForOpenBC_Tangent_Zmom (const amrex::Box& bxz,
                                    const OpenSide x_side,
                                    const OpenSide y_side,
                                    const amrex::Array4<      amrex::Real>& rho_w_rhs,
                                    const amrex::Array4<const amrex::Real>& w,
                                    const amrex::Array4<const amrex::Real>& rho_u,
                                    const amrex::Array4<const amrex::Real>& rho_v,
                                    const amrex::Array4<const amrex::Real>& Omega,
                                    const amrex::Array4<const amrex::Real>& ax,
                                    const amrex::Array4<const amrex::Real>& ay,
                                    const amrex::Array4<const amrex::Real>& az,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const int domhi_z);

/** Compute advection tendencies for scalar fields tangential to open BC.
 *
 * @param[in] bx box over which scalar tendencies are updated
 * @param[in] x_side lateral x boundary, if any
 * @param[in] y_side lateral y boundary, if any
 * @param[in] scalar_icomp first component of the supplied scalar view to read
 * @param[in] rhs_icomp first component of the supplied tendency view to write
 * @param[in] ncomp number of components to update
 * @param[out] rhs tendency view
 * @param[in] scalar supplied scalar view
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] detJ metric Jacobian
 * @param[in] cellSizeInv inverse grid spacing
 *
 * At a corner both sides are set, so neither direction is differenced across
 * its open boundary.
 */
void
AdvectionSrcForOpenBC_Tangent_Scalars (const amrex::Box& bx,
                                       const OpenSide x_side,
                                       const OpenSide y_side,
                                       const int scalar_icomp,
                                       const int rhs_icomp,
                                       const int ncomp,
                                       const amrex::Array4<      amrex::Real>& rhs,
                                       const amrex::Array4<const amrex::Real>& scalar,
                                       const amrex::Array4<const amrex::Real>& avg_xmom,
                                       const amrex::Array4<const amrex::Real>& avg_ymom,
                                       const amrex::Array4<const amrex::Real>& avg_zmom,
                                       const amrex::Array4<const amrex::Real>& detJ,
                                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

/** Compute advection tendencies for conserved variables tangential to open BC
 *  `x_side` and `y_side` say which lateral boundaries the cells being updated lie on;
 *  both are set at a corner, where neither direction may be differenced across.
 *  This native adapter maps conserved component `icomp + n` to primitive component
 *  `icomp + n - 1` before calling the component-explicit scalar operation. */
void
AdvectionSrcForOpenBC_Tangent_Cons (const amrex::Box& bx,
                                    const OpenSide x_side,
                                    const OpenSide y_side,
                                    const int& icomp,
                                    const int& ncomp,
                                    const amrex::Array4<      amrex::Real>& cell_rhs,
                                    const amrex::Array4<const amrex::Real>& cell_prim,
                                    const amrex::Array4<const amrex::Real>& avg_xmom,
                                    const amrex::Array4<const amrex::Real>& avg_ymom,
                                    const amrex::Array4<const amrex::Real>& avg_zmom,
                                    const amrex::Array4<const amrex::Real>& detJ,
                                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv);

/** Compute advection tendencies in the normal direction for variables tangent to open BC */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
AdvectionSrcForOpenBC_Tangent (const int& i,
                               const int& j,
                               const int& k,
                               const int& nprim,
                               const int& dir,
                               const amrex::Array4<const amrex::Real>& prim_tang_arr,
                               const amrex::Array4<const amrex::Real>& mom_norm_arr,
                               const amrex::Real& dxInv,
                               const bool do_lo=false);

/** Return the effective advection type used for the requested Runge-Kutta stage. */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
AdvType
EfficientAdvType (int nrk, AdvType adv_type)
{
    AdvType eff_adv_type;

    if (nrk == 0) {

        eff_adv_type = AdvType::Centered_2nd;

    } else if (nrk == 1) {

        if ( (adv_type == AdvType::Centered_2nd) ||
             (adv_type == AdvType::Upwind_3rd)   ||
             (adv_type == AdvType::Centered_4th) ||
             (adv_type == AdvType::Weno_3)       ||
             (adv_type == AdvType::Weno_3Z)      ||
             (adv_type == AdvType::Weno_3MZQ)    )
            {
                eff_adv_type = AdvType::Centered_2nd;

            } else if ( (adv_type == AdvType::Upwind_5th) ||
                        (adv_type == AdvType::Weno_5)     ||
                        (adv_type == AdvType::Weno_5Z)    )
            {
                eff_adv_type = AdvType::Upwind_3rd;

            } else { // (adv_type == AdvType::Centered_6th)

            eff_adv_type = AdvType::Centered_4th;
        }

    } else {

        eff_adv_type = adv_type;
    }

    return(eff_adv_type);
}
#endif
