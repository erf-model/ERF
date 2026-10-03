/**
 * \file ERF_MLTerrainPoisson_K.H
 *
 * Device helpers for MLTerrainPoisson, the multigrid form of the terrain-fitted
 * anelastic Poisson operator.  The stencil itself lives in ERF_TerrainPoisson_3D_K.H
 * and is shared with the GMRES operator; this file adds the pieces multigrid needs:
 * the accessor that lets the MLMG vectors (no ghost cells in a hidden direction)
 * feed the stencil, the model of the ghost-cell fill that the smoother uses to read
 * off the exact stencil coefficients, and the column (z-line) relaxation.
 */
#ifndef ERF_MLTERRAINPOISSON_K_H_
#define ERF_MLTERRAINPOISSON_K_H_

#include <AMReX_Array4.H>
#include <AMReX_Dim3.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_TerrainPoisson_3D_K.H"

/**
 * Cell-centered accessor for the MLMG vectors.
 *
 * In a hidden direction (a direction one cell wide, see LPInfo::setHiddenDirection)
 * the MLMG vectors have no ghost cells, and the +1/-1 neighbours the stencil asks for
 * are the cell itself; the accessor folds them back onto it.  Everywhere else it is a
 * plain Array4 read.
 */
struct TerrainStencilAccess
{
    amrex::Array4<amrex::Real const> a;
    int hidden_dir;   //!< -1 for none, 0 for x, 1 for y
    int hidden_idx;   //!< the one index value in the hidden direction

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k) const noexcept
    {
        if (hidden_dir == 0) { i = hidden_idx; }
        else if (hidden_dir == 1) { j = hidden_idx; }
        return a(i,j,k);
    }
};

/**
 * Sign a ghost value carries relative to the cell it mirrors.
 *
 * +1 for Neumann (ghost = interior), -1 for Dirichlet (ghost = -interior, the
 * second-order homogeneous fill), chosen per face of the solve domain; a box face
 * that lies inside the domain (coarse/fine or box-union boundary) uses cf_sign.
 */
struct TerrainGhostRule
{
    amrex::Dim3 dlo;                  //!< domain cell lower bounds
    amrex::Dim3 dhi;                  //!< domain cell upper bounds
    amrex::GpuArray<int,3> lo_sign;   //!< sign beyond the domain low face, per direction
    amrex::GpuArray<int,3> hi_sign;   //!< sign beyond the domain high face, per direction
    int cf_sign;                      //!< sign at a box face inside the domain

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int sign (int dir, int idx) const noexcept
    {
        int dl = (dir == 0) ? dlo.x : ((dir == 1) ? dlo.y : dlo.z);
        int dh = (dir == 0) ? dhi.x : ((dir == 1) ? dhi.y : dhi.z);
        if (idx < dl) { return lo_sign[dir]; }
        if (idx > dh) { return hi_sign[dir]; }
        return cf_sign;
    }
};

/**
 * Resolve where an uncovered ghost cell takes its value from.
 *
 * MLTerrainPoisson::applyBC fills every ghost cell that no box of the level covers
 * by mirroring it across the boundary it lies beyond: a domain face (sign from the
 * domain boundary condition) or, inside the domain, a face of the box union
 * (coarse/fine sign).  A cell beyond two boundaries is mirrored across both.  The
 * mirror source is always a valid or covered cell, so the fill is a single gather
 * once the sources are known.  This function follows the chain of mirrors for one
 * ghost cell and returns the source cell and the product of the signs.
 *
 * @param i, j, k the ghost cell
 * @param blo, bhi valid box bounds
 * @param mask 1 on uncovered ghost cells
 * @param rule signs per domain face and for interior (coarse/fine) faces
 * @param si, sj, sk the source cell (output)
 * @return sign the ghost value carries relative to the source
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real terrain_ghost_source (int i, int j, int k,
                                  amrex::Dim3 const& blo, amrex::Dim3 const& bhi,
                                  amrex::Array4<int const> const& mask,
                                  TerrainGhostRule const& rule,
                                  int& si, int& sj, int& sk) noexcept
{
    amrex::Real sgn = amrex::Real(1.0);
    for (int step = 0; step < 3; ++step)
    {
        bool in_i = (i >= blo.x && i <= bhi.x);
        bool in_j = (j >= blo.y && j <= bhi.y);
        bool in_k = (k >= blo.z && k <= bhi.z);
        if ((in_i && in_j && in_k) || mask(i,j,k) == 0) { break; }

        int d = -1;
        // A domain face first: the boundary condition decides
        if      (i < rule.dlo.x || i > rule.dhi.x) { d = 0; }
        else if (j < rule.dlo.y || j > rule.dhi.y) { d = 1; }
        else if (k < rule.dlo.z || k > rule.dhi.z) { d = 2; }
        else {
            // Inside the domain: the face of the box union, i.e. the direction whose
            // inward neighbour is valid or covered
            int ii = in_i ? i : ((i < blo.x) ? i+1 : i-1);
            int jj = in_j ? j : ((j < blo.y) ? j+1 : j-1);
            int kk = in_k ? k : ((k < blo.z) ? k+1 : k-1);
            if      (!in_i && mask(ii,j,k) == 0) { d = 0; }
            else if (!in_j && mask(i,jj,k) == 0) { d = 1; }
            else if (!in_k && mask(i,j,kk) == 0) { d = 2; }
            else if (!in_i) { d = 0; }   // beyond two union faces: mirror across both in turn
            else if (!in_j) { d = 1; }
            else            { d = 2; }
        }
        if (d == 0) { sgn *= static_cast<amrex::Real>(rule.sign(0,i)); i = (i < blo.x) ? i+1 : i-1; }
        if (d == 1) { sgn *= static_cast<amrex::Real>(rule.sign(1,j)); j = (j < blo.y) ? j+1 : j-1; }
        if (d == 2) { sgn *= static_cast<amrex::Real>(rule.sign(2,k)); k = (k < blo.z) ? k+1 : k-1; }
    }
    si = i; sj = j; sk = k;
    return sgn;
}

/**
 * The unit vector e_(i0,j0,k0) as the stencil sees it after the ghost-cell fill.
 *
 * A covered ghost cell holds the neighbouring box's data (or the periodic image);
 * an uncovered one holds sgn times its source cell (see terrain_ghost_source), the
 * tables MLTerrainPoisson::applyBC itself gathers from.  Reading this functor through
 * terrpoisson_adotx_value at cell (i,j,k) gives the exact matrix entry
 * A[(i,j,k),(i0,j0,k0)] including every boundary fold, which is how the smoother's
 * tridiagonal coefficients are obtained without re-deriving the stencil.
 */
struct TerrainUnitVector
{
    int i0, j0, k0;                      //!< the unknown whose column is probed
    amrex::Dim3 blo;                     //!< valid box lower bounds
    amrex::Dim3 bhi;                     //!< valid box upper bounds
    amrex::Array4<int const> mask;       //!< 1 on uncovered ghost cells, 0 elsewhere
    amrex::Array4<int const> src;        //!< source cell of an uncovered ghost cell (3 components)
    amrex::Array4<amrex::Real const> sgn;//!< sign of an uncovered ghost cell relative to its source
    amrex::Dim3 dlo;                     //!< domain lower bounds
    amrex::Dim3 dhi;                     //!< domain upper bounds
    amrex::Dim3 nper;                    //!< domain length in a periodic direction, 0 otherwise
    int hidden_dir;
    int hidden_idx;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k) const noexcept
    {
        if (hidden_dir == 0) { i = hidden_idx; }
        else if (hidden_dir == 1) { j = hidden_idx; }

        amrex::Real s = amrex::Real(1.0);
        bool inside = (i >= blo.x && i <= bhi.x && j >= blo.y && j <= bhi.y && k >= blo.z && k <= bhi.z);
        if (!inside && mask(i,j,k) == 1) {
            s = sgn(i,j,k);
            int ii = src(i,j,k,0);
            int jj = src(i,j,k,1);
            int kk = src(i,j,k,2);
            i = ii; j = jj; k = kk;
            inside = (i >= blo.x && i <= bhi.x && j >= blo.y && j <= bhi.y && k >= blo.z && k <= bhi.z);
        }
        if (!inside) {
            // Covered ghost cell: another box's data or a periodic image
            if (nper.x > 0) { if (i < dlo.x) { i += nper.x; } else if (i > dhi.x) { i -= nper.x; } }
            if (nper.y > 0) { if (j < dlo.y) { j += nper.y; } else if (j > dhi.y) { j -= nper.y; } }
            if (nper.z > 0) { if (k < dlo.z) { k += nper.z; } else if (k > dhi.z) { k -= nper.z; } }
        }
        return (i == i0 && j == j0 && k == k0) ? s : amrex::Real(0.0);
    }
};

/**
 * Whether a ghost cell lies across a coarse/fine (box-union) face: it is uncovered and
 * not beyond a non-periodic domain face.  The correction flux through such a face is
 * zero, as through the surface, because the momentum there is fixed by the coarse level.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool terrain_is_cf_cell (int i, int j, int k,
                         amrex::Array4<int const> const& mask,
                         amrex::Dim3 const& dlo, amrex::Dim3 const& dhi,
                         amrex::Dim3 const& nper) noexcept
{
    if (mask(i,j,k) != 1) { return false; }
    bool beyond = (nper.x == 0 && (i < dlo.x || i > dhi.x))
               || (nper.y == 0 && (j < dlo.y || j > dhi.y))
               || (nper.z == 0 && (k < dlo.z || k > dhi.z));
    return !beyond;
}

/**
 * The terrain Poisson operator at a cell of a box that touches a coarse/fine face:
 * assembled from the face fluxes, with the flux through every coarse/fine face set to
 * zero.  A Neumann mirror alone leaves the cross-term part of the flux on sloping
 * terrain, which makes an all-Neumann problem on a box union inconsistent.  Away from
 * such faces this equals terrpoisson_adotx_value up to round-off; boxes without
 * coarse/fine faces use that kernel so the finest level stays bitwise the GMRES operator.
 */
template <typename XT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real terrain_adotx_cf (int i, int j, int k, XT const& x,
                              amrex::Array4<int const> const& mask,
                              amrex::Dim3 const& dlo, amrex::Dim3 const& dhi, amrex::Dim3 const& nper,
                              int hidden_dir,
                              amrex::Array4<amrex::Real const> const& ax,
                              amrex::Array4<amrex::Real const> const& ay,
                              amrex::Array4<amrex::Real const> const& az,
                              amrex::Array4<amrex::Real const> const& dJ,
                              amrex::Array4<amrex::Real const> const& zp,
                              amrex::Real dxinv, amrex::Real dyinv, amrex::Real dzinv) noexcept
{
    using amrex::Real;
    Real px_lo = -terrpoisson_flux_x(i  , j, k, x, zp, dxinv);
    Real px_hi = -terrpoisson_flux_x(i+1, j, k, x, zp, dxinv);
    Real py_lo = -terrpoisson_flux_y(i, j  , k, x, zp, dyinv);
    Real py_hi = -terrpoisson_flux_y(i, j+1, k, x, zp, dyinv);
    Real pz_lo = -terrpoisson_flux_z(i, j, k  , x, zp, dxinv, dyinv);
    Real pz_hi = -terrpoisson_flux_z(i, j, k+1, x, zp, dxinv, dyinv);
    if (hidden_dir != 0) {
        if (terrain_is_cf_cell(i-1,j,k, mask, dlo, dhi, nper)) { px_lo = Real(0.0); }
        if (terrain_is_cf_cell(i+1,j,k, mask, dlo, dhi, nper)) { px_hi = Real(0.0); }
    }
    if (hidden_dir != 1) {
        if (terrain_is_cf_cell(i,j-1,k, mask, dlo, dhi, nper)) { py_lo = Real(0.0); }
        if (terrain_is_cf_cell(i,j+1,k, mask, dlo, dhi, nper)) { py_hi = Real(0.0); }
    }
    if (terrain_is_cf_cell(i,j,k-1, mask, dlo, dhi, nper)) { pz_lo = Real(0.0); }
    if (terrain_is_cf_cell(i,j,k+1, mask, dlo, dhi, nper)) { pz_hi = Real(0.0); }
    Real y =  (ax(i+1,j,k)*px_hi - ax(i,j,k)*px_lo) * dxinv
             +(ay(i,j+1,k)*py_hi - ay(i,j,k)*py_lo) * dyinv
             +(az(i,j,k+1)*pz_hi - az(i,j,k)*pz_lo) * dzinv;
    return y / dJ(i,j,k);
}

/**
 * Relax one column segment (i,j,klo..khi) of the terrain Poisson operator.
 *
 * Computes the residual of the current iterate in the column, solves the
 * tridiagonal system formed by the column's own entries (sub-, main and
 * super-diagonal in tri components 0, 1, 2) for the correction, and adds it to
 * sol.  Entries that couple the column to its neighbours stay on the right-hand
 * side, so with the (i+j) red-black ordering used by MLTerrainPoisson::Fsmooth
 * this is an exact block Gauss-Seidel sweep: the stencil never couples cells
 * whose (i,j) differ in both directions.
 *
 * @tparam XT Accessor type for the current iterate
 * @param i x-index of the column
 * @param j y-index of the column
 * @param klo first cell of the segment
 * @param khi last cell of the segment
 * @param sol iterate, updated in place on (i,j,klo..khi)
 * @param x read accessor of sol (ghost cells filled)
 * @param rhs right-hand side
 * @param tri tridiagonal coefficients (3 components)
 * @param ax, ay, az, dJ, zp metric terms of this multigrid level
 * @param dxinv, dyinv, dzinv inverse grid spacings of this multigrid level
 * @param cp scratch array, overwritten on the column
 * @param dp scratch array, overwritten on the column
 * @param has_cf whether the box touches a coarse/fine face (then the zero-flux form of the operator is used)
 * @param mask 1 on uncovered ghost cells
 * @param dlo, dhi domain bounds
 * @param nper domain length in each periodic direction, 0 otherwise
 * @param hidden_dir hidden direction (-1 for none)
 */
template <typename XT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void terrain_column_relax (int i, int j, int klo, int khi,
                           amrex::Array4<amrex::Real> const& sol,
                           XT const& x,
                           amrex::Array4<amrex::Real const> const& rhs,
                           amrex::Array4<amrex::Real const> const& tri,
                           amrex::Array4<amrex::Real const> const& ax,
                           amrex::Array4<amrex::Real const> const& ay,
                           amrex::Array4<amrex::Real const> const& az,
                           amrex::Array4<amrex::Real const> const& dJ,
                           amrex::Array4<amrex::Real const> const& zp,
                           amrex::Real dxinv, amrex::Real dyinv, amrex::Real dzinv,
                           amrex::Array4<amrex::Real> const& cp,
                           amrex::Array4<amrex::Real> const& dp,
                           bool has_cf,
                           amrex::Array4<int const> const& mask,
                           amrex::Dim3 const& dlo, amrex::Dim3 const& dhi, amrex::Dim3 const& nper,
                           int hidden_dir) noexcept
{
    using amrex::Real;

    // Residual of the current iterate, stored in dp
    for (int k = klo; k <= khi; ++k) {
        Real ax_ = has_cf ? terrain_adotx_cf(i, j, k, x, mask, dlo, dhi, nper, hidden_dir, ax, ay, az, dJ, zp, dxinv, dyinv, dzinv)
                          : terrpoisson_adotx_value(i, j, k, x, ax, ay, az, dJ, zp, dxinv, dyinv, dzinv);
        dp(i,j,k) = rhs(i,j,k) - ax_;
    }

    // Thomas algorithm on the column's tridiagonal block
    {
        Real b = tri(i,j,klo,1);
        cp(i,j,klo) = tri(i,j,klo,2) / b;
        dp(i,j,klo) = dp(i,j,klo) / b;
    }
    for (int k = klo+1; k <= khi; ++k) {
        Real a = tri(i,j,k,0);
        Real m = tri(i,j,k,1) - a * cp(i,j,k-1);
        cp(i,j,k) = tri(i,j,k,2) / m;
        dp(i,j,k) = (dp(i,j,k) - a * dp(i,j,k-1)) / m;
    }
    Real delta = dp(i,j,khi);
    sol(i,j,khi) += delta;
    for (int k = khi-1; k >= klo; --k) {
        delta = dp(i,j,k) - cp(i,j,k) * delta;
        sol(i,j,k) += delta;
    }
}

#endif
