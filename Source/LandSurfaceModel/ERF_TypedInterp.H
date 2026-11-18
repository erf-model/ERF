#ifndef ERF_TYPED_INTERP_H_
#define ERF_TYPED_INTERP_H_

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_MultiFab.H>

#include <limits>

/**
 * Conservative linear interpolation that does not mix cells of different types, in the spirit
 * of WRF's land/water-masked nest interpolation (erf-model/ERF#4390). Every cell is
 * interpolated; a cell's type decides only which neighbours its slopes and limiter may see.
 * The types are whatever the caller encodes: land and water (SLM's state), a valid value and
 * the lsm_undefined sentinel (Noah-MP), or both at once (SLM's fluxes).
 *
 * The coarse data carries nfield fields followed by nind type indicators (one value per type;
 * NaN, no type, outside the domain). With nind == 1 every field shares the one indicator; with
 * nind == nfield field n has indicator n. In each direction the slope of a coarse cell is taken
 * only if both of its neighbours are of the same type as the cell; otherwise it is zero, as it
 * is at a domain boundary. The min/max limiter likewise looks only at neighbours of the same
 * type. Where every neighbour is of the same type this is exactly AMReX's cell_cons_interp
 * (CellConservativeLinear without linear limiting), so the result changes only next to a type
 * boundary. Each coarse cell's mean is kept, and the indicator is copied (piecewise constant),
 * so a cell whose neighbours are all of other types (an isolated sentinel, say) comes through
 * exactly.
 */
class TypedConsLinear
    : public amrex::Interpolater
{
public:
    explicit TypedConsLinear (int nind = 1) : m_nind(nind) {}

    amrex::Box CoarseBox (const amrex::Box& fine, const amrex::IntVect& ratio) override
    {
        amrex::Box crse = amrex::coarsen(fine, ratio);
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
            if (ratio[dim] > 1) { crse.grow(dim, 1); }
        }
        return crse;
    }

    amrex::Box CoarseBox (const amrex::Box& fine, int ratio) override
    {
        return CoarseBox(fine, amrex::IntVect(ratio));
    }

    void interp (const amrex::FArrayBox& crse, int crse_comp,
                 amrex::FArrayBox& fine, int fine_comp, int ncomp,
                 const amrex::Box& fine_region, const amrex::IntVect& ratio,
                 const amrex::Geometry& crse_geom, const amrex::Geometry& /*fine_geom*/,
                 amrex::Vector<amrex::BCRec> const& bcr,
                 int /*actual_comp*/, int /*actual_state*/, amrex::RunOn runon) override;

private:
    int m_nind;
};

/**
 * Slopes of the field (component nu) of coarse cell (i,j,k), from neighbours of its own type only
 * (indicator in component nm). Follows mf_cell_cons_lin_interp_mcslope term for term; the slope
 * in direction d goes to component ns + d*nslope, the layout mf_cell_cons_lin_interp reads.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void typed_mcslope (int i, int j, int k, amrex::Array4<amrex::Real> const& slope,
                    int ns, int nslope,
                    amrex::Array4<amrex::Real const> const& u, int nu, int nm,
                    amrex::Box const& domain, amrex::IntVect const& ratio,
                    amrex::BCRec const& bc) noexcept
{
    using amrex::Real;
    // NaN (outside the domain) is never of the same type
    const Real m0 = u(i,j,k,nm);
    auto same = [&] (int ii, int jj, int kk) { return u(ii,jj,kk,nm) == m0; };

    Real sx = Real(0.);
    Real sy = Real(0.);
    Real sz = Real(0.);

    if (ratio[0] > 1 && same(i-1,j,k) && same(i+1,j,k)) {
        Real dc = amrex::mf_compute_slopes_x(i, j, k, u, nu, domain, bc);
        Real df = Real(2.0) * (u(i+1,j,k,nu) - u(i  ,j,k,nu));
        Real db = Real(2.0) * (u(i  ,j,k,nu) - u(i-1,j,k,nu));
        sx = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sx = std::copysign(Real(1.),dc)*amrex::min(sx,std::abs(dc));
    }
    if (ratio[1] > 1 && same(i,j-1,k) && same(i,j+1,k)) {
        Real dc = amrex::mf_compute_slopes_y(i, j, k, u, nu, domain, bc);
        Real df = Real(2.0) * (u(i,j+1,k,nu) - u(i,j  ,k,nu));
        Real db = Real(2.0) * (u(i,j  ,k,nu) - u(i,j-1,k,nu));
        sy = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sy = std::copysign(Real(1.),dc)*amrex::min(sy,std::abs(dc));
    }
    if (ratio[2] > 1 && same(i,j,k-1) && same(i,j,k+1)) {
        Real dc = amrex::mf_compute_slopes_z(i, j, k, u, nu, domain, bc);
        Real df = Real(2.0) * (u(i,j,k+1,nu) - u(i,j,k  ,nu));
        Real db = Real(2.0) * (u(i,j,k  ,nu) - u(i,j,k-1,nu));
        sz = (df*db >= Real(0.0)) ? amrex::min(std::abs(df),std::abs(db)) : Real(0.);
        sz = std::copysign(Real(1.),dc)*amrex::min(sz,std::abs(dc));
    }

    Real alpha = Real(1.0);
    if (sx != Real(0.0) || sy != Real(0.0) || sz != Real(0.0)) {
        Real dumax = std::abs(sx) * Real(ratio[0]-1)/Real(2*ratio[0])
            +        std::abs(sy) * Real(ratio[1]-1)/Real(2*ratio[1])
            +        std::abs(sz) * Real(ratio[2]-1)/Real(2*ratio[2]);
        Real umax = u(i,j,k,nu);
        Real umin = u(i,j,k,nu);
        const int ilim = ratio[0] > 1 ? 1 : 0;
        const int jlim = ratio[1] > 1 ? 1 : 0;
        const int klim = ratio[2] > 1 ? 1 : 0;
        for (int koff = -klim; koff <= klim; ++koff) {
        for (int joff = -jlim; joff <= jlim; ++joff) {
        for (int ioff = -ilim; ioff <= ilim; ++ioff) {
            if (same(i+ioff,j+joff,k+koff)) {
                umin = amrex::min(umin, u(i+ioff,j+joff,k+koff,nu));
                umax = amrex::max(umax, u(i+ioff,j+joff,k+koff,nu));
            }
        }}}
        // dumax > 0 here; this only keeps an if-converted divide below from dividing by zero
        dumax = dumax == Real(0.) ? Real(1.) : dumax;
        if (dumax * alpha > (umax - u(i,j,k,nu))) {
            alpha = (umax - u(i,j,k,nu)) / dumax;
        }
        if (dumax * alpha > (u(i,j,k,nu) - umin)) {
            alpha = (u(i,j,k,nu) - umin) / dumax;
        }
    }

    slope(i,j,k,ns         ) = sx * alpha;
    slope(i,j,k,ns+  nslope) = sy * alpha;
    slope(i,j,k,ns+2*nslope) = sz * alpha;
}

inline void
TypedConsLinear::interp (const amrex::FArrayBox& crse, int crse_comp,
                         amrex::FArrayBox& fine, int fine_comp, int ncomp,
                         const amrex::Box& fine_region, const amrex::IntVect& ratio,
                         const amrex::Geometry& crse_geom, const amrex::Geometry& /*fine_geom*/,
                         amrex::Vector<amrex::BCRec> const& bcr,
                         int /*actual_comp*/, int /*actual_state*/, amrex::RunOn runon)
{
    using namespace amrex;
    const int nind   = m_nind;
    const int nfield = ncomp - nind;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfield >= 1 && (nind == 1 || nind == nfield),
        "TypedConsLinear interpolates nfield fields followed by 1 or nfield type indicators");
    AMREX_ASSERT(fine.box().contains(fine_region));
    AMREX_ASSERT(!bcr.empty());

    const Box& cdomain = crse_geom.Domain();
    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();

    Box cslope_bx = CoarseBox(fine_region, ratio);
    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
        if (ratio[dim] > 1) { cslope_bx.grow(dim, -1); }
    }
    const bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    FArrayBox ccfab(cslope_bx, AMREX_SPACEDIM*nfield, run_on_gpu ? The_Async_Arena() : The_Cpu_Arena());
    Array4<Real> const& tmp = ccfab.array();
    Array4<Real const> const& ctmp = ccfab.const_array();

    // Every field shares the boundary condition of the first
    const BCRec bc = bcr[0];
    const int ind_comp = crse_comp + nfield;

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, nfield, i, j, k, n,
    {
        const int nm = ind_comp + (nind == 1 ? 0 : n);
        typed_mcslope(i, j, k, tmp, n, nfield, crsearr, crse_comp+n, nm, cdomain, ratio, bc);
    });

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
    {
        if (n < nfield) {
            mf_cell_cons_lin_interp(i, j, k, n, finearr, fine_comp, ctmp, crsearr, crse_comp,
                                    nfield, ratio);
        } else {
            finearr(i,j,k,fine_comp+n) = crsearr(amrex::coarsen(i,ratio[0]),
                                                 amrex::coarsen(j,ratio[1]),
                                                 amrex::coarsen(k,ratio[2]), crse_comp+n);
        }
    });
}

/**
 * TypedConsLinear on data already packed for it. Components 0 to nfield-1 of csrc are the
 * coarse fields and the nind (1 or nfield) after them their type indicators, one value per
 * type; ghost cells included, since the interpolation reads them. The indicators are set to NaN
 * (no type) outside the domain here. fdst has the same layout: its fields are interpolated, and
 * keep their values in the cells the interpolation does not fill (ghost cells outside the
 * domain); its indicator components are scratch. Every field takes the boundary condition bc.
 *
 * Packing many fields moves them in one ParallelCopy rather than one apiece; a caller that packs
 * them anyway saves interp_by_type's copies by building csrc and fdst directly.
 */
inline void
interp_by_type_packed (amrex::MultiFab& fdst, amrex::MultiFab& csrc, int nind,
                       const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                       const amrex::IntVect& ratio, const amrex::BCRec& bc)
{
    using namespace amrex;
    const int ntot   = csrc.nComp();
    const int nfield = ntot - nind;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fdst.nComp() == ntot,
        "interp_by_type_packed: csrc and fdst have different numbers of components");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nfield >= 1 && (nind == 1 || nind == nfield),
        "interp_by_type_packed: the indicators need 1 component or one per field");
    TypedConsLinear typed_interp(nind);

    csrc.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(), nfield, nind, cgeom);

    // InterpFromCoarseLevel reads a BCRec for every component it interpolates
    const Vector<BCRec> bcs(ntot, bc);
    InterpFromCoarseLevel(fdst, fdst.nGrowVect(), IntVect(0), csrc, 0, 0, ntot,
                          cgeom, fgeom, ratio, &typed_interp, bcs, 0);
}

/**
 * Interpolate the coarse fields src into the fine fields dst (same number of components) with
 * TypedConsLinear. The type of each coarse cell is given by indicator, one value
 * per type, on any BoxArray (it is copied onto src's cells, ghost cells included): with one
 * component it applies to every field, with src.nComp() components component n applies to field
 * n. Coarse cells it does not reach, and cells outside the domain, have no type and are never
 * anyone's neighbour. Fine cells the interpolation does not fill (ghost cells outside the domain)
 * keep their values in dst. Every field takes the boundary condition bcs[bcscomp].
 */
inline void
interp_by_type (amrex::MultiFab& dst, const amrex::MultiFab& src,
                const amrex::MultiFab& indicator,
                const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                const amrex::IntVect& ratio,
                const amrex::Vector<amrex::BCRec>& bcs, int bcscomp)
{
    using namespace amrex;
    const int nfield = src.nComp();
    const int nind   = indicator.nComp();
    const int ntot   = nfield + nind;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dst.nComp() == nfield,
        "interp_by_type: src and dst have different numbers of components");

    // The fields and then their indicators, interpolated together
    MultiFab csrc(src.boxArray(), src.DistributionMap(), ntot, src.nGrowVect());
    MultiFab::Copy(csrc, src, 0, 0, nfield, src.nGrowVect());
    csrc.setVal(std::numeric_limits<Real>::quiet_NaN(), nfield, nind, src.nGrowVect());
    csrc.ParallelCopy(indicator, 0, nfield, nind, indicator.nGrowVect(), src.nGrowVect(),
                      cgeom.periodicity());

    MultiFab fdst(dst.boxArray(), dst.DistributionMap(), ntot, dst.nGrowVect());
    MultiFab::Copy(fdst, dst, 0, 0, nfield, dst.nGrowVect());
    fdst.setVal(Real(0.0), nfield, nind, dst.nGrowVect());

    interp_by_type_packed(fdst, csrc, nind, cgeom, fgeom, ratio, bcs[bcscomp]);
    MultiFab::Copy(dst, fdst, 0, 0, nfield, dst.nGrowVect());
}

#endif
