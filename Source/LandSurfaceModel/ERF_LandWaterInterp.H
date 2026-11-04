#ifndef ERF_LAND_WATER_INTERP_H_
#define ERF_LAND_WATER_INTERP_H_

#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_MultiFab.H>

#include <limits>

/**
 * Conservative linear interpolation of a land-surface field that does not mix land and water
 * (erf-model/ERF#4390), in the spirit of WRF's land/water-masked nest interpolation.
 *
 * The coarse data carries two components: the field, and a land/water indicator (one value for
 * water, another for land; NaN outside the domain). In each direction the slope of a coarse cell
 * is taken only if both of its neighbours are of the same type as the cell; otherwise it is zero,
 * as it is at a domain boundary. The min/max limiter likewise looks only at neighbours of the
 * same type. Where every neighbour is of the same type this is exactly AMReX's cell_cons_interp
 * (CellConservativeLinear without linear limiting), so the result changes only next to a
 * coastline. Each coarse cell's mean is kept, and the indicator is copied (piecewise constant).
 */
class LandWaterConsLinear
    : public amrex::Interpolater
{
public:
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
};

/**
 * Slopes of the field (component nu) of coarse cell (i,j,k), from neighbours of its own type only
 * (indicator in component nm). Follows mf_cell_cons_lin_interp_mcslope term for term.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void land_water_mcslope (int i, int j, int k, amrex::Array4<amrex::Real> const& slope,
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

    slope(i,j,k,0) = sx * alpha;
    slope(i,j,k,1) = sy * alpha;
    slope(i,j,k,2) = sz * alpha;
}

inline void
LandWaterConsLinear::interp (const amrex::FArrayBox& crse, int crse_comp,
                             amrex::FArrayBox& fine, int fine_comp, int ncomp,
                             const amrex::Box& fine_region, const amrex::IntVect& ratio,
                             const amrex::Geometry& crse_geom, const amrex::Geometry& /*fine_geom*/,
                             amrex::Vector<amrex::BCRec> const& bcr,
                             int /*actual_comp*/, int /*actual_state*/, amrex::RunOn runon)
{
    using namespace amrex;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncomp == 2,
        "LandWaterConsLinear interpolates one field together with its land/water indicator");
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
    FArrayBox ccfab(cslope_bx, AMREX_SPACEDIM, run_on_gpu ? The_Async_Arena() : The_Cpu_Arena());
    Array4<Real> const& tmp = ccfab.array();
    Array4<Real const> const& ctmp = ccfab.const_array();

    const BCRec bc = bcr[0];
    const int nu = crse_comp;
    const int nm = crse_comp + 1;

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
    {
        land_water_mcslope(i, j, k, tmp, crsearr, nu, nm, cdomain, ratio, bc);
    });

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, fine_region, i, j, k,
    {
        mf_cell_cons_lin_interp(i, j, k, 0, finearr, fine_comp, ctmp, crsearr, nu, 1, ratio);
        finearr(i,j,k,fine_comp+1) = crsearr(amrex::coarsen(i,ratio[0]),
                                             amrex::coarsen(j,ratio[1]),
                                             amrex::coarsen(k,ratio[2]), nm);
    });
}

/**
 * Interpolate the coarse field src (one component) into the fine field dst with
 * LandWaterConsLinear. The land/water type of each coarse cell is component 0 of indicator, one
 * value per type, on any BoxArray (it is copied onto src's cells, ghost cells included). Coarse
 * cells it does not reach, and cells outside the domain, have no type and are never anyone's
 * neighbour. Fine cells the interpolation does not fill (ghost cells outside the domain) keep
 * their values in dst.
 */
inline void
interp_land_water (amrex::MultiFab& dst, const amrex::MultiFab& src,
                   const amrex::MultiFab& indicator,
                   const amrex::Geometry& cgeom, const amrex::Geometry& fgeom,
                   const amrex::IntVect& ratio,
                   const amrex::Vector<amrex::BCRec>& bcs, int bcscomp)
{
    using namespace amrex;
    static LandWaterConsLinear land_water_interp;
    constexpr Real no_type = std::numeric_limits<Real>::quiet_NaN();

    // The field and its indicator, interpolated together
    MultiFab csrc(src.boxArray(), src.DistributionMap(), 2, src.nGrowVect());
    MultiFab::Copy(csrc, src, 0, 0, 1, src.nGrowVect());
    csrc.setVal(no_type, 1, 1, src.nGrowVect());
    csrc.ParallelCopy(indicator, 0, 1, 1, indicator.nGrowVect(), src.nGrowVect(),
                      cgeom.periodicity());
    csrc.setDomainBndry(no_type, 1, 1, cgeom);

    MultiFab fdst(dst.boxArray(), dst.DistributionMap(), 2, dst.nGrowVect());
    MultiFab::Copy(fdst, dst, 0, 0, 1, dst.nGrowVect());
    fdst.setVal(Real(0.0), 1, 1, dst.nGrowVect());
    InterpFromCoarseLevel(fdst, dst.nGrowVect(), IntVect(0), csrc, 0, 0, 2,
                          cgeom, fgeom, ratio, &land_water_interp, bcs, bcscomp);
    MultiFab::Copy(dst, fdst, 0, 0, 1, dst.nGrowVect());
}

#endif
