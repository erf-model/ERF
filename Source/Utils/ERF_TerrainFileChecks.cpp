/**
 * \file ERF_TerrainFileChecks.cpp
 */
#include "ERF_TerrainFileChecks.H"

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <sstream>

using namespace amrex;

namespace erf_terrain_file {

namespace {

// x coordinate of table point (i,j) and the height there, in either layout
Real
x_at (bool is_usgs, int nx, const Vector<Real>& x, int i, int j)
{
    return is_usgs ? x[j*nx + i] : x[i];
}

Real
z_at (bool is_usgs, int nx, int ny, const Vector<Real>& z, int i, int j)
{
    return is_usgs ? z[j*nx + i] : z[i*ny + j];
}

// The problem with a strictly increasing list of n values v(0..n-1) spanning [lo, hi], or "".
// Each end may stop at most one spacing of the list (plus a rounding allowance) inside the domain.
template <typename F>
std::string
list_problem (const char* name, int n, F v, Real lo, Real hi)
{
    for (int m = 1; m < n; ++m) {
        if (!(v(m) > v(m-1))) {
            std::ostringstream os;
            os << "the " << name << " values must increase strictly, but " << name << "[" << m
               << "] = " << v(m) << " follows " << name << "[" << m-1 << "] = " << v(m-1);
            return os.str();
        }
    }
    const Real slack = Real(1.e-6) * (hi - lo);
    const Real gap_lo = v(0) - lo;
    const Real gap_hi = hi - v(n-1);
    if (gap_lo > (v(1) - v(0)) + slack || gap_hi > (v(n-1) - v(n-2)) + slack) {
        std::ostringstream os;
        os << "the " << name << " values " << v(0) << " .. " << v(n-1) << " do not cover the domain "
           << lo << " .. " << hi << ": each end may stop at most one table spacing inside the"
           << " domain, since heights beyond the table would be held at its edge value";
        return os.str();
    }
    return "";
}

} // namespace

std::string
TableProblem (bool is_usgs, int nx, int ny,
              const Vector<Real>& x, const Vector<Real>& y, const Vector<Real>& z,
              const GpuArray<Real,AMREX_SPACEDIM>& prob_lo,
              const GpuArray<Real,AMREX_SPACEDIM>& prob_hi)
{
    std::ostringstream os;
    if (nx < 2 || ny < 2) {
        os << "nx = " << nx << " and ny = " << ny << " must both be at least 2";
        return os.str();
    }
    const Long npts = Long(nx) * Long(ny);
    const Long nxv  = is_usgs ? npts : Long(nx);
    if (Long(x.size()) != nxv || Long(y.size()) != Long(ny) || Long(z.size()) != npts) {
        os << "expected " << nxv << " x values, " << ny << " y values and " << npts
           << " heights, found " << x.size() << ", " << y.size() << " and " << z.size();
        return os.str();
    }
    auto finite = [] (const Vector<Real>& v, const char* name) -> std::string {
        for (Long n = 0; n < Long(v.size()); ++n) {
            if (!std::isfinite(v[n])) {
                std::ostringstream o;
                o << name << " value " << n << " is not finite (" << v[n] << ")";
                return o.str();
            }
        }
        return "";
    };
    for (const auto& s : {finite(x,"x"), finite(y,"y"), finite(z,"height")}) {
        if (!s.empty()) { return s; }
    }

    std::string s = list_problem("y", ny, [&] (int m) { return y[m]; }, prob_lo[1], prob_hi[1]);
    if (!s.empty()) { return s; }
    if (is_usgs) {
        for (int j = 0; j < ny; ++j) {
            s = list_problem("x", nx, [&] (int m) { return x[j*nx + m]; }, prob_lo[0], prob_hi[0]);
            if (!s.empty()) {
                return s + " (row " + std::to_string(j) + " of the x-fastest USGS layout)";
            }
        }
    } else {
        s = list_problem("x", nx, [&] (int m) { return x[m]; }, prob_lo[0], prob_hi[0]);
        if (!s.empty()) { return s; }
    }
    return "";
}

std::string
TableSummary (const std::string& fname, bool is_usgs, int nx, int ny,
              const Vector<Real>& x, const Vector<Real>& y, const Vector<Real>& z)
{
    int imin = 0, jmin = 0, imax = 0, jmax = 0;
    Real zmin = z_at(is_usgs, nx, ny, z, 0, 0), zmax = zmin;
    Real sx = Real(-1.), sy = Real(-1.);
    Real sx_x = 0, sx_y = 0, sy_x = 0, sy_y = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const Real h = z_at(is_usgs, nx, ny, z, i, j);
            if (h < zmin) { zmin = h; imin = i; jmin = j; }
            if (h > zmax) { zmax = h; imax = i; jmax = j; }
            if (i+1 < nx) {
                const Real x0 = x_at(is_usgs, nx, x, i, j), x1 = x_at(is_usgs, nx, x, i+1, j);
                const Real s = std::abs(z_at(is_usgs, nx, ny, z, i+1, j) - h) / (x1 - x0);
                if (s > sx) { sx = s; sx_x = Real(0.5)*(x0 + x1); sx_y = y[j]; }
            }
            if (j+1 < ny) {
                const Real s = std::abs(z_at(is_usgs, nx, ny, z, i, j+1) - h) / (y[j+1] - y[j]);
                if (s > sy) { sy = s; sy_x = x_at(is_usgs, nx, x, i, j); sy_y = Real(0.5)*(y[j] + y[j+1]); }
            }
        }
    }

    std::ostringstream os;
    os << "Height file " << fname << ": " << nx << " x " << ny << " heights over x "
       << x_at(is_usgs, nx, x, 0, 0) << " .. " << x_at(is_usgs, nx, x, nx-1, 0)
       << ", y " << y[0] << " .. " << y[ny-1] << "\n";
    if (is_usgs) {
        os << "    read in USGS order: x varies fastest, height n is at (x[n], y[n / nx])\n";
    } else {
        os << "    read in ERF order: y varies fastest, height n is at (x[n / ny], y[n % ny])\n";
    }
    os << "    lowest  " << zmin << " m at (x, y) = (" << x_at(is_usgs, nx, x, imin, jmin)
       << ", " << y[jmin] << ")\n";
    os << "    highest " << zmax << " m at (x, y) = (" << x_at(is_usgs, nx, x, imax, jmax)
       << ", " << y[jmax] << ")\n";
    os << "    largest |dh/dx| " << sx << " at (" << sx_x << ", " << sx_y << "), largest |dh/dy| "
       << sy << " at (" << sy_x << ", " << sy_y << ")\n";
    return os.str();
}

Array<Real,4>
FaceBandSlopes (const MultiFab& z_phys_nd, const Geometry& geom, int band)
{
    AMREX_ALWAYS_ASSERT(band >= 1);
    const Box& dom = geom.Domain();
    const int klo = dom.smallEnd(2);
    const Real dx = geom.CellSize(0);
    const Real dy = geom.CellSize(1);

    // The bands of cells next to each lateral face, at the surface only
    Array<Box,4> bands;
    for (int f = 0; f < 4; ++f) {
        Box b(dom);
        b.setRange(2, klo);
        const int dir = f / 2;
        if (f % 2 == 0) {
            b.setBig(dir, amrex::min(dom.smallEnd(dir) + band - 1, dom.bigEnd(dir)));
        } else {
            b.setSmall(dir, amrex::max(dom.bigEnd(dir) - band + 1, dom.smallEnd(dir)));
        }
        bands[f] = b;
    }

    Array<Real,4> smax{Real(-1.), Real(-1.), Real(-1.), Real(-1.)};
    for (int f = 0; f < 4; ++f) {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        bool any = false;
        for (MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
            // z_phys_nd is nodal: its valid cells are the enclosed cells of the nodal box
            const Box bx = enclosedCells(mfi.validbox()) & bands[f];
            if (bx.isEmpty()) { continue; }
            any = true;
            const Array4<const Real> z_nd = z_phys_nd.const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return {SurfaceSlope(i, j, k, z_nd, dx, dy)};
            });
        }
        if (any) { smax[f] = amrex::get<0>(reduce_data.value(reduce_op)); }
    }
    ParallelDescriptor::ReduceRealMax(smax.data(), 4);
    return smax;
}

std::string
OpenFaceSlopeWarning (int lev, int band,
                      const Array<Real,4>& slopes,
                      const Array<bool,4>& open,
                      const Array<std::string,4>& bc_names)
{
    static const char* face_names[4] = {"xlo", "xhi", "ylo", "yhi"};
    std::ostringstream os;
    for (int f = 0; f < 4; ++f) {
        if (open[f] && slopes[f] > open_face_slope_warn) {
            os << "WARNING: at level " << lev << " the terrain slope |grad h| reaches " << slopes[f]
               << " within " << band << " cells of the " << face_names[f] << " face ("
               << bc_names[f] << "), above " << open_face_slope_warn << "\n";
        }
    }
    if (os.str().empty()) { return ""; }
    os << "    No steep-terrain option is on (erf.smag2d_slope_limiter, erf.smag2d_kh_cap,"
       << " erf.implicit_terrain_metric, erf.diffusive_dt_limit); consider one of them, or"
       << " flatten the terrain next to inflow and outflow faces\n";
    return os.str();
}

} // namespace erf_terrain_file
