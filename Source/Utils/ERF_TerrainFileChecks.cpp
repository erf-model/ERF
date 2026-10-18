/**
 * \file ERF_TerrainFileChecks.cpp
 */
#include "ERF_TerrainFileChecks.H"

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
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
// Each end may stop at most one spacing of the list inside the domain.  The allowance on top is
// 1e-3 of that spacing (coordinates written with a few digits) plus a few ulp of the coordinates
// (single precision at map offsets).
template <typename F>
std::string
list_problem (const char* name, int n, F v, Real lo, Real hi)
{
    // A single value makes the surface constant in this direction: nothing to bracket
    if (n == 1) { return ""; }
    for (int m = 1; m < n; ++m) {
        if (!(v(m) > v(m-1))) {
            std::ostringstream os;
            os << "the " << name << " values must increase strictly, but " << name << "[" << m
               << "] = " << v(m) << " follows " << name << "[" << m-1 << "] = " << v(m-1);
            return os.str();
        }
    }
    const Real ulps = Real(8.) * std::numeric_limits<Real>::epsilon()
                    * amrex::max(std::abs(lo), std::abs(hi));
    const Real d_lo = v(1) - v(0);
    const Real d_hi = v(n-1) - v(n-2);
    const Real gap_lo = v(0) - lo;
    const Real gap_hi = hi - v(n-1);
    if (gap_lo > d_lo + Real(1.e-3)*d_lo + ulps || gap_hi > d_hi + Real(1.e-3)*d_hi + ulps) {
        std::ostringstream os;
        os << "the " << name << " values " << v(0) << " .. " << v(n-1) << " do not cover the domain "
           << lo << " .. " << hi << ": each end may stop at most one table spacing inside the"
           << " domain, since the reader only extends the heights from the edge of the table"
           << " (linearly before its first point, and in the ERF format by holding the last"
           << " value after its last point)";
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
    const int nmin = is_usgs ? 2 : 1;
    if (nx < nmin || ny < nmin) {
        os << "nx = " << nx << " and ny = " << ny << " must both be at least " << nmin
           << (is_usgs ? " in the USGS format" : "");
        return os.str();
    }
    const Long npts = Long(nx) * Long(ny);
    const Long nxv  = is_usgs ? npts : Long(nx);
    if (Long(x.size()) != nxv || Long(y.size()) != Long(ny) || Long(z.size()) != npts) {
        os << "expected " << nxv << " x values, " << ny << " y values and " << npts
           << " heights, found " << x.size() << ", " << y.size() << " and " << z.size();
        return os.str();
    }
    // A file cannot hold nan or inf (the reader stops on them as "not a number"); this is a
    // guard for tables built otherwise
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
    if (sx < Real(0.)) {
        os << "    one x value: constant in x";
    } else {
        os << "    largest |dh/dx| " << sx << " at (" << sx_x << ", " << sx_y << ")";
    }
    if (sy < Real(0.)) {
        os << ", one y value: constant in y\n";
    } else {
        os << ", largest |dh/dy| " << sy << " at (" << sy_x << ", " << sy_y << ")\n";
    }
    return os.str();
}

namespace {

// Read and check the table on the I/O rank (aborts on any problem), then broadcast it
HeightTable
read_and_broadcast (const std::string& fname, bool is_usgs,
                    const GpuArray<Real,AMREX_SPACEDIM>& prob_lo,
                    const GpuArray<Real,AMREX_SPACEDIM>& prob_hi)
{
    HeightTable t;
    if (ParallelDescriptor::IOProcessor()) {
        Print() << "Reading terrain file: " << fname << std::endl;
        std::ifstream file(fname);
        if (!file.is_open()) {
            Abort("Error: Could not open the file " + fname + "\n");
        }
        if (file.peek() == std::ifstream::traits_type::eof()) {
            Abort("Error: The file  " + fname + " is empty.\n");
        }
        // nx and ny must be at least 1 to read the file (%nx, sizes) and nx*ny must fit an int;
        // the layout's own minimum is TableProblem's
        auto check_size = [&] (int nx, int ny, const std::string& where) {
            if (nx < 1 || ny < 1) {
                Abort("Terrain file " + fname + ": nx = " + std::to_string(nx) + " and ny = "
                      + std::to_string(ny) + " (" + where + ") must be positive");
            }
            if (Long(nx) * Long(ny) > Long(std::numeric_limits<int>::max())) {
                Abort("Terrain file " + fname + ": nx*ny is too large");
            }
        };

        if (is_usgs) {
            // Line 1: longitude and latitude of the origin; line 2: nx ny; then one
            // "x y height" per line with x varying fastest.  Read line by line, so that a
            // truncated or extra record is caught (blank lines are skipped).
            std::string line;
            int lineno = 0;
            auto next_line = [&] () -> bool {
                while (std::getline(file, line)) {
                    ++lineno;
                    if (line.find_first_not_of(" \t\r") != std::string::npos) { return true; }
                }
                return false;
            };
            auto rest_is_blank = [] (std::istringstream& ss) {
                std::string extra;
                return !(ss >> extra);
            };
            Real lon_min = 0, lat_min = 0;
            {
                std::istringstream ss;
                if (next_line()) { ss.str(line); }
                if (!(ss >> lon_min >> lat_min) || !rest_is_blank(ss)) {
                    Abort("Terrain file " + fname + ": the first line must give the longitude"
                          " and latitude of the origin");
                }
            }
            if (std::abs(lon_min) > Real(180.0)) {
                Error("The value of longitude for entry in the first line in " + fname
                      + " should not exceed 180. It is " + std::to_string(lon_min));
            }
            if (std::abs(lat_min) > Real(90.0)) {
                Error("The value of latitude for entry in the first line in " + fname
                      + " should not exceed 90. It is " + std::to_string(lat_min));
            }
            {
                std::istringstream ss;
                if (next_line()) { ss.str(line); }
                if (!(ss >> t.nx >> t.ny) || !rest_is_blank(ss)) {
                    Abort("Terrain file " + fname + ": the second line must give nx and ny");
                }
            }
            check_size(t.nx, t.ny, "line 2");
            const int npts = t.nx * t.ny;
            t.x.reserve(npts);
            t.z.reserve(npts);
            int counter = 0;
            while (next_line()) {
                std::istringstream ss(line);
                Real xv, yv, zv;
                if (!(ss >> xv >> yv >> zv) || !rest_is_blank(ss)) {
                    Abort("Terrain file " + fname + ": line " + std::to_string(lineno)
                          + " must hold exactly three numbers, x y height: '" + line + "'");
                }
                if (counter == npts) {
                    Abort("Terrain file " + fname + ": line " + std::to_string(lineno)
                          + " is beyond the nx*ny = " + std::to_string(npts) + " records nx and ny ask for");
                }
                t.x.push_back(xv);
                if (counter % t.nx == 0) {
                    t.y.push_back(yv);
                } else if (yv != t.y.back()) {
                    // Only the first y of each row is used: the rest must repeat it
                    Abort("Terrain file " + fname + ": y changes within row " + std::to_string(counter/t.nx)
                          + " (line " + std::to_string(lineno) + "); each row of nx points must share"
                          + " one y, with x varying fastest");
                }
                t.z.push_back(zv);
                ++counter;
            }
            if (counter != npts) {
                Abort("Terrain file " + fname + ": expected nx*ny = " + std::to_string(npts)
                      + " lines of x y height after nx and ny, read " + std::to_string(counter));
            }
        } else {
            int cnt = 1;
            t.nx = erf_get_single_value<int>(file, cnt, fname); cnt++;
            t.ny = erf_get_single_value<int>(file, cnt, fname); cnt++;
            check_size(t.nx, t.ny, "lines 1 and 2");
            Print() << "Expecting " << t.nx << " values of x,    " << t.ny << " values of y, and "
                    << t.nx*t.ny << " values of z" << std::endl;
            t.x.resize(t.nx);
            t.y.resize(t.ny);
            t.z.resize(t.nx * t.ny);
            for (auto& v : t.x) { v = erf_get_single_value<Real>(file, cnt, fname); cnt++; }
            for (auto& v : t.y) { v = erf_get_single_value<Real>(file, cnt, fname); cnt++; }
            for (auto& v : t.z) { v = erf_get_single_value<Real>(file, cnt, fname); cnt++; }
            if (!(file >> std::ws).eof()) {
                Abort("Terrain file " + fname + ": text after line " + std::to_string(cnt-1)
                      + ", the last of the nx + ny + nx*ny values that nx = "
                      + std::to_string(t.nx) + " and ny = " + std::to_string(t.ny) + " ask for");
            }
        }

        // Structure and coverage of the table, and a summary of it, so that a transposed or
        // mis-ordered file can be spotted at start-up
        const std::string problem = TableProblem(is_usgs, t.nx, t.ny, t.x, t.y, t.z, prob_lo, prob_hi);
        if (!problem.empty()) {
            Abort("Terrain file " + fname + ": " + problem);
        }
        Print() << TableSummary(fname, is_usgs, t.nx, t.ny, t.x, t.y, t.z);
    }

    const int io = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::Bcast(&t.nx, 1, io);
    ParallelDescriptor::Bcast(&t.ny, 1, io);
    const int npts = t.nx * t.ny;
    t.x.resize(is_usgs ? npts : t.nx);
    t.y.resize(t.ny);
    t.z.resize(npts);
    ParallelDescriptor::Bcast(t.x.data(), t.x.size(), io);
    ParallelDescriptor::Bcast(t.y.data(), t.y.size(), io);
    ParallelDescriptor::Bcast(t.z.data(), t.z.size(), io);
    return t;
}

} // namespace

const HeightTable&
ReadHeightTable (const std::string& fname, bool is_usgs,
                 const GpuArray<Real,AMREX_SPACEDIM>& prob_lo,
                 const GpuArray<Real,AMREX_SPACEDIM>& prob_hi)
{
    // One entry per file and layout, the same on every rank (every rank calls with the same
    // name), so the file is read, checked and broadcast once however often it is asked for
    static std::map<std::string, HeightTable> tables;
    const std::string key = (is_usgs ? "USGS:" : "ERF:") + fname;
    auto it = tables.find(key);
    if (it == tables.end()) {
        it = tables.emplace(key, read_and_broadcast(fname, is_usgs, prob_lo, prob_hi)).first;
    }
    return it->second;
}

Array<Real,4>
FaceBandSlopes (const MultiFab& z_phys_nd, const MultiFab& mf_mx, const MultiFab& mf_my,
                const Geometry& geom, int band)
{
    AMREX_ALWAYS_ASSERT(band >= 1);
    const Box& dom = geom.Domain();
    const int klo = dom.smallEnd(2);
    const Real dxinv = geom.InvCellSize(0);
    const Real dyinv = geom.InvCellSize(1);

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
            // the map factors live on the same boxes flattened to one k (ba2d)
            const Array4<const Real> mfx = mf_mx.const_array(mfi);
            const Array4<const Real> mfy = mf_my.const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return {CellTerrainSlope(i, j, k, z_nd, mfx(i,j,0), mfy(i,j,0), dxinv, dyinv)};
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
    os << "    Flattening the terrain within " << band << " cells of inflow and outflow faces avoids"
       << " this.  With an eddy-diffusivity closure, erf.smag2d_slope_limiter or"
       << " erf.smag2d_kh_cap (Smagorinsky2D), erf.implicit_terrain_metric (needs"
       << " erf.vert_implicit_fac > 0) or erf.diffusive_dt_limit limit the terrain-metric"
       << " diffusion on steep cells; none of them is on here\n";
    return os.str();
}

} // namespace erf_terrain_file
