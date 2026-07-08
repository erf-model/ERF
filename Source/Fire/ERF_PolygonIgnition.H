#ifndef ERF_POLYGON_IGNITION_H
#define ERF_POLYGON_IGNITION_H

/**
 * @file ERF_PolygonIgnition.H
 * @brief Polygon and polyline ignition initialization
 *
 * Polygon mode: interior set to phi = -dist_to_boundary (burned); exterior set
 * to phi = +dist_to_boundary (unburned). Boundary is determined by the winding
 * number test (point is inside if winding number != 0).
 *
 * Polyline mode: cells within half_width set to phi = -half_width (burned);
 * others set to distance from nearest segment. A polyline is an open path of
 * line segments (not closed).
 *
 * Both modes use min(phi, new_val) stamping, so ignitions merge with existing
 * fire state without overwriting it.
 *
 * Reference: Finney (1998) RMRS-RP-4: "An ignition polygon defines the initial
 * fire perimeter"
 */

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>

/**
 * @brief Read polygon vertices from a CSV file.
 *
 * File format:
 *   x  y
 *   (one vertex per line; # or ! lines are comments; blank lines ignored)
 *
 * File is read only on rank 0 (IOProcessor). Vertex arrays are broadcast to
 * all ranks via MPI_Bcast before returning.
 *
 * @param[in]  filename  Path to vertex CSV file
 * @param[out] xs        x-coordinates of vertices
 * @param[out] ys        y-coordinates of vertices
 */
inline
void read_polygon_vertices(const std::string& filename,
                           std::vector<amrex::Real>& xs,
                           std::vector<amrex::Real>& ys)
{
    using namespace amrex;

    xs.clear();
    ys.clear();

    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            amrex::Abort("Cannot open polygon vertex file: " + filename);
        }

        std::string line;
        int line_num = 0;

        while (std::getline(file, line)) {
            ++line_num;

            // Skip empty lines
            if (line.empty()) continue;

            // Skip comment lines
            if (line[0] == '#' || line[0] == '!') continue;

            // Convert commas to spaces
            for (auto& c : line) {
                if (c == ',') c = ' ';
            }

            // Parse line
            std::istringstream iss(line);
            Real x, y;

            if (!(iss >> x >> y)) {
                amrex::Print() << "[WARNING] Line " << line_num
                               << " in " << filename
                               << " could not be parsed; skipping\n";
                continue;
            }

            xs.push_back(x);
            ys.push_back(y);
        }

        file.close();

        if (xs.empty()) {
            amrex::Abort("Polygon vertex file '" + filename + "' is empty or contains only comments");
        }
    }

    // Broadcast vertex count
    int n_verts = static_cast<int>(xs.size());
    if (ParallelDescriptor::NProcs() > 1) {
        MPI_Bcast(&n_verts, 1, MPI_INT, ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);

        if (n_verts > 0) {
            if (!ParallelDescriptor::IOProcessor()) {
                xs.resize(n_verts);
                ys.resize(n_verts);
            }
            MPI_Bcast(xs.data(), n_verts, MPI_DOUBLE, ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
            MPI_Bcast(ys.data(), n_verts, MPI_DOUBLE, ParallelDescriptor::IOProcessorNumber(), MPI_COMM_WORLD);
        }
    }
}

/**
 * @brief Compute the winding number of a point relative to a polygon.
 *
 * Uses the ray-casting algorithm: count how many times a horizontal ray from
 * the point to infinity crosses the polygon boundary. Odd crossings mean the
 * point is inside.
 *
 * @param[in] px, py     Query point
 * @param[in] xs, ys     Polygon vertex arrays (device pointers)
 * @param[in] n_verts    Number of vertices
 * @return Winding number: 0 = outside, != 0 = inside
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int polygon_winding_number(amrex::Real px, amrex::Real py,
                           const amrex::Real* xs, const amrex::Real* ys,
                           int n_verts)
{
    int winding = 0;

    for (int i = 0; i < n_verts; ++i) {
        int j = (i + 1) % n_verts;  // Next vertex (wrap around)
        amrex::Real x1 = xs[i];
        amrex::Real y1 = ys[i];
        amrex::Real x2 = xs[j];
        amrex::Real y2 = ys[j];

        if (y1 <= py) {
            if (y2 > py) {
                // Upward crossing: check if crossing is to the right of point
                amrex::Real cross = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1);
                if (cross > 0.0) {
                    ++winding;
                }
            }
        } else {
            if (y2 <= py) {
                // Downward crossing: check if crossing is to the left of point
                amrex::Real cross = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1);
                if (cross < 0.0) {
                    --winding;
                }
            }
        }
    }

    return winding;
}

/**
 * @brief Compute squared distance from a point to a line segment.
 *
 * Returns the squared distance from point P to the closest point on segment AB.
 *
 * @param[in] px, py     Query point
 * @param[in] ax, ay     Segment start point
 * @param[in] bx, by     Segment end point
 * @return Squared distance
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real segment_dist_sq(amrex::Real px, amrex::Real py,
                            amrex::Real ax, amrex::Real ay,
                            amrex::Real bx, amrex::Real by)
{
    amrex::Real dx = bx - ax;
    amrex::Real dy = by - ay;
    amrex::Real len_sq = dx * dx + dy * dy;

    if (len_sq < 1.0e-14) {
        // Segment is degenerate (point)
        amrex::Real dpx = px - ax;
        amrex::Real dpy = py - ay;
        return dpx * dpx + dpy * dpy;
    }

    // Parameter t of closest point on segment
    amrex::Real t = ((px - ax) * dx + (py - ay) * dy) / len_sq;
    t = amrex::max(0.0, amrex::min(1.0, t));

    amrex::Real cx = ax + t * dx;
    amrex::Real cy = ay + t * dy;
    amrex::Real dpx = px - cx;
    amrex::Real dpy = py - cy;
    return dpx * dpx + dpy * dpy;
}

/**
 * @brief Initialize level-set field from a closed polygon.
 *
 * Sets phi to a signed distance function:
 * - phi < 0 (burned) inside the polygon
 * - phi > 0 (unburned) outside the polygon
 * - Distance to nearest boundary at magnitude |phi|
 *
 * Interior is determined by winding number != 0. Uses min(phi, new_val)
 * semantics to merge with existing fire state.
 *
 * @param[in,out] phi       Level-set field (on fire grid)
 * @param[in]     geom_fire Fire grid geometry
 * @param[in]     xs        x-coordinates of polygon vertices
 * @param[in]     ys        y-coordinates of polygon vertices
 */
inline
void init_phi_from_polygon(
    amrex::MultiFab& phi,
    const amrex::Geometry& geom_fire,
    const std::vector<amrex::Real>& xs,
    const std::vector<amrex::Real>& ys)
{
    using namespace amrex;

    if (xs.empty() || xs.size() != ys.size()) {
        amrex::Abort("Polygon vertices: xs and ys must have same non-zero size");
    }

    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();
    int n_verts = static_cast<int>(xs.size());

    // Copy vertices to GPU if necessary
    Gpu::DeviceVector<Real> d_xs(xs.size());
    Gpu::DeviceVector<Real> d_ys(ys.size());
    Gpu::copy(Gpu::hostToDevice, xs.data(), xs.data() + xs.size(), d_xs.data());
    Gpu::copy(Gpu::hostToDevice, ys.data(), ys.data() + ys.size(), d_ys.data());

    auto* p_xs = d_xs.data();
    auto* p_ys = d_ys.data();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];

            // Compute winding number
            int winding = polygon_winding_number(x, y, p_xs, p_ys, n_verts);
            bool inside = (winding != 0);

            // Compute distance to nearest boundary segment
            Real min_dist = 1.0e10;
            for (int seg = 0; seg < n_verts; ++seg) {
                int seg_next = (seg + 1) % n_verts;
                Real dist_sq = segment_dist_sq(x, y,
                                               p_xs[seg], p_ys[seg],
                                               p_xs[seg_next], p_ys[seg_next]);
                if (dist_sq < min_dist * min_dist) {
                    min_dist = std::sqrt(dist_sq);
                }
            }

            // Set phi: negative inside, positive outside
            Real new_phi = inside ? -min_dist : +min_dist;

            // Apply min() semantics to merge with existing fire
            p(i, j, k) = amrex::min(p(i, j, k), new_phi);
        });
    }
}

/**
 * @brief Initialize level-set field from a polyline (open line fire).
 *
 * Sets phi for cells near a polyline segment sequence:
 * - phi = -half_width for cells within half_width of any segment
 * - phi = distance from segment for cells farther away
 * - Uses min(phi, new_val) semantics to merge with existing fire
 *
 * A polyline is an open path (not closed), defined by a sequence of vertices
 * connected by line segments.
 *
 * @param[in,out] phi        Level-set field (on fire grid)
 * @param[in]     geom_fire  Fire grid geometry
 * @param[in]     xs         x-coordinates of polyline vertices
 * @param[in]     ys         y-coordinates of polyline vertices
 * @param[in]     half_width Half-width of burned zone around polyline [m]
 */
inline
void init_phi_from_polyline(
    amrex::MultiFab& phi,
    const amrex::Geometry& geom_fire,
    const std::vector<amrex::Real>& xs,
    const std::vector<amrex::Real>& ys,
    amrex::Real half_width)
{
    using namespace amrex;

    if (xs.empty() || xs.size() != ys.size()) {
        amrex::Abort("Polyline vertices: xs and ys must have same non-zero size");
    }

    const auto* prob_lo = geom_fire.ProbLo();
    const auto* dx      = geom_fire.CellSize();
    int n_verts = static_cast<int>(xs.size());

    // Copy vertices to GPU if necessary
    Gpu::DeviceVector<Real> d_xs(xs.size());
    Gpu::DeviceVector<Real> d_ys(ys.size());
    Gpu::copy(Gpu::hostToDevice, xs.data(), xs.data() + xs.size(), d_xs.data());
    Gpu::copy(Gpu::hostToDevice, ys.data(), ys.data() + ys.size(), d_ys.data());

    auto* p_xs = d_xs.data();
    auto* p_ys = d_ys.data();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = prob_lo[0] + (i + 0.5) * dx[0];
            Real y = prob_lo[1] + (j + 0.5) * dx[1];

            // Compute distance to nearest segment
            Real min_dist = 1.0e10;
            for (int seg = 0; seg < n_verts - 1; ++seg) {
                Real dist_sq = segment_dist_sq(x, y,
                                               p_xs[seg], p_ys[seg],
                                               p_xs[seg + 1], p_ys[seg + 1]);
                Real dist = std::sqrt(dist_sq);
                if (dist < min_dist) {
                    min_dist = dist;
                }
            }

            // Set phi: negative within half_width, positive otherwise
            Real new_phi = (min_dist <= half_width) ? -half_width : min_dist;

            // Apply min() semantics to merge with existing fire
            p(i, j, k) = amrex::min(p(i, j, k), new_phi);
        });
    }
}

#endif // ERF_POLYGON_IGNITION_H
