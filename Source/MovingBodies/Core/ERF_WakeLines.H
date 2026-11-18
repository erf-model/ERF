#ifndef ERF_WAKE_LINES_H
#define ERF_WAKE_LINES_H

/** Wake sampling lines behind one rotor and their running time average. */

#include <array>
#include <string>
#include <vector>

#include <AMReX_REAL.H>

namespace erf_actuator {

/**
 * Wake sampling lines behind one rotor: at each downstream distance x/D along the horizontal
 * projection of the rotor axis through the hub (a wake follows the wind at hub height, not a
 * tilted shaft), one lateral line (horizontal, normal to that direction) and one vertical
 * line, each of num_points points spanning +-half_width diameters about the centre. The
 * sampled velocities are written as they are and accumulated into a running time average.
 *
 * Point order: for each x/D, the lateral line then the vertical line, each from -half_width
 * to +half_width; the vertical line is clipped at z_min (the ground), so its lowest point may
 * sit above -half_width. Every point's offset s (in diameters) is recorded in the files. The
 * lateral line runs along (-axis_y, axis_x, 0), from the right of the downstream direction to
 * its left; the vertical line runs up. Positions in m, ERF domain coordinates; velocities in m/s.
 */
class WakeLines
{
public:
    /**
     * Build the sampling points.
     *
     * @param hub         the rotor's hub (m)
     * @param axis        the rotor axis (any non-vertical vector; its horizontal projection is used)
     * @param diameter    the rotor diameter (m, > 0)
     * @param lines_xD    downstream distances of the lines, in diameters (at least one)
     * @param half_width  half the length of each line, in diameters (> 0)
     * @param num_points  points per line (>= 2)
     * @param z_min       the ground under the hub (m, finite): the vertical line starts no lower
     */
    WakeLines (std::string name, std::string output_root,
               const std::array<amrex::Real,3>& hub, const std::array<amrex::Real,3>& axis,
               amrex::Real diameter, const std::vector<amrex::Real>& lines_xD,
               amrex::Real half_width, int num_points, amrex::Real z_min);

    /**
     * The same with one ground height per downstream distance (z_min[n] under the centre of
     * the lines at lines_xD[n], from terrain_heights), so that on terrain rising or falling
     * downstream each vertical line stops at its own ground. z_min must have lines_xD's size.
     */
    WakeLines (std::string name, std::string output_root,
               const std::array<amrex::Real,3>& hub, const std::array<amrex::Real,3>& axis,
               amrex::Real diameter, const std::vector<amrex::Real>& lines_xD,
               amrex::Real half_width, int num_points, const std::vector<amrex::Real>& z_min);

    const std::string& name () const { return m_name; }
    int num_points () const { return static_cast<int>(m_pos.size() / 3); }
    /** 3*n positions of the sampling points (m). */
    const std::vector<amrex::Real>& positions () const { return m_pos; }
    const std::array<amrex::Real,3>& hub () const { return m_hub; }
    /** The unit downstream direction: the rotor axis projected on the horizontal plane. */
    const std::array<amrex::Real,3>& axis () const { return m_axis; }
    amrex::Real diameter () const { return m_diameter; }
    int num_samples () const { return m_count; }

    /** Add one set of sampled velocities (3*n, m/s, finite) to the running average. */
    void accumulate (const std::vector<amrex::Real>& vel);
    /** The running mean velocities (3*n); zeros before the first sample. */
    std::vector<amrex::Real> average () const;

    /**
     * Write the sampled velocities to <output_root>_wake.csv, one row per point: truncate starts
     * the file afresh, otherwise the rows are appended; the header is written when open_log()
     * asks for one. Units: time s, x y z m, u v w m/s.
     */
    void write_instantaneous (double time, const std::vector<amrex::Real>& vel, bool truncate) const;
    /** Rewrite <output_root>_wake_avg.csv with the running average (the file always holds the current mean). */
    void write_average () const;

    /**
     * The lines' state for a checkpoint, <dir>/<name>_wake_avg.dat: the sampling geometry
     * (hub, axis, diameter, the points) and the running sums. The geometry is stored because a
     * body model may rebuild the hub and diameter from lower-precision data after a restart and
     * land a few ulps away; the stored points keep the average on the same positions.
     */
    void write_state (const std::string& dir) const;
    /**
     * Restore the state written by write_state(), geometry included; false when there is none.
     * Aborts, changing nothing, when the file is malformed: a negative count, a diameter that is
     * not positive, an axis that is not a horizontal unit vector, or a non-finite value.
     */
    bool read_state (const std::string& dir);

private:
    void write_rows (std::ofstream& out, const std::vector<amrex::Real>& vel, const std::string& prefix) const;

    std::string m_name, m_output_root;
    std::array<amrex::Real,3> m_hub, m_axis, m_lateral;
    amrex::Real m_diameter{0.0};
    std::vector<amrex::Real> m_xD;
    amrex::Real m_half_width{1.5};
    int m_npts{0};
    std::vector<amrex::Real> m_s;     //!< n offsets along the lines, in diameters
    std::vector<amrex::Real> m_pos;   //!< 3*n
    std::vector<double> m_sum;        //!< 3*n running sums, in double whatever Real is
    int m_count{0};
};

} // namespace erf_actuator

#endif
