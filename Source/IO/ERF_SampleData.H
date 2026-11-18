/**
 * \file ERF_SampleData.H
 */
#ifndef ERF_SAMPLEDATA_H
#define ERF_SAMPLEDATA_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_RealVect.H>

#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include <ERF_Container.H>
#include <ERF_DataStruct.H>
#include <ERF_LatLonMap.H>

#include <iostream>

/**
 * The value line and plane sampling write at cells covered by immersed terrain
 * (erf.sampling_covered_value)
 */
inline amrex::Real
sampling_covered_value ()
{
    amrex::Real value = amrex::Real(-9999.0);
    amrex::ParmParse pp("erf");
    pp.queryAdd("sampling_covered_value", value);
    return value;
}

/**
 * The number of levels that exist, counted from level 0: vars_new holds a slot for every
 * level up to amr.max_level, but the state of a level that has not been made is undefined
 */
inline int
sampling_num_levels (const amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
{
    int nlev = 0;
    while (nlev < static_cast<int>(vars_new.size()) &&
           static_cast<int>(vars_new[nlev].size()) > Vars::cons &&
           vars_new[nlev][Vars::cons].ok()) {
        ++nlev;
    }
    return nlev;
}

struct LineSampler
{
    LineSampler ()
    {
        amrex::ParmParse pp("erf");

        bool has_line_idx    = pp.contains("sample_line_lo") || pp.contains("sample_line_hi");
        bool has_line_real   = pp.contains("sample_line_lo_real") || pp.contains("sample_line_hi_real");
        bool has_line_latlon = pp.contains("sample_line_lo_latlon") || pp.contains("sample_line_hi_latlon");
        if (int(has_line_idx) + int(has_line_real) + int(has_line_latlon) > 1) {
            amrex::Abort("Specify only one of erf.sample_line_lo/hi, erf.sample_line_lo_real/hi_real"
                         " or erf.sample_line_lo_latlon/hi_latlon");
        }

        // Count number of lo and hi points define the line
        int n_line_lo  = pp.countval("sample_line_lo") / AMREX_SPACEDIM;
        int n_line_hi  = pp.countval("sample_line_hi") / AMREX_SPACEDIM;
        int n_line_lo_real = pp.countval("sample_line_lo_real") / AMREX_SPACEDIM;
        int n_line_hi_real = pp.countval("sample_line_hi_real") / AMREX_SPACEDIM;
        int n_line_lo_latlon = pp.countval("sample_line_lo_latlon") / 3;
        int n_line_hi_latlon = pp.countval("sample_line_hi_latlon") / 3;
        int n_line_dir = pp.countval("sample_line_dir");
        if (has_line_idx) {
            AMREX_ALWAYS_ASSERT( (n_line_lo==n_line_hi ) &&
                                 (n_line_lo==n_line_dir) );
        } else if (has_line_real) {
            AMREX_ALWAYS_ASSERT( (n_line_lo_real==n_line_hi_real ) &&
                                 (n_line_lo_real==n_line_dir) );
        } else if (has_line_latlon) {
            if (pp.countval("sample_line_lo_latlon") != 3*n_line_lo_latlon ||
                pp.countval("sample_line_hi_latlon") != 3*n_line_hi_latlon ||
                n_line_lo_latlon != n_line_hi_latlon || n_line_lo_latlon != n_line_dir) {
                amrex::Abort("erf.sample_line_lo_latlon and erf.sample_line_hi_latlon must each give"
                             " lat, long and z for every line of erf.sample_line_dir");
            }
        }

        // Parse the data
        int nline = 0;
        if (has_line_idx) {
            nline = n_line_lo;
            if (nline > 0) {
                // Parse lo
                amrex::Vector<int> idx_lo; idx_lo.resize(nline*AMREX_SPACEDIM);
                amrex::Vector<amrex::IntVect> iv_lo; iv_lo.resize(nline);
                pp.queryarr("sample_line_lo",idx_lo,0,nline*AMREX_SPACEDIM);
                for (int i(0); i < nline; i++) {
                    amrex::IntVect iv(idx_lo[AMREX_SPACEDIM*i+0],
                                      idx_lo[AMREX_SPACEDIM*i+1],
                                      idx_lo[AMREX_SPACEDIM*i+2]);
                    iv_lo[i] = iv;
                }

                // Parse hi
                amrex::Vector<int> idx_hi; idx_hi.resize(nline*AMREX_SPACEDIM);
                amrex::Vector<amrex::IntVect> iv_hi; iv_hi.resize(nline);
                pp.queryarr("sample_line_hi",idx_hi,0,nline*AMREX_SPACEDIM);
                for (int i(0); i < nline; i++) {
                    amrex::IntVect iv(idx_hi[AMREX_SPACEDIM*i+0],
                                      idx_hi[AMREX_SPACEDIM*i+1],
                                      idx_hi[AMREX_SPACEDIM*i+2]);
                    iv_hi[i] = iv;
                }

                // Construct vector of bounding boxes
                m_bnd_bx.resize(nline);
                for (int i = 0; i < nline; i++){
                    amrex::Box lbx(iv_lo[i],iv_hi[i]);
                    m_bnd_bx[i] = lbx;
                }

                // The indices are those of level 0; on a finer level the same indices are
                // somewhere else, so these lines are sampled on level 0 only
                m_max_level = 0;
                int amr_max_level = 0;
                amrex::ParmParse("amr").query("max_level", amr_max_level);
                if (amr_max_level > 0) {
                    amrex::Print() << "Line sampling: erf.sample_line_lo/hi are level-0 indices, so "
                                   << "these lines are sampled on level 0 only; give "
                                   << "erf.sample_line_lo_real/hi_real to sample the finest level "
                                   << "that holds a line" << std::endl;
                }
            }
        } else if (has_line_real) {
            m_use_real_bx = true;
            nline = n_line_lo_real;
            if (nline > 0) {
                // Parse lo
                amrex::Vector<amrex::Real> real_lo; real_lo.resize(nline*AMREX_SPACEDIM);
                amrex::Vector<amrex::Vector<amrex::Real>> rv_lo;
                pp.queryarr("sample_line_lo_real",real_lo,0,nline*AMREX_SPACEDIM);
                for (int i(0); i < nline; i++) {
                    amrex::Vector<amrex::Real> rv = {real_lo[AMREX_SPACEDIM*i+0],
                                                     real_lo[AMREX_SPACEDIM*i+1],
                                                     real_lo[AMREX_SPACEDIM*i+2]};
                    rv_lo.push_back(rv);
                }

                // Parse hi
                amrex::Vector<amrex::Real> real_hi; real_hi.resize(nline*AMREX_SPACEDIM);
                amrex::Vector<amrex::Vector<amrex::Real>> rv_hi;
                pp.queryarr("sample_line_hi_real",real_hi,0,nline*AMREX_SPACEDIM);
                for (int i(0); i < nline; i++) {
                    amrex::Vector<amrex::Real> rv = {real_hi[AMREX_SPACEDIM*i+0],
                                                     real_hi[AMREX_SPACEDIM*i+1],
                                                     real_hi[AMREX_SPACEDIM*i+2]};
                    rv_hi.push_back(rv);
                }

                // Construct vector of bounding real boxes
                m_bnd_rbx.resize(nline);
                for (int i = 0; i < nline; i++){
                    amrex::RealBox rbx(rv_lo[i].data(),rv_hi[i].data());
                    m_bnd_rbx[i] = rbx;
                }

                // Construct vector of bounding boxes
                m_bnd_bx.resize(nline);
            }
        } else if (has_line_latlon) {
            // The ends are turned into domain coordinates by resolve_latlon, once the run
            //    has its latitude/longitude arrays; until then the line is not along any
            //    direction, so its ends are checked there
            m_use_real_bx = true;
            m_use_latlon  = true;
            nline = n_line_lo_latlon;
            if (nline > 0) {
                pp.getarr("sample_line_lo_latlon", m_latlon_lo, 0, 3*nline);
                pp.getarr("sample_line_hi_latlon", m_latlon_hi, 0, 3*nline);
                m_bnd_rbx.resize(nline);
                m_bnd_bx.resize(nline);
            }
        }

        if (nline > 0) {
            // Parse directionality
            m_dir.resize(n_line_dir);
            pp.queryarr("sample_line_dir",m_dir,0,n_line_dir);

            // A line runs along one direction: its two ends must agree in the other two, or
            // the sampler would take a slab of cells and add them up at each point
            for (int iline(0); iline<nline; ++iline) {
                const int dir = m_dir[iline];
                if (dir < 0 || dir >= AMREX_SPACEDIM) {
                    amrex::Abort("erf.sample_line_dir must be 0, 1 or 2");
                }
                if (m_use_latlon) { continue; }
                for (int d(0); d<AMREX_SPACEDIM; ++d) {
                    const bool differ = has_line_idx
                        ? (m_bnd_bx[iline].smallEnd(d) != m_bnd_bx[iline].bigEnd(d))
                        : (m_bnd_rbx[iline].lo(d) != m_bnd_rbx[iline].hi(d));
                    const bool reversed = has_line_idx
                        ? (m_bnd_bx[iline].smallEnd(d) > m_bnd_bx[iline].bigEnd(d))
                        : (m_bnd_rbx[iline].lo(d) > m_bnd_rbx[iline].hi(d));
                    if (d != dir && differ) {
                        amrex::Abort("Line " + std::to_string(iline) + " of erf.sample_line_lo/hi"
                                     + std::string(has_line_idx ? "" : "_real")
                                     + " is not along erf.sample_line_dir: its ends differ in direction "
                                     + std::to_string(d));
                    }
                    if (d == dir && reversed) {
                        amrex::Abort("Line " + std::to_string(iline) + " of erf.sample_line_lo/hi"
                                     + std::string(has_line_idx ? "" : "_real")
                                     + " has its lo end above its hi end");
                    }
                }
            }

            // Parse names
            std::string name_base = "plt_line_";
            m_name.resize(nline);
            int n_names = pp.countval("sample_line_name");
            if (n_names > 0) {
                AMREX_ALWAYS_ASSERT( n_names==nline );
                pp.queryarr("sample_line_name",m_name,0,n_names);
            } else {
                for (int iline(0); iline<nline; ++iline) {
                    m_name[iline] = amrex::Concatenate(name_base, iline , 5);
                }
            }

            // Allocate space for level indicator
            m_lev.resize(n_line_dir,0);

            // Allocate space for MF pointers
            m_ls_mf.resize(nline);

            // Get requested vars.  They are stored, sampled and written in the order of
            // line_sampling_var_order(), whatever order they are requested in
            if (pp.countval("line_sampling_vars") > 0) {
                m_varnames.clear();
                amrex::Vector<std::string> requested_vars;
                pp.queryarr("line_sampling_vars",requested_vars);
                amrex::Print() << "Selected line sampling vars :";
                for (const auto& name : line_sampling_var_order()) {
                    if (containerHasElement(requested_vars, name)) {
                        m_varnames.push_back(name);
                        amrex::Print() << " " << name;
                    }
                }
                amrex::Print() << std::endl;
            }

            // Write outputs to text files, one file per variable, with all
            // times appended to the same file.  The files are opened at the first
            // write, once "conserved" has been expanded into its components
            pp.queryAdd("line_sampling_text_output",m_write_ascii);

            // The z of erf.sample_line_lo_real/hi_real of a vertical line is the height above
            // the ground under it instead of the height above sea level
            pp.queryAdd("sample_line_z_agl",m_z_agl);

            // Interpolate to the requested heights rather than take the cells that hold them:
            //    a vertical line at erf.sample_line_npts heights from its lo z to its hi z, a
            //    horizontal line at its z in every column it crosses
            pp.queryAdd("sample_line_interp",m_interp);
            if (m_interp) {
                if (!m_use_real_bx) {
                    amrex::Abort("erf.sample_line_interp needs lines given by heights, with "
                                 "erf.sample_line_lo_real/hi_real or erf.sample_line_lo_latlon/hi_latlon; "
                                 "a line given by erf.sample_line_lo/hi names cells, not heights");
                }
                m_npts.assign(nline, 0);
                const int n_npts = pp.countval("sample_line_npts");
                if (n_npts > 0) {
                    if (n_npts != nline) {
                        amrex::Abort("erf.sample_line_npts must give one value per line (0 for a horizontal line)");
                    }
                    pp.getarr("sample_line_npts", m_npts, 0, nline);
                } else {
                    // A vertical line whose two ends are at one height is sampled at that height
                    for (int iline(0); iline<nline; ++iline) {
                        if (m_dir[iline] == 2) { m_npts[iline] = 1; }
                    }
                }
                // The heights of lines given by latitude and longitude are checked once they are placed
                if (!m_use_latlon) {
                    for (int iline(0); iline<nline; ++iline) { check_interp_heights(iline); }
                }
            } else if (pp.contains("sample_line_npts")) {
                amrex::Abort("erf.sample_line_npts is used only with erf.sample_line_interp = true");
            }

            // The cells of a line placed by height are found once, so they do not follow a
            // mesh that moves; lines given by indices do not use heights
            if (m_use_real_bx && use_column_heights() &&
                SolverChoice::terrain_type == TerrainType::MovingFittedMesh) {
                amrex::Abort("Line sampling with erf.sample_line_lo_real/hi_real is not supported "
                             "with erf.terrain_type = MovingFittedMesh");
            }

            // Where the levels are not equally spaced, the z of the plotfile of a line is not
            // the height of its samples, so the height is written with them; an interpolated
            // line writes the height of each sample above sea level on any mesh
            if ((SolverChoice::mesh_type != MeshType::ConstantDz || m_interp) &&
                !containerHasElement(m_varnames, std::string("z"))) {
                const bool cons_first = !m_varnames.empty() && m_varnames[0] == "conserved";
                m_varnames.insert(m_varnames.begin() + (cons_first ? 1 : 0), "z");
            }
        }
    }

    /**
     * Sample the full-depth vertical column through each cell on level 0, with the
     * variables erf.sample_line writes: every conserved component, the cell-centred
     * velocities and the six stresses at cell centres.  Each column is written to its own
     * log, one row per sample: the time, then each variable over the column, bottom to top.
     *
     * @param[in] cells     a cell in each column; its k is ignored
     * @param[in] log_names the log file of each column
     * @param[in] domain    the level-0 domain
     */
    LineSampler (const amrex::Vector<amrex::IntVect>& cells,
                 const amrex::Vector<std::string>& log_names,
                 const amrex::Box& domain)
    {
        AMREX_ALWAYS_ASSERT(cells.size() == log_names.size());
        const int nline = static_cast<int>(cells.size());
        constexpr int dir = 2;

        m_bnd_bx.resize(nline);
        for (int iline(0); iline<nline; ++iline) {
            amrex::IntVect lo(cells[iline]); lo[dir] = domain.smallEnd(dir);
            amrex::IntVect hi(cells[iline]); hi[dir] = domain.bigEnd(dir);
            m_bnd_bx[iline] = amrex::Box(lo, hi);
        }
        m_dir.assign(nline, dir);
        m_lev.assign(nline, 0);
        m_ls_mf.resize(nline);
        m_name = log_names;
        m_max_level = 0;
        m_varnames = {"conserved", "x_velocity", "y_velocity", "z_velocity",
                      "sgs_tau11", "sgs_tau12", "sgs_tau13", "sgs_tau22", "sgs_tau23", "sgs_tau33"};

        open_logs(log_names);
    }

    /**
     * Sample single cells with every conserved component, as erf.sample_point does: each
     * point is written to its own log, one row per sample, the time and then each component.
     * A point given by a cell (erf.sample_point) is that cell of level 0.  A point given by
     * real coordinates (erf.sample_point_real, or erf.sample_point_latlon once placed) is a
     * line of one cell, sampled on the finest level that holds it; on a stretched or
     * terrain-fitted mesh its cell is found from the heights of its column, and with z_agl
     * its z is the height above the ground under it.
     *
     * @param[in] cells     the cell of each point given by indices, or empty
     * @param[in] points    the (x,y,z) of each point given by real coordinates, or empty
     * @param[in] log_names the log file of each point
     * @param[in] z_agl     whether the z of points are heights above the ground
     */
    LineSampler (const amrex::Vector<amrex::IntVect>& cells,
                 const amrex::Vector<amrex::RealVect>& points,
                 const amrex::Vector<std::string>& log_names,
                 bool z_agl)
    {
        AMREX_ALWAYS_ASSERT(cells.empty() || points.empty());
        const int nline = static_cast<int>(log_names.size());
        AMREX_ALWAYS_ASSERT(static_cast<int>(cells.size() + points.size()) == nline);

        m_bnd_bx.resize(nline);
        if (!cells.empty()) {
            for (int iline(0); iline<nline; ++iline) {
                m_bnd_bx[iline] = amrex::Box(cells[iline], cells[iline]);
            }
            m_max_level = 0;
        } else {
            m_use_real_bx = true;
            m_z_agl = z_agl;
            m_bnd_rbx.resize(nline);
            for (int iline(0); iline<nline; ++iline) {
                m_bnd_rbx[iline] = amrex::RealBox(points[iline].dataPtr(), points[iline].dataPtr());
            }
            // As for lines placed by height, the cells are found once
            if (use_column_heights() && SolverChoice::terrain_type == TerrainType::MovingFittedMesh) {
                amrex::Abort("erf.sample_point_real and erf.sample_point_latlon are not supported "
                             "with erf.terrain_type = MovingFittedMesh");
            }
        }
        // A point is a vertical line of one cell
        m_dir.assign(nline, 2);
        m_lev.assign(nline, 0);
        m_ls_mf.resize(nline);
        m_name = log_names;
        m_varnames = {"conserved"};

        open_logs(log_names);
    }

    //! Open one log per line, appended to, on the I/O rank
    void
    open_logs (const amrex::Vector<std::string>& log_names)
    {
        m_write_log = true;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            const int nline = static_cast<int>(log_names.size());
            m_logstream.resize(nline);
            for (int iline(0); iline<nline; ++iline) {
                m_logstream[iline] = std::make_unique<std::fstream>();
                m_logstream[iline]->open(log_names[iline].c_str(),std::ios::out|std::ios::app);
                if (!m_logstream[iline]->good()) {
                    amrex::FileOpenFailed(log_names[iline]);
                }
            }
        }
        amrex::ParallelDescriptor::Barrier("LineSampler::open_logs");
    }

    //! The variables that may be sampled, in the order they are stored and written
    static const amrex::Vector<std::string>&
    line_sampling_var_order ()
    {
        static const amrex::Vector<std::string> order = {
            "conserved", "z", "volfrac", "density", "x_velocity", "y_velocity", "z_velocity", "magvel", "theta",
            "qv", "qc", "pressure", "sgs_tke",
            "sgs_tau11", "sgs_tau12", "sgs_tau13", "sgs_tau22", "sgs_tau23", "sgs_tau33", "sgs_hfx3"};
        return order;
    }

    //! Replace "conserved" by one name per component of the conserved state, cons_0, cons_1, ...
    void
    expand_conserved (int ncons)
    {
        auto it = std::find(m_varnames.begin(), m_varnames.end(), std::string("conserved"));
        if (it == m_varnames.end()) { m_ncons = 0; return; }
        amrex::Vector<std::string> names;
        for (int n(0); n<ncons; ++n) { names.push_back(amrex::Concatenate("cons_", n, 1)); }
        it = m_varnames.erase(it);
        m_varnames.insert(it, names.begin(), names.end());
        m_ncons = ncons;
    }

    void
    open_ascii_streams ()
    {
        if (!m_datastream.empty() || !amrex::ParallelDescriptor::IOProcessor()) { return; }
        int nline = static_cast<int>(m_name.size());
        int nvar = static_cast<int>(m_varnames.size());
        m_datastream.resize(nline * nvar);
        int i = 0;
        for (int iline(0); iline<nline; ++iline) {
            for (int ivar(0); ivar<nvar; ++ivar) {
                std::string filename = m_name[iline] + "." + m_varnames[ivar];
                m_datastream[i] = std::make_unique<std::fstream>();
                m_datastream[i]->open(filename.c_str(),std::ios::out|std::ios::app);
                if (!m_datastream[i]->good()) {
                    amrex::FileOpenFailed(filename);
                }
                i++;
            }
        }
    }

    //! Whether the lines are given by latitude and longitude (erf.sample_line_lo/hi_latlon)
    [[nodiscard]] bool
    uses_latlon () const
    {
        return m_use_latlon;
    }

    /**
     * Turn the ends of the lines given by latitude and longitude into domain coordinates,
     * with the level-0 latitude/longitude map of the run, and check them.  A vertical line
     * gives one lat/lon at both ends.  The two ends of a horizontal line are in general not
     * on one row of cells: they must lie within one cell of each other across the line,
     * which is then sampled along the line of cells through its lo end.  Collective.
     */
    void
    resolve_latlon (const LatLonMap& latlon_map, const amrex::Geometry& geom0)
    {
        const int nline = static_cast<int>(m_bnd_rbx.size());
        for (int iline(0); iline<nline; ++iline) {
            const int dir = m_dir[iline];
            const std::string line = "Line " + std::to_string(iline) + " of erf.sample_line_lo/hi_latlon";
            const amrex::Real* req[2] = {&m_latlon_lo[3*iline], &m_latlon_hi[3*iline]};

            if (dir == 2) {
                if (req[0][0] != req[1][0] || req[0][1] != req[1][1]) {
                    amrex::Abort(line + " is vertical (erf.sample_line_dir = 2), so its two ends"
                                 " must give the same latitude and longitude");
                }
                if (req[0][2] > req[1][2]) {
                    amrex::Abort(line + " has its lo end above its hi end");
                }
            } else if (req[0][2] != req[1][2]) {
                amrex::Abort(line + " is horizontal, so its two ends must give the same z");
            }

            // A vertical line has one lat/lon
            const int nend = (dir == 2) ? 1 : 2;
            LatLonLocation where[2];
            for (int e(0); e<nend; ++e) {
                where[e] = latlon_map.locate_or_abort(req[e][0], req[e][1],
                                                      line + (e == 0 ? ", lo end," : ", hi end,"));
            }
            if (dir == 2) { where[1] = where[0]; }

            amrex::Real lo[AMREX_SPACEDIM] = {where[0].x, where[0].y, req[0][2]};
            amrex::Real hi[AMREX_SPACEDIM] = {where[1].x, where[1].y, req[1][2]};
            amrex::Real snap = amrex::Real(0.0);
            if (dir != 2) {
                const int across = 1 - dir;
                snap = hi[across] - lo[across];
                const amrex::Real gap = std::abs(snap) / geom0.CellSize(across);
                if (gap > amrex::Real(1.0)) {
                    amrex::Abort(line + " is not along erf.sample_line_dir: its ends are " +
                                 std::to_string(gap) + " cells apart in direction " +
                                 std::to_string(across) + ", and a line between two lat/lon"
                                 " points is sampled along one grid direction, so its ends may"
                                 " be at most one cell apart across it");
                }
                if (hi[dir] < lo[dir]) {
                    amrex::Abort(line + " has its lo end beyond its hi end in direction " +
                                 std::to_string(dir) + "; swap the two ends");
                }
                hi[across] = lo[across];
            }
            m_bnd_rbx[iline] = amrex::RealBox(lo, hi);
            if (m_interp) { check_interp_heights(iline); }

            amrex::Print() << line << ": lo end at lat=" << where[0].lat << " long=" << where[0].lon
                           << " (x=" << lo[0] << " y=" << lo[1] << ")";
            if (dir != 2) {
                amrex::Print() << ", hi end at lat=" << where[1].lat << " long=" << where[1].lon
                               << " (x=" << hi[0] << " y=" << hi[1] << "), sampled along the line of cells"
                               << " through the lo end, " << std::abs(snap) << " m from the hi end";
            }
            amrex::Print() << std::endl;
        }
    }

    //! Whether the sampler needs the ground under its lines (erf.sample_line_z_agl, or
    //! erf.sample_line_interp, which writes the samples below the ground as missing)
    [[nodiscard]] bool
    needs_ground () const
    {
        return m_use_real_bx && (m_z_agl || m_interp);
    }

    //! Abort unless the number of heights of line iline suits its ends (erf.sample_line_interp)
    void
    check_interp_heights (int iline) const
    {
        const std::string line = "Line " + std::to_string(iline) + " of erf.sample_line_interp";
        const int n = m_npts[iline];
        const amrex::RealBox& rb = m_bnd_rbx[iline];
        if (m_dir[iline] != 2) {
            if (n != 0) {
                amrex::Abort(line + " is horizontal, so it is sampled once in every column it crosses: "
                             "its erf.sample_line_npts must be 0");
            }
        } else if (n < 1) {
            amrex::Abort(line + " is vertical: give erf.sample_line_npts, the number of heights (at least 1)");
        } else if (n == 1 && rb.lo(2) != rb.hi(2)) {
            amrex::Abort(line + " has one height, so its lo z and hi z must be equal; give erf.sample_line_npts "
                         "of at least 2 for heights from its lo z to its hi z");
        } else if (n > 1 && !(rb.lo(2) < rb.hi(2))) {
            amrex::Abort(line + " has " + std::to_string(n) + " heights, so its lo z must be below its hi z");
        }
    }

    //! The requested heights of line iline (erf.sample_line_interp): erf.sample_line_npts evenly
    //! spaced from lo z to hi z for a vertical line, its one z for a horizontal line
    [[nodiscard]] amrex::Vector<amrex::Real>
    interp_heights (int iline) const
    {
        const amrex::RealBox& rb = m_bnd_rbx[iline];
        if (m_dir[iline] != 2) { return {rb.lo(2)}; }
        const int n = m_npts[iline];
        amrex::Vector<amrex::Real> h(n, rb.lo(2));
        for (int m(1); m<n; ++m) {
            h[m] = rb.lo(2) + (rb.hi(2) - rb.lo(2)) * amrex::Real(m) / amrex::Real(n - 1);
        }
        return h;
    }

    //! The samples of line iline, along its direction: its heights (vertical) or its columns (horizontal)
    [[nodiscard]] amrex::Box
    interp_sample_box (int iline, const amrex::Box& cols) const
    {
        if (m_dir[iline] != 2) { return cols; }
        amrex::IntVect lo(cols.smallEnd()); lo[2] = 0;
        amrex::IntVect hi(cols.smallEnd()); hi[2] = m_npts[iline] - 1;
        return amrex::Box(lo, hi);
    }

    /**
     * Components scomp .. scomp+ncomp-1 of src over the cells (or nodes) of bx, in the order of
     * a FArrayBox on bx, on the I/O rank, and on every rank with to_all.  Collective.
     */
    static amrex::Vector<amrex::Real>
    gather_box (const amrex::MultiFab& src, int scomp, int ncomp, const amrex::Box& bx, bool to_all)
    {
        const int io = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::MultiFab dst(amrex::BoxArray(bx), amrex::DistributionMapping(amrex::Vector<int>{io}),
                            ncomp, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
        dst.setVal(amrex::Real(0.0));
        dst.ParallelCopy(src, scomp, 0, ncomp);
        amrex::Gpu::streamSynchronize();
        amrex::Vector<amrex::Real> v(static_cast<std::size_t>(bx.numPts()) * ncomp, amrex::Real(0.0));
        if (amrex::ParallelDescriptor::IOProcessor()) {
            const amrex::Real* p = dst[0].dataPtr();
            std::copy(p, p + v.size(), v.begin());
        }
        if (to_all) { amrex::ParallelDescriptor::Bcast(v.data(), v.size(), io); }
        return v;
    }

    /**
     * Where the samples of a line with erf.sample_line_interp come from on one level, with the
     * grids that was found on
     */
    struct InterpPlan
    {
        amrex::BoxArray ba;                 //!< the grids of the level
        bool holds{false};                  //!< the level holds every cell the samples read
        amrex::Box cols;                    //!< the columns of the line on the level (k = the domain's bottom)
        amrex::Box cells;                   //!< the cells the samples read (empty if none)
        amrex::Vector<amrex::Real> zt;      //!< the height of each sample above sea level
        amrex::Vector<int> col, k0, k1;     //!< per sample: its column along the line, the cells it reads
        amrex::Vector<amrex::Real> w;       //!< the weight of k1
        amrex::Vector<char> missing;        //!< the sample is below the ground or above the top, or no fluid brackets it
    };

    /**
     * The samples of line iline on a level (erf.sample_line_interp).  Each sample is in a column of
     * cells of the level: the one that holds the line (vertical), or each one it crosses
     * (horizontal).  Its height is the requested one, above sea level or, with
     * erf.sample_line_z_agl, above the ground under the column.  The ground and the top of a
     * column are those of the level, the means of the four nodes under and over it, where the
     * level holds the column's bottom or top cell; elsewhere the ground is the level-0 ground,
     * bilinear at the column's centre (a level that stops above the ground).  Its value is
     * linear in height between the two cell centres of the column that bracket it; between the
     * ground and the lowest cell centre with fluid it is that cell's, and between the top cell
     * centre and the top the top cell's.  A sample below the ground or above the top is
     * missing.  The level holds the line if it holds every cell the samples read; level 0
     * always does.  Collective: every rank gets the same plan.
     *
     * @param[in] z_nd     the nodal heights of the level, or nullptr (the top is then prob_hi)
     * @param[in] volfrac  the volume fraction of the level with immersed terrain, or nullptr
     * @param[in] ground   the ground of the level (z_nd, or the immersed terrain surface), or nullptr
     * @param[in] ground0  the ground of level 0, or nullptr (the ground is then prob_lo)
     */
    InterpPlan
    make_interp_plan (int iline, int lev, const amrex::Geometry& geom, const amrex::Geometry& geom0,
                      const amrex::MultiFab& z_cc, const amrex::MultiFab* z_nd,
                      const amrex::MultiFab* volfrac, const amrex::MultiFab* ground,
                      const amrex::MultiFab* ground0) const
    {
        InterpPlan plan;
        plan.ba = z_cc.boxArray();
        plan.cols = lineIndexBoxXY(m_bnd_rbx[iline], geom);
        const amrex::Box& dom = geom.Domain();
        const int kdlo = dom.smallEnd(2);
        const int nz = dom.length(2);
        const int ncol = static_cast<int>(plan.cols.numPts());
        const int dir = m_dir[iline];
        const int along = (dir == 2) ? 0 : dir;
        const amrex::IntVect clo = plan.cols.smallEnd();

        // The heights of the cell centres of the columns, which cells the level holds, and their
        // volume fractions, from the bottom of the domain to its top
        amrex::Box colbox(plan.cols);
        colbox.setSmall(2, kdlo);
        colbox.setBig(2, dom.bigEnd(2));
        amrex::MultiFab t(z_cc.boxArray(), z_cc.DistributionMap(), 3, 0);
        amrex::MultiFab::Copy(t, z_cc, 0, 0, 1, 0);
        t.setVal(amrex::Real(1.0), 1, 1, 0);
        if (volfrac != nullptr) {
            amrex::MultiFab::Copy(t, *volfrac, 0, 2, 1, 0);
        } else {
            t.setVal(amrex::Real(1.0), 2, 1, 0);
        }
        const amrex::Vector<amrex::Real> c = gather_box(t, 0, 3, colbox, true);
        const auto ccol = static_cast<std::size_t>(colbox.numPts());
        auto at = [&] (int icol, int k, int n) {
            return c[static_cast<std::size_t>(icol) + static_cast<std::size_t>(ncol) * k + ccol * n];
        };
        auto held_at = [&] (int icol, int k) { return at(icol, k, 1) > amrex::Real(0.5); };

        // The mean of the four nodes at node level k of mf under or over each column of the level
        auto column_node_mean = [&] (const amrex::MultiFab& mf, int k) {
            amrex::IntVect nlo(clo); nlo[2] = k;
            amrex::IntVect nhi(plan.cols.bigEnd() + amrex::IntVect(1,1,0)); nhi[2] = k;
            const amrex::Box nbx(nlo, nhi, amrex::IntVect::TheNodeVector());
            const amrex::Vector<amrex::Real> g = gather_box(mf, 0, 1, nbx, true);
            const int nnx = nbx.length(0);
            amrex::Vector<amrex::Real> mean(ncol);
            for (int icol(0); icol<ncol; ++icol) {
                amrex::IntVect cell(clo);
                cell[along] += icol;
                const int i = cell[0] - nlo[0], j = cell[1] - nlo[1];
                mean[icol] = amrex::Real(0.25) * (g[i + nnx*j] + g[i+1 + nnx*j] + g[i + nnx*(j+1)] + g[i+1 + nnx*(j+1)]);
            }
            return mean;
        };

        // The ground under the centre of each column, bilinear between the level-0 nodes around
        // it: for the columns whose bottom the level does not hold
        amrex::Vector<amrex::Real> ground_col(ncol, geom0.ProbLo(2));
        if (ground0 != nullptr) {
            const amrex::Box& dom0 = geom0.Domain();
            amrex::Vector<int> i0(ncol), j0(ncol);
            amrex::Vector<amrex::Real> fx(ncol), fy(ncol);
            for (int icol(0); icol<ncol; ++icol) {
                amrex::IntVect cell(clo);
                cell[along] += icol;
                int* idx[2] = {&i0[icol], &j0[icol]};
                amrex::Real* f[2] = {&fx[icol], &fy[icol]};
                for (int d(0); d<2; ++d) {
                    const amrex::Real xc = geom.ProbLo(d) + (cell[d] + amrex::Real(0.5)) * geom.CellSize(d);
                    const amrex::Real s = (xc - geom0.ProbLo(d)) / geom0.CellSize(d);
                    *idx[d] = std::min(std::max(static_cast<int>(std::floor(s)), dom0.smallEnd(d)), dom0.bigEnd(d));
                    *f[d] = std::min(std::max(s - *idx[d], amrex::Real(0.0)), amrex::Real(1.0));
                }
            }
            const amrex::IntVect nlo(*std::min_element(i0.begin(), i0.end()), *std::min_element(j0.begin(), j0.end()),
                                     dom0.smallEnd(2));
            const amrex::IntVect nhi(*std::max_element(i0.begin(), i0.end()) + 1,
                                     *std::max_element(j0.begin(), j0.end()) + 1, dom0.smallEnd(2));
            const amrex::Box nbx(nlo, nhi, amrex::IntVect::TheNodeVector());
            const amrex::Vector<amrex::Real> g = gather_box(*ground0, 0, 1, nbx, true);
            const int nnx = nbx.length(0);
            auto node = [&] (int i, int j) { return g[(i - nlo[0]) + nnx * (j - nlo[1])]; };
            for (int icol(0); icol<ncol; ++icol) {
                const int i = i0[icol], j = j0[icol];
                const amrex::Real wx = fx[icol], wy = fy[icol];
                ground_col[icol] = (amrex::Real(1.0) - wy) * ((amrex::Real(1.0) - wx) * node(i,j)   + wx * node(i+1,j))
                                 +                    wy  * ((amrex::Real(1.0) - wx) * node(i,j+1) + wx * node(i+1,j+1));
            }
        }
        // the level's own ground where it holds the column's bottom cell
        if (ground != nullptr) {
            const amrex::Vector<amrex::Real> own = column_node_mean(*ground, kdlo);
            for (int icol(0); icol<ncol; ++icol) {
                if (held_at(icol, 0)) { ground_col[icol] = own[icol]; }
            }
        }
        // The top of each column, where the level holds its top cell; elsewhere no sample above
        // the level's highest cell is bracketed, so the level does not hold it
        amrex::Vector<amrex::Real> top_col(ncol, geom.ProbHi(2));
        if (z_nd != nullptr) {
            const amrex::Vector<amrex::Real> own = column_node_mean(*z_nd, dom.bigEnd(2) + 1);
            for (int icol(0); icol<ncol; ++icol) {
                if (held_at(icol, nz-1)) { top_col[icol] = own[icol]; }
            }
        }

        // The samples: the heights of a vertical line in its column, or one height in each column
        const amrex::Vector<amrex::Real> h = interp_heights(iline);
        const int nsamp = (dir == 2) ? static_cast<int>(h.size()) : ncol;
        plan.zt.resize(nsamp);
        plan.col.resize(nsamp);
        plan.k0.assign(nsamp, 0);
        plan.k1.assign(nsamp, 0);
        plan.w.assign(nsamp, amrex::Real(0.0));
        plan.missing.assign(nsamp, 0);
        plan.holds = true;
        int kmin = nz, kmax = -1;
        for (int s(0); s<nsamp; ++s) {
            const int icol = (dir == 2) ? 0 : s;
            const amrex::Real hs = (dir == 2) ? h[s] : h[0];
            plan.col[s] = icol;
            plan.zt[s] = hs + (m_z_agl ? ground_col[icol] : amrex::Real(0.0));
            const amrex::Real z = plan.zt[s];
            if (z < ground_col[icol] || z > top_col[icol]) { plan.missing[s] = 1; continue; }

            auto held  = [&] (int k) { return held_at(icol, k); };
            auto fluid = [&] (int k) { return held(k) && at(icol, k, 2) > amrex::Real(0.0); };
            auto zc    = [&] (int k) { return at(icol, k, 0); };
            int k0 = -1, k1 = -1;
            amrex::Real w = amrex::Real(0.0);
            // between two cell centres with fluid
            for (int k(0); k<nz-1 && k0 < 0; ++k) {
                if (fluid(k) && fluid(k+1) && zc(k) <= z && z < zc(k+1)) {
                    k0 = k; k1 = k+1; w = (z - zc(k)) / (zc(k+1) - zc(k));
                }
            }
            // between the top cell centre and the top
            if (k0 < 0 && fluid(nz-1) && z >= zc(nz-1)) { k0 = k1 = nz-1; }
            // between the ground and the lowest cell centre with fluid, if the level holds the
            // column up to it
            if (k0 < 0) {
                int kf = 0;
                while (kf < nz && held(kf) && !fluid(kf)) { ++kf; }
                if (kf < nz && fluid(kf) && z < zc(kf)) { k0 = k1 = kf; }
            }
            if (k0 < 0) {
                // Level 0 holds every cell, so there no cells with fluid bracket the height (a
                // column without fluid, or a height under an overhang of immersed terrain); a
                // finer level does not hold the cells around it
                if (lev == 0) { plan.missing[s] = 1; continue; }
                plan.holds = false;
                return plan;
            }
            plan.k0[s] = k0;
            plan.k1[s] = k1;
            plan.w[s] = w;
            kmin = std::min(kmin, k0);
            kmax = std::max(kmax, k1);
        }
        if (kmax >= kmin) {
            plan.cells = plan.cols;
            plan.cells.setSmall(2, kdlo + kmin);
            plan.cells.setBig(2, kdlo + kmax);
        }
        return plan;
    }

    //! The entry of level lev, or nullptr past the end
    static const amrex::MultiFab*
    level_ptr (const amrex::Vector<const amrex::MultiFab*>& v, int lev)
    {
        return lev < static_cast<int>(v.size()) ? v[lev] : nullptr;
    }

    //! make_interp_plan, kept until the grids of the level change
    const InterpPlan&
    cachedInterpPlan (int iline, int lev, const amrex::Geometry& geom, const amrex::Geometry& geom0,
                      const amrex::MultiFab& z_cc, const amrex::MultiFab* z_nd,
                      const amrex::MultiFab* volfrac, const amrex::MultiFab* ground,
                      const amrex::MultiFab* ground0)
    {
        if (lev >= static_cast<int>(m_interp_plan.size())) { m_interp_plan.resize(lev+1); }
        auto& plans = m_interp_plan[lev];
        const int nline = static_cast<int>(m_bnd_rbx.size());
        if (static_cast<int>(plans.size()) != nline) { plans.assign(nline, InterpPlan{}); }
        if (plans[iline].ba.empty() || !(plans[iline].ba == z_cc.boxArray())) {
            plans[iline] = make_interp_plan(iline, lev, geom, geom0, z_cc, z_nd, volfrac, ground, ground0);
        }
        return plans[iline];
    }

    /**
     * The samples of an interpolated line from the variables data of its level, as a MultiFab of
     * one box (interp_sample_box) on the I/O rank.  z is the height of each sample above sea
     * level; a missing sample has volume fraction 0 and erf.sampling_covered_value elsewhere.
     * Collective.
     */
    amrex::MultiFab
    interp_line_data (int iline, const InterpPlan& plan, const amrex::MultiFab& data, int z_comp, int vf_comp) const
    {
        const int ncomp = data.nComp();
        const int nsamp = static_cast<int>(plan.zt.size());
        amrex::Vector<amrex::Real> d;
        if (plan.cells.ok()) { d = gather_box(data, 0, ncomp, plan.cells, false); }

        const amrex::Box sbx = interp_sample_box(iline, plan.cols);
        const int io = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::MultiFab out(amrex::BoxArray(sbx), amrex::DistributionMapping(amrex::Vector<int>{io}), ncomp, 0);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            const int ncol = static_cast<int>(plan.cols.numPts());
            const auto ncell = static_cast<std::size_t>(plan.cells.ok() ? plan.cells.numPts() : 0);
            const int kc0 = plan.cells.ok() ? plan.cells.smallEnd(2) - plan.cols.smallEnd(2) : 0;
            auto cell = [&] (int icol, int k, int n) {
                return d[static_cast<std::size_t>(icol) + static_cast<std::size_t>(ncol) * (k - kc0) + ncell * n];
            };
            // in the order of the box's FArrayBox, whose samples run along its one long direction
            amrex::Vector<amrex::Real> v(static_cast<std::size_t>(nsamp) * ncomp);
            for (int n(0); n<ncomp; ++n) {
                for (int s(0); s<nsamp; ++s) {
                    amrex::Real val;
                    if (n == z_comp) {
                        val = plan.zt[s];
                    } else if (plan.missing[s] != 0) {
                        val = (n == vf_comp) ? amrex::Real(0.0) : m_covered_value;
                    } else {
                        const amrex::Real w = plan.w[s];
                        val = (amrex::Real(1.0) - w) * cell(plan.col[s], plan.k0[s], n);
                        if (w > amrex::Real(0.0)) { val += w * cell(plan.col[s], plan.k1[s], n); }
                    }
                    v[static_cast<std::size_t>(s) + static_cast<std::size_t>(nsamp) * n] = val;
                }
            }
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, v.begin(), v.end(), out[0].dataPtr());
        }
        return out;
    }

    //! Whether the z of a line given by real coordinates is found from the heights of a
    //! column: on a stretched or terrain-fitted mesh, or for a height above ground
    [[nodiscard]] bool
    use_column_heights () const
    {
        return m_z_agl || SolverChoice::mesh_type != MeshType::ConstantDz;
    }

    /**
     * The heights of the faces of the column through cell (i,j) on a level, from the nodal
     * heights z_nd: the lower face of each cell, its upper face, and 1 where the level holds
     * the cell (0 where it does not).  Each is indexed from the bottom of the domain.
     */
    static amrex::Vector<amrex::Gpu::HostVector<amrex::Real>>
    column_face_heights (const amrex::MultiFab& z_nd, const amrex::Geometry& geom,
                         const amrex::IntVect& cell)
    {
        const amrex::Box& dom = geom.Domain();
        amrex::IntVect lo(cell); lo[2] = dom.smallEnd(2);
        amrex::IntVect hi(cell); hi[2] = dom.bigEnd(2);
        const amrex::Box col(lo, hi);

        amrex::MultiFab f(amrex::convert(z_nd.boxArray(), amrex::IntVect(0)),
                          z_nd.DistributionMap(), 3, 0);
        f.setVal(amrex::Real(0.0));
        for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox() & col;
            if (!b.ok()) { continue; }
            const auto& out = f.array(mfi);
            const auto& z = z_nd.const_array(mfi);
            amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                out(i,j,k,0) = amrex::Real(0.25) * (z(i,j,k  ) + z(i+1,j,k  ) + z(i,j+1,k  ) + z(i+1,j+1,k  ));
                out(i,j,k,1) = amrex::Real(0.25) * (z(i,j,k+1) + z(i+1,j,k+1) + z(i,j+1,k+1) + z(i+1,j+1,k+1));
                out(i,j,k,2) = amrex::Real(1.0);
            });
        }
        return {amrex::sumToLine(f, 0, 1, col, 2), amrex::sumToLine(f, 1, 1, col, 2),
                amrex::sumToLine(f, 2, 1, col, 2)};
    }

    /**
     * The height of the ground under cell (i,j) of level 0: the mean of the four nodes at the
     * bottom of its column in ground_nd, which is z_phys_nd or, with immersed terrain, the
     * terrain surface (ERF::station_ground)
     */
    static amrex::Real
    column_ground (const amrex::MultiFab& ground_nd, const amrex::Geometry& geom,
                   const amrex::IntVect& cell)
    {
        const int i = cell[0], j = cell[1], k = geom.Domain().smallEnd(2);
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
        for (amrex::MFIter mfi(ground_nd); mfi.isValid(); ++mfi) {
            // The nodal box that holds the cell, so that a node shared by two boxes counts once
            const amrex::Box& vb = mfi.validbox();
            if (i < vb.smallEnd(0) || i >= vb.bigEnd(0) ||
                j < vb.smallEnd(1) || j >= vb.bigEnd(1) ||
                k < vb.smallEnd(2) || k >  vb.bigEnd(2)) { continue; }
            const auto& z = ground_nd.const_array(mfi);
            reduce_op.eval(amrex::Box(amrex::IntVect(i,j,k), amrex::IntVect(i,j,k)), reduce_data,
                [=] AMREX_GPU_DEVICE (int ii, int jj, int kk) -> amrex::GpuTuple<amrex::Real> {
                    return { amrex::Real(0.25) * (z(ii,jj,kk) + z(ii+1,jj,kk) + z(ii,jj+1,kk) + z(ii+1,jj+1,kk)) };
                });
        }
        amrex::Real ground = amrex::get<0>(reduce_data.value(reduce_op));
        amrex::ParallelAllReduce::Sum(ground, amrex::ParallelContext::CommunicatorSub());
        return ground;
    }

    //! The ground under line iline, found once: level 0 and the terrain do not change
    amrex::Real
    line_ground (int iline, const amrex::Geometry& geom0, const amrex::MultiFab* ground_nd)
    {
        if (!m_z_agl) { return amrex::Real(0.0); }
        if (m_ground.empty()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ground_nd != nullptr,
                "erf.sample_line_z_agl needs the height of the ground on level 0");
            const int nline = static_cast<int>(m_bnd_rbx.size());
            m_ground.resize(nline);
            for (int l(0); l<nline; ++l) {
                const amrex::Box ibx = lineIndexBoxXY(m_bnd_rbx[l], geom0);
                m_ground[l] = column_ground(*ground_nd, geom0, ibx.smallEnd());
            }
        }
        return m_ground[iline];
    }

    /**
     * lineIndexBox, found once per level and kept until the grids of the level change: the
     * mesh does not move, so only which cells the level holds can change
     */
    const amrex::Box&
    cachedLineIndexBox (int iline, int lev, const amrex::Geometry& geom,
                        const amrex::MultiFab& z_nd, amrex::Real ground)
    {
        if (lev >= static_cast<int>(m_box_cache.size())) {
            m_box_cache.resize(lev+1);
            m_box_cache_found.resize(lev+1);
            m_box_cache_ba.resize(lev+1);
        }
        const int nline = static_cast<int>(m_bnd_rbx.size());
        if (static_cast<int>(m_box_cache[lev].size()) != nline || !(m_box_cache_ba[lev] == z_nd.boxArray())) {
            m_box_cache[lev].assign(nline, amrex::Box());
            m_box_cache_found[lev].assign(nline, 0);
            m_box_cache_ba[lev] = z_nd.boxArray();
        }
        if (m_box_cache_found[lev][iline] == 0) {
            m_box_cache[lev][iline] = lineIndexBox(iline, geom, z_nd, ground);
            m_box_cache_found[lev][iline] = 1;
        }
        return m_box_cache[lev][iline];
    }

    /**
     * The cells of a line given by real coordinates on a level, from the face heights of the
     * column through its first point.  A vertical line runs from the cell that holds the lower
     * z to the cell that holds the upper z (each clamped to the column); a horizontal line is
     * the row of cells at the level of cells that holds z in that column, so it follows the
     * terrain-fitted mesh.  With erf.sample_line_z_agl, ground is added to z.  The box is
     * empty if the level does not hold these cells of the column.
     */
    amrex::Box
    lineIndexBox (int iline, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                  amrex::Real ground)
    {
        const amrex::RealBox& rb = m_bnd_rbx[iline];
        const amrex::Box ibx = lineIndexBoxXY(rb, geom);
        const amrex::IntVect cell = ibx.smallEnd();
        const auto h = column_face_heights(z_nd, geom, cell);
        const auto& zl = h[0];
        const auto& zu = h[1];
        const auto& held = h[2];
        const int n = static_cast<int>(zl.size());
        const int kdlo = geom.Domain().smallEnd(2);
        const amrex::Real offset = m_z_agl ? ground : amrex::Real(0.0);

        // The cell (from the bottom) that holds height z, or -1 where the level does not hold it
        auto find_cell = [&] (amrex::Real z) -> int {
            if (held[0] > amrex::Real(0.5) && z < zl[0]) { return 0; }
            if (held[n-1] > amrex::Real(0.5) && z >= zu[n-1]) { return n-1; }
            for (int m = 0; m < n; ++m) {
                if (held[m] > amrex::Real(0.5) && zl[m] <= z && z < zu[m]) { return m; }
            }
            return -1;
        };
        const int mlo = find_cell(rb.lo(2) + offset);
        const int mhi = (m_dir[iline] == 2) ? find_cell(rb.hi(2) + offset) : mlo;
        if (mlo < 0 || mhi < 0 || mhi < mlo) { return amrex::Box(); }
        for (int m = mlo; m <= mhi; ++m) {
            if (held[m] < amrex::Real(0.5)) { return amrex::Box(); }
        }
        amrex::IntVect lo(ibx.smallEnd()); lo[2] = kdlo + mlo;
        amrex::IntVect hi(ibx.bigEnd());   hi[2] = kdlo + mhi;
        return amrex::Box(lo, hi);
    }

    //! The cells of a real box in x and y on a level, with k = 0: its k comes from heights
    static amrex::Box
    lineIndexBoxXY (const amrex::RealBox& rb, const amrex::Geometry& geom)
    {
        amrex::IntVect lo, hi;
        for (int d(0); d<2; ++d) {
            lo[d] = static_cast<int>(std::floor((rb.lo(d) - geom.ProbLo(d))/geom.CellSize(d)));
            hi[d] = static_cast<int>(std::floor((rb.hi(d) - geom.ProbLo(d))/geom.CellSize(d)));
            lo[d] = std::max(lo[d], geom.Domain().smallEnd(d));
            hi[d] = std::min(hi[d], geom.Domain().bigEnd(d));
        }
        lo[2] = hi[2] = geom.Domain().smallEnd(2);
        return amrex::Box(lo, hi);
    }

    amrex::Box
    getIndexBox (const amrex::RealBox& real_box,
                 const amrex::Geometry& geom) {
        amrex::IntVect slice_lo, slice_hi;

        AMREX_D_TERM(slice_lo[0]=static_cast<int>(std::floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0)));,
                     slice_lo[1]=static_cast<int>(std::floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1)));,
                     slice_lo[2]=static_cast<int>(std::floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2))););

        AMREX_D_TERM(slice_hi[0]=static_cast<int>(std::floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0)));,
                     slice_hi[1]=static_cast<int>(std::floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1)));,
                     slice_hi[2]=static_cast<int>(std::floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2))););

        return amrex::Box(slice_lo, slice_hi) & geom.Domain();
    }

    void
    write_coords (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& z_phys_cc,
                  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& z_phys_nd,
                  amrex::Vector<amrex::Geometry>& geom,
                  const amrex::MultiFab* ground_nd = nullptr,
                  const amrex::Vector<const amrex::MultiFab*>& ground_lev = {},
                  const amrex::Vector<const amrex::MultiFab*>& volfrac_lev = {})
    {
        if (!m_write_ascii) return;

        // Only the levels that are sampled
        int nlev = static_cast<int>(z_phys_cc.size());
        if (m_max_level >= 0) { nlev = std::min(nlev, m_max_level + 1); }
        for (int lev(0); lev < nlev; ++lev) {
            // A level allowed by amr.max_level that has not been made has no heights, and
            // neither does any level above it
            if (!z_phys_cc[lev] || lev >= static_cast<int>(z_phys_nd.size()) || !z_phys_nd[lev]) { break; }

            // Written for each set of grids of the level: when the level is made, and again
            // when a regrid changes its grids
            if (lev >= static_cast<int>(m_zcc_ba.size())) { m_zcc_ba.resize(lev+1); }
            if (m_zcc_ba[lev] == z_phys_cc[lev]->boxArray()) { continue; }
            m_zcc_ba[lev] = z_phys_cc[lev]->boxArray();

            amrex::Print() << "Writing out line coordinates of level " << lev << " to text" << std::endl;

            // Write one text file per level
            std::ofstream outfile;
            if (amrex::ParallelDescriptor::IOProcessor()) {
                std::string fname = amrex::Concatenate("plt_line_lev", lev, 1);
                fname += ".zcc";
                outfile.open(fname);

                if (!outfile.is_open()) {
                    amrex::AllPrint() << "Could not open " << fname << std::endl;
                }
            }

            // Loop over each line
            int nline = static_cast<int>(m_ls_mf.size());
            for (int iline(0); iline<nline; ++iline) {
                int dir = m_dir[iline];
                if (m_interp) {
                    // The height of each sample above sea level, where the level holds the line
                    // (the plan the samples are then taken with)
                    const InterpPlan& plan = cachedInterpPlan(iline, lev, geom[lev], geom[0], *z_phys_cc[lev],
                                                              z_phys_nd[lev].get(), level_ptr(volfrac_lev, lev),
                                                              level_ptr(ground_lev, lev), ground_nd);
                    if (amrex::ParallelDescriptor::IOProcessor()) {
                        if (plan.holds) {
                            for (const auto& zval : plan.zt) { outfile << " " << zval; }
                        }
                        outfile << std::endl;
                    }
                    continue;
                }
                amrex::Box bnd_bx;
                if (m_use_real_bx && use_column_heights()) {
                    // As in get_sample_data: the cells of the line come from the heights
                    const amrex::Real ground = line_ground(iline, geom[0], ground_nd);
                    bnd_bx = cachedLineIndexBox(iline, lev, geom[lev], *z_phys_nd[lev], ground);
                    if (bnd_bx.isEmpty()) {
                        // This level does not hold the line
                        if (amrex::ParallelDescriptor::IOProcessor()) { outfile << std::endl; }
                        continue;
                    }
                } else {
                    bnd_bx = (m_use_real_bx) ? getIndexBox(m_bnd_rbx[iline], geom[lev])
                                             : m_bnd_bx[iline];
                }
                amrex::IntVect first_cell = bnd_bx.smallEnd();

                // Create multifab with "sampled" z_phys values
                amrex::MultiFab line_coords_mf = get_line_data(
                    *z_phys_cc[lev], dir, first_cell, bnd_bx
                );

                // Convert multifab to vector
                amrex::Gpu::HostVector<amrex::Real> vec = sumToLine(
                    line_coords_mf, 0, 1, bnd_bx, dir
                );

                // Append to file
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    for (const auto& zval : vec) {
                        outfile << " " << zval;
                    }
                    outfile << std::endl;
                }
            } // line loop

            outfile.close();
        } // level loop
    }

    void
    get_sample_data (amrex::Vector<amrex::Geometry>& geom,
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new,
                     const amrex::Vector<amrex::Vector<amrex::MultiFab*>>& tau_lev = {},
                     const amrex::Vector<amrex::MultiFab*>& hfx3_lev = {},
                     const amrex::Vector<int>& prognostic_tke_available = {},
                     const amrex::Vector<amrex::MultiFab*>& z_cc_lev = {},
                     const amrex::Vector<amrex::MultiFab*>& z_nd_lev = {},
                     const amrex::MultiFab* ground_nd = nullptr,
                     const amrex::Vector<const amrex::MultiFab*>& volfrac_lev = {},
                     const amrex::Vector<const amrex::MultiFab*>& ground_lev = {})
    {
        if (m_ncons < 0) { expand_conserved(vars_new[0][Vars::cons].nComp()); }

        int nlev  = sampling_num_levels(vars_new);
        if (m_max_level >= 0) { nlev = std::min(nlev, m_max_level + 1); }
        int nline = static_cast<int>(m_bnd_bx.size());
        int ncomp = static_cast<int>(m_varnames.size());

        int qv_comp = -1;

        // The stress tau_type on level ilev, or nullptr where it is not stored
        auto get_tau = [&] (int ilev, int tau_type) -> amrex::MultiFab* {
            if (ilev >= static_cast<int>(tau_lev.size()) ||
                tau_type >= static_cast<int>(tau_lev[ilev].size())) { return nullptr; }
            return tau_lev[ilev][tau_type];
        };

        // For a line given by real coordinates, its cells come from the heights of the column
        // through its first point, and a height above ground from the ground under it
        const bool by_height = m_use_real_bx && use_column_heights();
        if (by_height) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!z_nd_lev.empty() && z_nd_lev[0] != nullptr,
                "Line sampling by height needs the nodal heights of each level");
        }
        // An interpolated line is sampled on the finest level that holds every cell it reads
        auto interp_plan = [&] (int iline, int ilev) -> const InterpPlan& {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ilev < static_cast<int>(z_cc_lev.size()) && z_cc_lev[ilev] != nullptr,
                "Line sampling with erf.sample_line_interp needs the cell-centre heights of each level");
            const amrex::MultiFab* z_nd = ilev < static_cast<int>(z_nd_lev.size()) ? z_nd_lev[ilev] : nullptr;
            return cachedInterpPlan(iline, ilev, geom[ilev], geom[0], *z_cc_lev[ilev], z_nd,
                                    level_ptr(volfrac_lev, ilev), level_ptr(ground_lev, ilev), ground_nd);
        };
        auto line_box = [&] (int iline, int ilev) -> amrex::Box {
            if (m_interp) {
                const InterpPlan& plan = interp_plan(iline, ilev);
                return plan.holds ? plan.cols : amrex::Box();
            }
            if (by_height) {
                if (ilev >= static_cast<int>(z_nd_lev.size()) || z_nd_lev[ilev] == nullptr) { return amrex::Box(); }
                return cachedLineIndexBox(iline, ilev, geom[ilev], *z_nd_lev[ilev],
                                          line_ground(iline, geom[0], ground_nd));
            }
            return (m_use_real_bx) ? getIndexBox(m_bnd_rbx[iline], geom[ilev]) : m_bnd_bx[iline];
        };

        // Search each level, from the finest, to get the finest data possible.  The data of a
        // level are built once for all the lines that are sampled there
        amrex::Vector<int> done(nline, 0);
        for (int ilev(nlev-1); ilev>=0; --ilev) {
            amrex::Vector<int> lines;
            amrex::Vector<amrex::Box> boxes;
            for (int iline(0); iline<nline; ++iline) {
                if (done[iline] != 0) { continue; }
                const amrex::Box bnd_bx = line_box(iline, ilev);
                if (bnd_bx.isEmpty()) { continue; }
                lines.push_back(iline);
                boxes.push_back(bnd_bx);
            }
            if (lines.empty()) { continue; }

            // Construct CC velocities
            amrex::MultiFab mf_cc_vel;
            auto ba = vars_new[ilev][Vars::cons].boxArray();
            auto dm = vars_new[ilev][Vars::cons].DistributionMap();
            mf_cc_vel.define(ba, dm, AMREX_SPACEDIM, amrex::IntVect(1,1,1));
            average_face_to_cellcenter(mf_cc_vel,0,
                                       amrex::Array<const amrex::MultiFab*,3>{&vars_new[ilev][Vars::xvel],
                                                                              &vars_new[ilev][Vars::yvel],
                                                                              &vars_new[ilev][Vars::zvel]});

            // Construct MultiFab holding requested variables
            amrex::MultiFab mf_cc_data;
            mf_cc_data.define(ba, dm, ncomp, 1);

            int mf_comp = 0;
            int z_comp  = -1;
            int vf_comp = -1;

            // The raw conserved components, cons_0, cons_1, ..., come first
            if (m_ncons > 0) {
                amrex::MultiFab::Copy(mf_cc_data, vars_new[ilev][Vars::cons], 0, mf_comp, m_ncons, 0);
                mf_comp += m_ncons;
            }

            // The height of each cell centre
            if (containerHasElement(m_varnames, std::string("z"))) {
                amrex::MultiFab* z_cc = ilev < static_cast<int>(z_cc_lev.size()) ? z_cc_lev[ilev] : nullptr;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(z_cc != nullptr,
                    "Line sampling of z needs the cell-centre heights of each level");
                amrex::MultiFab::Copy(mf_cc_data, *z_cc, 0, mf_comp, 1, 0);
                z_comp = mf_comp;
                mf_comp += 1;
            }

            // The volume fraction of each cell: 1 without immersed terrain, 0 where covered
            const amrex::MultiFab* volfrac = ilev < static_cast<int>(volfrac_lev.size())
                ? volfrac_lev[ilev] : nullptr;
            if (containerHasElement(m_varnames, std::string("volfrac"))) {
                if (volfrac != nullptr) {
                    amrex::MultiFab::Copy(mf_cc_data, *volfrac, 0, mf_comp, 1, 0);
                } else {
                    mf_cc_data.setVal(1.0, mf_comp, 1, 0);
                }
                vf_comp = mf_comp;
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "density")) {
                amrex::MultiFab::Copy(mf_cc_data, vars_new[ilev][Vars::cons], Rho_comp, mf_comp, 1, 0);
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "x_velocity")) {
                amrex::MultiFab::Copy(mf_cc_data, mf_cc_vel, 0, mf_comp, 1, 0);
                mf_comp += 1;
            }
            if (containerHasElement(m_varnames, "y_velocity")) {
                amrex::MultiFab::Copy(mf_cc_data, mf_cc_vel, 1, mf_comp, 1, 0);
                mf_comp += 1;
            }
            if (containerHasElement(m_varnames, "z_velocity")) {
                amrex::MultiFab::Copy(mf_cc_data, mf_cc_vel, 2, mf_comp, 1, 0);
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "magvel")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const amrex::Box& tbx = mfi.tilebox();
                    auto const& dfab = mf_cc_data.array(mfi);
                    auto const& vfab = mf_cc_vel.array(mfi);

                    amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        dfab(i,j,k,mf_comp) = std::sqrt(vfab(i,j,k,0)*vfab(i,j,k,0)
                                                      + vfab(i,j,k,1)*vfab(i,j,k,1)
                                                      + vfab(i,j,k,2)*vfab(i,j,k,2)) ;
                    });
                }
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "theta")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const amrex::Box& tbx = mfi.tilebox();
                    auto const& dfab = mf_cc_data.array(mfi);
                    auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                    amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        dfab(i,j,k,mf_comp) = cfab(i,j,k,RhoTheta_comp) / cfab(i,j,k,Rho_comp);
                    });
                }
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "qv")) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vars_new[ilev][Vars::cons].nComp() > RhoQ1_comp,
                    "qv sampling requested but moisture components not present in state");
                if (qv_comp >= 0) AMREX_ALWAYS_ASSERT(qv_comp == mf_comp);
                qv_comp = mf_comp;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const amrex::Box& tbx = mfi.tilebox();
                    auto const& dfab = mf_cc_data.array(mfi);
                    auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                    amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        dfab(i,j,k,mf_comp) = cfab(i,j,k,RhoQ1_comp) / cfab(i,j,k,Rho_comp);
                    });
                }
                mf_comp += 1;
            }
            if (containerHasElement(m_varnames, "qc")) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vars_new[ilev][Vars::cons].nComp() > RhoQ2_comp,
                    "qc sampling requested but moisture components not present in state");
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const amrex::Box& tbx = mfi.tilebox();
                    auto const& dfab = mf_cc_data.array(mfi);
                    auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                    amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                    {
                        dfab(i,j,k,mf_comp) = cfab(i,j,k,RhoQ2_comp) / cfab(i,j,k,Rho_comp);
                    });
                }
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "pressure")) {
                if (vars_new[ilev][Vars::cons].nComp() > RhoQ1_comp) {
                    // moist pressure: use qv from dfab if already sampled, else compute inline
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        auto const& dfab = mf_cc_data.array(mfi);
                        auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            amrex::Real qv_val = (qv_comp >= 0) ? dfab(i,j,k,qv_comp)
                                                                 : cfab(i,j,k,RhoQ1_comp) / cfab(i,j,k,Rho_comp);
                            dfab(i,j,k,mf_comp) = getPgivenRTh(cfab(i,j,k,RhoTheta_comp), qv_val);
                        });
                    }
                } else {
                    // dry pressure
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        auto const& dfab = mf_cc_data.array(mfi);
                        auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            dfab(i,j,k,mf_comp) = getPgivenRTh(cfab(i,j,k,RhoTheta_comp));
                        });
                    }
                }
                mf_comp += 1;
            }

            if (containerHasElement(m_varnames, "sgs_tke")) {
                const bool tke_available =
                    ilev < static_cast<int>(prognostic_tke_available.size()) &&
                    prognostic_tke_available[ilev] != 0;
                if (!tke_available) {
                    mf_cc_data.setVal(amrex::Real(0.0), mf_comp, 1, 0);
                } else {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const amrex::Box& tbx = mfi.tilebox();
                    const auto& out = mf_cc_data.array(mfi);
                    const auto& state = vars_new[ilev][Vars::cons].const_array(mfi);
                    amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                        out(i,j,k,mf_comp) = state(i,j,k,RhoKE_comp) / state(i,j,k,Rho_comp);
                        });
                    }
                }
                mf_comp += 1;
            }

            // The stresses at cell centres: tau11, tau22 and tau33 are stored there, and the
            // others are averaged from the four edges around the cell
            auto sample_tau = [&] (const char* name, int tau_type,
                                   amrex::IntVect const& e1, amrex::IntVect const& e2) {
                if (!containerHasElement(m_varnames, std::string(name))) { return; }
                amrex::MultiFab* tau = get_tau(ilev, tau_type);
                if (tau != nullptr) {
                    const bool at_centre = (e1 == amrex::IntVect(0));
                    const int comp = mf_comp;
                    const int ax = e1[0], ay = e1[1], az = e1[2];
                    const int bx = e2[0], by = e2[1], bz = e2[2];
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        const auto& out = mf_cc_data.array(mfi);
                        const auto& t = tau->const_array(mfi);
                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            out(i,j,k,comp) = at_centre ? t(i,j,k) : amrex::Real(0.25) *
                                (t(i,j,k) + t(i+ax,j+ay,k+az) +
                                 t(i+bx,j+by,k+bz) + t(i+ax+bx,j+ay+by,k+az+bz));
                        });
                    }
                } else {
                    mf_cc_data.setVal(0.0, mf_comp, 1, 0);
                }
                mf_comp += 1;
            };
            const amrex::IntVect ex(1,0,0), ey(0,1,0), ez(0,0,1), e0(0);
            sample_tau("sgs_tau11", TauType::tau11, e0, e0);
            sample_tau("sgs_tau12", TauType::tau12, ex, ey);
            sample_tau("sgs_tau13", TauType::tau13, ex, ez);
            sample_tau("sgs_tau22", TauType::tau22, e0, e0);
            sample_tau("sgs_tau23", TauType::tau23, ey, ez);
            sample_tau("sgs_tau33", TauType::tau33, e0, e0);

            if (containerHasElement(m_varnames, "sgs_hfx3")) {
                amrex::MultiFab* hfx3 = ilev < static_cast<int>(hfx3_lev.size())
                    ? hfx3_lev[ilev] : nullptr;
                if (hfx3 != nullptr) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        const auto& out = mf_cc_data.array(mfi);
                        const auto& hfx = hfx3->const_array(mfi);
                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            out(i,j,k,mf_comp) = amrex::Real(0.5) * (hfx(i,j,k) + hfx(i,j,k+1));
                        });
                    }
                } else {
                    mf_cc_data.setVal(0.0, mf_comp, 1, 0);
                }
                mf_comp += 1;
            }

            // A cell covered by immersed terrain holds no fluid, so every variable but its
            // height and volume fraction is written there as erf.sampling_covered_value
            if (volfrac != nullptr) {
                const amrex::Real covered_value = m_covered_value;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const amrex::Box& tbx = mfi.tilebox();
                    const auto& out = mf_cc_data.array(mfi);
                    const auto& vf  = volfrac->const_array(mfi);
                    amrex::ParallelFor(tbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                        if (vf(i,j,k) == amrex::Real(0.0) && n != z_comp && n != vf_comp) {
                            out(i,j,k,n) = covered_value;
                        }
                    });
                }
            }

            for (int n(0); n<static_cast<int>(lines.size()); ++n) {
                const int iline = lines[n];
                const amrex::Box& bnd_bx = boxes[n];
                m_lev[iline] = ilev;
                if (m_interp) {
                    m_ls_mf[iline] = interp_line_data(iline, interp_plan(iline, ilev), mf_cc_data, z_comp, vf_comp);
                    done[iline] = 1;
                    continue;
                }
                m_ls_mf[iline] = get_line_data(mf_cc_data, m_dir[iline], bnd_bx.smallEnd(), bnd_bx);

                // A line is done once a level holds all of it
                if (bnd_bx == m_ls_mf[iline].boxArray().minimalBox()) { done[iline] = 1; }
            }
        } // ilev

        // Fill bnd_bx if use_rbx
        if (m_use_real_bx) {
            for (int iline(0); iline<nline; ++iline) {
                m_bnd_bx[iline] = m_ls_mf[iline].boxArray().minimalBox();
            }
        }
    }

    void
    write_sample_data (amrex::Vector<double>& time,
                       amrex::Vector<int>& level_steps,
                       amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Geometry>& geom)
    {
        if (m_write_log) {
            write_line_log(time);
        } else if (m_write_ascii) {
            write_line_ascii(time);
        } else {
            write_line_plotfile(time, level_steps, ref_ratio, geom);
        }
    }

    //! One row per line in its log: the time, then each variable along the line
    void
    write_line_log (amrex::Vector<double>& time)
    {
        // same as definitions in ERF.H
        constexpr int datwidth = 14;
        constexpr int datprecision = 6;

        int nline = static_cast<int>(m_ls_mf.size());
        int nvar  = static_cast<int>(m_varnames.size());
        for (int iline(0); iline<nline; ++iline) {
            int dir = m_dir[iline];
            double m_time = time[m_lev[iline]];
            amrex::Box m_dom = m_bnd_bx[iline];

            amrex::Vector<amrex::Gpu::HostVector<amrex::Real>> vecs(nvar);
            for (int ivar(0); ivar<nvar; ++ivar) {
                vecs[ivar] = sumToLine(m_ls_mf[iline], ivar, 1, m_dom, dir);
            }

            if (amrex::ParallelDescriptor::IOProcessor()) {
                std::ostream& fs = *m_logstream[iline];
                if (fs.good()) {
                    fs << std::setw(datwidth) << std::setprecision(datprecision) << m_time;
                    for (const auto& vec : vecs) {
                        for (const auto& val : vec) {
                            fs << std::setw(datwidth) << std::setprecision(datprecision) << val;
                        }
                    }
                    fs << std::endl;
                }
            }
        }
    }

    void
    write_line_ascii (amrex::Vector<double>& time)
    {
        // same as definitions in ERF.H
        constexpr int datwidth = 14;
        constexpr int datprecision = 6;
        constexpr int timeprecision = 13;

        open_ascii_streams();

        int nline = static_cast<int>(m_ls_mf.size());
        int nvar  = static_cast<int>(m_varnames.size());
        for (int iline(0); iline<nline; ++iline) {
            int dir = m_dir[iline];
            int lev = m_lev[iline];
            double m_time = time[lev];
            amrex::Box m_dom   = m_bnd_bx[iline];

            for (int ivar(0); ivar<nvar; ++ivar) {
                // Convert multifab to vector
                amrex::Gpu::HostVector<amrex::Real> vec = sumToLine(m_ls_mf[iline], ivar, 1, m_dom, dir);

                if (amrex::ParallelDescriptor::IOProcessor()) {
                    int ifile = iline*nvar + ivar;
                    std::ostream& fs = *m_datastream[ifile];
                    fs << std::setw(datwidth) << std::setprecision(timeprecision) << m_time
                       << std::setw(datwidth) << std::setprecision(datprecision);
                    for (const auto& val : vec) {
                        fs << " " << val;
                    }
                    fs << std::endl;
                }
            }
        }
    }

    void
    write_line_plotfile (amrex::Vector<double>& time,
                         amrex::Vector<int>& level_steps,
                         amrex::Vector<amrex::IntVect>& ref_ratio,
                         amrex::Vector<amrex::Geometry>& geom)
    {
        int nline = static_cast<int>(m_ls_mf.size());
        for (int iline(0); iline<nline; ++iline) {
            // Data members that can be used as-is
            int dir = m_dir[iline];
            int lev = m_lev[iline];
            double m_time = time[lev];
            amrex::Vector<int> m_level_steps = {level_steps[lev]};
            amrex::Vector<amrex::IntVect> m_ref_ratio  = {ref_ratio[lev]};

            // Create modified geometry object corresponding to the line
            auto plo = geom[lev].ProbLo();
            auto dx  = geom[lev].CellSize();
            amrex::Vector<amrex::Geometry> m_geom; m_geom.resize(1);
            amrex::Vector<int> is_per(AMREX_SPACEDIM,0);
            amrex::Box m_dom = m_bnd_bx[iline];
            amrex::RealBox m_rb;
            for (int d(0); d<AMREX_SPACEDIM; ++d) {
                amrex::Real offset = (d==dir) ? 0 : myhalf;
                amrex::Real lo = plo[d] + ( m_dom.smallEnd(d) - offset ) * dx[d];
                amrex::Real hi = plo[d] + ( m_dom.bigEnd(d)   + offset ) * dx[d];

                m_rb.setLo(d,lo);
                m_rb.setHi(d,hi);

                is_per[d] = geom[lev].isPeriodic(d);
            }
            // An interpolated line is at its requested heights (above the ground with
            // erf.sample_line_z_agl): its z cells are centred on them
            if (m_interp) {
                const amrex::Vector<amrex::Real> h = interp_heights(iline);
                const int nh = static_cast<int>(h.size());
                const amrex::Real dh = (nh > 1) ? (h[nh-1] - h[0]) / amrex::Real(nh - 1) : dx[2];
                m_rb.setLo(2, h[0] - myhalf * dh);
                m_rb.setHi(2, h[nh-1] + myhalf * dh);
                is_per[2] = 0;
            }
            m_geom[0].define(m_dom, &m_rb, geom[lev].Coord(), is_per.data());

            // Create plotfile name
            std::string name_line = m_name[iline];
            name_line += "_step_";
            std::string plotfilename = amrex::Concatenate(name_line, m_level_steps[0], 5);

            // Get the data
            amrex::Vector<const amrex::MultiFab*> mf = {&(m_ls_mf[iline])};

            // Write each line
            WriteMultiLevelPlotfile(plotfilename, 1, mf,
                                    m_varnames, m_geom, static_cast<amrex::Real>(m_time),
                                    m_level_steps, m_ref_ratio);
        }
    }

    amrex::Vector<int> m_dir;
    amrex::Vector<int> m_lev;
    amrex::Vector<amrex::Box> m_bnd_bx;
    amrex::Vector<amrex::RealBox> m_bnd_rbx;
    amrex::Vector<amrex::MultiFab> m_ls_mf;
    amrex::Vector<std::string> m_name;

    bool m_use_real_bx{false};
    bool m_use_latlon{false};   // erf.sample_line_lo/hi_latlon
    amrex::Vector<amrex::Real> m_latlon_lo, m_latlon_hi; // lat, long, z of each line's ends
    bool m_write_ascii{false};
    bool m_write_log{false};
    bool m_z_agl{false};   // erf.sample_line_z_agl
    bool m_interp{false};  // erf.sample_line_interp
    amrex::Vector<int> m_npts;   // erf.sample_line_npts: the heights of each vertical line (0 for a horizontal one)
    amrex::Real m_covered_value{sampling_covered_value()}; // erf.sampling_covered_value
    int  m_max_level{-1};  // finest level to sample; -1 for all
    int  m_ncons{-1};      // components "conserved" expands into; -1 until the first sample
    amrex::Vector<std::string> m_varnames {"magvel","theta"};
    amrex::Vector<std::unique_ptr<std::fstream> > m_datastream;
    amrex::Vector<std::unique_ptr<std::fstream> > m_logstream;

    // The ground under each line, and the cells of each line on each level with the grids
    // they were found on (lines given by real coordinates on a stretched or terrain-fitted mesh)
    amrex::Vector<amrex::Real> m_ground;
    amrex::Vector<amrex::Vector<amrex::Box>> m_box_cache;
    amrex::Vector<amrex::Vector<int>> m_box_cache_found;
    amrex::Vector<amrex::BoxArray> m_box_cache_ba;
    // The grids of each level when its line coordinates were last written (write_coords)
    amrex::Vector<amrex::BoxArray> m_zcc_ba;
    // The samples of each interpolated line on each level (erf.sample_line_interp)
    amrex::Vector<amrex::Vector<InterpPlan>> m_interp_plan;
};


struct PlaneSampler
{
    PlaneSampler ()
    {
        amrex::ParmParse pp("erf");

        // Count number of lo and hi points define the plane
        int n_plane_lo  = pp.countval("sample_plane_lo") / AMREX_SPACEDIM;
        int n_plane_hi  = pp.countval("sample_plane_hi") / AMREX_SPACEDIM;
        int n_plane_dir = pp.countval("sample_plane_dir");
        AMREX_ALWAYS_ASSERT( (n_plane_lo==n_plane_hi ) &&
                             (n_plane_lo==n_plane_dir) );

        // Parse the data
        if (n_plane_lo > 0) {
            // Parse lo
            amrex::Vector<amrex::Real> r_lo; r_lo.resize(n_plane_lo*AMREX_SPACEDIM);
            amrex::Vector<amrex::Vector<amrex::Real>> rv_lo;
            pp.queryarr("sample_plane_lo",r_lo,0,n_plane_lo*AMREX_SPACEDIM);
            for (int i(0); i < n_plane_lo; i++) {
                amrex::Vector<amrex::Real> rv = {r_lo[AMREX_SPACEDIM*i+0],
                                                 r_lo[AMREX_SPACEDIM*i+1],
                                                 r_lo[AMREX_SPACEDIM*i+2]};
                rv_lo.push_back(rv);
            }

            // Parse hi
            amrex::Vector<amrex::Real> r_hi; r_hi.resize(n_plane_hi*AMREX_SPACEDIM);
            amrex::Vector<amrex::Vector<amrex::Real>> rv_hi;
            pp.queryarr("sample_plane_hi",r_hi,0,n_plane_hi*AMREX_SPACEDIM);
            for (int i(0); i < n_plane_hi; i++) {
                amrex::Vector<amrex::Real> rv = {r_hi[AMREX_SPACEDIM*i+0],
                                                 r_hi[AMREX_SPACEDIM*i+1],
                                                 r_hi[AMREX_SPACEDIM*i+2]};
                rv_hi.push_back(rv);
            }

            // Construct vector of bounding real boxes
            m_bnd_rbx.resize(n_plane_lo);
            for (int i(0); i < n_plane_hi; i++){
                amrex::RealBox rbx(rv_lo[i].data(),rv_hi[i].data());
                m_bnd_rbx[i] = rbx;
            }

            // Parse directionality
            m_dir.resize(n_plane_dir);
            pp.queryarr("sample_plane_dir",m_dir,0,n_plane_dir);

            // Parse names
            std::string name_base = "plt_plane_";
            m_name.resize(n_plane_lo);
            int n_names = pp.countval("sample_plane_name");
            if (n_names > 0) {
                AMREX_ALWAYS_ASSERT( n_names==n_plane_lo );
                pp.queryarr("sample_plane_name",m_name,0,n_names);
            } else {
                for (int iplane(0); iplane<n_plane_lo; ++iplane) {
                    m_name[iplane] = amrex::Concatenate(name_base, iplane , 5);
                }
            }

            // Optional cap on the finest level written (-1 => all intersecting levels)
            pp.queryAdd("plane_sampling_max_level", m_max_level);

            // Allocate per-plane vectors of per-level MF pointers
            m_ps_mf.resize(n_plane_lo);

            // Get requested vars
            if (pp.countval("plane_sampling_vars") > 0) {
                m_varnames.clear();
                amrex::Vector<std::string> requested_vars;
                pp.queryarr("plane_sampling_vars",requested_vars);
                amrex::Print() << "Selected plane sampling vars :";
                if (containerHasElement(requested_vars, "density")) {
                    m_varnames.push_back("density");
                    amrex::Print() << " " << "density";
                }
                if (containerHasElement(requested_vars, "x_velocity")) {
                    m_varnames.push_back("x_velocity");
                    amrex::Print() << " " << "x_velocity";
                }
                if (containerHasElement(requested_vars, "y_velocity")) {
                    m_varnames.push_back("y_velocity");
                    amrex::Print() << " " << "y_velocity";
                }
                if (containerHasElement(requested_vars, "z_velocity")) {
                    m_varnames.push_back("z_velocity");
                    amrex::Print() << " " << "z_velocity";
                }
                if (containerHasElement(requested_vars, "magvel")) {
                    m_varnames.push_back("magvel");
                    amrex::Print() << " " << "magvel";
                }
                if (containerHasElement(requested_vars, "theta")) {
                    m_varnames.push_back("theta");
                    amrex::Print() << " " << "theta";
                }
                if (containerHasElement(requested_vars, "qv")) {
                    m_varnames.push_back("qv");
                    amrex::Print() << " " << "qv";
                }
                if (containerHasElement(requested_vars, "qc")) {
                    m_varnames.push_back("qc");
                    amrex::Print() << " " << "qc";
                }
                if (containerHasElement(requested_vars, "pressure")) {
                    m_varnames.push_back("pressure");
                    amrex::Print() << " " << "pressure";
                }
                if (containerHasElement(requested_vars, "volfrac")) {
                    m_varnames.push_back("volfrac");
                    amrex::Print() << " " << "volfrac";
                }
                amrex::Print() << std::endl;
            }
        }
    }

    // This must match what is in AMReX_MultiFabUtil.H
    amrex::Box
    getIndexBox (const amrex::RealBox& real_box,
                 const amrex::Geometry& geom) {
        amrex::IntVect slice_lo, slice_hi;

        AMREX_D_TERM(slice_lo[0]=static_cast<int>(std::floor((real_box.lo(0) - geom.ProbLo(0))/geom.CellSize(0)));,
                     slice_lo[1]=static_cast<int>(std::floor((real_box.lo(1) - geom.ProbLo(1))/geom.CellSize(1)));,
                     slice_lo[2]=static_cast<int>(std::floor((real_box.lo(2) - geom.ProbLo(2))/geom.CellSize(2))););

        AMREX_D_TERM(slice_hi[0]=static_cast<int>(std::floor((real_box.hi(0) - geom.ProbLo(0))/geom.CellSize(0)));,
                     slice_hi[1]=static_cast<int>(std::floor((real_box.hi(1) - geom.ProbLo(1))/geom.CellSize(1)));,
                     slice_hi[2]=static_cast<int>(std::floor((real_box.hi(2) - geom.ProbLo(2))/geom.CellSize(2))););

        return amrex::Box(slice_lo, slice_hi) & geom.Domain();
    }

    void
    get_sample_data (amrex::Vector<amrex::Geometry>& geom,
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new,
                     const amrex::Vector<const amrex::MultiFab*>& volfrac_lev = {})
    {
        int nlev   = sampling_num_levels(vars_new);
        int nplane = static_cast<int>(m_bnd_rbx.size());
        int ncomp  = static_cast<int>(m_varnames.size());
        bool interpolate = true;

        int qv_comp = -1;

        // Loop over each plane
        for (int iplane(0); iplane<nplane; ++iplane) {
            int dir = m_dir[iplane];
            amrex::RealBox bnd_rbx  = m_bnd_rbx[iplane];
            amrex::Real point = bnd_rbx.lo(dir);

            // Collect every level 0..lev_cap whose grids intersect the plane.
            int lev_cap = (m_max_level >= 0 && m_max_level < nlev-1) ? m_max_level : nlev-1;
            m_ps_mf[iplane].clear();

            for (int ilev(0); ilev<=lev_cap; ++ilev) {

                // Stop ascending once a level's grids no longer touch the flattened plane
                {
                    amrex::Box plane_bx = getIndexBox(bnd_rbx, geom[ilev]);
                    int k_l = static_cast<int>(std::floor((point - geom[ilev].ProbLo(dir))
                                                          / geom[ilev].CellSize(dir)));
                    plane_bx.setSmall(dir, k_l); plane_bx.setBig(dir, k_l);
                    if (!vars_new[ilev][Vars::cons].boxArray().intersects(plane_bx)) { break; }
                }

                // Construct CC velocities
                amrex::MultiFab mf_cc_vel;
                auto ba = vars_new[ilev][Vars::cons].boxArray();
                auto dm = vars_new[ilev][Vars::cons].DistributionMap();
                mf_cc_vel.define(ba, dm, AMREX_SPACEDIM, amrex::IntVect(1,1,1));
                average_face_to_cellcenter(mf_cc_vel,0,
                                           amrex::Array<const amrex::MultiFab*,3>{&vars_new[ilev][Vars::xvel],
                                                                                  &vars_new[ilev][Vars::yvel],
                                                                                  &vars_new[ilev][Vars::zvel]});

                // Construct MultiFab holding requested variables, and after them whether
                // each cell is covered by immersed terrain
                amrex::MultiFab mf_cc_data;
                mf_cc_data.define(ba, dm, ncomp+1, 1);

                int mf_comp = 0;

                if (containerHasElement(m_varnames, "density")) {
                    amrex::MultiFab::Copy(mf_cc_data, vars_new[ilev][Vars::cons], Rho_comp, mf_comp, 1, 0);
                    mf_comp += 1;
                }

                if (containerHasElement(m_varnames, "x_velocity")) {
                    amrex::MultiFab::Copy(mf_cc_data, mf_cc_vel, 0, mf_comp, 1, 0);
                    mf_comp += 1;
                }
                if (containerHasElement(m_varnames, "y_velocity")) {
                    amrex::MultiFab::Copy(mf_cc_data, mf_cc_vel, 1, mf_comp, 1, 0);
                    mf_comp += 1;
                }
                if (containerHasElement(m_varnames, "z_velocity")) {
                    amrex::MultiFab::Copy(mf_cc_data, mf_cc_vel, 2, mf_comp, 1, 0);
                    mf_comp += 1;
                }

                if (containerHasElement(m_varnames, "magvel")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        auto const& dfab = mf_cc_data.array(mfi);
                        auto const& vfab = mf_cc_vel.array(mfi);

                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            dfab(i,j,k,mf_comp) = std::sqrt(vfab(i,j,k,0)*vfab(i,j,k,0)
                                                           + vfab(i,j,k,1)*vfab(i,j,k,1)
                                                           + vfab(i,j,k,2)*vfab(i,j,k,2));
                        });
                    }
                    mf_comp += 1;
                }

                if (containerHasElement(m_varnames, "theta")) {
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        auto const& dfab = mf_cc_data.array(mfi);
                        auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            dfab(i,j,k,mf_comp) = cfab(i,j,k,RhoTheta_comp) / cfab(i,j,k,Rho_comp);
                        });
                    }
                    mf_comp += 1;
                }

                if (containerHasElement(m_varnames, "qv")) {
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vars_new[ilev][Vars::cons].nComp() > RhoQ1_comp,
                        "qv sampling requested but moisture components not present in state");
                    if (qv_comp >= 0) AMREX_ALWAYS_ASSERT(qv_comp == mf_comp);
                    qv_comp = mf_comp;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        auto const& dfab = mf_cc_data.array(mfi);
                        auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            dfab(i,j,k,mf_comp) = cfab(i,j,k,RhoQ1_comp) / cfab(i,j,k,Rho_comp);
                        });
                    }
                    mf_comp += 1;
                }
                if (containerHasElement(m_varnames, "qc")) {
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vars_new[ilev][Vars::cons].nComp() > RhoQ2_comp,
                        "qc sampling requested but moisture components not present in state");
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& tbx = mfi.tilebox();
                        auto const& dfab = mf_cc_data.array(mfi);
                        auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                        amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                            dfab(i,j,k,mf_comp) = cfab(i,j,k,RhoQ2_comp) / cfab(i,j,k,Rho_comp);
                        });
                    }
                    mf_comp += 1;
                }

                if (containerHasElement(m_varnames, "pressure")) {
                    if (vars_new[ilev][Vars::cons].nComp() > RhoQ1_comp) {
                        // moist pressure: use qv from dfab if already sampled, else compute inline
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                        for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                            const amrex::Box& tbx = mfi.tilebox();
                            auto const& dfab = mf_cc_data.array(mfi);
                            auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                            {
                                amrex::Real qv_val = (qv_comp >= 0) ? dfab(i,j,k,qv_comp)
                                                                     : cfab(i,j,k,RhoQ1_comp) / cfab(i,j,k,Rho_comp);
                                dfab(i,j,k,mf_comp) = getPgivenRTh(cfab(i,j,k,RhoTheta_comp), qv_val);
                            });
                        }
                    } else {
                        // dry pressure
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                        for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                            const amrex::Box& tbx = mfi.tilebox();
                            auto const& dfab = mf_cc_data.array(mfi);
                            auto const& cfab = vars_new[ilev][Vars::cons].array(mfi);

                            amrex::ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                            {
                                dfab(i,j,k,mf_comp) = getPgivenRTh(cfab(i,j,k,RhoTheta_comp));
                            });
                        }
                    }
                    mf_comp += 1;
                }

                // The volume fraction, 1 without immersed terrain, and 1 where a cell is
                // covered (0 elsewhere)
                const amrex::MultiFab* volfrac = ilev < static_cast<int>(volfrac_lev.size())
                    ? volfrac_lev[ilev] : nullptr;
                int vf_comp = -1;
                if (containerHasElement(m_varnames, "volfrac")) {
                    if (volfrac != nullptr) {
                        amrex::MultiFab::Copy(mf_cc_data, *volfrac, 0, mf_comp, 1, 1);
                    } else {
                        mf_cc_data.setVal(1.0, mf_comp, 1, 1);
                    }
                    vf_comp = mf_comp;
                    mf_comp += 1;
                }
                if (volfrac != nullptr) {
                    const int cov_comp = ncomp;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                    for (amrex::MFIter mfi(mf_cc_data, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                        const amrex::Box& gbx = mfi.growntilebox(1);
                        const auto& out = mf_cc_data.array(mfi);
                        const auto& vf  = volfrac->const_array(mfi);
                        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                            out(i,j,k,cov_comp) = (vf(i,j,k) == amrex::Real(0.0)) ? amrex::Real(1.0)
                                                                                  : amrex::Real(0.0);
                        });
                    }
                } else {
                    mf_cc_data.setVal(0.0, ncomp, 1, 1);
                }

                // The interpolation reads the cell on the other side of the plane, which may be
                // a ghost cell.  It is taken from the neighbouring grid, or periodic image, where
                // there is one; at the domain boundary or the edge of the level it is the cell
                // itself, so the plane there takes the value of the cell it lies in
                {
                    const int fdir = dir;
                    for (amrex::MFIter mfi(mf_cc_data); mfi.isValid(); ++mfi) {
                        const amrex::Box& vbx = mfi.validbox();
                        const int vlo = vbx.smallEnd(fdir);
                        const int vhi = vbx.bigEnd(fdir);
                        auto const& dfab = mf_cc_data.array(mfi);
                        amrex::ParallelFor(amrex::grow(vbx, fdir, 1), ncomp+1,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
                        {
                            amrex::IntVect iv(i,j,k);
                            if (iv[fdir] < vlo || iv[fdir] > vhi) {
                                iv[fdir] = amrex::min(amrex::max(iv[fdir], vlo), vhi);
                                dfab(i,j,k,n) = dfab(iv,n);
                            }
                        });
                    }
                    mf_cc_data.FillBoundary(geom[ilev].periodicity());
                }

                auto slice = get_slice_data(dir, point, mf_cc_data, geom[ilev],
                                            0, ncomp+1, interpolate, bnd_rbx);

                // Defensive: an empty return means nothing was actually sampled.
                if (!slice || slice->boxArray().size() == 0) { break; }

                // Broadcast the single sampled plane across 2^ilev cells in dir so
                // every level shares an isotropic refinement ratio. The dir extent
                // becomes [0, 2^ilev-1]; each cell holds the interpolated plane value
                // replicated across the shared level-0 slab thickness.
                const int rr_l    = 1 << ilev;
                const int k_slice = slice->boxArray().minimalBox().smallEnd(dir);
                const int l_dir   = dir;
                const amrex::Real covered_value = m_covered_value;

                amrex::BoxList bl;
                const amrex::BoxArray& slice_ba = slice->boxArray();
                for (int ib(0); ib<static_cast<int>(slice_ba.size()); ++ib) {
                    amrex::Box b = slice_ba[ib];
                    b.setSmall(dir, 0);
                    b.setBig(dir, rr_l - 1);
                    bl.push_back(b);
                }
                amrex::BoxArray out_ba(std::move(bl));
                auto out_mf = std::make_unique<amrex::MultiFab>(out_ba, slice->DistributionMap(),
                                                                ncomp, 0);

                for (amrex::MFIter mfi(*out_mf); mfi.isValid(); ++mfi) {
                    const amrex::Box& obx = mfi.validbox();
                    auto const& ofab = out_mf->array(mfi);
                    auto const& sfab = slice->array(mfi);
                    amrex::ParallelFor(obx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
                    {
                        int si = i, sj = j, sk = k;
                        if      (l_dir == 0) { si = k_slice; }
                        else if (l_dir == 1) { sj = k_slice; }
                        else                 { sk = k_slice; }
                        // A sample interpolated from a covered cell is written as
                        // erf.sampling_covered_value, except for the volume fraction
                        ofab(i,j,k,n) = (sfab(si,sj,sk,ncomp) > amrex::Real(0.0) && n != vf_comp)
                                      ? covered_value : sfab(si,sj,sk,n);
                    });
                }

                m_ps_mf[iplane].push_back(std::move(out_mf));

            } // ilev
        }// iplane
    }

    void
    write_sample_data (amrex::Vector<double>& time,
                       amrex::Vector<int>& level_steps,
                       amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Geometry>& geom)
    {
        amrex::ignore_unused(ref_ratio);
        int nplane = m_ps_mf.size();
        for (int iplane(0); iplane<nplane; ++iplane) {
            int nlev_c = static_cast<int>(m_ps_mf[iplane].size());
            if (nlev_c == 0) { continue; }

            int dir = m_dir[iplane];
            amrex::RealBox bnd_rbx = m_bnd_rbx[iplane];
            amrex::Real point = bnd_rbx.lo(dir);

            // One shared physical RealBox for every level: full in-plane extent of the
            // requested plane, and a dir-slab of level-0 thickness centered on the plane.
            // Level 0 supplies prob_lo/hi to the writer, so this object must be shared.
            amrex::Real dx0 = geom[0].CellSize(dir);
            amrex::RealBox shared_rbx = bnd_rbx;
            shared_rbx.setLo(dir, point - myhalf*dx0);
            shared_rbx.setHi(dir, point + myhalf*dx0);

            amrex::Vector<int> is_per(AMREX_SPACEDIM,0);
            for (int d(0); d<AMREX_SPACEDIM; ++d) { is_per[d] = geom[0].isPeriodic(d); }

            // Build the index chain by refining the level-0 plane box so that
            // Domain_l == refine(Domain_{l-1},2) holds exactly in every axis.
            amrex::Box B0 = getIndexBox(bnd_rbx, geom[0]);

            amrex::Vector<const amrex::MultiFab*> mf(nlev_c);
            amrex::Vector<amrex::Geometry>        m_geom(nlev_c);
            amrex::Vector<int>                    m_level_steps(nlev_c);
            amrex::Vector<amrex::IntVect>         m_ref_ratio(nlev_c-1, amrex::IntVect(2));

            for (int l(0); l<nlev_c; ++l) {
                const int rr_l = 1 << l;

                // The nested plane plotfile stores a single scalar refinement
                // ratio per level, so multi-level output requires factor-2
                // isotropic refinement.
                for (int d(0); d<AMREX_SPACEDIM; ++d) {
                    amrex::Real ratio = geom[0].CellSize(d) / geom[l].CellSize(d);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                        std::abs(ratio - static_cast<amrex::Real>(rr_l)) < amrex::Real(1.e-6) * rr_l,
                        "PlaneSampler multi-level output requires factor-2 isotropic refinement (amr.ref_ratio = 2 2 2)");
                }

                amrex::Box domain_l = amrex::refine(B0, rr_l);
                domain_l.setSmall(dir, 0);
                domain_l.setBig(dir, rr_l - 1);

                mf[l] = m_ps_mf[iplane][l].get();
                m_geom[l].define(domain_l, &shared_rbx, geom[l].Coord(), is_per.data());
                m_level_steps[l] = level_steps[l];
                AMREX_ASSERT(domain_l.contains(mf[l]->boxArray().minimalBox()));
            }

            // Create plotfile name
            std::string name_plane = m_name[iplane];
            name_plane += "_step_";
            std::string plotfilename = amrex::Concatenate(name_plane, m_level_steps[0], 5);

            // Write the nested multi-level plane plotfile
            WriteMultiLevelPlotfile(plotfilename, nlev_c, mf,
                                    m_varnames, m_geom, static_cast<amrex::Real>(time[0]),
                                    m_level_steps, m_ref_ratio);
        } // iplane
    }

    int m_max_level = -1;
    amrex::Vector<int> m_dir;
    amrex::Vector<amrex::RealBox> m_bnd_rbx;
    amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>> m_ps_mf;
    amrex::Vector<std::string> m_name;

    amrex::Vector<std::string> m_varnames {"magvel","theta"};
    amrex::Real m_covered_value{sampling_covered_value()}; // erf.sampling_covered_value
};
#endif
