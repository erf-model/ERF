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

#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include <ERF_Container.H>
#include <ERF_DataStruct.H>

#include <iostream>

struct LineSampler
{
    LineSampler ()
    {
        amrex::ParmParse pp("erf");

        bool has_line_idx  = pp.contains("sample_line_lo") || pp.contains("sample_line_hi");
        bool has_line_real = pp.contains("sample_line_lo_real") || pp.contains("sample_line_hi_real");
        if (has_line_idx && has_line_real) {
            amrex::Abort("Specify only one of erf.sample_line_lo/hi or erf.sample_line_lo_real/hi_real");
        }

        // Count number of lo and hi points define the line
        int n_line_lo  = pp.countval("sample_line_lo") / AMREX_SPACEDIM;
        int n_line_hi  = pp.countval("sample_line_hi") / AMREX_SPACEDIM;
        int n_line_lo_real = pp.countval("sample_line_lo_real") / AMREX_SPACEDIM;
        int n_line_hi_real = pp.countval("sample_line_hi_real") / AMREX_SPACEDIM;
        int n_line_dir = pp.countval("sample_line_dir");
        if (has_line_idx) {
            AMREX_ALWAYS_ASSERT( (n_line_lo==n_line_hi ) &&
                                 (n_line_lo==n_line_dir) );
        } else if (has_line_real) {
            AMREX_ALWAYS_ASSERT( (n_line_lo_real==n_line_hi_real ) &&
                                 (n_line_lo_real==n_line_dir) );
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

            // The cells of a line placed by height are found once, so they do not follow a
            // mesh that moves; lines given by indices do not use heights
            if (m_use_real_bx && use_column_heights() &&
                SolverChoice::terrain_type == TerrainType::MovingFittedMesh) {
                amrex::Abort("Line sampling with erf.sample_line_lo_real/hi_real is not supported "
                             "with erf.terrain_type = MovingFittedMesh");
            }

            // Where the levels are not equally spaced, the z of the plotfile of a line is not
            // the height of its samples, so the height is written with them
            if (SolverChoice::mesh_type != MeshType::ConstantDz &&
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

        m_write_log = true;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            m_logstream.resize(nline);
            for (int iline(0); iline<nline; ++iline) {
                m_logstream[iline] = std::make_unique<std::fstream>();
                m_logstream[iline]->open(log_names[iline].c_str(),std::ios::out|std::ios::app);
                if (!m_logstream[iline]->good()) {
                    amrex::FileOpenFailed(log_names[iline]);
                }
            }
        }
        amrex::ParallelDescriptor::Barrier("LineSampler::LineSampler");
    }

    //! The variables that may be sampled, in the order they are stored and written
    static const amrex::Vector<std::string>&
    line_sampling_var_order ()
    {
        static const amrex::Vector<std::string> order = {
            "conserved", "z", "density", "x_velocity", "y_velocity", "z_velocity", "magvel", "theta",
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

    //! Whether the sampler needs the ground under its lines (erf.sample_line_z_agl)
    [[nodiscard]] bool
    needs_ground () const
    {
        return m_use_real_bx && m_z_agl;
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
                  const amrex::MultiFab* ground_nd = nullptr)
    {
        if (!m_write_ascii) return;

        amrex::Print() << "Writing out line coordinates to text" << std::endl;

        // Only the levels that are sampled
        int nlev = static_cast<int>(z_phys_cc.size());
        if (m_max_level >= 0) { nlev = std::min(nlev, m_max_level + 1); }
        for (int lev(0); lev < nlev; ++lev) {
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
                     const amrex::MultiFab* ground_nd = nullptr)
    {
        if (m_ncons < 0) { expand_conserved(vars_new[0][Vars::cons].nComp()); }

        int nlev  = static_cast<int>(vars_new.size());
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
        auto line_box = [&] (int iline, int ilev) -> amrex::Box {
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

            for (int n(0); n<static_cast<int>(lines.size()); ++n) {
                const int iline = lines[n];
                const amrex::Box& bnd_bx = boxes[n];
                m_lev[iline] = ilev;
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
    bool m_write_ascii{false};
    bool m_write_log{false};
    bool m_z_agl{false};   // erf.sample_line_z_agl
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
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
    {
        int nlev   = static_cast<int>(vars_new.size());
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

                // Construct MultiFab holding requested variables
                amrex::MultiFab mf_cc_data;
                mf_cc_data.define(ba, dm, ncomp, 1);

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

                auto slice = get_slice_data(dir, point, mf_cc_data, geom[ilev],
                                            0, ncomp, interpolate, bnd_rbx);

                // Defensive: an empty return means nothing was actually sampled.
                if (!slice || slice->boxArray().size() == 0) { break; }

                // Broadcast the single sampled plane across 2^ilev cells in dir so
                // every level shares an isotropic refinement ratio. The dir extent
                // becomes [0, 2^ilev-1]; each cell holds the interpolated plane value
                // replicated across the shared level-0 slab thickness.
                const int rr_l    = 1 << ilev;
                const int k_slice = slice->boxArray().minimalBox().smallEnd(dir);
                const int l_dir   = dir;

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
                        ofab(i,j,k,n) = sfab(si,sj,sk,n);
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
};
#endif
