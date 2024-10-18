#ifndef ERF_SAMPLEDATA_H
#define ERF_SAMPLEDATA_H

#include <memory>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <ERF_IndexDefines.H>

struct LineSampler
{
    LineSampler ()
    {
        amrex::ParmParse pp("erf");

        // Count number of lo and hi points define the line
        int n_line_lo  = pp.countval("sample_line_lo") / AMREX_SPACEDIM;
        int n_line_hi  = pp.countval("sample_line_hi") / AMREX_SPACEDIM;
        int n_line_dir = pp.countval("sample_line_dir");
        AMREX_ALWAYS_ASSERT( (n_line_lo==n_line_hi ) &&
                             (n_line_lo==n_line_dir) );

        // Parse the data
        if (n_line_lo > 0) {
            // Parse lo
            amrex::Vector<int> idx_lo; idx_lo.resize(n_line_lo*AMREX_SPACEDIM);
            amrex::Vector<amrex::IntVect> iv_lo; iv_lo.resize(n_line_lo);
            pp.queryarr("sample_line_lo",idx_lo,0,n_line_lo*AMREX_SPACEDIM);
            for (int i(0); i < n_line_lo; i++) {
                amrex::IntVect iv(idx_lo[AMREX_SPACEDIM*i+0],
                                  idx_lo[AMREX_SPACEDIM*i+1],
                                  idx_lo[AMREX_SPACEDIM*i+2]);
                iv_lo[i] = iv;
            }

            // Parse hi
            amrex::Vector<int> idx_hi; idx_hi.resize(n_line_hi*AMREX_SPACEDIM);
            amrex::Vector<amrex::IntVect> iv_hi; iv_hi.resize(n_line_hi);
            pp.queryarr("sample_line_hi",idx_hi,0,n_line_hi*AMREX_SPACEDIM);
            for (int i(0); i < n_line_hi; i++) {
                amrex::IntVect iv(idx_hi[AMREX_SPACEDIM*i+0],
                                  idx_hi[AMREX_SPACEDIM*i+1],
                                  idx_hi[AMREX_SPACEDIM*i+2]);
                iv_hi[i] = iv;
            }

            // Construct vector of bounding boxes
            m_bnd_bx.resize(n_line_lo);
            for (int i = 0; i < n_line_hi; i++){
                amrex::Box lbx(iv_lo[i],iv_hi[i]);
                m_bnd_bx[i] = lbx;
            }

            // Parse directionality
            m_dir.resize(n_line_dir);
            pp.queryarr("sample_line_dir",m_dir,0,n_line_dir);

            // Allocate space for level indicator
            m_lev.resize(n_line_dir,0);

            // Allocate space for MF pointers
            m_ls_mf.resize(n_line_lo);
        }
    }

    void
    get_line_mfs (amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
    {
        int nlev  = vars_new.size();
        int nline = m_bnd_bx.size();

        // Loop over each line
        for (int iline(0); iline<nline; ++iline) {
            int dir = m_dir[iline];
            amrex::Box bnd_bx   = m_bnd_bx[iline];
            amrex::IntVect cell = bnd_bx.smallEnd();

            // Search each level to get the finest data possible
            for (int ilev(nlev-1); ilev>=0; --ilev) {
                amrex::MultiFab& mf = vars_new[ilev][Vars::cons];
                  m_lev[iline] = ilev;
                m_ls_mf[iline] = get_line_data(mf, dir, cell, bnd_bx);

                // We can stop if we got the entire line
                auto min_bnd_bx = m_ls_mf[iline].boxArray().minimalBox();
                if (bnd_bx == min_bnd_bx) { continue; }

            } // ilev
        }// iline
    }

    void
    write_line_mfs (amrex::Vector<amrex::Real>& time,
                    amrex::Vector<int>& level_steps,
                    amrex::Vector<amrex::IntVect>& ref_ratio,
                    amrex::Vector<amrex::Geometry>& geom)
    {
        std::string name_base = "plt_line_";
        amrex::Vector<std::string> varnames = {"r", "rT", "rKE", "rQKE", "rS"};

        int nline = m_ls_mf.size();
        for (int iline(0); iline<nline; ++iline) {
            // Data members that can be used as-is
            int dir = m_dir[iline];
            int lev = m_lev[iline];
            amrex::Real m_time = time[lev];
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
                int offset = (d==dir) ? 0 : 1;
                amrex::Real lo = plo[d] + ( m_dom.smallEnd(d)        ) * dx[d];
                amrex::Real hi = plo[d] + ( m_dom.bigEnd(d) + offset ) * dx[d];

                m_rb.setLo(d,lo);
                m_rb.setHi(d,hi);

                is_per[d] = geom[lev].isPeriodic(d);
            }
            m_geom[0].define(m_dom, &m_rb, geom[lev].Coord(), is_per.data());

            // Create plotfile name
            std::string name_line    = amrex::Concatenate(name_base, iline           , 5);
            name_line += "_step_";
            std::string plotfilename = amrex::Concatenate(name_line, m_level_steps[0], 5);

            // Get the data
            amrex::Vector<const amrex::MultiFab*> mf = {&(m_ls_mf[iline])};

            // Write each line
            WriteMultiLevelPlotfile(plotfilename, 1, mf,
                                    varnames, m_geom, m_time,
                                    m_level_steps, m_ref_ratio);
        }
    }

    amrex::Vector<int> m_dir;
    amrex::Vector<int> m_lev;
    amrex::Vector<amrex::Box> m_bnd_bx;
    amrex::Vector<amrex::MultiFab> m_ls_mf;
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

            // Allocate space for level indicator
            m_lev.resize(n_plane_dir,0);

            // Allocate space for MF pointers
            m_ps_mf.resize(n_plane_lo);
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
    get_plane_mfs (amrex::Vector<amrex::Geometry>& geom,
                   amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
    {
        int nlev   = vars_new.size();
        int nplane = m_bnd_rbx.size();
        int ncomp  = vars_new[0][Vars::cons].nComp();
        bool interpolate = true;

        // Loop over each plane
        for (int iplane(0); iplane<nplane; ++iplane) {
            int dir = m_dir[iplane];
            amrex::RealBox bnd_rbx  = m_bnd_rbx[iplane];
            amrex::Real point = bnd_rbx.lo(dir);

            // Search each level to get the finest data possible
            for (int ilev(nlev-1); ilev>=0; --ilev) {
                amrex::MultiFab& mf = vars_new[ilev][Vars::cons];
                  m_lev[iplane] = ilev;
                m_ps_mf[iplane] = get_slice_data(dir, point, mf, geom[ilev],
                                                 0, ncomp, interpolate, bnd_rbx);

                // We can stop if we got the entire plane
                auto min_bnd_bx   = m_ps_mf[iplane]->boxArray().minimalBox();
                amrex::Box bnd_bx = getIndexBox(bnd_rbx, geom[ilev]);
                if (bnd_bx == min_bnd_bx) { continue; }

            } // ilev
        }// iplane
    }

    void
    write_plane_mfs (amrex::Vector<amrex::Real>& time,
                     amrex::Vector<int>& level_steps,
                     amrex::Vector<amrex::IntVect>& ref_ratio,
                     amrex::Vector<amrex::Geometry>& geom)
    {
        std::string name_base = "plt_plane_";
        amrex::Vector<std::string> varnames = {"r", "rT", "rKE", "rQKE", "rS"};

        int nplane = m_ps_mf.size();
        for (int iplane(0); iplane<nplane; ++iplane) {
            // Data members that can be used as-is
            int dir = m_dir[iplane];
            int lev = m_lev[iplane];
            amrex::Real m_time = time[lev];
            amrex::Vector<int> m_level_steps = {level_steps[lev]};
            amrex::Vector<amrex::IntVect> m_ref_ratio  = {ref_ratio[lev]};

            // Create modified geometry object corresponding to the plane
            amrex::RealBox m_rb = m_bnd_rbx[iplane];
            amrex::Box m_dom    = getIndexBox(m_rb, geom[lev]);
            amrex::Real point   = m_rb.hi(dir);
            amrex::Vector<int> is_per(AMREX_SPACEDIM,0);
            for (int d(0); d<AMREX_SPACEDIM; ++d) {
                if (d==dir) {
                    m_dom.growHi(d,1);
                    m_rb.setHi(d,point+geom[lev].CellSize(d));
                }
                is_per[d] = geom[lev].isPeriodic(d);
            }
            amrex::Vector<amrex::Geometry> m_geom; m_geom.resize(1);
            m_geom[0].define(m_dom, &m_rb, geom[lev].Coord(), is_per.data());

            // Create plotfile name
            std::string name_plane   = amrex::Concatenate(name_base, iplane           , 5);
            name_plane += "_step_";
            std::string plotfilename = amrex::Concatenate(name_plane, m_level_steps[0], 5);

            // Get the data
            amrex::Vector<const amrex::MultiFab*> mf = {m_ps_mf[iplane].get()};

            // Write each plane
            WriteMultiLevelPlotfile(plotfilename, 1, mf,
                                    varnames, m_geom, m_time,
                                    m_level_steps, m_ref_ratio);
        } // iplane
    }

    amrex::Vector<int> m_dir;
    amrex::Vector<int> m_lev;
    amrex::Vector<amrex::RealBox> m_bnd_rbx;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_ps_mf;
};


class SampleData
{
public:
    explicit SampleData (bool do_line=false,
                         bool do_plane=false)
    {
        if(do_line)  m_ls = std::make_unique<LineSampler >();
        if(do_plane) m_ps = std::make_unique<PlaneSampler>();
    }

    void
    get_sample_data (amrex::Vector<amrex::Geometry>& geom,
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
    {
        if (m_ls) m_ls->get_line_mfs(vars_new);
        if (m_ps) m_ps->get_plane_mfs(geom,vars_new);
    }

    void
    write_sample_data (amrex::Vector<amrex::Real>& time,
                       amrex::Vector<int>& level_steps,
                       amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Geometry>& geom)
    {
        if (m_ls) m_ls->write_line_mfs(time, level_steps, ref_ratio, geom);
        if (m_ls) m_ps->write_plane_mfs(time, level_steps, ref_ratio, geom);
    }

private:

    // Structures for getting line MFs
    std::unique_ptr<LineSampler> m_ls = nullptr;

    // Structures for getting plane MFs
    std::unique_ptr<PlaneSampler> m_ps = nullptr;
};
#endif
