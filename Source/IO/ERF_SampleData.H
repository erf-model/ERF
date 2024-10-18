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
            for (int i = 0; i < n_line_lo; i++) {
                amrex::IntVect iv(idx_lo[AMREX_SPACEDIM*i+0],
                                  idx_lo[AMREX_SPACEDIM*i+1],
                                  idx_lo[AMREX_SPACEDIM*i+2]);
                iv_lo[i] = iv;
            }

            // Parse hi
            amrex::Vector<int> idx_hi; idx_hi.resize(n_line_hi*AMREX_SPACEDIM);
            amrex::Vector<amrex::IntVect> iv_hi; iv_hi.resize(n_line_hi);
            pp.queryarr("sample_line_hi",idx_hi,0,n_line_hi*AMREX_SPACEDIM);
            for (int i = 0; i < n_line_hi; i++) {
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
                if (bnd_bx == min_bnd_bx) continue;

            } // ilev
        }// iline
    }

    void
    write_line_mfs (amrex::Vector<amrex::Geometry>& geom)
    {
        int lev = m_lev[0];

        amrex::Real time = 0.0;
        amrex::Vector<int> level_steps  = {0};
        amrex::Vector<amrex::IntVect> ref_ratio = {amrex::IntVect(1)};
        amrex::Vector<std::string> varnames = {"r", "rT", "rKE", "rQKE", "rS"};
        std::string name_base = "plt_line_";

        std::string plotfilename = amrex::Concatenate(name_base, 0, 5);

        amrex::Vector<const amrex::MultiFab*> mfp = {&(m_ls_mf[0])};
        amrex::Vector<amrex::Geometry> temp = {geom[lev]};

        // HACK GEOM TO LINE
        auto plo = geom[lev].ProbLo();
        auto dx  = geom[lev].CellSize();

        amrex::Vector<amrex::Geometry> m_geom; m_geom.resize(1);
        amrex::Vector<int> is_per(3,0);
        amrex::Box m_dom = m_bnd_bx[0];
        amrex::RealBox m_rb;
        for (int i(0); i<3; ++i) {
            int offset = (i==2) ? 0 : 1;
            amrex::Real xlo = plo[i] + ( m_dom.smallEnd(i) ) * dx[i];
            amrex::Real xhi = plo[i] + ( m_dom.bigEnd(i) + offset ) * dx[i];

            m_rb.setLo(i,xlo);
            m_rb.setHi(i,xhi);

            is_per[i] = geom[lev].isPeriodic(i);
        }


        m_geom[0].define(m_dom, &m_rb, geom[lev].Coord(), is_per.data());
        WriteMultiLevelPlotfile(plotfilename, 1, mfp,
                                varnames, m_geom, time, level_steps, ref_ratio);


        /*
        WriteMultiLevelPlotfile(plotfilename, 1, mfp,
                                varnames, temp, time, level_steps, ref_ratio);
        */
    }

    amrex::Vector<int> m_dir;
    amrex::Vector<int> m_lev;
    amrex::Vector<amrex::Box> m_bnd_bx;
    amrex::Vector<amrex::MultiFab> m_ls_mf;
};

class SampleData
{
public:
    explicit SampleData ()
    {
        amrex::ParmParse pp("erf");
        if(pp.contains("sample_line")) m_ls = std::make_unique<LineSampler>();
    }

    void
    get_sample_data (amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
    {
        if (m_ls) m_ls->get_line_mfs(vars_new);
    }

    void
    write_sample_data (amrex::Vector<amrex::Geometry>& geom)
    {
        if (m_ls) m_ls->write_line_mfs(geom);
    }

private:

    // Geometry objects for all levels
    amrex::Vector<amrex::Geometry> m_geom;

    // Variables for IO
    amrex::Vector<amrex::Vector<std::string>> m_var_names;

    // Structures for getting line MFs
    std::unique_ptr<LineSampler> m_ls = nullptr;

    // Structures for getting plane MFs
    //std::unique_ptr<PlaneSampler> m_ps = nullptr;
};
#endif
