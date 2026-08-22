/**
 * \file ERF_SampleData.H
 */
#ifndef ERF_SAMPLEDATA_H
#define ERF_SAMPLEDATA_H

#include <memory>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include <ERF_Container.H>

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

            // Get requested vars
            if (pp.countval("line_sampling_vars") > 0) {
                m_varnames.clear();
                amrex::Vector<std::string> requested_vars;
                pp.queryarr("line_sampling_vars",requested_vars);
                amrex::Print() << "Selected line sampling vars :";
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

            // Write outputs to text files, one file per variable, with all
            // times appended to the same file
            pp.query("line_sampling_text_output",m_write_ascii);
            if (m_write_ascii && amrex::ParallelDescriptor::IOProcessor()) {
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
        }
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
                  amrex::Vector<amrex::Geometry>& geom)
    {
        if (!m_write_ascii) return;

        amrex::Print() << "Writing out line coordinates to text" << std::endl;

        for (int lev(0); lev < z_phys_cc.size(); ++lev) {
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
                amrex::Box bnd_bx = (m_use_real_bx)
                                  ? getIndexBox(m_bnd_rbx[iline], geom[lev])
                                  : m_bnd_bx[iline];
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
                     amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new)
    {
        int nlev  = static_cast<int>(vars_new.size());
        int nline = static_cast<int>(m_bnd_bx.size());
        int ncomp = static_cast<int>(m_varnames.size());

        int qv_comp = -1;

        // Loop over each line
        for (int iline(0); iline<nline; ++iline) {
            int dir = m_dir[iline];

            // Search each level to get the finest data possible
            for (int ilev(nlev-1); ilev>=0; --ilev) {
                amrex::Box bnd_bx = (m_use_real_bx)
                                  ? getIndexBox(m_bnd_rbx[iline], geom[ilev])
                                  : m_bnd_bx[iline];
                amrex::IntVect cell = bnd_bx.smallEnd();

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

                  m_lev[iline] = ilev;
                m_ls_mf[iline] = get_line_data(mf_cc_data, dir, cell, bnd_bx);

                // We can stop if we got the entire line
                auto min_bnd_bx = m_ls_mf[iline].boxArray().minimalBox();
                if (bnd_bx == min_bnd_bx) { break; }

            } // ilev

            // Fill bnd_bx if use_rbx
            if (m_use_real_bx) {
                m_bnd_bx[iline] = m_ls_mf[iline].boxArray().minimalBox();
            }

        }// iline
    }

    void
    write_sample_data (amrex::Vector<double>& time,
                       amrex::Vector<int>& level_steps,
                       amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Geometry>& geom)
    {
        if (m_write_ascii) {
            write_line_ascii(time);
        } else {
            write_line_plotfile(time, level_steps, ref_ratio, geom);
        }
    }

    void
    write_line_ascii (amrex::Vector<double>& time)
    {
        // same as definitions in ERF.H
        constexpr int datwidth = 14;
        constexpr int datprecision = 6;
        constexpr int timeprecision = 13;

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
    amrex::Vector<std::string> m_varnames {"magvel","theta"};
    amrex::Vector<std::unique_ptr<std::fstream> > m_datastream;
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
            pp.query("plane_sampling_max_level", m_max_level);

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
