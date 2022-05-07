#include <EOS.H>
#include <ERF.H>
#include "AMReX_Interp_3D_C.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

// get plotfile name
std::string
ERF::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

void
ERF::setPlotVariables ()
{
    ParmParse pp("amr");

    if (pp.contains("plot_vars"))
    {
        std::string nm;

        int nPltVars = pp.countval("plot_vars");

        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("plot_vars", nm, i);

            if (nm == "ALL") {
                // put all conserved state variables in the plot variables list
                plot_state_names = cons_names;
                // put all velocity components in the plot variables list
                plot_state_names.insert(plot_state_names.end(), velocity_names.begin(), velocity_names.end());
            } else if (nm == "NONE" || nm == "None") {
                // put no state variables in the plot variable list
                plot_state_names.clear();
            } else {
                // add the named variable to our list of plot variables
                // if it is not already in the list
                if (!containerHasElement(plot_state_names, nm)) {
                    plot_state_names.push_back(nm);
                }
            }
        }
    }
    else
    {
        //
        // The default is to add all state conserved and velocity variables to the plot variables list
        //
        plot_state_names = cons_names;
        plot_state_names.insert(plot_state_names.end(), velocity_names.begin(), velocity_names.end());
    }

    // Get state variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_plot_names;
    for (int i = 0; i < Cons::NumVars; ++i) {
        if (containerHasElement(plot_state_names, cons_names[i])) {
            tmp_plot_names.push_back(cons_names[i]);
        }
    }
    // check for velocity since it's not in cons_names
    // if we are asked for any velocity component, we will need them all
    if (containerHasElement(plot_state_names, "x_velocity") ||
        containerHasElement(plot_state_names, "y_velocity") ||
        containerHasElement(plot_state_names, "z_velocity")) {
        tmp_plot_names.push_back("x_velocity");
        tmp_plot_names.push_back("y_velocity");
        tmp_plot_names.push_back("z_velocity");
    }
    // check to see if we found all the requested vairables
    for (auto plot_name : plot_state_names) {
      if (!containerHasElement(tmp_plot_names, plot_name)) {
    amrex::Warning("\nWARNING: Requested to plot variable '" + plot_name + "' but it is not available");
      }
    }
    plot_state_names = tmp_plot_names;

    if (pp.contains("derive_plot_vars"))
    {
        std::string nm;

        int nDrvPltVars = pp.countval("derive_plot_vars");

        for (int i = 0; i < nDrvPltVars; i++)
        {
            pp.get("derive_plot_vars", nm, i);

            if (nm == "ALL") {
                // put all diagnostic variables in the plot variables list
                plot_deriv_names = derived_names;
            } else if (nm == "NONE" || nm == "None") {
                // put no diagnostic variables in the plot variable list
                plot_deriv_names.clear();
            } else {
                // add the named variable to our list of plot variables
                // if it is not already in the list
                if (!containerHasElement(plot_deriv_names, nm)) {
                    plot_deriv_names.push_back(nm);
                }
            }
        }
    }
    else
    {
        //
        // The default is to add none of the diagnostic variables to the plot variables list
        //
        plot_deriv_names.clear();
    }

    // Get derived variables in the same order as we define them,
    // since they may be in any order in the input list
    Vector<std::string> tmp_deriv_names;
    for (int i = 0; i < derived_names.size(); ++i) {
        if (containerHasElement(plot_deriv_names, derived_names[i])) {
            tmp_deriv_names.push_back(derived_names[i]);
        }
    }
    // check to see if we found all the requested vairables
    for (auto plot_name : plot_deriv_names) {
        if (!containerHasElement(tmp_deriv_names, plot_name)) {
            amrex::Warning("\nWARNING: Requested to plot derived variable '" + plot_name + "' but it is not available");
        }
    }

    plot_deriv_names = tmp_deriv_names;
}

// set plotfile variable names
Vector<std::string>
ERF::PlotFileVarNames () const
{
    Vector<std::string> names;

    names.insert(names.end(), plot_state_names.begin(), plot_state_names.end());
    names.insert(names.end(), plot_deriv_names.begin(), plot_deriv_names.end());

    return names;

}

// write plotfile to disk
void
ERF::WritePlotFile () const
{
    const Vector<std::string> varnames = PlotFileVarNames();
    const int ncomp_mf = varnames.size();

    if (ncomp_mf == 0)
        return;

    Vector<MultiFab> mf(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        mf[lev].define(grids[lev], dmap[lev], ncomp_mf, 0);
    }

#ifdef ERF_USE_TERRAIN
    Vector<MultiFab> mf_nd(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray nodal_grids(grids[lev]); nodal_grids.surroundingNodes();
        mf_nd[lev].define(nodal_grids, dmap[lev], ncomp_mf, 0);
        mf_nd[lev].setVal(0.);
    }
#endif

    for (int lev = 0; lev <= finest_level; ++lev) {
        int mf_comp = 0;

        // First, copy any of the conserved state variables into the output plotfile
        AMREX_ALWAYS_ASSERT(cons_names.size() == Cons::NumVars);
        for (int i = 0; i < Cons::NumVars; ++i) {
            if (containerHasElement(plot_state_names, cons_names[i])) {
                MultiFab::Copy(mf[lev],vars_new[lev][Vars::cons],i,mf_comp,1,0);
                mf_comp++;
            }
        }

        // Next, check for velocities and if desired, output them
        if (containerHasElement(plot_state_names, "x_velocity") ||
            containerHasElement(plot_state_names, "y_velocity") ||
            containerHasElement(plot_state_names, "z_velocity")) {

            average_face_to_cellcenter(mf[lev],mf_comp,
                Array<const MultiFab*,3>{&vars_new[lev][Vars::xvel],&vars_new[lev][Vars::yvel],&vars_new[lev][Vars::zvel]});
            mf_comp += AMREX_SPACEDIM;
        }

        // Finally, check for any derived quantities and compute them, inserting
        // them into our output multifab
        auto calculate_derived = [&](const std::string der_name,
                                     const decltype(derived::erf_dernull)& der_function)
        {
            if (containerHasElement(plot_deriv_names, der_name)) {
                MultiFab dmf(mf[lev], amrex::make_alias, mf_comp, 1);

                for (MFIter mfi(dmf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    auto& dfab = dmf[mfi];
                    auto& sfab = vars_new[lev][Vars::cons][mfi];
                    der_function(bx, dfab, 0, 1, sfab, Geom(lev), t_new[0], nullptr, lev);
                }

                mf_comp++;
            }
        };

        // Note: All derived variables must be computed in order of "derived_names" defined in ERF.H
        calculate_derived("pressure",    derived::erf_derpres);
        calculate_derived("soundspeed",  derived::erf_dersoundspeed);
        calculate_derived("temp",        derived::erf_dertemp);
        calculate_derived("theta",       derived::erf_dertheta);
        calculate_derived("KE",          derived::erf_derKE);
        calculate_derived("QKE",         derived::erf_derQKE);
        calculate_derived("scalar",      derived::erf_derscalar);

        if (containerHasElement(plot_deriv_names, "pres_hse"))
        {
#ifdef ERF_USE_TERRAIN
            MultiFab::Copy(mf[lev],pres_hse[lev],0,mf_comp,1,0);
#else
            auto d_pres_hse_lev = d_pres_hse[lev].dataPtr();
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                amrex::ParallelFor(bx, [=, ng_pres_hse=ng_pres_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i, j, k, mf_comp) = d_pres_hse_lev[k+ng_pres_hse];
                });
            }
#endif
            mf_comp += 1;
        }

        if (containerHasElement(plot_deriv_names, "dens_hse"))
        {
#ifdef ERF_USE_TERRAIN
            MultiFab::Copy(mf[lev],dens_hse[lev],0,mf_comp,1,0);
#else
            auto d_dens_hse_lev = d_dens_hse[lev].dataPtr();
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
                amrex::ParallelFor(bx, [=, ng_dens_hse=ng_dens_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    derdat(i, j, k, mf_comp) = d_dens_hse_lev[k+ng_dens_hse];
                });
            }
#endif
            mf_comp ++;
        }

        if (containerHasElement(plot_deriv_names, "pert_pres"))
        {
#ifndef ERF_USE_TERRAIN
            auto d_pres_hse_lev = d_pres_hse[lev].dataPtr();
#endif
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat = mf[lev].array(mfi);
#ifdef ERF_USE_TERRAIN
                const Array4<Real const>& p0_arr = pres_hse[lev].const_array(mfi);
#endif
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
                amrex::ParallelFor(bx, [=, ng_pres_hse=ng_pres_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    const Real rhotheta = S_arr(i,j,k,RhoTheta_comp);
#ifdef ERF_USE_TERRAIN
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rhotheta) - p0_arr(i,j,k);
#else
                    derdat(i, j, k, mf_comp) = getPgivenRTh(rhotheta) - d_pres_hse_lev[k+ng_pres_hse];
#endif
                });
            }
            mf_comp ++;
        }

        if (containerHasElement(plot_deriv_names, "pert_dens"))
        {
#ifndef ERF_USE_TERRAIN
            auto d_dens_hse_lev = d_dens_hse[lev].dataPtr();
#endif
            for ( amrex::MFIter mfi(mf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                const Array4<Real>& derdat  = mf[lev].array(mfi);
                const Array4<Real const>& S_arr = vars_new[lev][Vars::cons].const_array(mfi);
#ifdef ERF_USE_TERRAIN
                const Array4<Real const>& r0_arr = dens_hse[lev].const_array(mfi);
#endif
                amrex::ParallelFor(bx, [=, ng_dens_hse=ng_dens_hse] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef ERF_USE_TERRAIN
                    derdat(i, j, k, mf_comp) = S_arr(i,j,k,Rho_comp) - r0_arr(i,j,k);
#else
                    derdat(i, j, k, mf_comp) = S_arr(i,j,k,Rho_comp) - d_dens_hse_lev[k+ng_dens_hse];
#endif
                });
            }
            mf_comp ++;
        }

#ifdef ERF_USE_TERRAIN
        if (containerHasElement(plot_deriv_names, "z_phys"))
        {
            MultiFab::Copy(mf[lev],z_phys_cc[lev],0,mf_comp,1,0);
            mf_comp ++;
        }
    if (containerHasElement(plot_deriv_names, "detJ"))
        {
            MultiFab::Copy(mf[lev],detJ_cc[lev],0,mf_comp,1,0);
            mf_comp ++;
        }
#endif
    }

    const std::string& plotfilename = PlotFileName(istep[0]);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    if (finest_level == 0)
    {
        if (plotfile_type == "amrex") {
#ifdef ERF_USE_TERRAIN
            // We started with mf_nd holding 0 in every component; here we fill only the offset in z
            int lev = 0;
            MultiFab::Copy(mf_nd[lev],z_phys_nd[lev],0,2,1,0);
            Real dz = Geom()[lev].CellSizeArray()[2];
            for (MFIter mfi(mf_nd[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                Array4<      Real> mf_arr = mf_nd[lev].array(mfi);
                Array4<const Real>  z_arr = z_phys_nd[lev].const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    mf_arr(i,j,k,2) -= k * dz;
                });
            }
            WriteMultiLevelPlotfileWithTerrain(plotfilename, finest_level+1,
                                               GetVecOfConstPtrs(mf),
                                               GetVecOfConstPtrs(mf_nd),
                                               varnames,
                                               Geom(), t_new[0], istep, refRatio());
#else
            amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                           GetVecOfConstPtrs(mf),
                                           varnames,
                                           Geom(), t_new[0], istep, refRatio());
#endif
            writeJobInfo(plotfilename);
#ifdef ERF_USE_NETCDF
        } else {
             writeNCPlotFile(plotfilename, GetVecOfConstPtrs(mf), varnames, istep, t_new[0]);
#endif
        }

    } else {

        amrex::Vector<IntVect>   r2(finest_level);
        amrex::Vector<Geometry>  g2(finest_level+1);
        amrex::Vector<MultiFab> mf2(finest_level+1);

        mf2[0].define(grids[0], dmap[0], ncomp_mf, 0);

        // Copy level 0 as is
        MultiFab::Copy(mf2[0],mf[0],0,0,mf[0].nComp(),0);

        // Define a new multi-level array of Geometry's so that we pass the new "domain" at lev > 0
        Array<int,AMREX_SPACEDIM> periodicity =
                     {Geom()[0].isPeriodic(0),Geom()[0].isPeriodic(1),Geom()[0].isPeriodic(2)};
        g2[0].define(Geom()[0].Domain(),&(Geom()[0].ProbDomain()),0,periodicity.data());

        r2[0] = IntVect(1,1,ref_ratio[0][0]);
        for (int lev = 1; lev <= finest_level; ++lev) {
            if (lev > 1) {
                r2[lev-1][0] = 1;
                r2[lev-1][1] = 1;
                r2[lev-1][2] = r2[lev-2][2] * ref_ratio[lev-1][0];
            }

            mf2[lev].define(amrex::refine(grids[lev],r2[lev-1]), dmap[lev], ncomp_mf, 0);

            // Set the new problem domain
            Box d2(Geom()[lev].Domain());
            d2.refine(r2[lev-1]);

            g2[lev].define(d2,&(Geom()[lev].ProbDomain()),0,periodicity.data());
        }

        // Do piecewise interpolation of mf into mf2
        for (int lev = 1; lev <= finest_level; ++lev) {
            for (MFIter mfi(mf2[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                pcinterp_interp(bx,mf2[lev].array(mfi), 0, mf[lev].nComp(), mf[lev].const_array(mfi),0,r2[lev-1]);
            }
        }

        // Define an effective ref_ratio which is isotropic to be passed into WriteMultiLevelPlotfile
        amrex::Vector<IntVect> rr(finest_level);
        for (int lev = 0; lev < finest_level; ++lev) {
            rr[lev] = IntVect(ref_ratio[lev][0],ref_ratio[lev][1],ref_ratio[lev][0]);
        }

        if (plotfile_type == "amrex") {
            amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, GetVecOfConstPtrs(mf2), varnames,
                                           g2, t_new[0], istep, rr);
            writeJobInfo(plotfilename);
#ifdef ERF_USE_NETCDF
        } else {
             writeNCPlotFile(plotfilename, GetVecOfConstPtrs(mf2), varnames, istep, t_new[0]);
#endif
        }
    } // end multi-level
}

#ifdef ERF_USE_TERRAIN
void
ERF::WriteMultiLevelPlotfileWithTerrain (const std::string& plotfilename, int nlevels,
                                         const Vector<const MultiFab*>& mf,
                                         const Vector<const MultiFab*>& mf_nd,
                                         const Vector<std::string>& varnames,
                                         const Vector<Geometry>& geom, Real time,
                                         const Vector<int>& level_steps,
                                         const Vector<IntVect>& ref_ratio,
                                         const std::string &versionName,
                                         const std::string &levelPrefix,
                                         const std::string &mfPrefix,
                                         const Vector<std::string>& extra_dirs) const
{
    BL_PROFILE("WriteMultiLevelPlotfileWithTerrain()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::MyProc() == ParallelDescriptor::NProcs()-1) {
        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        auto f = [=]() {
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
            std::string HeaderFileName(plotfilename + "/Header");
            std::ofstream HeaderFile;
            HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) FileOpenFailed(HeaderFileName);
            WriteGenericPlotfileHeaderWithTerrain(HeaderFile, nlevels, boxArrays, varnames,
                                                  geom, time, level_steps, ref_ratio, versionName,
                                                  levelPrefix, mfPrefix);
        };

        if (AsyncOut::UseAsyncOut()) {
            AsyncOut::Submit(std::move(f));
        } else {
            f();
        }
    }

    std::string mf_nodal_prefix = "Nu_nd";
    for (int level = 0; level <= finest_level; ++level)
    {
        if (AsyncOut::UseAsyncOut()) {
            VisMF::AsyncWrite(*mf[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix),
                              true);
            VisMF::AsyncWrite(*mf_nd[level],
                              MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix),
                              true);
        } else {
            const MultiFab* data;
            std::unique_ptr<MultiFab> mf_tmp;
            if (mf[level]->nGrowVect() != 0) {
                mf_tmp = std::make_unique<MultiFab>(mf[level]->boxArray(),
                                                    mf[level]->DistributionMap(),
                                                    mf[level]->nComp(), 0, MFInfo(),
                                                    mf[level]->Factory());
                MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
                data = mf_tmp.get();
            } else {
                data = mf[level];
            }
            VisMF::Write(*data       , MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
            VisMF::Write(*mf_nd[level], MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mf_nodal_prefix));
        }
    }
}

void
ERF::WriteGenericPlotfileHeaderWithTerrain (std::ostream &HeaderFile,
                                            int nlevels,
                                            const Vector<BoxArray> &bArray,
                                            const Vector<std::string> &varnames,
                                            const Vector<Geometry> &geom,
                                            Real time,
                                            const Vector<int> &level_steps,
                                            const Vector<IntVect> &ref_ratio,
                                            const std::string &versionName,
                                            const std::string &levelPrefix,
                                            const std::string &mfPrefix) const
{
        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= geom.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        int finest_level(nlevels - 1);

        HeaderFile.precision(17);

        // ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
            HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbLo(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbHi(i) << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << geom[i].Domain() << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
        }
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
            }
            HeaderFile << '\n';
        }
        HeaderFile << (int) geom[0].Coord() << '\n';
        HeaderFile << "0\n";

        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
            HeaderFile << level_steps[level] << '\n';

            const IntVect& domain_lo = geom[level].Domain().smallEnd();
            for (int i = 0; i < bArray[level].size(); ++i)
            {
                // Need to shift because the RealBox ctor we call takes the
                // physical location of index (0,0,0).  This does not affect
                // the usual cases where the domain index starts with 0.
                const Box& b = amrex::shift(bArray[level][i], -domain_lo);
                RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
                }
            }

            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
        }
        HeaderFile << "1" << "\n";
        HeaderFile << "3" << "\n";
        HeaderFile << "amrexvec_nu_x" << "\n";
        HeaderFile << "amrexvec_nu_y" << "\n";
        HeaderFile << "amrexvec_nu_z" << "\n";
        std::string mf_nodal_prefix = "Nu_nd";
        for (int level = 0; level <= finest_level; ++level) {
            HeaderFile << MultiFabHeaderPath(level, levelPrefix, mf_nodal_prefix) << '\n';
        }
}
#endif
