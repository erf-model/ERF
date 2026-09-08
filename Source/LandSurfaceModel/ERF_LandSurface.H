#ifndef ERF_LANDSURFACE_H
#define ERF_LANDSURFACE_H

#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_NullSurf.H>
#include <ERF_SLM.H>
#include <ERF_MM5.H>
#include <ERF_LSM_Types.H>
#include <ERF_OceanSurf.H>

#if ERF_USE_NOAHMP
#include <ERF_NOAHMP.H>
#endif

class LandSurface {

public:

    LandSurface () { }

    ~LandSurface () = default;

    void
    ReSize (const int& nlev) { m_lsm_model.resize(nlev); }

    template<class NewSurfModel>
    void
    SetModel ()
    {
        for (int lev(0); lev<m_lsm_model.size(); ++lev) {
            m_lsm_model[lev] = std::make_unique<NewSurfModel>();
        }
    }

    void
    Define (const int& lev,
            SolverChoice& sc)
    {
        m_lsm_model[lev]->Define(sc);
    }

    void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::MultiFab& u_in,
          const amrex::MultiFab& v_in,
          const amrex::Geometry& geom,
          const amrex::Real& dt_advance,
          std::unique_ptr<amrex::MultiFab>& z_phys_nd)
    {
        m_lsm_model[lev]->Init(lev, cons_in, u_in, v_in, geom, dt_advance, z_phys_nd);
    }

    void
    Advance (const int& lev,
             amrex::MultiFab& cons_in,
             amrex::MultiFab& xvel_in,
             amrex::MultiFab& yvel_in,
             amrex::MultiFab* hfx3_out,
             amrex::MultiFab* qfx3_out,
             const amrex::Real& time,
             const amrex::Real& dt_advance,
             const int& nstep)
    {
        m_lsm_model[lev]->Advance_With_State(lev, cons_in, xvel_in, yvel_in,
                                             hfx3_out, qfx3_out, time, dt_advance, nstep);
    }

    void
    Advance (const int& lev, const amrex::Real& dt_advance, const amrex::Real &time, const amrex::Real &start_time)
    {
        m_lsm_model[lev]->Advance(dt_advance, time, start_time);
    }

    void
    Plot(const int& lev, const int& nstep)
    {
        m_lsm_model[lev]->Plot_Landfile(nstep);
    }

    void
    Update_Lsm_Vars_Lev (const int& lev, amrex::MultiFab& cons_in,
                         amrex::MultiFab& u_in, amrex::MultiFab& v_in)
    {
        m_lsm_model[lev]->Update_Lsm_Vars(cons_in, u_in, v_in);
    }

    void
    Update_State_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        m_lsm_model[lev]->Update_State_Vars(cons_in);
    }

    amrex::MultiFab*
    Get_Data_Ptr (const int& lev, const int& varIdx) { return m_lsm_model[lev]->Lsm_Data_Ptr(varIdx); }

    amrex::MultiFab*
    Get_Flux_Ptr (const int& lev, const int& varIdx) { return m_lsm_model[lev]->Lsm_Flux_Ptr(varIdx); }

    amrex::Geometry
    Get_Lsm_Geom (const int& lev ) { return m_lsm_model[lev]->Lsm_Geom(); }

    int
    Get_Data_Size () { return m_lsm_model[0]->Lsm_Data_Size(); }

    int
    Get_Flux_Size () { return m_lsm_model[0]->Lsm_Flux_Size(); }

    std::string
    Get_DataName (const int& varIdx) { return m_lsm_model[0]->Lsm_DataName(varIdx); }

    int
    Get_DataIdx (const int& lev, const std::string& varname) { return m_lsm_model[lev]->Lsm_DataIndex(varname); }

    std::string
    Get_FluxName (const int& varIdx) { return m_lsm_model[0]->Lsm_FluxName(varIdx); }

    int
    Get_FluxIdx (const int& lev, const std::string& varname) { return m_lsm_model[lev]->Lsm_FluxIndex(varname); }

    std::unordered_map<std::string,std::string>&
    Get_WRFInputNames() { return m_lsm_model[0]->Lsm_WRFInputNames(); }

    // Scalar LSM step counter (e.g., NoahMP itimestep) for checkpoint/restart.
    int
    Get_LSM_Step (const int& lev) const { return m_lsm_model[lev]->Get_LSM_Step(); }

    void
    Set_LSM_Step (const int& lev, int step) { m_lsm_model[lev]->Set_LSM_Step(step); }

    // Full LSM prognostic-state restart I/O (e.g., NoahMP soil/snow/canopy).
    // Write is const (called from the const WriteCheckpointFile); it only
    // serializes state, it does not mutate the ERF-side model wrapper.
    void
    Write_Lsm_Restart (const int& lev, const std::string& dir) const { m_lsm_model[lev]->Write_Lsm_Restart(dir); }

    void
    Read_Lsm_Restart (const int& lev, const std::string& dir) { m_lsm_model[lev]->Read_Lsm_Restart(dir); }

    void
    Plot_Lsm_Data (const PlotFileType plotfile_type, amrex::Real time,
                   const amrex::Vector<int> &level_steps,
                   const amrex::Vector<amrex::IntVect> &ref_ratio, const int finest_level)
    {
        int nlev = m_lsm_model.size();
        int nvar = this->Get_Data_Size();
        int nlev_valid = finest_level+1;
        amrex::Vector<std::string> varnames;

        // // Only write if we have valid pointers
        // if (this->Get_Data_Ptr(0,0)) {
        varnames.resize(nvar);
        m_lsm_geom_lev.resize(nlev);
        m_lsm_data_lev.resize(nlev);
        std::string plotfilename = amrex::Concatenate(plot_file_lsm, level_steps[0], 5);
        // for (int lev(0); lev<nlev; ++lev) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_lsm_geom_lev[lev] = this->Get_Lsm_Geom(lev);

            amrex::MultiFab* mf_lsm = this->Get_Data_Ptr(lev,0);
            amrex::IntVect ng(0,0,0);
            amrex::BoxArray ba = mf_lsm->boxArray();
            amrex::DistributionMapping dm = mf_lsm->DistributionMap();
            m_lsm_data_lev[lev].define(ba, dm, nvar, ng);
            for (int n(0); n<nvar; ++n) {
                mf_lsm = this->Get_Data_Ptr(lev,n);
                amrex::MultiFab::Copy(m_lsm_data_lev[lev],*(mf_lsm),0,n,1,0);
                if (lev==0) varnames[n] = this->Get_DataName(n);
            }
        }
        if (plotfile_type == PlotFileType::Amrex) {
            WriteMultiLevelPlotfile (plotfilename, nlev, GetVecOfConstPtrs(m_lsm_data_lev),
                                        varnames, m_lsm_geom_lev, time, level_steps, ref_ratio);
#ifdef ERF_USE_NETCDF
        } else if (plotfile_type == PlotFileType::Netcdf) {
            int lev   = 0;
            int l_which = 0;
            //writeNCPlotFile(lev, l_which, plotfilename, GetVecOfConstPtrs(m_lsm_data_lev), varnames, level_steps, time);
#endif
        } else {
            amrex::Abort("Dont know this plot_filetype");
        }

        m_lsm_geom_lev.clear();
        m_lsm_data_lev.clear();

        if (nlev > 0 && m_lsm_model[0]->Lsm_Type() == LSM_Type::SLM)
        {
            // for (int lev(0); lev<nlev_valid; ++lev) {
            amrex::Vector<amrex::MultiFab> fab_2d(finest_level+1);
            amrex::Vector<amrex::Geometry> geom_2d(finest_level+1);
            amrex::Vector<std::string> varnames;
            std::string plotfilename_2d = amrex::Concatenate(plot_file_lsm + "2D_", level_steps[0], 5);
            for (int lev = 0; lev <= finest_level; ++lev) {
                static_cast<SLM*>(m_lsm_model[lev].get())->writeSLM_Data(plotfile_type, time, plot_file_lsm, level_steps[0], lev, finest_level, fab_2d[lev], geom_2d[lev], varnames);
            }
            if (plotfile_type == PlotFileType::Amrex) {
                amrex::WriteMultiLevelPlotfile(plotfilename_2d, finest_level+1, GetVecOfConstPtrs(fab_2d), varnames, geom_2d, time, level_steps, ref_ratio);
            }
        }
    }

    void set_LSM_precip_input (const int& lev,
                               const amrex::MultiFab* precip_in)
    {
        m_lsm_model[lev]->set_precip_input(precip_in);
    }

    void set_LSM_terrain_inputs(const int& lev,
                                const amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab>>>& sst_in,
                                const amrex::Vector<amrex::Vector<std::unique_ptr<amrex::iMultiFab>>>& lmask_in)
    {
        m_lsm_model[lev]->set_terrain_inputs(sst_in[lev], lmask_in[lev]);
    }

    template<class SurfModelType>
    SurfModelType* get_model_lev(const int &lev) {
        return static_cast<SurfModelType*>(m_lsm_model[lev].get());
    }

    void WriteCheckpoint(const int &lev, const std::string &checkpointname) const
    {
        m_lsm_model[lev]->WriteCheckpoint(lev, checkpointname);
    }

    void ReadCheckpoint(const int &lev, const std::string &checkpointname)
    {
        m_lsm_model[lev]->ReadCheckpoint(lev, checkpointname);
    }

private:
    // lsm model at each level
    amrex::Vector<std::unique_ptr<NullSurf>> m_lsm_model;

    // plotfile prefix
    std::string plot_file_lsm {"plt_lsm_"};

    // Vector of geometry
    amrex::Vector<amrex::Geometry> m_lsm_geom_lev;

    // Vector of data pointers
    amrex::Vector<amrex::MultiFab> m_lsm_data_lev;
};
#endif
