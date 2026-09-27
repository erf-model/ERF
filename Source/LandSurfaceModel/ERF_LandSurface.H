#ifndef ERF_LANDSURFACE_H
#define ERF_LANDSURFACE_H

#include <typeinfo>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_NullSurf.H>
#include <ERF_SLM.H>

#if ERF_USE_NOAHMP
#include <ERF_NOAHMP.H>
#endif

/**
 * @brief Wrapper class for managing land surface models across AMR levels.
 */
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
          const amrex::Geometry& geom0,
          amrex::Vector<amrex::BCRec>& domain_bcs_type,
          amrex::IntVect& refRatio,
          const amrex::Real& dt_advance,
          std::unique_ptr<amrex::MultiFab>& z_phys_nd,
          amrex::Vector<amrex::Vector<std::string>>& nc_init_file)
    {
        m_lsm_model[lev]->Init(lev, cons_in, u_in, v_in, geom, geom0,
                               domain_bcs_type, refRatio,
                               dt_advance, z_phys_nd, nc_init_file);
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
        amrex::Vector<amrex::BCRec> domain_bcs_type;
        amrex::IntVect refRatio(1);
        amrex::Vector<amrex::Vector<std::string>> nc_init_file;
        Init(lev, cons_in, u_in, v_in, geom, geom, domain_bcs_type,
             refRatio, dt_advance, z_phys_nd, nc_init_file);
    }

    void
    Advance (const int& lev,
             amrex::MultiFab& cons_in,
             amrex::MultiFab& xvel_in,
             amrex::MultiFab& yvel_in,
             amrex::MultiFab* hfx3_out,
             amrex::MultiFab* qfx3_out,
             const SurfacePrecipAccumulationSources& precip_sources,
             const double& time_d,
             const double& dt_advance_d,
             const int& nstep,
             const bool updated_lev0)
    {
        amrex::Real time = static_cast<amrex::Real>(time_d);
        amrex::Real dt_advance = static_cast<amrex::Real>(dt_advance_d);
        m_lsm_model[lev]->Advance_With_State(lev, cons_in, xvel_in, yvel_in,
                                             hfx3_out, qfx3_out,
                                             precip_sources,
                                             time, dt_advance, nstep, updated_lev0);
    }

    void
    Advance (const int& lev, const double& dt_advance_d)
    {
        amrex::Real dt_advance = static_cast<amrex::Real>(dt_advance_d);
        m_lsm_model[lev]->Advance(dt_advance);
    }

    void Advance (const int& lev, const amrex::Real& dt,
                  const amrex::Real& time, const amrex::Real& start_time)
    {
        m_lsm_model[lev]->Advance(dt, time, start_time);
    }

    void
    Plot(const int& lev, const int& nstep)
    {
        m_lsm_model[lev]->Plot_Landfile(nstep);
    }

    void
    Update_Micro_Vars_Lev (const int& lev, amrex::MultiFab& cons_in)
    {
        m_lsm_model[lev]->Update_Micro_Vars(cons_in);
    }

    void Update_Lsm_Vars_Lev (const int& lev, amrex::MultiFab& cons_in,
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
    Get_Data_Ptr (const int& lev, const std::string& varname)
    {
        std::string name = varname;
        const int index = m_lsm_model[lev]->Lsm_DataIndex(name);
        return (index >= 0) ? m_lsm_model[lev]->Lsm_Data_Ptr(index) : nullptr;
    }

    bool
    Has_Model () const { return !m_lsm_model.empty() && m_lsm_model[0] != nullptr; }

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

    amrex::Vector<std::string>
    Get_DataNames ()
    {
        amrex::Vector<std::string> names;
        if (m_lsm_model.empty() || m_lsm_model[0] == nullptr) {
            return names;
        }
        const int nvars = Get_Data_Size();
        names.reserve(nvars);
        for (int var = 0; var < nvars; ++var) {
            names.push_back(Get_DataName(var));
        }
        return names;
    }

    int
    Get_DataIdx (const int& lev, const std::string& varname) { return m_lsm_model[lev]->Lsm_DataIndex(varname); }

    /**
     * Registers a source level's exposed state with a destination level.
     */
    void
    Set_Source_Ptrs (const int& lev, const int& source_lev)
    {
        AMREX_ALWAYS_ASSERT(lev > 0);
        AMREX_ALWAYS_ASSERT(source_lev >= 0);
        AMREX_ALWAYS_ASSERT(m_lsm_model[lev]->Lsm_Type() == m_lsm_model[source_lev]->Lsm_Type());

        const int ndata = m_lsm_model[lev]->Lsm_Data_Size();
        AMREX_ALWAYS_ASSERT(ndata == m_lsm_model[source_lev]->Lsm_Data_Size());
        for (int varIdx(0); varIdx<ndata; ++varIdx) {
            m_lsm_model[lev]->Lsm_Set_Source_Data_Ptr(
                varIdx, m_lsm_model[source_lev]->Lsm_Data_Ptr(varIdx));
        }

        const int nflux = m_lsm_model[lev]->Lsm_Flux_Size();
        AMREX_ALWAYS_ASSERT(nflux == m_lsm_model[source_lev]->Lsm_Flux_Size());
        for (int varIdx(0); varIdx<nflux; ++varIdx) {
            m_lsm_model[lev]->Lsm_Set_Source_Flux_Ptr(
                varIdx, m_lsm_model[source_lev]->Lsm_Flux_Ptr(varIdx));
        }

        const int nunmapped = m_lsm_model[lev]->Lsm_Unmapped_Size();
        AMREX_ALWAYS_ASSERT(nunmapped == m_lsm_model[source_lev]->Lsm_Unmapped_Size());
        for (int varIdx(0); varIdx<nunmapped; ++varIdx) {
            m_lsm_model[lev]->Lsm_Set_Source_Unmapped_Ptr(
                varIdx, m_lsm_model[source_lev]->Lsm_Unmapped_Ptr(varIdx));
        }
    }

    /**
     * Initializes a level from a source level using the requested state mode.
     */
    void
    Initialize_From_Source (const int& lev, const int& source_lev,
                            const amrex::IntVect& ref_ratio,
                            const LSMTransferMode mode)
    {
        Set_Source_Ptrs(lev, source_lev);
        m_lsm_model[lev]->Lsm_Interpolate_From_Source(
            m_lsm_model[source_lev]->Lsm_Atmos_Geom(),
            m_lsm_model[source_lev]->Lsm_Geom(), ref_ratio, mode);
    }

    /**
     * Rebuilds a land-surface level after ERF remakes its atmospheric grid.
     */
    void
    Remake_Level (const int& lev,
                  const amrex::MultiFab& cons_in,
                  const amrex::MultiFab& u_in,
                  const amrex::MultiFab& v_in,
                  const amrex::Geometry& geom,
                  const amrex::Geometry& geom0,
                  amrex::Vector<amrex::BCRec>& domain_bcs_type,
                  amrex::IntVect& ref_ratio,
                  const amrex::Real& dt_advance,
                  std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                  amrex::Vector<amrex::Vector<std::string>>& nc_init_file,
                  const int source_lev = -1,
                  const amrex::IntVect& source_ref_ratio = amrex::IntVect(1))
    {
        const bool has_source = source_lev >= 0;
        if (has_source) {
            Set_Source_Ptrs(lev, source_lev);
        }

        const amrex::Geometry& source_geom = has_source
            ? m_lsm_model[source_lev]->Lsm_Atmos_Geom() : geom;
        const amrex::Geometry& source_lsm_geom = has_source
            ? m_lsm_model[source_lev]->Lsm_Geom() : geom;
        m_lsm_model[lev]->Lsm_Regrid_Level(
            lev, cons_in, u_in, v_in, geom, geom0, domain_bcs_type,
            ref_ratio, dt_advance, z_phys_nd, nc_init_file,
            source_geom, source_lsm_geom, source_ref_ratio, has_source);
    }

    void
    Set_Lev0_Data_Ptr (const int& lev) {
        // Ensure the same type of model for copy in place of flux data
        if ( (lev > 0) && typeid(m_lsm_model[lev]) == typeid(m_lsm_model[0]) ) {
            int ndata = this->Get_Data_Size();
            for (int varIdx(0); varIdx<ndata; ++varIdx) {
                amrex::MultiFab* Lev0_Data_Ptr = this->Get_Data_Ptr(0, varIdx);
                m_lsm_model[lev]->Lsm_Set_Lev0_Data_Ptr(varIdx, Lev0_Data_Ptr);
            }
        } else {
            amrex::Abort("Differing types of land surface models across levels is not supported!");
        }
    }

    std::string
    Get_FluxName (const int& varIdx) { return m_lsm_model[0]->Lsm_FluxName(varIdx); }

    int
    Get_FluxIdx (const int& lev, const std::string& varname) { return m_lsm_model[lev]->Lsm_FluxIndex(varname); }

    void
    Set_Lev0_Flux_Ptr (const int& lev) {
        // Ensure the same type of model for copy in place of flux data
        if ( (lev > 0) && typeid(m_lsm_model[lev]) == typeid(m_lsm_model[0]) ) {
            int nflux = this->Get_Flux_Size();
            for (int varIdx(0); varIdx<nflux; ++varIdx) {
                amrex::MultiFab* Lev0_Flux_Ptr = this->Get_Flux_Ptr(0, varIdx);
                m_lsm_model[lev]->Lsm_Set_Lev0_Flux_Ptr(varIdx, Lev0_Flux_Ptr);
            }
        } else {
            amrex::Abort("Differing types of land surface models across levels is not supported!");
        }
    }

    bool
    Get_LSM_Update_Status (const int& lev) { return m_lsm_model[lev]->Get_Update_Status(); }

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
    Plot_Lsm_Data (const PlotFileType plotfile_type,
                   double time,
                   const amrex::Vector<int> &level_steps,
                   const amrex::Vector<amrex::IntVect> &ref_ratio, const int finest_level)
    {
        // Note we must use finest_level here rather than m_lsm_model.size(), which is
        //    (max_level+1) -- the models at levels that have not (yet) been created
        //    have been allocated but not Define'd or Init'ed, so their data pointers
        //    are still null
        int nlev = finest_level + 1;
        int nlev_valid = finest_level + 1;
        int nvar = this->Get_Data_Size();
        amrex::Vector<std::string> varnames;

        // Only write if we have valid pointers
        if (this->Get_Data_Ptr(0,0)) {
            varnames.resize(nvar);
            m_lsm_geom_lev.resize(nlev);
            m_lsm_data_lev.resize(nlev);
            std::string plotfilename = amrex::Concatenate(plot_file_lsm, level_steps[0], 5);
            for (int lev = 0; lev <= finest_level; ++lev) {
                m_lsm_geom_lev[lev] = this->Get_Lsm_Geom(lev);

                if (m_lsm_geom_lev[lev].Domain().isEmpty()) {
                    nlev_valid--;
                    continue;
                }

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
            if (plotfile_type == PlotFileType::Amrex
#ifdef ERF_USE_NETCDF
                || plotfile_type == PlotFileType::Netcdf
#endif
                ) {
#ifdef ERF_USE_NETCDF
                // NetCDF output of the land-surface fields is not implemented yet, so fall back
                // to the AMReX format -- which is what these plotfiles were before the plotfile
                // type was consulted here -- rather than silently writing nothing at all.
                warn_nc_lsm_plot_fallback(plotfile_type);
#endif
                WriteMultiLevelPlotfile (plotfilename, nlev_valid, GetVecOfConstPtrs(m_lsm_data_lev),
                                            varnames, m_lsm_geom_lev, static_cast<amrex::Real>(time),
                                            level_steps, ref_ratio);
            } else {
                amrex::Abort("Dont know this plotfile_type");
            }

            m_lsm_geom_lev.clear();
            m_lsm_data_lev.clear();
        }

        if (nlev_valid > 0 && m_lsm_model[0]->Lsm_Type() == LSM_Type::SLM)
        {
            amrex::Vector<amrex::MultiFab> fab_2d(finest_level+1);
            amrex::Vector<amrex::Geometry> geom_2d(finest_level+1);
            amrex::Vector<std::string> varnames_2d;
            std::string plotfilename_2d = amrex::Concatenate(plot_file_lsm + "2D_", level_steps[0], 5);
            for (int lev = 0; lev <= finest_level; ++lev) {
                static_cast<SLM*>(m_lsm_model[lev].get())->writeSLM_Data(plotfile_type, static_cast<amrex::Real>(time), plot_file_lsm, level_steps[0], lev, finest_level, fab_2d[lev], geom_2d[lev], varnames_2d);
            }
            // For a NetCDF run writeSLM_Data has already written these fields itself, one file
            // per level, so only the AMReX format is assembled here.
            if (plotfile_type == PlotFileType::Amrex) {
                amrex::WriteMultiLevelPlotfile(plotfilename_2d, finest_level+1, GetVecOfConstPtrs(fab_2d), varnames_2d, geom_2d, static_cast<amrex::Real>(time), level_steps, ref_ratio);
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
#ifdef ERF_USE_NETCDF
    /**
     * Tell the user once per run that the land-surface plotfile is written in the AMReX
     * format even though erf.plotfile_type is netcdf, because there is no NetCDF writer
     * for these fields.
     */
    void warn_nc_lsm_plot_fallback (const PlotFileType plotfile_type)
    {
        if (plotfile_type != PlotFileType::Netcdf || m_warned_nc_lsm_plot) { return; }
        m_warned_nc_lsm_plot = true;
        amrex::Print() << "NOTE: erf.plotfile_type = netcdf, but NetCDF output is not implemented "
                          "for the land-surface fields; writing " << plot_file_lsm
                       << "* in the AMReX plotfile format instead." << std::endl;
    }
#endif

    // lsm model at each level
    amrex::Vector<std::unique_ptr<NullSurf>> m_lsm_model;

    // plotfile prefix
    std::string plot_file_lsm {"plt_lsm_"};

    // Vector of geometry
    amrex::Vector<amrex::Geometry> m_lsm_geom_lev;

    // Vector of data pointers
    amrex::Vector<amrex::MultiFab> m_lsm_data_lev;

#ifdef ERF_USE_NETCDF
    // Whether we have already told the user that LSM plotfiles fall back to the AMReX format
    [[maybe_unused]] bool m_warned_nc_lsm_plot = false;
#endif
};
#endif
