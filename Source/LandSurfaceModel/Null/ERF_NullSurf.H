#ifndef ERF_NULLSURF_H
#define ERF_NULLSURF_H

#include "AMReX_BCRec.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_Geometry.H"
#include "ERF_DataStruct.H"
#include "ERF_SurfacePrecipitation.H"

class NullSurf {

 public:
    NullSurf () {}

    virtual ~NullSurf () = default;

    virtual
    void
    Define (SolverChoice& /*sc*/) { }

    virtual
    void Init (const int& /*lev*/,
               const amrex::MultiFab& /*cons_in*/,
               const amrex::Geometry& /*geom*/,
               const amrex::Geometry& /*geom0*/,
               amrex::Vector<amrex::BCRec>& /*domain_bcs_type*/,
               amrex::IntVect& /*refRatio*/,
               const amrex::Real& /*dt_advance*/,
               amrex::Vector<amrex::Vector<std::string>>& /*nc_init_file*/) { }

    virtual
    void
    Advance (const double& /*dt_advance*/) { }

    virtual
    void
    Advance_With_State (const int& /*lev*/,
                        amrex::MultiFab& /*cons_in*/,
                        amrex::MultiFab& /*xvel_in*/,
                        amrex::MultiFab& /*yvel_in*/,
                        amrex::MultiFab* /*hfx3_out*/,
                        amrex::MultiFab* /*qfx3_out*/,
                        const SurfacePrecipAccumulationSources& /*precip_sources*/,
                        const amrex::Real& /*elapsed_time*/,
                        const amrex::Real& /*dt_advance*/,
                        const int& /*nstep*/,
                        const bool /*updated_lev0*/) { }

    virtual
    void
    Plot_Landfile (const int& /*nstep*/) { }

    virtual
    void
    Update_Micro_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Update_State_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_State_to_Micro (const amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_Micro_to_State (amrex::MultiFab& /*cons_in*/) { }

    virtual
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    amrex::Geometry
    Lsm_Geom ( ) { return m_lsm_geom; }

    virtual
    int
    Lsm_Data_Size () { return NullSurf::m_lsm_size; }

    virtual
    int
    Lsm_Flux_Size () { return NullSurf::m_lsm_size; }

    virtual
    std::string
    Lsm_DataName (const int& /*varIdx*/) { return dataname; }

    virtual
    std::string
    Lsm_FluxName (const int& /*varIdx*/) { return fluxname; }

    virtual
    int
    Lsm_DataIndex (std::string /*varname*/) { return (-1); }

    virtual
    int
    Lsm_FluxIndex (std::string /*varname*/) { return (-1); }

    virtual
    std::unordered_map<std::string, std::string>&
    Lsm_WRFInputNames() { return wrfinput_map; }

    virtual
    void
    Lsm_Set_Lev0_Data_Ptr (const int& /*varIdx*/,
                           amrex::MultiFab* /*lev0_data_ptr*/)
    {}

    virtual
    void
    Lsm_Set_Lev0_Flux_Ptr (const int& /*varIdx*/,
                           amrex::MultiFab* /*lev0_flux_ptr*/)
    {}

    virtual
    bool
    Get_Update_Status () { return false; }

    // Scalar LSM step counter (NoahMP itimestep). Persisted in checkpoint so
    // the LSM substepping schedule survives restart. 0 means "uninitialized".
    virtual
    int
    Get_LSM_Step () const { return 0; }

    virtual
    void
    Set_LSM_Step (int /*step*/) { }

    // Full LSM prognostic-state restart I/O. No-op for models that carry no
    // persistent prognostic state beyond their checkpointed MultiFabs.
    virtual
    void
    Write_Lsm_Restart (const std::string& /*dir*/) const { }

    virtual
    void
    Read_Lsm_Restart (const std::string& /*dir*/) { }

 private:
    int m_lsm_size = 0;
    amrex::Geometry m_lsm_geom;
    std::string dataname = " ";
    std::string fluxname = " ";

    std::unordered_map<std::string,std::string> wrfinput_map;
};
#endif
