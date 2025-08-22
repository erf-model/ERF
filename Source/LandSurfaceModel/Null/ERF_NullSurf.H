#ifndef ERF_NULLSURF_H
#define ERF_NULLSURF_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>
#include <ERF_DataStruct.H>

#include <ERF_LSM_Types.H>

class NullSurf {

 public:
    NullSurf () {}

    virtual ~NullSurf () = default;

    virtual
    void
    Define (SolverChoice& /*sc*/) { }

    // u_in and v_in might not be needed in all LSM models (e.g MM5), maybe refactor?
    virtual
    void Init (const int& /*lev*/,
               const amrex::MultiFab& /*cons_in*/,
               const amrex::MultiFab& /*u_in*/,
               const amrex::MultiFab& /*v_in*/,
               const amrex::Geometry& /*geom*/,
               const amrex::Real& /*dt_advance*/, 
               std::unique_ptr<amrex::MultiFab>& /*z_phys_cc*/) { } 

    virtual
    void
    Advance (const amrex::Real& /*dt_advance*/) { }

    virtual
    void
    Advance_With_State (const int& /*lev*/,
                        amrex::MultiFab& /*cons_in*/,
                        amrex::MultiFab& /*xvel_in*/,
                        amrex::MultiFab& /*yvel_in*/,
                        amrex::MultiFab* /*hfx3_out*/,
                        amrex::MultiFab* /*qfx3_out*/,
                        const amrex::Real& /*dt_advance*/,
                        const int& /*nstep*/) { }

    virtual
    void
    Update_Lsm_Vars (amrex::MultiFab& /*cons_in*/,
                     amrex::MultiFab& /*u_in*/,
                     amrex::MultiFab& /*v_in*/) { }

    virtual
    void
    Update_State_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_State_to_Lsm (const amrex::MultiFab& /*cons_in*/,
                       const amrex::MultiFab& /*u_in*/,
                       const amrex::MultiFab& /*v_in*/) { }

    virtual
    void
    Copy_Lsm_to_State (amrex::MultiFab& /*cons_in*/) { }

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
    std::string
    Lsm_VarName (const int& /*varIdx*/)
    {
        return varname;
    }

    virtual
    void
    set_flux_inputs(const amrex::MultiFab* sw_lw_fluxes_in,
                    const amrex::MultiFab* zenith_in) { }

    virtual
    void
    set_precip_input(const amrex::MultiFab* precip_in) { }

    virtual
    void
    set_terrain_inputs(const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& sst_in,
                       const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& lmask_in) { }

    virtual int Lsm_Type()
    {
        return LSM_Type::None;
    }

    virtual
    int
    Lsm_VarIndex (std::string /*varname*/)
    {
        return (-1);
    }

    virtual
    std::unordered_map<std::string, std::string>&
    Lsm_WRFInputNames() {
        return wrfinput_map;
    }

 private:
    int m_lsm_size = 1;
    amrex::Geometry m_lsm_geom;
    std::string varname = " ";

    std::unordered_map<std::string,std::string> wrfinput_map;
};
#endif
