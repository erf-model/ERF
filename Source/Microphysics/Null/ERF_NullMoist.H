#ifndef ERF_NULLMOIST_H
#define ERF_NULLMOIST_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Geometry.H>
#include <ERF_DataStruct.H>

class NullMoist {

public:
    NullMoist () {}

    virtual ~NullMoist () = default;

    virtual
    void
    Define (SolverChoice& /*sc*/) { }

    virtual
    void Init (const amrex::MultiFab& /*cons_in*/,
               const amrex::BoxArray& /*grids*/,
               const amrex::Geometry& /*geom*/,
               const amrex::Real& /*dt_advance*/,
               std::unique_ptr<amrex::MultiFab>& /*z_phys_nd*/,
               std::unique_ptr<amrex::MultiFab>& /*detJ_cc*/) { }

    virtual
    void
    Advance (const amrex::Real& /*dt_advance*/,
             const SolverChoice& /*solverChoce*/) { }

    virtual
    void
    Update_Micro_Vars (amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Update_State_Vars (amrex::MultiFab& /*cons_in*/,
                       const amrex::MultiFab& /*z_phys_nd*/) { }

    virtual
    void
    Copy_State_to_Micro (const amrex::MultiFab& /*cons_in*/) { }

    virtual
    void
    Copy_Micro_to_State (amrex::MultiFab& /*cons_in*/) { }

    virtual
    amrex::MultiFab*
    Qmoist_Ptr (const int& /*varIdx*/ ) { return nullptr; }

    virtual
    int
    Qmoist_Size () { return NullMoist::m_qmoist_size; }

    virtual
    int
    Qstate_Moist_Size () { return NullMoist::m_qstate_moist_size; }

    virtual
    int
    Qstate_Moist_NumConc_Size () { return NullMoist::m_qstate_moist_numconc_size; }

    virtual
    int
    Qstate_NonMoist_Size () { return NullMoist::m_qstate_nonmoist_size; }

    virtual
    void
    Qmoist_Restart_Vars ( const SolverChoice&, std::vector<int>& a_idx, std::vector<std::string>& a_names) const
    {
        a_idx.clear();
        a_names.clear();
    }

    virtual
    void
    GetPlotVarNames (amrex::Vector<std::string>& a_vec) const
    {
        a_vec.resize(0);
    }

    virtual
    void
    GetPlotVar (const std::string& /* a_name */,
                amrex::MultiFab& /* a_mf */ ) const
    {
        amrex::Abort("NullMoist::GetPlotVar() should never be called.");
    }

    virtual
    void
    GetPlotVar (const std::string& a_name,
                amrex::MultiFab& a_mf,
                const int /*a_lev*/) const
    {
        GetPlotVar(a_name, a_mf);
    }

    virtual
    void
    SetCurrentLevel (const int&) {}

    virtual
    void
    InitLevel (const int, const amrex::MultiFab&) {}

    virtual
    int
    getDiagnosticsInterval () const { return 0; }

    virtual
    void
    Set_dzmin (const amrex::Real /*dz_min*/) { }

    virtual
    void
    // Current SatAdj/Kessler development paths do not trim by real width, but
    // tests still exercise Set_RealWidth(0/1/2) as a future-proof contract.
    // Keep the base hook as a no-op unless a concrete microphysics model
    // chooses to implement width-sensitive behavior.
    Set_RealWidth (const int /*real_width*/) { }


private:
    int m_qmoist_size = 0;
    int m_qstate_moist_size = 0;
    int m_qstate_moist_numconc_size = 0;
    int m_qstate_nonmoist_size = 0;
};

#endif
