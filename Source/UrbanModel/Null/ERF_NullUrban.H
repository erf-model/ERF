#ifndef ERF_NULLURBAN_H
#define ERF_NULLURBAN_H

#include <memory>
#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_DataStruct.H>

class NullUrban {

public:

    NullUrban () = default;

    virtual ~NullUrban () = default;

    /**
     * @brief Defines a no-op urban model.
     */
    virtual void Define (SolverChoice& /*sc*/) { }

    /**
     * @brief Initializes a no-op urban model.
     */
    virtual void Init (const amrex::MultiFab& /*cons_in*/,
                       const amrex::MultiFab& /*u_in*/,
                       const amrex::MultiFab& /*v_in*/,
                       const amrex::Geometry& /*geom*/,
                       const amrex::Real& /*dt_advance*/,
                       const amrex::MultiFab& /*z_phys_nd*/,
                       const amrex::iMultiFab& /*land_type*/,
                       const amrex::MultiFab& /*urb_frac*/) { }

    /**
     * @brief Advances a no-op urban model.
     */
    virtual void Advance (const amrex::Real& /*dt_advance*/,
                          const amrex::Real& /*time*/,
                          const amrex::Real& /*start_time*/,
                          const amrex::Geometry& /*geom*/,
                          const amrex::MultiFab* /*z_phys_nd*/,
                          amrex::MultiFab& /*eddyDiffs_lev*/,
                          const amrex::Real& /*calday*/) { }

    /**
     * @brief Updates no urban coupling fields.
     */
    virtual void Update_Urban_Vars (amrex::MultiFab& /*cons_in*/,
                                    amrex::MultiFab& /*u_in*/,
                                    amrex::MultiFab& /*v_in*/) { }

    /**
     * @brief Updates no atmospheric state fields.
     */
    virtual void Update_State_Vars (amrex::MultiFab& /*cons_in*/) { }

    /**
     * @brief Returns no urban data field.
     */
    virtual amrex::MultiFab* Urban_Data_Ptr (const int& /*varIdx*/) { return nullptr; }

    /**
     * @brief Returns no urban flux field.
     */
    virtual amrex::MultiFab* Urban_Flux_Ptr (const int& /*varIdx*/) { return nullptr; }

    /**
     * @brief Returns the null-model geometry.
     */
    virtual amrex::Geometry Urban_Geom () { return m_urban_geom; }

    /**
     * @brief Returns the number of null-model fields.
     */
    virtual int Urban_Data_Size () { return m_urban_size; }

    /**
     * @brief Returns a null-model field name.
     */
    virtual std::string Urban_DataName (const int& /*varIdx*/) { return varname; }

    /**
     * @brief Returns no field index.
     */
    virtual int Urban_DataIndex (const std::string& /*varname*/) { return -1; }

    /**
     * @brief Accepts no flux input.
     */
    virtual void set_flux_inputs (const amrex::Real& /*declin*/) { }

    /**
     * @brief Accepts no terrain input.
     */
    virtual void set_terrain_inputs (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& /*tsk_in*/,
        const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& /*lmask_in*/,
        const amrex::Vector<std::unique_ptr<amrex::iMultiFab>>& /*land_type_in*/,
        const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& /*urb_frac_in*/,
        amrex::MultiFab* /*lat_in*/, amrex::MultiFab* /*lon_in*/) { }

    /**
     * @brief Writes no checkpoint state.
     */
    virtual void WriteCheckpoint (const int& /*lev*/, const std::string& /*checkpointname*/) const { }

    /**
     * @brief Reads no checkpoint state.
     */
    virtual void ReadCheckpoint (const int& /*lev*/, const std::string& /*checkpointname*/) { }

private:
    int m_urban_size = 1;
    amrex::Geometry m_urban_geom;
    std::string varname = " ";
};

#endif
