#ifndef ERF_RADIATION_INTERFACE_H
#define ERF_RADIATION_INTERFACE_H

#include <memory>
#include <iostream>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

// True when every box of ba spans dom in z. The radiation column sweeps need a whole column
// inside ONE box, so this tests each box rather than the level's bounding box: a level with a
// full-column box and a box aloft elsewhere has a bounding box reaching both domain ends while
// not every box holds a column. The one definition of the test, shared by every radiation model
// and by ERF::rad_level_needs_interpolation, so they cannot disagree about a level.
inline bool
rad_boxes_span_column (const amrex::BoxArray& ba, const amrex::Box& dom)
{
    for (int ib = 0; ib < ba.size(); ++ib) {
        if (ba[ib].smallEnd(2) != dom.smallEnd(2) || ba[ib].bigEnd(2) != dom.bigEnd(2)) {
            return false;
        }
    }
    return true;
}

// Basic interface class for wrapping different radiation models
class IRadiation
{
public:
    virtual ~IRadiation () = default;

    virtual
    void
    Init (const amrex::Geometry& geom,
          const amrex::BoxArray& ba,
          amrex::MultiFab* cons_in) = 0;

    virtual
    void
    Run (int& level,
         int& step,
         double& time,
         const double& dt,
         const amrex::BoxArray& ba,
         amrex::Geometry& geom,
         amrex::MultiFab* cons_in,
         amrex::iMultiFab* lmask,
         amrex::MultiFab*  t_surf,
         const amrex::Vector<const amrex::MultiFab*>& lsm_input_ptrs,
         const amrex::Vector<amrex::MultiFab*>& lsm_output_ptrs,
         amrex::MultiFab* qheating_rates,
         amrex::MultiFab* rad_fluxes,
         amrex::MultiFab* z_phys,
         amrex::MultiFab* z_sfc,
         amrex::MultiFab* lat,
         amrex::MultiFab* lon,
         const bool force_update,
         amrex::Real& declin,
         amrex::Real& calday) = 0;

    /**
     * @brief Reports whether the last radiation call updated its outputs.
     * @return True when the radiation kernel wrote output fields.
     */
    virtual bool radiation_updated () const { return true; }

    virtual
    void
    WriteDataLog (const double& time) = 0;

    // Get names input varnames for lsm
    virtual
    amrex::Vector<std::string>
    get_lsm_input_varnames ()
    {
        amrex::Vector<std::string> empty = {" "};
        return empty;
    }

    // Get names output varnames for lsm
    virtual
    amrex::Vector<std::string>
    get_lsm_output_varnames ()
    {
        amrex::Vector<std::string> empty = {" "};
        return empty;
    }

    void
    setupDataLog ()
    {
        amrex::ParmParse pp("erf.rad");
        if (pp.contains("datalog")) {
            pp.queryAdd("datalog", datalogname);

            setRecordDataInfo(datalogname);
        }
    }

    void
    setDataLogFrequency (const int nstep) {
        datalog_int = nstep;
    }

    bool
    hasDatalog () { return (datalog != nullptr && datalog->good()); }

protected:
    std::unique_ptr<std::fstream> datalog = nullptr;
    std::string datalogname;

    int datalog_int = -1; // datalog frequency in steps

    // Formatting of the datalog columns, shared by every radiation model
    static constexpr int datalog_width          = 14;
    static constexpr int datalog_precision      = 9;
    static constexpr int datalog_time_precision = 13;

private:
    void
    setRecordDataInfo (const std::string& filename)
    {
        datalog = std::make_unique<std::fstream>();
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            datalog->open(filename.c_str(),std::ios::out|std::ios::app);
            if (!datalog->good()) {
                amrex::FileOpenFailed(filename);
            }
            datalogname = filename;
        }
        amrex::ParallelDescriptor::Barrier("ERF::Radiation::setRecordDataInfo");
    }
};

#endif
