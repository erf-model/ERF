#ifndef ERF_RADIATION_INTERFACE_H
#define ERF_RADIATION_INTERFACE_H

#include <memory>
#include <iostream>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

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
         amrex::Real& time,
         const amrex::Real& dt,
         const amrex::BoxArray& ba,
         amrex::Geometry& geom,
         amrex::MultiFab* cons_in,
         amrex::iMultiFab* lmask,
         amrex::MultiFab*  t_surf,
         amrex::MultiFab* lsm_fluxes,
         amrex::MultiFab* lsm_zenith,
         amrex::Vector<amrex::MultiFab*>& lsm_input_ptrs,
         amrex::Vector<amrex::MultiFab*>& lsm_output_ptrs,
         amrex::MultiFab* qheating_rates,
         amrex::MultiFab* rad_fluxes,
         amrex::MultiFab* z_phys,
         amrex::MultiFab* lat,
         amrex::MultiFab* lon) = 0;

    virtual
    void
    WriteDataLog (const amrex::Real& time) = 0;

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
            pp.query("datalog", datalogname);

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
