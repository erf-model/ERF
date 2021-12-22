/** \addtogroup Utilities
 * @{
 *
 * \file IOManager.H
 */

#ifndef IOMANAGER_H
#define IOMANAGER_H

#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>
#include <unordered_map>
#include <set>
#include <memory>

#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_buildInfo.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_Amr.H>


// forwarded class
class ERF;

/**
 * \brief Input/Output manager
 *
 *  This class contains infrastructure to output solution fields during the
 *  simulation as well as code that reads a checkpoint file during a restarted
 *  simulation. By default the \ref eqsys register the necessary fields during
 *  initialization and these fields are output into plot-file and checkpoint
 *  file formats based on output frequencies managed by SimTime. The user can
 *  request additional fields be output by setting appropriate parameters in the
 *  input file. The class also provides the ability to override output of the
 *  default fields and output a subset of those fields.
 */
class IOManager {
public:
    explicit IOManager(ERF&);

    ~IOManager();

    //! Write all user-requested fields to disk
    void writePlotFile(const std::string& dir, std::ostream& os, amrex::VisMF::How how);

    //! write small plotfile
    void writeSmallPlotFile(const std::string& dir, std::ostream& os, amrex::VisMF::How how);

    //! Write checkpoint file
    void checkPoint(const std::string& dir, std::ostream& os, amrex::VisMF::How how, bool dump_old_default);

    //! restart using checkpoint file
    void restart(amrex::Amr& papa, std::istream& is, bool bReadSpecial);

    //! set state in checkpoint
    void set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint);

    //! setup the plot variables
    void setPlotVariables();

    //! write job information
    void writeJobInfo(const std::string& dir);

    //! Create 1D vertical column output for coupling
    void createNCColumnFile(const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc);

    //! Write a timestep to 1D vertical column output for coupling
    void writeToNCColumnFile(const std::string& colfile_name, const amrex::Real xloc, const amrex::Real yloc);

    //! write build information
    static void writeBuildInfo(std::ostream& os);

#ifdef ERF_USE_NETCDF
    //! Write all user-requested fields to disk
    void writeNCPlotFile(const std::string& dir, std::ostream& os);

    //! Write all necessary fields for restart
    void NCWriteCheckpointFile (const std::string& dir, std::ostream& os, const bool& dump_old);

    //! restart using nc checkpoint file
    void ncrestart(amrex::Amr& papa, std::istream& is, bool bReadSpecial);
#endif

private:
    // ERF class that iomanager belongs to
    ERF& erf;
};

#endif /* IOMANAGER_H */
/** @}*/
