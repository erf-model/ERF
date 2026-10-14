/**
 * \file ERF_WriteBndryPlanes.H
 */
#ifndef ERF_WRITEBNDRYPLANES_H
#define ERF_WRITEBNDRYPLANES_H

#include "AMReX_Gpu.H"
#include "AMReX_AmrCore.H"
#include <AMReX_BndryRegister.H>


/** Interface for writing boundary planes
 *
 *  This class performs the necessary file operations to write boundary planes
 *
 */
class WriteBndryPlanes
{
public:
    explicit WriteBndryPlanes (amrex::Vector<amrex::BoxArray>& grids,
                               amrex::Vector<amrex::Geometry>& geom);

    void write_planes (int t_step, double time,
                       amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_new,
                       bool is_moist);

    /** Decide which rows of time.dat a run starting at step start_step (at time start_time,
     *  with start_time the absolute time written in time.dat) continues from: a restart keeps
     *  the rows the run before it wrote up to that step, a fresh start keeps none. The file is
     *  rewritten only when the run writes its first plane, so a run that stops before then
     *  leaves it as it was. Returns whether the kept rows end with a plane at start_step itself.
     */
    bool start_series (int start_step, double start_time, bool restarting);

    /** How many leading rows of a time.dat with these steps and times a run starting at
     *  start_step and start_time keeps: on a restart the rows up to and including the restart
     *  that increase in step and time and are not later than the restart; on a fresh start none.
     */
    static int rows_to_keep (const amrex::Vector<int>& steps, const amrex::Vector<double>& times,
                             int start_step, double start_time, bool restarting);

private:

    //! IO output box region
    amrex::Box target_box;

    //! Geometry objects for all levels
    amrex::Vector<amrex::Geometry>& m_geom;

    //! File name for IO
    std::string m_filename{""};

    //! File name for Native time file
    std::string m_time_file{""};

    //! Variables for IO
    amrex::Vector<std::string> m_var_names;

    //! Rows of time.dat the run continues from, written back before its first plane (I/O rank)
    amrex::Vector<std::string> m_kept_rows;

    //! Whether time.dat still has to be rewritten to m_kept_rows (I/O rank)
    bool m_rewrite_time_file = false;

    //! Timestep and times to be stored in time.dat
    amrex::Vector<double> m_in_times;
    amrex::Vector<int> m_in_timesteps;

    // This defaults to level 0, but if we can fit the whole target_box
    // inside a higher level, we will use the finest data possible
    static int bndry_lev;

    //! controls extents on native bndry output
    int m_in_rad        = 1;
    const int m_out_rad = 1;
    const int m_extent_rad = 0;
};

#endif /* ERF_BOUNDARYPLANE_H */
