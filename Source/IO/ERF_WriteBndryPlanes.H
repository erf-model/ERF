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

    /** Prepare time.dat for a run that starts at step start_step: a restart keeps the rows the
     *  run before it wrote up to that step, a fresh start keeps none. Returns whether the series
     *  holds a plane already, in which case the start-up plane must not be written.
     */
    bool start_series (int start_step, bool restarting);

    /** How many leading rows of a time.dat with these steps a run starting at start_step keeps:
     *  on a restart the increasing rows up to and including start_step, on a fresh start none.
     */
    static int rows_to_keep (const amrex::Vector<int>& steps, int start_step, bool restarting);

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
