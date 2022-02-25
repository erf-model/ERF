#ifndef ERF_BOUNDARYPLANE_H
#define ERF_BOUNDARYPLANE_H

#include "AMReX_Gpu.H"
#include "AMReX_AmrCore.H"
#include <AMReX_BndryRegister.H>

/** Collection of data structures and operations for reading data
 *
 *  This class contains the inlet data structures and operations to
 *  read and interpolate inflow data.
 */
class ReadBndryPlanes
{
    using PlaneVector = amrex::Vector<amrex::FArrayBox>;

public:
    explicit ReadBndryPlanes(const amrex::Geometry& geom);

    void read_time_file();

    void read_input_files(amrex::Real time, amrex::Real dt);

    void populate_data(const int, const amrex::Real, amrex::MultiFab&, amrex::MultiFab&) const;

    void define_level_data();

    void read_file(const int idx, amrex::Vector<std::unique_ptr<PlaneVector>>& data_to_fill);

    void read_data_native(amrex::FArrayBox& xlo_plane, amrex::FArrayBox& ylo_plane,
                          amrex::FArrayBox& xhi_plane, amrex::FArrayBox& yhi_plane);
#if 0
    void read_data_native(
        const amrex::OrientationIter oit,
        amrex::BndryRegister& bndry_n,
        amrex::BndryRegister& bndry_np1,
        amrex::BndryRegister& bndry_np2,
        const int lev,
        const amrex::MultiFab*,
        const amrex::Real time,
        const amrex::Vector<amrex::Real>&);
#endif

    void interpolate(const amrex::Real);
    bool is_populated(amrex::Orientation) const;
    const amrex::FArrayBox&
    interpolate_data(const amrex::Orientation ori, const int lev) const
    {
        return (*m_data_interp[ori])[lev];
    }

    amrex::Real tinterp() const { return m_tinterp; }

private:

    //! The times for which we currently have data
    amrex::Real m_tn;
    amrex::Real m_tnp1;
    amrex::Real m_tnp2;

    //! Data at time m_tn
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_n;

    //! Data at time m_tnp1
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_np1;

    //! Data at time m_tnp2
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_np2;

    //! Data interpolated to the time requested
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_interp;

    //! Time for plane at interpolation
    amrex::Real m_tinterp{-1.0};

    //! Geometry at level 0
    amrex::Geometry m_geom;

    //! File name for IO
    std::string m_filename{""};

    //! File name for file holding timesteps and times
    std::string m_time_file{""};

    //! The timesteps / times that we read from time.dat
    amrex::Vector<amrex::Real> m_in_times;
    amrex::Vector<int> m_in_timesteps;

    //! Variables to be read in
    amrex::Vector<std::string> m_var_names;

    //! controls extents on native bndry output
    const int m_in_rad = 1;
    const int m_out_rad = 1;
    const int m_extent_rad = 0;

    int last_file_read;
};

#endif /* ERF_BOUNDARYPLANE_H */
