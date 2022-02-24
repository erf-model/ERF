#ifndef ERF_BOUNDARYPLANE_H
#define ERF_BOUNDARYPLANE_H

#include "AMReX_Gpu.H"
#include "AMReX_AmrCore.H"
#include <AMReX_BndryRegister.H>

struct BufferData
{
    amrex::Gpu::ManagedVector<amrex::Real> data;
    amrex::Vector<size_t> start{0, 0, 0, 0};
    amrex::Vector<size_t> count{0, 0, 0, 0};
};

/** Collection of data structures and operations for reading data
 *  \ingroup we_abl
 *
 *  This class contains the inlet data structures and operations to
 *  read and interpolate inflow data.
 */
class ReadBndryPlanes
{
    using PlaneVector = amrex::Vector<amrex::FArrayBox>;

public:
    explicit ReadBndryPlanes(const amrex::Geometry& geom); 

    void read_header();

    void read_input_files(amrex::Real time, amrex::Real dt);

    void populate_data(const int, const amrex::Real, amrex::MultiFab&, amrex::MultiFab&) const;

    bool is_initialized() const { return m_is_initialized; }

    bool box_intersects_boundary(
        const amrex::Box&, const int, const amrex::Orientation) const;

    amrex::Real m_time;

    void define_level_data(
        const amrex::Orientation, const amrex::Box&, const size_t);

    void read_file(const int idx);

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

    int& component(const int id) { return m_components[id]; }

    int component(const int id) const { return m_components.at(id); }

    int nlevels(const amrex::Orientation ori) const
    {
        return (*m_data_interp[ori]).size();
    }

    amrex::Real tn() const { return m_tn; }
    amrex::Real tnp1() const { return m_tnp1; }
    amrex::Real tnp2() const { return m_tnp2; }
    amrex::Real tinterp() const { return m_tinterp; }

private:
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_n;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_np1;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_np2;
    amrex::Vector<std::unique_ptr<PlaneVector>> m_data_interp;

    //! Time for plane at n
    amrex::Real m_tn{-1.0};

    //! Time for plane at n + 1
    amrex::Real m_tnp1{-1.0};

    //! Time for plane at n + 2
    amrex::Real m_tnp2{-1.0};

    //! Time for plane at interpolation
    amrex::Real m_tinterp{-1.0};

    //! Map of `{variableId : component}`
    std::unordered_map<int, int> m_components;

    //! Plane names
    amrex::Vector<std::string> m_plane_names{
        {"xlo", "ylo", "zlo", "xhi", "yhi", "zhi"}};

    //! IO boundary planes
    amrex::Vector<std::string> m_planes;

    //! IO output box region
    amrex::Box target_box;

    //! Geometry at level 0
    amrex::Geometry m_geom;

    //! File name for IO
    std::string m_filename{""};

    //! File name for Native time file
    std::string m_time_file{""};

    //! Variables for IO
    amrex::Vector<std::string> m_var_names;

    //! NetCDF input level time
    amrex::Vector<amrex::Real> m_in_times;
    amrex::Vector<int> m_in_timesteps;

    //! Flag indicating if this capability is available
    bool m_is_initialized{false};

    //! controls extents on native bndry output
    const int m_in_rad = 1;
    const int m_out_rad = 1;
    const int m_extent_rad = 0;

    int last_file_read;
};

#endif /* ERF_BOUNDARYPLANE_H */
