#ifndef MOSTAverage_H
#define MOSTAverage_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>
#include <IndexDefines.H>

class MOSTAverage {
public:
    explicit MOSTAverage(const amrex::Vector<amrex::Geometry>& geom,
                         const amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd);
    MOSTAverage() = default;
    ~MOSTAverage()
    {
        for (int lev(0); lev<m_maxlev; ++lev){
            for (int iavg(0); iavg<m_navg; ++iavg) delete m_averages[lev][iavg];
            delete m_k_indx[lev];
        }
    };

    // Compute ncells per plane
    void set_plane_normalization();

    // Populate a 2D iMF with a uniform k value
    void set_uniform_k_indices();

    // Driver for the different average policies
    void compute_averages(int lev);

    // Fill averages for policy::plane
    void compute_plane_averages(int lev);

    // Fill averages for policy::point
    void compute_point_averages(int lev);

    // Write k indices for each average
    void write_k_indices(int lev);

    // Write averages on 2D mf
    void write_averages(int lev);

    // Get pointer to the 2D mf of averages
    const amrex::MultiFab* get_average(int lev, int comp) const { return m_averages[lev][comp]; };

protected:

    // Passed through constructor
    //--------------------------------------------
    const amrex::Vector<amrex::Geometry> m_geom;                     // Geometry at each level
    amrex::Vector<amrex::Vector<const amrex::MultiFab*>>  m_fields;  // Ptr to fields to be averaged
    amrex::Vector<amrex::MultiFab*> m_z_phys_nd;                     // Ptr to terrain heigh coords

    // General vars for multiple or all policies
    //--------------------------------------------
    int m_nvar{3};                                                   // 3 field for U/V/T
    int m_navg{4};                                                   // 4 averages for U/V/T/Umag
    int m_maxlev{0};                                                 // Total number of levels
    int m_policy{0};                                                 // Policy for type of averaging
    amrex::Real m_zref;                                              // Height for MOST BC
    std::string m_pp_prefix {"erf"};                                 // ParmParse prefix
    amrex::Vector<amrex::iMultiFab*> m_k_indx;                       // Ptr to 2D imf to hold k indices (maxlev)
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_averages;      // Ptr to 2D mf to hold averages (maxlev,navg)

    // Vars for planar average policy
    //--------------------------------------------
    amrex::Vector<amrex::Vector<int>> m_ncell_plane;                 // Number of cells in plane (maxlev,navg)
    amrex::Vector<amrex::Vector<amrex::Real>> m_plane_average;       // Plane avgs (maxlev,navg)

    // Vars for point/region average policy
    //--------------------------------------------
    int m_k_ind{-1};                                                 // k_index for point/region average
    int m_radius{0};                                                 // Radius around k_index

    // Time average w/ exponential filter fun
    //--------------------------------------------
    bool m_t_avg{false};                                             // Flag to do moving average in time
    amrex::Vector<int> m_t_init;                                     // Flag to specify if averages are initialized
    amrex::Real m_time_window{1.0e-16};                              // Width of the exp filter function
    amrex::Real m_fact_new, m_fact_old;                              // Time average factors for new and old means
};
#endif
