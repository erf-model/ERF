#ifndef MOSTAverage_H
#define MOSTAverage_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

class MOSTAverage {
public:
    explicit MOSTAverage(amrex::Vector<const amrex::MultiFab*> fields,
                         amrex::Vector<amrex::MultiFab*> averages,
                         amrex::Vector<amrex::iMultiFab*> k_indx,
                         amrex::Vector<amrex::Real> z_ref, 
                         const amrex::MultiFab* z_nd,
                         const amrex::Geometry geom,
                         int policy,
                         bool update_k);
     MOSTAverage() = default;
    ~MOSTAverage() = default;

    // Populate a 2D iMF with the k indices for averaging
    void compute_plane_k_indices();

    // Driver for the different average policies
    void compute_averages();

    // Fill plane storage with averages
    void compute_plane_averages();

    // Get number of components for a mf in m_fields
    int ncomp()      const  { return m_ncomps[0]; }
    int ncomp(int i) const  { return m_ncomps[i]; }

    // Get pointer to mf in m_fields
    const amrex::MultiFab* get_field()      const { return m_fields[0]; }
    const amrex::MultiFab* get_field(int i) const { return m_fields[i]; }

    // Get pointer to the 2D mf of averages
    const amrex::MultiFab* get_average()      const { return m_averages[0]; };
    const amrex::MultiFab* get_average(int i) const { return m_averages[i]; };

    // Get plane averages long dir()
    const amrex::Vector<amrex::Real>& get_plane_average()      const { return m_plane_average[0]; }
    const amrex::Vector<amrex::Real>& get_plane_average(int i) const { return m_plane_average[i]; }

    // Write k indices for each average
    void write_k_indices();
    
    // Write averages on 2D mf
    void write_averages();
        
protected:
    amrex::Real m_dz;  // Mesh spacing in axis direction
    amrex::Real m_zlo; // Bottom of domain in axis direction

    amrex::Vector<int> m_ncomps;      // Number of components (mf)
    amrex::Vector<int> m_ncell_plane; // Number of cells in plane (mf)

    amrex::Vector<amrex::Vector<amrex::Real>> m_plane_average; // Plane avgs (mf,2)

    amrex::Vector<const amrex::MultiFab*> m_fields;  // Ptr to fields to be averaged
    amrex::Vector<amrex::MultiFab*> m_averages;      // Ptr to 2D mf to hold averages
    amrex::Vector<amrex::iMultiFab*> m_k_indx;       // Ptr to 2D imf to hold k indices
    amrex::Vector<amrex::Real> m_z_ref;              // Vector of reals for height above MOST BC
    amrex::Vector<amrex::Real> c_interp;             // Interpolation factors
    const amrex::MultiFab* m_z_nd;                   // Ptr to terrain heigh coords
    amrex::Geometry m_geom;                          // Geometry
    const int m_policy;                              // Policy for averaging
    bool m_update_k;                                 // Flag to update 2D imf with k indices
};
#endif
