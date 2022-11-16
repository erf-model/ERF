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
    explicit MOSTAverage(amrex::Real& zref,
                         const amrex::Vector<amrex::Geometry>& geom,
                         const amrex::Vector<amrex::Vector<amrex::MultiFab>>& vars_old,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Theta_prim,
                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& z_phys_nd);
     MOSTAverage() = default;
    ~MOSTAverage()
     {
       for (int lev(0); lev<m_geom.size(); ++lev)
       {
         for (int iavg(0); iavg<navg; ++iavg) delete m_averages[lev][iavg];
         delete m_k_indx[lev];
       }
     };

     // Driver for the different average policies
     void compute_averages();

     // Get pointer to the 2D mf of averages
     const amrex::MultiFab* get_average(int lev, int comp) const { return m_averages[lev][comp]; };

protected:

    int nvar{3};                                                // 3 field for U/V/T
    int navg{4};                                                // 4 averages for U/V/T/Umag
    amrex::Real m_zref;                                         // Vector of reals for height above MOST BC
    const amrex::Vector<amrex::Geometry> m_geom;                      // Geometry at each level
    amrex::Vector<amrex::Vector<const amrex::MultiFab*>>  m_fields;   // Ptr to fields to be averaged
    amrex::Vector<amrex::Vector<amrex::MultiFab*>>  m_averages; // Ptr to 2D mf to hold averages
    amrex::Vector<amrex::MultiFab*> m_z_phys_nd;// Ptr to terrain heigh coords
    amrex::Vector<amrex::iMultiFab*> m_k_indx;                  // Ptr to 2D imf to hold k indices

};
#endif
