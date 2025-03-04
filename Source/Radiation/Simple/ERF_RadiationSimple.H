#ifndef ERF_RADIATION_SIMPLE_H
#define ERF_RADIATION_SIMPLE_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>

#include <ERF_DataStruct.H>

class RadiationSimple
{
public:
    RadiationSimple(const int& lev,
                    SolverChoice& sc) {
        // Check if we have a valid moisture model
        if (sc.moisture_type != MoistureType::None) { moist = true; }

        // Check if we have a moisture model with ice
        if (sc.moisture_type == MoistureType::SAM)  { ice = true; }
    };

    ~RadiationSimple() = default;

    void Init(const amrex::Geometry& geom,
              const amrex::BoxArray& ba,
              amrex::MultiFab* cons_in);

    void Run(const amrex::Real &dt,
             amrex::Geometry& geom,
             amrex::MultiFab* cons_in,
             amrex::MultiFab* qheating_rates,
             amrex::MultiFab* z_nd);

    const amrex::MultiFab *get_radlwdn() { return radlwdn; }
    const amrex::MultiFab *get_radqrlw() { return radqrlw; }

private:
    bool moist = false;
    bool ice = false;

    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;

    amrex::MultiFab deltaq, flux;
    amrex::MultiFab *radlwdn = nullptr, *radqrlw = nullptr;

};

#endif