#ifndef ERF_PERTURBATION_WIND_H_
#define ERF_PERTURBATION_WIND_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <ERF_NumericalConstants.H>

struct PerturbationWindSample {
    amrex::Real u = zero;
    amrex::Real v = zero;
    bool complete = false;
};

// Sample horizontal winds at one height AGL and average over the columns in
// region. Interpolate between cell centers; use the nearest center within the
// bottom/top half-cell. A height outside the available cells is not covered.
// All ranks participate. complete is true only when every column is covered.
PerturbationWindSample sample_perturbation_wind (
    const amrex::Box& region, amrex::Real height_agl,
    const amrex::Geometry& geom, const amrex::MultiFab& state,
    const amrex::MultiFab& xvel, const amrex::MultiFab& yvel,
    const amrex::MultiFab* z_cc, const amrex::MultiFab* z_nd);

#endif
