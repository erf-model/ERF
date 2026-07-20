#ifndef ERF_PROB_H
#define ERF_PROB_H

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>

struct ProbParm {
    // Gaussian terrain parameters
    amrex::Real hill_height    = 100.0;  // [m] positive = hill, negative = pit
    amrex::Real hill_halfwidth = 500.0;  // sigma [m]
    amrex::Real hill_cx        = -1.0;   // x-centre [-1 = domain centre]
    amrex::Real hill_cy        = -1.0;   // y-centre [-1 = domain centre]
};

class Problem {
public:
    Problem() {
        amrex::ParmParse pp("prob");
        pp.query("hill_height",    parms.hill_height);
        pp.query("hill_halfwidth", parms.hill_halfwidth);
        pp.query("hill_cx",        parms.hill_cx);
        pp.query("hill_cy",        parms.hill_cy);
    }

    void init_terrain_surface(
        const amrex::Geometry& geom,
        amrex::FArrayBox&      terrain_fab,
        const amrex::Real      /*time*/) const
    {
        const auto& dx  = geom.CellSizeArray();
        const auto& plo = geom.ProbLoArray();
        const auto& phi = geom.ProbHiArray();

        // Default centre to domain centre if not specified
        const amrex::Real cx = (parms.hill_cx < 0)
            ? 0.5*(plo[0]+phi[0]) : parms.hill_cx;
        const amrex::Real cy = (parms.hill_cy < 0)
            ? 0.5*(plo[1]+phi[1]) : parms.hill_cy;

        const amrex::Real H     = parms.hill_height;
        const amrex::Real sigma = parms.hill_halfwidth;
        const amrex::Real inv2s2 = 1.0 / (2.0*sigma*sigma);

        const amrex::Box& bx = terrain_fab.box();
        auto z_arr = terrain_fab.array();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int /*k*/) noexcept {
            // Nodal x,y coordinates
            amrex::Real x = plo[0] + i * dx[0];
            amrex::Real y = plo[1] + j * dx[1];
            amrex::Real r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy);
            z_arr(i,j,0) = H * std::exp(-r2 * inv2s2);
        });
    }

    ProbParm parms;
};

#endif
